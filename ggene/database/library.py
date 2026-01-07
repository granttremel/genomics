"""
GenomeLibrary: Unified management of local genomic data libraries.

Handles creation, indexing, metadata tracking, and retrieval of genomic datasets
including sequences (FASTA), features (BED/GTF), and arrays (NPZ).

Example usage:
    lib = GenomeLibrary("my_ncrna_project")

    # Streaming write - feeds data as it's generated
    with lib.writer("snrnas", "features") as w:
        for feature in generate_features():
            w.write(feature)  # accepts dicts or UFeatures

    # Retrieve as stream
    stream = lib.get_stream("snrnas")
    for feature in stream.stream("1", start=1000, end=50000):
        print(feature)

    # Quick sample
    samples = lib.sample("snrnas", n=5)
"""

import os
import json
import gzip
import subprocess
import shutil
import tempfile
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import (
    Optional, Dict, Any, List, Iterator, Union,
    Callable, Literal, TextIO, TYPE_CHECKING
)
from contextlib import contextmanager
import logging

from ggene.config import DEFAULT_LIBRARY

if TYPE_CHECKING:
    from ggene.database.annotations import TabularStream
    from ggene.database.ufeature import UFeature

logger = logging.getLogger(__name__)


# =============================================================================
# Data Types
# =============================================================================

DatasetType = Literal["features", "sequences", "arrays"]

@dataclass
class DatasetMeta:
    """Metadata for a single dataset in the library."""
    name: str
    path: str  # Relative to project root
    dtype: DatasetType

    # Stream reconstruction info
    stream_class: str = "TabularStream"
    columns: List[str] = field(default_factory=list)
    chr_format: str = "{chrstr}"
    feature_type: str = ""
    delimiter: str = "\t"

    # Provenance
    source: str = ""  # e.g., "ensembl_v110", "encode_ENCSR123ABC"
    organism: str = ""  # e.g., "homo_sapiens", "saccharomyces_cerevisiae"
    assembly: str = ""  # e.g., "GRCh38", "R64"
    description: str = ""

    # Tracking
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    modified: str = ""
    record_count: int = 0
    indexed: bool = False
    compressed: bool = False

    # Custom fields
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "DatasetMeta":
        # Handle extra fields gracefully
        known_fields = {f.name for f in cls.__dataclass_fields__.values()}
        extra = d.pop("extra", {})

        # Move unknown fields to extra
        for key in list(d.keys()):
            if key not in known_fields:
                extra[key] = d.pop(key)

        d["extra"] = extra
        return cls(**d)


@dataclass
class LibraryManifest:
    """Master manifest for a library project."""
    name: str
    version: str = "1.0"
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    description: str = ""
    datasets: Dict[str, DatasetMeta] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "created": self.created,
            "description": self.description,
            "datasets": {k: v.to_dict() for k, v in self.datasets.items()}
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "LibraryManifest":
        datasets = {
            k: DatasetMeta.from_dict(v)
            for k, v in d.get("datasets", {}).items()
        }
        return cls(
            name=d["name"],
            version=d.get("version", "1.0"),
            created=d.get("created", ""),
            description=d.get("description", ""),
            datasets=datasets
        )


# =============================================================================
# Library Writer - Streaming writes
# =============================================================================

class LibraryWriter:
    """
    Streaming writer for library datasets.

    Accepts dicts or UFeatures and writes them incrementally to avoid
    holding large datasets in memory. Handles sorting, compression,
    and indexing on close.

    Usage:
        with lib.writer("my_dataset", "features") as w:
            w.write({"chrom": "1", "start": 100, "end": 200, "name": "feat1"})
            w.write(some_ufeature)
    """

    # Default column order for BED-like output
    CORE_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]

    def __init__(
        self,
        output_path: Path,
        dtype: DatasetType,
        columns: Optional[List[str]] = None,
        header: bool = True,
        chr_format: str = "{chrstr}",
    ):
        self.output_path = Path(output_path)
        self.dtype = dtype
        self.columns = columns
        self.header = header
        self.chr_format = chr_format

        self._temp_file: Optional[Path] = None
        self._handle: Optional[TextIO] = None
        self._record_count = 0
        self._columns_written = False
        self._inferred_columns: List[str] = []

    def open(self):
        """Open the writer for streaming."""
        # Write to temp file first, then sort/compress on close
        self._temp_file = Path(tempfile.mktemp(suffix=".tsv"))
        self._handle = open(self._temp_file, "w")
        self._record_count = 0
        return self

    def write(self, record: Union[Dict[str, Any], "UFeature"]):
        """
        Write a single record.

        Args:
            record: Dict with field values, or UFeature object
        """
        
        if self._handle is None:
            raise RuntimeError("Writer not opened. Use 'with' context or call open()")
        # Convert UFeature to dict if needed
        if hasattr(record, "to_dict"):
            data = record.to_dict()
        elif hasattr(record, "__dict__"):
            data = self._ufeature_to_dict(record)
        else:
            data = record

        # Infer columns from first record if not specified
        if not self._columns_written:
            self._infer_and_write_header(data)
        
        # Write the record
        line = self._format_record(data)
        self._handle.write(line + "\n")
        self._record_count += 1

    def write_many(self, records: Iterator[Union[Dict[str, Any], "UFeature"]]):
        """Write multiple records from an iterator."""
        for record in records:
            self.write(record)

    def _ufeature_to_dict(self, uf: "UFeature") -> Dict[str, Any]:
        """Convert a UFeature to a dict."""
        data = {
            "chrom": uf.chrom,
            "start": uf.start,
            "end": uf.end,
            "feature_type": uf.feature_type,
            "name": getattr(uf, "name", ""),
            "score": getattr(uf, "score", "."),
            "strand": getattr(uf, "strand", "."),
        }
        # Add attributes
        if hasattr(uf, "attributes"):
            data.update(uf.attributes)
        return data

    def _infer_and_write_header(self, first_record: Dict[str, Any]):
        """Infer column order and optionally write header."""
        if self.columns:
            self._inferred_columns = self.columns
        else:
            # Put core columns first, then extras in sorted order
            core = [c for c in self.CORE_COLUMNS if c in first_record]
            extra = sorted(k for k in first_record.keys() if k not in self.CORE_COLUMNS)
            self._inferred_columns = core + extra

        if self.header and self.dtype == "features":
            header_line = "#" + "\t".join(self._inferred_columns)
            self._handle.write(header_line + "\n")

        self._columns_written = True

    def _format_record(self, data: Dict[str, Any]) -> str:
        """Format a record as a TSV line."""
        values = []
        for col in self._inferred_columns:
            val = data.get(col, ".")
            if val is None:
                val = "."
            elif isinstance(val, float):
                val = f"{val:.6g}"
            else:
                val = str(val)
            values.append(val)
        return "\t".join(values)

    def close(self) -> tuple[Path, int, List[str]]:
        """
        Close writer and finalize the file.

        Returns:
            Tuple of (final_path, record_count, columns)
        """
        if self._handle:
            self._handle.close()
            self._handle = None

        if self._temp_file and self._temp_file.exists():
            # Move temp file to final location
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(self._temp_file), str(self.output_path))

        return self.output_path, self._record_count, self._inferred_columns

    def __enter__(self):
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class FASTAWriter:
    """
    Streaming writer for FASTA sequences.

    Usage:
        with lib.fasta_writer("my_seqs") as w:
            w.write("seq1", "ATGCATGC...")
            w.write("seq2", "GCTAGCTA...", description="some gene")

    For indexed access to descriptions, consider encoding metadata in the
    sequence ID using delimiters, since pysam only exposes the first word
    of headers as the reference name:
        w.write("gene1|chr1|1000-2000|+", sequence)
    """

    def __init__(self, output_path: Path, line_width: int = 80):
        self.output_path = Path(output_path)
        self.line_width = line_width
        self._temp_file: Optional[Path] = None
        self._handle: Optional[TextIO] = None
        self._record_count = 0
        self._opened = False

    @property
    def is_open(self):
        return self._opened

    def open(self):
        if self._opened:
            return self
        self._temp_file = Path(tempfile.mktemp(suffix='.fa'))
        self._handle = open(self._temp_file, "w")
        self._record_count = 0
        self._opened = True
        return self

    def write(self, seq_id: str, sequence: str, description: str = ""):
        """Write a single FASTA record."""
        if self._handle is None:
            raise RuntimeError("Writer not opened")

        # Write header
        header = f">{seq_id}"
        if description:
            header += f" {description}"
        self._handle.write(header + "\n")

        # Write sequence with line wrapping
        for i in range(0, len(sequence), self.line_width):
            self._handle.write(sequence[i:i + self.line_width] + "\n")

        self._record_count += 1

    def write_feature(self, feature: "UFeature", sequence: str):
        """Write a feature as a FASTA record with genomic coordinates in header."""
        coords = f"{feature.chrom}:{feature.start}-{feature.end}"
        strand = getattr(feature, "strand", "+")
        name = getattr(feature, "name", feature.chrom)
        desc = f"{coords}({strand})"
        self.write(name, sequence, description=desc)

    def close(self) -> tuple[Path, int]:
        if self._handle:
            self._handle.close()
            self._handle = None

        if self._temp_file and self._temp_file.exists():
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(self._temp_file), str(self.output_path))

        return self.output_path, self._record_count

    def __enter__(self):
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


# =============================================================================
# Library FASTA Reader
# =============================================================================

@dataclass
class LibrarySequence:
    """A sequence from a library FASTA with parsed metadata."""
    seq_id: str
    sequence: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[str] = None
    description: str = ""

    @classmethod
    def parse_header(cls, header: str) -> Dict[str, Any]:
        """
        Parse a library FASTA header.

        Supports formats:
            >name chrom:start-end(strand) extra...
            >name|chrom|start|end|strand
        """
        parts = header.split(None, 1)
        seq_id = parts[0]
        desc = parts[1] if len(parts) > 1 else ""

        result = {"seq_id": seq_id, "description": desc}

        # Try to parse coordinates from description
        import re
        coord_match = re.search(r'(\w+):(\d+)-(\d+)\(([+-])\)', desc)
        if coord_match:
            result["chrom"] = coord_match.group(1)
            result["start"] = int(coord_match.group(2))
            result["end"] = int(coord_match.group(3))
            result["strand"] = coord_match.group(4)
        else:
            # Try pipe-delimited format in seq_id
            if '|' in seq_id:
                pipe_parts = seq_id.split('|')
                if len(pipe_parts) >= 4:
                    result["seq_id"] = pipe_parts[0]
                    result["chrom"] = pipe_parts[1]
                    try:
                        result["start"] = int(pipe_parts[2])
                        result["end"] = int(pipe_parts[3])
                    except ValueError:
                        pass
                    if len(pipe_parts) >= 5:
                        result["strand"] = pipe_parts[4]

        return result


class LibraryFASTAReader:
    """
    Reader for library FASTA files with metadata parsing.

    Unlike FASTAStream (for genomic DNA), this is designed for
    discontinuous sequence collections like ncRNA libraries.

    Usage:
        reader = LibraryFASTAReader("library.fa.gz")
        for seq in reader.iter_sequences():
            print(seq.seq_id, seq.chrom, seq.start)

        # Or fetch by name
        seq = reader.fetch("U1-snRNA")
    """

    def __init__(self, filepath: Union[str, Path]):
        self.filepath = Path(filepath)
        self._fasta = None
        self._header_cache: Dict[str, str] = {}

        if not self.filepath.exists():
            raise FileNotFoundError(f"FASTA not found: {filepath}")

        try:
            import pysam
            self._fasta = pysam.FastaFile(str(self.filepath))
        except Exception as e:
            logger.error(f"Failed to open FASTA: {e}")
            raise

    @property
    def references(self) -> List[str]:
        """List all sequence IDs."""
        return list(self._fasta.references)

    def fetch(self, seq_id: str) -> Optional[LibrarySequence]:
        """Fetch a sequence by ID."""
        if seq_id not in self._fasta.references:
            return None

        sequence = self._fasta.fetch(seq_id)
        header = self._get_header(seq_id)
        parsed = LibrarySequence.parse_header(header)

        return LibrarySequence(
            seq_id=parsed.get("seq_id", seq_id),
            sequence=sequence,
            chrom=parsed.get("chrom"),
            start=parsed.get("start"),
            end=parsed.get("end"),
            strand=parsed.get("strand"),
            description=parsed.get("description", ""),
        )

    def _get_header(self, seq_id: str) -> str:
        """
        Get full header line for a sequence.

        pysam doesn't expose descriptions directly, so we scan the file.
        Results are cached.
        """
        if seq_id in self._header_cache:
            return self._header_cache[seq_id]

        # Scan file for headers (only do this once)
        if not self._header_cache:
            self._scan_headers()

        return self._header_cache.get(seq_id, seq_id)

    def _scan_headers(self):
        """Scan file to extract full headers."""
        opener = gzip.open if str(self.filepath).endswith('.gz') else open
        mode = 'rt' if str(self.filepath).endswith('.gz') else 'r'

        with opener(self.filepath, mode) as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    seq_id = header.split()[0]
                    self._header_cache[seq_id] = header

    def iter_sequences(self, limit: Optional[int] = None) -> Iterator[LibrarySequence]:
        """Iterate over all sequences."""
        count = 0
        for ref in self._fasta.references:
            if limit and count >= limit:
                break
            seq = self.fetch(ref)
            if seq:
                yield seq
                count += 1

    def __len__(self):
        return len(self._fasta.references)

    def __repr__(self):
        return f"LibraryFASTAReader('{self.filepath}', n={len(self)})"


# =============================================================================
# Indexing Utilities
# =============================================================================

class IndexingError(Exception):
    """Raised when indexing operations fail."""
    pass


def _check_tool(tool: str) -> bool:
    """Check if a command-line tool is available."""
    return shutil.which(tool) is not None


def sort_bed(input_path: Path, output_path: Optional[Path] = None) -> Path:
    """
    Sort a BED/GTF file by chromosome and position.

    Args:
        input_path: Input file path
        output_path: Output path (default: overwrites input)

    Returns:
        Path to sorted file
    """
    if output_path is None:
        output_path = input_path

    temp_sorted = Path(tempfile.mktemp(suffix=".sorted"))

    try:
        # Use sort with natural chromosome ordering
        result = subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(input_path), "-o", str(temp_sorted)],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            raise IndexingError(f"sort failed: {result.stderr}")

        shutil.move(str(temp_sorted), str(output_path))
        return output_path

    except FileNotFoundError:
        raise IndexingError("'sort' command not found")
    finally:
        if temp_sorted.exists():
            temp_sorted.unlink()


def bgzip_file(filepath: Path, keep_original: bool = False, force: bool = True) -> Path:
    """
    Compress a file with bgzip.

    Args:
        filepath: File to compress
        keep_original: If True, keep the uncompressed file
        force: If True, overwrite existing .gz file without prompting

    Returns:
        Path to compressed file (.gz)
    """
    if not _check_tool("bgzip"):
        raise IndexingError("bgzip not found. Install samtools/htslib.")

    args = ["bgzip"]
    if force:
        args.append("-f")
    if keep_original:
        args.append("-k")
    args.append(str(filepath))

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise IndexingError(f"bgzip failed: {result.stderr}")

    return Path(str(filepath) + ".gz")

def trim_header(filepath:Path, output_path = None):
    
    if not output_path:
        tmp = True
        output_path = Path(f"/tmp/{filepath.stem}_tmp")
    
    with open(filepath, 'r') as f_r:
        with open(output_path, 'w') as f_w:
            for i, line in enumerate(f_r):
                if i==0:
                    continue
                f_w.write(line)
    
    if tmp:
        os.rename(output_path, filepath)
        output_path.unlink()
        output_path = filepath
        
    return output_path

def scrub_header(filepath:Path, output_path = None, header_start = "C"):
    
    if not output_path:
        tmp = True
        output_path = Path(f"/tmp/{filepath.stem}_tmp")
    
    with open(filepath, 'r') as f_r:
        with open(output_path, 'w') as f_w:
            for i, line in enumerate(f_r):
                if line.startswith(header_start):
                    continue
                else:
                    f_w.write(line)

    if tmp:
        os.rename(output_path, filepath)
        try:
            output_path.unlink()
        except:
            pass
        output_path = filepath
        
    return output_path
    

def gzip_file(filepath:Path, keep_original:bool = False, force: bool = True):
    
    if '.gz' in filepath.suffixes:
        decompress = True
        outpath = Path(str(filepath).removesuffix(".gz"))
    else:
        decompress = False
        outpath = Path(str(filepath) + ".gz")
        
    args = ["gzip"]
    # if force:
    #     args.append("-f")
    if keep_original:
        args.append("-k")
    if decompress:
        args.append("-d")
    args.append(str(filepath))
    
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise IndexingError(f"bgzip failed: {result.stderr}")
    
    return outpath


def tabix_index(filepath: Path, preset: str = "bed") -> Path:
    """
    Create a tabix index for a bgzipped file.

    Args:
        filepath: Path to bgzipped file
        preset: File format preset ('bed', 'gff', 'vcf')

    Returns:
        Path to index file (.tbi)
    """
    if not _check_tool("tabix"):
        raise IndexingError("tabix not found. Install samtools/htslib.")

    result = subprocess.run(
        ["tabix", "-p", preset, str(filepath)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise IndexingError(f"tabix failed: {result.stderr}")

    return Path(str(filepath) + ".tbi")


def index_fasta(filepath: Path) -> Path:
    """
    Create a FASTA index (.fai).

    Args:
        filepath: Path to FASTA file

    Returns:
        Path to index file (.fai)
    """
    if not _check_tool("samtools"):
        raise IndexingError("samtools not found. Install samtools.")

    result = subprocess.run(
        ["samtools", "faidx", str(filepath)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise IndexingError(f"samtools faidx failed: {result.stderr}")

    return Path(str(filepath) + ".fai")


def prepare_features_file(filepath: Path, preset: str = "bed", do_trim = False, do_scrub = False, header_start = "C") -> Path:
    """
    Full pipeline: sort, compress, and index a features file.

    Args:
        filepath: Input BED/GTF file
        preset: tabix preset

    Returns:
        Path to final compressed and indexed file
    """
    logger.info(f"Preparing features file: {filepath}")
    
    newpath = filepath
    if '.gz' in filepath.suffixes:
        newpath = gzip_file(newpath)
    
    if do_scrub:
        newpath = scrub_header(newpath, header_start = header_start)
    elif do_trim:
        newpath = trim_header(newpath)
        
    # Sort
    logger.info("  Sorting...")
    newpath = sort_bed(newpath)

    # Compress
    logger.info("  Compressing with bgzip...")
    gz_path = bgzip_file(newpath)

    # Index
    logger.info("  Creating tabix index...")
    tabix_index(gz_path, preset=preset)

    logger.info(f"  Done: {gz_path}")
    return gz_path


def prepare_fasta_file(filepath: Path) -> Path:
    """
    Compress and index a FASTA file.

    Args:
        filepath: Input FASTA file

    Returns:
        Path to final compressed and indexed file
    """
    logger.info(f"Preparing FASTA file: {filepath}")

    # Compress with bgzip
    logger.info("  Compressing with bgzip...")
    gz_path = bgzip_file(filepath)
    
    # Index (samtools faidx works on bgzipped FASTA)
    logger.info("  Creating FASTA index...")
    index_fasta(gz_path)

    logger.info(f"  Done: {gz_path}")
    return gz_path


# =============================================================================
# Main GenomeLibrary Class
# =============================================================================

class GenomeLibrary:
    """
    Manages a local library of genomic datasets.

    Provides unified interface for:
    - Creating and organizing datasets by project
    - Streaming writes to avoid memory issues
    - Automatic sorting, compression, and indexing
    - Metadata tracking via manifest
    - Retrieval as TabularStream or direct queries

    Directory structure:
        library_root/
        ├── project_name/
        │   ├── manifest.json
        │   ├── sequences/
        │   │   └── dataset.fa.gz
        │   ├── features/
        │   │   └── dataset.bed.gz
        │   └── arrays/
        │       └── dataset.npz

    Example:
        lib = GenomeLibrary("ncrna_analysis")

        # Write data
        with lib.writer("snrnas", "features", organism="human") as w:
            for f in generate_snrnas():
                w.write(f)

        # Read back
        stream = lib.get_stream("snrnas")
        for feature in stream.stream("chr1"):
            print(feature)
    """

    SUBDIRS = ["sequences", "features", "arrays"]

    def __init__(
        self,
        project: str,
        library_root: Optional[Path] = None,
        create: bool = True
    ):
        """
        Initialize or open a library project.

        Args:
            project: Project name (creates subdirectory)
            library_root: Root library directory (default: config DEFAULT_LIBRARY)
            create: If True, create project directory if it doesn't exist
        """
        self.library_root = Path(library_root) if library_root else DEFAULT_LIBRARY
        self.project = project
        self.project_dir = self.library_root / project
        self.manifest_path = self.project_dir / "manifest.json"

        if create:
            self._init_project()

        self.manifest = self._load_manifest()

    def _init_project(self):
        """Initialize project directory structure."""
        if not self.project_dir.exists():
            self.project_dir.mkdir(parents=True)
            for subdir in self.SUBDIRS:
                (self.project_dir / subdir).mkdir()

            # Create empty manifest
            manifest = LibraryManifest(name=self.project)
            self._save_manifest(manifest)
            logger.info(f"Created new library project: {self.project_dir}")

    def _load_manifest(self) -> LibraryManifest:
        """Load manifest from disk."""
        if self.manifest_path.exists():
            with open(self.manifest_path) as f:
                data = json.load(f)
            return LibraryManifest.from_dict(data)
        return LibraryManifest(name=self.project)

    def _save_manifest(self, manifest: Optional[LibraryManifest] = None):
        """Save manifest to disk."""
        if manifest is None:
            manifest = self.manifest
        with open(self.manifest_path, "w") as f:
            json.dump(manifest.to_dict(), f, indent=2)

    def _get_dataset_path(self, name: str, dtype: DatasetType) -> Path:
        """Get the path for a dataset file."""
        ext_map = {
            "features": ".bed",
            "sequences": ".fa",
            "arrays": ".npz"
        }
        ext = ext_map.get(dtype, ".dat")
        return self.project_dir / dtype / f"{name}{ext}"

    # -------------------------------------------------------------------------
    # Writing
    # -------------------------------------------------------------------------

    @contextmanager
    def writer(
        self,
        name: str,
        dtype: DatasetType = "features",
        columns: Optional[List[str]] = None,
        chr_format: str = "{chrstr}",
        feature_type: str = "",
        source: str = "",
        organism: str = "Hsa",
        assembly: str = "",
        description: str = "",
        index: bool = True,
        overwrite: bool = False,
        **extra_meta
    ) -> Iterator[LibraryWriter]:
        """
        Context manager for streaming writes to a dataset.

        Args:
            name: Dataset name
            dtype: Type of data ("features", "sequences", "arrays")
            columns: Column names (inferred from first record if not provided)
            chr_format: Chromosome format string for TabularStream
            feature_type: Feature type for TabularStream
            source: Data source (e.g., "ensembl_v110")
            organism: Organism name
            assembly: Genome assembly
            description: Human-readable description
            index: If True, sort/compress/index after writing
            overwrite: If True, overwrite existing dataset
            **extra_meta: Additional metadata fields

        Yields:
            LibraryWriter for streaming writes

        Example:
            with lib.writer("my_data", "features", organism="human") as w:
                w.write({"chrom": "1", "start": 100, "end": 200})
        """
        if name in self.manifest.datasets and not overwrite:
            raise ValueError(f"Dataset '{name}' already exists. Use overwrite=True.")

        output_path = self._get_dataset_path(name, dtype)

        writer = LibraryWriter(
            output_path=output_path,
            dtype=dtype,
            columns=columns,
            chr_format=chr_format,
        )

        try:
            writer.open()
            yield writer
        finally:
            final_path, record_count, inferred_columns = writer.close()
            
            

            # Process and index
            if index and dtype == "features" and record_count > 0:
                try:
                    final_path = prepare_features_file(final_path)
                    indexed = True
                    compressed = True
                except IndexingError as e:
                    logger.warning(f"Indexing failed: {e}")
                    indexed = False
                    compressed = False
            else:
                indexed = False
                compressed = False

            # Update manifest
            rel_path = str(final_path.relative_to(self.project_dir))
            meta = DatasetMeta(
                name=name,
                path=rel_path,
                dtype=dtype,
                columns=inferred_columns,
                chr_format=chr_format,
                feature_type=feature_type,
                source=source,
                organism=organism,
                assembly=assembly,
                description=description,
                record_count=record_count,
                indexed=indexed,
                compressed=compressed,
                extra=extra_meta,
            )
            self.manifest.datasets[name] = meta
            self._save_manifest()

            logger.info(f"Wrote dataset '{name}': {record_count} records")

    # -------------------------------------------------------------------------
    # Simple Writer Factory (for batch/dict patterns)
    # -------------------------------------------------------------------------

    def create_fasta_writer(
        self,
        name: str,
        description: str = "",
        source: str = "",
        organism: str = "Hsa",
        assembly: str = "GRCh38",
        overwrite: bool = False,
        **extra_meta
    ) -> FASTAWriter:
        """
        Create a FASTAWriter directly (not as context manager).

        Use this when you need multiple writers in a dict pattern.
        Call finalize_fasta_writer() when done writing.

        Example:
            writers = {}
            for rtype in rna_types:
                writers[rtype] = lib.create_fasta_writer(f"ncrna_{rtype}")

            for feature in features:
                writers[feature.type].write(...)

            for name, writer in writers.items():
                lib.finalize_fasta_writer(name, writer)
        """
        if name in self.manifest.datasets and not overwrite:
            raise ValueError(f"Dataset '{name}' already exists. Use overwrite=True.")

        output_path = self._get_dataset_path(name, "sequences")

        # Store metadata for finalization
        self._pending_meta = getattr(self, '_pending_meta', {})
        self._pending_meta[name] = {
            "source": source,
            "organism": organism,
            "assembly": assembly,
            "description": description,
            "extra": extra_meta,
        }

        writer = FASTAWriter(output_path=output_path)
        writer.open()
        return writer

    def finalize_fasta_writer(
        self,
        name: str,
        writer: FASTAWriter,
        index: bool = True
    ):
        """
        Finalize a FASTAWriter created with create_fasta_writer().

        Closes the writer, compresses, indexes, and updates manifest.
        """
        
        final_path, record_count = writer.close()
        
        # Get stored metadata
        meta_info = getattr(self, '_pending_meta', {}).pop(name, {})

        # Compress and index
        if index and record_count > 0:
            try:
                final_path = prepare_fasta_file(final_path)
                indexed = True
                compressed = True
            except IndexingError as e:
                logger.warning(f"Indexing failed: {e}")
                indexed = False
                compressed = False
        else:
            indexed = False
            compressed = False

        # Update manifest
        
        rel_path = str(final_path.relative_to(self.project_dir))
        meta = DatasetMeta(
            name=name,
            path=rel_path,
            dtype="sequences",
            source=meta_info.get("source", ""),
            organism=meta_info.get("organism", ""),
            assembly=meta_info.get("assembly", ""),
            description=meta_info.get("description", ""),
            record_count=record_count,
            indexed=indexed,
            compressed=compressed,
            extra=meta_info.get("extra", {}),
        )
        self.manifest.datasets[name] = meta
        self._save_manifest()

        logger.info(f"Finalized FASTA dataset '{name}': {record_count} sequences")

    def create_features_writer(
        self,
        name: str,
        columns: Optional[List[str]] = None,
        chr_format: str = "{chrstr}",
        feature_type: str = "",
        source: str = "",
        organism: str = "Hsa",
        assembly: str = "",
        description: str = "",
        overwrite: bool = False,
        **extra_meta
    ) -> LibraryWriter:
        """
        Create a LibraryWriter directly (not as context manager).

        Use this when you need multiple writers in a dict pattern.
        Call finalize_features_writer() when done writing.

        Example:
            writers = {}
            for rtype in rna_types:
                writers[rtype] = lib.create_features_writer(f"ncrna_{rtype}")

            for feature in features:
                writers[feature.type].write(...)

            for name, writer in writers.items():
                lib.finalize_features_writer(name, writer)
        """
        if name in self.manifest.datasets and not overwrite:
            raise ValueError(f"Dataset '{name}' already exists. Use overwrite=True.")

        output_path = self._get_dataset_path(name, "features")

        # Store metadata for finalization
        self._pending_meta = getattr(self, '_pending_meta', {})
        self._pending_meta[name] = {
            "chr_format": chr_format,
            "feature_type": feature_type,
            "source": source,
            "organism": organism,
            "assembly": assembly,
            "description": description,
            "columns": columns,
            "extra": extra_meta,
        }

        writer = LibraryWriter(
            output_path=output_path,
            dtype="features",
            columns=columns,
            chr_format=chr_format,
        )
        writer.open()
        return writer

    def finalize_features_writer(
        self,
        name: str,
        writer: LibraryWriter,
        index: bool = True
    ):
        """
        Finalize a LibraryWriter created with create_features_writer().

        Closes the writer, sorts, compresses, indexes, and updates manifest.
        """
        final_path, record_count, inferred_columns = writer.close()

        # logger.debug(f"feature writer inferred columns: {inferred_columns}")
        # logger.debug(f"feature_writer path {final_path}, rcord count {record_count}")

        # Get stored metadata
        meta_info = getattr(self, '_pending_meta', {}).pop(name, {})

        # Process and index
        if index and record_count > 0:
            try:
                final_path = prepare_features_file(final_path)
                indexed = True
                compressed = True
            except IndexingError as e:
                logger.warning(f"Indexing failed: {e}")
                indexed = False
                compressed = False
        else:
            indexed = False
            compressed = False

        # Update manifest
        rel_path = str(final_path.relative_to(self.project_dir))
        meta = DatasetMeta(
            name=name,
            path=rel_path,
            dtype="features",
            columns=inferred_columns,
            chr_format=meta_info.get("chr_format", "{chrstr}"),
            feature_type=meta_info.get("feature_type", ""),
            source=meta_info.get("source", ""),
            organism=meta_info.get("organism", ""),
            assembly=meta_info.get("assembly", ""),
            description=meta_info.get("description", ""),
            record_count=record_count,
            indexed=indexed,
            compressed=compressed,
            extra=meta_info.get("extra", {}),
        )
        self.manifest.datasets[name] = meta
        self._save_manifest()

        logger.info(f"Finalized features dataset '{name}': {record_count} records")

    @contextmanager
    def fasta_writer(
        self,
        name: str,
        source: str = "",
        organism: str = "Hsa",
        assembly: str = "",
        description: str = "",
        index: bool = True,
        overwrite: bool = False,
        **extra_meta
    ) -> Iterator[FASTAWriter]:
        """
        Context manager for streaming FASTA writes.

        Args:
            name: Dataset name
            source: Data source
            organism: Organism name
            assembly: Genome assembly
            description: Human-readable description
            index: If True, compress and index after writing
            overwrite: If True, overwrite existing dataset
            **extra_meta: Additional metadata fields

        Yields:
            FASTAWriter for streaming writes
        """
        if name in self.manifest.datasets and not overwrite:
            raise ValueError(f"Dataset '{name}' already exists. Use overwrite=True.")

        output_path = self._get_dataset_path(name, "sequences")

        writer = FASTAWriter(output_path=output_path)

        try:
            writer.open()
            yield writer
        finally:
            final_path, record_count = writer.close()

            # Compress and index
            if index and record_count > 0:
                try:
                    final_path = prepare_fasta_file(final_path)
                    indexed = True
                    compressed = True
                except IndexingError as e:
                    logger.warning(f"Indexing failed: {e}")
                    indexed = False
                    compressed = False
            else:
                indexed = False
                compressed = False

            # Update manifest
            rel_path = str(final_path.relative_to(self.project_dir))
            meta = DatasetMeta(
                name=name,
                path=rel_path,
                dtype="sequences",
                source=source,
                organism=organism,
                assembly=assembly,
                description=description,
                record_count=record_count,
                indexed=indexed,
                compressed=compressed,
                extra=extra_meta,
            )
            self.manifest.datasets[name] = meta
            self._save_manifest()

            logger.info(f"Wrote FASTA dataset '{name}': {record_count} sequences")

    # -------------------------------------------------------------------------
    # Reading
    # -------------------------------------------------------------------------

    def get_stream(self, name: str) -> "TabularStream":
        """
        Get a TabularStream for a features dataset.

        Args:
            name: Dataset name

        Returns:
            Configured TabularStream ready for querying
        """
        from ggene.database.annotations import TabularStream, ColumnSpec

        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found in library")

        meta = self.manifest.datasets[name]
        if meta.dtype != "features":
            raise ValueError(f"Dataset '{name}' is type '{meta.dtype}', not 'features'")

        filepath = self.project_dir / meta.path

        # Build column specs - convert all to ColumnSpec to avoid parse_line bug
        # with string columns in the secondary loop
        columns = []
        for col in meta.columns:
            if col in TabularStream.core_columns:
                # Get the actual ColumnSpec from core_columns
                columns.append(TabularStream.core_columns[col])
            else:
                # Custom column - create a ColumnSpec
                columns.append(ColumnSpec(col, str, "", target="attribute"))

        # Create a dynamic stream class with stored configuration
        class LibraryStream(TabularStream):
            comment_chars = ('#',)  # Our files use # for header

        LibraryStream.columns = columns
        LibraryStream.chr_format = meta.chr_format
        LibraryStream.feature_type = meta.feature_type

        return LibraryStream(str(filepath), source_name=name)

    def get_fasta(self, name: str) -> LibraryFASTAReader:
        """
        Get a LibraryFASTAReader for a sequences dataset.

        Unlike FASTAStream (for genomic DNA), this is designed for
        library collections with discontinuous sequences.

        Args:
            name: Dataset name

        Returns:
            LibraryFASTAReader with header parsing support
        """
        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found in library")

        meta = self.manifest.datasets[name]
        if meta.dtype != "sequences":
            raise ValueError(f"Dataset '{name}' is type '{meta.dtype}', not 'sequences'")

        filepath = self.project_dir / meta.path
        return LibraryFASTAReader(filepath)

    def get_genomic_fasta(self, name: str):
        """
        Get a FASTAStream for genomic DNA access (coordinate-based queries).

        Use get_fasta() instead for library sequences.
        """
        from ggene.database.sequences import FASTAStream

        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found in library")

        meta = self.manifest.datasets[name]
        filepath = self.project_dir / meta.path
        return FASTAStream(str(filepath))

    def sample(
        self,
        name: str,
        n: int = 5,
        chrom: Optional[str] = None
    ) -> List[Any]:
        """
        Get a quick sample of records from a dataset.

        Useful for inspecting data without loading everything.

        Args:
            name: Dataset name
            n: Number of records to return
            chrom: Optional chromosome filter

        Returns:
            List of records (UFeatures for features, sequences for FASTA)
        """
        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found in library")

        meta = self.manifest.datasets[name]

        if meta.dtype == "features":
            stream = self.get_stream(name)
            results = []
            for feature in stream.stream(chrom=chrom, start=1):
                results.append(feature)
                if len(results) >= n:
                    break
            return results

        elif meta.dtype == "sequences":
            # For FASTA, read first n records from raw file
            filepath = self.project_dir / meta.path
            results = []

            opener = gzip.open if filepath.suffix == ".gz" else open
            mode = "rt" if filepath.suffix == ".gz" else "r"

            with opener(filepath, mode) as f:
                current_id = None
                current_seq = []

                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_id and len(results) < n:
                            results.append((current_id, "".join(current_seq)))
                        current_id = line[1:].split()[0]
                        current_seq = []
                        if len(results) >= n:
                            break
                    else:
                        current_seq.append(line)

                # Don't forget last sequence
                if current_id and len(results) < n:
                    results.append((current_id, "".join(current_seq)))

            return results

        else:
            raise ValueError(f"Sampling not supported for dtype '{meta.dtype}'")

    # -------------------------------------------------------------------------
    # Metadata & Info
    # -------------------------------------------------------------------------

    def list_datasets(self) -> List[str]:
        """List all dataset names in this project."""
        return list(self.manifest.datasets.keys())

    def get_meta(self, name: str) -> DatasetMeta:
        """Get metadata for a dataset."""
        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found")
        return self.manifest.datasets[name]

    def update_meta(self, name: str, **updates):
        """
        Update metadata for a dataset.

        Args:
            name: Dataset name
            **updates: Fields to update
        """
        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found")

        meta = self.manifest.datasets[name]

        for key, value in updates.items():
            if hasattr(meta, key):
                setattr(meta, key, value)
            else:
                meta.extra[key] = value

        meta.modified = datetime.now().isoformat()
        self._save_manifest()

    def describe(self, name: Optional[str] = None) -> str:
        """
        Get a human-readable description of a dataset or the whole library.

        Args:
            name: Dataset name, or None for library summary

        Returns:
            Formatted description string
        """
        if name:
            meta = self.get_meta(name)
            lines = [
                f"Dataset: {meta.name}",
                f"  Type: {meta.dtype}",
                f"  Path: {meta.path}",
                f"  Records: {meta.record_count:,}",
                f"  Organism: {meta.organism or 'not specified'}",
                f"  Assembly: {meta.assembly or 'not specified'}",
                f"  Source: {meta.source or 'not specified'}",
                f"  Created: {meta.created}",
                f"  Indexed: {meta.indexed}",
            ]
            if meta.description:
                lines.append(f"  Description: {meta.description}")
            if meta.columns:
                lines.append(f"  Columns: {', '.join(meta.columns)}")
            if meta.extra:
                lines.append(f"  Extra: {meta.extra}")
            return "\n".join(lines)

        else:
            lines = [
                f"Library: {self.manifest.name}",
                f"  Location: {self.project_dir}",
                f"  Datasets: {len(self.manifest.datasets)}",
                ""
            ]
            for ds_name, meta in self.manifest.datasets.items():
                lines.append(f"  [{meta.dtype}] {ds_name}: {meta.record_count:,} records")
            return "\n".join(lines)

    def delete_dataset(self, name: str, confirm: bool = False):
        """
        Delete a dataset and its files.

        Args:
            name: Dataset name
            confirm: Must be True to actually delete
        """
        if not confirm:
            raise ValueError("Must pass confirm=True to delete")

        if name not in self.manifest.datasets:
            raise KeyError(f"Dataset '{name}' not found")

        meta = self.manifest.datasets[name]
        filepath = self.project_dir / meta.path

        # Delete main file and any index files
        for suffix in ["", ".tbi", ".fai", ".gzi"]:
            p = Path(str(filepath) + suffix)
            if p.exists():
                p.unlink()
                logger.info(f"Deleted: {p}")

        del self.manifest.datasets[name]
        self._save_manifest()
        logger.info(f"Deleted dataset '{name}' from manifest")

    def __repr__(self):
        return f"GenomeLibrary('{self.project}', datasets={len(self.manifest.datasets)})"


# =============================================================================
# Convenience functions
# =============================================================================

def list_projects(library_root: Optional[Path] = None) -> List[str]:
    """List all projects in the library root."""
    root = Path(library_root) if library_root else DEFAULT_LIBRARY
    if not root.exists():
        return []

    projects = []
    for p in root.iterdir():
        if p.is_dir() and (p / "manifest.json").exists():
            projects.append(p.name)
    return projects


def open_library(project: str, library_root: Optional[Path] = None) -> GenomeLibrary:
    """Open an existing library project (does not create if missing)."""
    return GenomeLibrary(project, library_root=library_root, create=False)

# examples



def build_rna_library(gm, chromes=[], project_name="ncrna_library"):
    """
    Build an ncRNA library using the GenomeLibrary system.

    Creates indexed BED files for each RNA type with sequences,
    merged across all chromosomes into single files.

    Args:
        gm: GenomeManager instance
        chromes: List of chromosomes to process (empty = all)
        project_name: Name for the library project

    Returns:
        GenomeLibrary instance with populated datasets
    """
    # from ggene.database.library import GenomeLibrary

    from ggene.genome import ncrna
    rna_types = ncrna.rna_types
    lib = GenomeLibrary(project_name)

    gstr = gm.annotations.streams.get("genes")

    # Create a writer for each RNA type using factory pattern
    # (avoids contextmanager issues with dict-based batch writing)
    writers = {}
    writer_names = {}
    for rtp in rna_types:
        name = f"ncrna_{rtp}"
        writers[rtp] = lib.create_features_writer(
            name=name,
            columns=["chrom", "start", "end", "gene_id", "name", "strand", "seq"],
            feature_type=rtp,
            organism="homo_sapiens",
            source="ensembl",
            description=f"{rtp} genes with sequences",
            overwrite=True,
        )
        writer_names[rtp] = name

    ns = 0
    try:
        for chrom in gm.iter_chromes():
            if chromes and chrom not in chromes:
                continue

            print(f"Processing chromosome {chrom}")

            for f in gstr.stream(chrom, start=1):
                if f.gene_biotype not in rna_types:
                    continue

                fseq = gm.get_sequence(f.chrom, f.start, f.end)

                # Use gene_id as fallback if name is empty/None
                feat_name = f.name if f.name else f.gene_id
                if not feat_name:
                    feat_name = f"{f.chrom}_{f.start}_{f.end}"

                # Write directly to the appropriate writer
                writers[f.gene_biotype].write({
                    "chrom": f.chrom,
                    "start": f.start,
                    "end": f.end,
                    "gene_id": f.gene_id,
                    "name": feat_name,
                    "strand": f.strand,
                    "seq": fseq,
                })

                ns += 1
                if ns % 1000 == 0:
                    print(f"  collected {ns} sequences")

    finally:
        # Finalize all writers (closes, sorts, compresses, indexes, updates manifest)
        for rtp, writer in writers.items():
            print(f"Finalizing {writer_names[rtp]}...")
            lib.finalize_features_writer(writer_names[rtp], writer)

    print(f"\nDone! Built library with {ns} total sequences")
    print(lib.describe())

    return lib


def build_rna_fasta_library(gm, chromes=[], project_name="ncrna_sequences", overwrite = False):
    """
    Build an ncRNA FASTA library using the GenomeLibrary system.

    Creates indexed FASTA files for each RNA type.

    Args:
        gm: GenomeManager instance
        chromes: List of chromosomes to process (empty = all)
        project_name: Name for the library project

    Returns:
        GenomeLibrary instance with populated datasets
    """
    from ggene.genome import ncrna
    # from ggene.database.library import GenomeLibrary

    rna_types = ncrna.rna_types
    lib = GenomeLibrary(project_name)

    gstr = gm.annotations.streams.get("genes")

    # Create a FASTA writer for each RNA type using factory pattern
    # (avoids contextmanager issues with dict-based batch writing)
    writers = {}
    writer_names = {}
    for rtp in rna_types:
        name = f"ncrna_{rtp}"
        writers[rtp] = lib.create_fasta_writer(
            name=name,
            organism="homo_sapiens",
            source="ensembl",
            description=f"{rtp} sequences",
            overwrite=overwrite,
        )
        writer_names[rtp] = name

    ns = 0
    try:
        for chrom in gm.iter_chromes():
            if chromes and chrom not in chromes:
                continue

            print(f"Processing chromosome {chrom}")

            for f in gstr.stream(chrom, start=1):
                if f.gene_biotype not in rna_types:
                    continue

                fseq = gm.get_sequence(f.chrom, f.start, f.end)
                desc = f"{f.chrom}:{f.start}-{f.end}({f.strand}) {f.gene_id}"

                # Use gene_id as fallback if name is empty/None
                seq_id = f.name if f.name else f.gene_id
                if not seq_id:
                    seq_id = f"{f.chrom}_{f.start}_{f.end}"

                writer = writers[f.gene_biotype]
                writer.write(seq_id, fseq, description=desc)

                ns += 1
                if ns % 1000 == 0:
                    print(f"  collected {ns} sequences")

    finally:
        # Finalize all writers (closes, compresses, indexes, updates manifest)
        for rtp, writer in writers.items():
            print(f"Finalizing {writer_names[rtp]}...")
            lib.finalize_fasta_writer(writer_names[rtp], writer)

    print(f"\nDone! Built FASTA library with {ns} total sequences")
    print(lib.describe())

    return lib



