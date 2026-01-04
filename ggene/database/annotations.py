"""Unified streaming interface for genomic annotations.

This module provides a common interface to stream and merge annotations from
multiple genomic databases (GTF, VCF, BED, JASPAR, etc.) into a unified
Feature-based representation.

IMPORTANT: All data files should be bgzipped and tabix-indexed for performance.
Without indexing, queries will be extremely slow on large files.
"""

import heapq
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Iterator, Dict, Any, List, Optional, Tuple, Union, Callable, TypeVar
import gzip
import json
from pathlib import Path
import logging
import pysam
import os
import warnings
import gzip

from ggene.database.ufeature import UFeature

logger = logging.getLogger(__name__)

T = TypeVar('T')

logger.setLevel("CRITICAL")
# logger.setLevel(logging.DEBUG)

chr_lens = {'1': 248937043, '10': 133778498, '11': 135075908, '12': 133238549, '13': 114346637, '14': 106879812, '15': 101979093, '16': 90222678, '17': 83240391, '18': 80247514, '19': 58599303, '2': 242175634, '20': 64327972, '21': 46691226, '22': 50799123, '3': 198228376, '4': 190195978, '5': 181472430, '6': 170745977, '7': 159233377, '8': 145066516, '9': 138320835, 'MT': 16023, 'X': 156027877, 'Y': 57214397}

@dataclass
class ColumnSpec:
    """Specification for parsing a column from tabular data.

    Attributes:
        name: Column/attribute name
        type_: Type converter function (int, float, str, etc.)
        default: Default value if parsing fails or value is missing
        target: Where to store - 'core' for UFeature fields, 'attribute' for attributes dict
        formatter: Optional custom parser function (receives raw string, returns parsed value)
    """
    name: str
    type_: Callable[[str], T] = str
    default: Any = None
    target: str = "attribute"  # 'core' or 'attribute'
    formatter: Optional[Callable[[str], Any]] = None

    def parse(self, value: str) -> Any:
        """Parse a string value according to this spec."""
        if value is None or value == "" or value == ".":
            return self.default

        try:
            if self.formatter:
                return self.formatter(value)
            return self.type_(value)
        except (ValueError, TypeError):
            return self.default


@dataclass
class DerivedSpec:
    """Specification for a derived/computed attribute.

    The compute function receives the UFeature after initial construction
    and returns the computed value.

    Attributes:
        name: Attribute name to store the computed value
        compute: Function that takes UFeature and returns the derived value
        default: Default value if computation fails
    """
    name: str
    compute: Callable[['UFeature'], Any]
    default: Any = None

    def evaluate(self, feature: 'UFeature') -> Any:
        """Compute the derived value for a feature."""
        try:
            return self.compute(feature)
        except Exception:
            return self.default

class AnnotationStream(ABC):
    """Abstract base class for annotation streams."""
    
    def __init__(self, source_name: str):
        self.source_name = source_name
        self._current_chrom = None
        self._buffer = []
        
    @abstractmethod
    def parse_line(self, line: str) -> Optional[UFeature]:
        """Parse a single line into a UFeature."""
        pass
    
    @abstractmethod
    def stream(self, chrom: Optional[str] = None, 
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UFeature]:
        """Stream features, optionally filtered by region."""
        pass
    
    def query_range(self, chrom: str, start: int, end: int) -> List[UFeature]:
        """Query features in a specific range."""
        features = []
        for feature in self.stream(chrom, start, end):
            if feature.overlaps(start, end):
                features.append(feature)
        return features


class TabularStream(AnnotationStream):
    """Base class for tabular file streams using ColumnSpec pattern.

    Subclasses define columns and derived attributes as class attributes,
    and this base class handles the parsing automatically.

    Example:
        class MyStream(TabularStream):
            columns = [
                ColumnSpec("chrom", str, "", "core"),
                ColumnSpec("start", int, -1, "core"),
                ColumnSpec("end", int, -1, "core"),
                ColumnSpec("score", float, None, "attribute"),
            ]
            derived = [
                DerivedSpec("length", lambda uf: uf.end - uf.start),
            ]
    """
    
    core_columns: Dict[str, ColumnSpec] = {
        "chrom":ColumnSpec("chrom", str, "", "core", formatter = lambda s: TabularStream.format_chrom(s)[1]),
        "start":ColumnSpec("start", int, -1, "core"),
        "end":ColumnSpec("end", int, -1, "core"),
        "name":ColumnSpec("name", str, "", "core"),
        "score":ColumnSpec("score", float, None, "core"),
        "strand":ColumnSpec("strand", str, "", "core"),
        "id":ColumnSpec("id", str, "", "core"),
    }
    
    # Subclasses override these
    columns: List[Union[ColumnSpec, str]] = []
    derived: List[DerivedSpec] = []
    feature_type: str = "region"
    comment_chars: Tuple[str, ...] = ('#', 'track')
    delimiter: str = '\t'
    chr_format: str = "{chrstr}"  # or "chr{chrstr}"

    fe_delimiter:str = ';'
    fe_kv_delimiter:str = " "

    def __init__(self, filepath: str, source_name: str = "Tabular"):
        super().__init__(source_name)
        self.filepath = Path(filepath).absolute()

        # Setup tabix for indexed access
        self.tabix = None
        if self.filepath.suffix == '.gz':
            index_file = Path(str(self.filepath) + '.tbi')
            if index_file.exists():
                try:
                    self.tabix = pysam.TabixFile(str(self.filepath))
                    logger.info(f"Using indexed access for {self.filepath}")
                except Exception as e:
                    logger.warning(f"Failed to open tabix index: {e}")

        if not self.tabix and self.filepath.suffix == '.gz':
            logger.warning(f"No index found for {self.filepath}.")
            
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UFeature]:
        """Stream features using indexed access."""
        if not self.tabix:
            logger.warning(f"No tabix index for {self.filepath}")
            return
        
        # clsname = type(self).__name__
        # if clsname == "VariantStream":
        #     print("entered stream")
        chr_formatted, chr_raw = self.format_chrom(chrom) if chrom else ("", "")

        query_start = int(start) if start else 1
        query_end = int(end) if end else chr_lens.get(chr_raw, 999999999)

        try:
            for line in self.tabix.fetch(chr_formatted, start=query_start, end=query_end):
                feature = self.parse_line(line)
                if self.validate_feature(feature):
                    
                    # Apply region filter
                    if chrom and feature.chrom != chr_raw:
                        continue
                    if start and feature.end < start:
                        continue
                    if end and feature.start > end:
                        continue
                    
                    
                    yield feature
        except Exception as e:
            logger.error(f"Error fetching {chr_formatted}:{query_start}-{query_end}: {e}")

    def validate_feature(self, feature):
        return bool(feature)

    def preparse_line(self, line:str):
        parts = line.strip().split(self.delimiter)
        return parts, {}

    def parse_line(self, line: str) -> Optional[UFeature]:
        """Parse a line using column specs."""
        # Skip comments
        if any(line.startswith(c) for c in self.comment_chars):
            return None
        
        parts, fe_dict = self.preparse_line(line)
        
        uf_cols = []
        
        # Build data dict from parts
        data = {}
        for i, col in enumerate(self.columns):
            if isinstance(col, str):
                col = self.core_columns.get(col)
            if col is None:
                continue
            
            uf_cols.append(col)
            try:
                if i < len(parts):
                    v = col.parse(parts[i])
                    data[col.name] = v
                else:
                    data[col.name] = ""
                    break
            except Exception as e:
                print(f"exception on column {col}: {str(e)}")
        
        for j, col in enumerate(self.columns[i:]):
            
            if col and col.name in fe_dict:
                data[col.name] = fe_dict[col.name]
                uf_cols.append(col)
            
        # Add feature_type and source
        # Add feature_type if not already set from columns
        if 'feature_type' not in data:
            data['feature_type'] = self.feature_type
        data['source'] = self.source_name

        # Create UFeature with flat data and compute derived attributes
        return UFeature.from_parsed(data, columns=uf_cols, derived=self.derived, raw_line = line)


    def format_group_entry(self, group_entry:str):
        
        fparts = group_entry.split(self.fe_delimiter)
        
        feparts = {}
        
        for p in fparts:
            kv = p.strip().split(self.fe_kv_delimiter)
            if len(kv) > 2:
                k = kv[0].strip()
                v = " ".join(kv[1:])
            elif len(kv) == 2:
                k, v = kv
            else:
                k = v = ""
            
            if k and v:
                
                k = k.strip()
                v = v.strip().strip('"')
                
                if k in feparts:
                    vlist = feparts[k]
                    if not isinstance(vlist, list):
                        vlist = [vlist]
                    vlist.append(v)
                    feparts[k] = vlist
                else:
                    feparts[k] = v
            else:
                # logger.debug(f"skipped part {p} with kv {kv}")
                pass
        
        return feparts
    
    @classmethod
    def format_chrom(cls, chrom: str) -> Tuple[str, str]:
        """Format chromosome for tabix query and internal use."""
        chr_raw = str(chrom).removeprefix("chr")
        chr_formatted = cls.chr_format.format(chrstr=chr_raw)
        return chr_formatted, chr_raw
    
    @classmethod
    def get_headers(cls):
        hdrs = []
        for c in cls.columns + cls.derived:
            if isinstance(c, str):
                c = cls.core_columns.get(c, "")
            if not c:
                continue
            hdrs.append(c.name)
        # hdrs= [c.name for c in cls.columns] + [c.name for c in cls.derived]
        return hdrs


class MiniStream(TabularStream):
    
    delimiter = "\t"
    columns = [
        "chrom","start","end",
        ColumnSpec("id", str, "", target="core"),
        ColumnSpec("type", str),
        "name", "strand", 
        ColumnSpec("mean_exp", float),
        ColumnSpec("total_exp", float),
        ColumnSpec("num_experiments", int),
        ColumnSpec("num_expressed", int),
    ]
    
    chr_format = "chr{chrstr}"
    
    def __init__(self, file_path, source_name, feature_type, filter = None):
        super().__init__(source_name)
        self.file_path =Path(file_path)
        self.feature_type = feature_type
        self._data = {}
        self.load_data()
        self._filter = filter
            
    def stream(self, chrom: Optional[str] = None, 
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UFeature]:
        
        if not chrom:
            return []
        
        clsname = type(self).__name__
        
        chrom = chrom.removeprefix("chr")
        
        query_start = int(start) if start else 1
        query_end = int(end) if end else chr_lens.get(chrom, 999999999)
        
        feats = self._data.get(chrom, [])
        
        for feature in feats:               
            if self.validate_feature(feature):
                # Apply region filter
                # logger.debug(f"feature passed validation: {feature.name}")
                
                if clsname == "VariantStream":
                    pass
                
                if chrom and feature.chrom != chrom:
                    continue
                if query_start and feature.end < query_start:
                    continue
                if query_end and feature.start > query_end:
                    continue
                logger.debug(f"feature passed all challenges: {feature.name}")
                yield feature
    
    def validate_feature(self, feature):        
        if self._filter:
            return self._filter(feature)
        else:
            return bool(feature)
    
    def load_data(self):
        
        if '.gz' in self.file_path.suffixes:
            with gzip.open(self.file_path) as f:
                data = f.readlines()
            data = [d.decode("utf-8") for d in data]
        else:
            with open(self.file_path, "r") as f:
                data = f.readlines()
        
        feats = {}
        for i,line in enumerate(data):
            if i==0: continue
            
            f = self.parse_line(line)
            if not f.chrom in feats:
                feats[f.chrom] = []
            feats[f.chrom].append(f)
        ddict = {}
        for chrom, fs in feats.items():
            fs_srt = list(sorted(fs, key = lambda f:f.start))
            ddict[chrom] = fs_srt
        
        logger.setLevel("DEBUG")
        logger.debug(fs_srt[0])
        self._data = ddict

def get_gs_feature_type(f):
    
    fts = []
    
    if "pseudo" in f.gene_biotype:
        return "pseudogene"
    elif "lncRNA" in f.gene_biotype:
        if f.feature_type == 'gene':
            ft = 'lncRNA'
        else:
            ft = "lncRNA" +'_'+ f.feature_type
    elif "RNA" in f.gene_biotype:
        return "ncRNA"
    else:
        ft = f.feature_type
    
    fts.append(ft)
    
    if "nonsense" in f.gene_biotype:
        fts.append("NMD")
    
    return "_".join(fts)
    
def get_gs_feature_name(f):
    
    if f.transcript_name:
        parname = f.transcript_name
        if f.feature_type == "exon":
            return parname + f'-ex{f.exon_number}'
        elif f.feature_type == "CDS":
            return parname + f"-CDS{f.exon_number}"
        elif f.feature_type == "three_prime_utr":
            return parname + '-3pUTR'
        elif f.feature_type == "five_prime_utr":
            return parname + '-5pUTR'
        elif f.feature_type == "stop_codon":
            return parname + '-stop'
        elif f.feature_type == "start_codon":
            return parname + '-start'
        else:
            return parname
    elif f.gene_name:
        return f.gene_name

def get_gs_gene_id(f):
    
    fid = ""
    
    if f.gene_id:
        fid = f.gene_id
    elif f.transcript_id:
        fid = f.transcript_id
    else:
        fid = f"{f.feature_type}:chr{f.chrom}:{f.start}-{f.end}"
    
    return fid

class GeneStream(TabularStream):
    
    feature_type = ""
    chr_format = "{chrstr}"
    fe_delimiter = ";"
    fe_kv_delimiter = " "
    
    columns = [
        "chrom","source",
        # "feature_type",
        ColumnSpec("feature_type", str, ""),
        "start","end", None, "strand", ColumnSpec("frame", int, default = 0),
        ColumnSpec("gene_id", str, ""),
        ColumnSpec("gene_name", str, ""),
        ColumnSpec("gene_source", str, ""),
        ColumnSpec("gene_biotype", str, ""),
        ColumnSpec("gene_version", str, ""),
        
        ColumnSpec("transcript_id", str, ""),
        ColumnSpec("ccds_id", str, ""),
        ColumnSpec("transcript_name", str, ""),
        ColumnSpec("transcript_source", str, ""),
        ColumnSpec("transcript_biotype", str, ""),
        ColumnSpec("transcript_version", str, ""),
        
        ColumnSpec("exon_id", str, ""),
        ColumnSpec("exon_number", int, ""),
        ColumnSpec("exon_version", str, ""),
        
        ColumnSpec("protein_id", str, ""),
        ColumnSpec("protein_version", str, ""),
        
        ColumnSpec("tag", str, ""),
    ]
    # other cols
    # ColumnSpec("transcript_support_level", str, ""),
    
    derived = [
        DerivedSpec("feature_type", get_gs_feature_type, default = "gene"),
        DerivedSpec("id", get_gs_gene_id, default = "[no_id]"),
        DerivedSpec("name", get_gs_feature_name, default = "?"),
    ]
    
    def __init__(self, filepath: str):
        super().__init__(filepath, source_name="GTF")
    
    def preparse_line(self, line:str):
        
        parts = line.split(self.delimiter)
        
        fpart = parts[-1]
        parts = parts[:-1]
        
        fe_dict = self.format_group_entry(fpart)
        
        return parts, fe_dict

def get_genotype(f):
    
    gts = {"0":f.ref, "1":f.alt[0], "2":f.alt[1]}
    
    gt = f.GT
    
    for k, gtb in gts.items():
        gt = gt.replace(k, gtb)
    
    return gt
    

class VariantStream(TabularStream):
    """
    chr1	1078047	.	G	A	5.67	PASS	AC=1;AF=0.5;AN=2;DP=6;FS=0;MQ=20.88;MQRankSum=0.72;QD=4.11;ReadPosRankSum=-1.38;SOR=2.303;FractionInformativeReads=1	GT:AD:AF:DP:F1R2:F2R1:GQ:PL:GP:PRI:SB:MB	0/1:4,2:0.333:6:4,1:0,1:4:38,0,2:5.6683,2.4607,7.9203:0,34.77,37.78:2,2,0,2:2,2,1,1
    """
    feature_type: str = "variant"
    chr_format = "chr{chrstr}"
    
    fe_delimiter = ";"
    fe_kv_delimiter = "="
    
    columns = [
        "chrom", "start", None,
        ColumnSpec("ref", str, ""),
        ColumnSpec("alt", str, "", formatter = lambda s: s.split(',') if ',' in s else (s,s)),
        ColumnSpec("qual", float, ""),
        ColumnSpec("depth_result", str, ""),
        
        ColumnSpec("AD", str, ""), # Allelic depths        
        ColumnSpec("AC", str, ""), # Allele count
        ColumnSpec("AF", str, ""), # Allele frequency
        ColumnSpec("AN", str, ""), # Allele number
        ColumnSpec("GQ", str, ""), # genotype quality
        ColumnSpec("GT", str, ""), # Genotype
        ColumnSpec("PS", str, ""), # Physical phasing ID information
        ColumnSpec("SQ", str, ""), # somatic quality
        ColumnSpec("DP", str, ""), # read depth, approximate, unfiltered
        
        # ColumnSpec("F1R2", str, ""), # Count of reads in F1R2 pair orientation supporting each allele
        # ColumnSpec("F2R1", str, ""), # same ^ 
        # ColumnSpec("GP", str, ""), # Phred-scaled posterior p for genotypes
        # ColumnSpec("MB", str, ""), # Per-sample component statistics
        # ColumnSpec("PL", str, ""), # Normalized, phred-scaled likelihoods for genotypes
        # ColumnSpec("PRI", str, ""), # Phred-scaled prior probabilities for genotypes
        # ColumnSpec("SB", str, ""), # per-sample component statistics
        # ColumnSpec("FS", float, ""),  # phred-scaled p-value, Fisher's exact test (strand bias)
        # ColumnSpec("MQ", float, ""), # RMS mapping quality
        # ColumnSpec("MQRankSum", float, ""), # Z-score for wilcoxon rank sum test of alt vs ref read mapping qualities
        # ColumnSpec("ReadPosRankSum", float, ""),  # same but for read position bias
        # ColumnSpec("SOR", float, ""), # symmetric odds ratio of 2x2 contigency table (strand bias)
        # ColumnSpec("FractionInformativeReads", float, ""),
    ]
    
    derived = [
        DerivedSpec("end", lambda uf:uf.start + 1),
        DerivedSpec("delta", lambda uf: max(len(a) for a in uf.alt)-len(uf.ref)),
        DerivedSpec("genotype", get_genotype),
    ]
    
    def __init__(self, filepath):
        super().__init__(filepath, source_name = "Variants")
    
    def parse_line(self, line:str):
        f = super().parse_line(line)
        f.attributes['_raw_line'] = line
        return f
    
    def preparse_line(self, line:str):
        parts = line.split(self.delimiter)
        
        fe0, fe1ks, fe1vs = parts[-3:]
        parts = parts[:-3]
        
        fedict = self.format_group_entry(fe0)
        
        for k, v in zip(fe1ks.split(":"), fe1vs.split(":")):
            fedict[k.strip()] = v.strip()
        
        return parts, fedict

class RepeatMaskerStream(TabularStream):
    
    feature_type = "repeat"
    chr_format = "chr{chrstr}"
    
    columns = [
        "chrom", "start", "end",
        ColumnSpec("name", str, None),
        ColumnSpec("repeat_type", str, None),
        ColumnSpec("strand", str, formatter = lambda s:"-" if s == "C" else s)
    ]
    
    def __init__(self, filepath: str):
        super().__init__(filepath, source_name="RepeatMasker")

class SimpleRepeatStream(RepeatMaskerStream):
    
    feature_type = "simple_repeat"
    
    columns = RepeatMaskerStream.columns + [
        ColumnSpec("motif", str, "",
            formatter=lambda s: s.strip("()n") if s.endswith(")n") else s),
    ]
    
    derived = [
        DerivedSpec("motif", lambda uf: uf.name.strip("()n")),
        DerivedSpec("repeat_count",
            lambda uf: (uf.end - uf.start) / len(uf.motif),
            default=None
        ),
    ]
    
    def __init__(self, filepath: str):
        super().__init__(filepath)
        
    def validate_feature(self, feature):
        res = bool(feature) and feature.repeat_type == "Simple_repeat"
        return res
        
class DfamStream(TabularStream):
    
    feature_type = "dfam_hit"
    chr_format = "chr{chrstr}"
    
    columns = [
        "chrom", "start", "end",
        ColumnSpec("name", str, None),
        ColumnSpec("bits", float, None),
        "strand",
        ColumnSpec("id", str, None),
        ColumnSpec("e-value", float, None),
        ColumnSpec("bias", float, None),
        ColumnSpec("start_hmm", int, None),
        ColumnSpec("end_hmm", int, None),
        ColumnSpec("start_ali", int, None),
        ColumnSpec("end_ali", int, None),
        ColumnSpec("sq-len", float, None),
        ColumnSpec("kimura_div", float, None),
    ]
    
    def __init__(self, filepath: str):
        super().__init__(filepath, source_name="Dfam")
        

class ClinvarStream(TabularStream):
    """Stream ClinVar clinical variant annotations.

    Uses the new TabularStream/ColumnSpec pattern for clean column definitions.
    """
    feature_type = "clinical_variant"
    chr_format = "{chrstr}"

    # Define columns with types and targets
    columns = [
        "chrom", "start", "end", "id",
        ColumnSpec("ref", str, ""),
        ColumnSpec("alt", str, ""),
        ColumnSpec("clnsig", str, ""),
        ColumnSpec("clndn", str, ""),  # Clinical disease name
        ColumnSpec("loc", str, ""),
        ColumnSpec("clnvc", str, ""),  # Variant class
        ColumnSpec("af_esp", float, None, formatter = lambda s:s.lstrip("C=")),
        ColumnSpec("af_exac", float, None, formatter = lambda s:s.lstrip("C=")),
        ColumnSpec("af_tgp", float, None, formatter = lambda s:s.lstrip("C=")),
        ColumnSpec("clndnincl", str, ""),   # empty ?
        ColumnSpec("mc", str, ""),  # Molecular consequence
    ]

    # Example derived attributes
    derived = [
        DerivedSpec('gene_name',
            lambda uf: uf.loc.split(":")[0]
        ),
        DerivedSpec("is_pathogenic",
            lambda uf: "pathogenic" in uf.get("clnsig", "").lower()
        ),
    ]

    def __init__(self, filepath: str):
        super().__init__(filepath, source_name="ClinVar")

def get_experiment_stream(file_path, source_name, feature_type, filter = None, gm = None):
    
    mstream = MiniStream(file_path, source_name, feature_type, filter = filter)
    if gm:
        gm.annotations.add_source(source_name, mstream)
    return mstream

class UGenomeAnnotations:
    """Unified interface for all genomic annotations and sequences.
    
    Merges multiple annotation streams and sequence data efficiently.
    """
    
    def __init__(self, fasta_path: Optional[str] = None, vcf_path: Optional[str] = None):
        self.streams = {}  # Annotation streams (tabix-indexed)
        self.motif_streams = {}  # Motif streams (sequence-scanning)
        self.indices = {}  # For fast range queries
        self.sequence_stream = None

        # Initialize sequence streaming if FASTA provided
        if fasta_path:
            self.setup_sequence_stream(fasta_path, vcf_path)
        
    def setup_sequence_stream(self, fasta_path: str, vcf_path: Optional[str] = None,
                            min_qual: float = 5.0):
        """Setup sequence streaming with optional variant integration.
        
        Args:
            fasta_path: Path to FASTA file
            vcf_path: Optional path to VCF file for variants
            min_qual: Minimum quality threshold for variants
        """
        from .sequences import FASTAStream, SequenceStreamWithVariants
        
        try:
            fasta_stream = FASTAStream(fasta_path)
            if vcf_path:
                self.sequence_stream = SequenceStreamWithVariants(
                    fasta_stream, vcf_path, min_qual
                )
                logger.info(f"Initialized sequence stream with variants from {vcf_path}")
            else:
                self.sequence_stream = fasta_stream
                logger.info(f"Initialized sequence stream from {fasta_path}")
                
        except Exception as e:
            logger.error(f"Failed to setup sequence stream: {e}")
    
    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        """Get reference sequence for a region."""
        # logger.debug("get_sequence called on UGenomeAnnotations")
        if self.sequence_stream:
            return self.sequence_stream.get_sequence(chrom, start, end)
        return ""
    
    def get_personal_sequence(self, chrom: str, start: int, end: int) -> str:
        """Get personalized sequence with variants."""
        if self.sequence_stream:
            return self.sequence_stream.get_personal_sequence(chrom, start, end)
        return ""
    
    def get_aligned_sequences(self, chrom: str, start: int, end: int) -> Tuple[str, str, List[Tuple[int, int]]]:
        """Get aligned sequences with gaps for visualization.
        
        Returns:
            Tuple of (aligned_ref, aligned_pers, variant_list)
            where variant_list contains tuples of (position, delta)
        """
        if self.sequence_stream and hasattr(self.sequence_stream, 'get_sequence_with_alignment'):
            return self.sequence_stream.get_sequence_with_alignment(chrom, start, end)
        
        # Fallback to unaligned with empty variant list
        ref = self.get_sequence(chrom, start, end)
        pers = self.get_personal_sequence(chrom, start, end)
        return ref, pers, []
    
    def add_source(self, name: str, stream: AnnotationStream):
        """Add an annotation source."""
        self.streams[name] = stream
        logger.info(f"Added annotation source: {name}")
    
    def add_genes(self, filepath:str, name:str = "genes"):
        self.add_source(name, GeneStream(filepath))
        
    def add_variants(self, filepath:str, name:str = "variants"):
        self.add_source(name, VariantStream(filepath))
    
    def add_repeatmasker(self, filepath: str, name: str = "repeats"):
        """Add RepeatMasker annotation source."""
        self.add_source(name, RepeatMaskerStream(filepath))
    
    def add_dfam(self, filepath, name="dfam"):
        self.add_source(name, DfamStream(filepath))
    
    def add_clinvar(self, filepath, name = "clinvar"):
        self.add_source(name, ClinvarStream(filepath))
    
    def add_motifs(self, motif_stream, name: str = "motifs"):
        """Register a motif stream for sequence scanning.

        The motif stream is bound to this object's sequence getter so it can
        fetch sequences during scanning. Motif streams are included in stream_all().

        Args:
            motif_stream: A MotifStream subclass instance (JasparStream, PatternStream, etc.)
            name: Name for this motif source
        """
        from ggene.database.motifs import MotifStream

        if not isinstance(motif_stream, MotifStream):
            raise TypeError(f"Expected MotifStream, got {type(motif_stream)}")

        if not self.sequence_stream:
            logger.warning(f"No sequence stream configured. Motif stream '{name}' "
                          "will fail until setup_sequence_stream() is called.")

        # Bind our sequence getter to the motif stream
        motif_stream.bind_sequence_getter(self.get_sequence)
        self.motif_streams[name] = motif_stream
        logger.info(f"Added motif source: {name}")
    
    def stream_all(self, chrom: Optional[str] = None,
                   start: Optional[int] = None,
                   end: Optional[int] = None) -> Iterator[UFeature]:
        """Stream all annotations and motifs merged by position.

        Uses heapq.merge for efficient sorted merging. Annotation streams
        query indexed files while motif streams scan sequences in chunks.
        """
        iterators = []

        # Add annotation stream iterators
        for name, stream in self.streams.items():
            try:
                iterator = stream.stream(chrom, start, end)
                iterators.append(iterator)
            except Exception as e:
                logger.warning(f"Failed to create annotation stream for {name}: {e}")

        # Add motif stream iterators (they'll scan sequences lazily)
        # for name, stream in self.motif_streams.items():
        #     try:
        #         iterator = stream.stream(chrom, start, end)
        #         iterators.append(iterator)
        #     except Exception as e:
        #         logger.warning(f"Failed to create motif stream for {name}: {e}")

        # Merge all iterators by position
        for feature in heapq.merge(*iterators):
            yield feature

    
    def query_point(self, chrom: str, position: int) -> List[UFeature]:
        """Query all features at a specific position."""
        features = []
        for name, stream in self.streams.items():
            try:
                for feature in stream.query_range(chrom, position, position):
                    features.append(feature)
            except Exception as e:
                logger.warning(f"Query failed for {name}: {e}")
        return sorted(features)
    
    def query_range(self, chrom: str, start: int, end: int) -> List[UFeature]:
        """Query all features in a range."""
        features = []
        for name, stream in self.streams.items():
            try:
                features.extend(stream.query_range(chrom, start, end))
            except Exception as e:
                logger.warning(f"Query failed for {name}: {e}")
        return sorted(features)
        
    def stream_all_motifs(self, chrom: Optional[str] = None,
                   start: Optional[int] = None,
                   end: Optional[int] = None) -> Iterator[UFeature]:
        iterators = []
        # Add motif stream iterators (they'll scan sequences lazily)
        for name, stream in self.motif_streams.items():
            try:
                iterator = stream.stream(chrom, start, end)
                iterators.append(iterator)
            except Exception as e:
                logger.warning(f"Failed to create motif stream for {name}: {e}")
                
        for feature in heapq.merge(*iterators):
            yield feature
            
    def query_motifs(self, chrom:str, start: int, end:int):
        
        seq = self.get_sequence(chrom, start, end)
        
        motifs = []
        
        for name, mstream in self.motif_streams.items():
            mtfs = mstream.scan_sequence(seq, chrom, start)
            motifs.extend(mtfs)
        
        return motifs
    
    def stream_by_types(self, feature_types: List[str],
                       chrom: Optional[str] = None,
                       start: Optional[int] = None,
                       end: Optional[int] = None) -> Iterator[UFeature]:
        """Stream only specific feature types."""
        
        # maybe should make this one not stream_all dependent (for performance)
        for feature in self.stream_all(chrom, start, end):
            if feature.feature_type in feature_types:
                yield feature
    
    def get_summary(self, chrom: str, start: int, end: int) -> Dict[str, int]:
        """Get summary statistics for a region."""
        summary = {}
        for feature in self.query_range(chrom, start, end):
            key = f"{feature.source}:{feature.feature_type}"
            summary[key] = summary.get(key, 0) + 1
        return summary
    
    def to_bedgraph(self, chrom: str, start: int, end: int,
                    feature_type: str, window: int = 100) -> List[Tuple[int, int, float]]:
        """Convert to bedGraph format for visualization."""
        # Count features in windows
        scores = {}
        for pos in range(start, end, window):
            window_end = min(pos + window, end)
            count = 0
            for feature in self.query_range(chrom, pos, window_end):
                if feature.feature_type == feature_type:
                    count += 1
            scores[pos] = count / window
            
        # Convert to bedGraph format
        bedgraph = []
        for pos in sorted(scores.keys()):
            bedgraph.append((pos, min(pos + window, end), scores[pos]))
        
        return bedgraph


class CachedUnifiedAnnotations(UGenomeAnnotations):
    """Cached version for better performance."""
    
    def __init__(self, cache_size: int = 10000):
        super().__init__()
        self.cache = {}
        self.cache_size = cache_size
        
    def query_range(self, chrom: str, start: int, end: int) -> List[UFeature]:
        """Query with caching."""
        cache_key = (chrom, start, end)
        
        if cache_key in self.cache:
            return self.cache[cache_key]
            
        # Query and cache
        features = super().query_range(chrom, start, end)
        
        # Simple LRU-ish cache
        if len(self.cache) >= self.cache_size:
            # Remove oldest entry
            oldest = next(iter(self.cache))
            del self.cache[oldest]
            
        self.cache[cache_key] = features
        return features


# Example usage
def example_usage():
    """Demonstrate the unified annotation system."""
    
    # Create unified annotation manager
    annotations = CachedUnifiedAnnotations()
    
    # Add various sources
    annotations.add_gtf("/path/to/genes.gtf.gz", "genes")
    annotations.add_vcf("/path/to/variants.vcf.gz", "variants")
    annotations.add_bed("/path/to/peaks.bed", "chip_peaks", "peak")
    annotations.add_repeatmasker("/path/to/repeats.out", "repeats")
    
    # Stream all annotations for a region
    print("All features in region:")
    for feature in annotations.stream_all("chr1", 1000000, 1001000):
        print(f"  {feature.source}:{feature.feature_type} at {feature.start}-{feature.end}")
    
    # Query specific position
    print("\nFeatures at position chr1:1000500:")
    for feature in annotations.query_point("chr1", 1000500):
        print(f"  {feature.name} ({feature.feature_type})")
    
    # Get summary
    summary = annotations.get_summary("chr1", 1000000, 2000000)
    print(f"\nSummary: {summary}")
    
    # Stream only genes and variants
    print("\nGenes and variants only:")
    for feature in annotations.stream_by_types(["gene", "variant"], "chr1", 1000000, 1100000):
        print(f"  {feature.feature_type}: {feature.name}")


if __name__ == "__main__":
    example_usage()