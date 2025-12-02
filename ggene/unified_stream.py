"""Unified streaming interface for genomic annotations.

This module provides a common interface to stream and merge annotations from
multiple genomic databases (GTF, VCF, BED, JASPAR, etc.) into a unified
Feature-based representation.

IMPORTANT: All data files should be bgzipped and tabix-indexed for performance.
Without indexing, queries will be extremely slow on large files.
"""

import heapq
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Iterator, Dict, Any, List, Optional, Tuple, Union
import gzip
import json
from pathlib import Path
import logging
import pysam
import os
import warnings

logger = logging.getLogger(__name__)
logger.setLevel("CRITICAL")
# logger.setLevel(logging.DEBUG)


@dataclass
class UnifiedFeature:
    """Common representation for all genomic features.
    
    This is the lingua franca for all annotation sources.
    """
    chrom: str
    start: int  # 1-based, inclusive
    end: int    # 1-based, inclusive
    feature_type: str  # gene, variant, motif, repeat, etc.
    source: str  # GTF, VCF, JASPAR, etc.
    score: Optional[float] = None
    strand: Optional[str] = ""
    name: Optional[str] = ""
    id: Optional[str] = ""
    attributes: Dict[str, Any] = field(default_factory=dict)
    
    def __lt__(self, other):
        """For heap-based merging."""
        return (self.chrom, self.start, self.end) < (other.chrom, other.start, other.end)
    
    def overlaps(self, start: int, end: int) -> bool:
        """Check if feature overlaps with given range."""
        return not (self.end < start or self.start > end)
    
    def get(self, att, default=None):
        if hasattr(self, att):
            return getattr(self, att)
        else:
            return self.attributes.get(att,default)
    
    def __getitem__(self, att):
        if isinstance(att, str):
            return self.get(att)
        else:
            raise IndexError
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'type': self.feature_type,
            'source': self.source,
            'score': self.score,
            'strand': self.strand,
            'name': self.name,
            'id': self.id,
            'attributes': self.attributes
        }


class AnnotationStream(ABC):
    """Abstract base class for annotation streams."""
    
    def __init__(self, source_name: str):
        self.source_name = source_name
        self._current_chrom = None
        self._buffer = []
        
    @abstractmethod
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Parse a single line into a UnifiedFeature."""
        pass
    
    @abstractmethod
    def stream(self, chrom: Optional[str] = None, 
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream features, optionally filtered by region."""
        pass
    
    def query_range(self, chrom: str, start: int, end: int) -> List[UnifiedFeature]:
        """Query features in a specific range."""
        features = []
        for feature in self.stream(chrom, start, end):
            if feature.overlaps(start, end):
                features.append(feature)
        return features


class GTFStream(AnnotationStream):
    """Stream GTF/GFF annotations using indexed access."""
    
    max_indices = {'1': 248937043, '10': 133778498, '11': 135075908, '12': 133238549, '13': 114346637, '14': 106879812, '15': 101979093, '16': 90222678, '17': 83240391, '18': 80247514, '19': 58599303, '2': 242175634, '20': 64327972, '21': 46691226, '22': 50799123, '3': 198228376, '4': 190195978, '5': 181472430, '6': 170745977, '7': 159233377, '8': 145066516, '9': 138320835, 'MT': 16023, 'X': 156027877, 'Y': 57214397}
    
    def __init__(self, filepath: str):
        super().__init__("GTF")
        self.filepath = Path(filepath)
        
        # Check if file is indexed
        self.tabix = None
        if self.filepath.suffix == '.gz':
            # Check for index file
            index_file = Path(str(self.filepath) + '.tbi')
            if index_file.exists():
                try:
                    self.tabix = pysam.TabixFile(str(self.filepath))
                    # self.tabix.
                    logger.info(f"Using indexed access for {self.filepath}")
                except Exception as e:
                    logger.warning(f"Failed to open tabix index: {e}")
        
        if not self.tabix:
            logger.warning(f"No index found for {self.filepath}. Performance will be slow!")
            logger.info(f"Create index with: bgzip -c {self.filepath} > {self.filepath}.gz && tabix -p gff {self.filepath}.gz")
        
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Parse GTF line."""
        if line.startswith('#'):
            return None
            
        parts = line.strip().split('\t')
        if len(parts) < 9:
            return None
            
        # Parse attributes
        attr_dict = {}
        for attr in parts[8].split(';'):
            if attr.strip():
                if ' ' in attr:
                    key, val = attr.strip().split(' ', 1)
                    attr_dict[key] = val.strip('"')
        
        return UnifiedFeature(
            chrom=parts[0],
            start=int(parts[3]),
            end=int(parts[4]),
            feature_type=parts[2],
            source=self.source_name,
            score=float(parts[5]) if parts[5] != '.' else None,
            strand=parts[6] if parts[6] != '.' else None,
            name=attr_dict.get('gene_name'),
            id=attr_dict.get('gene_id'),
            attributes=attr_dict
        )
    
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream GTF features using indexed access when available."""
        
        # Use indexed access if available
        if self.tabix and chrom:
            try:
                # Tabix uses 0-based half-open intervals
                query_start = (start - 1) if start else 0
                query_end = end if end else 999999999
                
                for line in self.tabix.fetch(chrom, query_start, query_end):
                    feature = self._parse_line(line)
                    if feature:
                        yield feature
                return
            except Exception as e:
                logger.debug(f"Tabix fetch failed, falling back to linear scan: {e}")
        
        # Fallback to linear scan (slow!)
        open_func = gzip.open if self.filepath.suffix == '.gz' else open
        with open_func(self.filepath, 'rt') as f:
            for line in f:
                feature = self._parse_line(line)
                if feature:
                    # Filter by region if specified
                    if chrom and feature.chrom != chrom:
                        continue
                    if start and feature.end < start:
                        continue
                    if end and feature.start > end:
                        continue
                    yield feature


class VCFStream(AnnotationStream):
    """Stream VCF variant annotations using cyvcf2 for efficient indexed access."""
    
    def __init__(self, filepath: str):
        super().__init__("VCF")
        self.filepath = Path(filepath)
        
        # Use cyvcf2.VCF which handles indexing internally
        try:
            from cyvcf2 import VCF
            self.vcf = VCF(str(self.filepath))
            
            # Check if index exists and set it
            index_files = [
                Path(str(self.filepath) + '.tbi'),
                Path(str(self.filepath) + '.csi')
            ]
            for index_file in index_files:
                if index_file.exists():
                    try:
                        self.vcf.set_index(str(index_file))
                        logger.info(f"Using indexed VCF access for {self.filepath}")
                    except:
                        pass  # VCF might auto-detect the index
                    break
            else:
                logger.warning(f"No index found for {self.filepath}. Create with: tabix -p vcf {self.filepath}")
                
        except ImportError:
            logger.error("cyvcf2 not installed. Install with: pip install cyvcf2")
            self.vcf = None
        except Exception as e:
            logger.warning(f"Failed to open VCF with cyvcf2: {e}")
            self.vcf = None
        
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Parse VCF line."""
        if line.startswith('#'):
            return None
            
        parts = line.strip().split('\t')
        if len(parts) < 8:
            return None
        
        # Parse INFO field
        info_dict = {}
        if parts[7] != '.':
            for item in parts[7].split(';'):
                if '=' in item:
                    key, val = item.split('=', 1)
                    info_dict[key] = val
                else:
                    info_dict[item] = True
        
        # Determine variant type
        ref = parts[3]
        alts = parts[4].split(',')
        
        variant_types = []
        for alt in alts:
            if len(ref) == len(alt) == 1:
                variant_types.append('SNP')
            elif len(ref) > len(alt):
                variant_types.append('DEL')
            elif len(ref) < len(alt):
                variant_types.append('INS')
            else:
                variant_types.append('COMPLEX')
        
        return UnifiedFeature(
            chrom=parts[0],
            start=int(parts[1]),
            end=int(parts[1]) + len(ref) - 1,
            feature_type='variant',
            source=self.source_name,
            score=float(parts[5]) if parts[5] != '.' else None,
            name=parts[2] if parts[2] != '.' else None,
            id=parts[2] if parts[2] != '.' else None,
            attributes={
                'ref': ref,
                'alt': alts,
                'qual': parts[5],
                'filter': parts[6],
                'info': info_dict,
                'variant_types': variant_types,
                'genotype': ""
            }
        )
    
    def _parse_variant(self, variant) -> Optional[UnifiedFeature]:
        """Parse cyvcf2 Variant object to UnifiedFeature."""
        # Parse INFO field
        info_dict = dict(variant.INFO)
        
        # Determine variant type
        ref = variant.REF
        alts = variant.ALT if variant.ALT else []
        
        variant_types = []
        for alt in alts:
            if alt is None:
                continue
            if len(ref) == len(alt) == 1:
                variant_types.append('SNP')
            elif len(ref) > len(alt):
                variant_types.append('DEL')
            elif len(ref) < len(alt):
                variant_types.append('INS')
            else:
                variant_types.append('COMPLEX')
        
        if variant.num_het > 0:
            zyg = "heterozygous"
        elif variant.num_hom_alt > 0:
            zyg = "homozygous alt"
        elif variant.num_hom_ref > 0:
            zyg = "homozygous ref"
        else:
            zyg = ""
        
        return UnifiedFeature(
            chrom=variant.CHROM,
            start=variant.POS,  # VCF is 1-based
            end=variant.POS + len(ref) - 1,
            feature_type='variant',
            source=self.source_name,
            score=variant.QUAL if variant.QUAL else None,
            name=variant.ID if variant.ID else None,
            id=variant.ID if variant.ID else None,
            attributes={
                'ref': ref,
                'alt': alts,
                'qual': variant.QUAL,
                'filter': variant.FILTER,
                'info': info_dict,
                'variant_types': variant_types,
                'genotype': variant.gt_bases[0],
                'zygosity':zyg
            }
        )
    
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream VCF features using cyvcf2's efficient indexed access."""
        
        # Use cyvcf2's efficient querying
        if self.vcf and chrom:
            chrstr = str(chrom)
            
            # Try both with and without chr prefix
            for chrom_format in [f"chr{chrstr}" if not chrstr.startswith('chr') else chrstr,
                                chrstr.lstrip('chr') if chrstr.startswith('chr') else chrstr]:
                try:
                    # Create region string for cyvcf2
                    if start and end:
                        region = f"{chrom_format}:{start}-{end}"
                    elif start:
                        region = f"{chrom_format}:{start}-"
                    else:
                        region = chrom_format
                    
                    
                    with warnings.catch_warnings(action="ignore"):
                        # Query using cyvcf2's efficient indexed access
                        found_any = False
                        for variant in self.vcf(region):
                            feature = self._parse_variant(variant)
                            if feature:
                                found_any = True
                                yield feature
                    
                    if found_any:
                        return  # Successfully found variants, stop trying
                except Exception as e:
                    logger.debug(f"VCF query failed for {region}: {e}")
                    continue  # Try next format
        
        # Fallback to using pysam if cyvcf2 isn't available
        # This is much slower but works as a backup
        if not self.vcf:
            logger.warning("Falling back to slow VCF parsing without cyvcf2")
            open_func = gzip.open if self.filepath.suffix == '.gz' else open
            
            with open_func(self.filepath, 'rt') as f:
                for line in f:
                    feature = self._parse_line(line)
                    if feature:
                        if chrom and feature.chrom != chrom:
                            continue
                        if start and feature.end < start:
                            continue
                        if end and feature.start > end:
                            continue
                        yield feature


class BEDStream(AnnotationStream):
    """Stream BED format annotations using indexed access."""
    
    def __init__(self, filepath: str, feature_type: str = "region"):
        super().__init__("BED")
        self.filepath = Path(filepath)
        self.feature_type = feature_type
        
        # Check if file is indexed
        self.tabix = None
        if self.filepath.suffix == '.gz':
            index_file = Path(str(self.filepath) + '.tbi')
            if index_file.exists():
                try:
                    print("doing tabix")
                    self.tabix = pysam.TabixFile(str(self.filepath))
                    logger.info(f"Using indexed access for {self.filepath}")
                except Exception as e:
                    logger.warning(f"Failed to open tabix index: {e}")
        
        if not self.tabix and self.filepath.suffix == '.gz':
            logger.warning(f"No index found for {self.filepath}. Performance will be slow!")
            logger.info(f"Create index with: tabix -p bed {self.filepath}")
        
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Parse BED line."""
        if line.startswith('#') or line.startswith('track'):
            return None
            
        parts = line.strip().split('\t')
        if len(parts) < 3:
            return None
        
        return UnifiedFeature(
            chrom=parts[0].removeprefix("chr"),
            start=int(parts[1]) + 1,  # BED is 0-based
            end=int(parts[2]),
            feature_type=self.feature_type,
            source=self.source_name,
            name=parts[3] if len(parts) > 3 else None,
            score = 0,
            strand=parts[5] if len(parts) > 5 else None,
            attributes={
                'motif': parts[4],
                'type':parts[5]
                # 'thickStart': int(parts[6]) if len(parts) > 6 else None,
                # 'thickEnd': int(parts[7]) if len(parts) > 7 else None,
                # 'itemRgb': parts[8] if len(parts) > 8 else None
            }
        )
    
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream BED features."""
        with open(self.filepath, 'r') as f:
            for line in f:
                feature = self._parse_line(line)
                if feature:
                    if chrom and feature.chrom != chrom:
                        continue
                    if start and feature.end < start:
                        continue
                    if end and feature.start > end:
                        continue
                    yield feature


class JASPARStream(AnnotationStream):
    """Stream JASPAR motif predictions."""
    
    def __init__(self, jaspar_file: str, genome_sequence_getter):
        super().__init__("JASPAR")
        self.jaspar_file = jaspar_file
        self.get_sequence = genome_sequence_getter
        self.motifs = self._load_jaspar_motifs()
        
    def _load_jaspar_motifs(self) -> Dict:
        """Load JASPAR PWMs."""
        # This would parse JASPAR format files
        # For now, returning empty dict
        return {}
    
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Not used for JASPAR - we scan sequences."""
        pass
    
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream motif predictions by scanning sequence."""
        if not chrom or not start or not end:
            return
            
        # Get sequence for region
        seq = self.get_sequence(chrom, start, end)
        
        # Scan with each motif
        for motif_id, pwm in self.motifs.items():
            hits = pwm.scan_sequence(seq)
            for hit_start, score, matched_seq in hits:
                yield UnifiedFeature(
                    chrom=chrom,
                    start=start + hit_start,
                    end=start + hit_start + len(matched_seq) - 1,
                    feature_type='tf_binding',
                    source=self.source_name,
                    score=score,
                    name=motif_id,
                    id=f"{motif_id}_{chrom}_{start + hit_start}",
                    attributes={
                        'sequence': matched_seq,
                        'motif_id': motif_id
                    }
                )


class RepeatMaskerStream(AnnotationStream):
    """Stream RepeatMasker annotations."""
    
    def __init__(self, filepath: str):
        super().__init__("RepeatMasker")
        self.filepath = Path(filepath)
        
    def _parse_line(self, line: str) -> Optional[UnifiedFeature]:
        """Parse RepeatMasker .out format."""
        parts = line.strip().split()
        if len(parts) < 15 or parts[0].isdigit() == False:
            return None
            
        return UnifiedFeature(
            chrom=parts[4],
            start=int(parts[5]),
            end=int(parts[6]),
            feature_type='repeat',
            source=self.source_name,
            name=parts[9],  # Repeat name
            strand=parts[8],
            score=float(parts[0]),  # SW score
            attributes={
                'repeat_class': parts[10],
                'repeat_family': parts[11] if len(parts) > 11 else None,
                'divergence': parts[1],
                'deletion': parts[2],
                'insertion': parts[3]
            }
        )
    
    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream RepeatMasker features."""
        with open(self.filepath, 'r') as f:
            # Skip header lines
            for _ in range(3):
                next(f)
                
            for line in f:
                feature = self._parse_line(line)
                if feature:
                    if chrom and feature.chrom != chrom:
                        continue
                    if start and feature.end < start:
                        continue
                    if end and feature.start > end:
                        continue
                    yield feature


class UnifiedGenomeAnnotations:
    """Unified interface for all genomic annotations and sequences.
    
    Merges multiple annotation streams and sequence data efficiently.
    """
    
    def __init__(self, fasta_path: Optional[str] = None, vcf_path: Optional[str] = None):
        self.streams = {}
        self.indices = {}  # For fast range queries
        self.sequence_stream = None
        
        # Initialize sequence streaming if FASTA provided
        if fasta_path:
            self.setup_sequence_stream(fasta_path, vcf_path)
            # print("annotations set up sequence stream")
        
    def setup_sequence_stream(self, fasta_path: str, vcf_path: Optional[str] = None,
                            min_qual: float = 5.0):
        """Setup sequence streaming with optional variant integration.
        
        Args:
            fasta_path: Path to FASTA file
            vcf_path: Optional path to VCF file for variants
            min_qual: Minimum quality threshold for variants
        """
        from .sequence_stream import FASTAStream, SequenceStreamWithVariants
        
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
        
    def add_gtf(self, filepath: str, name: str = "genes"):
        """Add GTF annotation source."""
        self.add_source(name, GTFStream(filepath))
        
    def add_vcf(self, filepath: str, name: str = "variants"):
        """Add VCF annotation source."""
        self.add_source(name, VCFStream(filepath))
        
    def add_bed(self, filepath: str, name: str, feature_type: str = "region"):
        """Add BED annotation source."""
        self.add_source(name, BEDStream(filepath, feature_type))
        
    def add_repeatmasker(self, filepath: str, name: str = "repeats"):
        """Add RepeatMasker annotation source."""
        self.add_source(name, RepeatMaskerStream(filepath))
    
    def stream_all(self, chrom: Optional[str] = None,
                   start: Optional[int] = None,
                   end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream all annotations merged by position.
        
        Uses heapq.merge for efficient sorted merging.
        """
        # Create iterators for each source
        iterators = []
        for name, stream in self.streams.items():
            try:
                iterator = stream.stream(chrom, start, end)
                # Wrap each iterator to handle exhaustion
                iterators.append(iterator)
            except Exception as e:
                logger.warning(f"Failed to create stream for {name}: {e}")
                
        # Merge all iterators by position
        for feature in heapq.merge(*iterators):
            yield feature
    
    def query_point(self, chrom: str, position: int) -> List[UnifiedFeature]:
        """Query all features at a specific position."""
        features = []
        for name, stream in self.streams.items():
            try:
                for feature in stream.query_range(chrom, position, position):
                    features.append(feature)
            except Exception as e:
                logger.warning(f"Query failed for {name}: {e}")
        return sorted(features)
    
    def query_range(self, chrom: str, start: int, end: int) -> List[UnifiedFeature]:
        """Query all features in a range."""
        features = []
        for name, stream in self.streams.items():
            try:
                features.extend(stream.query_range(chrom, start, end))
            except Exception as e:
                logger.warning(f"Query failed for {name}: {e}")
        return sorted(features)
    
    def stream_by_types(self, feature_types: List[str],
                       chrom: Optional[str] = None,
                       start: Optional[int] = None,
                       end: Optional[int] = None) -> Iterator[UnifiedFeature]:
        """Stream only specific feature types."""
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


class CachedUnifiedAnnotations(UnifiedGenomeAnnotations):
    """Cached version for better performance."""
    
    def __init__(self, cache_size: int = 10000):
        super().__init__()
        self.cache = {}
        self.cache_size = cache_size
        
    def query_range(self, chrom: str, start: int, end: int) -> List[UnifiedFeature]:
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