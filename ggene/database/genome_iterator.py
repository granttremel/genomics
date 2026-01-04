"""Enhanced Genome Iterator with unified streaming and multi-coordinate tracking.

This module provides an iterator over genomic positions with:
- Preloading of sequences and variants
- Three coordinate systems (reference, alternate, display)
- Integration with UGenomeAnnotations
- Efficient buffering for large indels
- Integrated motif detection
"""

import logging
from typing import Dict, List, Optional, Tuple, Union, Iterator, Any, TYPE_CHECKING
from dataclasses import dataclass, field
import numpy as np
import reprlib
import asyncio
from concurrent.futures import ThreadPoolExecutor

from ggene.seqs.bio import reverse_complement
from ggene.database.annotations import UFeature

if TYPE_CHECKING:
    from .genome_manager import GenomeManager
    from .annotations import UGenomeAnnotations

logger = logging.getLogger(__name__)
logger.setLevel("CRITICAL")
# logger.setLevel(logging.DEBUG)

@dataclass
class BaseWindow:
    
    pass

class BaseIterator:
    
    def __init__(self):
        pass
    
    def get_window(self)->BaseWindow:
        pass
    
    def get_window_at(self, position)->BaseWindow:
        pass
    
    def step(self, delta: int) -> None:
        pass
    
    def update(self, **kwargs):
        pass
    
    def cleanup(self):
        pass
    
    def __iter__(self):
        pass

@dataclass
class CoordinateState:
    """Tracks position in three coordinate systems."""
    ref_pos: int      # Position in reference genome
    alt_pos: int      # Position in alternate genome (with variants)
    display_pos: int  # Position for display (monotonically increasing)
    
    def advance_snp(self):
        """Advance all coordinates by 1 (SNP or match)."""
        self.ref_pos += 1
        self.alt_pos += 1
        self.display_pos += 1
    
    def advance_deletion(self, del_length: int):
        """Advance coordinates for a deletion."""
        self.ref_pos += del_length  # Reference advances
        # alt_pos stays same (bases deleted)
        self.display_pos += del_length  # Display shows gaps
    
    def advance_insertion(self, ins_length: int):
        """Advance coordinates for an insertion."""
        # ref_pos stays same (no ref bases)
        self.alt_pos += ins_length  # Alt advances
        self.display_pos += ins_length  # Display shows insertions


@dataclass 
class GenomeWindow:
    """Represents a window of genome with aligned sequences and features."""
    chrom: str
    start_ref: int          # Start in reference coordinates
    end_ref: int            # End in reference coordinates  
    start_alt: int          # Start in alternate coordinates
    end_alt: int            # End in alternate coordinates
    start_display: int      # Start in display coordinates
    end_display: int        # End in display coordinates
    ref_seq: str            # Reference sequence (may have gaps)
    alt_seq: str            # Alternate sequence (may have gaps)
    features: List[Any]     # Features in this window
    variant_deltas: List[Tuple[int, int]]  # List of (position, delta) tuples for coordinate tracking
    variant_features: List[Any] = None  # Actual variant feature objects
    motifs: List[Any] = None  # Detected motifs in this window
    
    @property
    def has_gaps(self) -> bool:
        """Check if this window has gaps from indels."""
        return '-' in self.ref_seq or '-' in self.alt_seq
    
    @property
    def ref_length(self) -> int:
        """Length of reference sequence (excluding gaps)."""
        return len(self.ref_seq.replace('-', ''))
    
    @property
    def alt_length(self) -> int:
        """Length of alternate sequence (excluding gaps)."""
        return len(self.alt_seq.replace('-', ''))
    
    @property
    def display_length(self) -> int:
        """Length for display (including gaps)."""
        return max(len(self.ref_seq), len(self.alt_seq))


class UGenomeIterator:
    """Enhanced genome iterator using unified streaming system."""
    
    # Buffer sizes
    PRELOAD_BUFFER = 500  # Preload 500bp upstream/downstream for variants
    SEQUENCE_CACHE_SIZE = 10000  # Cache 10kb of sequence
    
    def __init__(self,
                 genome_manager,
                 chrom: Union[str, int],
                 start: int,
                 end: Optional[int] = None,
                 window_size: int = 240,
                 stride: Optional[int] = None,
                 integrate_variants: bool = True,
                 track_features: bool = True,
                 feature_types: Optional[List[str]] = None,
                 detect_motifs: bool = True):
        """Initialize the unified genome iterator.
        
        Args:
            genome_manager: GenomeManager with UGenomeAnnotations
            chrom: Chromosome to iterate
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            window_size: Size of windows to return
            stride: Step size (defaults to window_size)
            integrate_variants: Whether to apply variants
            track_features: Whether to track features
            feature_types: Specific feature types to track
            detect_motifs: Whether to detect motifs in sequences
        """
        self.gm:'GenomeManager' = genome_manager
        self.annotations:'UGenomeAnnotations' = genome_manager.annotations
        self.chrom = str(chrom)
        self.start = start
        self.end = end
        self.window_size = window_size
        self.stride = stride or window_size
        self.integrate_variants = integrate_variants
        self.track_features = track_features
        self.feature_types = feature_types
        self.detect_motifs = detect_motifs
        
        # Coordinate tracking
        self.coords = CoordinateState(
            ref_pos=start,
            alt_pos=start,
            display_pos=start
        )
        
        # Buffered sequences and variants
        # self._buffer_start = start - self.PRELOAD_BUFFER
        self._buffer_start = start - self.SEQUENCE_CACHE_SIZE//2
        self._buffer_end = start + window_size + self.SEQUENCE_CACHE_SIZE//2
        self._ref_buffer = ""
        self._alt_buffer = ""
        self._variant_list = []  # List of (position, delta) tuples
        self._cumulative_delta = 0  # Track cumulative variant delta
        
        # Cache for features
        self._feature_cache = {}
        
        # Motif detection
        self._motif_cache = {}  # Cache detected motifs by position
        # self._motif_detector = None
        # if detect_motifs:
        #     self._motif_detector = self.gm.motif_detector

            # self._setup_motif_detection()

        # Async preloading infrastructure
        self._executor = ThreadPoolExecutor(max_workers=1)
        self._preload_future = None
        self._next_buffer_start = None
        self._next_buffer_end = None
        self._is_preloading = False

        # Preload initial buffer
        self._preload_buffer()
    
    def _scan_buffer_for_motifs(self, seq: str, buffer_start: int) -> Dict[int, List]:
        """Scan a buffer sequence for motifs and cache results.
        
        Args:
            seq: Sequence to scan
            buffer_start: Genomic start position of the buffer
            
        Returns:
            Dictionary of position -> motif list
        """
        if not seq:
            return {}
        
        motif_results = {}
        
        strand = "+"
        
        for src, mtf_strm in self.gm.annotations.motif_streams.items():
            
            mtfs = mtf_strm.scan_sequence(seq, self.chrom, buffer_start, strand = strand)
            
            for f in mtfs:
                
                genomic_start = f.start
                
                if not genomic_start in motif_results:
                    motif_results[genomic_start] = []
                
                motif_results[genomic_start].append(f)
        
        return motif_results
    #     for motif_name, instances in all_insts.items():
    #         if instances:
    #             # for motif_start, motif_end, score, is_rc, mtf_cls in instances:
    #             for motif in instances:
                    
    #                 motif_start = motif.get("start", 0)
    #                 motif_end = motif.get("end", 0)
    #                 score = motif.get("score", 0)
    #                 is_rc = motif.get("is_rc", False)
    #                 mtf_cls = motif.get("class", "")
                    
    #                 genomic_start = buffer_start + motif_start
                    
    #                 if not genomic_start in motif_results:
    #                     motif_results[genomic_start] = []
                    
    #                 mseq = seq[motif_start:motif_end]
    #                 mtfstrand = strand
    #                 if is_rc:
    #                     # mseq = reverse_complement(mseq)
    #                     mtfstrand = '-' if mtfstrand == '+' else '+'
                    
    #                 motif_results[genomic_start].append(UFeature(
    #                     chrom=self.chrom,
    #                     start=buffer_start + motif_start,
    #                     end=buffer_start + motif_end - 1,
    #                     feature_type="motif",
    #                     source="MotifDetector",
    #                     score=score,
    #                     strand=mtfstrand,
    #                     name=motif_name,
    #                     attributes={
    #                         'sequence': mseq,
    #                         'is_rc':is_rc,
    #                         'class':mtf_cls,
    #                         'caller':'iterator',
    #                     }
    #                 ))

    #     return motif_results
    
    def _preload_buffer(self) -> None:
        """Preload sequence and variant buffer, and scan for motifs (synchronous version for initial load)."""
        # Ensure buffer boundaries
        buffer_start = max(1, self._buffer_start)
        buffer_end = self._buffer_end

        if self.integrate_variants and self.annotations.sequence_stream:
            # Get aligned sequences with variant information
            aligned_ref, aligned_alt, variant_list = self.annotations.get_aligned_sequences(
                self.chrom, buffer_start, buffer_end
            )
            self._ref_buffer = aligned_ref
            self._alt_buffer = aligned_alt
            self._variant_list = variant_list

            # Calculate cumulative delta up to current position
            self._cumulative_delta = sum(
                delta for pos, delta in variant_list
                if pos < self.coords.ref_pos
            )
        else:
            # No variants - just get reference
            ref_seq = self.annotations.get_sequence(self.chrom, buffer_start, buffer_end)
            self._ref_buffer = ref_seq
            self._alt_buffer = ref_seq
            self._variant_list = []
            self._cumulative_delta = 0

        # Scan for motifs if enabled
        if self.detect_motifs and self._ref_buffer:
            # Get clean sequence without gaps for motif scanning
            clean_seq = self._ref_buffer.replace('-', '')
            if clean_seq:
                motif_results = self._scan_buffer_for_motifs(clean_seq, buffer_start)
                # motif_results = {}
                # Update cache with new motif results
                self._motif_cache.update(motif_results)
                logger.debug(f"Found {len(motif_results)} motif positions in buffer")

        logger.debug(f"Preloaded buffer {buffer_start}-{buffer_end}: "
                    f"{len(self._ref_buffer)} ref bases, {len(self._variant_list)} variants")

    def _load_buffer_async_worker(self, buffer_start: int, buffer_end: int) -> Tuple[str, str, List, int]:
        """Worker function that runs in thread pool to load buffer data."""
        buffer_start = max(1, buffer_start)

        if self.integrate_variants and self.annotations.sequence_stream:
            # Get aligned sequences with variant information
            aligned_ref, aligned_alt, variant_list = self.annotations.get_aligned_sequences(
                self.chrom, buffer_start, buffer_end
            )

            # Calculate cumulative delta
            cumulative_delta = sum(
                delta for pos, delta in variant_list
                if pos < self.coords.ref_pos
            )
        else:
            # No variants - just get reference
            ref_seq = self.annotations.get_sequence(self.chrom, buffer_start, buffer_end)
            aligned_ref = ref_seq
            aligned_alt = ref_seq
            variant_list = []
            cumulative_delta = 0

        # Scan for motifs if enabled
        motif_cache = {}
        if self.detect_motifs and aligned_ref:
            # Get clean sequence without gaps for motif scanning
            clean_seq = aligned_ref.replace('-', '')
            if clean_seq:
                motif_cache = self._scan_buffer_for_motifs(clean_seq, buffer_start)
                # motif_cache = {}
                logger.debug(f"Found {len(motif_cache)} motif positions in async buffer")

        logger.debug(f"Async preloaded buffer {buffer_start}-{buffer_end}: "
                    f"{len(aligned_ref)} ref bases, {len(variant_list)} variants")

        return aligned_ref, aligned_alt, variant_list, cumulative_delta, motif_cache

    def _start_async_preload(self, buffer_start: int, buffer_end: int) -> None:
        """Start asynchronously preloading the next buffer in the background."""
        if self._is_preloading:
            logger.debug("Already preloading, skipping")
            return

        self._is_preloading = True
        self._next_buffer_start = buffer_start
        self._next_buffer_end = buffer_end

        logger.debug(f"Starting async preload for {buffer_start}-{buffer_end}")
        self._preload_future = self._executor.submit(
            self._load_buffer_async_worker, buffer_start, buffer_end
        )

    def _wait_for_preload(self) -> None:
        """Wait for async preload to complete and update buffers."""
        if not self._is_preloading or self._preload_future is None:
            return

        try:
            # Wait for the preload to complete and get results
            aligned_ref, aligned_alt, variant_list, cumulative_delta, motif_cache = self._preload_future.result()

            # Update the buffer state
            self._buffer_start = self._next_buffer_start
            self._buffer_end = self._next_buffer_end
            self._ref_buffer = aligned_ref
            self._alt_buffer = aligned_alt
            self._variant_list = variant_list
            self._cumulative_delta = cumulative_delta

            # Update motif cache
            self._motif_cache.update(motif_cache)

            logger.debug(f"Async preload completed and applied")
        finally:
            self._is_preloading = False
            self._preload_future = None
            self._next_buffer_start = None
            self._next_buffer_end = None
    
    def _get_features_in_window(self, start: int, end: int) -> List[Any]:
        """Get features in a genomic window."""
        if not self.track_features:
            return []
        
        # Check cache
        cache_key = (start, end)
        if cache_key in self._feature_cache:
            return self._feature_cache[cache_key]
        

        streams = list(self.annotations.streams.keys())
        logger.debug(f"available streams: {streams}")
        
        # Query features
        # features = self.annotations.query_range(self.chrom, start, end)
        # features = [f for f in self.annotations.stream_all(self.chrom, start, end)]
        features = [f for f in self.annotations.stream_by_types(self.feature_types,self.chrom, start, end)]
        
        fts = set()
        for f in features:
            fts.add(f.feature_type)
        logger.debug(f"feature types in window: {list(fts)}")
        
        # Filter by type if specified
        if self.feature_types:
            features = [f for f in features if f.feature_type in self.feature_types]
        
        # Cache result
        self._feature_cache[cache_key] = features
        return features
    
    def _get_motifs_in_window(self, start:int, end:int):
        if not self.detect_motifs:
            logger.debug("skipping motifs")
            return []
        
        # logger.debug("detecting motifs")
        
        # Check cache
        # cache_key = (start, end)
        # if cache_key in self._feature_cache:
        #     return self._feature_cache[cache_key]
        
        # Query motifs
        motifs = self.annotations.query_motifs(self.chrom, start, end)
        
        # logger.debug(f"identified {len(motifs)} motifs")
        
        # Filter by type if specified
        if self.feature_types:
            motifs = [f for f in motifs if f.feature_type in self.feature_types]
        
        # Cache result
        # self._feature_cache[cache_key] = features
        return motifs
        
    def _calculate_display_coordinate(self, ref_pos: int, alt_pos: int) -> int:
        """Calculate display coordinate based on variant deltas.
        
        Display coordinate shows both insertions and deletions:
        - For insertions: display = ref + sum(positive deltas)
        - For deletions: display = alt + sum(|negative deltas|)
        """
        # Get all variants up to this position
        positive_delta = sum(
            delta for pos, delta in self._variant_list
            if pos <= ref_pos and delta > 0
        )
        negative_delta = sum(
            abs(delta) for pos, delta in self._variant_list
            if pos <= ref_pos and delta < 0
        )
        
        # Display coordinate includes space for both insertions and deletions
        return ref_pos + positive_delta + negative_delta
    
    def _extract_window(self, window_start_ref: int, window_end_ref: int) -> GenomeWindow:
        """Extract a window from the buffer."""
        # Check if we need to reload buffer
        logger.debug(f"buffer offsets before: {self._buffer_start}, {self._buffer_end}")

        # Check if the window extends beyond current buffer (either direction)
        if window_end_ref > self._buffer_end or window_start_ref < self._buffer_start:
            # If we have a preload in progress for this region, wait for it
            if self._is_preloading and self._next_buffer_start is not None:
                if window_start_ref >= self._next_buffer_start and window_end_ref <= self._next_buffer_end:
                    logger.debug("Using async preloaded buffer")
                    self._wait_for_preload()
                else:
                    # Wrong region, need to reload synchronously
                    direction = "backward" if window_start_ref < self._buffer_start else "forward"
                    logger.debug(f"Preload was for wrong region, reloading synchronously ({direction})")
                    # self._buffer_start = window_start_ref - self.PRELOAD_BUFFER
                    self._buffer_start = window_start_ref - self.SEQUENCE_CACHE_SIZE//2
                    self._buffer_end = window_end_ref + self.SEQUENCE_CACHE_SIZE//2
                    self._preload_buffer()
            else:
                # No preload in progress, load synchronously
                direction = "backward" if window_start_ref < self._buffer_start else "forward"
                logger.debug(f"No async preload available, loading synchronously ({direction})")
                # self._buffer_start = window_start_ref - self.PRELOAD_BUFFER
                self._buffer_start = window_start_ref - self.SEQUENCE_CACHE_SIZE//2
                self._buffer_end = window_end_ref + self.SEQUENCE_CACHE_SIZE//2
                self._preload_buffer()

        # Predictively start loading the next buffer if we're getting close to the end
        # Start async preload when we're 75% through the current buffer
        buffer_usage_fwd = (window_end_ref - self._buffer_start) / (self._buffer_end - self._buffer_start)
        buffer_usage_rev = (window_start_ref - self._buffer_start) / (self._buffer_end - self._buffer_start)
        if (buffer_usage_fwd > 0.9 or buffer_usage_rev < 0.1) and not self._is_preloading:
            # next_buffer_start = window_end_ref - self.PRELOAD_BUFFER
            next_buffer_start = window_start_ref - self.SEQUENCE_CACHE_SIZE//2
            next_buffer_end = window_end_ref + self.SEQUENCE_CACHE_SIZE//2
            logger.debug(f"Predictively starting async preload (buffer usage: {buffer_usage_fwd:.1%})")
            self._start_async_preload(next_buffer_start, next_buffer_end)

        logger.debug(f"window and consts: {window_start_ref}, {window_end_ref}, {self.PRELOAD_BUFFER}, {self.SEQUENCE_CACHE_SIZE}")
        
        # Calculate position in buffer
        buffer_offset_start = window_start_ref - max(1, self._buffer_start)
        buffer_offset_end = window_end_ref - max(1, self._buffer_start)
        logger.debug(f"buffer offsets after: {buffer_offset_start}, {buffer_offset_end}, {self._buffer_start}, {self._buffer_end}")
        
        # Account for gaps from variants
        # Count variants before and in window
        variants_before = [(p, d) for p, d in self._variant_list if p < window_start_ref]
        variants_in_window = [(p, d) for p, d in self._variant_list 
                              if window_start_ref <= p <= window_end_ref]
        
        # Calculate cumulative delta before window
        delta_before = sum(d for _, d in variants_before)
        
        # Adjust buffer positions for gaps
        # This is approximate - may need refinement based on actual gap positions
        adjusted_start = buffer_offset_start
        adjusted_end = buffer_offset_end + sum(abs(d) for _, d in variants_in_window)
        
        # Extract sequences (with gaps)
        ref_seq = self._ref_buffer[adjusted_start:adjusted_end] if self._ref_buffer else ""
        alt_seq = self._alt_buffer[adjusted_start:adjusted_end] if self._alt_buffer else ""
        
        logger.debug(f"extracted seqs with {len(ref_seq)} bases from {adjusted_start}-{adjusted_end}")
        if ref_seq:
            logger.debug(f"ref: {ref_seq[:32]}")
        logger.debug(f"has buffer?: {bool(self._ref_buffer)}, {bool(self._alt_buffer)}")
        
        # Get features
        features = self._get_features_in_window(window_start_ref, window_end_ref)
        
        # Get variant features separately (these are UFeature objects)
        variant_features = []
        if self.annotations and hasattr(self.annotations, 'streams') and 'variants' in self.annotations.streams:
            # Query the VCF stream for variant features in this window
            try:
                variant_stream = self.annotations.streams['variants']
                for var_feature in variant_stream.stream(self.chrom, window_start_ref, window_end_ref):
                    if var_feature.overlaps(window_start_ref, window_end_ref):
                        variant_features.append(var_feature)
            except Exception as e:
                logger.debug(f"Failed to get variant features: {e}")
        
        window_motifs = []
        # Get motifs in this window from cache
        # window_motifs = []
        if self.detect_motifs and self._motif_cache:
            window_motifs = self._get_motifs_in_window(window_start_ref, window_end_ref)
            
            for pos in range(window_start_ref, window_end_ref + 1):
                if pos in self._motif_cache:
                    # Filter motifs that actually overlap with the window
                    for motif in self._motif_cache[pos]:
                        if motif['start'] < window_end_ref and motif['end'] > window_start_ref:
                            # Only add if not already in list (avoid duplicates)
                            if motif not in window_motifs:
                                window_motifs.append(motif)
            logger.debug(f"found {len(window_motifs)} motifs in _extract_window")
        #     # Sort motifs by position for display
        #     window_motifs.sort(key=lambda m: m['start'])
        #     logger.debug(f"Found {len(window_motifs)} motifs in window {window_start_ref}-{window_end_ref}")
        
        # Calculate coordinate positions
        alt_start = window_start_ref + delta_before
        alt_end = alt_start + len(alt_seq.replace('-', ''))
        display_start = self._calculate_display_coordinate(window_start_ref, alt_start)
        display_end = display_start + max(len(ref_seq), len(alt_seq))
        
        return GenomeWindow(
            chrom=self.chrom,
            start_ref=window_start_ref,
            end_ref=window_end_ref,
            start_alt=alt_start,
            end_alt=alt_end,
            start_display=display_start,
            end_display=display_end,
            ref_seq=ref_seq,
            alt_seq=alt_seq,
            features=features,
            variant_deltas=variants_in_window,
            variant_features=variant_features,
            motifs=window_motifs
        )
    
    def __iter__(self) -> Iterator[GenomeWindow]:
        """Iterate over genome windows."""
        current_ref_pos = self.start
        
        while True:
            # Check if we've reached the end
            if self.end and current_ref_pos > self.end:
                break
            
            # Calculate window boundaries
            window_end_ref = min(current_ref_pos + self.window_size - 1,
                                self.end) if self.end else current_ref_pos + self.window_size - 1
            
            # Extract and yield window
            window = self._extract_window(current_ref_pos, window_end_ref)
            yield window
            
            # Advance position
            current_ref_pos += self.stride
            
            # Update coordinate state
            self.coords.ref_pos = current_ref_pos
            # Update alt_pos based on variants encountered
            variants_passed = [(p, d) for p, d in self._variant_list if p < current_ref_pos]
            self.coords.alt_pos = current_ref_pos + sum(d for _, d in variants_passed)
            self.coords.display_pos = self._calculate_display_coordinate(
                current_ref_pos, self.coords.alt_pos
            )
    
    def get_window(self):
        return self._extract_window(self.start, self.start + self.window_size)
    
    def get_window_at(self, position):
        return self._extract_window(position, position + self.window_size)

    def step(self, delta: int) -> None:
        """Move the iterator position by delta bases.

        This is optimized for sequential navigation and cooperates with
        the async buffer preloading system.

        Args:
            delta: Number of bases to move (positive = forward, negative = backward)
        """
        new_position = max(1, self.start + delta)
        self.start = new_position

        # Update coordinate state
        self.coords.ref_pos = new_position

        # Update alt_pos based on variants encountered
        variants_passed = [(p, d) for p, d in self._variant_list if p < new_position]
        self.coords.alt_pos = new_position + sum(d for _, d in variants_passed)
        self.coords.display_pos = self._calculate_display_coordinate(
            new_position, self.coords.alt_pos
        )

        logger.debug(f"Stepped by {delta} to position {new_position}")

    def update(self, chrom = None, position = None, window_size = None, stride = None, **kwargs):
        
        if chrom:
            self.update_chromosome(chrom, new_position = position)
        
        if position:
            self.update_position(position)
        
        if window_size:
            self.window_size = window_size
        
        if stride:
            self.stride = stride
            
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
        
            
    def update_position(self, new_position: int) -> None:
        """Update the iterator to a new position without jumping.

        This is similar to step() but takes an absolute position instead of a delta.
        Optimized for sequential navigation and preserves async preloading benefits.

        Args:
            new_position: New absolute position (1-based)
        """
        new_position = max(1, new_position)
        self.start = new_position

        # Update coordinate state
        self.coords.ref_pos = new_position

        # Update alt_pos based on variants encountered
        variants_passed = [(p, d) for p, d in self._variant_list if p < new_position]
        self.coords.alt_pos = new_position + sum(d for _, d in variants_passed)
        self.coords.display_pos = self._calculate_display_coordinate(
            new_position, self.coords.alt_pos
        )

        logger.debug(f"Updated position to {new_position}")

    def update_chromosome(self, new_chrom: Union[str, int], new_position: int = None) -> None:
        """Update the iterator to a new chromosome and optionally a new position.

        This requires reloading buffers since we're changing chromosomes.

        Args:
            new_chrom: New chromosome to navigate to
            new_position: Optional new position (defaults to current position or 1)
        """
        # Cancel any ongoing preload
        if self._is_preloading:
            logger.debug("Canceling async preload due to chromosome change")
            self._is_preloading = False
            self._preload_future = None

        self.chrom = str(new_chrom)
        self.start = max(1, new_position if new_position is not None else self.start)

        # Reset coordinate state
        self.coords.ref_pos = self.start
        self.coords.alt_pos = self.start
        self.coords.display_pos = self.start

        # Clear caches
        self._feature_cache = {}
        self._motif_cache = {}

        # Reload buffer for new chromosome
        self._buffer_start = self.start - self.PRELOAD_BUFFER
        self._buffer_end = self.start + self.SEQUENCE_CACHE_SIZE
        self._preload_buffer()

        logger.debug(f"Updated chromosome to {new_chrom} at position {self.start}")


    def jump_to(self, position: int, coord_system: str = 'ref') -> None:
        """Jump to a specific position in the genome.
        
        Args:
            position: Target position
            coord_system: Which coordinate system ('ref', 'alt', or 'display')
        """
        if coord_system == 'ref':
            self.coords.ref_pos = position
            # Recalculate other coordinates
            variants_before = [(p, d) for p, d in self._variant_list if p < position]
            delta = sum(d for _, d in variants_before)
            self.coords.alt_pos = position + delta
            self.coords.display_pos = self._calculate_display_coordinate(
                position, self.coords.alt_pos
            )
        elif coord_system == 'alt':
            # This is more complex - need to find corresponding ref position
            # For now, approximate
            self.coords.alt_pos = position
            self.coords.ref_pos = position - self._cumulative_delta
            self.coords.display_pos = self._calculate_display_coordinate(
                self.coords.ref_pos, position
            )
        elif coord_system == 'display':
            # Even more complex - would need reverse mapping
            raise NotImplementedError("Jump to display coordinate not yet implemented")
        else:
            raise ValueError(f"Unknown coordinate system: {coord_system}")
        
        # Reload buffer if needed
        if self.coords.ref_pos < self._buffer_start or self.coords.ref_pos > self._buffer_end:
            # Cancel any ongoing preload since we're jumping
            if self._is_preloading:
                logger.debug("Canceling async preload due to jump")
                # Don't wait for it, just mark as not preloading
                self._is_preloading = False
                self._preload_future = None

            self._buffer_start = self.coords.ref_pos - self.PRELOAD_BUFFER
            self._buffer_end = self.coords.ref_pos + self.SEQUENCE_CACHE_SIZE
            self._preload_buffer()
    
    def get_current_coordinates(self) -> Dict[str, int]:
        """Get current position in all coordinate systems."""
        return {
            'reference': self.coords.ref_pos,
            'alternate': self.coords.alt_pos,
            'display': self.coords.display_pos,
            'cumulative_delta': self._cumulative_delta
        }

    def cleanup(self):
        """Clean up resources (thread pool, etc.)."""
        if self._executor:
            logger.debug("Shutting down thread pool executor")
            self._executor.shutdown(wait=False)
            self._executor = None

    def __del__(self):
        """Cleanup when iterator is destroyed."""
        self.cleanup()

    def __repr__(self):
        safe_attrs = {k: v for k, v in self.__dict__.items()
                    if isinstance(v, (int, float, str, bool, type(None)))}
        return f"{self.__class__.__name__}({reprlib.repr(safe_attrs)})"