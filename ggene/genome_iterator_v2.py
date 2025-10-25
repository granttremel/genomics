"""Enhanced Genome Iterator with unified streaming and multi-coordinate tracking.

This module provides an iterator over genomic positions with:
- Preloading of sequences and variants
- Three coordinate systems (reference, alternate, display)
- Integration with UnifiedGenomeAnnotations
- Efficient buffering for large indels
- Integrated motif detection
"""

import logging
from typing import Dict, List, Optional, Tuple, Union, Iterator, Any
from dataclasses import dataclass, field
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel("CRITICAL")

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


class UnifiedGenomeIterator:
    """Enhanced genome iterator using unified streaming system."""
    
    # Buffer sizes
    PRELOAD_BUFFER = 500  # Preload 500bp upstream/downstream for variants
    SEQUENCE_CACHE_SIZE = 10000  # Cache 10kb of sequence
    
    def __init__(self,
                 genome_manager,
                 chrom: Union[str, int],
                 start: int,
                 end: Optional[int] = None,
                 window_size: int = 100,
                 stride: Optional[int] = None,
                 integrate_variants: bool = True,
                 track_features: bool = True,
                 feature_types: Optional[List[str]] = None,
                 detect_motifs: bool = True):
        """Initialize the unified genome iterator.
        
        Args:
            genome_manager: GenomeManager with UnifiedGenomeAnnotations
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
        self.gm = genome_manager
        self.annotations = genome_manager.annotations
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
        self._buffer_start = start - self.PRELOAD_BUFFER
        self._buffer_end = start + self.SEQUENCE_CACHE_SIZE
        self._ref_buffer = ""
        self._alt_buffer = ""
        self._variant_list = []  # List of (position, delta) tuples
        self._cumulative_delta = 0  # Track cumulative variant delta
        
        # Cache for features
        self._feature_cache = {}
        
        # Motif detection
        self._motif_cache = {}  # Cache detected motifs by position
        self._motif_detector = None
        if detect_motifs:
            self._setup_motif_detection()
        
        # Preload initial buffer
        self._preload_buffer()
    
    def _setup_motif_detection(self) -> None:
        """Setup motif detector with default motifs."""
        try:
            from ggene.motifs.motif import MotifDetector
            from ggene.motifs.pattern import PatternMotif
            
            self._motif_detector = MotifDetector()
            
            # Add common splice site motifs (compile with IGNORECASE)
            import re
            splice_donor = PatternMotif(
                "splice_donor", 
                re.compile(r"GT[AG]AGT", re.IGNORECASE),  # GT(A/G)AGT
                lambda x: 1.0
            )
            splice_acceptor = PatternMotif(
                "splice_acceptor",
                re.compile(r"[CT]{10,}[ATGC]CAG", re.IGNORECASE),  # Pyrimidine tract + CAG
                lambda x: 1.0
            )
            
            # Add TATA box motif
            tata_box = PatternMotif(
                "TATA_box",
                re.compile(r"TATA[AT]A[AT]", re.IGNORECASE),
                lambda x: 1.0
            )
            
            # Add CpG island pattern (simplified)
            cpg_pattern = PatternMotif(
                "CpG_island",
                re.compile(r"(CG){5,}", re.IGNORECASE),  # At least 5 CG dinucleotides
                lambda x: len(x) / 2  # Score by number of CG repeats
            )
            
            # Add poly-A signal
            polya_signal = PatternMotif(
                "polyA_signal",
                re.compile(r"AATAAA|ATTAAA", re.IGNORECASE),
                lambda x: 1.0
            )
            
            # Add motifs to detector
            self._motif_detector.add_motif(splice_donor)
            self._motif_detector.add_motif(splice_acceptor)
            self._motif_detector.add_motif(tata_box)
            self._motif_detector.add_motif(cpg_pattern)
            self._motif_detector.add_motif(polya_signal)
            
            logger.debug(f"Initialized motif detector with {len(self._motif_detector.motifs)} motifs")
            
        except ImportError:
            logger.warning("Motif modules not available, disabling motif detection")
            self.detect_motifs = False
            self._motif_detector = None
    
    def _scan_buffer_for_motifs(self, seq: str, buffer_start: int) -> Dict[int, List]:
        """Scan a buffer sequence for motifs and cache results.
        
        Args:
            seq: Sequence to scan
            buffer_start: Genomic start position of the buffer
            
        Returns:
            Dictionary of position -> motif list
        """
        if not self._motif_detector or not seq:
            return {}
        
        motif_results = {}
        
        # Helper function to get reverse complement
        def reverse_complement_seq(s):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                         'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                         'N': 'N', 'n': 'n'}
            return ''.join(complement.get(base, base) for base in s[::-1])
        
        # Scan forward strand
        for motif_name, motif in self._motif_detector.motifs.items():
            try:
                instances = motif.find_instances(seq)
                
                for start, end, score in instances:
                    # Convert to genomic coordinates
                    genomic_start = buffer_start + start
                    genomic_end = buffer_start + end
                    
                    # Create motif feature object
                    motif_feature = {
                        'type': 'motif',
                        'name': motif_name,
                        'start': genomic_start,
                        'end': genomic_end,
                        'score': score,
                        'sequence': seq[start:end],
                        'strand': '+'  # Forward strand
                    }
                    
                    # Add to cache by position
                    if genomic_start not in motif_results:
                        motif_results[genomic_start] = []
                    motif_results[genomic_start].append(motif_feature)
                    
            except Exception as e:
                logger.debug(f"Failed to scan for motif {motif_name}: {e}")
        
        # Scan reverse complement strand
        rev_seq = reverse_complement_seq(seq)
        for motif_name, motif in self._motif_detector.motifs.items():
            try:
                instances = motif.find_instances(rev_seq)
                
                for start, end, score in instances:
                    # Convert positions - reverse strand positions need adjustment
                    # Position on reverse strand maps to original sequence
                    rev_start = len(seq) - end
                    rev_end = len(seq) - start
                    
                    # Convert to genomic coordinates
                    genomic_start = buffer_start + rev_start
                    genomic_end = buffer_start + rev_end
                    
                    # Create motif feature object
                    motif_feature = {
                        'type': 'motif',
                        'name': motif_name,
                        'start': genomic_start,
                        'end': genomic_end,
                        'score': score,
                        'sequence': seq[rev_start:rev_end],  # Original sequence at this position
                        'strand': '-'  # Reverse strand
                    }
                    
                    # Add to cache by position
                    if genomic_start not in motif_results:
                        motif_results[genomic_start] = []
                    motif_results[genomic_start].append(motif_feature)
                    
            except Exception as e:
                logger.debug(f"Failed to scan reverse strand for motif {motif_name}: {e}")
        
        return motif_results
    
    def _preload_buffer(self) -> None:
        """Preload sequence and variant buffer, and scan for motifs."""
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
                # Update cache with new motif results
                self._motif_cache.update(motif_results)
                logger.debug(f"Found {len(motif_results)} motif positions in buffer")
        
        logger.debug(f"Preloaded buffer {buffer_start}-{buffer_end}: "
                    f"{len(self._ref_buffer)} ref bases, {len(self._variant_list)} variants")
    
    def _get_features_in_window(self, start: int, end: int) -> List[Any]:
        """Get features in a genomic window."""
        if not self.track_features:
            return []
        
        # Check cache
        cache_key = (start, end)
        if cache_key in self._feature_cache:
            return self._feature_cache[cache_key]
        
        # Query features
        features = self.annotations.query_range(self.chrom, start, end)
        
        # Filter by type if specified
        if self.feature_types:
            features = [f for f in features if f.feature_type in self.feature_types]
        
        # Cache result
        self._feature_cache[cache_key] = features
        return features
    
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
        if window_end_ref > self._buffer_end:
            self._buffer_start = window_start_ref - self.PRELOAD_BUFFER
            self._buffer_end = window_end_ref + self.SEQUENCE_CACHE_SIZE
            self._preload_buffer()
        
        # Calculate position in buffer
        buffer_offset_start = window_start_ref - max(1, self._buffer_start)
        buffer_offset_end = window_end_ref - max(1, self._buffer_start)
        
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
        
        # Get features
        features = self._get_features_in_window(window_start_ref, window_end_ref)
        
        # Get variant features separately (these are UnifiedFeature objects)
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
        
        # Get motifs in this window from cache
        window_motifs = []
        if self.detect_motifs and self._motif_cache:
            for pos in range(window_start_ref, window_end_ref + 1):
                if pos in self._motif_cache:
                    # Filter motifs that actually overlap with the window
                    for motif in self._motif_cache[pos]:
                        if motif['start'] < window_end_ref and motif['end'] > window_start_ref:
                            # Only add if not already in list (avoid duplicates)
                            if motif not in window_motifs:
                                window_motifs.append(motif)
            
            # Sort motifs by position for display
            window_motifs.sort(key=lambda m: m['start'])
            logger.debug(f"Found {len(window_motifs)} motifs in window {window_start_ref}-{window_end_ref}")
        
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
            motifs=window_motifs if window_motifs else None
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
    
    def get_window_at(self, position):
        
        return self._extract_window(position, position + self.window_size)
        
    
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