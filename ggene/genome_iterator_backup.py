import logging
from typing import Dict, List, Optional, Tuple, Union, Iterator, Any, Set, TYPE_CHECKING
from dataclasses import dataclass
import numpy as np

if TYPE_CHECKING:
    from .genomemanager import GenomeManager

logger = logging.getLogger(__name__)

all_features = ["gene","transcript","exon","CDS","intron","five_prime_utr","start_codon","stop_codon","three_prime_utr","variant"]

@dataclass
class GenomicPosition:
    """Represents a single position in the genome with associated features and variants."""
    chrom: Union[str, int]
    position: int
    reference_base: str
    personal_base: str
    features: List[Any]
    variants: List[Any]
    quality_scores: List[float]
    
    @property
    def has_variant(self) -> bool:
        """Check if this position has a variant."""
        return self.reference_base != self.personal_base
    
    @property
    def feature_types(self) -> Set[str]:
        """Get unique feature types at this position."""
        return {f.type for f in self.features if hasattr(f, 'type')}
    
    def get_features_by_type(self, feature_type: str) -> List[Any]:
        """Get all features of a specific type at this position."""
        return [f for f in self.features if hasattr(f, 'type') and f.type == feature_type]


class GenomeIterator:
    """Iterator over genomic positions with feature and variant integration."""
    
    def __init__(self, genome_manager: 'GenomeManager', 
                 chrom: Union[str, int],
                 start: int,
                 end: Optional[int] = None,
                 stride: int = 1,
                 integrate_variants: bool = True,
                 min_variant_quality: float = 5.0,
                 track_features: bool = True,
                 feature_types: Optional[List[str]] = all_features,
                 window_size: Optional[int] = None):
        """Initialize the genome iterator.
        
        Args:
            genome_manager: GenomeManager instance with loaded data
            chrom: Chromosome to iterate over
            start: Starting position (1-based)
            end: Ending position (1-based, inclusive). If None, iterate to chromosome end
            stride: Number of bases to advance each iteration
            integrate_variants: Whether to apply personal variants to sequence
            min_variant_quality: Minimum quality threshold for variants
            track_features: Whether to track features at each position
            feature_types: List of feature types to track (None = all types)
            window_size: If set, return windows of this size instead of single bases
        """
        self.gm = genome_manager
        self.chrom = chrom
        self.start = start
        self.end = end
        self.stride = stride
        self.integrate_variants = integrate_variants
        self.min_variant_quality = min_variant_quality
        self.track_features = track_features
        self.feature_types = feature_types
        self.window_size = window_size or 1
        
        # Current position
        self.current_pos = start
        
        # Cache for efficiency
        self._variant_cache = {}
        self._feature_cache = {}
        self._sequence_cache = {}
        self._cache_window = 10000  # Cache 10kb windows
        
        # Preload variants if integrating
        if self.integrate_variants:
            self._preload_variants()
        
        # Preload features if tracking
        if self.track_features:
            self._preload_features()
    
    def _preload_variants(self) -> None:
        """Preload variants for the region."""
        try:
            # Calculate region to preload
            cache_start = self.current_pos
            cache_end = min(self.current_pos + self._cache_window, 
                          self.end) if self.end else self.current_pos + self._cache_window
            
            # Query variants
            variants = []
            for var in self.gm.vcf(self.gm._make_index(self.chrom, cache_start, cache_end)):
                if var.QUAL >= self.min_variant_quality:
                    variants.append(var)
            
            # Build position-based cache
            self._variant_cache = {}
            for var in variants:
                pos = var.POS
                if pos not in self._variant_cache:
                    self._variant_cache[pos] = []
                self._variant_cache[pos].append(var)
                
            # logger.debug(f"Preloaded {len(variants)} variants for region {cache_start}-{cache_end}")
            
        except Exception as e:
            logger.warning(f"Failed to preload variants: {e}")
            self._variant_cache = {}
    
    def _preload_features(self) -> None:
        """Preload features for the region."""
        try:
            # Calculate region to preload
            cache_start = self.current_pos
            cache_end = min(self.current_pos + self._cache_window,
                          self.end) if self.end else self.current_pos + self._cache_window
            
            # Query features from gene map
            features = list(self.gm.gene_map.fetch(
                self.chrom, cache_start, cache_end, 
                features=tuple(self.feature_types) if self.feature_types else tuple()
            ))
            
            # Build position-based cache (features that overlap each position)
            self._feature_cache = {}
            for pos in range(cache_start, cache_end + 1):
                overlapping = [f for f in features 
                             if f.get('start', 0) <= pos <= f.get('end', 0)]
                if overlapping:
                    self._feature_cache[pos] = overlapping
                    
            # logger.debug(f"Preloaded {len(features)} features for region {cache_start}-{cache_end}")
            
        except Exception as e:
            logger.warning(f"Failed to preload features: {e}")
            self._feature_cache = {}
    
    def get_features_at_position(self, position: int, 
                                 exclude_types: Optional[List[str]] = None) -> List[Any]:
        """Get features at a specific genomic position.
        
        Args:
            position: Genomic position (1-based)
            exclude_types: List of feature types to exclude
            
        Returns:
            List of features at the position
        """
        # Use cache if available
        if position in self._feature_cache:
            features = self._feature_cache[position]
            feats = ', '.join(features.keys())
            print(f'using cached features: {feats}')
        else:
            # Direct query if not in cache
            features = list(self.gm.gene_map.get_feature(self.chrom, position))
        
        # Filter out excluded types if specified
        if exclude_types:
            features = [f for f in features 
                       if f.get('feature') not in exclude_types]
        
        return features
    
    def _get_sequence_window(self, start: int, end: int) -> Tuple[str, str]:
        """Get reference and personal sequence for a window.
        
        Returns:
            Tuple of (reference_sequence, personal_sequence)
        """
        # Check cache
        cache_key = (start, end)
        if cache_key in self._sequence_cache:
            return self._sequence_cache[cache_key]
        
        # Get reference sequence
        ref_seq = self.gm.get_sequence(self.chrom, start, end)
        if not ref_seq:
            return "", ""
        
        # Get personal sequence if variants are being integrated
        if self.integrate_variants:
            # Apply variants to the reference sequence
            personal_seq = self._apply_variants_to_window(ref_seq, start, end)
        else:
            personal_seq = ref_seq
        
        # Cache the result
        self._sequence_cache[cache_key] = (ref_seq, personal_seq)
        
        return ref_seq, personal_seq
    
    def _apply_variants_to_window(self, ref_seq: str, start: int, end: int) -> str:
        """Apply variants to a reference sequence window.
        
        Args:
            ref_seq: Reference sequence
            start: Genomic start position
            end: Genomic end position
            
        Returns:
            Personalized sequence with variants applied
        """
        # Convert to list for easier manipulation
        seq_list = list(ref_seq)
        
        # Get variants in this window
        variants_in_window = []
        for pos in range(start, end + 1):
            if pos in self._variant_cache:
                variants_in_window.extend(self._variant_cache[pos])
        
        # Sort by position
        variants_in_window.sort(key=lambda v: v.POS)
        
        # Track position offset due to indels
        offset = 0
        
        for variant in variants_in_window:
            # Skip if variant doesn't have alt allele
            if not variant.ALT:
                continue
                
            # Get the first alt allele (assuming diploid, homozygous or we take first)
            alt_allele = variant.ALT[0]
            ref_allele = variant.REF
            
            # Calculate position in sequence (0-based)
            seq_pos = variant.POS - start + offset
            
            # Validate position
            if seq_pos < 0 or seq_pos >= len(seq_list):
                continue
            
            # Apply variant
            try:
                if len(ref_allele) == len(alt_allele):
                    # SNP - simple substitution
                    for i in range(len(ref_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                            
                elif len(ref_allele) > len(alt_allele):
                    # Deletion
                    del_length = len(ref_allele) - len(alt_allele)
                    # Replace ref with alt
                    for i in range(len(alt_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                    # Delete extra bases
                    for _ in range(del_length):
                        if seq_pos + len(alt_allele) < len(seq_list):
                            del seq_list[seq_pos + len(alt_allele)]
                    offset -= del_length
                    
                else:
                    # Insertion
                    ins_length = len(alt_allele) - len(ref_allele)
                    # Replace ref bases
                    for i in range(len(ref_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                    # Insert additional bases
                    for i in range(len(ref_allele), len(alt_allele)):
                        seq_list.insert(seq_pos + i, alt_allele[i])
                    offset += ins_length
                    
            except Exception as e:
                logger.debug(f"Failed to apply variant at {variant.POS}: {e}")
                
        return ''.join(seq_list)
    
    def __iter__(self) -> Iterator[Union[GenomicPosition, Tuple[str, str, List[Any]]]]:
        """Iterate over genomic positions."""
        return self
    
    def __next__(self) -> Union[GenomicPosition, Tuple[str, str, List[Any]]]:
        """Get the next genomic position or window."""
        # Check if we've reached the end
        if self.end and self.current_pos > self.end:
            raise StopIteration
        
        # Check if we need to reload caches
        if self.current_pos % self._cache_window == 0:
            if self.integrate_variants:
                self._preload_variants()
            if self.track_features:
                self._preload_features()
        
        # Get sequence window
        window_end = self.current_pos + self.window_size - 1
 
        ref_seq, personal_seq = self._get_sequence_window(self.current_pos, window_end)

        
        if not ref_seq:
            # Can't get sequence, advance and try again
            self.current_pos += self.stride
            if self.end and self.current_pos > self.end:
                raise StopIteration
            return self.__next__()
        
        # Get all features that overlap with the window
        features = []
        if self.track_features:
            if self.window_size == 1:
                features = self.get_features_at_position(self.current_pos)
            else:
                # For windows, get all features that overlap with any part of the window
                window_features = list(self.gm.gene_map.fetch(
                    self.chrom, self.current_pos, window_end,
                    features=tuple(self.feature_types) if self.feature_types else tuple()
                ))
                features = window_features
        
        # Get variants in the window
        variants = []
        quality_scores = []
        if self.integrate_variants:
            if self.window_size == 1:
                variants = self._variant_cache.get(self.current_pos, [])
                quality_scores = [v.QUAL for v in variants]
            else:
                # For windows, collect all variants and add them as features
                for pos in range(self.current_pos, min(window_end + 1, self.current_pos + self.window_size)):
                    if pos in self._variant_cache:
                        for var in self._variant_cache[pos]:
                            variants.append(var)
                            # Add variant as a feature for display
                            var_end = var.POS + len(var.REF) - 1
                            features.append({
                                'feature': 'variant',
                                'start': var.POS,
                                'end': var_end,
                                'ref': var.REF,
                                'alt': var.ALT[0] if var.ALT else None,
                                'qual': var.QUAL,
                                'info': {'ref': var.REF, 'alt': var.ALT[0] if var.ALT else None}
                            })
                quality_scores = [v.QUAL for v in variants]
        
        # Prepare result
        if self.window_size == 1:
            # Single base mode
            result = GenomicPosition(
                chrom=self.chrom,
                position=self.current_pos,
                reference_base=ref_seq[0] if ref_seq else 'N',
                personal_base=personal_seq[0] if personal_seq else 'N',
                features=features,
                variants=variants,
                quality_scores=quality_scores
            )
        else:
            # Window mode - return tuple of (ref_seq, personal_seq, features)
            result = (ref_seq, personal_seq, features)
        
        # Advance position
        self.current_pos += self.stride
        
        return result
    
    def jump_to(self, position: int) -> None:
        """Jump to a specific genomic position.
        
        Args:
            position: Genomic position to jump to (1-based)
        """
        self.current_pos = position
        
        # Reload caches for new position
        if self.integrate_variants:
            self._preload_variants()
        if self.track_features:
            self._preload_features()
        
        # Clear sequence cache
        self._sequence_cache = {}
    
    def scan_for_feature(self, feature_type: str, 
                        max_distance: Optional[int] = None) -> Optional[int]:
        """Scan forward for the next occurrence of a feature type.
        
        Args:
            feature_type: Type of feature to search for
            max_distance: Maximum distance to scan (None = unlimited)
            
        Returns:
            Position of next feature or None if not found
        """
        scan_start = self.current_pos
        scan_end = self.end if self.end else scan_start + (max_distance or 1000000)
        
        if max_distance:
            scan_end = min(scan_end, scan_start + max_distance)
        
        # Query features in the scan region
        features = list(self.gm.gene_map.fetch(
            self.chrom, scan_start, scan_end,
            features=(feature_type,)
        ))
        
        if features:
            # Find the closest feature
            closest = min(features, key=lambda f: f.get('start', float('inf')))
            return closest.get('start')
        
        return None
    
    def extract_feature_sequence(self, 
                                feature_type: str,
                                upstream: int = 0,
                                downstream: int = 0,
                                break_on_types: Optional[List[str]] = None) -> Optional[str]:
        """Extract sequence until the end of a feature of specified type.
        
        Args:
            feature_type: Type of feature to extract
            upstream: Bases to include before feature start
            downstream: Bases to include after feature end
            break_on_types: Stop extraction if these feature types are encountered
            
        Returns:
            Extracted sequence or None if feature not found
        """
        sequences = []
        feature_started = False
        feature_ended = False
        
        for item in self:
            if self.window_size == 1:
                # Single base mode
                pos_data = item
                
                # Check if we're in the target feature
                if feature_type in pos_data.feature_types:
                    feature_started = True
                    sequences.append(pos_data.personal_base)
                elif feature_started:
                    # We've left the feature
                    feature_ended = True
                    
                    # Add downstream if requested
                    if downstream > 0:
                        for _ in range(downstream):
                            sequences.append(pos_data.personal_base)
                            try:
                                pos_data = next(self)
                            except StopIteration:
                                break
                    break
                    
                # Check for break conditions
                if break_on_types:
                    if any(bt in pos_data.feature_types for bt in break_on_types):
                        if feature_started:
                            break
            else:
                # Window mode
                ref_seq, personal_seq, features = item
                
                # Check features in window
                has_target = any(f.get('feature') == feature_type for f in features)
                
                if has_target:
                    feature_started = True
                    sequences.append(personal_seq)
                elif feature_started:
                    feature_ended = True
                    break
        
        if sequences:
            return ''.join(sequences)
        return None
    
    def collect_variants_in_window(self, window_size: int) -> List[Any]:
        """Collect all variants within a window from current position.
        
        Args:
            window_size: Size of window to examine
            
        Returns:
            List of variants in the window
        """
        window_start = self.current_pos
        window_end = self.current_pos + window_size - 1
        
        variants = []
        for pos in range(window_start, window_end + 1):
            if pos in self._variant_cache:
                variants.extend(self._variant_cache[pos])
        
        return variants
    
    def get_gc_content_profile(self, window_size: int = 100, 
                              num_windows: int = 10) -> List[float]:
        """Calculate GC content in sliding windows.
        
        Args:
            window_size: Size of each window
            num_windows: Number of windows to calculate
            
        Returns:
            List of GC content values for each window
        """
        gc_contents = []
        
        for _ in range(num_windows):
            window_end = self.current_pos + window_size - 1
            ref_seq, personal_seq = self._get_sequence_window(self.current_pos, window_end)
            
            if personal_seq:
                gc_count = personal_seq.upper().count('G') + personal_seq.upper().count('C')
                gc_content = gc_count / len(personal_seq) if len(personal_seq) > 0 else 0
                gc_contents.append(gc_content)
            
            self.current_pos += self.stride
            
            if self.end and self.current_pos > self.end:
                break
        
        return gc_contents


class FeatureExtractor:
    """Extract sequences for ML analysis with feature context."""
    
    def __init__(self, genome_manager: 'GenomeManager'):
        """Initialize the feature extractor.
        
        Args:
            genome_manager: GenomeManager instance
        """
        self.gm = genome_manager
    
    def extract_around_features(self, 
                               feature_list: List[Any],
                               upstream: int = 1000,
                               downstream: int = 1000,
                               integrate_variants: bool = True) -> List[Dict[str, Any]]:
        """Extract sequences around a list of features.
        
        Args:
            feature_list: List of features to extract sequences around
            upstream: Bases to include upstream
            downstream: Bases to include downstream
            integrate_variants: Whether to integrate personal variants
            
        Returns:
            List of dictionaries with sequences and metadata
        """
        results = []
        
        for feature in feature_list:
            if not hasattr(feature, 'chrom') or not hasattr(feature, 'start'):
                continue
            
            # Create iterator for this region
            iterator = GenomeIterator(
                self.gm,
                feature.chrom,
                max(1, feature.start - upstream),
                feature.end + downstream,
                stride=1,
                integrate_variants=integrate_variants,
                window_size=feature.end - feature.start + upstream + downstream + 1
            )
            
            # Get the sequence
            try:
                ref_seq, personal_seq, features = next(iterator)
                
                result = {
                    'feature_id': getattr(feature, 'sfid', f"{feature.type}_{feature.start}"),
                    'feature_type': feature.type,
                    'chrom': feature.chrom,
                    'start': feature.start - upstream,
                    'end': feature.end + downstream,
                    'reference_sequence': ref_seq,
                    'personal_sequence': personal_seq,
                    'upstream_size': upstream,
                    'downstream_size': downstream,
                    'has_variants': ref_seq != personal_seq,
                    'feature_start_in_seq': upstream,
                    'feature_end_in_seq': upstream + (feature.end - feature.start)
                }
                
                # Add variant information
                if integrate_variants:
                    variants = iterator.collect_variants_in_window(len(personal_seq))
                    result['num_variants'] = len(variants)
                    result['variant_positions'] = [v.POS - feature.start + upstream for v in variants]
                
                # Calculate GC content
                if personal_seq:
                    gc_count = personal_seq.upper().count('G') + personal_seq.upper().count('C')
                    result['gc_content'] = gc_count / len(personal_seq)
                
                results.append(result)
                
            except StopIteration:
                logger.warning(f"Could not extract sequence for feature {feature}")
                continue
        
        return results
    
    def prepare_for_ml(self, 
                       sequences: List[Dict[str, Any]],
                       sequence_length: int = 1000,
                       one_hot_encode: bool = True) -> Dict[str, Any]:
        """Prepare sequences for ML model input.
        
        Args:
            sequences: List of sequence dictionaries from extract_around_features
            sequence_length: Fixed length for ML input (pad or truncate)
            one_hot_encode: Whether to one-hot encode sequences
            
        Returns:
            Dictionary with prepared data for ML
        """
        import numpy as np
        
        # DNA to index mapping
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        
        prepared_sequences = []
        metadata = []
        
        for seq_dict in sequences:
            seq = seq_dict['personal_sequence']
            
            # Pad or truncate to fixed length
            if len(seq) > sequence_length:
                # Take center portion
                start = (len(seq) - sequence_length) // 2
                seq = seq[start:start + sequence_length]
            elif len(seq) < sequence_length:
                # Pad with N's
                pad_left = (sequence_length - len(seq)) // 2
                pad_right = sequence_length - len(seq) - pad_left
                seq = 'N' * pad_left + seq + 'N' * pad_right
            
            if one_hot_encode:
                # Convert to one-hot encoding
                encoded = np.zeros((sequence_length, 5))
                for i, base in enumerate(seq.upper()):
                    if base in base_to_idx:
                        encoded[i, base_to_idx[base]] = 1
                    else:
                        encoded[i, 4] = 1  # Unknown base
                prepared_sequences.append(encoded)
            else:
                # Convert to indices
                indices = [base_to_idx.get(b.upper(), 4) for b in seq]
                prepared_sequences.append(indices)
            
            # Store metadata
            meta = {
                'feature_id': seq_dict['feature_id'],
                'feature_type': seq_dict['feature_type'],
                'chrom': seq_dict['chrom'],
                'original_length': len(seq_dict['personal_sequence']),
                'has_variants': seq_dict['has_variants'],
                'gc_content': seq_dict.get('gc_content', 0)
            }
            metadata.append(meta)
        
        return {
            'sequences': np.array(prepared_sequences),
            'metadata': metadata,
            'sequence_length': sequence_length,
            'encoding': 'one_hot' if one_hot_encode else 'index'
        }