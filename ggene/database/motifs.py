"""Motif detection streams that integrate with UGenomeAnnotations.

These streams scan sequences for motifs (PWMs, regex patterns, HMMs) and yield
UFeature objects that can be merged with annotation streams via heapq.merge.

Key design:
- MotifStream is the base class with chunked scanning
- Sequence data is provided via a callable injected by UGenomeAnnotations
- Features are buffered and yielded in sorted order for heapq.merge compatibility
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Iterator, List, Tuple, Dict, Any, Callable

from ggene.database.ufeature import UFeature


class MotifStream(ABC):
    """Base class for motif streams that scan sequences in chunks.

    Subclasses implement `scan_sequence()` to detect motifs in a sequence chunk.
    The base class handles chunking, caching, and sorted iteration.

    Attributes:
        chunk_size: Size of sequence chunks to scan (default 2048)
        overlap: Overlap between chunks to catch boundary-spanning motifs
        feature_type: Type string for generated features
        source_name: Source string for generated features
    """
    chunk_size: int = 2048
    overlap: int = 64  # Must be >= max motif length
    feature_type: str = "motif"
    source_name: str = "Motif"

    def __init__(self):
        self._get_sequence: Optional[Callable[[str, int, int], str]] = None
        self._loaded = False

    def bind_sequence_getter(self, getter: Callable[[str, int, int], str]):
        """Called by UGenomeAnnotations when registering this stream.

        Args:
            getter: Callable(chrom, start, end) -> sequence string
        """
        self._get_sequence = getter

    @abstractmethod
    def load(self):
        """Load motif definitions (PWMs, patterns, etc). Called once before scanning."""
        pass

    @abstractmethod
    def scan_sequence(self, seq: str, chrom: str, offset: int, strand: str = "+") -> List[UFeature]:
        """Scan a sequence chunk and return detected motif features.

        Args:
            seq: The sequence string to scan
            chrom: Chromosome name (for feature coordinates)
            offset: Genomic position of seq[0] (1-based)
            strand: Strand being scanned ('+' or '-')

        Returns:
            List of UFeature objects with genomic coordinates (offset-adjusted)
        """
        pass

    def stream(self, chrom: Optional[str] = None,
               start: Optional[int] = None,
               end: Optional[int] = None) -> Iterator[UFeature]:
        """Stream motif features, scanning sequences in chunks as needed.

        This generator scans lazily - chunks are processed as features are consumed.
        Features are yielded in sorted genomic order for heapq.merge compatibility.
        """
        if not self._get_sequence:
            raise RuntimeError("MotifStream not bound to sequence getter. "
                              "Register with UGenomeAnnotations first.")

        if not self._loaded:
            self.load()
            self._loaded = True

        if not chrom or start is None or end is None:
            return

        # Generator state
        scanned_up_to = start
        buffer: List[UFeature] = []

        def scan_next_chunk():
            nonlocal scanned_up_to, buffer

            chunk_start = scanned_up_to
            chunk_end = min(chunk_start + self.chunk_size, end + self.overlap)

            if chunk_start >= end:
                return False  # Nothing more to scan

            try:
                seq = self._get_sequence(chrom, chunk_start, chunk_end)
            except Exception:
                return False

            if not seq:
                return False

            # Scan forward strand
            features = self.scan_sequence(seq, chrom, chunk_start, "+")

            # TODO: Optionally scan reverse complement for reverse strand motifs
            # rev_seq = reverse_complement(seq)
            # features.extend(self.scan_sequence(rev_seq, chrom, chunk_start, "-"))

            # Filter to features that start within the valid scan range
            # (avoid duplicates from overlapping chunks)
            scan_end = chunk_end - self.overlap if chunk_end < end + self.overlap else chunk_end
            for f in features:
                if chunk_start <= f.start < scan_end and f.end <= end:
                    buffer.append(f)

            buffer.sort()  # Keep sorted for heapq.merge

            # Advance, with overlap to catch boundary-spanning motifs
            scanned_up_to = chunk_end - self.overlap
            return True

        # Initial scan
        scan_next_chunk()

        while True:
            # If buffer empty and more to scan, scan ahead
            while not buffer and scanned_up_to < end:
                if not scan_next_chunk():
                    break

            if not buffer:
                return  # Done

            # Yield next feature
            yield buffer.pop(0)

            # Pre-scan if buffer getting low (reduces latency)
            if len(buffer) < 5 and scanned_up_to < end:
                scan_next_chunk()


class JasparStream(MotifStream):
    """Stream for JASPAR PWM motif detection."""

    feature_type = "tf_binding"
    source_name = "JASPAR"

    def __init__(self, filepath: str, min_score: float = 0.85):
        super().__init__()
        self.filepath = Path(filepath).absolute()
        self.min_score = min_score
        self.motifs = None  # JasparLibrary

    def load(self):
        """Load JASPAR PWM library from directory."""
        from ggene.motifs import jaspar
        self.motifs = jaspar.Jaspar.from_path(self.filepath)

    def scan_sequence(self, seq: str, chrom: str, offset: int, strand: str = "+") -> List[UFeature]:
        """Scan sequence with all loaded PWMs."""
        if not self.motifs:
            return []

        features = []
        instances = self.motifs.find_instances(seq, min_score=self.min_score)

        for (mtf_id, mtf_name), inst_data in instances.items():
            for positions in inst_data:
                # inst_data is list of (start_indices,) tuples from numpy nonzero
                if isinstance(positions, tuple):
                    # Handle numpy nonzero output
                    for rel_start in positions:
                        mtf_len = self.motifs.jaspars[mtf_id].num_pos
                        features.append(UFeature({
                            'chrom': chrom,
                            'start': offset + int(rel_start),
                            'end': offset + int(rel_start) + mtf_len - 1,
                            'feature_type': self.feature_type,
                            'source': self.source_name,
                            'name': mtf_name,
                            'id': mtf_id,
                            'strand': strand,
                            'score': self.min_score,  # Could track actual score
                            'motif_id': mtf_id,
                        }))
                else:
                    # Handle (start, end, score) tuple format
                    rel_start, rel_end, score = positions
                    features.append(UFeature({
                        'chrom': chrom,
                        'start': offset + int(rel_start),
                        'end': offset + int(rel_end),
                        'feature_type': self.feature_type,
                        'source': self.source_name,
                        'name': mtf_name,
                        'id': mtf_id,
                        'strand': strand,
                        'score': float(score),
                        'motif_id': mtf_id,
                    }))

        return features


class PatternStream(MotifStream):
    """Stream for regex pattern motif detection."""

    feature_type = "pattern"
    source_name = "Pattern"

    def __init__(self, patterns: Optional[Dict[str, str]] = None):
        """Initialize with pattern dict: {name: regex_pattern}."""
        super().__init__()
        self._pattern_defs = patterns or {}
        self._compiled = {}

    def add_pattern(self, name: str, pattern: str):
        """Add a regex pattern."""
        self._pattern_defs[name] = pattern

    def load(self):
        """Compile regex patterns."""
        import re
        for name, pattern in self._pattern_defs.items():
            self._compiled[name] = re.compile(pattern, re.IGNORECASE)

    def scan_sequence(self, seq: str, chrom: str, offset: int, strand: str = "+") -> List[UFeature]:
        """Find all pattern matches in sequence."""
        features = []

        for name, regex in self._compiled.items():
            for match in regex.finditer(seq):
                features.append(UFeature({
                    'chrom': chrom,
                    'start': offset + match.start(),
                    'end': offset + match.end() - 1,
                    'feature_type': self.feature_type,
                    'source': self.source_name,
                    'name': name,
                    'strand': strand,
                    'matched_seq': match.group(),
                }))

        return features


class HMMStream(MotifStream):
    """Stream for HMM-based motif detection (placeholder)."""

    feature_type = "hmm_hit"
    source_name = "HMM"

    def __init__(self, model_path: Optional[str] = None):
        super().__init__()
        self.model_path = model_path

    def load(self):
        """Load HMM model."""
        # TODO: Implement HMM loading (pyhmmer, etc.)
        pass

    def scan_sequence(self, seq: str, chrom: str, offset: int, strand: str = "+") -> List[UFeature]:
        """Scan with HMM model."""
        # TODO: Implement HMM scanning
        return []
