"""Streaming motif detection with scale-aware processing."""

from abc import ABC, abstractmethod
from collections import deque
import numpy as np
from typing import Iterator, Tuple, List, Optional

from motifs.repeat import RepeatMotif

class StreamingMotif(ABC):
    """Base class for streaming motif detection."""
    
    def __init__(self, name: str, window_size: int, stride: int = 1):
        self.name = name
        self.window_size = window_size
        self.stride = stride
        self.buffer = deque(maxlen=window_size)
        self.position = 0
        
    @abstractmethod
    def process_window(self, seq: str) -> Optional[Tuple[int, int, float]]:
        """Process current window, return (start, end, score) if motif found."""
        pass
    
    def stream_sequence(self, seq_iterator: Iterator[str]) -> Iterator[Tuple[int, int, float]]:
        """Stream through sequence, yielding motif occurrences."""
        for base in seq_iterator:
            self.buffer.append(base)
            self.position += 1
            
            # Only process when we have a full window and at stride intervals
            if len(self.buffer) == self.window_size and self.position % self.stride == 0:
                window_seq = ''.join(self.buffer)
                result = self.process_window(window_seq)
                if result:
                    # Adjust position to account for buffer
                    start, end, score = result
                    yield (self.position - self.window_size + start, 
                          self.position - self.window_size + end, 
                          score)


class MultiScaleMotifDetector:
    """Manages motifs at different scales efficiently."""
    
    def __init__(self):
        self.fast_motifs = []  # Check every base
        self.medium_motifs = []  # Check every 10 bases
        self.slow_motifs = []  # Check every 100 bases
        self.chunk_motifs = []  # Check on full chunks
        
    def add_motif(self, motif, scale='auto'):
        """Add motif at appropriate scale."""
        if scale == 'auto':
            scale = self._determine_scale(motif)
            
        if scale == 'fast':
            self.fast_motifs.append(motif)
        elif scale == 'medium':
            self.medium_motifs.append(motif)
        elif scale == 'slow':
            self.slow_motifs.append(motif)
        else:
            self.chunk_motifs.append(motif)
    
    def _determine_scale(self, motif):
        """Automatically determine processing scale based on motif properties."""
        if hasattr(motif, 'pattern') and len(motif.pattern) < 10:
            return 'fast'  # Short patterns need frequent checking
        elif hasattr(motif, 'window_size') and motif.window_size > 1000:
            return 'chunk'  # Large-scale features
        elif isinstance(motif, RepeatMotif):
            return 'medium'  # Repeats can skip some positions
        else:
            return 'medium'
    
    def process_stream(self, seq_iterator, chunk_size=10000):
        """Process sequence stream at multiple scales."""
        position = 0
        chunk_buffer = []
        
        for base in seq_iterator:
            position += 1
            chunk_buffer.append(base)
            
            # Fast motifs - every base
            if self.fast_motifs:
                for motif in self.fast_motifs:
                    # Process with small window
                    pass
            
            # Medium motifs - every 10 bases
            if position % 10 == 0 and self.medium_motifs:
                for motif in self.medium_motifs:
                    # Process with medium window
                    pass
            
            # Slow motifs - every 100 bases
            if position % 100 == 0 and self.slow_motifs:
                for motif in self.slow_motifs:
                    # Process with large window
                    pass
            
            # Chunk motifs - when buffer is full
            if len(chunk_buffer) >= chunk_size:
                chunk_seq = ''.join(chunk_buffer)
                for motif in self.chunk_motifs:
                    # Process full chunk
                    pass
                chunk_buffer = chunk_buffer[-1000:]  # Keep overlap


class AdaptiveMotifScanner:
    """Dynamically adjusts scanning based on sequence content."""
    
    def __init__(self):
        self.motifs = []
        self.active_regions = {}  # Track regions where motifs are dense
        
    def scan_with_adaptation(self, seq):
        """Scan more carefully in interesting regions."""
        # Quick initial scan
        hotspots = self._find_hotspots(seq)
        
        # Detailed scan in hotspots
        for start, end in hotspots:
            self._detailed_scan(seq[start:end], offset=start)
    
    def _find_hotspots(self, seq, window=1000):
        """Identify regions likely to contain motifs."""
        hotspots = []
        
        # Look for characteristics like:
        # - High GC content (regulatory regions)
        # - Low complexity (repeats)
        # - Conserved patterns
        
        for i in range(0, len(seq), window):
            window_seq = seq[i:i+window]
            if self._is_interesting(window_seq):
                hotspots.append((i, min(i+window, len(seq))))
        
        return hotspots