"""Metric-based motif detection (GC content, etc)."""

import numpy as np
from typing import List
from .motif import BaseMotif, MatchResult


class MetricMotif(BaseMotif):
    """Compute sequence metrics (GC content, etc).

    Unlike pattern motifs, metric motifs compute a value for a window
    rather than finding discrete instances.
    """

    def __init__(self, name: str, scoring_function, window_size: int = 1):
        super().__init__(name, motif_id=name)
        self.score_func = scoring_function
        self._window_size = window_size

    @property
    def length(self) -> int:
        return self._window_size

    def scan(self, seq: str) -> np.ndarray:
        """Compute metric at each position."""
        if self._window_size == 1:
            # Per-base metric
            return np.array([self.score_func(seq)])
        else:
            # Sliding window
            scores = []
            for i in range(len(seq) - self._window_size + 1):
                window = seq[i:i + self._window_size]
                scores.append(self.score_func(window))
            return np.array(scores) if scores else np.array([])

    def find_instances(self, seq: str, threshold: float = 0.0,
                      chrom: str = "", offset: int = 0) -> List[MatchResult]:
        """Metric motifs don't have discrete instances."""
        return []

    def __call__(self, seq: str) -> float:
        return self.score_func(seq)


def gc_content(seq: str) -> float:
    """Compute GC content of sequence."""
    if not seq:
        return 0.0
    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq)


gcpattern = MetricMotif('gc_content', gc_content)
