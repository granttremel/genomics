
from ggene import aliases_dna, aliases_rna, to_rna, reverse_complement, is_rna, is_dna
from abc import ABC, abstractmethod
import numpy as np
import re

from .motifio import MotifIO

def calculate_significance(self, seq, motif_score):
      """Calculate p-value for motif occurrence"""
      # Shuffle sequence to create null distribution
      # Or use pre-computed background model
      pass

class MotifDetector:

    def __init__(self):
        self.motifs={}
        self.io = MotifIO()
    
    def add_motif(self, motif):
        self.motifs[motif.name]=motif
    
    def identify(self, seq, motifs=[]):
        if not motifs:
            motifs=self.motifs.keys()
            
        for m in motifs:
            pass
    
    def score(self, seq, motifs=[]):
        pass

class IndexedMotifDetector:
    """Pre-index genome for fast motif searches"""

    def __init__(self):
        self.suffix_array = None
        self.bwt_index = None

    def build_index(self, genome_sequence):
        """Build suffix array or BWT for O(m) searches"""
        pass

class BaseMotif(ABC):
    
    def __init__(self, name):
        self.name=name
        
    @abstractmethod
    def score(self, seq, return_positions=False):
        
        pass
    
    @abstractmethod
    def find_instances(self, seq, threshold=None):
        """Return list of (start, end, score) tuples"""
        pass
    
    def to_features(self, seq, window_size=None):
        """Convert to ML-ready features"""
        if window_size:
            return self._windowed_features(seq, window_size)
        return self._sequence_features(seq)



