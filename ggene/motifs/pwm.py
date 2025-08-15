
"""Position Weight Matrix implementation for motif detection."""

import numpy as np
from typing import List, Tuple, Optional, Dict
import math

class PWM:
    """Position Weight Matrix for sequence motifs.
    
    PWMs are powerful because they:
    1. Capture position-specific nucleotide preferences
    2. Account for evolutionary variation
    3. Provide probabilistic scoring
    4. Can be learned from real binding data
    """
    
    def __init__(self, counts: np.ndarray, pseudocount: float = 0.1, 
                 background: Dict[str, float] = None):
        """Initialize PWM from count matrix.
        
        Args:
            counts: 4xN matrix where rows are [A,C,G,T] and columns are positions
            pseudocount: Small value to avoid log(0)
            background: Background nucleotide frequencies
        """
        self.bases = ['A', 'C', 'G', 'T']
        self.base_to_idx = {b: i for i, b in enumerate(self.bases)}
        
        # Add pseudocounts to avoid zeros
        self.counts = counts + pseudocount
        self.length = counts.shape[1]
        
        # Calculate probability matrix
        self.probs = self.counts / self.counts.sum(axis=0)
        
        # Background frequencies (default: uniform)
        if background is None:
            background = {b: 0.25 for b in self.bases}
        self.background = background
        
        # Calculate log-odds matrix for efficient scoring
        self.log_odds = self._calculate_log_odds()
        
        # Information content per position (for visualization)
        self.ic = self._calculate_information_content()
        
        # Consensus sequence
        self.consensus = self._get_consensus()
        
        # Score thresholds
        self.max_score = self._calculate_max_score()
        self.min_score = self._calculate_min_score()
    
    def _calculate_log_odds(self) -> np.ndarray:
        """Calculate log-odds scoring matrix."""
        log_odds = np.zeros_like(self.probs)
        for i, base in enumerate(self.bases):
            log_odds[i, :] = np.log2(self.probs[i, :] / self.background[base])
        return log_odds
    
    def _calculate_information_content(self) -> np.ndarray:
        """Calculate information content at each position."""
        ic = np.zeros(self.length)
        for pos in range(self.length):
            # Shannon entropy
            entropy = -sum(p * np.log2(p) if p > 0 else 0 
                          for p in self.probs[:, pos])
            # Information content = max_entropy - entropy
            ic[pos] = 2 - entropy  # 2 bits is max for 4 bases
        return ic
    
    def _get_consensus(self) -> str:
        """Get consensus sequence using IUPAC codes."""
        iupac = {
            'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
            'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 
            'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H',
            'AGT': 'D', 'CGT': 'B', 'ACGT': 'N'
        }
        
        consensus = []
        for pos in range(self.length):
            # Get bases with >25% frequency
            significant = [self.bases[i] for i in range(4) 
                          if self.probs[i, pos] > 0.25]
            key = ''.join(sorted(significant)) if significant else 'N'
            consensus.append(iupac.get(key, 'N'))
        
        return ''.join(consensus)
    
    def _calculate_max_score(self) -> float:
        """Calculate maximum possible score."""
        return sum(self.log_odds.max(axis=0))
    
    def _calculate_min_score(self) -> float:
        """Calculate minimum possible score."""
        return sum(self.log_odds.min(axis=0))
    
    def score_sequence(self, seq: str) -> float:
        """Score a sequence using log-odds."""
        if len(seq) != self.length:
            raise ValueError(f"Sequence length {len(seq)} doesn't match PWM length {self.length}")
        
        score = 0.0
        for pos, base in enumerate(seq.upper()):
            if base in self.base_to_idx:
                score += self.log_odds[self.base_to_idx[base], pos]
            else:
                # Handle ambiguous bases by averaging
                score += self.log_odds[:, pos].mean()
        
        return score
    
    def scan_sequence(self, seq: str, threshold: float = 0.8) -> List[Tuple[int, float, str]]:
        """Scan sequence for motif occurrences.
        
        Args:
            seq: Sequence to scan
            threshold: Score threshold as fraction of max score
            
        Returns:
            List of (position, score, matched_sequence) tuples
        """
        min_score = self.min_score + threshold * (self.max_score - self.min_score)
        hits = []
        
        for i in range(len(seq) - self.length + 1):
            subseq = seq[i:i + self.length]
            score = self.score_sequence(subseq)
            
            if score >= min_score:
                # Calculate p-value (simplified)
                relative_score = (score - self.min_score) / (self.max_score - self.min_score)
                hits.append((i, score, subseq, relative_score))
        
        return hits