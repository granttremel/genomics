"""Position Weight Matrix implementation for motif detection.

This module provides PWM-based motif detection with:
- Log-odds scoring for efficient computation
- Vectorized batch scanning across multiple PWMs
- Integration with MotifLibrary interface for streaming

Classes:
    PWM: Single position weight matrix (extends BaseMotif)
    PWMLibrary: Collection of PWMs with vectorized batch scanning (extends MotifLibrary)
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Iterator

from ggene.seqs import bio
from ggene.motifs.motif import BaseMotif, MotifLibrary, MatchResult


class PWM(BaseMotif):
    """Position Weight Matrix for sequence motifs.

    PWMs capture position-specific nucleotide preferences learned from
    binding site alignments. This implementation uses log-odds scoring
    against background frequencies for efficient and statistically
    meaningful scoring.

    Attributes:
        counts: Raw count matrix (4 x length), rows are [A, C, G, T]
        probs: Probability matrix (counts normalized per position)
        log_odds: Log-odds scoring matrix against background
        ic: Information content per position (bits)
        consensus: Consensus sequence using IUPAC codes
        max_score: Maximum possible score
        min_score: Minimum possible score
    """

    bases = ['A', 'C', 'G', 'T']
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    def __init__(self, name: str, counts: np.ndarray,
                 pseudocount: float = 0.1,
                 background: Dict[str, float] = None,
                 motif_id: str = ""):
        """Initialize PWM from count matrix.

        Args:
            name: Human-readable motif name
            counts: 4xN matrix where rows are [A,C,G,T] and columns are positions
            pseudocount: Small value to avoid log(0)
            background: Background nucleotide frequencies (default: uniform 0.25)
            motif_id: Unique identifier (defaults to name)
        """
        super().__init__(name, motif_id=motif_id or name)

        # Store raw counts with pseudocounts
        self._counts = counts.astype(float) + pseudocount
        self._length = counts.shape[1]

        # Calculate probability matrix (normalized per position)
        self.probs = self._counts / self._counts.sum(axis=0)

        # Background frequencies (default: uniform)
        if background is None:
            background = {b: 0.25 for b in self.bases}
        self.background = background

        # Calculate log-odds matrix for efficient scoring
        self.log_odds = self._calculate_log_odds()

        # Information content per position (for visualization)
        self.ic = self._calculate_information_content()

        # Score bounds for threshold normalization
        self.max_score = float(np.sum(self.log_odds.max(axis=0)))
        self.min_score = float(np.sum(self.log_odds.min(axis=0)))

        # Consensus sequence (computed on demand)
        self._consensus = None

    @property
    def length(self) -> int:
        """Motif length (number of positions)."""
        return self._length

    @property
    def counts(self) -> np.ndarray:
        """Raw count matrix."""
        return self._counts

    @property
    def consensus(self) -> str:
        """Consensus sequence using IUPAC ambiguity codes."""
        if self._consensus is None:
            self._consensus = self._compute_consensus()
        return self._consensus

    def _calculate_log_odds(self) -> np.ndarray:
        """Calculate log-odds scoring matrix.

        Log-odds = log2(P(base|motif) / P(base|background))

        This transformation:
        - Makes scores additive across positions
        - Accounts for background nucleotide composition
        - Gives positive scores for enriched bases, negative for depleted
        """
        log_odds = np.zeros_like(self.probs)
        for i, base in enumerate(self.bases):
            log_odds[i, :] = np.log2(self.probs[i, :] / self.background[base])
        return log_odds

    def _calculate_information_content(self) -> np.ndarray:
        """Calculate information content at each position (bits).

        IC = 2 - entropy, where max entropy for 4 bases is 2 bits.
        Higher IC means more conserved positions.
        """
        ic = np.zeros(self._length)
        for pos in range(self._length):
            # Shannon entropy
            probs = self.probs[:, pos]
            entropy = -np.sum(probs * np.log2(probs + 1e-10))
            ic[pos] = 2 - entropy
        return ic

    def _compute_consensus(self, min_p_ratio: float = 2.0) -> str:
        """Compute consensus sequence with IUPAC ambiguity codes.

        At each position, bases with probability less than 1/(n*min_p_ratio)
        are eliminated, where n is current number of candidates. The remaining
        bases are combined into an IUPAC code.

        Args:
            min_p_ratio: Threshold ratio for eliminating low-probability bases
        """
        cons = []

        for pos in range(self._length):
            p_bases = dict(zip(self.bases, self.probs[:, pos]))
            n_bases = 4

            while n_bases > 1:
                min_base, min_p = min(p_bases.items(), key=lambda x: x[1])
                threshold = 1 / n_bases / min_p_ratio

                if min_p < threshold:
                    del p_bases[min_base]
                    n_bases -= 1
                else:
                    break

            # Convert remaining bases to IUPAC code
            remaining = list(p_bases.keys())
            cons.append(bio.get_alias(*remaining))

        return "".join(cons)

    def score_sequence(self, seq: str) -> float:
        """Score a sequence using log-odds.

        Args:
            seq: Sequence of exactly self.length bases

        Returns:
            Log-odds score (sum across positions)
        """
        if len(seq) != self._length:
            raise ValueError(f"Sequence length {len(seq)} doesn't match PWM length {self._length}")

        score = 0.0
        for pos, base in enumerate(seq.upper()):
            if base in self.base_to_idx:
                score += self.log_odds[self.base_to_idx[base], pos]
            else:
                # Handle ambiguous bases by averaging
                score += self.log_odds[:, pos].mean()

        return score

    def scan(self, seq: str) -> np.ndarray:
        """Score each position in sequence using log-odds.

        Args:
            seq: DNA sequence to scan

        Returns:
            Array of normalized scores (0-1 range), length = len(seq) - length + 1
        """
        seq = seq.upper()
        n_positions = len(seq) - self._length + 1

        if n_positions <= 0:
            return np.array([])

        scores = np.zeros(n_positions)

        for i in range(n_positions):
            subseq = seq[i:i + self._length]
            raw_score = self.score_sequence(subseq)
            # Normalize to 0-1 range
            scores[i] = (raw_score - self.min_score) / (self.max_score - self.min_score)

        return scores

    def find_instances(self, seq: str, threshold: float = 0.8,
                      chrom: str = "", offset: int = 0) -> List[MatchResult]:
        """Find all matches above threshold.

        Args:
            seq: Sequence to scan
            threshold: Minimum normalized score (0-1 range, default 0.8)
            chrom: Chromosome for genomic coordinates
            offset: Genomic offset for coordinates

        Returns:
            List of MatchResult objects
        """
        scores = self.scan(seq)
        results = []

        for i, score in enumerate(scores):
            if score >= threshold:
                end = i + self._length
                results.append(MatchResult(
                    start=offset + i if offset else i,
                    end=offset + end if offset else end,
                    score=float(score),
                    name=self.name,
                    motif_id=self.id,
                    strand="+",
                    seq=seq[i:end],
                    chrom=chrom,
                ))

        return results

    @property
    def total_ic(self) -> float:
        """Total information content (sum across positions)."""
        return float(np.sum(self.ic))

    @property
    def entropy(self) -> float:
        """Average entropy per position."""
        return float(2 - np.mean(self.ic))

    def __repr__(self) -> str:
        return f"PWM({self.name!r}, len={self._length}, IC={self.total_ic:.2f})"


class PWMLibrary(MotifLibrary):
    """Collection of PWMs with vectorized batch scanning.

    This library provides efficient batch scanning by:
    1. Assembling all PWMs into a 3D tensor (num_pwms x 4 x max_length)
    2. Converting sequences to one-hot encoding
    3. Using vectorized numpy operations for parallel scoring

    The 3D tensor approach provides ~10-100x speedup over individual scans
    when scanning with many PWMs.
    """

    bases = "ACGT"

    def __init__(self):
        self._pwms: Dict[str, PWM] = {}
        self._pwm_ids: List[str] = []

        # Cached assembled tensor (invalidated on add)
        self._assembled: Optional[Tuple[np.ndarray, np.ndarray]] = None

    @property
    def motif_ids(self) -> List[str]:
        """List of PWM IDs in insertion order."""
        return self._pwm_ids.copy()

    @property
    def max_length(self) -> int:
        """Maximum PWM length."""
        if not self._pwms:
            return 0
        return max(pwm.length for pwm in self._pwms.values())

    def add(self, motif: BaseMotif) -> None:
        """Add a PWM to the library."""
        if not isinstance(motif, PWM):
            raise TypeError(f"Expected PWM, got {type(motif)}")
        self._pwms[motif.id] = motif
        if motif.id not in self._pwm_ids:
            self._pwm_ids.append(motif.id)
        self._assembled = None  # Invalidate cache

    def add_pwm(self, pwm: PWM) -> None:
        """Convenience alias for add()."""
        self.add(pwm)

    def get_motif(self, motif_id: str) -> PWM:
        """Get PWM by ID."""
        return self._pwms[motif_id]

    def get_pwm(self, pwm_id: str) -> PWM:
        """Convenience alias for get_motif()."""
        return self.get_motif(pwm_id)

    def _assemble_pwms(self) -> Tuple[np.ndarray, np.ndarray]:
        """Assemble all PWMs into a 3D tensor for vectorized scanning.

        Returns:
            Tuple of:
            - big_pwm: (num_pwms, 4, max_length) tensor of log-odds matrices
            - pwm_lens: (num_pwms,) array of actual lengths
        """
        if self._assembled is not None:
            return self._assembled

        num_pwms = len(self._pwms)
        max_len = self.max_length

        big_pwm = np.zeros((num_pwms, 4, max_len))
        pwm_lens = np.zeros(num_pwms)

        for i, pwm_id in enumerate(self._pwm_ids):
            pwm = self._pwms[pwm_id]
            big_pwm[i, :, :pwm.length] = pwm.log_odds
            pwm_lens[i] = pwm.length

        self._assembled = (big_pwm, pwm_lens)
        return self._assembled

    @staticmethod
    def seq_to_onehot(seq: str, pad_length: int = 0) -> np.ndarray:
        """Convert sequence to one-hot encoding.

        Args:
            seq: DNA sequence
            pad_length: Additional padding at end (for sliding window)

        Returns:
            (4, len(seq) + pad_length) one-hot array, rows are [A, C, G, T]
        """
        base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        oh = np.zeros((4, len(seq) + pad_length))

        for i, b in enumerate(seq.upper()):
            if b in base_map:
                oh[base_map[b], i] = 1

        return oh

    def scan_all(self, seq: str) -> Dict[str, np.ndarray]:
        """Scan sequence with all PWMs using vectorized operations.

        Uses 3D tensor multiplication for efficient batch scoring:
        - Assembles all PWMs into (num_pwms, 4, max_length) tensor
        - One-hot encodes sequence once
        - Scores all PWMs at each position in parallel

        Args:
            seq: Sequence to scan

        Returns:
            Dict mapping pwm_id -> normalized score array
        """
        if not self._pwms:
            return {}

        big_pwm, pwm_lens = self._assemble_pwms()
        pwm_len = big_pwm.shape[2]
        seq_len = len(seq)

        # One-hot encode sequence with padding for sliding window
        seq_oh = self.seq_to_onehot(seq, pad_length=pwm_len)

        # Raw scores array: (num_pwms, seq_len)
        raw_scores = np.zeros((len(self._pwms), seq_len))

        # Slide across sequence positions
        for pos in range(seq_len):
            # Mask for PWMs that fit at this position
            valid_mask = (pwm_lens <= seq_len - pos)

            # Extract subsequence one-hot
            subseq = seq_oh[:, pos:pos + pwm_len]
            subseq_exp = subseq[None, :, :]  # (1, 4, max_len)

            # Vectorized score: element-wise multiply and sum
            scores = np.sum(big_pwm * subseq_exp, axis=(1, 2))

            # Zero out scores for PWMs that don't fit
            raw_scores[:, pos] = scores * valid_mask

        # Normalize scores to 0-1 range per PWM
        results = {}
        for i, pwm_id in enumerate(self._pwm_ids):
            pwm = self._pwms[pwm_id]
            n_valid = seq_len - pwm.length + 1

            if n_valid <= 0:
                results[pwm_id] = np.array([])
                continue

            # Normalize by score range
            raw = raw_scores[i, :n_valid]
            normalized = (raw - pwm.min_score) / (pwm.max_score - pwm.min_score)
            results[pwm_id] = normalized

        return results

    def find_all_instances(self, seq: str, threshold: float = 0.0,
                          chrom: str = "", offset: int = 0) -> Iterator[MatchResult]:
        """Find instances of all PWMs above threshold.

        Args:
            seq: Sequence to scan
            threshold: Minimum normalized score (0-1 range)
            chrom: Chromosome for genomic coordinates
            offset: Genomic offset for coordinates

        Yields:
            MatchResult objects sorted by position
        """
        all_scores = self.scan_all(seq)
        results = []

        for pwm_id, scores in all_scores.items():
            pwm = self._pwms[pwm_id]

            # Find positions above threshold
            hits = np.where(scores >= threshold)[0]

            for i in hits:
                end = i + pwm.length
                results.append(MatchResult(
                    start=offset + i if offset else i,
                    end=offset + end if offset else end,
                    score=float(scores[i]),
                    name=pwm.name,
                    motif_id=pwm_id,
                    strand="+",
                    seq=seq[i:end],
                    chrom=chrom,
                ))

        # Sort by position for heapq.merge compatibility
        results.sort(key=lambda r: (r.start, r.end))
        yield from results

    def __repr__(self) -> str:
        return f"PWMLibrary({len(self._pwms)} PWMs)"
