"""
Sequence processing utilities for genome browser.

This module handles sequence transformations, variant detection, codon
identification, and translation.
"""

from typing import List, Tuple, Dict, Optional, Set
import re
# from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, to_rna, reverse_complement

from ggene.seqs.vocab import infer_vocab
from ggene.seqs import bio
from ggene.seqs.bio import CODON_TABLE, COMPLEMENT_MAP, to_rna, reverse_complement

class SequenceProcessor:
    """Handles sequence transformations and analysis."""

    def __init__(self):
        """Initialize the sequence processor."""
        self.start_codons = {'ATG', 'AUG'}
        self.stop_codons = {'TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA'}

    def detect_variants(self, ref_seq: str, alt_seq: str) -> List[int]:
        """Identify variant positions between two aligned sequences.

        Args:
            ref_seq: Reference sequence
            alt_seq: Alternative/personal sequence

        Returns:
            List of positions (0-based) where sequences differ
        """
        variants = []
        for i, (ref_base, alt_base) in enumerate(zip(ref_seq, alt_seq)):
            if ref_base != alt_base and alt_base != '-' and ref_base != '-':
                # Skip N bases and gaps
                if ref_base != 'N' and alt_base != 'N':
                    variants.append(i)
        return variants

    def find_codons(self, seq: str, frame: int = 0) -> Dict[str, List[int]]:
        """Find start and stop codons in a sequence.

        Args:
            seq: DNA or RNA sequence
            frame: Reading frame (0, 1, or 2)

        Returns:
            Dictionary with 'start' and 'stop' keys containing position lists
        """
        seq_upper = seq.upper()
        start_positions = []
        stop_positions = []

        # Determine if sequence is RNA
        is_rna = 'U' in seq_upper

        # Scan for codons in the specified frame
        for i in range(frame, len(seq) - 2, 3):
            codon = seq_upper[i:i+3]

            if codon in self.start_codons:
                start_positions.append(i)
            elif codon in self.stop_codons:
                stop_positions.append(i)

        return {
            'start': start_positions,
            'stop': stop_positions
        }

    def find_all_codons(self, seq: str) -> Dict[int, Dict[str, List[int]]]:
        """Find codons in all three reading frames.

        Args:
            seq: DNA or RNA sequence

        Returns:
            Dictionary mapping frame (0,1,2) to codon positions
        """
        results = {}
        for frame in range(3):
            results[frame] = self.find_codons(seq, frame)
        return results

    def translate_dna(self, dna_seq: str, to_single_letter: bool = True) -> str:
        """Translate DNA sequence to amino acids.

        Args:
            dna_seq: DNA sequence (should be multiple of 3)
            to_single_letter: If True, use single-letter AA codes

        Returns:
            Amino acid sequence
        """
        # Remove gaps
        dna_seq = dna_seq.replace('-', '')

        # Determine codon table based on sequence type
        seq_upper = dna_seq.upper()
        vc = infer_vocab(seq_upper)
        codon_table = bio.get_codon_table(vc)

        amino_acids = []
        for i in range(0, len(dna_seq) - 2, 3):
            codon = seq_upper[i:i+3]
            if 'N' in codon:
                amino_acids.append('X')  # Unknown
            else:
                amino_acids.append(codon_table.get(codon, 'X'))

        return ''.join(amino_acids)

    def find_orfs(self, seq: str, min_length: int = 100) -> List[Tuple[int, int]]:
        """Find open reading frames (ORFs) in a sequence.

        Args:
            seq: DNA or RNA sequence
            min_length: Minimum ORF length in nucleotides

        Returns:
            List of (start, end) tuples for ORFs
        """
        orfs = []
        seq_upper = seq.upper()

        for frame in range(3):
            # Find start codons in this frame
            start_codons = []
            for i in range(frame, len(seq) - 2, 3):
                codon = seq_upper[i:i+3]
                if codon in self.start_codons:
                    start_codons.append(i)

            # For each start codon, find the next stop codon
            for start in start_codons:
                for i in range(start + 3, len(seq) - 2, 3):
                    codon = seq_upper[i:i+3]
                    if codon in self.stop_codons:
                        if i + 3 - start >= min_length:
                            orfs.append((start, i + 3))
                        break

        return orfs

    def find_motifs(self, seq: str, motifs: Dict[str, str]) -> Dict[str, List[int]]:
        """Find motif occurrences in a sequence.

        Args:
            seq: DNA or RNA sequence
            motifs: Dictionary mapping motif names to regex patterns

        Returns:
            Dictionary mapping motif names to lists of match positions
        """
        results = {}
        seq_upper = seq.upper()

        for motif_name, pattern in motifs.items():
            positions = []
            for match in re.finditer(pattern, seq_upper):
                positions.append(match.start())
            if positions:
                results[motif_name] = positions

        return results

    def align_sequences(self, seq1: str, seq2: str, gap_penalty: int = -1,
                        mismatch_penalty: int = -1, match_score: int = 2) -> Tuple[str, str]:
        """Simple pairwise sequence alignment.

        Args:
            seq1: First sequence
            seq2: Second sequence
            gap_penalty: Penalty for gaps
            mismatch_penalty: Penalty for mismatches
            match_score: Score for matches

        Returns:
            Tuple of aligned sequences with gaps
        """
        # This is a placeholder for simple alignment
        # In practice, you might want to use more sophisticated algorithms
        # or integrate with the existing seqs.py alignment functions
        if len(seq1) == len(seq2):
            return seq1, seq2

        # Simple gap padding for different lengths
        if len(seq1) < len(seq2):
            seq1 = seq1 + '-' * (len(seq2) - len(seq1))
        else:
            seq2 = seq2 + '-' * (len(seq1) - len(seq2))

        return seq1, seq2

    def calculate_gc_content(self, seq: str, window_size: int = None) -> float:
        """Calculate GC content of a sequence.

        Args:
            seq: DNA sequence
            window_size: If provided, calculate windowed GC content

        Returns:
            GC content as a fraction (0-1)
        """
        seq_upper = seq.upper()

        if window_size is None:
            # Overall GC content
            gc_count = seq_upper.count('G') + seq_upper.count('C')
            total = len(seq_upper.replace('N', '').replace('-', ''))
            return gc_count / total if total > 0 else 0.0

        # Windowed GC content
        gc_values = []
        for i in range(len(seq) - window_size + 1):
            window = seq_upper[i:i+window_size]
            gc_count = window.count('G') + window.count('C')
            total = len(window.replace('N', '').replace('-', ''))
            gc_values.append(gc_count / total if total > 0 else 0.0)

        return gc_values

    def find_repeats(self, seq: str, min_length: int = 3,
                     max_length: int = 10) -> List[Tuple[str, List[int]]]:
        """Find tandem repeats in a sequence.

        Args:
            seq: DNA or RNA sequence
            min_length: Minimum repeat unit length
            max_length: Maximum repeat unit length

        Returns:
            List of (repeat_unit, positions) tuples
        """
        repeats = []
        seq_upper = seq.upper()

        for unit_len in range(min_length, max_length + 1):
            for i in range(len(seq) - unit_len * 2 + 1):
                unit = seq_upper[i:i+unit_len]
                # Check if unit repeats at least once
                if seq_upper[i+unit_len:i+unit_len*2] == unit:
                    # Count total repetitions
                    count = 2
                    while (i + unit_len * count < len(seq) and
                           seq_upper[i+unit_len*count:i+unit_len*(count+1)] == unit):
                        count += 1

                    if count >= 2:  # At least 2 repetitions
                        repeats.append((unit, [i, i + unit_len * count]))

        return repeats