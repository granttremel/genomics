"""
Dyad class for RNA/DNA palindrome structures.

A dyad is a palindromic structure separated by a loop region.
"""

from typing import Optional, Tuple


class Dyad:

    def __init__(self,
                 stem_length: int,
                 stem_start: int,
                 loop_length: int,
                 reverse_stem_start: int,
                 ref_pos = 0,
                 sequence: Optional[str] = None):
        self.stem_length = stem_length
        self.stem_start = stem_start
        self.loop_length = loop_length
        self.reverse_stem_start = reverse_stem_start
        self.sequence = sequence
        self._dyad_seq = None  # Cache for extracted sequence

        # Validate consistency
        if self.reverse_stem_start != self.stem_start + self.stem_length + self.loop_length:
            raise ValueError(f"Inconsistent dyad positions: reverse_stem should be at "
                           f"{self.stem_start + self.stem_length + self.loop_length}, "
                           f"but is at {self.reverse_stem_start}")

    @classmethod
    def from_tuple(cls, dyad_tuple: Tuple[int, int, int, int], sequence: Optional[str] = None):
        """Create a Dyad from the old tuple format (stem_length, stem_start, loop_length, reverse_stem_start)."""
        return cls(*dyad_tuple, sequence=sequence)

    def to_tuple(self) -> Tuple[int, int, int, int]:
        """Convert back to tuple format for compatibility."""
        return (self.stem_length, self.stem_start, self.loop_length, self.reverse_stem_start)

    @property
    def end_position(self) -> int:
        """End position of the entire dyad in the sequence."""
        return self.reverse_stem_start + self.stem_length

    @property
    def total_length(self) -> int:
        """Total length of the dyad structure."""
        return 2 * self.stem_length + self.loop_length

    @property
    def loop_start(self) -> int:
        """Start position of the loop region."""
        return self.stem_start + self.stem_length

    @property
    def center_position(self) -> float:
        """Center position of the dyad."""
        return self.stem_start + self.total_length / 2

    def extract_sequence(self, seq: Optional[str] = None) -> str:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        # Cache the result if using stored sequence
        if seq is self.sequence and self._dyad_seq is not None:
            return self._dyad_seq

        dyad_seq = seq[self.stem_start:self.end_position]

        if seq is self.sequence:
            self._dyad_seq = dyad_seq

        return dyad_seq

    def extract_stems(self, seq: Optional[str] = None) -> Tuple[str, str]:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        first_stem = seq[self.stem_start:self.loop_start]
        second_stem = seq[self.reverse_stem_start:self.end_position]
        return first_stem, second_stem

    def extract_loop(self, seq: Optional[str] = None) -> str:
        """Extract the loop sequence."""
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        return seq[self.loop_start:self.reverse_stem_start]

    def test_dyad(self, seq: Optional[str] = None,
                  complement_map: Optional[dict] = None,
                  err_tol: Optional[int] = None) -> int:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        if complement_map is None:
            complement_map = {'A':'U','U':'A','C':'G','G':'C','N':'N'}

        first_stem, second_stem = self.extract_stems(seq)

        # Reverse complement the first stem
        first_rc = self._reverse_complement(first_stem, complement_map)

        return self._compare_sequences(first_rc, second_stem, err_tol)

    def check_subdyad(self, other: 'Dyad') -> bool:
        # Ensure self is the larger dyad
        if other.stem_length > self.stem_length:
            return other.check_subdyad(self)
        elif other.stem_length == self.stem_length:
            return self == other

        # Check containment
        return (other.stem_start >= self.stem_start and
                other.stem_start + other.stem_length <= self.stem_start + self.stem_length and
                other.reverse_stem_start >= self.reverse_stem_start and
                other.reverse_stem_start + other.stem_length <= self.reverse_stem_start + self.stem_length)

    def check_mutual_subdyad(self, other: 'Dyad') -> Tuple[bool, Optional['Dyad']]:
        # Ensure self is the larger dyad
        if other.stem_length > self.stem_length:
            return other.check_mutual_subdyad(self)

        # Check if they share the same center
        if self.center_position == other.center_position:
            # Create the minimal common superdyad
            stem_start = min(self.stem_start, other.stem_start)
            loop_length = min(self.loop_length, other.loop_length)
            reverse_stem_start = min(self.reverse_stem_start, other.reverse_stem_start)
            stem_length = reverse_stem_start - loop_length - stem_start

            superdyad = Dyad(stem_length, stem_start, loop_length, reverse_stem_start,
                           sequence=self.sequence)
            return True, superdyad

        return False, None

    def expand(self, seq: Optional[str] = None,
               complement_map: Optional[dict] = None) -> Tuple[Optional['Dyad'], int]:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        # Check bounds
        if self.stem_start == 0 or self.end_position >= len(seq):
            return None, 1

        expanded = Dyad(
            self.stem_length + 1,
            self.stem_start - 1,
            self.loop_length,
            self.reverse_stem_start,
            sequence=seq
        )

        err = expanded.test_dyad(seq, complement_map)
        return expanded, err

    def contract(self, seq: Optional[str] = None,
                 complement_map: Optional[dict] = None) -> Tuple[Optional['Dyad'], int]:
        if seq is None:
            seq = self.sequence

        # Check minimum size and loop length
        if self.stem_length <= 1 or self.loop_length <= 2:
            return None, 1

        contracted = Dyad(
            self.stem_length + 1,
            self.stem_start + 1,
            self.loop_length - 2,
            self.reverse_stem_start - 1,
            sequence=seq
        )

        err = contracted.test_dyad(seq, complement_map)
        return contracted, err

    def shift(self, offset: int) -> 'Dyad':
        return Dyad(
            self.stem_length,
            self.stem_start + offset,
            self.loop_length,
            self.reverse_stem_start + offset,
            sequence=self.sequence
        )

    @staticmethod
    def _reverse_complement(seq: str, complement_map: dict) -> str:
        """Helper to get reverse complement of a sequence."""
        complement = ''.join(complement_map.get(base, base) for base in seq)
        return complement[::-1]

    @staticmethod
    def _compare_sequences(s1: str, s2: str, err_tol: Optional[int] = None) -> int:
        """Helper to compare two sequences and count mismatches."""
        if len(s1) != len(s2):
            return len(s1)

        if err_tol is None:
            err_tol = len(s1)

        nerr = 0
        for a, b in zip(s1, s2):
            if a != b:
                nerr += 1
            if nerr > err_tol:
                return nerr
        return nerr

    def __eq__(self, other: 'Dyad') -> bool:
        """Check equality based on all position attributes."""
        if not isinstance(other, Dyad):
            return False
        return (self.stem_length == other.stem_length and
                self.stem_start == other.stem_start and
                self.loop_length == other.loop_length and
                self.reverse_stem_start == other.reverse_stem_start)

    def __repr__(self) -> str:
        """String representation of the dyad."""
        seq_info = f", has_seq={self.sequence is not None}"
        return (f"Dyad(stem={self.stem_length}, start={self.stem_start}, "
                f"loop={self.loop_length}, rev_start={self.reverse_stem_start}{seq_info})")

    def __str__(self) -> str:
        """Human-readable string representation."""
        if self.sequence:
            stems = self.extract_stems()
            loop = self.extract_loop()
            return f"Dyad: {stems[0]}-{loop}-{stems[1]} @ pos {self.stem_start}"
        return f"Dyad: {self.stem_length}bp stems, {self.loop_length}bp loop @ pos {self.stem_start}"
    
    
    
