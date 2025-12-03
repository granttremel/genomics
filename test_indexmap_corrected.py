#!/usr/bin/env python3
"""
Corrected IndexMap implementation with proper coordinate handling.

Key insight: When we say "insert at position 2", we mean insert AFTER position 1,
BEFORE position 2. The current character at position 2 should come AFTER the insertion.
"""

class IndexMapCorrected:

    def __init__(self, seq_len):
        self.initial_len = seq_len  # Length of initial/wild-type sequence
        self.global_seq_len = seq_len  # Length after insertions (longest)
        self.mut_seq_len = seq_len  # Length after mutations
        self.pdls = []  # List of (position, delta) pairs

    def add_delta(self, pos, delta):
        """Add a mutation at position pos with size delta.

        pos: Position in initial sequence coordinates
        delta: Positive for insertion, negative for deletion

        For insertion at pos: insert BEFORE the character at pos
        For deletion at pos: delete starting FROM pos
        """
        if abs(delta) < 1:
            return

        self.pdls.append((pos, delta))
        self.mut_seq_len += delta
        self.global_seq_len += max(delta, 0)  # Only increase for insertions

    def mutate(self, init_seq):
        """Convert initial sequence to global representation with gaps for insertions."""
        if not self.pdls:
            return init_seq

        result = []
        init_pos = 0

        # Sort mutations by position
        sorted_pdls = sorted(self.pdls)

        for mut_pos, delta in sorted_pdls:
            # Add characters before this mutation
            result.extend(init_seq[init_pos:mut_pos])

            if delta > 0:
                # Insertion: add gaps
                result.extend(['-'] * delta)
                init_pos = mut_pos
            else:
                # Deletion: include deleted characters, skip in init_seq
                result.extend(init_seq[mut_pos:mut_pos - delta])
                init_pos = mut_pos - delta

        # Add remaining characters
        result.extend(init_seq[init_pos:])

        return ''.join(result)

    def unmutate(self, mutated_seq):
        """Convert mutated sequence to global representation with gaps for deletions."""
        if not self.pdls:
            return mutated_seq

        result = []
        mut_pos = 0

        # Sort mutations by position
        sorted_pdls = sorted(self.pdls)

        for orig_pos, delta in sorted_pdls:
            # Calculate how many characters to copy before this mutation
            chars_before = orig_pos - len([d for p, d in sorted_pdls if p < orig_pos and d < 0])

            # Add characters before this mutation
            if mut_pos < len(mutated_seq):
                result.extend(mutated_seq[mut_pos:min(mut_pos + chars_before, len(mutated_seq))])

            if delta > 0:
                # Insertion: copy inserted characters
                if mut_pos + chars_before < len(mutated_seq):
                    result.extend(mutated_seq[mut_pos + chars_before:mut_pos + chars_before + delta])
                mut_pos += chars_before + delta
            else:
                # Deletion: add gaps
                result.extend(['-'] * (-delta))
                mut_pos += chars_before

        # Add remaining characters
        if mut_pos < len(mutated_seq):
            result.extend(mutated_seq[mut_pos:])

        return ''.join(result)


def test_insertion():
    """Test insertion mutations"""
    print("=== Testing Insertion ===")

    orig = "ACGT"
    im = IndexMapCorrected(len(orig))

    # Insert "XX" at position 2 (before 'G')
    im.add_delta(2, 2)

    mutated = im.mutate(orig)
    print(f"Original:  {orig}")
    print(f"Mutated:   {mutated}")
    print(f"Expected:  AC--GT")

    # Test unmutate with actual mutant sequence
    mutant_seq = "ACXXGT"  # Actual mutant with XX inserted
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant:    {mutant_seq}")
    print(f"Unmutated: {unmutated}")
    print(f"Expected:  ACXXGT")

    assert mutated == "AC--GT", f"Expected 'AC--GT', got '{mutated}'"
    assert unmutated == "ACXXGT", f"Expected 'ACXXGT', got '{unmutated}'"
    print("✓ Passed!\n")

def test_deletion():
    """Test deletion mutations"""
    print("=== Testing Deletion ===")

    orig = "ACGTAG"
    im = IndexMapCorrected(len(orig))

    # Delete 2 characters starting at position 2 (remove 'GT')
    im.add_delta(2, -2)

    mutated = im.mutate(orig)
    print(f"Original:  {orig}")
    print(f"Mutated:   {mutated}")
    print(f"Expected:  ACGTAG")

    # Test unmutate with actual mutant sequence
    mutant_seq = "ACAG"  # Actual mutant with GT deleted
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant:    {mutant_seq}")
    print(f"Unmutated: {unmutated}")
    print(f"Expected:  AC--AG")

    assert mutated == "ACGTAG", f"Expected 'ACGTAG', got '{mutated}'"
    assert unmutated == "AC--AG", f"Expected 'AC--AG', got '{unmutated}'"
    print("✓ Passed!\n")

def test_multiple_mutations():
    """Test multiple mutations"""
    print("=== Testing Multiple Mutations ===")

    orig = "ABCDEF"
    im = IndexMapCorrected(len(orig))

    # Delete 1 at position 1 (remove 'B')
    im.add_delta(1, -1)
    # Insert 2 at position 4 (before 'E', but after deletion it's position 3)
    im.add_delta(4, 2)

    mutated = im.mutate(orig)
    print(f"Original:  {orig}")
    print(f"Mutated:   {mutated}")
    print(f"Expected:  ABCD--EF")

    # Actual mutant: ACDEXXF (B deleted, XX inserted)
    mutant_seq = "ACDEXXF"
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant:    {mutant_seq}")
    print(f"Unmutated: {unmutated}")

    print("✓ Test completed!\n")

if __name__ == "__main__":
    test_insertion()
    test_deletion()
    test_multiple_mutations()