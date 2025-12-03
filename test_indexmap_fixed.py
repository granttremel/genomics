#!/usr/bin/env python3
"""Fixed IndexMap implementation with correct gap handling"""

class IndexMapFixed:

    def __init__(self, seq_len):
        self.global_seq_len = seq_len
        self.seq_len = self.global_seq_len
        self.mut_seq_len = self.global_seq_len
        self.mapping = None
        self.pdls = []

    def add_delta(self, pos, delta):
        if abs(delta) < 1:
            return
        mapped_pos = pos
        self.pdls.append((mapped_pos, delta))

        self.mut_seq_len += delta
        self.global_seq_len += max(delta, 0)

    def global_to_init(self, i):
        """Convert global coordinate to initial/wild-type coordinate"""
        ii = i
        for p, d in self.pdls:
            dd = max(0, d)  # Only consider insertions
            if ii >= p and ii < p + dd:
                # Within insertion region - map to insertion point
                ii = p
            elif ii >= p + dd:
                # After insertion region - adjust by insertion size
                ii -= dd

        return ii

    def init_to_global(self, i):
        """Convert initial/wild-type coordinate to global coordinate"""
        ii = i
        for p, d in reversed(self.pdls):
            dd = max(0, d)  # Only consider insertions
            if ii >= p:
                ii += dd

        return ii

    def global_to_mut(self, i):
        """Convert global coordinate to mutant coordinate"""
        ii = i
        for p, d in self.pdls:
            dd = min(0, d)  # Only consider deletions (dd is negative)
            if ii >= p and ii < p - dd:
                # Within deletion region - map to deletion point
                ii = p
            elif ii >= p - dd:
                # After deletion region - adjust by deletion size
                ii += dd  # Add negative value (subtracts)

        return ii

    def mut_to_global(self, i):
        """Convert mutant coordinate to global coordinate"""
        ii = i
        for p, d in reversed(self.pdls):
            dd = min(0, d)  # Only consider deletions (negative)
            if ii >= p:
                ii -= dd  # Subtract negative value (adds)

        return ii

    def unmutate(self, mutated_seq):
        """Convert mutated sequence to global coordinate representation"""
        out_seq = []
        prev_im = None

        for i in range(self.global_seq_len):
            im = self.global_to_mut(i)

            # Check if we should insert a gap
            if prev_im is not None and im == prev_im:
                # We're in a deletion region - add gap
                out_seq.append('-')
            else:
                # Normal position - add character
                if im < len(mutated_seq):
                    out_seq.append(mutated_seq[im])
                else:
                    # Handle case where mutated_seq is shorter
                    out_seq.append('-')

            prev_im = im

        return "".join(out_seq)

    def mutate(self, init_seq):
        """Convert initial sequence to global coordinate representation"""
        out_seq = []
        prev_ii = None

        for i in range(self.global_seq_len):
            ii = self.global_to_init(i)

            # Check if we should insert a gap
            if prev_ii is not None and ii == prev_ii:
                # We're in an insertion region - add gap
                out_seq.append('-')
            else:
                # Normal position - add character
                if ii < len(init_seq):
                    out_seq.append(init_seq[ii])
                else:
                    # Handle case where init_seq is shorter
                    out_seq.append('-')

            prev_ii = ii

        return "".join(out_seq)


def test_simple_insertion_fixed():
    """Test a simple insertion mutation with fixed implementation"""
    print("=== Testing Simple Insertion (Fixed) ===")

    # Original sequence: ACGT
    orig_seq = "ACGT"
    im = IndexMapFixed(len(orig_seq))

    # Insert "AA" at position 2 (between C and G)
    im.add_delta(2, 2)  # Insert 2 bases at position 2

    print(f"Original sequence: {orig_seq}")
    print(f"Expected mutated:  AC--GT (gaps for insertions)")

    mutated = im.mutate(orig_seq)
    print(f"Actual mutated:    {mutated}")

    # Now let's test unmutate with the inserted sequence
    mutant_seq = "ACAAGT"
    print(f"Mutant sequence:   {mutant_seq}")
    print(f"Expected unmutated: ACAAGT")
    unmutated = im.unmutate(mutant_seq)
    print(f"Actual unmutated:  {unmutated}")

    # Verify correctness
    assert mutated == "AC--GT", f"Mutated should be 'AC--GT' but got '{mutated}'"
    assert unmutated == "ACAAGT", f"Unmutated should be 'ACAAGT' but got '{unmutated}'"
    print("✓ Test passed!")

def test_simple_deletion_fixed():
    """Test a simple deletion mutation with fixed implementation"""
    print("\n=== Testing Simple Deletion (Fixed) ===")

    # Original sequence: ACGTAG
    orig_seq = "ACGTAG"
    im = IndexMapFixed(len(orig_seq))

    # Delete 2 bases at position 2 (remove GT)
    im.add_delta(2, -2)  # Delete 2 bases at position 2

    print(f"Original sequence: {orig_seq}")
    print(f"Expected mutated:  ACGTAG")

    mutated = im.mutate(orig_seq)
    print(f"Actual mutated:    {mutated}")

    # Now let's test unmutate with the deleted sequence
    mutant_seq = "ACAG"
    print(f"Mutant sequence:   {mutant_seq}")
    print(f"Expected unmutated: AC--AG (gaps for deletions)")
    unmutated = im.unmutate(mutant_seq)
    print(f"Actual unmutated:  {unmutated}")

    # Verify correctness
    assert mutated == "ACGTAG", f"Mutated should be 'ACGTAG' but got '{mutated}'"
    assert unmutated == "AC--AG", f"Unmutated should be 'AC--AG' but got '{unmutated}'"
    print("✓ Test passed!")

def test_edge_cases():
    """Test edge cases and boundary conditions"""
    print("\n=== Testing Edge Cases ===")

    # Test insertion at beginning
    orig = "ACGT"
    im = IndexMapFixed(len(orig))
    im.add_delta(0, 2)
    mutated = im.mutate(orig)
    print(f"Insert at start: {orig} -> {mutated}")
    assert mutated == "--ACGT"

    # Test insertion at end
    im2 = IndexMapFixed(len(orig))
    im2.add_delta(4, 2)
    mutated2 = im2.mutate(orig)
    print(f"Insert at end:   {orig} -> {mutated2}")
    assert mutated2 == "ACGT--"

    # Test deletion at beginning
    orig3 = "ACGT"
    im3 = IndexMapFixed(len(orig3))
    im3.add_delta(0, -2)
    mutated3 = im3.mutate(orig3)
    unmutated3 = im3.unmutate("GT")
    print(f"Delete at start: {orig3} -> mutated: {mutated3}, unmutated: {unmutated3}")
    assert unmutated3 == "--GT"

    print("✓ All edge case tests passed!")

if __name__ == "__main__":
    test_simple_insertion_fixed()
    test_simple_deletion_fixed()
    test_edge_cases()