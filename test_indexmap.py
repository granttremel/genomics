#!/usr/bin/env python3
"""Test script to debug IndexMap off-by-one errors"""

from ggene.seqs.align import IndexMap

def test_simple_insertion():
    """Test a simple insertion mutation"""
    print("=== Testing Simple Insertion ===")

    # Original sequence: ACGT
    orig_seq = "ACGT"
    im = IndexMap(len(orig_seq))

    # Insert "AA" at position 2 (between C and G)
    im.add_delta(2, 2)  # Insert 2 bases at position 2

    print(f"Original sequence: {orig_seq}")
    print(f"Expected mutated:  AC--GT (gaps for insertions)")
    print(f"Expected unmutated: ACAAGT")

    mutated = im.mutate(orig_seq)
    print(f"Actual mutated:    {mutated}")

    # Now let's test unmutate with the inserted sequence
    mutant_seq = "ACAAGT"
    unmutated = im.unmutate(mutant_seq)
    print(f"Actual unmutated:  {unmutated}")

    # Test coordinate mappings
    print("\nCoordinate mappings (global -> init, global -> mut):")
    for i in range(im.global_seq_len):
        init_i = im.global_to_init(i)
        mut_i = im.global_to_mut(i)
        print(f"  Global {i} -> Init {init_i}, Mut {mut_i}")

    return mutated, unmutated

def test_simple_deletion():
    """Test a simple deletion mutation"""
    print("\n=== Testing Simple Deletion ===")

    # Original sequence: ACGTAG
    orig_seq = "ACGTAG"
    im = IndexMap(len(orig_seq))

    # Delete 2 bases at position 2 (remove GT)
    im.add_delta(2, -2)  # Delete 2 bases at position 2

    print(f"Original sequence: {orig_seq}")
    print(f"Expected mutated:  ACGTAG")
    print(f"Expected unmutated: AC--AG (gaps for deletions)")

    mutated = im.mutate(orig_seq)
    print(f"Actual mutated:    {mutated}")

    # Now let's test unmutate with the deleted sequence
    mutant_seq = "ACAG"
    unmutated = im.unmutate(mutant_seq)
    print(f"Actual unmutated:  {unmutated}")

    # Test coordinate mappings
    print("\nCoordinate mappings (global -> init, global -> mut):")
    for i in range(im.global_seq_len):
        init_i = im.global_to_init(i)
        mut_i = im.global_to_mut(i)
        print(f"  Global {i} -> Init {init_i}, Mut {mut_i}")

    return mutated, unmutated

def test_mixed_mutations():
    """Test combination of insertions and deletions"""
    print("\n=== Testing Mixed Mutations ===")

    # Original sequence: ACGTAG
    orig_seq = "ACGTAG"
    im = IndexMap(len(orig_seq))

    # First: Delete 1 base at position 1 (remove C)
    im.add_delta(1, -1)
    # Then: Insert 2 bases at position 3 (after G)
    im.add_delta(3, 2)

    print(f"Original sequence: {orig_seq}")
    print(f"Mutations: Delete at pos 1, Insert 2 at pos 3")

    mutated = im.mutate(orig_seq)
    print(f"Mutated sequence:  {mutated}")

    # Create the actual mutant sequence
    # Original: ACGTAG -> Delete C -> AGTAG -> Insert XX at pos 3 -> AGTXXAG
    mutant_seq = "AGTTTAG"  # Using TT as insertion
    unmutated = im.unmutate(mutant_seq)
    print(f"Unmutated sequence: {unmutated}")

    # Test coordinate mappings
    print("\nCoordinate mappings:")
    print("  Global -> Init, Mut")
    for i in range(im.global_seq_len):
        init_i = im.global_to_init(i)
        mut_i = im.global_to_mut(i)
        print(f"  {i:6} -> {init_i:4}, {mut_i:3}")

    return mutated, unmutated

def test_roundtrip():
    """Test that mutate -> unmutate preserves information"""
    print("\n=== Testing Roundtrip ===")

    orig_seq = "ACGTACGT"
    im = IndexMap(len(orig_seq))

    # Add some mutations
    im.add_delta(2, 3)   # Insert 3 at position 2
    im.add_delta(5, -2)  # Delete 2 at position 5

    print(f"Original: {orig_seq}")

    # Apply mutations
    mutated = im.mutate(orig_seq)
    print(f"After mutate(): {mutated}")

    # Try to recover - but we need the actual mutant sequence
    # Original: ACGTACGT
    # After insert 3 at pos 2: AC---GTACGT
    # After delete 2 at pos 5: AC---GTACGT (positions shift)
    # The actual mutant would be: ACXXXGTGT (where XXX is insertion)

    # Let's manually create what the mutant should be
    # Insert AAA at position 2: ACAAAGTACGT
    # Then delete 2 at position 5 (in original coords): ACAAAGT
    expected_mutant = "ACAAAGT"

    unmutated = im.unmutate(expected_mutant)
    print(f"After unmutate(): {unmutated}")

    return mutated, unmutated

if __name__ == "__main__":
    test_simple_insertion()
    test_simple_deletion()
    test_mixed_mutations()
    test_roundtrip()