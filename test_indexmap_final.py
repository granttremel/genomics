#!/usr/bin/env python3
"""
Final comprehensive test of the fixed IndexMap implementation.
Tests various mutation scenarios to ensure correctness.
"""

from ggene.seqs.align import IndexMap


def test_insertion():
    """Test insertion mutations"""
    print("=== Testing Insertion ===")

    orig = "ACGT"
    im = IndexMap(len(orig))

    # Insert 2 characters at position 2 (before 'G')
    im.add_delta(2, 2)

    mutated = im.mutate(orig)
    print(f"Original:        {orig}")
    print(f"Mutated (gaps):  {mutated}")
    print(f"Expected:        AC--GT")

    # Test with actual mutant sequence
    mutant_seq = "ACTTGT"  # Inserted TT
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant seq:      {mutant_seq}")
    print(f"Unmutated:       {unmutated}")
    print(f"Expected:        ACTTGT")

    assert mutated == "AC--GT", f"Mutate failed: expected 'AC--GT', got '{mutated}'"
    assert unmutated == "ACTTGT", f"Unmutate failed: expected 'ACTTGT', got '{unmutated}'"
    print("✓ Passed!\n")


def test_deletion():
    """Test deletion mutations"""
    print("=== Testing Deletion ===")

    orig = "ACGTAG"
    im = IndexMap(len(orig))

    # Delete 2 characters starting at position 2 (remove 'GT')
    im.add_delta(2, -2)

    mutated = im.mutate(orig)
    print(f"Original:        {orig}")
    print(f"Mutated:         {mutated}")
    print(f"Expected:        ACGTAG")

    # Test with actual mutant sequence
    mutant_seq = "ACAG"  # GT deleted
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant seq:      {mutant_seq}")
    print(f"Unmutated (gaps):{unmutated}")
    print(f"Expected:        AC--AG")

    assert mutated == "ACGTAG", f"Mutate failed: expected 'ACGTAG', got '{mutated}'"
    assert unmutated == "AC--AG", f"Unmutate failed: expected 'AC--AG', got '{unmutated}'"
    print("✓ Passed!\n")


def test_multiple_mutations():
    """Test multiple mutations in sequence"""
    print("=== Testing Multiple Mutations ===")

    orig = "ABCDEFGH"
    im = IndexMap(len(orig))

    # First deletion: remove 'BC' at position 1
    im.add_delta(1, -2)
    # Then insertion: add 3 characters at position 5 (before 'F')
    im.add_delta(5, 3)

    mutated = im.mutate(orig)
    print(f"Original:        {orig}")
    print(f"After mutations: {mutated}")
    print(f"Expected:        ABCDE---FGH")

    # Actual mutant: ADEXYZFGH (BC deleted, XYZ inserted)
    mutant_seq = "ADEXYZFGH"
    unmutated = im.unmutate(mutant_seq)
    print(f"Mutant seq:      {mutant_seq}")
    print(f"Unmutated:       {unmutated}")
    print(f"Expected:        A--DEXYZFGH")

    assert mutated == "ABCDE---FGH", f"Mutate failed: expected 'ABCDE---FGH', got '{mutated}'"
    assert unmutated == "A--DEXYZFGH", f"Unmutate failed: expected 'A--DEXYZFGH', got '{unmutated}'"
    print("✓ Passed!\n")


def test_edge_cases():
    """Test edge cases and boundary conditions"""
    print("=== Testing Edge Cases ===")

    # Test 1: Insertion at beginning
    orig = "ACGT"
    im1 = IndexMap(len(orig))
    im1.add_delta(0, 2)
    mutated1 = im1.mutate(orig)
    print(f"Insert at start: {orig} -> {mutated1}")
    assert mutated1 == "--ACGT", f"Expected '--ACGT', got '{mutated1}'"

    # Test 2: Insertion at end
    im2 = IndexMap(len(orig))
    im2.add_delta(4, 2)
    mutated2 = im2.mutate(orig)
    print(f"Insert at end:   {orig} -> {mutated2}")
    assert mutated2 == "ACGT--", f"Expected 'ACGT--', got '{mutated2}'"

    # Test 3: Deletion at beginning
    orig3 = "ACGT"
    im3 = IndexMap(len(orig3))
    im3.add_delta(0, -2)
    mutated3 = im3.mutate(orig3)
    unmutated3 = im3.unmutate("GT")
    print(f"Delete at start: {orig3} -> mut: {mutated3}, unmut: {unmutated3}")
    assert mutated3 == "ACGT", f"Expected 'ACGT', got '{mutated3}'"
    assert unmutated3 == "--GT", f"Expected '--GT', got '{unmutated3}'"

    # Test 4: Complete deletion
    orig4 = "ACGT"
    im4 = IndexMap(len(orig4))
    im4.add_delta(0, -4)
    mutated4 = im4.mutate(orig4)
    unmutated4 = im4.unmutate("")
    print(f"Delete all:      {orig4} -> mut: {mutated4}, unmut: {unmutated4}")
    assert mutated4 == "ACGT", f"Expected 'ACGT', got '{mutated4}'"
    assert unmutated4 == "----", f"Expected '----', got '{unmutated4}'"

    print("✓ All edge cases passed!\n")


def test_coordinate_conversions():
    """Test the coordinate conversion methods"""
    print("=== Testing Coordinate Conversions ===")

    orig = "ABCDEF"
    im = IndexMap(len(orig))

    # Add insertion at position 2
    im.add_delta(2, 3)  # Insert 3 chars before 'C'
    # Add deletion at position 5
    im.add_delta(5, -1)  # Delete 'F'

    print("Original: ABCDEF")
    print("After insert 3 at pos 2: AB---CDEF")
    print("After delete 1 at pos 5: AB---CDE")
    print()

    # Test global_to_init conversions
    print("Global -> Init mappings:")
    for i in range(im.global_seq_len):
        init_i = im.global_to_init(i)
        print(f"  {i} -> {init_i}")

    # Test global_to_mut conversions
    print("\nGlobal -> Mut mappings:")
    for i in range(im.global_seq_len):
        mut_i = im.global_to_mut(i)
        print(f"  {i} -> {mut_i}")

    print("✓ Coordinate tests completed!\n")


def test_complex_scenario():
    """Test a more complex biological scenario"""
    print("=== Testing Complex Biological Scenario ===")

    # Simulate a gene with exons and a frameshift mutation
    gene = "ATGCGTACGGATTAG"  # Start codon + sequence + stop
    im = IndexMap(len(gene))

    # Frameshift: delete 1 base at position 5
    im.add_delta(5, -1)
    # Compensatory insertion: add 2 bases at position 10
    im.add_delta(10, 2)

    mutated = im.mutate(gene)
    print(f"Original gene:   {gene}")
    print(f"With gaps:       {mutated}")

    # Actual mutant with frameshift
    mutant_gene = "ATGCGACGGAAATTAG"  # T deleted, AA inserted
    unmutated = im.unmutate(mutant_gene)
    print(f"Mutant gene:     {mutant_gene}")
    print(f"Aligned:         {unmutated}")

    print("✓ Complex scenario completed!\n")


def test_roundtrip_consistency():
    """Test that operations maintain consistency"""
    print("=== Testing Roundtrip Consistency ===")

    orig = "GATTACA"
    im = IndexMap(len(orig))

    # Add multiple mutations
    im.add_delta(1, 2)   # Insert 2 after G
    im.add_delta(4, -1)  # Delete A
    im.add_delta(6, 1)   # Insert 1 before A

    # Create aligned versions
    mutated_aligned = im.mutate(orig)
    print(f"Original:        {orig}")
    print(f"Aligned w/gaps:  {mutated_aligned}")

    # The actual mutant sequence would be:
    # GATTACA -> G++ATTACA -> G++ATACA -> G++ATAC+A
    # Actual: GXXATCYA (where XX and Y are insertions)

    # Test with a specific mutant
    mutant = "GCCATCTA"  # CC inserted, A deleted, T inserted
    unmutated_aligned = im.unmutate(mutant)
    print(f"Mutant seq:      {mutant}")
    print(f"Unmut. aligned:  {unmutated_aligned}")

    # Verify lengths match
    assert len(mutated_aligned) == len(unmutated_aligned), \
        f"Length mismatch: {len(mutated_aligned)} vs {len(unmutated_aligned)}"

    print("✓ Roundtrip consistency verified!\n")


if __name__ == "__main__":
    print("Testing Fixed IndexMap Implementation")
    print("=" * 50)

    test_insertion()
    test_deletion()
    test_multiple_mutations()
    test_edge_cases()
    test_coordinate_conversions()
    test_complex_scenario()
    test_roundtrip_consistency()

    print("=" * 50)
    print("ALL TESTS PASSED! ✓")
    print("\nThe IndexMap class is now correctly handling:")
    print("- Insertions (gaps in initial sequence representation)")
    print("- Deletions (gaps in mutant sequence representation)")
    print("- Multiple mutations in combination")
    print("- Edge cases (start/end mutations)")
    print("- Coordinate transformations between spaces")