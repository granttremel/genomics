"""
Test script to verify the refactored correlation functions.
"""

import sys
sys.path.insert(0, '/home/grant/genomics-prj')

from ggene.seqs.process import (
    correlate, correlate_v2,
    correlate_longest_subseq, correlate_longest_subseq_v2,
    convolve_generator, correlate_general
)

# Test sequences
seq1 = "AUGCAUGCAUGCUUUAAAGGG"
seq2 = "AUGCAUGCAUGCUUUAAAGGG"
seq3 = "GCAUGCAUGCUUUAAA"

print("=" * 80)
print("Test 1: Convolve generator basic functionality")
print("=" * 80)
test_seq1 = "ABCDEF"
test_seq2 = "12345"
print(f"seq1: {test_seq1}")
print(f"seq2: {test_seq2}")
print("\nGenerated pairs (no scale):")
for shift, s1, s2, overlap in convolve_generator(test_seq1, test_seq2):
    print(f"  shift={shift:2d}, overlap={overlap}, s1='{s1}', s2='{s2}'")

print("\nGenerated pairs (scale=3):")
for shift, s1, s2, overlap in convolve_generator(test_seq1, test_seq2, scale=3):
    print(f"  shift={shift:2d}, overlap={overlap}, s1='{s1}', s2='{s2}'")

print()
print("=" * 80)
print("Test 2: Compare old correlate vs new correlate_v2")
print("=" * 80)
print(f"seq1: {seq1}")
print(f"seq2: {seq2}")
print()

# Old version
direct_old, rc_old = correlate(seq1, seq2, scale=10)
print(f"Old correlate (first 10 values):")
print(f"  Direct: {[f'{v:.3f}' for v in direct_old[:10]]}")
print(f"  RC:     {[f'{v:.3f}' for v in rc_old[:10]]}")
print()

# New version
direct_new, rc_new = correlate_v2(seq1, seq2, scale=10)
print(f"New correlate_v2 (first 10 values):")
print(f"  Direct: {[f'{v:.3f}' for v in direct_new[:10]]}")
print(f"  RC:     {[f'{v:.3f}' for v in rc_new[:10]]}")
print()

# Check if results match
matches_direct = sum(1 for a, b in zip(direct_old, direct_new) if abs(a - b) < 1e-6)
matches_rc = sum(1 for a, b in zip(rc_old, rc_new) if abs(a - b) < 1e-6)
print(f"Matching values: Direct={matches_direct}/{len(direct_old)}, RC={matches_rc}/{len(rc_old)}")
print()

print("=" * 80)
print("Test 3: Compare old correlate_longest_subseq vs new v2")
print("=" * 80)
print(f"seq1: {seq1}")
print(f"seq3: {seq3}")
print()

# Old version
runs_old, inds_old, shifts_old = correlate_longest_subseq(seq1, seq3)
print(f"Old correlate_longest_subseq (first 10 values):")
print(f"  Runs:   {runs_old[:10]}")
print(f"  Inds:   {inds_old[:10]}")
print(f"  Shifts: {shifts_old[:10]}")
print()

# New version
runs_new, inds_new, shifts_new = correlate_longest_subseq_v2(seq1, seq3)
print(f"New correlate_longest_subseq_v2 (first 10 values):")
print(f"  Runs:   {runs_new[:10]}")
print(f"  Inds:   {inds_new[:10]}")
print(f"  Shifts: {shifts_new[:10]}")
print()

print("=" * 80)
print("Test 4: Custom analysis function using correlate_general")
print("=" * 80)

def h_bond_analyzer(shift, seq1_subseq, seq2_subseq, overlap_len):
    """
    Example: count potential H-bonds (A-U and G-C pairs).
    Returns (h_bonds, gc_count, au_count)
    """
    h_bonds = 0
    gc_count = 0
    au_count = 0

    for s1_base, s2_base in zip(seq1_subseq, seq2_subseq):
        if (s1_base, s2_base) in [('A', 'U'), ('U', 'A')]:
            h_bonds += 2  # A-U has 2 H-bonds
            au_count += 1
        elif (s1_base, s2_base) in [('G', 'C'), ('C', 'G')]:
            h_bonds += 3  # G-C has 3 H-bonds
            gc_count += 1

    return h_bonds, gc_count, au_count

h_bonds, gc_counts, au_counts = correlate_general(seq1, seq2, h_bond_analyzer, scale=15)

print(f"Custom H-bond analysis (first 10 values):")
print(f"  H-bonds:   {h_bonds[:10]}")
print(f"  GC pairs:  {gc_counts[:10]}")
print(f"  AU pairs:  {au_counts[:10]}")
print()

# Find maximum H-bonding shift
max_h_bond_idx = h_bonds.index(max(h_bonds))
# Shift values go from -max_shift to +max_shift
max_shift = len(seq1) // 2
shift_at_max = max_h_bond_idx - max_shift

print(f"Maximum H-bonds: {max(h_bonds)} at shift={shift_at_max}")
print(f"  GC pairs: {gc_counts[max_h_bond_idx]}, AU pairs: {au_counts[max_h_bond_idx]}")
print()

print("All tests completed!")
