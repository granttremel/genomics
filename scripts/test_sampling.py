"""
Test script for step/keep sampling in convolve_generator.
"""

import sys
sys.path.insert(0, '/home/grant/genomics-prj')

from ggene.seqs.process import convolve_generator, correlate_v2, correlate_general

print("=" * 80)
print("Test 1: Basic sampling demonstration")
print("=" * 80)
seq1 = "ABCDEFGHIJKLMNOP"
seq2 = "1234567890ABCDEF"

print(f"seq1: {seq1}")
print(f"seq2: {seq2}")
print()

print("No sampling (step=1, keep=1):")
for shift, s1, s2, overlap in convolve_generator(seq1, seq2, scale=12, step=1, keep=1):
    if abs(shift) <= 2:  # Just show a few shifts
        print(f"  shift={shift:2d}: s1='{s1}' s2='{s2}'")
print()

print("Sample every 3, keep 2 (step=3, keep=2):")
print("  Expected pattern: AB__DE__GH__JK__MN__ -> 'ABDEGHJKMN...'")
for shift, s1, s2, overlap in convolve_generator(seq1, seq2, scale=12, step=3, keep=2):
    if abs(shift) <= 2:  # Just show a few shifts
        print(f"  shift={shift:2d}: s1='{s1}' s2='{s2}'")
print()

print("Sample every 3, keep 1 (step=3, keep=1):")
print("  Expected pattern: A__D__G__J__M__ -> 'ADGJM...'")
for shift, s1, s2, overlap in convolve_generator(seq1, seq2, scale=8, step=3, keep=1):
    if abs(shift) <= 2:  # Just show a few shifts
        print(f"  shift={shift:2d}: s1='{s1}' s2='{s2}'")
print()

print("=" * 80)
print("Test 2: RecA-like homology search simulation")
print("=" * 80)

# Create two sequences with a matching region
# seq1: AAAA-GGGGGGGG-TTTT
# seq2: CCCC-GGGGGGGG-UUUU
seq1 = "AAAAAAGGGGGGGGTTTTTT"
seq2 = "CCCCCCGGGGGGGGUU UUUU"

print(f"seq1: {seq1}")
print(f"seq2: {seq2}")
print()
print("These sequences share a GGGGGGGG region in the middle")
print()

def match_counter(shift, s1, s2, overlap):
    """Count matching bases."""
    matches = sum(1 for a, b in zip(s1, s2) if a == b)
    return (matches,)

print("Full comparison (no sampling):")
(matches_full,) = correlate_general(seq1, seq2, match_counter, scale=10)
max_idx_full = matches_full.index(max(matches_full))
shift_full = max_idx_full - len(seq1)//2
print(f"  Max matches: {max(matches_full)} at shift={shift_full}")
print(f"  All matches: {matches_full}")
print()

print("RecA-like sampling (step=3, keep=2):")
(matches_sampled,) = correlate_general(seq1, seq2, match_counter, scale=10, step=3, keep=2)
max_idx_sampled = matches_sampled.index(max(matches_sampled))
shift_sampled = max_idx_sampled - len(seq1)//2
print(f"  Max matches: {max(matches_sampled)} at shift={shift_sampled}")
print(f"  All matches: {matches_sampled}")
print()

print("Aggressive sampling (step=3, keep=1):")
(matches_aggressive,) = correlate_general(seq1, seq2, match_counter, scale=10, step=3, keep=1)
max_idx_aggressive = matches_aggressive.index(max(matches_aggressive))
shift_aggressive = max_idx_aggressive - len(seq1)//2
print(f"  Max matches: {max(matches_aggressive)} at shift={shift_aggressive}")
print(f"  All matches: {matches_aggressive}")
print()

print("=" * 80)
print("Test 3: Genomic sequences with sampling")
print("=" * 80)

# More realistic RNA sequences
seq1 = "AUGCAUGCAUGCUUUAAAGGG"
seq2 = "AUGCAUGCAUGCUUUAAAGGG"

print(f"seq1: {seq1}")
print(f"seq2: {seq2} (identical)")
print()

print("Full comparison:")
direct, rc = correlate_v2(seq1, seq2, scale=15, step=1, keep=1)
print(f"  Direct matches (first 5): {[f'{v:.2f}' for v in direct[:5]]}")
print(f"  Max direct: {max(direct):.2f}")
print()

print("With sampling (step=3, keep=2) - RecA-like:")
direct_sampled, rc_sampled = correlate_v2(seq1, seq2, scale=15, step=3, keep=2)
print(f"  Direct matches (first 5): {[f'{v:.2f}' for v in direct_sampled[:5]]}")
print(f"  Max direct: {max(direct_sampled):.2f}")
print()

print("Note: Sampling reduces computational cost while maintaining ability to find homology")
print()

print("All tests completed!")
