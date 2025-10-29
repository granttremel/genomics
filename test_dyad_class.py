#!/usr/bin/env python3
"""
Test script for the new Dyad class.
"""

import sys
sys.path.insert(0, '/home/grant/genomics-prj')

from ggene.motifs.dyad_class import Dyad
from ggene.motifs import dyad

# Test sequence with a known dyad
test_seq = "GGGCAUUUGGGCCCCAAAAUGGCCCAAA"
#           -----         -----       (GGGCA and reverse complement UGCCC)

print("Testing Dyad class")
print("=" * 50)
print(f"Test sequence: {test_seq}\n")

# Create a Dyad instance from the old tuple format
# (stem_length=5, stem_start=0, loop_length=8, reverse_stem_start=13)
old_tuple = (5, 0, 8, 13)
d1 = Dyad.from_tuple(old_tuple, sequence=test_seq)

print(f"Created dyad from tuple: {d1}")
print(f"Back to tuple: {d1.to_tuple()}")
print(f"Total length: {d1.total_length}")
print(f"Center position: {d1.center_position}")
print()

# Extract sequences
print("Extracted sequences:")
print(f"Full dyad: {d1.extract_sequence()}")
stems = d1.extract_stems()
print(f"Stem 1: {stems[0]}")
print(f"Loop: {d1.extract_loop()}")
print(f"Stem 2: {stems[1]}")
print()

# Test validity
err = d1.test_dyad()
print(f"Test dyad validity: {err} errors (0 = perfect palindrome)")
print()

# Try to expand
print("Testing expand:")
expanded, err = d1.expand()
if expanded:
    print(f"Expanded: {expanded}")
    print(f"Expanded sequence: {expanded.extract_sequence()}")
else:
    print(f"Could not expand (error count: {err})")
print()

# Create another dyad for subdyad testing
d2 = Dyad(3, 2, 10, 15, sequence=test_seq)
print(f"Created second dyad: {d2}")
print(f"d2 sequence: {d2.extract_sequence()}")
print(f"Is d2 subdyad of d1? {d1.check_subdyad(d2)}")
print()

# Test mutual subdyad
has_super, superdyad = d1.check_mutual_subdyad(d2)
print(f"Do d1 and d2 have mutual superdyad? {has_super}")
if superdyad:
    print(f"Superdyad: {superdyad}")
print()

# Test shifting
shifted = d1.shift(5)
print(f"Shifted dyad by 5 positions: {shifted}")
print()

# Compare with existing functions
print("Comparing with existing dyad module functions:")
print("=" * 50)

# Test the old test_dyad function vs new method
old_err = dyad.test_dyad(test_seq, old_tuple)
new_err = d1.test_dyad()
print(f"Old test_dyad: {old_err} errors")
print(f"New test_dyad: {new_err} errors")
print(f"Match: {old_err == new_err}")
print()

# Test subdyad checking
old_subdyad = dyad.check_subdyad(old_tuple, d2.to_tuple())
new_subdyad = d1.check_subdyad(d2)
print(f"Old check_subdyad: {old_subdyad}")
print(f"New check_subdyad: {new_subdyad}")
print(f"Match: {old_subdyad == new_subdyad}")
print()

print("Test completed!")