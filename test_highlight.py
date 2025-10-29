#!/usr/bin/env python3
"""
Test the improved highlight_sequences function.
"""

import sys
sys.path.insert(0, '/home/grant/genomics-prj')

from ggene.motifs.utils import highlight_sequences

# Test sequence
test_seq = "AUGCAUGCGGCUAAUGCUAUGCCGAUGGCAUGCUAAA"
print(f"Original sequence:\n{test_seq}\n")

# Test case 1: Non-overlapping sequences
print("Test 1: Non-overlapping subsequences")
subseqs1 = ["AUGC", "GGCU", "CGAU"]
result1 = highlight_sequences(test_seq, subseqs1, min_len=4)
print(result1)
print()

# Test case 2: Overlapping sequences
print("Test 2: Overlapping subsequences")
subseqs2 = ["AUGCAU", "GCAUGC", "CAUGCG"]
result2 = highlight_sequences(test_seq, subseqs2, min_len=5, debug=True)
print(result2)
print()

# Test case 3: Complex overlaps
print("Test 3: Multiple overlapping patterns")
subseqs3 = ["AUGC", "GCAU", "CAUG", "AUGCA", "UGCAU"]
result3 = highlight_sequences(test_seq, subseqs3, min_len=4)
print(result3)
print()

# Test case 4: Show color palette
print("Test 4: Color palette demonstration")
print("256 Color ANSI Palette Sample:")
for i in range(16, 256, 8):
    row = []
    for j in range(8):
        color_num = i + j
        if color_num <= 255:
            row.append(f"\x1b[38;5;{color_num}m{color_num:3d}\x1b[0m")
    print(" ".join(row))
print()

# Explanation of the color ranges
print("Color ranges in 256-color mode:")
print("  0-15:   Standard terminal colors (0 = black)")
print("  16-231: 6x6x6 RGB color cube")
print("  232-255: Grayscale ramp (232 = darkest, 255 = lightest)")
print()
print("Our function uses colors 20-230 to ensure good visibility!")