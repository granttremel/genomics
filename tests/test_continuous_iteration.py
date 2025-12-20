#!/usr/bin/env python3
"""Test continuous iteration through large insertions."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene.database.genome_manager import GenomeManager
from ggene.genome_iterator import GenomeIterator

def test_continuous_through_insertion():
    """Test that we iterate continuously through the large insertion."""
    
    print("Loading genome data...")
    gm = GenomeManager()
    
    # Test at the large insertion position
    chrom = '1'
    start_pos = 204195220  # Start just before the insertion
    window_size = 40
    stride = 20
    
    print(f"\nTesting continuous iteration through insertion at 1:204195223")
    print(f"Window size: {window_size}, Stride: {stride}")
    print("-" * 80)
    
    # Create iterator
    iterator = GenomeIterator(
        gm, 
        chrom, 
        start_pos,
        start_pos + 600,  # Cover the full 455bp insertion
        window_size=window_size,
        stride=stride,
        integrate_variants=True,
        track_features=False  # Simplify output
    )
    
    iteration = 0
    prev_end_seq = ""
    continuous = True
    
    print(f"\nIterating from position {start_pos}:")
    print("=" * 80)
    
    for ref_seq, pers_seq, _ in iterator:
        iteration += 1
        
        # Check continuity - the start of this window should overlap with end of previous
        if iteration > 1:
            overlap_size = window_size - stride
            expected_start = prev_end_seq[-overlap_size:] if len(prev_end_seq) >= overlap_size else prev_end_seq
            actual_start = pers_seq[:len(expected_start)]
            
            if expected_start and expected_start != actual_start:
                continuous = False
                print(f"\n❌ DISCONTINUITY at iteration {iteration}!")
                print(f"   Expected start: ...{expected_start[-20:]}")
                print(f"   Actual start:   {actual_start[:20]}...")
        
        # Count gaps to detect when we're in the insertion
        ref_gaps = ref_seq.count('-')
        pers_gaps = pers_seq.count('-')
        
        # Display information
        print(f"\nIteration {iteration} (pos ~{iterator.current_pos}):")
        print(f"  Ref:  {ref_seq[:30]}{'...' if len(ref_seq) > 30 else ''}")
        print(f"  Pers: {pers_seq[:30]}{'...' if len(pers_seq) > 30 else ''}")
        print(f"  Gaps: ref={ref_gaps}, pers={pers_gaps}")
        
        # Check if we're in the insertion
        if ref_gaps > window_size * 0.5:  # More than half gaps means we're in the insertion
            print(f"  📍 Inside large insertion!")
        
        prev_end_seq = pers_seq
        
        # Stop after enough iterations to see the pattern
        if iteration >= 15:
            print("\n... (stopping after 15 iterations)")
            break
    
    print("\n" + "=" * 80)
    if continuous:
        print("✅ Iteration is continuous - no sequence was skipped!")
    else:
        print("❌ Discontinuity detected - sequence was skipped!")
    
    return continuous

def test_position_advancement():
    """Test that position advances correctly through insertions."""
    
    print("\n" + "=" * 80)
    print("Testing position advancement...")
    
    gm = GenomeManager()
    
    chrom = '1'
    start_pos = 204195220
    window_size = 40
    stride = 20
    
    iterator = GenomeIterator(
        gm, 
        chrom, 
        start_pos,
        start_pos + 600,
        window_size=window_size,
        stride=stride,
        integrate_variants=True
    )
    
    positions = []
    for i, (ref_seq, pers_seq, _) in enumerate(iterator):
        positions.append(iterator.current_pos)
        
        # Check position advancement
        if i > 0:
            advancement = positions[i] - positions[i-1]
            ref_gaps = ref_seq[:stride].count('-')
            
            print(f"Iteration {i+1}: pos={iterator.current_pos}, "
                  f"advanced by {advancement}, gaps in stride={ref_gaps}")
        
        if i >= 10:
            break
    
    print("\nPosition sequence:", positions[:10])

if __name__ == "__main__":
    continuous = test_continuous_through_insertion()
    test_position_advancement()
    
    if continuous:
        print("\n✅ All tests passed!")
    else:
        print("\n❌ Tests failed - check implementation")