#!/usr/bin/env python3
"""Test the aligned coordinate system in GenomeIterator."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene.genomemanager import GenomeManager
from ggene.genome_iterator import GenomeIterator

def test_large_insertion():
    """Test the iterator with the large insertion at 1:204195223."""
    
    print("Loading genome data...")
    # Use default paths from GenomeManager
    gm = GenomeManager()
    
    # Test at the large insertion position
    chrom = '1'
    position = 204195210  # Start a bit before the insertion
    window_size = 80
    
    print(f"\nTesting at chromosome {chrom}, position {position} (insertion at 1:204195223)")
    print(f"Window size: {window_size}")
    print("-" * 80)
    
    # Create iterator with variants
    iterator = GenomeIterator(
        gm, 
        chrom, 
        position,
        position + window_size + 500,  # Get extra to account for insertion
        window_size=window_size,
        integrate_variants=True,
        track_features=True
    )
    
    # Get first window
    try:
        ref_seq, pers_seq, features = next(iterator)
        
        print(f"\nAligned sequences (length: {len(ref_seq)}):")
        print(f"Reference:  {ref_seq[:80]}")
        print(f"Personal:   {pers_seq[:80]}")
        
        # Count gaps
        ref_gaps = ref_seq.count('-')
        pers_gaps = pers_seq.count('-')
        print(f"\nGaps in reference: {ref_gaps}")
        print(f"Gaps in personal: {pers_gaps}")
        
        # Check if we have the expected window size
        if len(ref_seq) == window_size and len(pers_seq) == window_size:
            print(f"✓ Window size is correct: {window_size} bases")
        else:
            print(f"✗ Window size mismatch: ref={len(ref_seq)}, pers={len(pers_seq)}, expected={window_size}")
        
        # Check for variants
        variant_count = 0
        for i in range(min(len(ref_seq), len(pers_seq))):
            if ref_seq[i] != pers_seq[i]:
                variant_count += 1
        
        print(f"\nVariant positions: {variant_count}")
        
        # Show variant features
        variant_features = [f for f in features if f.get('feature') == 'variant']
        if variant_features:
            print(f"\nVariants in window:")
            for var in variant_features[:5]:  # Show first 5
                print(f"  - Position {var.get('start')}: {var.get('ref')} → {var.get('alt')}")
                ref_len = len(var.get('ref', ''))
                alt_len = len(var.get('alt', ''))
                if alt_len > ref_len:
                    print(f"    (Insertion of {alt_len - ref_len} bp)")
                elif ref_len > alt_len:
                    print(f"    (Deletion of {ref_len - alt_len} bp)")
        
        # Test moving forward
        print("\n" + "=" * 80)
        print("Testing iterator advancement...")
        
        for i in range(3):
            try:
                ref_seq, pers_seq, features = next(iterator)
                print(f"\nIteration {i+2}: Window size = {len(ref_seq)}")
                if len(ref_seq) != window_size:
                    print(f"  Warning: Size mismatch!")
            except StopIteration:
                print(f"\nReached end of sequence at iteration {i+2}")
                break
                
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

def test_small_indels():
    """Test with smaller indels to verify basic functionality."""
    
    print("\n" + "=" * 80)
    print("Testing with smaller indels...")
    
    # Create a simple test case
    gm = GenomeManager()
    
    # Find a region with small indels
    chrom = '1'
    position = 1000000
    window_size = 60
    
    iterator = GenomeIterator(
        gm,
        chrom,
        position,
        position + 200,
        window_size=window_size,
        integrate_variants=True
    )
    
    try:
        ref_seq, pers_seq, _ = next(iterator)
        print(f"\nWindow at {chrom}:{position}")
        print(f"Reference:  {ref_seq}")
        print(f"Personal:   {pers_seq}")
        print(f"Length check: ref={len(ref_seq)}, pers={len(pers_seq)}, expected={window_size}")
        
    except Exception as e:
        print(f"Error in small indel test: {e}")

if __name__ == "__main__":
    test_large_insertion()
    test_small_indels()
    print("\nTest completed!")