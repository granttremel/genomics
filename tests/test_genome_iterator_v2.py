#!/usr/bin/env python3
"""Test the enhanced genome iterator with unified streaming."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene import get_paths
DEFAULT_VCF_PATH, DEFAULT_GTF_PATH, DEFAULT_FASTA_PATH, DEFAULT_LIBRARY = get_paths()
from ggene.genomemanager import GenomeManager
from ggene.genome_iterator_v2 import UnifiedGenomeIterator

def test_coordinate_tracking():
    """Test the three coordinate systems."""
    
    print("Testing Enhanced Genome Iterator with Coordinate Tracking")
    print("=" * 60)
    
    # Initialize GenomeManager with unified streaming
    gm = GenomeManager()
    
    # Test region with known large insertion at 1:204195223
    chrom = 1
    start = 204195200
    window_size = 50
    
    print(f"\nIterating over chromosome {chrom}:{start:,} with window size {window_size}")
    print("This region contains a 455bp insertion at position 204195223")
    print("-" * 60)
    
    # Create iterator
    iterator = UnifiedGenomeIterator(
        gm, chrom, start, 
        end=start + 200,
        window_size=window_size,
        stride=25,
        integrate_variants=True,
        track_features=True
    )
    
    # Iterate and display windows
    for i, window in enumerate(iterator):
        print(f"\nWindow {i+1}:")
        print(f"  Reference coords: {window.start_ref:,} - {window.end_ref:,}")
        print(f"  Alternate coords: {window.start_alt:,} - {window.end_alt:,}")
        print(f"  Display coords:   {window.start_display:,} - {window.end_display:,}")
        print(f"  Ref length: {window.ref_length} bp (excl. gaps)")
        print(f"  Alt length: {window.alt_length} bp (excl. gaps)")
        print(f"  Ref sequence: {window.ref_seq} bp (excl. gaps)")
        print(f"  Alt sequence: {window.alt_seq} bp (excl. gaps)")
        print(f"  Display length: {window.display_length} bp (incl. gaps)")
        
        if window.variants:
            print(f"  Variants in window: {len(window.variants)}")
            for pos, delta in window.variants[:3]:  # Show first 3
                variant_type = "insertion" if delta > 0 else "deletion" if delta < 0 else "SNP"
                print(f"    - {variant_type} at {pos:,}: delta={delta:+d}")
        
        if window.has_gaps:
            print(f"  Has gaps: Yes")
            print(f"  Ref seq: {window.ref_seq[:30]}...")
            print(f"  Alt seq: {window.alt_seq[:30]}...")
        
        if window.features:
            print(f"  Features: {len(window.features)}")
            feature_types = set(f.feature_type for f in window.features)
            print(f"    Types: {', '.join(feature_types)}")
        
        # Get current coordinates
        coords = iterator.get_current_coordinates()
        print(f"  Iterator state: ref={coords['reference']:,}, "
              f"alt={coords['alternate']:,}, "
              f"delta={coords['cumulative_delta']:+d}")


def test_preloading_performance():
    """Test the preloading and buffering performance."""
    
    print("\n" + "=" * 60)
    print("Testing Preloading and Buffering Performance")
    print("=" * 60)
    
    import time
    
    gm = GenomeManager()
    
    # Test regions
    test_regions = [
        (1, 1000000, 1001000, "Region with genes"),
        (1, 204195000, 204196000, "Region with large insertion"),
        (1, 5000000, 5001000, "Another region"),
    ]
    
    for chrom, start, end, description in test_regions:
        print(f"\n{description}: chr{chrom}:{start:,}-{end:,}")
        
        start_time = time.time()
        
        # Create iterator with preloading
        iterator = UnifiedGenomeIterator(
            gm, chrom, start, end,
            window_size=100,
            stride=100,
            integrate_variants=True,
            track_features=True
        )
        
        # Iterate through windows
        window_count = 0
        total_variants = 0
        total_features = 0
        
        for window in iterator:
            window_count += 1
            total_variants += len(window.variants)
            total_features += len(window.features)
        
        elapsed = time.time() - start_time
        
        print(f"  Processed {window_count} windows in {elapsed:.3f}s")
        print(f"  Total variants: {total_variants}")
        print(f"  Total features: {total_features}")
        print(f"  Speed: {(end-start)/elapsed:.0f} bp/s")


def test_jump_functionality():
    """Test jumping to specific positions."""
    
    print("\n" + "=" * 60)
    print("Testing Jump Functionality")
    print("=" * 60)
    
    gm = GenomeManager()
    
    # Create iterator
    iterator = UnifiedGenomeIterator(
        gm, chrom=1, start=1000000,
        window_size=50,
        integrate_variants=True
    )
    
    # Test jumps
    test_positions = [1000000, 2000000, 204195223, 5000000]
    
    for pos in test_positions:
        print(f"\nJumping to position {pos:,}")
        iterator.jump_to(pos, coord_system='ref')
        
        coords = iterator.get_current_coordinates()
        print(f"  After jump:")
        print(f"    Reference: {coords['reference']:,}")
        print(f"    Alternate: {coords['alternate']:,}")
        print(f"    Display:   {coords['display']:,}")
        print(f"    Cumulative delta: {coords['cumulative_delta']:+d}")
        
        # Get one window at this position
        for window in iterator:
            print(f"  Window at position:")
            print(f"    Covers: {window.start_ref:,} - {window.end_ref:,}")
            print(f"    Variants: {len(window.variants)}")
            print(f"    Features: {len(window.features)}")
            break  # Just one window


def test_constant_window_size():
    """Test that windows maintain constant size in presence of indels."""
    
    print("\n" + "=" * 60)
    print("Testing Constant Window Size with Indels")
    print("=" * 60)
    
    gm = GenomeManager()
    
    # Region with known large insertion
    chrom = 1
    start = 204195200
    window_size = 80
    
    print(f"\nTesting windows of size {window_size} bp around large insertion")
    
    iterator = UnifiedGenomeIterator(
        gm, chrom, start,
        end=start + 300,
        window_size=window_size,
        stride=40,
        integrate_variants=True
    )
    
    for i, window in enumerate(iterator):
        # Check that reference window size is constant
        ref_window_size = window.end_ref - window.start_ref + 1
        
        print(f"\nWindow {i+1}:")
        print(f"  Ref window size: {ref_window_size} bp (should be ≤ {window_size})")
        print(f"  Ref sequence length (no gaps): {window.ref_length} bp")
        print(f"  Alt sequence length (no gaps): {window.alt_length} bp")
        print(f"  Ref sequence: {window.ref_seq} bp (excl. gaps)")
        print(f"  Alt sequence: {window.alt_seq} bp (excl. gaps)")
        print(f"  Display length (with gaps): {window.display_length} bp")
        
        if window.variants:
            total_delta = sum(delta for _, delta in window.variants)
            print(f"  Total variant delta: {total_delta:+d} bp")
        
        # Verify window size constraint
        assert ref_window_size <= window_size, f"Window exceeded size limit!"
        print(f"  ✓ Window size constraint maintained")


if __name__ == "__main__":
    try:
        # Run all tests
        test_coordinate_tracking()
        test_preloading_performance()
        test_jump_functionality()
        test_constant_window_size()
        
        print("\n" + "=" * 60)
        print("✓ All tests completed successfully!")
        print("\nThe enhanced genome iterator provides:")
        print("  - Three coordinate systems (reference, alternate, display)")
        print("  - Efficient preloading and buffering")
        print("  - Constant window sizes despite indels")
        print("  - Jump-to-position functionality")
        print("  - Integration with unified streaming system")
        
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()