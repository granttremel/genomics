#!/usr/bin/env python3
"""
Test the improved genome iterator implementation.
"""

import logging
from ggene.database.genomemanager import GenomeManager
from ggene.genome_iterator import GenomeIterator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_improved_sequence_extraction():
    """Test that the improved sequence extraction works correctly."""
    print("\n=== Testing Improved Sequence Extraction ===")
    
    gm = GenomeManager()
    
    # Test basic sequence retrieval
    chrom = 1
    start = 1000000
    end = 1000100
    
    # Get reference sequence directly
    ref_seq_direct = gm.get_sequence(chrom, start, end)
    print(f"Direct reference sequence length: {len(ref_seq_direct) if ref_seq_direct else 0}")
    
    # Get sequence through iterator
    iterator = GenomeIterator(
        gm,
        chrom=chrom,
        start=start,
        end=end,
        window_size=end - start + 1,
        integrate_variants=False
    )
    
    try:
        ref_seq_iter, _, _ = next(iterator)
        print(f"Iterator reference sequence length: {len(ref_seq_iter)}")
        
        # Check they match
        if ref_seq_direct and ref_seq_iter:
            match = ref_seq_direct == ref_seq_iter
            print(f"Sequences match: {match}")
            if not match:
                print(f"  Direct: {ref_seq_direct[:50]}...")
                print(f"  Iterator: {ref_seq_iter[:50]}...")
    except StopIteration:
        print("Iterator did not return sequence")
    
    # Test with variant integration
    print("\n--- Testing with variant integration ---")
    iterator_var = GenomeIterator(
        gm,
        chrom=chrom,
        start=start,
        end=end,
        window_size=end - start + 1,
        integrate_variants=True
    )
    
    try:
        ref_seq, personal_seq, _ = next(iterator_var)
        print(f"Personal sequence length: {len(personal_seq)}")
        # Count differences
        if ref_seq and personal_seq:
            differences = sum(1 for r, p in zip(ref_seq, personal_seq) if r != p)
            print(f"Number of differences: {differences}")
            
            # Show first few differences
            for i, (r, p) in enumerate(zip(ref_seq, personal_seq)):
                if r != p:
                    print(f"  Position {start + i}: {r} -> {p}")
                    if i >= 5:  # Show only first 5
                        break
    except StopIteration:
        print("Iterator did not return sequence")


def test_feature_retrieval():
    """Test the improved feature retrieval with exclusion."""
    print("\n=== Testing Feature Retrieval ===")
    
    gm = GenomeManager()
    
    # Find a position with features
    chrom = 1
    test_position = 1040000
    
    # Get all features at position
    all_features = gm.get_features_at_position(chrom, test_position)
    print(f"All features at position {test_position}: {len(all_features)}")
    
    if all_features:
        feature_types = set(f.get('feature', 'unknown') for f in all_features)
        print(f"Feature types: {', '.join(feature_types)}")
    
    # Get features excluding certain types
    filtered_features = gm.get_features_at_position(
        chrom, test_position, 
        exclude_types=['intron', 'transcript']
    )
    print(f"Filtered features: {len(filtered_features)}")
    
    if filtered_features:
        feature_types = set(f.get('feature', 'unknown') for f in filtered_features)
        print(f"Remaining feature types: {', '.join(feature_types)}")


def test_variant_application():
    """Test the variant application in window."""
    print("\n=== Testing Variant Application ===")
    
    gm = GenomeManager()
    
    # Find a region with variants
    chrom = 1
    start = 1000000
    end = 1010000
    
    # Count variants in region
    variant_count = 0
    variant_positions = []
    
    try:
        for var in gm.vcf(gm._make_index(chrom, start, end)):
            if var.QUAL > 10:  # Quality filter
                variant_count += 1
                variant_positions.append(var.POS)
                if variant_count >= 10:  # Get first 10
                    break
    except Exception as e:
        logger.warning(f"Error counting variants: {e}")
    
    print(f"Found {variant_count} quality variants in region")
    
    if variant_positions:
        # Test iterator on region with known variants
        test_start = variant_positions[0] - 10
        test_end = variant_positions[0] + 10
        
        iterator = GenomeIterator(
            gm,
            chrom=chrom,
            start=test_start,
            end=test_end,
            stride=1,
            integrate_variants=True
        )
        
        print(f"\nChecking positions {test_start} to {test_end}:")
        variants_found = 0
        
        for pos_data in iterator:
            if pos_data.has_variant:
                variants_found += 1
                print(f"  Position {pos_data.position}: "
                      f"{pos_data.reference_base} -> {pos_data.personal_base}")
        
        print(f"Iterator found {variants_found} variant positions")


def test_feature_exclusion_in_iterator():
    """Test using feature exclusion in the iterator."""
    print("\n=== Testing Feature Exclusion in Iterator ===")
    
    gm = GenomeManager()
    
    # Create iterator with specific feature types
    iterator = GenomeIterator(
        gm,
        chrom=1,
        start=1000000,
        end=1001000,
        stride=100,
        track_features=True,
        feature_types=['gene', 'exon', 'CDS']  # Only track these types
    )
    
    feature_summary = {}
    
    for pos_data in iterator:
        for feature_type in pos_data.feature_types:
            feature_summary[feature_type] = feature_summary.get(feature_type, 0) + 1
    
    print(f"Feature type summary:")
    for ftype, count in feature_summary.items():
        print(f"  {ftype}: {count} occurrences")


def main():
    """Run all tests."""
    tests = [
        test_improved_sequence_extraction,
        test_feature_retrieval,
        test_variant_application,
        test_feature_exclusion_in_iterator
    ]
    
    print("=" * 60)
    print("IMPROVED GENOME ITERATOR TESTS")
    print("=" * 60)
    
    for test in tests:
        try:
            test()
        except Exception as e:
            logger.error(f"Test failed: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("Tests completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()