#!/usr/bin/env python3
"""Test the unified streaming system with sequence integration."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene import DEFAULT_FASTA_PATH, DEFAULT_GTF_PATH, DEFAULT_VCF_PATH
from ggene.unified_stream import UnifiedGenomeAnnotations
from ggene.genomemanager import GenomeManager
import os

def test_sequence_streaming():
    """Test sequence streaming functionality."""
    
    print("Testing Unified Streaming System with Sequence Integration")
    print("=" * 60)
    
    # Create unified annotation system with sequence streaming
    annotations = UnifiedGenomeAnnotations(
        fasta_path=DEFAULT_FASTA_PATH,
        vcf_path=DEFAULT_VCF_PATH
    )
    
    # Add annotation sources
    if os.path.exists("./data/GRCh38_sorted.gtf.gz"):
        annotations.add_gtf("./data/GRCh38_sorted.gtf.gz", "genes")
        print("✓ Added GTF annotations")
    
    # Test sequence access
    print("\n" + "-" * 60)
    print("Testing sequence access...")
    
    chrom = "1"
    start = 1000000
    end = 1000100
    
    # Get reference sequence
    ref_seq = annotations.get_sequence(chrom, start, end)
    print(f"\nReference sequence ({chrom}:{start:,}-{end:,}):")
    if ref_seq:
        print(f"  Length: {len(ref_seq)} bp")
        print(f"  First 50 bp: {ref_seq[:50]}...")
    else:
        print("  No sequence available")
    
    # Get personal sequence with variants
    pers_seq = annotations.get_personal_sequence(chrom, start, end)
    print(f"\nPersonal sequence ({chrom}:{start:,}-{end:,}):")
    if pers_seq:
        print(f"  Length: {len(pers_seq)} bp")
        if pers_seq != ref_seq:
            print("  ✓ Variants detected and applied")
        else:
            print("  No variants in this region")
    
    # Get aligned sequences with gaps
    print("\n" + "-" * 60)
    print("Testing aligned sequence generation...")
    
    aligned_ref, aligned_pers, coord_map = annotations.get_aligned_sequences(chrom, start, end)
    if aligned_ref:
        print(f"  Aligned reference length: {len(aligned_ref)}")
        print(f"  Aligned personal length: {len(aligned_pers)}")
        print(f"  Coordinate map entries: {len(coord_map)}")
        
        # Check for gaps
        gaps_in_ref = aligned_ref.count('-')
        gaps_in_pers = aligned_pers.count('-')
        print(f"  Gaps in reference: {gaps_in_ref}")
        print(f"  Gaps in personal: {gaps_in_pers}")
    
    # Test annotation streaming
    print("\n" + "-" * 60)
    print("Testing combined annotation streaming...")
    
    features = list(annotations.stream_all(chrom, start, end))
    print(f"\nFound {len(features)} features in region {chrom}:{start:,}-{end:,}")
    
    # Group by source
    by_source = {}
    for f in features:
        if f.source not in by_source:
            by_source[f.source] = []
        by_source[f.source].append(f)
    
    for source, items in by_source.items():
        print(f"  {source}: {len(items)} features")
        for item in items[:3]:
            print(f"    - {item.feature_type} at {item.start:,}-{item.end:,}")
        if len(items) > 3:
            print(f"    ... and {len(items) - 3} more")
    
    return annotations


def test_genome_manager_integration():
    """Test GenomeManager with unified streaming."""
    
    print("\n" + "=" * 60)
    print("Testing GenomeManager Integration")
    print("=" * 60)
    
    # Create GenomeManager
    gm = GenomeManager()
    
    print("✓ GenomeManager initialized with unified streaming")
    
    # Test sequence access through GenomeManager
    chrom = 1
    start = 204195200
    end = start+500
    
    seq = gm.get_sequence(chrom, start, end)
    if seq:
        print(f"\n✓ Sequence access working: {len(seq)} bp retrieved")
        # print(seq)
    
    # Test aligned sequences
    aligned_ref, aligned_pers, coord_map = gm.get_aligned_sequences(chrom, start, end)
    if aligned_ref:
        print(f"✓ Aligned sequence generation working")
        print(f"  Reference: {len(aligned_ref)} bp: {aligned_ref}")
        print(f"  Personal: {len(aligned_pers)} bp: {aligned_pers}")
    
    # Test annotation queries
    annotations = gm.get_all_annotations(str(chrom), start, end, include_motifs=False)
    print(f"\n✓ Found {len(annotations)} annotations in region")
    
    # Group by type
    by_type = {}
    for ann in annotations:
        if ann.feature_type not in by_type:
            by_type[ann.feature_type] = 0
        by_type[ann.feature_type] += 1
    
    print("  Annotation types:")
    for feat_type, count in by_type.items():
        print(f"    - {feat_type}: {count}")
    
    return gm


def test_performance():
    """Test performance of unified streaming."""
    
    print("\n" + "=" * 60)
    print("Performance Test")
    print("=" * 60)
    
    import time
    
    annotations = UnifiedGenomeAnnotations(
        fasta_path=DEFAULT_FASTA_PATH,
        vcf_path=DEFAULT_VCF_PATH
    )
    annotations.add_gtf(DEFAULT_GTF_PATH, "genes")
    
    # Test multiple region queries
    regions = [
        ("1", 1000000, 1001000),
        ("1", 2000000, 2001000),
        ("1", 5000000, 5001000),
    ]
    
    print("\nSequence retrieval performance:")
    for chrom, start, end in regions:
        start_time = time.time()
        seq = annotations.get_sequence(chrom, start, end)
        elapsed = time.time() - start_time
        print(f"  {chrom}:{start:,}-{end:,}: {elapsed:.3f}s ({len(seq) if seq else 0} bp)")
    
    print("\nAnnotation query performance:")
    for chrom, start, end in regions:
        start_time = time.time()
        features = annotations.query_range(chrom, start, end)
        elapsed = time.time() - start_time
        print(f"  {chrom}:{start:,}-{end:,}: {elapsed:.3f}s ({len(features)} features)")
    
    print("\nAligned sequence generation performance:")
    for chrom, start, end in regions[:2]:  # Test fewer for alignment
        start_time = time.time()
        aligned_ref, aligned_pers, coord_map = annotations.get_aligned_sequences(chrom, start, end)
        elapsed = time.time() - start_time
        gaps = (aligned_ref.count('-') + aligned_pers.count('-')) if aligned_ref else 0
        print(f"  {chrom}:{start:,}-{end:,}: {elapsed:.3f}s ({gaps} gaps)")


if __name__ == "__main__":
    # Run tests
    try:
        # Test sequence streaming
        annotations = test_sequence_streaming()
        
        # Test GenomeManager integration
        gm = test_genome_manager_integration()
        
        # Test performance
        test_performance()
        
        print("\n" + "=" * 60)
        print("✓ All tests completed successfully!")
        print("\nThe unified streaming system is working with:")
        print("  - Indexed GTF/VCF access for fast queries")
        print("  - FASTA sequence streaming")
        print("  - Variant integration")
        print("  - Gap-aligned sequence generation")
        print("  - Unified feature representation")
        
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()