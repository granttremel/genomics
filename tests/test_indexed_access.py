#!/usr/bin/env python3
"""Test indexed access performance for unified streams."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

import time
from ggene.database.annotations import GTFStream, VCFStream
from ggene.database.index_utils import ensure_indexed, check_all_indexes

def test_indexed_vs_linear():
    """Compare performance of indexed vs linear access."""
    
    print("Testing indexed access performance...")
    print("=" * 60)
    
    # Check index status
    files = [
        "./data/GRCh38_sorted.gtf.gz",
        "./data/gt_vcf.gz"
    ]
    
    print("Checking index status:")
    status = check_all_indexes(files)
    for filepath, stat in status.items():
        print(f"  {filepath}: {stat}")
    
    # Test GTF access
    print("\n" + "-" * 60)
    print("Testing GTF stream...")
    
    gtf_stream = GTFStream("./data/GRCh38_sorted.gtf.gz")
    
    # Test a specific region query
    chrom = "1"
    start = 1000000
    end = 1100000
    
    print(f"\nQuerying region {chrom}:{start:,}-{end:,}")
    
    start_time = time.time()
    features = list(gtf_stream.stream(chrom, start, end))
    elapsed = time.time() - start_time
    
    print(f"Found {len(features)} features in {elapsed:.3f} seconds")
    
    if features:
        print("\nFirst 5 features:")
        for f in features[:5]:
            print(f"  {f.feature_type} at {f.start:,}-{f.end:,}: {f.name}")
    
    # Test VCF access
    print("\n" + "-" * 60)
    print("Testing VCF stream...")
    
    vcf_stream = VCFStream("./data/gt_vcf.gz")
    
    start_time = time.time()
    variants = list(vcf_stream.stream(chrom, start, end))
    elapsed = time.time() - start_time
    
    print(f"Found {len(variants)} variants in {elapsed:.3f} seconds")
    
    if variants:
        print("\nFirst 5 variants:")
        for v in variants[:5]:
            ref = v.attributes.get('ref', '?')
            alt = v.attributes.get('alt', ['?'])[0] if v.attributes.get('alt') else '?'
            print(f"  {v.start:,}: {ref} → {alt}")
    
    print("\n" + "=" * 60)
    print("Performance comparison:")
    print("  With indexed access: < 1 second")
    print("  Without indexed access: > 10 seconds (for large files)")
    print("\nMake sure your files are indexed for optimal performance!")


def test_random_access():
    """Test multiple random access queries."""
    
    print("\n" + "=" * 60)
    print("Testing random access performance...")
    
    gtf_stream = GTFStream("./data/GRCh38_sorted.gtf.gz")
    
    # Test multiple random regions
    regions = [
        ("1", 1000000, 1001000),
        ("2", 5000000, 5001000),
        ("X", 1000000, 1001000),
        ("1", 100000000, 100001000),
    ]
    
    total_time = 0
    total_features = 0
    
    for chrom, start, end in regions:
        start_time = time.time()
        features = list(gtf_stream.stream(chrom, start, end))
        elapsed = time.time() - start_time
        total_time += elapsed
        total_features += len(features)
        
        print(f"  {chrom}:{start:,}-{end:,}: {len(features)} features in {elapsed:.3f}s")
    
    print(f"\nTotal: {total_features} features in {total_time:.3f} seconds")
    print(f"Average: {total_time/len(regions):.3f} seconds per query")
    
    if total_time < 1:
        print("✓ Excellent! Indexed access is working properly.")
    else:
        print("⚠ Queries are slow. Check if files are properly indexed.")


if __name__ == "__main__":
    test_indexed_vs_linear()
    test_random_access()
    
    print("\n" + "=" * 60)
    print("Test complete!")
    print("\nTip: If queries are slow, create indexes with:")
    print("  tabix -p gff your_file.gtf.gz")
    print("  tabix -p vcf your_file.vcf.gz")
    print("  tabix -p bed your_file.bed.gz")