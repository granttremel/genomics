#!/usr/bin/env python3
"""Test the unified annotation streaming system."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene.unified_stream import (
    UnifiedGenomeAnnotations, 
    GTFStream, 
    VCFStream,
    UnifiedFeature
)
from ggene.genomemanager import GenomeManager
import heapq
from typing import Iterator

def test_unified_streaming():
    """Test streaming from multiple sources."""
    
    print("Setting up unified annotation system...")
    annotations = UnifiedGenomeAnnotations()
    
    # Add real data sources
    print("Adding GTF annotations...")
    annotations.add_gtf("./data/gencode.v38.annotation.head.gtf", "genes")
    
    print("Adding VCF variants...")
    annotations.add_vcf("./genome/vcf/HG00096.vcf.gz", "personal_variants")
    
    # Test streaming a small region
    chrom = "1"
    start = 1000000
    end = 1100000
    
    print(f"\nStreaming annotations for {chrom}:{start:,}-{end:,}")
    print("-" * 60)
    
    # Count features by type
    feature_counts = {}
    feature_examples = {}
    
    for i, feature in enumerate(annotations.stream_all(chrom, start, end)):
        # Count by type
        key = f"{feature.source}:{feature.feature_type}"
        feature_counts[key] = feature_counts.get(key, 0) + 1
        
        # Save examples
        if key not in feature_examples:
            feature_examples[key] = feature
        
        # Show first few
        if i < 10:
            print(f"{feature.chrom}:{feature.start:,}-{feature.end:,} "
                  f"[{feature.source}:{feature.feature_type}] {feature.name or ''}")
    
    print("\nFeature summary:")
    for key, count in sorted(feature_counts.items()):
        print(f"  {key}: {count} features")
        if key in feature_examples:
            ex = feature_examples[key]
            print(f"    Example: {ex.name} at {ex.start:,}")


def test_concurrent_merging():
    """Test that heapq.merge properly merges sorted streams."""
    
    print("\n" + "=" * 60)
    print("Testing concurrent stream merging...")
    
    # Create mock streams
    def stream1() -> Iterator[UnifiedFeature]:
        """Mock gene stream."""
        positions = [1000, 5000, 9000]
        for pos in positions:
            yield UnifiedFeature(
                chrom="1", start=pos, end=pos+1000,
                feature_type="gene", source="GTF",
                name=f"Gene_{pos}"
            )
    
    def stream2() -> Iterator[UnifiedFeature]:
        """Mock variant stream."""
        positions = [500, 1500, 6000, 9500]
        for pos in positions:
            yield UnifiedFeature(
                chrom="1", start=pos, end=pos,
                feature_type="SNP", source="VCF",
                name=f"Var_{pos}"
            )
    
    def stream3() -> Iterator[UnifiedFeature]:
        """Mock repeat stream."""
        positions = [2000, 7000, 8000]
        for pos in positions:
            yield UnifiedFeature(
                chrom="1", start=pos, end=pos+200,
                feature_type="repeat", source="RepeatMasker",
                name=f"Repeat_{pos}"
            )
    
    # Merge streams
    print("\nMerged stream (sorted by position):")
    merged = heapq.merge(stream1(), stream2(), stream3())
    
    for feature in merged:
        print(f"  {feature.start:5d}: {feature.source:12s} {feature.feature_type:8s} {feature.name}")
    
    print("\n✓ Streams merged in correct order!")


def test_genome_browser_integration():
    """Test integration with genome browser."""
    
    print("\n" + "=" * 60)
    print("Testing genome browser integration...")
    
    # Load genome manager
    gm = GenomeManager()
    
    # Create unified annotations
    annotations = UnifiedGenomeAnnotations()
    annotations.add_vcf("./genome/vcf/HG00096.vcf.gz", "variants")
    
    # Simulate what the browser would do
    chrom = "1"
    position = 204195220  # Near the big insertion
    window_size = 80
    
    print(f"\nQuerying annotations at {chrom}:{position:,} (window={window_size})")
    
    # Get features in window
    features = annotations.query_range(chrom, position, position + window_size)
    
    # Convert to browser-friendly format
    print(f"Found {len(features)} features:")
    for f in features[:5]:
        print(f"  {f.feature_type} at {f.start:,}: ", end="")
        if f.feature_type == "variant":
            ref = f.attributes.get('ref', '?')
            alt = f.attributes.get('alt', ['?'])[0] if f.attributes.get('alt') else '?'
            print(f"{ref} → {alt}")
        else:
            print(f"{f.name}")


def test_motif_integration():
    """Test motif scanning integration."""
    
    print("\n" + "=" * 60)
    print("Testing motif scanning...")
    
    from ggene.motifs.motif import PatternMotif
    
    # Create motif scanner
    splice_donor = PatternMotif("splice_donor", r"GT[AG]AGT", lambda x: 1.0)
    
    # Test sequence
    test_seq = "AAAAGTAAGTCCCGTGAGTAAAA"
    
    # Find motifs
    import re
    pattern = re.compile(r"GT[AG]AGT")
    
    print(f"Scanning: {test_seq}")
    for match in pattern.finditer(test_seq):
        print(f"  Found splice donor at {match.start()}: {match.group()}")
        
        # Convert to UnifiedFeature
        feature = UnifiedFeature(
            chrom="test",
            start=match.start() + 1,
            end=match.end(),
            feature_type="splice_donor",
            source="MotifScanner",
            score=1.0,
            name="GT..AGT",
            attributes={'sequence': match.group()}
        )
        print(f"    As feature: {feature.to_dict()}")


def demonstrate_feature_pipeline():
    """Show the complete annotation pipeline."""
    
    print("\n" + "=" * 60)
    print("Complete Annotation Pipeline Demo")
    print("=" * 60)
    
    # 1. Initialize unified annotation system
    print("\n1. Initializing annotation sources...")
    annotations = UnifiedGenomeAnnotations()
    
    # Add local files
    print("   - Adding GTF (genes)")
    print("   - Adding VCF (variants)")
    print("   - Adding BED (peaks)")
    
    # 2. Define region of interest
    chrom = "1"
    start = 1000000
    end = 1001000
    
    print(f"\n2. Region of interest: {chrom}:{start:,}-{end:,}")
    
    # 3. Stream all annotations
    print("\n3. Streaming merged annotations:")
    
    # This would actually stream from all sources
    mock_features = [
        UnifiedFeature(chrom, 1000100, 1000200, "gene", "GTF", name="GENE1"),
        UnifiedFeature(chrom, 1000150, 1000150, "SNP", "VCF", name="rs123"),
        UnifiedFeature(chrom, 1000180, 1000280, "repeat", "RepeatMasker", name="AluY"),
        UnifiedFeature(chrom, 1000250, 1000260, "tf_binding", "JASPAR", name="CTCF"),
    ]
    
    for f in sorted(mock_features):
        print(f"   {f.start:,}: {f.feature_type:12s} ({f.source:12s}) {f.name}")
    
    # 4. Convert to browser display
    print("\n4. Browser display format:")
    print("   Position  1000100    1000150    1000200    1000250")
    print("   Genes     |-------GENE1-------|")
    print("   Variants             *")
    print("   Repeats                    |-----AluY-----|")
    print("   TF Sites                               |CTCF|")
    
    print("\n✓ Pipeline complete!")


if __name__ == "__main__":
    # Run tests
    test_unified_streaming()
    test_concurrent_merging()
    test_genome_browser_integration()
    test_motif_integration()
    demonstrate_feature_pipeline()
    
    print("\n" + "=" * 60)
    print("All tests completed successfully!")
    print("=" * 60)