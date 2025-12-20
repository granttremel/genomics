#!/usr/bin/env python3
"""Test variant display in genome browser."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene import get_paths
DEFAULT_VCF_PATH, DEFAULT_GTF_PATH, DEFAULT_FASTA_PATH, DEFAULT_LIBRARY = get_paths()
from ggene.database.genomemanager import GenomeManager
from ggene.database.genome_iterator import UGenomeIterator
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

def test_variant_features():
    """Test that variants are properly retrieved and displayed."""
    
    print("Testing Variant Feature Display")
    print("=" * 60)
    
    # Initialize GenomeManager
    gm = GenomeManager()
    
    # Test region with known variants
    chrom = 1
    start = 204195200
    end = start + 100
    
    print(f"\nTesting region chr{chrom}:{start:,}-{end:,}")
    print("-" * 60)
    
    # Create iterator
    iterator = UGenomeIterator(
        gm, chrom, start, end,
        window_size=50,
        integrate_variants=True,
        track_features=True
    )
    
    # Get a window
    window = iterator.get_window_at(start)
    
    print(f"\nWindow information:")
    print(f"  Reference sequence: {window.ref_seq[:50]}...")
    print(f"  Alternate sequence: {window.alt_seq[:50]}...")
    print(f"  Has gaps: {window.has_gaps}")
    
    # Check variant deltas (coordinate tracking)
    if window.variant_deltas:
        print(f"\n  Variant deltas (for coordinates): {len(window.variant_deltas)}")
        for pos, delta in window.variant_deltas[:3]:
            print(f"    Position {pos:,}: delta={delta:+d}")
    else:
        print(f"\n  No variant deltas found")
    
    # Check variant features (actual variant objects)
    if window.variant_features:
        print(f"\n  Variant features (for display): {len(window.variant_features)}")
        for var in window.variant_features[:3]:
            attrs = var.attributes
            ref = attrs.get('ref', '?')
            alt = attrs.get('alt', ['?'])[0] if isinstance(attrs.get('alt'), list) else attrs.get('alt', '?')
            qual = attrs.get('qual', 0)
            print(f"    {var.start:,}: {ref} → {alt} (Q={qual})")
    else:
        print(f"\n  No variant features found")
    
    # Also check if variants come through as regular features
    variant_in_features = [f for f in window.features if f.feature_type == 'variant']
    if variant_in_features:
        print(f"\n  Variants in feature list: {len(variant_in_features)}")
        for var in variant_in_features[:3]:
            print(f"    {var.start:,}: {var.feature_type}")
    
    # Test direct annotation query
    print("\n" + "-" * 60)
    print("Testing direct annotation query...")
    
    annotations = gm.get_all_annotations(str(chrom), start, end, include_motifs=False)
    variants_in_annotations = [a for a in annotations if a.feature_type == 'variant']
    
    print(f"  Total annotations: {len(annotations)}")
    print(f"  Variant annotations: {len(variants_in_annotations)}")
    
    if variants_in_annotations:
        for var in variants_in_annotations[:3]:
            attrs = var.attributes
            ref = attrs.get('ref', '?')
            alt = attrs.get('alt', ['?'])[0] if isinstance(attrs.get('alt'), list) else attrs.get('alt', '?')
            print(f"    {var.start:,}: {ref} → {alt}")
    
    # Test VCF stream directly
    print("\n" + "-" * 60)
    print("Testing VCF stream directly...")
    
    if hasattr(gm.annotations, 'streams') and 'variants' in gm.annotations.streams:
        vcf_stream = gm.annotations.streams['variants']
        print(f"  VCF stream available: {vcf_stream}")
        
        count = 0
        for var in vcf_stream.stream(str(chrom), start, end):
            count += 1
            if count <= 3:
                attrs = var.attributes
                ref = attrs.get('ref', '?')
                alt = attrs.get('alt', ['?'])[0] if isinstance(attrs.get('alt'), list) else attrs.get('alt', '?')
                print(f"    {var.start:,}: {ref} → {alt}")
        
        print(f"  Total variants from stream: {count}")


if __name__ == "__main__":
    test_variant_features()