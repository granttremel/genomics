#!/usr/bin/env python3
"""
Example usage of the GenomeIterator for genomic analysis with variant integration.
"""

import logging
import sys
from ggene.database.genome_manager import GenomeManager
from ggene.genome_iterator import GenomeIterator, FeatureExtractor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def example_basic_iteration():
    """Example 1: Basic iteration over genomic positions."""
    print("\n=== Example 1: Basic Genome Iteration ===")
    
    # Initialize GenomeManager
    gm = GenomeManager()
    
    # Create an iterator for chromosome 1, starting at position 1000000
    iterator = GenomeIterator(
        gm, 
        chrom=1,
        start=1000000,
        end=1000100,  # Iterate 100 bases
        stride=1,
        integrate_variants=True,
        track_features=True
    )
    
    # Iterate and collect information
    variant_positions = []
    feature_positions = []
    
    for position in iterator:
        # Check if this position has a variant
        if position.has_variant:
            variant_positions.append(position.position)
            print(f"Variant at position {position.position}: "
                  f"{position.reference_base} -> {position.personal_base}")
        
        # Check if this position has features
        if position.features:
            feature_positions.append(position.position)
            feature_types = ', '.join(position.feature_types)
            print(f"Features at position {position.position}: {feature_types}")
    
    print(f"\nFound {len(variant_positions)} positions with variants")
    print(f"Found {len(feature_positions)} positions with features")


def example_window_iteration():
    """Example 2: Iterate using windows for efficient sequence extraction."""
    print("\n=== Example 2: Window-based Iteration ===")
    
    gm = GenomeManager()
    
    # Create iterator with 100bp windows
    iterator = GenomeIterator(
        gm,
        chrom=1,
        start=1000000,
        end=1001000,
        stride=100,  # Move 100bp at a time
        window_size=100,  # Get 100bp windows
        integrate_variants=True
    )
    
    windows_with_variants = 0
    
    for ref_seq, personal_seq, features in iterator:
        # Check if this window has variants
        if ref_seq != personal_seq:
            windows_with_variants += 1
            num_differences = sum(1 for r, p in zip(ref_seq, personal_seq) if r != p)
            print(f"Window at position {iterator.current_pos - iterator.stride}: "
                  f"{num_differences} differences from reference")
    
    print(f"\n{windows_with_variants} windows contain variants")


def example_feature_extraction():
    """Example 3: Extract sequences around specific features."""
    print("\n=== Example 3: Feature Extraction for ML ===")
    
    gm = GenomeManager()
    
    # Find and assemble a gene (example with BRCA1 or any gene you have)
    gene_name = "KISS1"  # You can change this to any gene in your dataset
    
    # Try to load or assemble the gene
    gene = gm.get_gene(gene_name, chrom=1)  # Adjust chromosome as needed
    
    if not gene:
        print(f"Could not find gene {gene_name}, using example coordinates")
        # Create a mock feature for demonstration
        from ggene.genome.features import Feature
        gene = Feature({
            'feature': 'gene',
            'chrom': 1,
            'start': 1000000,
            'end': 1010000,
            'gene_name': 'EXAMPLE_GENE'
        })
    
    # Create feature extractor
    extractor = FeatureExtractor(gm)
    
    # Extract sequences around the gene with 2kb flanking regions
    sequences = extractor.extract_around_features(
        [gene],
        upstream=2000,
        downstream=2000,
        integrate_variants=True
    )
    
    for seq_data in sequences:
        print(f"\nExtracted sequence for {seq_data['feature_id']}:")
        print(f"  Length: {len(seq_data['personal_sequence'])} bp")
        print(f"  GC content: {seq_data.get('gc_content', 0):.2%}")
        print(f"  Has variants: {seq_data['has_variants']}")
        print(f"  Number of variants: {seq_data.get('num_variants', 0)}")
    
    # Prepare for ML models
    if sequences:
        ml_data = extractor.prepare_for_ml(
            sequences,
            sequence_length=1000,  # Standard input size for many models
            one_hot_encode=True
        )
        
        print(f"\nPrepared ML data:")
        print(f"  Shape: {ml_data['sequences'].shape}")
        print(f"  Encoding: {ml_data['encoding']}")
        print(f"  Number of sequences: {len(ml_data['metadata'])}")


def example_motif_search():
    """Example 4: Search for motifs while iterating."""
    print("\n=== Example 4: Motif Search with Iteration ===")
    
    gm = GenomeManager()
    
    # Define a motif to search for (e.g., TATA box)
    motif = "TATAAA"
    
    # Create iterator
    iterator = GenomeIterator(
        gm,
        chrom=1,
        start=1000000,
        end=1100000,
        stride=1,
        window_size=len(motif),
        integrate_variants=True,
        track_features=True
    )
    
    motif_positions = []
    
    for ref_seq, personal_seq, features in iterator:
        # Check if motif is present in personal sequence
        if personal_seq.upper() == motif:
            pos = iterator.current_pos - iterator.stride
            motif_positions.append(pos)
            
            # Check what features overlap this motif
            feature_types = set()
            for f in features:
                if 'feature' in f:
                    feature_types.add(f['feature'])
            
            print(f"Found motif at position {pos}")
            if feature_types:
                print(f"  Overlapping features: {', '.join(feature_types)}")
    
    print(f"\nFound {len(motif_positions)} occurrences of motif '{motif}'")


def example_variant_density():
    """Example 5: Calculate variant density across regions."""
    print("\n=== Example 5: Variant Density Analysis ===")
    
    gm = GenomeManager()
    
    # Analyze variant density in 10kb windows
    window_size = 10000
    iterator = GenomeIterator(
        gm,
        chrom=1,
        start=1000000,
        end=2000000,  # Analyze 1Mb region
        stride=window_size,
        window_size=window_size,
        integrate_variants=True
    )
    
    densities = []
    
    for ref_seq, personal_seq, _ in iterator:
        if ref_seq and personal_seq:
            # Count differences
            num_variants = sum(1 for r, p in zip(ref_seq, personal_seq) if r != p)
            density = num_variants / len(ref_seq) * 1000  # Variants per kb
            densities.append(density)
            
            pos = iterator.current_pos - iterator.stride
            print(f"Window {pos}-{pos + window_size}: "
                  f"{density:.2f} variants/kb")
    
    if densities:
        import numpy as np
        print(f"\nVariant density statistics:")
        print(f"  Mean: {np.mean(densities):.2f} variants/kb")
        print(f"  Median: {np.median(densities):.2f} variants/kb")
        print(f"  Max: {np.max(densities):.2f} variants/kb")
        print(f"  Min: {np.min(densities):.2f} variants/kb")


def example_extract_feature_sequence():
    """Example 6: Extract a complete feature sequence."""
    print("\n=== Example 6: Extract Complete Feature Sequence ===")
    
    gm = GenomeManager()
    
    # Look for an exon and extract its sequence
    iterator = GenomeIterator(
        gm,
        chrom=1,
        start=1000000,
        end=2000000,
        stride=1,
        integrate_variants=True,
        track_features=True
    )
    
    # Scan for next exon
    exon_position = iterator.scan_for_feature('exon', max_distance=100000)
    
    if exon_position:
        print(f"Found exon at position {exon_position}")
        
        # Jump to just before the exon
        iterator.jump_to(exon_position - 10)
        
        # Extract the exon sequence
        exon_seq = iterator.extract_feature_sequence(
            'exon',
            upstream=10,
            downstream=10
        )
        
        if exon_seq:
            print(f"Extracted exon sequence:")
            print(f"  Length: {len(exon_seq)} bp")
            print(f"  First 50 bases: {exon_seq[:50]}...")
            print(f"  Last 50 bases: ...{exon_seq[-50:]}")
            
            # Calculate GC content
            gc_count = exon_seq.upper().count('G') + exon_seq.upper().count('C')
            gc_content = gc_count / len(exon_seq) if len(exon_seq) > 0 else 0
            print(f"  GC content: {gc_content:.2%}")
    else:
        print("No exon found in the search region")


def main():
    """Run all examples."""
    
    examples = [
        ("Basic Iteration", example_basic_iteration),
        ("Window Iteration", example_window_iteration),
        ("Feature Extraction", example_feature_extraction),
        ("Motif Search", example_motif_search),
        ("Variant Density", example_variant_density),
        ("Extract Feature Sequence", example_extract_feature_sequence)
    ]
    
    print("=" * 60)
    print("GENOME ITERATOR EXAMPLES")
    print("=" * 60)
    
    for name, func in examples:
        try:
            func()
        except Exception as e:
            logger.error(f"Error in {name}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    # You can run specific examples or all of them
    if len(sys.argv) > 1:
        example_num = int(sys.argv[1])
        examples = {
            1: example_basic_iteration,
            2: example_window_iteration,
            3: example_feature_extraction,
            4: example_motif_search,
            5: example_variant_density,
            6: example_extract_feature_sequence
        }
        if example_num in examples:
            examples[example_num]()
        else:
            print(f"Unknown example number: {example_num}")
            print(f"Available examples: 1-{len(examples)}")
    else:
        main()