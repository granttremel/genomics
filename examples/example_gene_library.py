#!/usr/bin/env python3
"""Example script demonstrating gene library save/load functionality."""

from ggene.database.genome_manager import GenomeManager
from ggene.genome.features import Gene
import json
import os

def main():
    # Initialize GenomeManager
    gm = GenomeManager()
    
    # Example 1: Load a single gene from JSON
    print("=== Loading KISS1 gene from JSON ===")
    kiss1_gene = gm.load_gene("analysis/KISS1.json")
    if kiss1_gene:
        print(f"Loaded {kiss1_gene.name}")
        print(f"  Transcripts: {len(kiss1_gene.transcripts)}")
        print(f"  Exons: {len(kiss1_gene.exons)}")
        print(f"  Variants: {len(kiss1_gene.variants)}")
        
        # Test the linking
        for variant in kiss1_gene.variants[:3]:
            parents = variant.get_parent_ids()
            print(f"  Variant {variant.sfid} parents: {parents}")
    
    # Example 2: Get sequences with personal variants
    if kiss1_gene and kiss1_gene.exons:
        print("\n=== Comparing reference vs personal sequences ===")
        exon = list(kiss1_gene.exons.values())[0]
        comparison = gm.get_sequence_comparison(exon, upstream=10, downstream=10)
        
        if comparison:
            print(f"Exon {exon.exon_number}:")
            print(f"  Reference: {comparison['reference_sequence'][:50]}...")
            print(f"  Personal:  {comparison['personal_sequence'][:50]}...")
            print(f"  Variants in region: {comparison['num_variants']}")
            print(f"  Sequences differ: {comparison['sequences_differ']}")
    
    # Example 3: Save a gene back to JSON
    if kiss1_gene:
        print("\n=== Saving gene to new JSON file ===")
        output_file = "analysis/KISS1_reloaded.json"
        gm.save_gene(kiss1_gene, output_file, include_ordered_features=True)
        print(f"Saved to {output_file}")
    
    # Example 4: Load multiple genes from a directory
    print("\n=== Loading gene library ===")
    gene_library = gm.load_gene_library("analysis", "*.json")
    print(f"Loaded {len(gene_library)} genes:")
    for gene_name, gene in gene_library.items():
        print(f"  - {gene_name}: {len(gene.variants)} variants")
    
    # Example 5: Use get_or_load_gene for smart loading
    print("\n=== Smart gene loading ===")
    # This will load from cache since we already loaded it
    gene = gm.get_or_load_gene("KISS1")
    print(f"Got {gene.name} from cache")
    
    # This will try to find JSON file or assemble if chromosome provided
    # gene = gm.get_or_load_gene("BRCA1", chrom="17")

if __name__ == "__main__":
    main()