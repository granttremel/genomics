#!/usr/bin/env python3
"""
Test the enhanced genome browser features.
"""

import logging
from ggene.genomemanager import GenomeManager

# Keep it quiet
logging.basicConfig(level=logging.WARNING)


def test_browser_features():
    """Test the enhanced browser features."""
    print("=" * 60)
    print("ENHANCED GENOME BROWSER TEST")
    print("=" * 60)
    print()
    print("New Features:")
    print("-" * 40)
    print("1. Compressed feature labels (exon 1/3)")
    print("2. Start/stop codon highlighting in sequence")
    print("3. Amino acid translation with variant effects")
    print("4. Press 'a' to toggle amino acid display")
    print("-" * 40)
    print()
    
    print("Loading genome data...")
    gm = GenomeManager()
    
    # Find a good test position - look for a gene with known features
    test_cases = [
        ("Look for exons with variants", 1, 204500000),  # Near KISS1
        ("CDS region with potential variants", 1, 204510000),
        ("Start codon region", 1, 204500100),
    ]
    
    print("\nTest positions:")
    for i, (desc, chrom, pos) in enumerate(test_cases, 1):
        print(f"{i}. {desc} - chr{chrom}:{pos:,}")
    print("0. Exit")
    
    try:
        choice = input("\nSelect test case (0-3): ")
        choice = int(choice)
        
        if choice == 0:
            return
        elif 1 <= choice <= len(test_cases):
            desc, chrom, pos = test_cases[choice-1]
            print(f"\nLaunching browser: {desc}")
            print("\nRemember:")
            print("- Press 'a' to toggle amino acid translation")
            print("- Press 'f' to toggle feature display")
            print("- Look for compressed labels like 'exon 1/3' or 'txs(5)'")
            print("- Green = start codons, Red = stop codons in sequence")
            print()
            
            gm.browse(chrom, pos, window_size=90)
        else:
            print("Invalid choice")
            
    except (KeyboardInterrupt, ValueError):
        print("\nExiting...")


def find_coding_region_with_variants():
    """Helper to find a CDS region with variants for testing."""
    print("\nSearching for CDS regions with variants...")
    
    gm = GenomeManager()
    
    # Look through a region likely to have genes
    chrom = 1
    start = 1000000
    end = 2000000
    
    # Find CDS features
    cds_features = list(gm.gene_map.fetch(
        chrom, start, end, 
        features=('CDS',)
    ))
    
    if cds_features:
        print(f"Found {len(cds_features)} CDS features")
        
        # Check first few for variants
        for cds in cds_features[:10]:
            cds_start = cds.get('start', 0)
            cds_end = cds.get('end', 0)
            
            # Check for variants
            variant_count = 0
            try:
                for var in gm.vcf(gm._make_index(chrom, cds_start, cds_end)):
                    if var.QUAL > 20:
                        variant_count += 1
            except:
                pass
            
            if variant_count > 0:
                gene_name = cds.get('info', {}).get('gene_name', 'Unknown')
                print(f"\nFound CDS with {variant_count} variants:")
                print(f"  Gene: {gene_name}")
                print(f"  Position: chr{chrom}:{cds_start:,}-{cds_end:,}")
                print(f"\nTry: gm.browse({chrom}, {cds_start})")
                return chrom, cds_start
    
    print("No CDS with variants found in test region")
    return None, None


def main():
    """Run tests."""
    test_browser_features()
    
    # Optional: find a good test region
    print("\n" + "=" * 60)
    response = input("Search for CDS regions with variants? (y/n): ")
    if response.lower() == 'y':
        find_coding_region_with_variants()


if __name__ == "__main__":
    main()