#!/usr/bin/env python3
"""
Demo script for the interactive genome browser.
Shows how to use the browser to explore your personal genome.
"""

import logging
from ggene.database.genome_manager import GenomeManager

# Suppress debug output for cleaner display
logging.basicConfig(level=logging.WARNING)


def main():
    print("=" * 60)
    print("INTERACTIVE GENOME BROWSER DEMO")
    print("=" * 60)
    print()
    print("This browser lets you explore your personal genome")
    print("compared to the reference, with variant highlighting")
    print("and feature annotations.")
    print()
    print("Loading genome data...")
    
    # Initialize GenomeManager
    gm = GenomeManager()
    
    print("\nBrowser Controls:")
    print("-" * 40)
    print("→ / l       : Move forward")
    print("← / h       : Move backward")
    print("↑ / k       : Jump forward (10x)")
    print("↓ / j       : Jump backward (10x)")
    print("Space/Enter : Move forward")
    print("g           : Go to position/gene")
    print("f           : Toggle features")
    print("i           : Show position info")
    print("?           : Help")
    print("q           : Quit")
    print("-" * 40)
    
    # Example starting positions
    examples = [
        ("KISS1", None),     # Jump to KISS1 gene
        ("HTR2A", None),     # Jump to HTR2A gene
        (1, 1000000),        # Chromosome 1, position 1M
        (2, 50000000),       # Chromosome 2, position 50M
    ]
    
    print("\nExample starting points:")
    for i, (chrom_or_gene, pos) in enumerate(examples, 1):
        if pos:
            print(f"{i}. Chromosome {chrom_or_gene}, position {pos:,}")
        else:
            print(f"{i}. Gene {chrom_or_gene}")
    print(f"{len(examples)+1}. Custom position")
    print("0. Exit")
    
    try:
        choice = input("\nSelect option (0-{}): ".format(len(examples)+1))
        choice = int(choice)
        
        if choice == 0:
            print("Exiting...")
            return
        elif 1 <= choice <= len(examples):
            chrom_or_gene, pos = examples[choice-1]
            
            if pos:
                # Direct position
                print(f"\nStarting browser at chr{chrom_or_gene}:{pos:,}")
                gm.browse(chrom_or_gene, pos)
            else:
                # Gene name - find position
                gene_name = chrom_or_gene
                print(f"\nSearching for gene {gene_name}...")
                chrom, gene_info = gm.gene_map.find_gene(gene_name)
                
                if chrom and gene_info:
                    gene_data = gene_info[0]
                    start_pos = gene_data['start']
                    print(f"Found {gene_name} on chromosome {chrom} at position {start_pos:,}")
                    gm.browse(chrom, start_pos)
                else:
                    print(f"Gene {gene_name} not found")
        else:
            # Custom position
            try:
                chrom_input = input("Enter chromosome (e.g., 1, 2, X): ")
                pos_input = input("Enter position: ")
                pos = int(pos_input.replace(',', ''))
                
                print(f"\nStarting browser at chr{chrom_input}:{pos:,}")
                gm.browse(chrom_input, pos)
            except ValueError:
                print("Invalid input")
                
    except (KeyboardInterrupt, ValueError):
        print("\nExiting...")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()