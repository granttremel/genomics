#!/usr/bin/env python3
"""
Launch the interactive genome browser.
"""

import sys
import argparse
import logging
from ggene.genomemanager import GenomeManager
from ggene.genome_browser import browse_genome

logging.basicConfig(level=logging.WARNING)  # Keep it quiet for browsing
import os
os.environ['QT_AUTO_SCREEN_SCALE_FACTOR'] = '5'

def main():
    parser = argparse.ArgumentParser(description='Interactive Genome Browser')
    parser.add_argument('--chrom', '-c', type=str, default='1',
                       help='Chromosome to browse (default: 1)')
    parser.add_argument('--gui', type=int, nargs=1, default=1,
                       help='Add flag to use gui, otherwise use terminal')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--window', '-w', type=int, default=240,
                       help='Window size in base pairs (default: 80)')
    parser.add_argument('--gene', '-g', type=str,
                       help='Jump to a specific gene (overrides position)')
    parser.add_argument('--debug', '-d', action="store_true",
                       help='debug')
    
    args = parser.parse_args()
    
    print("Loading genome data...")
    gm = GenomeManager()
    
    # If gene is specified, find its position
    if args.gene:
        gene_name = args.gene.upper()
        chrom, gene_info = gm.gene_map.find_gene(gene_name)
        
        if chrom and gene_info:
            gene_data = gene_info[0]
            args.chrom = chrom
            args.position = gene_data['start']
            print(f"Found {gene_name} on chromosome {chrom} at position {args.position:,}")
        else:
            print(f"Gene '{gene_name}' not found, using default position")
    # gui = int(args.gui[0])==1
    # Start the browser
    
    
    browse_genome(gm, args.chrom, args.position, args.window, debug = args.debug, use_gui=False)


if __name__ == "__main__":
    main()