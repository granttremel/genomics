#!/usr/bin/env python3
"""
Launch the interactive genome browser.
"""

import sys
import argparse
import logging
from ggene import other_paths
from ggene.database.genomemanager import GenomeManager
from ggene.browser.genome_browser import browse_genome as browse_genome_v1
from ggene.browser.genome_browser_v2 import browse_genome as browse_genome_v2

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
    parser.add_argument('--random', '-r', action="store_true",
                       help='start at random position')
    parser.add_argument('--debug', '-d', action="store_true",
                       help='debug')
    parser.add_argument('--old', '-o', action="store_true",
                       help='do version 1')
    
    parser.add_argument('--artist-type', type=str, default = "line", help = "normal, seq, line, zoom ... ?")
    
    args = parser.parse_args()
    
    print(f"args: {args}")
    
    print("Loading genome data...")
    gm = GenomeManager(**other_paths)
    
    # If gene is specified, find its position
    chrom = ""
    pos = -1
    if args.position:
        pos = args.position
        chrom = args.chrom
    
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
    elif args.random:
        
        import random
        chrom = str(random.choice(range(23)))
        max_ind = gm.gene_map.max_indices.get(chrom)
        pos = int(random.random() * (max_ind - 1e6) + 1e6)
    
    artist_defaults = {"show_ruler":True}
    
    # gui = int(args.gui[0])==1
    # Start the browser
    
    window_sz = args.window
    if args.artist_type == "line":
        window_sz *= 4
      
    if args.old:
        browse_genome_v1(gm, chrom, pos, args.window, debug = args.debug, use_gui=False, artist_type = args.artist_type, **artist_defaults)
    else:
        browse_genome_v2(gm, chrom, pos, args.window, debug = args.debug, use_gui=False)


if __name__ == "__main__":
    main()