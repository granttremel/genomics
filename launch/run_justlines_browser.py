        
from typing import Tuple, Optional
import argparse
import random

from ggene.draw import scalar_plot

from ggene.seqs import lambdas

from ggene.database.genome_manager import GenomeManager
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.colors import FColors

from ggene.browser.builds.justlines_browser import JustLinesBrowser, JustLinesBrowserState

def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")

def main():
    
    
    parser = argparse.ArgumentParser(description='Interactive Genome Browser')
    parser.add_argument('--chrom', '-c', type=str, default='1',
                       help='Chromosome to browse (default: 1)')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--window', '-w', type=int, default=240, dest = "window_size",
                       help='Starting position (default: 1000000)')
    parser.add_argument('--stride', '-s', type=int, default=20,
                       help='stride in base pairs (default: 80)')
    parser.add_argument('--debug', '-d', action="store_true",
                       help='debug', dest = "debug")
    parser.add_argument('--random', '-r', action="store_true",
                       help='random', dest = "random")
    
    args = parser.parse_args()
    
    if args.random:
        args.chrom, args.position = GenomeManager.get_random_location(margin = 10e6)
    
    # window_size = 256
    max_disp = 256
    display_width = max_disp
    display_height = 32
    
    brws = JustLinesBrowser(display_width = display_width, display_height = display_height, **args.__dict__)
    brws.start(display_width = display_width, **args.__dict__)
    
    # write_view(brws._rendered_view, "with_align")
    
    brws.to_yaml("./data/browser_data/test_browser_1.yaml")
    
    
    
    
if __name__=="__main__":
    main()

