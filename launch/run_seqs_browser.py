        
from typing import Tuple, Optional
import argparse
import random

from ggene.draw import scalar_plot

from ggene.seqs import lambdas

from ggene.database.genome_manager import GenomeManager
from ggene.draw.colors import visible_len, visible_slice
from ggene.display.artists import *
from ggene.display.colors import FColors

from ggene.browser.builds.seqs_browser import SeqsBrowser, SeqsBrowserState

def write_view(view, fname):
    
    with open(f"./data/browser_views/{fname}.txt", "w+") as f:    
        for line in view:
            f.write(FColors.scrub_codes(line) + "\n")

def load_genome():
    gm = GenomeManager()
    return gm


def get_repeat_stats():
    gm = load_genome()
    
    
    all_tes = {}
    
    for chr in gm.iter_chromes():
        
        tes = {}
        print(f"starting chrome {chr}")
        
        for feat in gm.annotations.stream_by_types("dfam_hit", chr, start = 10e6):
            
            if not feat.name in tes:
                tes[feat.name] = 0
            
            tes[feat.name] += 1
            
            # if len(tes)> 10:
            #     break
            
        plot_repeat_stats(tes)
        input(f"that was chr{chr}. enter to keep goin")
        
        all_tes = merge_te_dicts(all_tes, tes)
    
    print("all te stats!!")
    plot_repeat_stats(all_tes)

def merge_te_dicts(te1, te2):
    
    all_tes = list(set(list(te1.keys()) + list(te2.keys())))
    te_out = {te:0 for te in all_tes}
    
    for te, ct in te1.items():
        te_out[te] += ct
    for te, ct in te2.items():
        te_out[te] += ct
    
    return te_out

def plot_repeat_stats(tes):

    tes_srt = sorted(tes.items(), key = lambda k:-k[1])
    te_names = [te for te, ct in tes_srt]
    te_cts = [ct for te, ct in tes_srt]
    
    max_ct = te_cts[0]
    
    bits = 64
    sctxt = scalar_plot.scalar_to_text_nbh(te_cts, minval = 0, bit_depth = 8*bits, fg_color = 130, bg_color = 234)
    
    nr = 0
    for (lbl, ct), row in zip(tes_srt, sctxt):
        
        print("{:<16}{} ({})".format(lbl, row, ct))
        nr += 1
        
        if nr%256 == 0:
            input("enter for more")
            
        if ct < max_ct / bits:
            break
    
def main():
    
    parser = argparse.ArgumentParser(description='Interactive Genome Browser')
    parser.add_argument('--chrom', '-c', type=str, default='1',
                       help='Chromosome to browse (default: 1)')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--window', '-w', type=int, default=240, dest = "window_size",
                       help='Starting position (default: 1000000)')
    parser.add_argument('--stride', '-s', type=int, default=20,
                       help='Window size in base pairs (default: 80)')
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
    
    brws = SeqsBrowser(display_width = display_width, display_height = display_height, **args.__dict__)
    brws.start(display_width = display_width, **args.__dict__)
    
    # write_view(brws._rendered_view, "with_align")
    
    brws.to_yaml("./data/browser_data/test_browser_1.yaml")
    
    
    
    
if __name__=="__main__":
    main()

