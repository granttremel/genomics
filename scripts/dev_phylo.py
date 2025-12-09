
from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np
import random

from ggene import DATA_DIR, other_paths
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align, repeats, compare
from ggene.seqs.repeats import Repeat, RepeatLibrary, Motif, MotifNetwork
from ggene.seqs.vocab import VOCAB

def load_genome():
    gm = GenomeManager(**other_paths)
    rpts = gm.annotations.streams.get("repeats")
    return gm, rpts

def get_repeats(gm, rpts, chr, name = "", type = "", skip = 1):
    
    out_rpts = []
    
    num_rpts = 0
    for i, rpt in enumerate(rpts.stream(chr, start = 1e6)):
        
        if name and rpt.get("name") != name:
            continue
        if type and rpt.get("type") != type:
            continue
        
        num_rpts += 1
        
        if num_rpts % skip != 0:
            continue
        
        loc = repeats.Loc.from_feature(rpt)
        
        new_rpt = repeats.RepeatSeq.get_instance(gm, rpt, context_sz = 128)
        out_rpts.append(new_rpt)
    
    return out_rpts

def organize_repeats(gm, rpts, chr, name = "", type = "", skip = 1, max_num = None):
    
    rpts = get_repeats(gm, rpts, chr, name=name, type=type, skip=skip)
    rpts = list(sorted(rpts, key = lambda rpt:-len(rpt.repeat)))
    
    print(f"identified {len(rpts)} repeats")
    if max_num:
        rpts = rpts[:max_num]
    
    return rpts

def build_phylo(seqs, refa = "", refb = ""):
    if not refa:
        refa = seqs.pop(0)
    if not refb:
        refb = seqs.pop(-1)
    
    phylo = compare.Phylo(refa, refb)
    
    for seq in seqs:
        phylo.insert(seq)
        
    return phylo

def print_phylo(phn, tab = 0, max_len = 32):
    
    if isinstance(phn, compare.Phylo):
        phn = phn.root
    
    if isinstance(phn, compare.PhyloLeaf):
        phn.print()
        
    else:
        phn.print()
        print()
        for ch in [phn.childa, phn.childb]:
            print_phylo(ch, tab = tab+1, max_len=max_len)

def main():
    
    gm, rpts = load_genome()    
    test_chr = "19"
    # rpt_name = "LTR38C"
    rpt_name = "AluJr"
    rpt_type = ""
    
    my_rpts = organize_repeats(gm, rpts, test_chr, name=rpt_name, type=rpt_type, skip = 50, max_num = 16)
    
    seqs = [rpt.repeat for rpt in my_rpts]

    phylo = build_phylo(seqs)
    phylo.print(show_length = 64)
    phylo.recompute_consensus()
    phylo.print(show_length = 64)
    
    
    
if __name__=="__main__":
    main()


