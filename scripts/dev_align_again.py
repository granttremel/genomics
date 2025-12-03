
import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

from tabulate import tabulate

import numpy as np

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import bio, process, find, align


def load_genome():
    return GenomeManager()

def get_test_seqs(ind, with_context = False):
    """
    CCCCCTGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGG
    CCCC-TGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGGG
    
      TGGGGAGAAACAAGGGA AGGGAAGAAAAGGGACCGACCC
    TGTGGGGCATTT GGGGGTGTGGG  GAAAAA GACC AGAGTA
    
    TGGGGAAGCC  GGGCCGCTGTGAGCCTGGGAGA      CAGGAG ATGCTATCAGGAGGG
    TGG  AAATGTGGGGCATTTGGGGGTGTGGGGAAAAAGACCAG AGTATGGGA CAG
    
    CAGACCCCAGCGTGGGGAAGCCGGGCCGCTGTGAGCCTGGGAGACA
    AGCTCCCCAGGCTGG AAATGTGGGGCATTTGGGGGTGTGGGGAAAA
    
    TTTT GAAATGACCT AATTATCTAAGA
    TTTTTGAAAACAC TGAATT TGTAAAA
    
    
    """
    
    
    if ind == -1:
        seqa, seqb = "CCCCCTGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGG", "CCCCTGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGGG"
        if with_context:
            seqa = "ATA"
    elif ind == 0:
        seqa, seqb = "TGGGGAGAAACAAGGGAAGGGAAGAAAAGGGACCGACCC", "TGTGGGGCATTTGGGGGTGTGGGGAAAAAGACCAGAGTA"
        if with_context:
            seqa = "GGC" + seqa + "CTC"
            seqb = "AAA" + seqb + "TGG"
    
    elif ind == 1:
        seqa, seqb = "TGGGGAAGCCGGGCCGCTGTGAGCCTGGGAGACAGGAGATGCTATCAGGAGGG", "TGGAAATGTGGGGCATTTGGGGGTGTGGGGAAAAAGACCAGAGTATGGGACAG"
        if with_context:
            seqa = "GCG" + seqa + "TTG"
            seqb = "GGC" + seqb + "ACA"
    
    elif ind == 2:
        seqa, seqb = "CAGACCCCAGCGTGGGGAAGCCGGGCCGCTGTGAGCCTGGGAGACA", "AGCTCCCCAGGCTGGAAATGTGGGGCATTTGGGGGTGTGGGGAAAA"
    
    elif ind == 3:
        seqa, seqb = "TTTTGAAATGACCTAATTATCTAAGA", "TTTTTGAAAACACTGAATTTGTAAAA"
    
    return seqa, seqb


def plot_runs(seqa, seqb, max_err):
    runs, inds, shifts, errs = process.correlate_longest_subseq_err(seqa, seqb, fill = None, max_err = max_err)
    ScalarPlot(runs, add_range = True, minval = 0, xmin = 0, xmax = len(runs), ruler = True, num_labels = 2, ticks = 0, minor_ticks = len(runs)).show()
    
    tops, topdatas = process.extract_max_runs(seqa, seqb, runs, inds, shifts, 3)
    
    nr = 0
    for ssa, ssb in tops:
        print(f"run{nr}",ssa)
        print(f"run{nr}",ssb)
        print()
        nr += 1
    

def plot_corr(seqa, seqb):
    corr, rccorr = process.correlate(seqa, seqb, fill = None)  
    ScalarPlot(corr, add_range = True, minval = 0, xmin = 0, xmax = len(corr), ruler = True, num_labels = 2, ticks = 0, minor_ticks = len(corr)).show()    
    print(seqa)
    print(seqb)
    ScalarPlot(rccorr, add_range = True, minval = 0, xmin = 0, xmax = len(corr), ruler = True, num_labels = 2, ticks = 0, minor_ticks = len(corr)).show()
    print(seqa)
    print(bio.reverse_complement(seqb))
    print()
    return corr, rccorr

def test_extract_corr(seqa, seqb, topk):
    
    corr, rccorr = plot_corr(seqa, seqb)
    seqs_corr = process.extract_top_correlated(seqa, seqb, corr, topk)
    
    entrs = []
    
    for ssa, ssb in seqs_corr:
        cons = bio.merge(ssa, ssb)
        cons_entr = bio.consensus_entropy(cons)
        entrs.append(cons_entr)
        print(ssa)
        print(ssb)
        print(cons.replace("N", " "), f"entropy = {cons_entr:0.3f}")
        print()
    
    print("entropy vs rank")
    ScalarPlot(entrs, add_range = True, minval = 0).show()
    

def test_align(seqa, seqb):
    
    
    seqa, seqb, checkpts = align.get_align_checkpoints(seqa, seqb, err_ratio = 0, num_checkpoints = 3)
    print(checkpts)
    
    plot_runs(seqa, seqb, 20)
    plot_runs(seqa, bio.reverse_complement(seqb), 20)
    
    diffs = align.align_sequences_bad(seqa, seqb)
    seqafill, seqbfill = align.fill_aligned(seqa, seqb, diffs)
    print("aligned (?)")
    print(seqafill)
    print(seqbfill)
    print()
    
    
    pass

def main():
    
    # gm = load_genome()
    
    ind = 3
    seqa, seqb = get_test_seqs(ind)
    
    topk = 10
    # test_extract_corr(seqa, seqb, topk)
    
    test_align(seqa, seqb)
    test_align(seqb, seqa)
    
    
    
    

    


if __name__=="__main__":
    main()
