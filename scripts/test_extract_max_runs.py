



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

"""
TGCCTGCCGGCTGCTAATTTATTAACTCCCAGTGAATCATGTCCTGTGAAGGGACTGAATATTAGTGGCAATTTATGTTGATGATTTGTATTTTGAATAGT
TTGAATACATAGAACATTAAGCTTGTATACATTTTGAAAATAGTATTTTAATATTCTACTGTGTCATAGTTACAATGATTGGATATATGTTGAATTTATAT
GTACTTTGAGTTGTTCTATGTTTATGGTCTTTAGCATTCTAACGTGCAATTGTATATCTGTTAAGTCTTTTTTTTTTCGAGATTAGACTGATTTATTGAGG
CGTCTGTTTGATGCCACATTAAGTGGCCCAGGCTTTGTGTAGGGTTGAGGTTAAAGCAGGAAGAAGGGTGGTGAGAGGCGGGGCACCAGGGTTAGGTTGGA
ATACCTGGGGGTGCTCTGAGGCTCCCCAAGTTTCCCTGGTCTTGGCCGGCTGTGCTGCTGGCCTGGGCATCTGATGGGCTTGCAAGGGTGGTCCAGGGGCT
AGGGCAGGGACTTTGGAGTGACGCCGTTGGCTTTGAATCCAGACTCCTACACTTGGTAGCTGTGAAC
"""

def get_test_seqs():
    
    seqa = "AGTTTCCCTGGTCTTGGCCGGCTGTGCTGCTGGCCTGGGCATTGCAAAGTCTGATG"
    seqb = "AGTTTCCCTTATTCGGTCTTGGCCGGCTGTGCTGCTGGCCTGGTGAGCATCTGATG"
    return seqa, seqb


def get_run_metrics_fwd(seqa, seqb):
    
    buffer = 0
    topk = 1
    max_err = 5
    
    runs, inds, shifts, errs = process.correlate_longest_subseq_err(seqa, seqb, fill = 0, max_err = max_err)
    
    print("forward runs")
    ScalarPlot(runs, add_range = True, minval = 0).show_chunks()
    tops, topdatas = process.extract_max_runs(seqa, seqb, runs, inds, shifts, topk, buffer = buffer)
    draw.highlight_run(seqa, seqb, tops, topdatas, suppress = False)
    
    print("run hist")
    frhist, bins = np.histogram(runs, bins = max(runs), range = (0, max(runs)), density = True)
    ScalarPlot(frhist, add_range = True, minval = 0).show()
    
    cons = bio.merge(*tops[0])
    cons_entr = bio.consensus_entropy(cons)
    print(f"consensus: {cons}, entropy = {cons_entr:0.3f}")
    
def get_run_metrics_rev(seqa, seqb):
    
    buffer = 0
    topk = 1
    max_err = 5
    
    ##### reverse comp #####
    rc_seqb = bio.reverse_complement(seqb)
    rcruns, rcinds, rcshifts, rcerrs = process.correlate_longest_subseq_err(seqa, rc_seqb, max_err = max_err)
    
    print("rc runs")
    ScalarPlot(rcruns, add_range = True, minval = 0).show_chunks()
    rctops, rctopdatas = process.extract_max_runs(seqa, rc_seqb, rcruns, rcinds, rcshifts, topk, buffer = buffer)
    draw.highlight_run(seqa, rc_seqb, rctops, rctopdatas, suppress = False)
    
    print("rc run hist")
    rcrhist, bins = np.histogram(rcruns, bins = max(rcruns), range = (0, max(rcruns)), density = True)
    ScalarPlot(rcrhist, add_range = True, minval = 0).show()
    
    cons = bio.merge(*rctops[0])
    cons_entr = bio.consensus_entropy(cons)
    print(f"consensus: {cons}, entropy = {cons_entr:0.3f}")

def main():
    
    seqa, seqb = get_test_seqs()
    get_run_metrics_fwd(seqa, seqb)
    
    rc_seqb = bio.reverse_complement(seqb)
    get_run_metrics_rev(seqa, rc_seqb)
    
    



if __name__=="__main__":
    main()
