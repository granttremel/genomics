
from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np
import random

from Bio import Align

from ggene import DATA_DIR
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align, repeats
from ggene.seqs.repeats import Repeat
from ggene.seqs.vocab import VOCAB

from Bio.Align import substitution_matrices


def load_genome():
    gm = GenomeManager()
    rpts = unified_stream.BEDStream(DATA_DIR / "repeatmasker/repeats.sorted.bed.gz", feature_type = "repeats")
    gm.annotations.add_source("repeats", rpts)
    return gm, rpts


def extract_repeats(gm, chr, repeat_name = "", repeat_type = "", motif = "", start = None, end = None):
    
    rpts = []
    for rpt in gm.annotations.streams.get("repeats").stream(chr, start=start, end=end):
        
        if repeat_name and rpt.name != repeat_name:
            continue
        if repeat_type and rpt.attributes.get("type") != repeat_type:
            continue
        if motif and rpt.attributes.get("motif") != motif:
            continue
        
        seq = gm.annotations.get_sequence(chr, rpt.start, rpt.end)
        rpt.attributes["seq"] = seq
        rpts.append(rpt)
    
    return rpts

def gather_repeats(gm:GenomeManager, rpt_proto, chrs = [], context = 128, length_range = 0):
    
    rpt_proto, ind, arg_min = bio.get_least_seq_index(rpt_proto)
    rpt_test = 2*rpt_proto
    rpt_len = len(rpt_proto)
    
    rpts = gm.annotations.streams.get("repeats")
    
    out_rpts = []
    
    if not chrs:
        chrs = [str(i) for i in range(24)] + ["X","Y"]
    
    for chr in chrs:
        testchr = str(chr)
        
        for rpt in rpts.stream(testchr, start = None):
            new_seq = ""
            mtf = rpt.attributes.get("motif","")
            if abs(rpt_len - len(mtf)) <= length_range:
                rmtf, ind, mtf_arg_min = bio.get_least_seq_index(mtf)
                # is_rc = False
                # if mtf in rpt_test:
                #     new_seq = gm.annotations.get_sequence(testchr, rpt.start - context, rpt.end + context)
                # elif bio.reverse_complement(mtf) in rpt_test:
                #     new_seq = bio.reverse_complement(gm.annotations.get_sequence(testchr, rpt.start - context, rpt.end + context))
                #     is_rc = True
                # else:
                #     continue
                # phase = len(mtf) - mtf_arg_min
                # phase = get_phase(rpt_proto, new_seq[context:len(new_seq) - context])
                # new_rpt = Repeat(rpt, mtf, rpt.chrom, rpt.start, rpt.end, is_rc, new_seq)
                # out_rpts.append(new_rpt)
    
    for rpt in out_rpts:
        for g in gm.annotations.stream_by_types(["gene"], rpt.chr, start = rpt.start, end = rpt.end):
            gn = g.get("info",{}).get("gene_name","")
            if gn:
                rpt.parent_gene_name = gn
                break
    
    return rpt_proto, out_rpts

def test_rpt_lib(gm, rpts, test_chr = "19", max_num = 1e5, topk = 24):
    
    lib = repeats.RepeatLibrary(motifs_only = True)
    
    nrpt = 0
    for rpt_feat in rpts.stream(test_chr, start = 1e6, end=None):
        inst = lib.get_instance(gm, rpt_feat, context_sz = 128)
        if inst is not None:
            nrpt += 1
        
        if nrpt > max_num:
            break
    
    tops = lib.get_top(topk = topk)
    
    for i in range(topk):
        for j in range(i, topk):
            
            
            rpti = tops[i]
            rptj = tops[j]
            
            if abs(len(rpti.motif) - len(rptj.motif)) > 1:
                continue
            
            ham = rpti.get_hamming_distance(rptj)
            # algn = align.align_sequences(rpti.motif, rptj.motif)[0]
            
            print(f"***** {rpti.motif}->{rptj.motif} HD: {ham} *****")
    
    return lib


def main():
    
    gm, rpts = load_genome()
    
    # test_rpt_align(double_second = False)
    
    max_rpts = 1e4
    lib = test_rpt_lib(gm, rpts, max_num=max_rpts, topk = 12)
    print(lib.motifs, len(lib.motifs), len(list(set(lib.motifs))))
    neighbs = lib.build_neighborhood()
    
    print(neighbs)
    
    # names = substitution_matrices.load()
    # # print(names)
    # for name in names:
    #     print(name)
    #     print(substitution_matrices.load(name))
        
    
    # seqa = "AAAT"
    # seqb = "AT"
    
    # test_ham(seqa, seqb)
    # test_ham(seqb, seqa)
        
    
    
    pass


if __name__=="__main__":
    main()

