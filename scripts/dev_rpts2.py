
from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np
import random

from ggene import DATA_DIR
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align, repeats
from ggene.seqs.repeats import Repeat, RepeatLibrary, Motif, MotifNetwork
from ggene.seqs.vocab import VOCAB


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
            
            print(f"***** {rpti.motif}->{rptj.motif} HD: {ham} *****")
    
    return lib

def load_motifs(chr = "1"):
    
    motifs = []
    with open(f'./data/motifs_chr{chr}.txt', "r") as f:
        for line in f:
            motifs.append(line.strip())
    
    motifs = list(sorted(motifs, key = lambda k: bio.get_seq_index_abs(k)))
    
    return motifs

def build_motif_network(motifs, max_len = 6):
    
    curr_len = 2
    
    mtfs = MotifNetwork()
    
    for i, motif in enumerate(motifs):
        new_motif = Motif(motif)
        mtfs.add_to_tree(new_motif)
        
        if new_motif.motif_len > curr_len:
            curr_len += 1
        
        if curr_len > max_len:
            break
    
    mtfs.join_peers()
    return mtfs

def test_motif_network(gm, rpts):
    test_chr = "1"
    max_rpts = 1e4
    lib = test_rpt_lib(gm, rpts, test_chr = test_chr, max_num=max_rpts, topk = 12)
    
    with open(f"./data/motifs_chr{test_chr}.txt","w") as f:
        for mtf in lib.motifs:
            f.write(mtf + "\n")
    
    motifs = load_motifs(chr = test_chr)
    mtfs = build_motif_network(motifs)

    mtfs.print()
    mtfs.print_stats()
    
    return mtfs

def verify_rpt_rotation(motifa, motifb):
    
    max_err = len(motifb)
    
    rmotifbs = [rotate(motifb, r) for r in range(len(motifb))]
    
    res = repeats.hamming_distance_re(motifa, rmotifbs, max_err = max_err, min_err = 0)
    errs = [r[1] for r in res]
    # sctxt = draw.scalar_to_text_nb(errs, bit_depth = 8)
    
    # print(motifa, sctxt[0])
    for mtfb, (r, s, i, b )in zip(rmotifbs,res):
        print(mtfb, s+i+b, s, i, b)
    
    return errs

def verify_rpt_corotation(motifa, motifb):
    
    max_err = len(motifb)
    
    rcmotifb = rotate(bio.reverse_complement(motifb), r=1)
    rmotifbs = [rotate(rcmotifb, r) for r in range(len(motifb))]
    rcmotifbs = [rotate(rcmotifb, -r) for r in range(len(motifb))]
    
    res = repeats.hamming_distance_re(motifa, rmotifbs, max_err = max_err, min_err = 0)
    rcres = repeats.hamming_distance_re(motifa, rcmotifbs, max_err = max_err, min_err = 0)
    errs = [r[1] for r in res]
    rcerrs = [r[1] for r in rcres]
    
    for mtfb, (r, s, i, b), (rr,ss,ii,bb) in zip(rmotifbs, res, rcres):
        print(mtfb, s+i+b, ss+ii+bb)
    
    return errs, rcerrs

def rotate(motif, r = 1):
    r = r % len(motif)
    return motif[r:] + motif[:r]

def rpt_rotation_mc(num_rounds, len_motifa, len_motifb):
    
    for n in range(num_rounds):
        
        mtfa = "".join(["ATGC"[random.randint(0,3)] for n in range(len_motifa)])
        ins_pos = random.randint(0, len_motifa)
        mtfb = mtfa[:ins_pos] + random.choice("ATGC") + mtfa[ins_pos:]
        
        mtfa_min, _, _, _ = bio.get_least_seq_index(mtfa)
        mtfb_min, _, _, _ = bio.get_least_seq_index(mtfb)
        rcmtfb = bio.reverse_complement(mtfb_min)
        # rcmtfb_min = rotate(rcmtfb, r = 1)
        rcmtfb_min = rcmtfb
        
        errs, rcerrs = verify_rpt_corotation(mtfa_min, mtfb_min)
        revdiff_errs = [r-rc for r,rc in zip(errs, rcerrs)]
        
        # errs, rcerrs = verify_rpt_rotation(mtfa_min, mtfb_min)
        # rcerrs = verify_rpt_rotation(mtfa_min, rcmtfb_min)
        # sum_errs = [r+rc for r,rc in zip(errs, rcerrs)]
        # revdiff_errs = [r-rc for r,rc in zip(errs, reversed(rcerrs))]
        
        argmin_errs = get_all_argmin(errs)
        argmin_revdiff = get_all_argmin(revdiff_errs)
        
        print(mtfa)
        print(mtfb)
        ScalarPlot(errs, bit_depth = 8, minval = 0, maxval = len(mtfb_min), fg_color = 59).show()
        ScalarPlot(rcerrs, bit_depth = 8, minval = 0, maxval = len(mtfb_min), flip = True, fg_color = 131).show()
        print()
        ScalarPlot(revdiff_errs, bit_depth = 16, center = 0, mode = "mid").show()
        print()
        print(f"argmin: {argmin_errs}, argmin revdiff: {argmin_revdiff}")
        print()

def get_all_argmin(vals):
    minval = min(vals)
    argmins = []
    for i, v in enumerate(vals):
        if v==minval:
            argmins.append(i)
    return argmins

def main():
    
    gm, rpts = load_genome()
    
    # mtfs = test_motif_network(gm, rpts)
    
    # while True:
        
    #     mtfa = random.choice(mtfs.motifs)
    #     if len(mtfa.children) < 1:
    #         continue
        
    #     mtfbs = list(mtfa.children)
    #     for mtfb in mtfbs:
    #         verify_rpt_rotation(str(mtfa), str(mtfb))
    #         print()
    #     input()
    
    num_rounds = 1
    len_motifa = 9
    len_motifb = 6
    rpt_rotation_mc(num_rounds, len_motifa, len_motifb)
    
    
    # seq = "GGAAAAGGT"
    # rcseq = bio.reverse_complement(seq)
    # rcseq_r = rotate(rcseq, 1)
    # print(seq)
    # print(rcseq)
    # print(rcseq_r)
    
    # ScalarPlot([-0.5,-0.25,0, 0.25,1], )
    
    
    """
    
    
    
    
    
    
    """
    
    pass
    
    




if __name__=="__main__":
    main()

