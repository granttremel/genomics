from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import re
import numpy as np
import random

from ggene import DATA_DIR, other_paths
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align, repeats, compare, lambdas
from ggene.seqs.repeats import Repeat, RepeatLibrary, Motif, MotifNetwork
from ggene.seqs.vocab import VOCAB

def load_genome():
    gm = GenomeManager(**other_paths)
    rpts = gm.annotations.streams.get("repeats")
    # rpts = unified_stream.BEDStream(DATA_DIR / "repeatmasker/repeats.sorted.bed.gz", feature_type = "repeats")
    # gm.annotations.add_source("repeats", rpts)
    return gm, rpts

def get_repeat_types(gm, rpts, chr = "1", name_prefixes = []):
    
    rpt_names_pre = {pre:{} for pre in name_prefixes}
    rpt_names = {}
    rpt_types = {}
    for rpt in rpts.stream(chr,start=1e6):
        rpt_name = rpt.get("name")
        rpt_type = rpt.get("type",{})
        
        if rpt_type == "Simple_repeat":
            continue
        
        done = False
        for pre in name_prefixes:
            if rpt_name.startswith(pre):
                done = True
                suff = rpt_name.removeprefix(pre)
                if not suff in rpt_names_pre[pre]:
                    rpt_names_pre[pre][suff] = 0
                rpt_names_pre[pre][suff] += 1
        if not done:
            if not rpt_name in rpt_names:
                rpt_names[rpt_name]= 0
            rpt_names[rpt_name] += 1
        
        if not rpt_type in rpt_types:
            rpt_types[rpt_type]= 0
        rpt_types[rpt_type] += 1
    
    return rpt_names_pre, rpt_names, rpt_types

def display_repeats(gm, rpts):
    
    pres = [
        "Alu", "L1", "L", "U", "LTR", "MIR", "MLT", "MST", "SVA", "CR1", "HUERS",
        "Charlie", "Eulor", "Tigger", "Mam", "Kanga",
        "FLAM", "UCON", "MER", "HERV", "HERV", "ERV", "PABL", "THE1",
    ]
    
    pres = list(sorted(pres, key = lambda k:-len(k)))
    
    prefixed, names, types = get_repeat_types(gm, rpts, name_prefixes = pres)
    
    pf_srt = {pf:sorted(pfdict.items(), key = lambda k:-k[1]) for pf, pfdict in prefixed.items()}
    names_srt = sorted(names.items(), key = lambda k:-k[1])
    types_srt = sorted(types.items(), key = lambda k:-k[1])
    
    print("repeat names:")
    for n, v in names_srt:
        if v < 5:
            continue
        print(n,v)
    print()
    
    print("prefixed repeat names:")
    for pf, pfdict in pf_srt.items():
        print(pf, sum([v for k,v in pfdict]))
        print(", ".join([k for k,v in pfdict]))
        # for n, v in pfdict:
        #     if v < 5:
        #         continue
        #     print("  ", n,v)
        print()
    print()
    
    print("repeat types:")
    for n, v in types_srt:
        print(n,v)
    print()
    
    pass

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
        # print(rpt)
        out_rpts.append(new_rpt)
    
    return out_rpts

def organize_repeats(gm, rpts, chr, name = "", type = "", skip = 1, max_num = None):
    
    rpts = get_repeats(gm, rpts, chr, name=name, type=type, skip=skip)
    rpts = list(sorted(rpts, key = lambda rpt:-len(rpt.repeat)))
    
    print(f"identified {len(rpts)} repeats")
    if max_num:
        rpts = rpts[:max_num]
    
    return rpts

def show_organized_repeats(org_rpts):
    
    data = []
    for irpt in range(len(org_rpts)):
        row = []
        for jrpt in range(len(org_rpts)):
            if irpt == jrpt:
                row.append(None)
                continue
            
            algn = align.align_sequences(org_rpts[irpt].repeat, org_rpts[jrpt].repeat)[0]
            row.append(algn.score)
        data.append(row)
        row_nn = [r for r in row if not r is None]
        rmean = np.mean(row_nn)
        
        print(f"distance stats for repeat {irpt}: mean={rmean:0.1f}, sd = {np.std(row_nn):0.3f}, min={min(row_nn):0.1f}, max={max(row_nn):0.1f}")
    
    draw.heatmap(data, center = None, row_space = 0, col_space = 0)
    
    all_seqs = [rpt.repeat for rpt in org_rpts]
    sf = lambda sa, sb: align.score_sequences(sa, sb)
    seqs_org, inds_org = compare.organize_sequences(all_seqs, score_func = sf)
    
    for s in seqs_org:
        print(s)
    
    sa = all_seqs[0]
    for sb in all_seqs[1:]:
        algn = align.align_sequences(sa, sb)[0]
        algn.print(show_consensus = True)
    
    data_reord = []
    for ni in inds_org:
        nd = data[ni]
        ndr = [nd[nii] for nii in inds_org]
        data_reord.append(ndr)
    
    draw.heatmap(data_reord, center = None, row_space = 0, col_space = 0)

def test_phylo(seqs):
    
    ref = seqs.pop(0)
    antiref = seqs.pop(-1)
    
    phylo = compare.Phylo()
    phylo.add_reference(ref, antiref = antiref)
    
    for seq in seqs:
        phylo.insert_sequence(seq)
        
    print(phylo.tree)

def test_multiple_algn(seqs, ref = ""):
    
    if not ref:
        ref = seqs.pop(0)
    
    cons = [ref]
    algns = []
    for sb in seqs:
        malgn = align.align_sequences(ref, sb)[0]
        algns.append(malgn)
        ref = malgn.consensus
        cons.append(ref)
    
    return algns, cons


def get_long_repeats(gm, rpts, chr, min_seq_len = None, min_rpt_len = None, rpt_type = "Simple_repeat"):
    
    out_rpts = []
    
    num_rpts = 0
    for i, rpt in enumerate(rpts.stream(chr, start = 1e6)):
        
        if rpt_type and rpt.get("type") != rpt_type:
            continue
        
        num_rpts += 1
        
        loc = repeats.Loc.from_feature(rpt)
        
        rpt_motif = rpt.get("attributes",{}).get("motif","")
        
        if min_seq_len and loc.end - loc.start < min_seq_len:
            continue
        if min_rpt_len and len(rpt_motif) < min_rpt_len:
            continue
        
        new_rpt = repeats.RepeatSeq.get_instance(gm, rpt, context_sz = 128)
        out_rpts.append(new_rpt)
    
    return out_rpts

def get_test_seqs(max_len = None):
    
    sa, sb, sc = "GGCTGGGTGTGGTGGCTCAGGTCTGTAATCCCAGCACTTTGGGAAGCTGAGGCGGCAGGAGCACTTGGAGCCAAGAGTTTGAGACCAGCCTGGACAACAAAGCAAGACCCTGTCGCTACAAAACATTGAACAAAAATTAGCTGGCTGTGGTGGCAAGTGCCTGTAATCCCAGCTACTTGGGAGACGGAGGTGGGAGGGTGTCTTGAGCCCAGGAGTTGAAGGCTGCAGTGAGCTATATGAGTTGAAGGCTGCAGTGAGCTATGATCACACCACTGAACTCCAGCCAAGGCCACAGAACAAGATCCTGTCTCTAAAAAATAATAATAATAATAAATAAA","GGCTGAGCATGCTGACTCATGCTTTTAATCTTAGCACTTTGGGAGGCCGAGGCAGGAGGATCACTTGAGGCCAGGAATTCGAGACCAGCCTGGGCAACATAGTGGAACCCCATCTCTATGAAAAATAAATTAAAAATTTTAAAAATTAGCCAAACATATTGGTGTGCTCCTATAGTCCCAGATCCTTGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAAGTTGAGGATGCAGTGAGCTACGATTGCACCACTGCACTCCAGCCTGGATGACAGAGTAAGACCTTGTCTCCAAAAGAAGAAATAAATAAATAAATAAA","GGCCAGGCACAGTGACTGGTGCCTGTAATCCTAGCTACTTGGGAGGCTGAAACTGGAGGATCGCTTGAGGCCAGGAGTTGGAGACCAGCCTGGGCAACATAGTTGAGAGCCTGTCCCCTGTCTCTACAAACTGAAAATAAAAAATTAGCCGGGCATGGTGGCCTGCACCTGTAGTCCCAGCCACTCAGTAGGCCAAGGTAGGAGGATTGCTTGAGTCCAGGAGATCGAGGTTGCAGTGAGACGTGATCGCACCACAGCACTCTAGTCTGGGCAACAGAGCAAGAGCCTGTCTTTAAACACGGAAAGAA"

    if max_len:
        sa, sb, sc = sa[:max_len], sb[:max_len], sc[:max_len]
    
    return sa, sb, sc

def cons_to_alias(cons, aliases):
    
    aliases = aliases
    
    outseq = []
    
    ali_dict = {}
    
    for ali in aliases:
        refs = bio.ALIASES_FULL.get(ali)
        for r in refs:
            if r not in ali_dict:
                ali_dict[r] = ali
    
    for b in cons:
        outseq.append(ali_dict.get(b, 'N'))
        
    return "".join(outseq)

def cons_to_SW(cons):
    return cons_to_alias(cons, "SW")

def cons_to_RY(cons):
    return cons_to_alias(cons, "RY")

def cons_to_MK(cons):
    return cons_to_alias(cons, "MK")

def main():
    
    gm, rpts = load_genome()
    
    test_chr = "19"
    # rpt_name = "LTR38C"
    rpt_name = "AluJr"
    rpt_type = ""
    skip = 50
    
    # display_repeats(gm, rpts)
    
    # my_rpts = get_repeats(gm, rpts, test_chr, name = rpt_name, type = rpt_type, skip = skip)
    # my_rpts = organize_repeats(gm, rpts, test_chr, name=rpt_name, type=rpt_type, skip = 50, max_num = 8)
    
    # show_organized_repeats(my_rpts)
    
    # seqs = [rpt.repeat for rpt in my_rpts]
    # test_phylo(seqs)
    
    # malgn, cons = test_multiple_algn(seqs)
    
    # for ma in malgn:
    #     ma.print()
    
    chunksz = 4*64*256
    start = 10e6
    length = 3*256*chunksz
    
    long_rpts = get_long_repeats(gm, rpts, test_chr, min_rpt_len = 10)
    
    for rpt in long_rpts:
        mtf = rpt.raw_repeat.get("attributes",{}).get("motif","")
        
        rpt_seq = rpt.upstream[-16:] + rpt.repeat + rpt.downstream[:16]
        
        folded = repeats.fold_repeat(mtf, rpt_seq)
        cons = bio.merge_all(folded)
        cons_sw = cons_to_SW(cons)
        
        print(mtf)
        print(mtf + mtf)
        print(bio.reverse_complement(mtf+mtf))
        print(rpt.repeat)
        # print(cons)
        
        print()
        
        cons_re = find.consensus_to_re(cons)
        print(cons_re)
        
        ptrn = re.compile(cons_re)
        
        seq_fn = lambda seq, feats: lambdas._seq_pattern(seq, feats, ptrn)
        
        # gm.display_chromosomal_quantity(test_chr, seq_fn, chunksz = chunksz, start=start, length=length)
        
        input()

    
    

        
    
    
    
"""
2:1472140
AGTGCTCTTCCTGTGATCCAAATGCTCATGGCATGCCCTTTCCATGACGTGAATGCTCATGGCGTGCCCCTCCTGTGTTCTGAATGCTTATAGAGTGCCCTTCACATGATCTGAGTGCTCATGGTGTGTCCTTTCTGTGATCCAAATACTCATGCTATGTCCTTCCCTTGATCCAAATGCTCATGGTGTGTCCTTCTGATGATCTGAGTGCTCATAGTGTGCTCTTCCCACGATATGAGT

"""

    
if __name__=="__main__":
    main()




