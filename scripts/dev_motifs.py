
import random
import string
import numpy as np
import matplotlib.pyplot as plt
import re

from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, splice_donor, splice_branch, splice_acceptor
from ggene import seqs
from ggene import draw
from ggene.draw import scalar_to_text_8b, scalar_to_text_16b, scalar_to_text_nb
from ggene.seqs import reverse, reverse_complement, COMPLEMENT_MAP
from ggene.translate import Ribosome
from ggene.motifs import dyad
from ggene.motifs import utils
from ggene.draw import Colors
from ggene.motifs.dyad import Dyad, find_all_dyads, frequency_rank_dyads
from ggene.genomemanager import GenomeManager


SCALE = " ▁▂▃▄▅▆▇█"

def load_genome():
    return GenomeManager()


def get_natural_seq(gm, seq_len = 240, seq_center = 1124574):
    c = seq_center
    seq = gm.get_sequence(1, c-seq_len//2, c+seq_len//2)
    return seq

def get_hammerhead():
    
    seq = "YYRRGCCGUUACCURCAGCUGAUGAGCUCCAARAAGAGCGAAACCXRXYAGGUCCUGYAGUAYUGGCYXRXXXXXX"
    subs = {"Y":"C","R":"G","C":"C","G":"G","U":"U","A":"A"}
    out = []
    xbases = "ACGUACGUA"
    
    for s in seq:
        
        ssub = subs.get(s, s)
        if ssub == "X":
            # ssub = random.sample(bases, 1)[0]
            ssub = xbases[0]
            xbases = xbases[1:]
        
        out.append(ssub)
    return "".join(out)

def get_test_sequence():
    seq = list(range(-100, 0, 14)) + list(range(0, -50, 10)) + list(range(-50,25,9)) + list(range(25,100,18)) + list(range(100, 50, -19)) + [50]
    return seq

def get_random_sequence(seq_len, minval = 0, maxval = 100, var = 1, bias = 0, nint = 1):
    
    rand_seq_d = [(2*random.random()-1)*var+bias for i in range(seq_len)]
    
    for n in range(nint):
        rand_seq = [sum(rand_seq_d[:i]) for i in range(seq_len)]
        rand_seq_d = rand_seq
    
    _max= max(rand_seq)
    _min = min(rand_seq)
    _ran = _max - _min
    
    tgt_ran = maxval - minval
    
    rand_seq_sc = [tgt_ran*(r-_min)/_ran+minval for r in rand_seq]
    return rand_seq_sc

def get_color_scheme(name):
    """
    returns bg, fg
    """
    if name == "gray":
        return 244, 236
    elif name == "blue":
        return 17, 38
    elif name == "foggy":
        return 36, 67
    elif name == "dusty":
        return 188, 138
    elif name == "ruddy":
        return 179, 131
    elif name == "icy":
        return 146, 225
    elif name == "vscode":
        return 234, 131
    elif name == "test":
        return 234, 65
    else:
        return 0,1



def test_convolve2(seq1, seq2, do_rc = False):
    
    bc, fc = get_color_scheme("test")
    
    runs, inds = seqs.convolve_longest_subseq(seq1, seq2)
    minval = 0
    meanrun = np.mean(runs)
    maxrun = max(runs)
    print(f"{len(runs)} runs identified, mean {meanrun:0.2f}, max {maxrun}")
    
    # print("autocorrelation response")
    # res = scalar_to_text_16b(runs, fg_color = fc, bg_color = bc)
    # for r in res:
    #     print(r)
    # print(seq1)
    # print()
    
    print_longest_subseqs(seq1, seq2, runs, inds, extra = 6, topk = 3, do_rc = do_rc)
    
    return seq1, runs, inds

def print_longest_subseqs(seq, seq2, runs, inds, min_run = 4, extra = 3, topk = 10, do_rc = False):
    
    seq_len = len(seq)
    
    high = 142
    low = 8
    
    srt = sorted(zip(runs, inds), key = lambda a:-a[0])
    top_runs = srt[:topk]
    
    seq1s = []
    seq2s=[]
    
    for r, (st, sh) in top_runs:
        if st < 0 or r < min_run:
            continue
        
        st2, sh2 = st, sh
        
        mn1 = max(0, st-extra)
        mx1 = min(seq_len, st+r+extra)
        mn2 = max(0, st2-sh2-extra)
        mx2 = min(seq_len, st2+r-sh2+extra)
        
        print(f"start = {st}, run length = {r}, shift = {sh}")
        res = draw.highlight_subsequences([seq[mn1:st], seq[st:st+r], seq[st+r:mx1]],[low, high, low], delimit= " ")
        print(res)
        res=draw.highlight_subsequences([seq2[mn2:st2-sh2], seq2[st2-sh2:st2+r-sh2], seq2[st2+r-sh2:mx2]],[low, high, low], delimit= " ")
        print(res)
        print()
        
        seq1s.append(seq[st:st+r])
    seq1s = list(set(seq1s))
    for eachseq in seq1s:
        subseqs = [eachseq]
        if do_rc:
            subseqs.append(reverse_complement(eachseq))
        draw.highlight_sequences(seq, subseqs)

def reverse_complement_index(seq_len, run, ind):
    return (seq_len - run - ind[0], -ind[1])
    

def test_flip(bd = 16):
    
    bc, fc = get_color_scheme("test")
    _bc = bc
    
    # effects = [1, 2, 3, 4, 5, 6, 9, 51, 52, 53, 90, 100] * 10
    effects = range(0, 10)
    
    for eff in effects:
        # effect = f"\x1b[{eff}m"
        effect = ""
        
        fc = random.randint(20, 230)
        
        seq = get_random_sequence(256, var = 16, bias = 0.1)
        res = draw.scalar_to_text_nb(seq, fg_color = fc, bg_color = bc, flip= False, bit_depth = bd, effect = effect)
        for r in res:
            print(r)
        
        # bc = random.randint(20, 230)
        if eff == effects[-1]:
            bc = _bc
            
        rseq=list(reversed(seq))
        res = draw.scalar_to_text_nb(rseq, fg_color = fc, bg_color = bc, flip= True, bit_depth = bd, effect = effect)
        for r in res:
            print(r)



def scan_for_repeats(gm, chr, start, step):
    
    old_vocab = "ATGC"
    
    while True:
        
        
        seq = gm.get_sequence(chr, start, start+step)
        seq = seqs.convert_to_vocab(seq, seqs.vocab, from_vocab = old_vocab)
        rcseq = reverse_complement(seq)
        
        print(f"***************{chr}:{start}-{start+step}***************")
        test_convolve2(seq, seq)
        test_convolve2(seq, rcseq, do_rc = True)
        
        res = input("keep going?")
        if 'n' in res.lower():
            break
        print()
        
        start+=step

    
        
def main():
    pass    
    
    seqs.set_vocab("ATGC")
    # seqs.set_vocab("■□▼▽")
    print(seqs.vocab, seqs.COMPLEMENT_MAP)

    bc, fc = get_color_scheme("test")

    # gm = load_genome()
    
    chr = 2
    seq_center = 1124574
    seq_len = 1024
    
    input()
    
    # scan_for_repeats(gm, chr, seq_center, 1024)
    
    # seq = get_natural_seq(gm, seq_len = seq_len)
    # seq = get_hammerhead()
    # seq = "AUGCGCGCGCAUGGCGCUA"
    # rcseq = reverse_complement(seq)


    
    
if __name__=="__main__":
    main()
