
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
from ggene.draw import Color
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
    

def find_nice_colors(seq_len = 256):
    
    # seq = list(range(-100, 0, 14)) + list(range(0, -50, 10)) + list(range(-50,25,9)) + list(range(25,100,18)) + list(range(100, 50, -19)) + [50]

    n = 0
    while True:
        seq = get_random_sequence(seq_len, var = 4)
        offset = 24
        # fc = 170
        bc = 234
        # bc = random.randint(20, 230-offset)
        fc = random.randint(20,230)
        # fc = min(bc+offset, 230)
        
        if n == 0:
            fc = 236
            bc = 244
        
        res = scalar_to_text_nb(seq, bg_color = bc, fg_color = fc, bit_depth = 16)
            
        print(f"bg color = {bc}, fg color = {fc}")
        for r in res:
            print(r)
        print()
        n += 1
        
        if n%5 == 0:
            res=input()
            if 'n' in res.lower():
                break
    
    return bc

def looks_cool(seq_len = 128, num_bias = 20, bias_range = 1.0):
    
    import os
    with open("./test.txt", "w") as f:
        f.write("")
    
    for i in range(num_bias):
        b = bias_range*(i/num_bias - 0.5)
        rseq = get_random_sequence(seq_len, bias = b)
        print()
        
        res = scalar_to_text_16b(rseq)
        for r in res:
            print(r)
        
        rclean = draw.clean_scalar_text(res)
        rclfull = "\n".join(rclean)
        
        with open("./test.txt", "a") as f:
            f.write(rclfull)
            # for r in rclean:
            #     f.write(r + "\n")
            f.write("\n\n")

def test_nb(seq_len = 256):
    
    bc, fc = get_color_scheme("test")
    
    rseq = get_random_sequence(seq_len, var = 4)
    
    test_bds = [8, 16, 24, 32, 40, 48, 56]
    for bd in test_bds:
        res = scalar_to_text_nb(rseq, bit_depth = bd, fg_color = fc, bg_color = bc)
        print(f"{bd}-BIT:")
        for r in res:
            print(r, "-")


def convolve(seq1, seq2, comparison_func = None):
    
    if not comparison_func:
        cmp = lambda x, y:x==y
    else:
        cmp = comparison_func
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]
    
    start = seq_len // 2
    ips = []
    comp_ips = []
    n = 0
    for t in range(-start, start + 1):
        
        sslen = seq_len - abs(t)
        s1start = max(t, 0)
        seq1t = seq1[s1start:s1start + sslen]
        s2start = max(-t, 0)
        seq2t = seq2[s2start:s2start + sslen]
        
        summ = 0
        csumm = 0
        for sa, sb in zip(seq1t, seq2t):
            
            csb = seqs.COMPLEMENT_MAP.get(sb)
            
            if cmp(sa, sb):
                summ += 1/sslen
            elif cmp(sa, csb):
                csumm += 1/sslen
        
        if t == 0:
            ips.append(0)
        else:
            ips.append(summ)
        comp_ips.append(csumm)
        n+=1
    
    return ips, comp_ips

def convolve_longest_subseq(seq1, seq2):
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]
    
    start = seq_len // 2
    runs = []
    inds = []
    found = set()
    
    for i,t in enumerate(range(-start, start+1)):
        
        sslen = seq_len - abs(t)
        s1start = max(t, 0)
        seq1t = seq1[s1start:s1start + sslen]
        s2start = max(-t, 0)
        seq2t = seq2[s2start:s2start + sslen]
        
        max_run = 0
        max_run_end = -1
        max_run_shift = -1
        run = 0
        for j, (sa, sb) in enumerate(zip(seq1t, seq2t)):
            
            if sa==sb:
                run += 1
            else:
                run = 0
                
            if run > max_run:
                max_run = run
                max_run_end = j + s1start
                max_run_shift = t
        
        if t == 0:
            runs.append(0)
        else:
            runs.append(max_run)
        s1sp = max_run_end - max_run + 1
        s2sp = s1sp - max_run_shift
        if (s1sp, s2sp) in found:
            inds.append((-1, -1))
            continue
        else:
            found.add((s1sp, s2sp))
            found.add((s2sp, s1sp))
            
            newind = (max_run_end - max_run + 1, max_run_shift)
            inds.append(newind)
    
    return runs, inds

def convolve_generator(seq1, seq2):
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]
    
    start = seq_len // 2
    
    for t in range(-start, start):
        
        sslen = seq_len - abs(t)
        seq1t = seq1[max(t, 0):max(t, 0) + sslen]
        seq2t = seq2[max(-t, 0):max(-t, 0) + sslen]
        yield seq1t, seq2t



def test_convolve(seq, cmp = None):
    
    bc, fc = get_color_scheme("test")
    
    summ, csumm = convolve(seq, seq, comparison_func = cmp)
    minval = 0
    
    print("autocorrelation response")
    res = scalar_to_text_16b(summ, minval = minval, fg_color = fc, bg_color = bc)
    for r in res:
        print(r)
    # print(seq)
    
    # print()
    # print("reverse complement autocorrelation response")
    res = scalar_to_text_16b(csumm,minval = minval, fg_color = fc - 6, bg_color = bc)
    for r in res:
        print(r)                                                         
    print(seq)
    return seq, summ, csumm

def test_convolve2(seq1, seq2, do_rc = False):
    
    bc, fc = get_color_scheme("test")
    
    runs, inds = convolve_longest_subseq(seq1, seq2)
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

def test_ansi():
    fc = "\x1b[38;5;120m"
    bc = "\x1b[48;5;90m"
    test_ints = [1, 2, 3, 4, 6, 7, 9, 21, 53]
    # test_codes = ["G"]
    # for ti in test_codes:
    for ti in test_ints:
        t = f"\x1b[{ti}m"
        print(ti,"hey deeev", t, "hi hellow how are you", draw.RESET, "thanks i m doing well, you?")
        print(ti,"hey deeev", fc, bc, t, "hi hellow how are you", draw.RESET, "thanks i m doing well, you?")
        print()

def test_hor():
    
    bc, fc = get_color_scheme("test")
    for i in range(10):
        seq = get_random_sequence(16)
        res = draw.scalar_to_text_nbh(seq, fg_color = fc, bg_color = bc,bit_depth = 80, flip = True, effect = "")
        for r in res:
            print(r)

def test_overflow():
    
    bc, fc = get_color_scheme("test")
    _bc = bc
    bd = 32
    
    effects = range(0, 5)
    
    for eff in effects:
        
        fc = random.randint(20, 230)
        seq = get_random_sequence(256, var = 256, nint=1)
        res = draw.scalar_to_text_nb(seq, maxval = 75, fg_color = fc, bg_color = bc, bit_depth = bd)
        for r in res:
            print(r)
        # bc = random.randint(20, 230)
        # fc = bc
        bc = random.randint(20, 230)
        # if eff == effects[-1]:
        #     bc = _bc
        rseq=list(reversed(seq))
        # rseq = seq
        res = draw.scalar_to_text_nb(rseq, maxval = 75, fg_color = fc, bg_color = bc, flip= True, bit_depth = bd)
        for r in res:
            print(r)
            
        # fc = bc
        # bc = random.randint(20, 230)

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
        

def test_border(ntests = 5, bd = 16):
    
    bc, fc = get_color_scheme("test")
    _bc = bc
    
    for eff in range(ntests):
        
        fc = random.randint(20, 230)
        
        seq = get_random_sequence(256, var = 16, bias = 0.1)
        res = draw.scalar_to_text_nb(seq, fg_color = fc, bg_color = bc, flip= False, bit_depth = bd, effect = "")
        for r in res:
            print(r)
        print()
        
        rseq = seq
        res = draw.scalar_to_text_nb(rseq, fg_color = fc, bg_color = bc, flip= False, bit_depth = bd, effect="border")
        for r in res:
            print(r)
        print()

def test_mid(ntests = 5):
    bc, fc = get_color_scheme("test")
    _bc = bc
    
    for eff in range(ntests):
        
        fc = random.randint(20, 230)
        
        seq = get_random_sequence(256, var = 16, bias = 0.1)
        res = draw.scalar_to_text_nb(seq, fg_color = fc, bg_color = bc, flip= False, bit_depth=16, effect = "")
        for r in res:
            print(r)
        print()
        
        rseq = seq
        res = draw.scalar_to_text_mid(rseq, fg_color = fc, bg_color = bc)
        for r in res:
            print(r)
        print()

    
        
def reveal_escapes(sctext):
    for r in sctext:
        print(repr(r))
    print()
        
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
    
    # scan_for_repeats(gm, chr, seq_center, 1024)
    
    # seq = get_natural_seq(gm, seq_len = seq_len)
    # seq = get_hammerhead()
    # seq = "AUGCGCGCGCAUGGCGCUA"
    # rcseq = reverse_complement(seq)
    
    
    # out = test_convolve(seq)
    # test_convolve2(seq, seq)
    # test_convolve2(seq, rcseq, do_rc = True)
    
    # test_flip(bd=32)
    # test_ansi()
    # test_hor()
    # test_overflow()
    # test_border(5)
    
    # test_mid(5)
    
    bc, fc = get_color_scheme("test")
    
    seqa = get_random_sequence(16)
    seqb = get_random_sequence(16)
    seqa, seqboff, off, ebd = draw.plot_adjacent(seqa, seqb)
    
    print(f"offset: {off}, effective bit depth: {ebd} vs {2*16}")
    print(" ".join(str(a) for a in seqa))
    print(" ".join(str(b) for b in seqboff))
    print(" ".join(str(a+b) for a,b in zip(seqa,seqboff)))
    
    resa = draw.scalar_to_text_mid(seqa, fg_color = fc, bg_color = bc)
    for r in resa:
        print(r)
    resb = draw.scalar_to_text_mid(seqb, fg_color = fc, bg_color = bc)
    for r in resb:
        print(r)
    
    # qseq = draw.quantize(seq, 16)
    # mseq = draw.quantize(seq, 16, mid=True)
    # print(" ".join(format(s, "0.1f") for s in seq))
    # print(" ".join(str(q) for q in qseq))
    # print(" ".join(str(m) for m in mseq))
    # pyrpur = {"AG":"R", "UC":"Y"}
    # yrcmp = lambda x, y: pyrpur.get(x) == pyrpur.get(y)
    
    # summ, csumm = test_convolve(gm, cmp = yrcmp, seq_len = seq_len)
    
    # print("alternate comparison")
    # res = scalar_to_text_16b(summ, minval = minval, fg_color = fc, bg_color = bc)
    # for r in res:
    #     print(r)
    
    # find_nice_colors()
    
    # looks_cool()
    
    # with open("./test.txt") as f:
    #     for i,line in enumerate(f):
    #         print(i,line)

    # test_nb()
    
    # find_nice_colors()    
    
    
if __name__=="__main__":
    main()
