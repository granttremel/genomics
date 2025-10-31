
import random
from statistics import correlation
import string
import numpy as np
import matplotlib.pyplot as plt
import re

from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, splice_donor, splice_branch, splice_acceptor
from ggene import seqs
from ggene import draw
from ggene.draw import scalar_to_text_8b, scalar_to_text_16b, scalar_to_text_mid, scalar_to_text_nb
from ggene.seqs import convolve, convolve_longest_subseq, reverse, reverse_complement, COMPLEMENT_MAP
from ggene.translate import Ribosome
from ggene.motifs import dyad
from ggene.motifs import utils
from ggene.draw import Color, get_color_scheme
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

def get_random_sequence(seq_len, minval = -100, maxval = 100, var = 1, bias = 0, nint = 1):
    
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

def autocorrelate_sequence(seq):
    
    bc, fc = get_color_scheme("test")
    
    ip, rcip = convolve(seq, seq)
    
    print("Autocorrelation:")
    res = draw.scalar_to_text_mid(ip, fg_color = fc, bg_color = bc)
    for r in res:
        print(r)
    print("reverse complement:")
    res = draw.scalar_to_text_mid(rcip, fg_color = fc, bg_color = bc)
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

def test_arrows(arrows = []):
    
    if not arrows:
        all_arrows = [str(chr(i)) for i in range(0x2190, 0x21FF+1)]
        all_arrows += [str(chr(i)) for i in range(0x27F0, 0x27FF+1)]
        all_arrows += [str(chr(i)) for i in range(0x2900, 0x297F+1)]
        all_arrows += [str(chr(i)) for i in range(0x2794, 0x27BE+1)]
    else:
        all_arrows = arrows
    
    for i,(hex_code, arr) in enumerate(all_arrows):
        delimit = 10*" "
        arrstr = []
        arrstr.extend(str(i))
        arrstr.append(   hex(ord(arr)).upper()            )      
        arrstr.append(f"({hex(ord(arr)).upper()}, '{arr}'), #")
        arrstr.extend((arr + "-"*5,"-"*5 + arr))
        print(delimit.join(arrstr))
        print()
    pass

def view_seq(gm, chr, c, len):
    
    bc, fc = get_color_scheme("test")
    
    
    for i in range(-3, 4):
        
        shift = 5*i
        seq = gm.get_sequence(chr, c-len//2 + shift, c+len//2+shift)
        ip, rcip = seqs.convolve(seq, seq, fill = 0.25)
        
        rule_palindrome(seq, 6)
        rule_palindrome(seqs.complement(seq), 6)
        
        res = draw.scalar_to_text_16b(ip, fg_color = fc, bg_color = bc)
        for r in res:
            print(r)
        
        res = draw.scalar_to_text_16b(rcip, fg_color = fc, bg_color = bc, flip=True)
        for r in res:
            print(r)
        
        print(" ".join([format(s, "0.1f") if s>0.05 else "!0.0!" for s in ip]))
        
def rule_palindrome(seq, dlen):
    
    c = len(seq)//2
    nrules = (len(seq)//dlen)//2 + 2
    
    
    # cols = [random.randint(20, 230),random.randint(20, 230)]
    cols = [142, 84]
    colors = {}
    subseqs =[]
    starts = {}
    for n in range(2,nrules):
        minidx = c - int(dlen*(n - 0.5))
        maxidx = c + int(dlen*(n - 0.5))
        # print(minidx,c, maxidx)
        s1 = seq[minidx: minidx+dlen]
        s2 = seq[maxidx-dlen: maxidx]
        colors[s1] = cols[n%2]
        colors[s2] = cols[n%2]
        subseqs.append(s1)
        subseqs.append(s2)
        starts[s1] = [minidx]
        starts[s2] = [maxidx-dlen]
    
    draw.highlight_sequences(seq, subseqs, start_pos= starts, colors=colors, show_key=False)
    
        
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

    gm = load_genome()
    
    # chr = 2
    # seq_center = 1124574
    # seq_len = 1024
    
    # weird place
    chr = 2
    seq_center = 1012063
    # also seq_center = 1013230
    seq_len = 256
    
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
    
    view_seq(gm, chr, seq_center, seq_len)
    
    # for i in range(10):
    
    #     seqa = get_random_sequence(64)
    #     seqb = get_random_sequence(64)
    #     resa, qnt = draw.scalar_to_text_mid(seqa, fg_color = fc, bg_color = bc)
    #     for r in resa:
    #         print(r)
    #     print()
    #     print(" ".join([format(q, "g") for q in qnt]))
        # resb = draw.scalar_to_text_mid(seqb, fg_color = fc, bg_color = bc)
        # for r in resb:
        #     print(r)
        # print()

        
        # seqa, seqboff, off, ebd = draw.plot_adjacent(seqa, seqb)

        # print(f"offset: {off}, effective bit depth: {ebd} vs {2*16}")
        # print(" ".join(str(a) for a in seqa))
        # print(" ".join(str(b) for b in seqboff))
        # print(" ".join(str(a+b) for a,b in zip(seqa,seqboff)))

        # resa = draw.scalar_to_text_mid(seqa, fg_color = fc, bg_color = bc)
        # for r in resa:
        #     print(r)
        # print()
        # resb = draw.scalar_to_text_mid(seqb, fg_color = fc, bg_color = bc)
        # for r in resb:
        #     print(r)
        # print()

    # print(draw._arrows)
    # print(draw._arrows2)
    # print(draw._arrows3)
    
    # print("⇇=")
    
    # test_arrows()
    
    good_arrows = [
        (0X27F5, '⟵'), #          ⟵-----          -----⟵
        (0X27F6, '⟶'), #          ⟶-----          -----⟶
        (0X21FD, '⇽'), #          ⇽-----          -----⇽
        (0X21FE, '⇾'), #          ⇾-----          -----⇾
        (0X21E6, '⇦'), #          ⇦-----          -----⇦
        (0X21E8, '⇨'), #          ⇨-----          -----⇨
        (0X21E0, '⇠'), #          ⇠-----          -----⇠
        (0X21E2, '⇢'), #          ⇢-----          -----⇢
        (0X21D0, '⇐'), #          ⇐-----          -----⇐
        (0X21D2, '⇒'), #          ⇒-----          -----⇒
        (0X21C0, '⇀'), #          ⇀-----          -----⇀
        (0X21C1, '⇁'), #          ⇁-----          -----⇁
        (0X21C0, '⇀'), #          ⇀-----          -----⇀
        (0X21BD, '↽'), #          ↽-----          -----↽
        (0X21BC, '↼'), #          ↼-----          -----↼
        (0X2190, '←'), #          ←-----          -----←
        (0X2192, '→'), #          →-----          -----→
        
    ]
    
    might_be_useful = [
            (0X293A, '⤺'), #          ⤺-----          -----⤺
            (0X293A, '⤺'), #          ⤺-----          -----⤺
            (0X293B, '⤻'), #          ⤻-----          -----⤻
            (0X293C, '⤼'), #          ⤼-----          -----⤼
            (0X293D, '⤽'), #          ⤽-----          -----⤽

            (0X21B6, '↶'), #          ↶-----          -----↶
            (0X21B7, '↷'), #          ↷-----          -----↷
            (0X219C, '↜'), #          ↜-----          -----↜
            (0X219D, '↝'), #          ↝-----          -----↝
    ]
    
    # print("good")
    # test_arrows(good_arrows)
    # print()
    # print("idk")
    # test_arrows(might_be_useful)
    
    # kinda going with:
    # (0X21C0, '⇀'), #          ⇀-----          -----⇀
    # (0X21BD, '↽'), #          ↽-----          -----↽
    
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
