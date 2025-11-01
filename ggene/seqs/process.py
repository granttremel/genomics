

from typing import Optional, Dict, List, Tuple, Any
from .bio import reverse_complement, get_aliases
from .vocab import VOCAB


def normalize_sequence(seq):
    """
    would be nice. define a seq as "normal", so that norm_seq(seq) = seq, norm_seq(rc(seq)) = seq
    A..A -> U..U
    A..U -> A..U # not possible, also U..A
    A..G -> C..U
    A..C -> G..U
    G..U -> A..C
    each pair AB -> {0,1}, 0 meaning do nothing, 1 meaning do rc
    """
    
    b0 = seq[0]
    bl = seq[-1]
    
    inds = {"A":0,"U":2,"G":1,"C":3}
    py = "UC"
    pu = "AG"
    
    v0 = inds[b0]
    v1 = inds[bl]
    
    # if reverse_complement(b0) == bl:
    #     b0 = seq[1]
    #     bl = seq[-1]
    
    if b0 in py and bl in pu or b0 in py and bl in pu:
        pass
    
    if seq[0] in 'UC':
        return reverse_complement(seq)
    else:
        return seq
    
    # if (v0 - v1) % 2 == 1:
    #     return reverse_complement(seq)
    # else:
    #     return seq


def convolve(seq1, seq2, comparison_func = None, fill = 0, scale = None):
    
    if not comparison_func:
        cmp = lambda x, y:x==y
    else:
        cmp = comparison_func
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]
    rcseq2 = reverse_complement(seq2)
    
    max_shift = seq_len//2
    
    ips = []
    comp_ips = []
    n = 0
    for t in range(-max_shift, max_shift + 1):
        
        sslen = seq_len - abs(t)
        if scale is None:
            sc = sslen
        else:
            sc = scale
        
        s1start = max(t, 0)
        s1c = s1start + sslen//2
        seq1t = seq1[s1c - sc//2:s1c+sc//2]
        s2start = max(-t, 0)
        s2c = s2start + sslen//2
        seq2t = seq2[s2c - sc//2:s2c+sc//2]
        rcseq2t = rcseq2[s2c - sc//2:s2c+sc//2]
        
        summ = 0
        csumm = 0
        for sa, sb, rcsb in zip(seq1t, seq2t, rcseq2t):
            
            if not sa in VOCAB or not sb in VOCAB:
                continue
            
            if cmp(sa, sb):
                summ += 1/sc
            elif cmp(sa, rcsb):
                csumm += 1/sc
        
        if t == 0:
            ips.append(fill)
            comp_ips.append(fill)
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
    