
import numpy as np
import math
import itertools

from typing import Optional, Dict, List, Tuple, Any
from .bio import reverse_complement, get_aliases, get_minimal_alias
from .vocab import VOCAB
from .process import correlate_longest_subseq
import random

def consensus_to_re(seq, err = 0):
    
    ia, iar = get_aliases()
    outstr = []
    for s in seq:
        if s in VOCAB:
            outstr.append(s)
        elif s in ia:
            outstr.append(f'[{ia.get(s)}]')
        else:
            outstr.append(".")
    
    return "".join(outstr)


def reverse_complement_re(pattern):
    plist = list(pattern)
    rcp = []
    
    in_grp = False
    grp = []
    for c in plist:
        if c == '[':
            in_grp = True
        elif c == ']':
            in_grp = False
            new_grp_bs = [v for v in VOCAB if not v in grp]
            rcp.append(f"[{"".join(new_grp_bs)}]")
            grp = []
        else:
            if in_grp:
                grp.append(c)
            else:
                rcp.append(c)
    return "".join(rcp)


def compare_sequences(s1: str, s2: str, err_tol: Optional[int] = None, allow_alias = False) -> int:
    """Helper to compare two sequences and count mismatches."""
    if len(s1) < len(s2):
        return len(s1)
    
    if err_tol is None:
        err_tol = len(s1)

    ia, iar = get_aliases(VOCAB)
    
    nerr = 0
    for a, b in zip(s1, s2):
        if a == b:
            pass
        elif allow_alias and a in ia and b in ia[a]:
            pass
        elif allow_alias and b in ia and a in ia[b]:
            pass
        else:
            nerr += 1
        if nerr > err_tol:
            return nerr
    return nerr


def find_subsequence(seq, subseq, frame_start = -1):
    starts = []
    # ends = []
    frame_size = len(subseq)
    
    # Find all occurrences of this subsequence
    start = 0
    while True:
        pos = seq.find(subseq, start)
        if pos == -1:
            break
        if not frame_start < 0:
            fr = (pos - frame_start) % frame_size
            if fr > 0:
                start = pos+1
                continue
        
        starts.append(pos)

        start = pos + 1  # Continue searching from next position
    return starts

def find_subsequences(seq, subseqs, do_rc = False, err_tol = 0):
    max_ss = max([len(ss) for ss in subseqs])
    subseq_pos = {s:[] for s in subseqs}
    for i in range(len(seq)):
        test_str = seq[i:i+max_ss]
        rctest_str = reverse_complement(test_str)
        
        for ss in subseqs:
            err = compare_sequences(test_str, ss, err_tol = err_tol)
            if err <= err_tol:
                subseq_pos[ss].append(i)
            elif do_rc:
                err = compare_sequences(rctest_str, ss, err_tol = err_tol)
                if err <= err_tol:
                    subseq_pos[ss].append(i - len(ss))
    
    return subseq_pos

def find_subsequences_fuzzy(seq, subseqs, err_tol = 1):
    
    max_ss = max([len(ss) for ss in subseqs])
    min_ss = min([len(ss) for ss in subseqs])
    subseq_pos = {s:[] for s in subseqs}
    templates = {s:set([s]) for s in subseqs}
    for i in range(len(seq)-min_ss):
        test_str = seq[i:i+max_ss]
        
        for ss in subseqs:
            err = compare_sequences(test_str, ss, err_tol = err_tol)
            if err == 0:
                subseq_pos[ss].append(i)
            elif err <= err_tol:
                subseq_pos[ss].append(i)
                templates[ss].add(test_str)
    
    return templates, subseq_pos


def frequency_rank(seqs, proc = [], min_len = 3, topk = None, do_rc = False, err_tol = 0):
    max_n = 0
    max_seq = ""
    
    # seqs = sorted(seqs, key = lambda a:-len(a))
    # if proc:
    #     proc = sorted(proc, )
    
    freqs = {}
    procout = {}
    for i in range(len(seqs)):
        seq = seqs[i]
        
        if len(seq) < min_len:
            continue
        
        key = seq
        if seq in freqs:
            key = seq
        elif do_rc and reverse_complement(seq) in freqs:
            key = seq
        else:
            for f in freqs:
                if seq in f:
                    key = f
                    break
        
        if not key in freqs:
            freqs[key] = 0
            procout[key] = []
        n = freqs.get(key, 0)
        freqs[key] = freqs[key] + 1
        
        if proc:
            procout[key].append(proc[i])
            print(f'{i}th sequence:')
            print(seq)
            proc[i].print()
        
        if n > max_n:
            max_n = n
            max_seq = key
    
    if topk is not None:
        top_seqs = sorted(freqs.keys(), key = lambda k:-freqs[k])[:topk]
        freqs = {k:freqs[k] for k in top_seqs}
        procout= {k:procout[k] for k in top_seqs}
    
    return freqs, procout, max_seq

def merge_seqs(seq1, seq2):
    
    ia, iar = get_aliases(VOCAB)
    
    outseq = []
    for s1, s2 in zip(seq1, seq2):
        if s1==s2:
            outseq.append(s1)
        else:
            smerge = get_minimal_alias(s1, s2)
            outseq.append(smerge)
            
    return "".join(outseq)

def group_templates(all_observed_seqs):
    
    err_tol = len(all_observed_seqs[0])//10
    
    oseq1 = all_observed_seqs[0]
    groups = {oseq1:[oseq1]}
    for seq in all_observed_seqs[1:]:
        grouped = False
        for group in groups:
            err = compare_sequences(seq, group, err_tol = err_tol, allow_alias = True)
            if err <= err_tol:
                if err > 0:
                    merge = merge_seqs(group, seq)
                    groups[merge] = groups[group]
                    del groups[group]
                    group = merge
                groups[group].append(seq)
                grouped = True
                break
        if not grouped:
            groups[seq] = []
    
    return groups

def compile_templates(observed_seqs):
    """
    template: ACTCACT->list of dict:base->number occurrences
    """
    if not observed_seqs:
        return []
    nobs = len(observed_seqs)
    temps = [{} for s in observed_seqs[0]]
    for oseq in observed_seqs:
        for i in range(len(oseq)):
            if len(temps) <= i:
                temps.append({})
            if not oseq[i] in temps[i]:
                temps[i][oseq[i]] = 0
            temps[i][oseq[i]] += 1/nobs
    return temps

def get_template_entropy(template):
    return [sum([-p*np.log(p) for p in pos.values() if p > 0]) for pos in template]

def sample_template(template):
    
    outseq = []
    for bdict in template:
        rs = random.random()
        vv = 0
        for k, v in bdict.items():
            vv += v
            if vv > rs:
                outseq.append(k)
                break
        
    return "".join(outseq)

def reduce_ratios(ns, factor):
    ns = list(set(ns))
    gcd = math.gcd(*ns)
    if gcd == 0:
        return 1
    else:
        return factor//gcd
    
def make_template_str(temp_list, factor = 16, reduce = True):
    if reduce:
        ns = itertools.chain.from_iterable([[int(v*factor) for k,v in ti.items()] for ti in temp_list])
        factor = reduce_ratios(ns, factor)
    outstrs = []
    for i in range(len(temp_list)):
        base_str = []
        for k, n in temp_list[i].items():
            nn = int(n*factor)
            if n==0:
                continue
            elif n==1:
                nstr = ""
            elif nn == 1:
                nstr = "1."
            else:
                nstr = str(int(n*factor))
                if nstr == "0":
                    nstr = "0."
            base_str.append(f"{nstr}{k}")
        if len(base_str) > 1:
            outstrs.append(f"[{" ".join(base_str)}]")
        else:
            outstrs.append(base_str[0])
    return "".join(outstrs)

def combine_templates(reference, templates):
    ref = [r for r in reference]
    
    for i in range(len(ref)):
        alltmp = {b:0 for b in VOCAB}
        for t in templates:
            alltmp[t[i]] += 1
        ref[i] = alltmp
    
    for d in ref:
        for b in VOCAB:
            if d[b] == 0:
                del d[b]

    return ref

def merge_template_dicts(temp, temp_new):
    
    for f, ts in temp_new.items():
        if f in temp:
            for tts in ts:
                temp[f].add(tts)
        else:
            temp[f] = ts
    return temp

def locate_repeats(seqa, seqb, zscore = 2, scale = None, do_rc = False):
    
    if do_rc:
        seqb = reverse_complement(seqb)

    runs, inds, shifts = correlate_longest_subseq(seqa, seqb, scale = scale)
    seq_len = len(seqa)
    
    cmean = np.mean(runs)
    csd = np.std(runs)
    
    maxthr = cmean + zscore*csd
    
    rpts = []
    rptinds = []
    
    for i,(v, j, sh) in enumerate(zip(runs, inds, shifts)):
        if i == seq_len//2:
            continue
        
        if v >= maxthr:
            rpseq, start, end = expand_repeat(seqa, j, v, sh)
            rplen = end - start
            rpts.append(rpseq)
            rptinds.append((start, rplen))
        
    return rpts, rptinds

def expand_repeat(seq, index, length, shift):
    seq_len = len(seq)
    
    done = False
    i = index+length
    while not done:
        b1 = seq[i]
        if i+shift >= seq_len:
            b2 = seq[i-shift]
        else:
            b2 = seq[i+shift]
        
        if b1 == b2:
            pass
        else:
            done = True
            break
        
        i += 1
        if i >= seq_len:
            break
    ipos = min(i, seq_len)
    
    done = False
    i = index
    while not done:
        
        b1 = seq[i]
        if i+shift >= seq_len:
            b2 = seq[i-shift]
        else:
            b2 = seq[i+shift]
        if b1 == b2:
            pass
        else:
            done = True
            break
        
        i -= 1
        if i <= 0:
            break
    ineg = max(i,0)
    return seq[ineg:ipos], ineg, ipos

def trim_repeats(rpts, inds, seq_len):
    
    rpts = sorted(rpts, key = lambda rpt:len(rpt))
    
    outrpts = []
    outrptinds = []
    i = 0
    j = len(rpts) - 1
    while i < len(rpts):
        rpt1 = rpts[i]
        st1, l1= inds[i]
        outrpts.append(rpt1)
        outrptinds.append((st1, l1))
        while j > i:
            rpt2 = rpts[j]
            st2, l2 = inds[j]
            
            if len(rpt2) <  len(rpt1):
                continue
            
            starts = find_subsequence(rpt2, rpt1)
            
            if starts:
                for start in starts:
                    newind = (st2+start, len(rpt1))
                    if newind in outrptinds:
                        continue
                    outrpts.append(rpt2[start:start+len(rpt1)])
                    outrptinds.append(newind)
                    
                del rpts[j]
            j -= 1
        i +=1 
    return outrpts, outrptinds