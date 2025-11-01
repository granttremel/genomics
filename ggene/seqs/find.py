

from typing import Optional, Dict, List, Tuple, Any
from .bio import reverse_complement, get_aliases
from .vocab import VOCAB


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


def compare_sequences(s1: str, s2: str, err_tol: Optional[int] = None) -> int:
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
        elif a in ia and b in ia[a]:
            pass
        elif b in ia and a in ia[b]:
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


