
import random

import numpy as np
import itertools

from typing import Optional, Dict, List, Tuple, Any
from .bio import reverse_complement, get_aliases
from .vocab import VOCAB

def rotate_seq(seq, nr = 1):
    nr = nr % len(seq)
    return seq[nr:] + seq[:nr]

def permute_seq(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

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

def bin_data(data, bin_size = 2, bin_func = None):
    
    if not bin_func:
        bin_func = np.mean
    
    outdata = []
    for i in range(0, len(data), bin_size):
        d = bin_func(data[i:i+bin_size])
        outdata.append(d)
    return outdata
    
def correlate_old(seq1, seq2, comparison_func = None, score_func = None, fill = 0, scale = None):
    
    if not comparison_func:
        cmp = lambda x, y:x==y
    else:
        cmp = comparison_func
    
    if not score_func:
        sf = lambda a, b, sc: int(a==b)/sc
    else:
        sf = score_func
    
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
        seq1t = seq1[s1c-sc//2:s1c+sc//2]
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
                summ += sf(sa, sb, sc)
            elif cmp(sa, rcsb):
                csumm += sf(sa, rcsb, sc)
        
        if t == 0 and fill is not None:
            ips.append(fill)
            comp_ips.append(fill/2)
        else:
            ips.append(summ)
            comp_ips.append(csumm)
        n+=1
    
    return ips, comp_ips

def correlate_longest_subseq_old(seq1, seq2, comparison_func = None, scale = None, fill = None, err_tol = 0):
    
    if not comparison_func:
        cmp = lambda x, y:x==y
    else:
        cmp = comparison_func
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]
    
    start = seq_len // 2
    runs = []
    inds = []
    shifts = []
    found = set()
    
    for i,t in enumerate(range(-start, start+1)):
        
        sslen = seq_len - abs(t)
        
        if scale is None:
            sc = sslen//2
        else:
            sc = scale//2
        
        s1start = max(t, 0)
        s1c = s1start + sslen//2
        seq1t = seq1[s1c-sc:s1c+sc]
        s2start = max(-t, 0)
        s2c = s2start + sslen//2
        seq2t = seq2[s2c-sc:s2c+sc]
        
        max_run = 0
        max_run_end = -1
        max_run_shift = -1
        run = 0
        err = 0
        for j, (sa, sb) in enumerate(zip(seq1t, seq2t)):
            
            if cmp(sa, sb):
                run += 1
            else:
                # run = 0
                err += 1
            
            if err > err_tol:
                run = 0
                err = 0
                
            if run > max_run:
                max_run = run
                max_run_end = j + s1start
                max_run_shift = t
        
        if t == 0 and fill is not None:
            runs.append(fill)
        else:
            runs.append(max_run)
        s1sp = max_run_end - max_run + 1
        s2sp = s1sp - max_run_shift
        if (s1sp, s2sp) in found:
            inds.append(-1)
            shifts.append(-1)
            continue
        else:
            found.add((s1sp, s2sp))
            found.add((s2sp, s1sp))
            
            inds.append(max_run_end - max_run + 1)
            shifts.append(max_run_shift)
    
    return runs, inds, shifts

def convolve_generator(seq1, seq2, scale=None, step=1, keep=None, shift_step = 1):
    if keep is None:
        keep = 1

    if keep > step:
        raise ValueError(f"keep ({keep}) must be <= step ({step})")

    if step < 1 or keep < 1:
        raise ValueError(f"step ({step}) and keep ({keep}) must be >= 1")

    def sample_sequence(seq, step, keep):
        """Sample sequence by keeping `keep` elements every `step` positions."""
        if step == 1:
            return seq  # No sampling needed

        result = []
        for i in range(0, len(seq), step):
            result.append(seq[i:i+keep])

        # Join back into same type as input
        if isinstance(seq, str):
            return ''.join(result)
        else:
            return list(itertools.chain.from_iterable(result))

    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[:seq_len]

    max_shift = seq_len // 2

    for shift in range(-max_shift, max_shift + 1, shift_step):
        # Calculate overlap length at this shift
        overlap_len = seq_len - abs(shift)

        # Calculate start positions in each sequence
        s1_start = max(shift, 0)
        s2_start = max(-shift, 0)

        # Determine how much to extract before sampling
        if scale is None:
            # Use the full overlap
            target_sampled_len = overlap_len
        else:
            # We want `scale` bases after sampling
            target_sampled_len = scale

        # Calculate how many bases we need to extract to get target_sampled_len after sampling
        # With step/keep, we keep (keep/step) fraction of bases
        # So we need: extracted_len * (keep/step) >= target_sampled_len
        # Therefore: extracted_len >= target_sampled_len * (step/keep)
        if step == keep:
            # No actual sampling happening (keep all)
            sc = min(target_sampled_len, overlap_len)
        else:
            # Need to extract more to account for sampling
            import math
            needed_before_sampling = math.ceil(target_sampled_len * step / keep)
            sc = min(needed_before_sampling, overlap_len)

        # Calculate center positions of the overlap
        s1_center = s1_start + overlap_len // 2
        s2_center = s2_start + overlap_len // 2

        # Extract centered subsequences of length `sc`
        seq1_subseq = seq1[s1_center - sc//2 : s1_center + sc//2 + (sc % 2)]
        seq2_subseq = seq2[s2_center - sc//2 : s2_center + sc//2 + (sc % 2)]

        # Apply step/keep sampling
        seq1_sampled = sample_sequence(seq1_subseq, step, keep)
        seq2_sampled = sample_sequence(seq2_subseq, step, keep)

        # Apply final scale limit (in case we extracted more than needed)
        if scale is not None:
            seq1_sampled = seq1_sampled[:scale]
            seq2_sampled = seq2_sampled[:scale]

        yield shift, seq1_sampled, seq2_sampled, overlap_len


def correlate_general(seq1, seq2, analysis_func, scale=None, step=1, keep=None, shift_step = 1, fill_at_zero=None):
    # Initialize result lists (we don't know how many yet)
    results = None

    for shift, seq1_subseq, seq2_subseq, overlap_len in convolve_generator(seq1, seq2, scale, step, keep, shift_step):
        # Handle fill at zero shift
        if shift == 0 and fill_at_zero is not None:
            # print("filling at zero")
            values = fill_at_zero
        else:
            # Call the analysis function
            values = analysis_func(shift, seq1_subseq, seq2_subseq, overlap_len)

        # Ensure values is a tuple
        if not isinstance(values, tuple):
            values = (values,)

        # Initialize results on first iteration
        if results is None:
            results = tuple([] for _ in values)

        # Append each value to its corresponding list
        for result_list, value in zip(results, values):
            result_list.append(value)

    # Return tuple of lists (or empty tuple if no iterations)
    return results if results is not None else tuple()


def correlate(seq1, seq2, comparison_func=None, score_func=None, fill=0, scale=None, step=1, keep=None, shift_step = 1):
    if not comparison_func:
        cmp = lambda x, y: x == y
    else:
        cmp = comparison_func

    if not score_func:
        sf = lambda a, b, sc: int(a == b)*len(seq1)/sc
    else:
        sf = score_func

    rcseq2 = reverse_complement(seq2)

    def analyze_match(shift, seq1_subseq, seq2_subseq, overlap_len):
        """Count direct and reverse-complement matches."""
        # Get the corresponding rc subsequence
        s2_start = max(-shift, 0)
        s2_center = s2_start + overlap_len // 2

        if scale is None:
            sc = overlap_len
        else:
            sc = min(scale, overlap_len)

        rcseq2_subseq = rcseq2[s2_center - sc//2 : s2_center + sc//2 + (sc % 2)]

        direct_sum = 0
        rc_sum = 0

        for sa, sb, rcsb in zip(seq1_subseq, seq2_subseq, rcseq2_subseq):
            if sa not in VOCAB or sb not in VOCAB:
                continue

            if cmp(sa, sb):
                score = sf(sa, sb, sc)
                direct_sum += score
                
            elif cmp(sa, rcsb):
                rc_sum += sf(sa, rcsb, sc)

        return direct_sum, rc_sum

    fill_tuple = (fill, fill) if fill is not None else None
    return correlate_general(seq1, seq2, analyze_match, scale=scale, step=step, keep=keep,
                            fill_at_zero=fill_tuple, shift_step=shift_step)


def correlate_longest_subseq(seq1, seq2, comparison_func=None, scale=None, step=1, keep=None, fill=None, err_tol = 0):

    if not comparison_func:
        cmp = lambda x, y: x == y
    else:
        cmp = comparison_func

    found = set()  # Track (s1_pos, s2_pos) pairs to avoid duplicates

    def analyze_longest_run(shift, seq1_subseq, seq2_subseq, overlap_len):
        """Find the longest matching run in the aligned subsequences."""
        max_run = 0
        max_run_end_local = -1
        run = 0

        # Calculate the actual start position in the full sequences
        s1_start = max(shift, 0)
        s1_center = s1_start + overlap_len // 2

        if scale is None:
            sc = overlap_len
        else:
            sc = min(scale, overlap_len)

        # Offset from center
        offset_from_center = sc // 2
        
        # space1 = " "*abs(shift)
        space1 = ""
        space2 = ""
        if shift < 0:
            space2, space1 = space1, space2
        
        err = 0
        for j, (sa, sb) in enumerate(zip(seq1_subseq, seq2_subseq)):
            if cmp(sa, sb):
                run += 1
            else:
                run += 1
                err += 1
            
            if err > err_tol:
                run = 0
                err = 0

            if run > max_run:
                max_run = run
                max_run_end_local = j
        
        # Convert local index to global index
        if max_run > 0:
            max_run_end_global = s1_center - offset_from_center + max_run_end_local
            max_run_start_global = max_run_end_global - max_run + 1
            s2_start_global = max_run_start_global - shift

            # Check if we've seen this run before
            if (max_run_start_global, s2_start_global) in found:
                return max_run, -1, -1

            found.add((max_run_start_global, s2_start_global))
            found.add((s2_start_global, max_run_start_global))

            return max_run, max_run_start_global, shift
        else:
            return 0, -1, -1

    fill_tuple = (fill, -1, -1) if fill is not None else None
    return correlate_general(seq1, seq2, analyze_longest_run, scale=scale, step=step, keep=keep,
                            fill_at_zero=fill_tuple)

def correlate_longest_subseq_err(seq1, seq2, max_err, comparison_func=None, scale=None, step=1, keep=None, fill=None):

    if not comparison_func:
        cmp = lambda x, y: x == y
    else:
        cmp = comparison_func

    found = set()  # Track (s1_pos, s2_pos) pairs to avoid duplicates

    def analyze_longest_run(shift, seq1_subseq, seq2_subseq, overlap_len):
        """Find the longest matching run in the aligned subsequences."""
        max_run = 0
        max_run_end_local = -1
        max_run_err = 0
        run = 0

        # Calculate the actual start position in the full sequences
        s1_start = max(shift, 0)
        s1_center = s1_start + overlap_len // 2

        if scale is None:
            sc = overlap_len
        else:
            sc = min(scale, overlap_len)

        # Offset from center
        offset_from_center = sc // 2
        
        space1 = ""
        space2 = ""
        if shift < 0:
            space2, space1 = space1, space2
        
        err = 0
        for j, (sa, sb) in enumerate(zip(seq1_subseq, seq2_subseq)):
            if cmp(sa, sb):
                run += 1
            else:
                run += 1
                err += 1
            
            if err > max_err:
                run = 0
                err = 0

            if run > max_run:
                max_run = run
                max_run_end_local = j
                max_run_err = err
        
        # Convert local index to global index
        if max_run > 0:
            max_run_end_global = s1_center - offset_from_center + max_run_end_local
            max_run_start_global = max_run_end_global - max_run + 1
            s2_start_global = max_run_start_global - shift

            # Check if we've seen this run before
            if (max_run_start_global, s2_start_global) in found:
                return max_run, -1, -1, max_run_err

            found.add((max_run_start_global, s2_start_global))
            found.add((s2_start_global, max_run_start_global))

            return max_run, max_run_start_global, shift, max_run_err
        else:
            return 0, -1, -1, -1

    fill_tuple = (fill, -1, -1) if fill is not None else None
    return correlate_general(seq1, seq2, analyze_longest_run, scale=scale, step=step, keep=keep,
                            fill_at_zero=fill_tuple)

def correlate_composition(seq, comp_type = "GC", scale=None, step = 1, keep = None, shift_step = 1, fill = 0.5):
    
    def calc_comp(shift, subseq1, subseq2, overlap_len):
        return sum([subseq1.count(b) for b in comp_type]) / len(subseq1)
    
    fill_tuple = (fill,) if fill is not None else None
    res = correlate_general(seq, seq, calc_comp, scale=scale, step=step, keep=keep, shift_step = shift_step, fill_at_zero = fill_tuple)
    return res

def correlate_gc(seq, scale = None, step = 1, keep = None, shift_step = 1, fill = 0.5):
    return correlate_composition(seq, comp_type = "GC", scale=scale, step=step, keep=keep, shift_step = shift_step, fill=fill)

def correlate_ag(seq, scale = None, step = 1, keep = None, shift_step = 1, fill = 0.5):
    return correlate_composition(seq, comp_type = "AG", scale=scale, step=step, keep=keep, shift_step = shift_step, fill=fill)

def correlate_ac(seq, scale = None, step = 1, keep = None, shift_step = 1, fill = 0.5):
    return correlate_composition(seq, comp_type = "AC", scale=scale, step=step, keep=keep, shift_step = shift_step, fill=fill)
    
def extract_run(seqa, seqb, run_len, run_ind, shift, buffer = 0):
    
    amin = run_ind - buffer
    amax = run_ind + run_len + buffer
    apre = ""
    apost = ""
    if amin < 0:
        apre = " "*abs(amin)
        amin = 0
    if amax > len(seqa):
        apost = " "*(amax - len(seqa))
        amax = len(seqa)
        
    subseqa = apre + seqa[amin:amax] + apost
    
    bmin = run_ind - shift - buffer
    bmax = run_ind - shift + run_len + buffer
    
    bpre = ""
    bpost = ""
    if bmin < 0:
        bpre = " "*abs(bmin)
        bmin = 0
    if bmax > len(seqb):
        bpost = " "*(bmax - len(seqb))
        bmax = len(seqb)
        
    subseqb = bpre + seqb[bmin:bmax] + bpost
    
    return subseqa, subseqb

def extract_max_runs(seqa, seqb, runs, inds, shifts, topk, buffer = 0, do_seqb_rc = False):
    
    topk_runs = sorted([(a, b, c) for a, b,c in zip(runs, inds, shifts)], key = lambda k:-k[0])[:topk]
    
    seqs = []
    
    for run, ind, shift in topk_runs:
        subseqa, subseqb = extract_run(seqa, seqb, run, ind, shift, buffer=buffer)
        seqs.append((subseqa, subseqb))
    
    return seqs, topk_runs

def extract_top_correlated(seqa, seqb, corr, topk, scale = None):
    
    seq_len = len(corr)
    half_len = seq_len//2
    if not scale:
        scale = half_len
    
    topk_corrs = sorted([(a - half_len, b) for a, b in enumerate(corr)], key = lambda k:-k[1])[:topk]
    
    seqs = []
    
    for shift, corr in topk_corrs:
        
        overlap_len = seq_len - abs(shift)
        s1_start = max(shift, 0)
        s2_start = max(-shift, 0)
        s1_center = s1_start + overlap_len // 2
        s2_center = s2_start + overlap_len // 2
        sc = min(overlap_len, scale)
        subseqa = seqa[s1_center - sc//2 : s1_center + sc//2 + (sc % 2)]
        subseqb = seqb[s2_center - sc//2 : s2_center + sc//2 + (sc % 2)]
        
        seqs.append((subseqa, subseqb, corr))
    
    return seqs
    
    