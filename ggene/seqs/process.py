
import random

import numpy as np
import itertools

from typing import Optional, Dict, List, Tuple, Any
from .bio import reverse_complement, ALIASES, ORDER, ALIASES_REV
from .vocab import VOCAB

def rotate_seq(seq, nr = 1):
    nr = nr % len(seq)
    return seq[nr:] + seq[:nr]

def permute_seq(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

def pad_sequence(seq:str, to_length, char = "N"):
    
    rfill = lfill = (to_length - len(seq))//2
    
    if rfill + lfill + len(seq) < to_length:
        rfill += 1
        
    return char * lfill + seq + char*rfill

def tile_sequence(seq, to_length):
    
    rfill = lfill = (to_length - len(seq))//2
    
    if rfill + lfill + len(seq) < to_length:
        rfill += 1
        
    left = lfill//len(seq) + 1
    right = rfill//len(seq) + 1
    
    leftseq = (seq*left)[-lfill:]
    rightseq = (seq*right)[:rfill]
    
    return leftseq + seq + rightseq

def stretch_sequence(seq, step, keep, frame = 0, to_length = None, fill_char = "N"):
    
    fill = step - keep
    out_seq = []
    
    frame = frame % step
    
    for i in range(frame, len(seq), keep):
        out_seq.append(seq[i:i+keep])
        out_seq.append(fill_char * fill)
    
    out_seq = "".join(out_seq)
    if to_length:
        rtrim = ltrim = ( -to_length + len(out_seq))//2
        if not rtrim + len(out_seq) + ltrim == to_length:
            ltrim += 0
        
        out_seq = out_seq[ltrim:len(out_seq) - rtrim]
    
    return out_seq

def crop_sequence(seq, to_length, mode = "center"):
    
    seq_len = len(seq)
    rhalf = (seq_len - to_length)//2
    
    if rhalf < 0:
        return seq
    
    if mode == "center":
        offset = rhalf
    elif mode == "left":
        offset = 0
    elif mode == "right":
        offset = 2*rhalf
    elif mode == "random":
        offset = random.randint(0, 2*rhalf)
    
    return seq[offset:offset+to_length]

def get_comp_entropy(seq):
    
    ent = 0
    for b in "ATGC":
        p = seq.count(b) / len(seq)
        if p > 0:
            ent += -p*np.log2(p)
    
    return ent

def get_ngram_entropy(seq, n):
    
    ngs = {}
    
    for s in range(len(seq) - n):
        
        ng = seq[s:s+n]
        if not ng in ngs:
            ngs[ng] = 0
        ngs[ng] += 1

    # num_ngs = 4**n
    num_ngs = (len(seq) - n)
    ent = 0
    for ng, num_ng in ngs.items():
        
        p_ng = num_ng / num_ngs
        if p_ng > 0:
            ent += -np.log2(p_ng) * p_ng
    
    return ent

def get_ngrams_entropy(seq, max_n):
    
    ns = list(range(1, max_n + 1))
    nds = {n:{} for n in ns}
    
    for s in range(len(seq)):
        
        for n in ns:
            if s+n > len(seq):
                continue
            
            nss = seq[s:s+n]
            if not nss in nds[n]:
                nds[n][nss] = 0
            nds[n][nss] += 1
            
    ent = 0
    for n in ns:
        num_ngs = len(seq) - n
        for ng, num_ng in nds[n].items():
            
            if num_ng > 0:
                
                p = num_ng/num_ngs
                ent += -np.log(p)*p
        
    return ent

def sum_tree(tree, depth = 0):
    
    if not tree:
        return 0, 0
    
    sum_ng = 0
    num_unique = 0
    for b, (ng, subtree) in tree.items():
        sum_ng += ng
        v_ng, v_unique = sum_tree(subtree, depth=depth+1)
        sum_ng += ng + v_ng
        num_unique += 1 + v_unique
        
    return sum_ng, num_unique

def diff_trees(tree1, tree2, depth = 0):
    
    if not tree1 and not tree2:
        return 0, 0
    
    sum_ng = 0
    num_unique = 0    
    
    for b, (ng, subtree) in tree1.items():
        
        if b in tree2:
            ng2, subtree2 = tree2[b]
        else:
            ng2 = 0
            subtree2 = {}
        
        v_ng, v_unique = diff_trees(subtree, subtree2, depth=depth+1)
        
        sum_ng += v_ng + ng - ng2
        num_unique += v_unique + 1 - int(ng2 > 0)
        
    return sum_ng, num_unique        

def place_sample(samp, tree):
    
    new_branch = tree
    
    for b in samp:
        branch = new_branch
        if b in branch:
            nb, new_branch = branch[b]
        else:
            nb, new_branch = 0, {}
            branch[b] = (0, {})
        
        # nb += 1
        branch[b] = (nb, new_branch)
    
    nb, last_branch = branch[b]
    branch[b] = (nb+1, last_branch)
    
    return tree

def get_LZ_complexity(seq, cond = {}):
    
    tree = {}
    
    for i in range(1, len(seq)+1):
        for j in range(0, i):
            sseq = seq[j:i]        
            tree = place_sample(sseq, tree)
    
    if cond:
        num_subseqs, num_unique = diff_trees(tree, cond)
    else:
        num_subseqs, num_unique = sum_tree(tree)
    return tree, num_unique


def get_inter_matrix(seqa, seqb, cmp_func = None, score_func = None):
    
    if not cmp_func:
        cmp_func = lambda a, b: a==b
    
    if not score_func:
        score_func = lambda a, b: 1
    
    mat = np.zeros((len(seqa), len(seqb)))
    rcseqb = reverse_complement(seqb)
    
    for i in range(len(seqa)):
        
        bi = seqa[i]
        
        for j in range(len(seqb)):
            
            bj = seqb[j]
            rcbj = rcseqb[j]
            # rcbj = bio.complement(bj)
            
            if cmp_func(bi, bj):
                mat[i, j] = score_func(bi, bj)
            elif cmp_func(bi, rcbj):
                mat[i,j] = -score_func(bi, rcbj)
            
    return mat

def get_combined_inter_matrix(seqa, seqb):

    fmat = get_full_inter_matrix(seqa, seqb)
    inds = np.array([0,1,2,3])[:, None, None, None]

    return np.sum(fmat*inds, axis=0)


# Precomputed lookup tables for get_full_inter_matrix
# Base indices: N=0, A=1, T=2, G=3, C=4
_BASE_TO_IDX = {'N': 0, 'A': 1, 'T': 2, 'G': 3, 'C': 4}
_IDX_TO_BASE = {0: 'N', 1: 'A', 2: 'T', 3: 'G', 4: 'C'}

def _build_inter_matrix_luts():
    """Build lookup tables for fast inter-matrix computation."""
    # alias_to_ind mapping
    alias_to_ind = {
        "A": 0, "T": 0, "G": 0, "C": 0, "N": 0,
        "R": 1, "Y": -1, "S": 2, "W": -2, "M": 3, "K": -3,
    }
    base_to_ind = {"N": 0, "A": 1, "T": 2, "G": 3, "C": 4}

    # LUTs: [bi_idx, bj_idx] -> (alias_bucket, value)
    alias_bucket_lut = np.zeros((5, 5), dtype=np.int8)
    alias_value_lut = np.zeros((5, 5), dtype=np.float32)

    for bi_idx in range(5):
        bi = _IDX_TO_BASE[bi_idx]
        for bj_idx in range(5):
            bj = _IDX_TO_BASE[bj_idx]

            # Compute alias
            ali = ALIASES_REV.get(bi + bj, ALIASES_REV.get(bj + bi, bi))
            find = alias_to_ind[ali]

            # Compute value
            if find:
                v = 1.0 if find > 0 else -1.0
            else:
                v = float(base_to_ind.get(bi, 0))

            alias_bucket_lut[bi_idx, bj_idx] = abs(find)
            alias_value_lut[bi_idx, bj_idx] = v

    return alias_bucket_lut, alias_value_lut

_ALIAS_BUCKET_LUT, _ALIAS_VALUE_LUT = _build_inter_matrix_luts()


def get_full_inter_matrix(seqa, seqb):
    """
    Compute mutation-type inter-matrix between two sequences.

    Uses precomputed lookup tables and numpy vectorization for performance.
    Returns a (4, 2, len(seqa), len(seqb)) matrix where:
        - axis 0: alias bucket (0=same/N, 1=Y/R, 2=S/W, 3=M/K)
        - axis 1: strand (0=forward, 1=reverse complement)
        - axes 2,3: position in seqa, seqb
    """
    len_a = len(seqa)
    len_b = len(seqb)

    # Convert sequences to index arrays
    seqa_idx = np.array([_BASE_TO_IDX.get(b, 0) for b in seqa], dtype=np.int8)
    seqb_idx = np.array([_BASE_TO_IDX.get(b, 0) for b in seqb], dtype=np.int8)

    # Compute reverse complement indices
    rcseqb = reverse_complement(seqb)
    rcseqb_idx = np.array([_BASE_TO_IDX.get(b, 0) for b in rcseqb], dtype=np.int8)

    # Create 2D index grids for all pairs: (len_a, len_b)
    a_grid = seqa_idx[:, np.newaxis]  # (len_a, 1)
    b_grid = seqb_idx[np.newaxis, :]  # (1, len_b)
    rcb_grid = rcseqb_idx[np.newaxis, :]  # (1, len_b)

    # Look up alias buckets and values for all pairs at once
    fwd_buckets = _ALIAS_BUCKET_LUT[a_grid, b_grid]  # (len_a, len_b)
    fwd_values = _ALIAS_VALUE_LUT[a_grid, b_grid]
    rc_buckets = _ALIAS_BUCKET_LUT[a_grid, rcb_grid]
    rc_values = _ALIAS_VALUE_LUT[a_grid, rcb_grid]

    # Build output using one-hot encoding for bucket dimension
    # Create masks for each bucket and multiply by values
    mat = np.zeros((4, 2, len_a, len_b), dtype=np.float32)

    for bucket in range(4):
        fwd_mask = (fwd_buckets == bucket)
        rc_mask = (rc_buckets == bucket)
        mat[bucket, 0] = np.where(fwd_mask, fwd_values, 0)
        mat[bucket, 1] = np.where(rc_mask, rc_values, 0)

    return mat


def get_full_inter_matrix_slow(seqa, seqb):
    """Original implementation for reference/validation."""
    base_to_ind = {
        "N":0,
        "A":1,
        "T":2,
        "G":3,
        "C":4,
    }

    alias_to_ind = {
        "A":0, "T":0, "G":0, "C":0, "N":0,
        "Y":1,
        "R":-1,
        "S":2,
        "W":-2,
        "M":3,
        "K":-3,
    }

    mat = np.zeros((4, 2, len(seqa), len(seqb)))

    rcseqb = reverse_complement(seqb)

    for i in range(len(seqa)):
        bi = seqa[i]

        for j in range(len(seqb)):
            bj = seqb[j]
            rcbj = rcseqb[j]

            ali = ALIASES_REV.get(bi+bj, ALIASES_REV.get(bj+bi, bi))
            rcali = ALIASES_REV.get(bi+rcbj, ALIASES_REV.get(rcbj+bi, bi))

            find = alias_to_ind[ali]
            rcfind = alias_to_ind[rcali]

            v = find/abs(find) if find else base_to_ind.get(bi,0)
            rcv = rcfind/abs(rcfind) if rcfind else base_to_ind.get(bi, 0)

            mat[abs(find), 0, i, j] = v
            mat[abs(find), 1, i, j] = rcv

    return mat

def shear_matrix(mat, axis = 0, fill = 0):
    
    if axis == 1:
        mat = mat.T
        
    r,c = mat.shape
    bigmat = np.zeros((r+c, c)) + fill
    
    for i in range(r):
        for j in range(c):
            bigmat[c-j-1+i, j] = mat[i, j]
        
    if axis == 1:
        bigmat = bigmat.T
    
    return bigmat

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

def _cf_std(a, b):
    return a==b

def _sf_std(a, b, sc):
    return int(a==b)/sc

def _ff_std(seqa, seqb):
    return 1.0

def _ff_entropy(seqa, seqb):
    
    compa = {b:seqa.count(b)/len(seqa) for b in "ATGC"}
    compb = {b:seqb.count(b)/len(seqb) for b in "ATGC"}
    
    pxy = {(b1, b2):(seqa.count(b1) + seqb.count(b2))/(len(seqa) + len(seqb)) for b1, b2 in itertools.product("ATGC","ATGC")}
    px = compa
    
    # diffs = {}
    # sumdiffs = 0
    # for b in compa.keys():
        
    #     diff = abs(compa[b] - compb[b])
    #     diffs[b] = diff
    #     sumdiffs += diff
    
    # p=diff[b]/sumdiffs is excess p_b in a that is not explained by b, and vice versa 
    # sumdiffs = 1
    
    # entr = sum([-np.log(d/sumdiffs) * d/sumdiffs if d else 0 for d in diffs.values()])
    entr = sum([np.log(pxy[(b1, b2)]/px[b1]) * pxy[(b1, b2)] for b1, b2 in pxy.keys()])
    return entr
    

class CorrFuncs:
    
    @classmethod
    def cf_standard(cls):
        return _cf_std
    
    @classmethod
    def sf_standard(cls):
        return _sf_std
    
    @classmethod
    def ff_standard(cls):
        return _ff_std
    
    @classmethod
    def ff_entropy(cls):
        return _ff_entropy

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


def correlate(seq1, seq2, comparison_func=None, score_func=None, factor_func = None, fill=0, scale=None, step=1, keep=None, shift_step = 1):
    if not comparison_func:
        cmp = CorrFuncs.cf_standard()
    else:
        cmp = comparison_func

    if not score_func:
        # sf = lambda a, b, sc: int(a == b)/sc
        
        sf = CorrFuncs.sf_standard()
    else:
        sf = score_func

    if not factor_func:
        # ff = lambda seqa, seqb: 1.0
        ff = CorrFuncs.ff_standard()
    else:
        ff = factor_func

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
            
            # rcsb = reverse_complement(sb)
            if sa not in VOCAB or sb not in VOCAB:
                continue
            
            score = 0
            if cmp(sa, sb):
                score = sf(sa, sb, sc)
            
            rcscore=0
            if cmp(sa, rcsb):
                rcscore = sf(sa, rcsb, sc)
            
            direct_sum += score
            rc_sum += rcscore
            
        
        fact = ff(seq1_subseq, seq2_subseq)
        rcfact = ff(seq1_subseq, seq2_subseq)
        
        direct_sum *= fact
        rc_sum *= rcfact
        
        # if shift == 0:
        #     print(f"{direct_sum:0.3f}, {rc_sum:0.3f}")
        #     print(seq2_subseq, "seqb")
        #     print(seq1_subseq, "seqa")
        #     print(rcseq2_subseq, "rcseqb")

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
################ deprecated ###############
        
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
