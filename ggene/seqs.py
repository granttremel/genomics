
from typing import Optional, Dict, List, Tuple, Any

_comp_strs = "AĀEĒGḠIĪOŌUŪYȲ"
vocab = "AUGC"
_vc = vocab

COMPLEMENT_MAP = {'A':'U','C':'G','G':'C','U':'A'}
_cm = COMPLEMENT_MAP

iupac_aliases = {
    'R': 'AG',   # puRine
    'Y': 'UC',   # pYrimidine  
    'S': 'GC',   # Strong (3 H bonds)
    'W': 'AU',   # Weak (2 H bonds)
    'K': 'UG',   # Keto
    'M': 'AC',   # aMino
    'B': 'UGC',  # not A
    'D': 'AUG',  # not C
    'H': 'AUC',  # not G
    'V': 'AGC',  # not T
    'N': 'AUGC', # aNy
}

iupac_aliases_rev = {
    'A':'RWMDHVN',
    'U':'YWKBDHN',
    'G':'RSKBDVN',
    'C':'YSMBHVN',
}

aliases = {}

def get_iupac_aliases(_vocab = None):
    if vocab:
        ia = {}
        iar = {}
        for k, v in iupac_aliases.items():
            ia[convert_to_vocab(k, _vocab, vocab)] = convert_to_vocab(v, _vocab, vocab)
        for k, v in iupac_aliases_rev.items():
            iar[convert_to_vocab(k, _vocab, vocab)] = convert_to_vocab(v, _vocab, vocab)
        return ia, iar
    else:
        return iupac_aliases, iupac_aliases_rev

def set_vocab(new_vocab, new_cm=None):
    global vocab, COMPLEMENT_MAP
    
    if not new_cm:
        new_cm = make_complement_map(new_vocab)
    
    vocab = new_vocab
    COMPLEMENT_MAP = new_cm

def reset_vocab():
    global vocab, COMPLEMENT_MAP
    vocab = _vc
    COMPLEMENT_MAP = _cm

def convert_to_vocab(seq, to_vocab, from_vocab = None):
    if not from_vocab:
        from_vocab = vocab
    vcdict = {f:t for f,t in zip(from_vocab, to_vocab)}
    return "".join(vcdict[s] for s in seq)

def make_complement_map(vocab):
    cm = {}
    for i in range(0, len(vocab)//2+1, 2):
        cm[vocab[i]] = vocab[i+1] 
        cm[vocab[i+1]] = vocab[i] 
    return cm
    
def reverse_complement(seq):
    return "".join(_reverse(_complement(seq)))

def _complement(seq):
    return [COMPLEMENT_MAP.get(s,s) for s in seq]

def complement(seq):
    return "".join(_complement(seq))

def _reverse(seq):
    return [s for s in reversed(seq)]

def reverse(seq):
    return "".join(_reverse(seq))

def consensus_to_re(seq, err = 0):
    
    ia, iar = get_iupac_aliases()
    outstr = []
    for s in seq:
        if s in vocab:
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
            new_grp_bs = [v for v in vocab if not v in grp]
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

    nerr = 0
    for a, b in zip(s1, s2):
        if a == b:
            pass
        elif a in iupac_aliases and b in iupac_aliases[a]:
            pass
        elif b in iupac_aliases and a in iupac_aliases[b]:
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
    
    # if not max_shift:
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
            
            if not sa in vocab or not sb in vocab:
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

def convolve(seq1, seq2, comparison_func = None, sequence_length = None, fill = 0):
    
    if not comparison_func:
        cmp = lambda x, y:x==y
    else:
        cmp = comparison_func
    
    if not sequence_length:
        sequence_length = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:sequence_length], seq2[:sequence_length]
    rcseq2 = reverse_complement(seq2)
    
    start = sequence_length // 2
    ips = []
    comp_ips = []
    n = 0
    for t in range(-start, start + 1):
        
        sslen = sequence_length - abs(t)
        s1start = max(t, 0)
        seq1t = seq1[s1start:s1start + sslen]
        s2start = max(-t, 0)
        seq2t = seq2[s2start:s2start + sslen]
        rcseq2t = rcseq2[s2start:s2start + sslen]
        
        summ = 0
        csumm = 0
        for sa, sb, rcsb in zip(seq1t, seq2t, rcseq2t):
            
            if cmp(sa, sb):
                summ += 1/sslen
            elif cmp(sa, rcsb):
                csumm += 1/sslen
        
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


def align_sequences_with_gaps(seqa: str, seqb: str, features: List[Dict[str, Any]], window_start = 0):
    """Align reference and personal sequences by inserting gaps at indel positions.
    
    Returns:
        Tuple of (aligned_ref, aligned_pers, variant_positions)
    """
    # Find variant features to identify indels
    variant_features = [f for f in features if f.feature_type == 'variant']
    
    # Create a list of indel events in the window
    indels = []
    
    for var in variant_features:
        var_start = var.start
        ref = var.get('ref', '')
        alt = var.get('alt', '')
        
        # Calculate position in our window (0-based)
        window_pos = var_start - window_start
        
        # Only process if the variant starts within our window
        if 0 <= window_pos < len(seqa):
            if len(ref) != len(alt):
                indels.append({
                    'pos': window_pos,
                    'ref_len': len(ref),
                    'alt_len': len(alt),
                    'ref': ref,
                    'alt': alt
                })
    
    # Sort indels by position
    indels.sort(key=lambda x: x['pos'])
    
    # Build aligned sequences with gaps
    aligned_a = []
    aligned_b = []
    variant_positions = []  # Track where variants are for coloring
    
    idx_a = 0
    idx_b = 0
    
    for i in range(len(seqa)):
        # Check if we're at an indel position
        current_indel = None
        for indel in indels:
            if indel['pos'] == i:
                current_indel = indel
                break
        
        if current_indel:
            ref_len = current_indel['ref_len']
            alt_len = current_indel['alt_len']
            
            if alt_len > ref_len:
                # Insertion - add gaps to reference
                for j in range(ref_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_a += 1
                    idx_b += 1
                
                # Add the inserted bases with gaps in reference
                for j in range(alt_len - ref_len):
                    aligned_a.append('-')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_b += 1
                    
            elif ref_len > alt_len:
                # Deletion - add gaps to personal  
                for j in range(alt_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_a += 1
                    idx_b += 1
                
                # Add the deleted bases with gaps in personal
                for j in range(ref_len - alt_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append('-')
                    variant_positions.append(True)
                    idx_a += 1
                    
            else:
                # SNP - just add both
                aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                variant_positions.append(seqa[idx_a] != seqb[idx_b])
                idx_a += 1
                idx_b += 1
        else:
            # No indel at this position
            if idx_a < len(seqa) and idx_b < len(seqb):
                aligned_a.append(seqa[idx_a])
                aligned_b.append(seqb[idx_b])
                variant_positions.append(seqa[idx_a] != seqb[idx_b])
                idx_a += 1
                idx_b += 1
            elif idx_a < len(seqa):
                aligned_a.append(seqa[idx_a])
                aligned_b.append('-')
                variant_positions.append(True)
                idx_a += 1
            elif idx_b < len(seqb):
                aligned_a.append('-')
                aligned_b.append(seqb[idx_b])
                variant_positions.append(True)
                idx_b += 1
    
    return ''.join(aligned_a), ''.join(aligned_b), variant_positions

