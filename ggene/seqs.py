
from typing import Optional

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