
from .vocab import VOCAB

def get_seq_index(seq):
    nbs = 4
    seq_len = len(seq)
    base_idx = {VOCAB[i]:i for i in range(nbs)}
    return sum([base_idx.get(s)*nbs**(seq_len-i-1) for i,s in enumerate(seq)])

def index_to_seq(ind, seq_len = 4):
    nbs = 4
    if ind == 0:
        return "A"*seq_len
    inds = [(ind//nbs**k)%nbs for k in range(seq_len-1, -1, -1)]
    return "".join(VOCAB[i] for i in inds)

def get_seq_index_abs(seq):
    seq = "A" + seq
    seq_len = len(seq)
    fullind = get_seq_index(seq)
    return fullind + sum([4**n for n in range(1, seq_len-1)])

def index_to_seq_abs(ind):
    seq_len = 0
    ip = ind+1
    while ip > 0:
        seq_len += 1
        ip -= 4**seq_len
    ipp = ip + 4**seq_len - 1
    seq = index_to_seq(ipp, seq_len = seq_len+1)
    if seq:
        return seq[1:]
    
    else:
        return ""
    
def permute(seq, n):
    return seq[n:] + seq[:n]
    
def find_min_permutation(seq):
    
    min_ind = get_seq_index(seq)
    min_perm = 0
    min_seq = seq
    
    for i in range(len(seq)):
        ns = permute(seq, i)
        ind = get_seq_index(ns)
        if ind < min_ind:
            min_ind = ind
            min_perm = i
            min_seq = ns
        elif ind == min_ind:
            if ns != min_seq:
                print(f"maybe tiebreaker... {ns} vs {min_seq}")
    return min_perm, min_ind, min_seq