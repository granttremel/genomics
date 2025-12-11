
import math

from .vocab import VOCAB


def totient(n):
    
    if n==1:
        return 1
    
    phi = int(n>1 and n)
    for p in range(2, int(n**0.5+1)):
        if not n%p: # n,p coprime
            phi -= phi//p
            while not n%p:
                n//=p
    
    if n>1: phi -= phi//n
    return phi

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

def find_min_permutation(seq):
    
    min_ind = get_seq_index(seq)
    min_perm = 0
    min_seq = seq
    
    for i in range(len(seq)):
        ns = rotate(seq, i)
        ind = get_seq_index(ns)
        if ind < min_ind:
            min_ind = ind
            min_perm = i
            min_seq = ns
        elif ind == min_ind:
            if ns != min_seq:
                print(f"maybe tiebreaker... {ns} vs {min_seq}")
    return min_perm, min_ind, min_seq


def get_all_sequences(n, k, ab = None):
    
    
    
    pass


def rotate(seq, n):
    return seq[n:] + seq[:n]

def count_rotations(n, k):
    """
    unique necklaces length n with k symbols, up to rotation. or, unique repeats (over quotient/modulo space)
    """
    N = 0
    for d in range(1, int(n**0.5)+1):
        
        if n%d==0:
            dd = n//d
            
            N += totient(d)*k**(dd)
            if d!=dd:
                N += totient(dd)*k**(d)
    
    return N//n
    
def to_canonical_r(seq, abt = "ATGC"):
    
    symbol_idx = {a:i for a,i in enumerate(abt)}
    # idx_symbol = {i:a for a,i in symbol_idx.items()}
    s = [symbol_idx[b] for b in seq]
    
    n = len(seq)
    f = [-1] * (2*n)
    k = 0
    for j in range(1, 2*n):
        i = f[j-k-1]
        while i != 1 and s[j%n] != s[(k+i+1) % n]:
            if s[j%n] < s[(k+i+1)%n]:
                k = j - i - 1
            i = f[i]
        if i==-1 and s[j%n] != s[(k+i+1)%n]:
            if s[j%n] < s[(k+i+1)%n]:
                k=j
            f[j-k] = -1
        else:
            f[j-k] = i+1
    
    return k, seq[k:] + seq[:k]