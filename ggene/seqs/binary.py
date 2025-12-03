
from typing import List, Dict, Optional, Any, Union
from enum import StrEnum
import random

from ggene.seqs.bio import ALIASES, ALIASES_REV, ORDER

base_inds = {b:i for i, b in enumerate("AGCT")}

def get_onehot_fwd():
    
    ohdict = {"-":[0,0,0,0]}
    
    for bs,alias in ALIASES_REV.items():
        oh = [0,0,0,0]
        for b in bs:
            b_ind = base_inds.get(b)
            oh[b_ind] = 1
        ohdict[alias] = oh
    
    return ohdict

def get_onehot_rev():
    return {tuple(v):k for k,v in get_onehot_fwd().items()}

ONE_HOT = get_onehot_fwd()
ONE_HOT_REV = get_onehot_rev()

class BinaryBase:
    
    def __init__(self, base):
        self.base = base
        self.oh = ONE_HOT.get(base, ONE_HOT['-'])
    
    @classmethod
    def from_onehot(cls, onehot):
        base = ONE_HOT_REV.get(tuple(onehot))
        return cls(base)
    
    def complement(self):
        return BinBase.from_onehot(self.oh[::-1])
    
    def logical_and(self, other:'BinBase'):
        oh_out = []
        
        for v1, v2 in zip(self.oh, other.oh):
            oh_out.append(v1 and v2)
        
        new_bb = BinBase.from_onehot(oh_out)
        
        return new_bb
        
    def logical_or(self, other:'BinBase'):
        oh_out = []
        
        for v1, v2 in zip(self.oh, other.oh):
            oh_out.append(v1 or v2)
        
        return BinBase.from_onehot(oh_out)
    
    def logical_xor(self, other:'BinBase'):
        oh_out = []
        
        for v1, v2 in zip(self.oh, other.oh):
            oh_out.append(v1 ^ v2)
        
        return BinBase.from_onehot(oh_out)
    
    def __repr__(self):
        return f"BinaryBase({self.base})"

    def __str__(self):
        return self.base
        
class BinarySequence:
    
    def __init__(self,seq):
        if isinstance(seq, list):
            seq = "".join(map(str, seq))
        self.seq:str = seq
        self.bbs:List[BinBase] = [BinBase(b) for b in seq]
    
    @property
    def multiplicity(self):
        return sum([len(ALIASES.get(b, "")) for b in self.seq])
    
    @property
    def bit_length(self):
        return len(self)
    
    ########## binary ###############
    
    def __and__(self, other:'BinSeq'):
        return self.logical_and(other)
    
    def __or__(self, other:'BinSeq'):
        return self.logical_or(other)
    
    def __xor__(self, other:'BinSeq'):
        return self.logical_xor(other)
    
    def __rshift__(self, n):
        return self.rshift(n)
    
    def __lshift__(self, n):
        return self.lshift(n)
    
    def logical_and(self, other:'BinSeq'):
        outseq = []
        
        for b1, b2 in zip(self, other):
            outseq.append(b1.logical_and(b2))
        
        return BinSeq(outseq)
        
    def logical_or(self, other:'BinSeq'):
        outseq = []
        
        for b1, b2 in zip(self, other):
            outseq.append(b1.logical_or(b2))
        
        return BinSeq(outseq)

    def logical_xor(self, other:'BinSeq'):
        outseq = []
        
        for b1, b2 in zip(self, other):
            outseq.append(b1.logical_xor(b2))
        
        return BinSeq(outseq)
    
    def rshift(self, n):
        return self.subseq(0, len(self) - n)
    
    def lshift(self, n):
        return BinSeq(self.seq + '-'*n)
    
    def rtrunc(self, pos):
        return self.subseq(len(self) - pos, len(self))
    
    def concat(self, other:'BinSeq'):
        return BinSeq(self.seq + other.seq)

    def double(self):
        newseq = "".join([b + b for b in self.seq])
        return BinSeq(newseq)

    def repeat(self, n):
        return BinSeq(self.seq * n)

    def lfill(self, b, to_len):
        return BinSeq(b*(to_len - len(self)) + self.seq)
    
    def rotate(self, n):
        if n < 0:
            n = len(self) + n
        elif n == 0:
            return BinSeq(self.seq)
        
        return BinSeq(self.seq[n:] + self.seq[:n])
    
    ############## other stuff ######################3
    
    def __getitem__(self, ind):
        return self.bbs[ind]

    def __getitems__(self, indexer):
        subseq = self.seq[indexer]
        newseq = BinSeq(subseq)
        return newseq

    def subseq(self, start = 0, stop = -1, step = 1):
        subseq = self.seq[slice(start, stop, step)]
        return BinSeq(subseq)
    
    def sample(self):
        sampled = [random.choice(ALIASES.get(b,"-")) for b in self.seq]
        return BinSeq(sampled)

    def reverse(self):
        return BinSeq(self.seq[::-1])

    def complement(self):
        return BinSeq([b.complement() for b in self.bbs])
    
    def reverse_complement(self):
        return BinSeq([b.complement() for b in self.bbs[::-1]])
    
    def to_palindrome(self):
        return BinSeq(self.seq + self.reverse().seq)
    
    def to_dyad(self):
        return BinSeq(self.seq + self.reverse_complement().seq)
    
    def __len__(self):
        return len(self.bbs)
    
    def __repr__(self):
        return f"BinSeq({self.seq})"
    
    def __str__(self):
        return self.seq

BinBase = BinaryBase
BinSeq = BinarySequence


def get_const(order = 1, parity = 0, nbits = 64):
    
    if order == 0:
        return BinSeq("N"*nbits)
    
    out = (BinSeq("N")<<1).repeat(nbits//2)
    
    for ord in range(order - 1):
        out = out.rtrunc(nbits//2).double()
        
    if parity:
        out = out >> 2**(order-1)
    
    return out.lfill("-", nbits).rtrunc(nbits)

def get_consts(nbits = 64, max_n = 5):
    
    _consts = {}
    for n in range(0,max_n):
        for p in [0,1]:
            if p and not n:
                continue
            _consts[(n, p)]=get_const(n, p, nbits = nbits)
    
    return _consts

"""
in binary, consts go like:
11111111 (0, 0)
01010101 (1, 0)
10101010 (1, 1)
00110011 (2, 0)
11001100 (2, 1)
and so on

here (ternary?), consts should have more four parities
properties:
- c[o, p] & c[o, p+1] = 0
- c[o, p] | c[o, p+1] = 2^n-1


so...
NNNNNNNN (0,0)

0001000100010001
-A-A-A-A (1,0)

0010001000100010
-G-G-G-G (1,1)
-C-C-C-C (1,2)
-T-T-T-T (1,3)

0001001000010010 = AGAG
0001010000010100 = ACAC
0001100000011000 = ATAT
...
GCGC
GTGT
CTCT

(2,0)=AAGGAAGG
(2,1)=AACCAACC
(2,2)=AATTAATT
(2,3)=GGCCGGCC
(2,4)=GGTTGGTT
(2,5)=CCTTCCTT

(2,0)=AAGGAAGG
(2,1)=AACCAACC
(2,2)=AATTAATT
(2,3)=GGCCGGCC
(2,4)=GGTTGGTT
(2,5)=CCTTCCTT



"""


BINARY_MAP1 = {
    "G": [0b00],
    "T": [0b10],
    "A": [0b01],
    "C": [0b11],
}
"""
properties:
- two bases are complementary iff they sum to 3. the comp of one base is its three's complement.
- purines even, pyrimidines odd
- G is zero and T is identity? so A+-A = G, and A*A^-1 = T. meaningless
- generally symmetric wrt GC and AT. I put G first because it feels more important than A (I like it better)
- also can swap A with T and G with C, though you lose the evenness property
- all bases have const bit depth = 2, no waste
"""

BINARY_MAP2 = {
    "" : [0b00],
    "G": [0b01],
    "T": [0b10],
    "A": [0b11],
    "C": [0b100],
}
"""
properties:
- two bases are complementary iff they sum to 5. the complement of one base is its five's complement.
- purines odd, pyrimidines even
- requires bit_depth=3 to represent, wastes a bunch of space. maybe use to represent methylation or smth? or other exotic (like deaminocytosine,  )
- G is identity..? 
   G T A C
G  G T A C
T  T C G A
A  A G C T
C  C A T G
- yeah not really meaningful?
   G T A C
G  R Y R Y
T  Y Y R R
A  R R Y Y
C  Y R Y R

   G T A C
G  S W W S
T  W S S W
A  W S S W
C  S W W S

   G T A C
G  K K M M
T  K M K M
A  M K M K
C  M M K K
- nah

"""

BINARY_MAP3 = {
    "T":[-0b10],
    "C":[-0b01],
    "" :[0b00],
    "G": [0b01],
    "A": [0b10],
}
"""
properties:
- two bases comp iff sum to zero
- purines positive, pyrimidines negative
- constant width though can't be concatenated into larger bin due to negatives
"""

BINARY_MAP4 = {
    "" : [0b00, 0b00],
    "G": [0b00, 0b01],
    "C": [0b00, 0b10], 
    "A": [0b01, 0b00],
    "T": [0b10, 0b00],
}
"""
properties:
- two bases comp iff sum has one component equal to 2
"""
BINARY_MAP5 = {
    "" : [0b00, 0b00],
    "G": [0b00, 0b01],
    "A": [0b00, 0b10],
    "C": [0b01, 0b00], 
    "T": [0b10, 0b00],
}
"""
- two bases comp iff sum contains identical entries
"""

BINARY_MAP10 = {
    "" : [0b0, 0b0],
    "G": [0b0, 0b1],
    "A": [0b1, 0b0],
    "C": [0b0,-0b1],
    "T": [-0b1,0b0],
}
"""
    A
  C   G
    T
- two bases complementary if they sum to zero, i.e. the null string. 
- G+A = [0b1, 0b1]. so what?
- represent a sequence as two binaries
"""

BINARY_MAP11 = {
    "":  [0b0, 0b0, 0b0, 0b0],
    "G": [0b0, 0b0, 0b0, 0b1],
    "A": [0b0, 0b0, 0b1, 0b0],
    "C": [0b0, 0b1, 0b0, 0b0],
    "T": [0b1, 0b0, 0b0, 0b0],
}
"""
one-hot encoding
properties:
- no algebraic properties to be spoken of
- no ambiguity possible
- one sequence requires four parallel binaries

"""

MAP_INDEX = 0
MAPS = [{"":[0b0]}, BINARY_MAP1, BINARY_MAP2, BINARY_MAP3, BINARY_MAP4, BINARY_MAP5, BINARY_MAP10, BINARY_MAP11]
BINARY_MAP = MAPS[MAP_INDEX]
# MAP_TYPES = [type(m.get("G")) for m in MAPS]
MAP_LENS = [len(m.get("G",[])) for m in MAPS]
MAP_WIDTHS = [max([max(mv) for mv in m.values()]).bit_length() for m in MAPS]

def set_binary_map(ind):
    global MAP_INDEX, BINARY_MAP
    MAP_INDEX = ind
    BINARY_MAP = MAPS[ind]
    return BINARY_MAP

def get_binary_data():
    return MAP_LENS[MAP_INDEX], MAP_WIDTHS[MAP_INDEX]

def to_binary(seq):
    blen, bwidth = get_binary_data()
    bins = [0b0 for i in range(blen)]
    seq_len = len(seq)
    for ib in range(seq_len):
        b = seq[ib]
        base_bin = BINARY_MAP.get(b)
        for nb in range(blen):
            bins[nb] += base_bin[nb] << bwidth*ib
    return bins

def to_sequence(binseq):
    
    
    pass
