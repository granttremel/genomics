

import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

from scipy.stats import poisson
from tabulate import tabulate

import numpy as np

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import bio, process, find, align, heal, vocab
from ggene.seqs.bio import ALIASES, ALIASES_REV

subs_map = {
    'A': 'TGC',  # not A
    'C': 'ATG',  # not C
    'G': 'ATC',  # not G
    'T': 'AGC',  # not T
}

def get_random_seq(seq_len = 256):
    return bio.get_random_sequence(seq_len)

def get_test_seqs(ind, with_context = False):
    
    if ind == -1:
        seqa, seqb = "CCCCCTGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGG", "CCCCTGGACACCCCCTTCCTTTGTTGGGCCCGGGGGAGGGTGGGCGGGGG"
    
    elif ind == 0:
        seqa, seqb = "ATGC", "TGCC"
    
    elif ind == 1:
        seqa, seqb = "TGGGGAGAAACAAGGGAAGGGAAGAAAAGGGACCGACCC", "TGTGGGGCATTTGGGGGTGTGGGGAAAAAGACCAGAGTA"
    
    elif ind == 2:
        seqa, seqb = "TGGGGAAGCCGGGCCGCTGTGAGCCTGGGAGACAGGAGATGCTATCAGGAGGG", "TGGAAATGTGGGGCATTTGGGGGTGTGGGGAAAAAGACCAGAGTATGGGACAG"
    
    elif ind == 3:
        seqa, seqb = "CAGACCCCAGCGTGGGGAAGCCGGGCCGCTGTGAGCCTGGGAGACA", "AGCTCCCCAGGCTGGAAATGTGGGGCATTTGGGGGTGTGGGGAAAA"
    
    elif ind == 4:
        seqa, seqb = "TTTTGAAATGACCTAATTATCTAAGA", "TTTTTGAAAACACTGAATTTGTAAAA"
    
    return seqa, seqb

def load_genome():
    
    return GenomeManager()


def get_mutated_seq(seq, num_subs = 0, num_ins = 0, ins_len = 0, num_dels = 0, del_len = 0, buffer = 32):
    
    changes = []
    
    seq_len = len(seq)
    seq_list = list(seq)
    
    for nsub in range(num_subs):
        
        spos = random.randint(buffer, seq_len - 1 - buffer)
        old_base = seq[spos]
        new_base = random.choice(subs_map[old_base])
        seq_list[spos] = new_base
        changes.append((spos, old_base, new_base))
        
    if num_dels > 0:
        
        for ndels in range(num_dels):
            
            dlen = poisson.rvs(del_len)
            dpos = random.randint(buffer, seq_len - 1 - buffer)
            del_seq = seq[dpos:dpos + dlen]
            
            del seq_list[dpos:dpos + dlen]
            changes.append((dpos, del_seq, ""))
            
            seq_len = len(seq_list)
            
    if num_ins > 0:
        for nins in range(num_ins):
            
            ilen = poisson.rvs(ins_len)
            ipos = random.randint(buffer, seq_len - 1 - buffer)
            
            ins_seq = vocab.get_random_sequence(ilen)
            seq_list.insert(ipos, ins_seq)
            changes.append((ipos, "", ins_seq))
            
    return "".join(seq_list), changes

def get_random_mutations(init_seq_len, num_subs, num_ins, num_dels, ins_len, del_len, buffer):
    
    mqueue = ["s"]*num_subs + ["i"]*num_ins + ["d"]*num_dels
    random.shuffle(mqueue)
    num_mutations = len(mqueue)
    mutations = []
    seq_len = init_seq_len - 2*buffer
    
    mi = 0
    while len(mutations) < num_mutations:
        
        mtype = mqueue[mi]
        mpos = random.randint(0, seq_len-1)
        mlen = 0
        
        if mtype == "i":
            mlen = poisson.rvs(ins_len)
        elif mtype == "d":
            mlen = -poisson.rvs(del_len)
            if mpos > seq_len + mlen:
                mpos = seq_len + mlen - 1
            
            if mpos < 0:
                continue
        
        mutations.append((mtype, mpos + buffer, mlen))
        seq_len += mlen
        mi += 1
    
    return mutations

def mutate_sequence_stepwise(seq, mutations, buffer = 16):
    
    seq_list = list(seq)
    seq_len = len(seq_list)
    changes = []
    
    for (mtype, mpos, mlen) in mutations:
        
        if mpos < buffer or mpos < seq_len - buffer:
            mpos = (mpos - buffer) % (seq_len - 2*buffer) + buffer
        
        if mtype == "s":
            old_base = seq[mpos]
            new_base = random.choice(subs_map[old_base])
            seq_list[mpos] = new_base
            changes.append((mpos, old_base, new_base))
        elif mtype == "i":
            ins_seq = vocab.get_random_sequence(mlen)
            for b in reversed(ins_seq):
                seq_list.insert(mpos, b)
            changes.append((mpos, "", ins_seq))
        elif mtype == "d":
            if mpos > seq_len + mlen - buffer:
                mpos = seq_len - buffer + mlen - 1
            old_seq = "".join(seq_list[mpos:mpos-mlen])
            del seq_list[mpos:mpos-mlen]
            changes.append((mpos, old_seq, ""))
        
        seq_len = len(seq_list)
    
    return "".join(seq_list), changes


def test_align():
    
    # seqa, seqb = get_test_seqs(0)
    
    base = 8
    seq_len = 2**base
    num_subs = base
    num_ins = base//2
    num_dels = base//2
    ins_len = base
    del_len = base
    
    buffer = seq_len // 16
    seqa = get_random_seq(seq_len)
    muts = get_random_mutations(len(seqa), num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    seqb, changes = mutate_sequence_stepwise(seqa, muts, buffer=buffer)
    
    print(seqa)
    print(seqb)
    
    score = align.score_sequences(seqa, seqb)
    print(score)
    
    alignment = align.align_sequences(seqa, seqb)[0]
    print(alignment)
    
    # print(dir(alignment))
    
    

def main():
    
    # gm = load_genome()
    
    test_align()
    
        
    pass

if __name__=="__main__":
    main()
