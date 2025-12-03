
import random
import numpy as np
from scipy.stats import poisson

import string

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import vocab, bio, process, find
from ggene.seqs.align import IndexMap

subs_map = {
    'A': 'TGC',  # not A
    'C': 'ATG',  # not C
    'G': 'ATC',  # not G
    'T': 'AGC',  # not T
}

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

def build_index_map(seq_len, changes):
    
    ind_map = IndexMap(seq_len)
    
    for p, old, new in changes:
        delta = len(new) - len(old)
        ind_map.add_delta(p, delta)
    
    return ind_map

def test_mutations_ins():    
    seq_len = 32
    v=4
    init_seq = bio.get_random_sequence(seq_len)
    num_subs = 0
    num_ins = 2
    num_dels = 0
    ins_len = 3
    del_len = 3
    buffer = 4
    
    new_seq, changes, ind_map = test_mutations(init_seq, num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    delta = len(changes[0][2]) - len(changes[0][1])
    
    eval_mutations(init_seq, new_seq, ind_map, fname_base = "in")
    
    return new_seq, changes, ind_map

def test_mutations_del():    
    seq_len = 32
    v=4
    init_seq = bio.get_random_sequence(seq_len)
    num_subs = 0
    num_ins = 0
    num_dels = 2
    ins_len = 3
    del_len = 3
    buffer = 4
    
    new_seq, changes, ind_map = test_mutations(init_seq, num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    delta = len(changes[0][2]) - len(changes[0][1])
    
    eval_mutations(init_seq, new_seq, ind_map, fname_base = "del")
    
    return new_seq, changes, ind_map

def test_mutations_indel():    
    seq_len = 32
    v=4
    init_seq = bio.get_random_sequence(seq_len)
    num_subs = 0
    num_ins = 1
    num_dels = 1
    ins_len = 3
    del_len = 3
    buffer = 4
    
    new_seq, changes, ind_map = test_mutations(init_seq, num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    delta = len(changes[0][2]) - len(changes[0][1])
    
    eval_mutations(init_seq, new_seq, ind_map, fname_base = "indel")
    
    return new_seq, changes, ind_map

def eval_mutations(init_seq, new_seq, ind_map, fname_base = "indel"):
    
    print(f"global max: {ind_map.global_seq_len}")
    print(init_seq, len(init_seq))
    print(new_seq, len(new_seq))
    
    print("unmutated mutant")
    unmut = ind_map.unmutate(new_seq)
    print(unmut, unmut.replace('-','') == init_seq)
    
    print("mutated WT")
    wtmut = ind_map.mutate(init_seq)
    print(wtmut, wtmut.replace('-','') == new_seq)
    
    
    print(ind_map.pdls)
    ind_map.plot(fname = f"{fname_base}_test.png")
    ind_map.plot(fname = f"{fname_base}_test2.png", xdata = "m", ydata = ["i"])
    
    pass

def test_mutations(init_seq, num_subs, num_ins, num_dels, ins_len, del_len, buffer):
    
    seq_len = len(init_seq)
    
    muts = get_random_mutations(seq_len, num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    new_seq, changes = mutate_sequence_stepwise(init_seq, muts, buffer=buffer)
    
    for p, old, new in changes:
        oldstr = old
        newstr = new
        if not old:
            oldstr = f"[{'-'*len(new)}]"
        if not new:
            newstr = f"[{'-'*len(old)}]"
        print(f"{p}: {oldstr} -> {newstr}")
    
    ind_map = build_index_map(seq_len, changes)
    
    return new_seq, changes, ind_map
    
    

def main():
    
    # test_mutations_ins()
    # test_mutations_del()
    test_mutations_indel()

    pass

if __name__=="__main__":
    main()

