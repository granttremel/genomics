
import random
import numpy as np
from scipy.stats import poisson

import string

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import vocab, bio, process, find, align
from ggene.seqs.align import IndexMap

subs_map = {
    'A': 'TGC',  # not A
    'C': 'ATG',  # not C
    'G': 'ATC',  # not G
    'T': 'AGC',  # not T
}

def get_random_seq(seq_len = 256):
    return bio.get_random_sequence(seq_len)

def get_argtopk(data, topk = 3):
    
    d = [(i, dat) for i, dat in enumerate(data)]
    d = list(sorted(d, key = lambda k:-k[1]))
    topd = d[:topk]
    top_arg = [d[0] for d in topd]
    top_data = [d[1] for d in topd]
    return top_arg, top_data

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

def build_index_map(seq_len, changes):
    
    ind_map = IndexMap(seq_len)
    
    for p, old, new in changes:
        delta = len(new) - len(old)
        ind_map.add_delta(p, delta)
    
    return ind_map

def show_mutations(orig_seq, new_seq, changes):
    
    dels = [old for p, old, new in changes if old and not new]
    ins = [new for p, old, new in changes if new and not old]
    
    # draw.highlight_sequences(orig_seq, dels, min_len = 0, show_key = False, suppress = False)
    # draw.highlight_sequences(new_seq, ins, min_len = 0, show_key = False, suppress = False)
    
    pass

def test_mutated_seq():
    
    seq_len = 256
    test_seq = get_random_seq(seq_len = seq_len)
    v = 4
    
    len_delta = 10
    conf = 3 # sigma
    corr_scale = conf*conf*(v-1) + 2*len_delta
    
    mean = np.sqrt(corr_scale) / v
    sd = np.sqrt(v-1) / v
    
    score_func = lambda a, b, sc: int(a==b)/np.sqrt(sc) 
    
    m1, ch1 = get_mutated_seq(test_seq, num_subs = 5)
    m2, ch2 = get_mutated_seq(test_seq, num_ins = 1, ins_len = len_delta)
    m3, ch3 = get_mutated_seq(test_seq, num_dels = 1, del_len = len_delta)
    
    m2_pos = ch2[0][0]
    m2_len = len(ch2[0][2])
    m3_pos = ch3[0][0]
    m3_len = len(ch3[0][1])
    
    res0, _ = process.correlate(test_seq, test_seq, fill = None, score_func = score_func, scale = corr_scale)
    res1, _ = process.correlate(test_seq, m1, fill = None, score_func = score_func, scale = corr_scale)
    res2, _ = process.correlate(test_seq, m2, fill = None, score_func = score_func, scale = corr_scale)
    res3, _ = process.correlate(m3, test_seq, fill = None, score_func = score_func, scale = corr_scale)
    
    res0_std = [(r-mean)/sd for r in res0]
    res1_std = [(r-mean)/sd for r in res1]
    res2_std = [(r-mean)/sd for r in res2]
    res3_std = [(r-mean)/sd for r in res3]
    
    init_pk = (np.sqrt(corr_scale) - mean)/sd
    res1_pk = res1_std[len(m1)//2]
    res2_pk1 = res2_std[seq_len//2]
    res2_pk2 = res2_std[seq_len//2 - m2_len]
    res3_pk1 = res3_std[seq_len//2 - m3_len//2]
    res3_pk2 = res3_std[seq_len//2 - m3_len//2 - m3_len]
    
    print(f"correlation scale {corr_scale}")
    print("initial sequence:")
    print(test_seq)
    ScalarPlot(res0_std, add_range = True, ruler = True, xmin = 0, xmax = seq_len).show()
    print(f"peak: {init_pk:0.3f}")
    print()
    
    print("with substitutions:")
    print(m1)
    print(ch1)
    
    ScalarPlot(res1_std, add_range = True, ruler = True, xmin = 0, xmax = len(res1)).show()
    print(f"peak: {res1_pk:0.3f}")
    print()
    
    print(f"with insertion length {m2_len} at position {m2_pos}:")
    print(m2)
    print(ch2)
    ScalarPlot(res2_std, add_range = True, ruler = True, xmin = 0, xmax = len(res2)).show()
    print(f"peaks: {res2_pk1:0.3f}, {res2_pk2:0.3f}")
    print()
    
    print(f"with deletion length {m3_len} at position {m3_pos}:")
    print(m3)
    print(ch3)
    ScalarPlot(res3_std, add_range = True, ruler = True, xmin = 0, xmax = len(res3)).show()
    print(f"peaks: {res3_pk1:0.3f}, {res3_pk2:0.3f}")
    print()

def test_mutation_metrics():
    
    seq_len = 256
    v=4
    init_seq = get_random_seq(seq_len)
    
    num_mutations = 20
    num_subs = num_mutations//2
    num_ins = num_subs//2
    num_dels = num_ins
    ins_len = 10
    del_len = 10
    buffer = 16
    muts = get_random_mutations(seq_len, num_subs, num_ins, num_dels, ins_len, del_len, buffer)
    
    new_seq, changes = mutate_sequence_stepwise(init_seq, muts, buffer=buffer)
    
    show_mutations(init_seq, new_seq, changes)
    ind_map = build_index_map(seq_len, changes)
    inds = list(range(seq_len))
    mapped_inds = [ind_map(i) for i in inds]
    print(" ".join(map(str, inds)))
    print(" ".join(map(str, mapped_inds)))
    
    for c in changes:
        print(c)
    
    mean = np.sqrt(seq_len) / v
    sd = np.sqrt(v-1) / v
    score_func = lambda a, b, sc: int(a==b)/np.sqrt(sc) 
    res0, _ = process.correlate_v2(init_seq, new_seq, score_func = score_func, fill = None)
    
    res0_std = [(r-mean)/sd for r in res0]
    ScalarPlot(res0, add_range = True, ruler = True, xmin = 0, xmax = len(new_seq)).show()
    
    runs, inds, shifts = process.correlate_longest_subseq_v2(init_seq, new_seq, fill = None)
    longs = [i for i, r in enumerate(runs) if r > 8]
    ScalarPlot(runs, minval = 10, add_range = True, fg_color = 47, bit_depth = 8).show()
    ScalarPlot(runs, minval = 0, maxval = 10, add_range = True).show()
    
    
    # primes = [3, 5, 7, 11, 13, 19, 23, 27, 31, 37]
    # fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    # step_seq = fibs
    # for ist in range(2, len(step_seq)):
    #     step = step_seq[ist]
    #     # for keep in range(step//3, 2*step//3):
    #     for keep in step_seq[ist-2:ist]:
    #         runs, inds, shifts = process.correlate_longest_subseq_v2(init_seq, new_seq, fill = None, step = step, keep = keep)
    #         found_longs = [str(longs.index(i)) for i, r in enumerate(runs) if r > 10]
    #         num_long = len(longs)
    #         print(f"step = {step}, keep = {keep}, r = {keep/step:0.3f} longs = {",".join(found_longs)}")
    #     print()
    # pass


def sequence_similarity(seqa, seqb):
    
    N = len(seqa)
    res = 0
    for a,b in zip(seqa, seqb):
        if a==b:
            res += 1
            
    return res/np.sqrt(N)

def run_seq_similarity_montecarlo(num_tests = 200, seq_len = 256, vocab_len = 4):
    
    sims = []
    
    if vocab_len != 4:
        new_vocab = string.ascii_uppercase[:vocab_len]
        vocab.set_vocab(new_vocab)
    
    ref_seq = vocab.get_random_sequence(seq_len)
    print(f"ref seq: {ref_seq}")
    
    for n in range(num_tests):
        
        test_seq = vocab.get_random_sequence(seq_len)
        sim = sequence_similarity(ref_seq, test_seq)
        sims.append(sim)
    
    mean = np.mean(sims)
    sd = np.std(sims)
    
    print(f"sequence similarity with {num_tests} trials and sequence length {seq_len}:")
    print(f"mean = {mean:0.3f}, sd = {sd:0.3f}")

def test_seqsim_montecarlo():
    
    N = 256
    v = 4

    run_seq_similarity_montecarlo(num_tests = 2000, seq_len = N, vocab_len = v)
    
    meansim = np.sqrt(N) / v
    varsim = (v-1)/v/v
    sdsim = np.sqrt(varsim)
    
    print(f"sequence similarity mean = {meansim:0.3f}, var = {varsim:0.3f}, sd = {sdsim:0.3f}")

def test_convolve_step():
    test_seq = "AAAAAAACTGACTGTGACTGACTTTTTTT"

    for shift, subseq1, subseq2, overlap in process.convolve_generator(test_seq, test_seq, step = 3):
        if shift > 0:
            s1 = " "*abs(shift) + subseq1
            s2 = subseq2
        else:
            s2 = " "*abs(shift) + subseq2
            s1 = subseq1
            
        print(s1)
        print(s2)
        print()

def test_checkpoints():
    
    seq_len = 128
    buffer = 32
    max_err = 5
    
    num_ins = 3
    num_dels = 3
    topk = num_ins + num_dels + 1
    
    seq = get_random_seq(seq_len = seq_len)
    muts = get_random_mutations(seq_len, num_subs = 20, num_ins = num_ins, ins_len = 10, num_dels = num_dels, del_len = 5, buffer = buffer)
    mseq, changes = mutate_sequence_stepwise(seq, muts, buffer=buffer)
    
    new_len = min(len(seq), len(mseq))
    seq = seq[:new_len]
    mseq = mseq[:new_len]
    
    runs, inds, shifts, errs = process.correlate_longest_subseq_err(seq, mseq, fill = None, max_err = max_err)
    ScalarPlot(runs, add_range = True, minval = 0).show()
    
    tops, topdatas = process.extract_max_runs(seq, mseq, runs, inds, shifts, topk, buffer = buffer)
    print(seq)
    print(mseq)
    print()
    
    checkpoints = align.get_align_checkpoints(seq, mseq, num_checkpoints = topk, err_ratio = 0.2)
    
    lmxa = lmxb = None
    for mina, maxa, minb, maxb, entr in checkpoints:
        
        if lmxa is None:
            pass
        else:
            ssaa = seq[lmxa:mina]
            ssbb = mseq[lmxb:minb]
            print("between")
            print(ssaa)
            print(ssbb)
            print()
        
        ssa = seq[mina:maxa]
        ssb = mseq[minb:maxb]
        
        print(f"checkpoint {mina}-{maxa}, {minb}-{maxb}, entropy = {entr:0.3f}")
        print(ssa)
        print(ssb)
        print()
        
        lmxa = maxa
        lmxb = maxb
    
    print("the rest")
    print(seq[lmxa:])
    print(mseq[lmxb:])
        
    
    # print(checkpoints)
    
    # diffs = align.align_sequences_chkpt(seq, mseq, checkpoints)
    
    # seq_algn, mseq_algn = align.fill_aligned(seq, mseq, diffs)
    
    # print(seq_algn)
    # print(mseq_algn)
    

def main():
    
    # test_mutation_metrics()
    # test_mutated_seq()
    test_checkpoints()

    pass

if __name__=="__main__":
    main()



