

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import is_consensus, reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons, hamming
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE, ALIASES, ALIASES_FULL

import numpy as np

from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.motifs.motif import MotifDetector

from ggene.ngrams.base import Base, Sequence
from ggene.ngrams.ngramtree import NGramTree

import random
from ggene.genomemanager import GenomeManager

    
def load_genome():
    return GenomeManager()

def load_test_seq(gm, ind):
    test_seq_locs = [
        ("X", 122741424, 122742161),
        ("X", 122710280, 122710690),
    ]
    test_seq = gm.get_sequence(*test_seq_locs[ind])
    return test_seq

def get_ngrams(seq, n):
    
    ngs = {}
    
    for i in range(len(seq) - n):
        ng = seq[i:i+n]
        if not ng in ngs:
            ngs[ng] = 0
        ngs[ng] += 1
        
    return ngs

def get_digrams(seq):
    return get_ngrams(seq, 2)

def get_trigrams(seq):
    return get_ngrams(seq, 3)

def get_4grams(seq):
    return get_ngrams(seq, 4)

def get_ngram_transitions(seq, n, ns, la = 0):
    """
    n: length of reference n-gram
    ns: length of subsequent n-gram
    la: lookahead
    
    e.g.:
    n=3, ns=1, la=0
    ACTGAG -> (ACT, G) -> ACTG
    
    n=3, ns=2, la=1
    ACTGAG -> (ACT, AG) -> ACTNAG
    
    n=2, ns=2, la=2
    ACTGAG -> (AC, AG) -> ACNNAG
    
    idk
    """
    if not la:
        la = n
    
    txs = {}
    
    for i in range(len(seq)-n-ns-la):
        
        ng = seq[i:i+n].isolate()
        tx = seq[i+la:i+la+ns].isolate()
        
        if not ng in txs:
            txs[ng] = {}
        if not tx in txs:
            txs[tx] = {}
        if not tx in txs[ng]:
            txs[ng][tx] = 0
        
        txs[ng][tx] +=1
    return txs

def order_seq_dist(seq_dist):
    
    ordered = [(k, v) for k, v in seq_dist.items()]
    ordered = sorted(ordered, key = lambda kv:bio.get_seq_index(kv[0]))
    return [k for k, v in ordered], [v for k, v in ordered]

def order_tx_dist(tx_dist):
    ordered = [(k, v) for k, v in tx_dist.items()]
    ordered = sorted(ordered, key = lambda kv:bio.get_seq_index("".join(kv[0])))
    return [k for k, v in ordered], [v for k, v in ordered]

def get_nested_tx_dist(tx_dist, n, ns):
    tx_nst = {}
    for (tx,ntx), f in tx_dist.items():
        if not tx in tx_nst:
            tx_nst[tx] = {}
        tx_nst[tx][ntx] = f
    return tx_nst

def normalize_tx_dist(tx_dist, n, ns):
    
    tx_dist_norm = {}
    tx_dist_freqs = {}
    
    for nga in tx_dist:
        fnorm = sum(tx_dist[nga].values())
        tx_dist_norm[nga] = {ngs:ct/fnorm for ngs, ct in tx_dist[nga].items()}
        tx_dist_freqs[nga] = fnorm
    return tx_dist_norm, tx_dist_freqs

def draw_tx_dist(tx_dist, n, ns, thresh = None):
    
    next_seqs = [bio.index_to_seq(i, seq_len = ns) for i in range(4**ns)]
    # print(next_seqs)
    
    lbls = ["".join([next_seqs[nb][nns] for nb in range(4**ns)]) for nns in range(ns)]
    lblmarg = ["", ""]
    
    cell_frm = "{:<5}{:<10}"
    
    for tx1 in tx_dist:
        data = [tx_dist[tx1].get(ns,0) for ns in next_seqs]
        if thresh and max(data)< thresh:
            continue
        
        bars = draw.scalar_to_text_nb(data, minval = 0, bit_depth= 16, add_range = True)
        
        lmarg = [f"{tx1}",""] + lblmarg
        for m, b in zip(lmarg, bars+lbls):
            print(cell_frm.format(m, b))
        print()

def generate_sequence(tx_dist, n, ns, seq_len, seed = ""):
    
    if not seed:
        seed = "".join(ORDER[random.randint(0, 3)] for nn in range(n))
    
    gen_seq = seed
    
    while len(gen_seq) < seq_len:
        
        ng = gen_seq[-1-n:-1]
        next_dist = tx_dist.get(ng, {})
        rs = random.random()
        vv = 0
        nb = ""
        for k, v in next_dist.items():
            vv+=v
            if vv > rs:
                nb = k
                break
        
        if not nb:
            nb = ORDER[random.randint(0,3)]
            print("flop")
        gen_seq = gen_seq + nb
    
    return gen_seq
    
def calc_entropy(data):
    return sum([-d*np.log(d) for d in data if d > 0])

def find_plausible_mutations(tx_dist, txfqs, n, ns, min_obs =10, max_dist = 1.0):
    
    mdict = {}
    
    ref_ord = sorted(tx_dist.keys(), key= lambda k:-txfqs[k])
    
    for sref in ref_ord:
        msrdict = {}
        
        srdict = tx_dist[sref]
        
        if txfqs[sref] < min_obs:
            continue
        
        p_marg = 0
        sq_ord = sorted(tx_dist[sref].keys(), key = lambda k:-tx_dist[sref][k])
        for isq in range(len(sq_ord)):
            ssq = sq_ord[isq]
            
            for issq2 in range(isq, len(sq_ord)):
                ssq2 = sq_ord[issq2]
                if hamming(ssq, ssq2) < max_dist:
                    cons = bio.merge(ssq, ssq2)
                    # msrdict[ssq] = srdict[ssq]
                    if not cons in msrdict:
                        msrdict[cons] = {}
                    msrdict[cons][ssq] = srdict[ssq]
                    msrdict[cons][ssq2] = srdict[ssq2]
                    
                    
            p_marg += srdict[ssq]
            if p_marg > 0.8:
                break
        
        if msrdict:
            mdict[sref] = msrdict
            
    return mdict

def construct_consensus(seq, mutation_dict, n, ns):
    
    
    
    pass

def unmutate(seq, mutation_dict, n, ns):
    
    outseq = [seq[:n]]
    i=0
    
    while i < len(seq):
        
        sref = seq[i:i+n]
        sref_m = outseq[-1]
        
        ssq = seq[i+n:i+n+ns]
        
        
        if sref in mutation_dict:
            srdict = mutation_dict[sref]
            
            done = False
            for cons, sqdict in srdict.items():
                
                if bio.is_consensus(ssq, cons):
                    maxsq = ""
                    maxpsq = 0
                    for _ssq, _psq in sqdict.items():
                        if _psq > maxpsq:
                            maxpsq = _psq
                            maxsq = _ssq
                    outseq.append(maxsq)
                    i += n
                    done = True
                else:
                    continue
            if not done:
                i+=1
            
        else:
            outseq.append(ssq)
            i += n
    
    return "".join(outseq)
                    
def test_seq_cls():
    gm = load_genome()
    test_seq = load_test_seq(gm, 0)
    rc_test_seq = reverse_complement(test_seq)
    
    seq = Sequence.from_string(test_seq)
    
    for b in seq[50:60]:
        print(repr(b),str(b), b.alias)
    
    subseq = seq[20:30]
    print(subseq)
    print(str(subseq))
    print(repr(subseq))
    
    subseq[3].update(A=2)
    
    print(subseq.crepr())
    
    print(subseq.to_consensus())
    
    rcseq = seq.RC
    print(rcseq)
    print(type(seq), type(subseq), type(rcseq))
    b1 = Base(fdict = {"A":2,"T":4,"C":1})
    b2 = Base(fdict = {"A":2,"T":4,"C":1})
    b3 = Base(fdict = {"A":3,"T":9})
    b4 = Base(fdict = {"A":1,"T":3,"C":3,"G":6})
    b5 = Base(fdict = {"G":12,"C":4})
    
    bs = [b1, b2, b3, b4, b5]
    
    for i in range(len(bs)):
        bi = bs[i]
        print(f"base {i}: {repr(bi)} entropy = {bi.entropy:0.3f}")
    print()
    
        
    for i in range(len(bs)):
        bi = bs[i]
        for j in range(i+1, len(bs)):
            bj = bs[j]
            print(f"base {i} vs {j}")
            print(f"joint entropy: {bi.joint_entropy(bj):0.3f}")    
            print(f"mutual information: {bi.mutual_information(bj):0.3f}")
            print(f"KL divergence ({i},{j}): {bi.KL_div(bj)}")
            # print(f"KL divergence ({j},{i}): {bj.KL_div(bi)}")
            print(f"symmetric KL divergence: {bj.KL_div_symm(bi)}")
            print(f"JS divergence: {bj.JS_div(bi)}")
            print()
    
    b12 = b1.combine(b2).combine(b3).combine(b3).combine(b3)
    print(repr(b12))
    
    b12.trim(p_min = 0.06, keep_obs = False)
    print(repr(b12))
        
    print(b12.mode)
    
    print("".join([b12.sample() for i in range(20)]))
    
    seq2 = Sequence(bases = bs)
    
    print(seq2.crepr())
    print(seq2.RC.crepr())
    
    print(seq2.consensus)
    print(seq2.RC.consensus)
        
    seq3 = seq2.merge(seq2.RC)
    print(seq3.crepr())
    print(seq3.consensus)
    
    for i in range(20):
        seqa = Sequence(bases = [ORDER[random.randint(0,3)] for j in range(3)])
        seqaa = Sequence(bases = [ORDER[random.randint(0,3)] for j in range(3)])
        
        seqb = Sequence(bases = [ORDER[random.randint(0,3)] for j in range(3)])
        seqbb = Sequence(bases = [ORDER[random.randint(0,3)] for j in range(3)])
        
        seqa = seqa.merge(seqaa).merge(seqbb)
        seqb = seqb.merge(seqbb)
        
        seqa = seqa.merge(seqb)
        # seqb = seqb.merge(seqa)
        
        seqa = seqa.merge(seqaa).merge(seqaa)
        seqa.trim(0.1)
        
        # seqb = seqa.copy()
        
        print(seqa.crepr(), seqa.consensus)
        print(seqb.crepr(), seqb.consensus)
        print()
        print(f"alias: {seqa.hamming(seqb, mode = "alias")}")
        print(f"mode: {seqa.hamming(seqb, mode = "mode")}")
        print(f"kl_div: {seqa.hamming(seqb, mode = "kl_div"):0.3f}")
        print(f"js_div: {seqa.hamming(seqb, mode = "js_div"):0.3f}")
        print()
    
    print(hash(seq[1:4]))
    print(hash(seq[2:5]))
    print(hash(seq[4:7]))
    
    mydict = {seq[1:4]:1}
    print(mydict)
    
    print(id(seq))
    seq = seq.copy()

    print(id(seq))
    ss1 = seq[2:5]
    ss2 = seq.slice(2,5)
    ss3 = seq.slice(2,5, keep_ref = True)
    
    print(id(seq), id(ss1), id(ss2), id(ss3),id(None))
    print(id(ss1.bases[0].parent), id(ss2.bases[0].parent), id(ss3.bases[0].parent))
    
    print()

def main():
    
    bg, fg = draw.get_color_scheme("test")

    gm = load_genome()

    # test_seq = load_test_seq(gm, 0)
    # seq = Sequence.from_string(test_seq)

    # n = 3
    # ns = 3
    # la = 0

    # ngt = NGramTree(n, ns, la=la)
    
    # ngt.from_sequence(seq)
    
    # ngs = ngt.ngs
    # ngtx = ngt.txs
    
    # ngt.print()
    
    cons = Sequence.from_string("TRT")
    
    print(Sequence.from_string("TAT").fit_consensus(cons))
    print(Sequence.from_string("TTT").fit_consensus(cons))
    print(Sequence.from_string("TGT").fit_consensus(cons))
    print(Sequence.from_string("TCT").fit_consensus(cons))
    
    # ngt.generalize_ngram(cons)
    
    # print("after generalization:")
    
    # ngt.print()
    

        
        
if __name__=="__main__":
    main()
