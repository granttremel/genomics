

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import is_consensus, reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons, hamming
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE

import numpy as np

from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.motifs.motif import MotifDetector

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
        
        ng = seq[i:i+n]
        tx = seq[i+la:i+la+ns]
        
        if not (ng, tx) in txs:
            txs[(ng, tx)] = 0
        txs[(ng, tx)] +=1
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

def find_plausible_mutations(tx_dist, txfqs, n, ns):
    
    mdict = {}
    
    ref_ord = sorted(tx_dist.keys(), key= lambda k:-txfqs[k])
    
    for sref in ref_ord:
        msrdict = {}
        
        srdict = tx_dist[sref]
        
        if txfqs[sref] < 10:
            continue
        
        p_marg = 0
        sq_ord = sorted(tx_dist[sref].keys(), key = lambda k:-tx_dist[sref][k])
        for isq in range(len(sq_ord)):
            ssq = sq_ord[isq]
            
            for issq2 in range(isq, len(sq_ord)):
                ssq2 = sq_ord[issq2]
                if hamming(ssq, ssq2) == 1:
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
                    
            

def main():
    
    bg, fg = draw.get_color_scheme("test")

    gm = load_genome()

    seq = load_test_seq(gm, 0)
    rcseq = reverse_complement(seq)
    
    for n in range(len(seq)//256 + 1):
        print(seq[n*256:(n+1)*256])
    print()
    
    for n in range(len(seq)//256 + 1):
        print(rcseq[n*256:(n+1)*256])
    print()
    # for n in [2,3,4]:
    #     tgs = get_ngrams(seq, n)
    #     tris, freqs = order_seq_dist(tgs)
    #     for t, f, in zip(tris, freqs):
    #         print(t, f)
    #     print()
    
    ng = 3
    ns = 3
    la = ng
    
    num_cats = 4**ns
    logcats = np.log(num_cats)
    
    tgtx1s = get_ngram_transitions(seq, ng, ns, la = la)
    txnst = get_nested_tx_dist(tgtx1s, ng, ns)
    txnst, txfqs = normalize_tx_dist(txnst, ng, ns)
    num_obs = sum(txfqs.values())
    
    ref_ord = sorted(txnst.keys(), key= lambda k:-txfqs[k])
    
    for sref in ref_ord:
        
        if txfqs[sref] < 10:
            continue
        
        p_obs = txfqs[sref] / num_obs
        ntr = calc_entropy(txnst[sref].values())
        scaled_ntr = ntr/logcats
        cond_ntr = -p_obs*np.log(p_obs)*ntr
        
        print(f"{sref} {p_obs:0.3f}, {ntr:0.3f}, {cond_ntr:0.3f}")
        
        sq_ord = sorted(txnst[sref].keys(), key = lambda k:-txnst[sref][k])
        
        for ssq in sq_ord:
            # ssq = bio.index_to_seq(isq, seq_len = ns)
            if not ssq in txnst[sref]:
                continue
            comb = sref + ssq
            mk = ""
            if comb[ns:ng+ns] == sref:
                mk = "*"
            
            print(f"  {sref} {ssq}{mk}: {txnst[sref][ssq]:0.3f}")
        print()
    
    mdict = find_plausible_mutations(txnst, txfqs, ng, ns)
    
    for sref in mdict:
        print(sref)
        for cons in mdict[sref]:
            pcons = sum(mdict[sref][cons].values())
            print(f"  {cons}: {pcons:0.3f}")
            for ssq in mdict[sref][cons]:
                print(f"   ⟶ {ssq}: {mdict[sref][cons][ssq]:0.3f}")
    print()
    
    seq_unmut = unmutate(seq, mdict, ng, ns)
    
    print(seq[:256])
    print(seq_unmut[:256])
    
    
    """
    TAT -> ATW (A=0.44, T=0.13)
    
    
    
    
    TTA -> TAK (T=0.9, G=0.07)
    TTA -> TRT (A=0.9, G=0.03)
    
    """
    
    
    
    # new_seq = generate_sequence(txnst, ng, ns, 256, "ATG")
    # print(new_seq)
    
    
    # print(txnst)
    
    # draw_tx_dist(txnst, ng, ns, thresh = 20)
    
    # txs, freqs = order_tx_dist(tgtx1s)
    # for tx, f in zip(txs, freqs):
    #     next_ng = tx[0][1:] + tx[1]
    #     print(f"{tx[0]} ⟶ {next_ng} : {f}")


if __name__=="__main__":
    main()
