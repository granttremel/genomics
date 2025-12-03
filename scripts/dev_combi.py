

import random

from tabulate import tabulate

from ggene.seqs import combi

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.vocab import VOCAB
from ggene.seqs.bio import is_consensus, reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons, hamming, transition, transvert
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE, ALIASES, ALIASES_FULL

def get_random_seq(seq_len):
    return "".join([VOCAB[random.randint(0,3)] for i in range(seq_len)])

def test_random_perm(seq_len = 32):
    
    seq = get_random_seq(seq_len)
    
    min_ind = float("inf")
    min_perm = 0
    
    inds = []
    for i in range(len(seq)):
        ns = combi.permute(seq, i)
        ind = combi.get_seq_index(ns)
        if ind < min_ind:
            min_ind = ind
            min_perm = i
        inds.append(ind)
    
    print(f"permutation with minimum index: {min_perm} with {min_ind}")
    print(seq)
    sctxt = draw.scalar_to_text_nb(inds, minval = 0, maxval = 4**seq_len, add_range = True, range_fstr = ".2e")
    for r in sctxt:
        print(r)
    print()
    sctxt_prm = draw.scalar_to_text_nb(combi.permute(inds, min_perm), minval = 0, maxval = 4**seq_len, add_range = True, range_fstr = ".2e")
    print(combi.permute(seq, min_perm))
    for r in sctxt_prm:
        print(r)

def find_min_perm(seq):
    
    min_ind = combi.get_seq_index(seq)
    min_perm = 0
    min_seq = seq
    
    for i in range(len(seq)):
        ns = combi.permute(seq, i)
        ind = combi.get_seq_index(ns)
        if ind < min_ind:
            min_ind = ind
            min_perm = i
            min_seq = ns
        elif ind == min_ind:
            if ns != min_seq:
                print(f"maybe tiebreaker... {ns} vs {min_seq}")
    return min_perm, min_ind, min_seq

def generate_members(seq):
    txs = [lambda a:a, reverse_complement, transition, transvert]
    tx_names = ["e","rc","tn","tv"]
    txstrs = []
    mems = []
    for i in range(len(seq)):
        for itx in range(len(txs)):
            ns = combi.permute(txs[itx](seq), i)
            mems.append(ns)
            txstrs.append(f"{tx_names[itx]} p{i}")
    return mems, txstrs

def count_nrepeats(n=3):
    
    all_seqs = [combi.index_to_seq(i, seq_len = n) for i in range(4**n)]
    
    # random.shuffle(all_seqs)
    
    fams = {}
    txs = {}
    
    min_tab = []
    
    for s in all_seqs:
        rcs = reverse_complement(s)
        tns = transition(s)
        tvs = transvert(s)
        
        perm, min_ind, min_seq = find_min_perm(s)
        rcperm, rcmin_ind, rcmin_seq = find_min_perm(rcs)
        tnperm, tnmin_ind, tnmin_seq = find_min_perm(tns)
        tvperm, tvmin_ind, tvmin_seq = find_min_perm(tvs)
        
        all_seqs = [min_seq, rcmin_seq, tnmin_seq, tvmin_seq]
        all_mins = [min_ind, rcmin_ind, tnmin_ind, tvmin_ind]
        min_tab.append([s] + all_mins)
        # if s == "ACG":
            # print(s, "|", min_seq, min_ind, rcmin_seq, rcmin_ind, tnmin_seq, tnmin_ind, tvmin_seq, tvmin_ind)
        # print(s, ", ".join([str(m) for m in all_mins]))
        all_min = min(all_mins)
        
        # if len(list(set(all_mins))) < len(all_mins) and len(list(set(all_seqs))) == len(all_seqs):
            # print(f"needs tiebreaker {min_seq} vs {rcmin_seq} vs {tnmin_seq} vs {tvmin_seq}")
        
        if all_mins.index(all_min) == 0:
            proto = min_seq
            tx = ('e', perm)
        elif all_mins.index(all_min) == 1:
            proto = rcmin_seq
            tx = ('rc', rcperm)
        elif all_mins.index(all_min) == 2:
            proto = tnmin_seq
            tx = ('tn', tnperm)
        elif all_mins.index(all_min) == 3:
            proto = tvmin_seq
            tx = ('tv', tvperm)
        else:
            # tiebreaker? idk
            print(f"needs tiebreaker {min_seq} vs {rcmin_seq} vs {tnmin_seq} vs {tvmin_seq}")
            proto = min_seq
            tx = ('?',0)
        
        # if s == "ACG":
        #     print(s, "|", proto, tx)
        
        if not proto in fams:
            # print(f"adding proto {proto} to dict for seq {s}")
            fams[proto] = []
        
        fams[proto].append(s)
        txs[s] = tx
        
    # print(tabulate(min_tab, headers = ["Seq", "E", "RC","TN","TV"]))
        
    return fams, txs
        
def test_count_repeats(n = 3, suppress = False):
    fams, txs = count_nrepeats(n=n)
    
    # protos_srt = list(sorted(fams.keys(), key = lambda s:combi.get_seq_index(s)))
    protos_srt = fams
    lines = []
    
    lines.append(f"Repeat with length 3 produces {4**n} sequences with {len(fams)} non-degenerate families")
    for f in protos_srt:
        lines.append(f"prototype: {f} ({combi.get_seq_index(f)})")
        memstrs = [f"{m} ({combi.get_seq_index(m)})" for m in fams[f]]
        lines.append(f"{len(fams[f])} members: {", ".join(memstrs)}")
        lines.append("\n")
    
    if not suppress:
        for line in lines:
            print(line)
        
    return fams, txs

def test_txs(test_str, nperms = 3):
    
    headers = ["Identity","Rev Comp","Transition","Transversion"]
    test_strs = [combi.permute(test_str, i) for i in range(nperms)]
    rows = []
    for ts in test_strs:
        row = [ts, reverse_complement(ts), transition(ts), transvert(ts)]
        row_strs = [f"{r} ({combi.get_seq_index(r)})" for r in row]
        rows.append(row_strs)
    print(tabulate(rows, headers = headers))
    print()
    
def main():
    
    bg, fg = draw.get_color_scheme("test")

    for n in range(2, 7):
        fams, txs = test_count_repeats(n=n, suppress = True)

        inds = list(sorted([combi.get_seq_index(f) for f in fams]))
        print(f"{n}: {", ".join([str(ind) for ind in inds])}")
    
    # test_random_perm(9)
    
    # test_seq = "ACG"
    # print()
    # test_txs(test_seq)
    # res = find_min_perm(test_seq)
    # print(res)
    # res = find_min_perm(reverse_complement(test_seq))
    # print(res)
    # res = find_min_perm(transition(test_seq))
    # print(res)
    # res = find_min_perm(transvert(test_seq))
    # print(res)
    
    # ACG (14) reduces to ATC (7), not ATG (9)
    
    # test_txs("AAA", nperms = 1)
    # test_txs("AAT", nperms = 3)
    # test_txs("AAG", nperms = 3)
    # test_txs("AAC", nperms = 3)
    # test_txs("ATG", nperms = 3)
    # test_txs("AGT", nperms = 3)
    

if __name__=="__main__":
    main()
