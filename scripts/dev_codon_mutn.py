
from math import e
from subprocess import list2cmdline
import numpy as np
import random
from matplotlib import pyplot as plt
import os
from pathlib import Path

from tabulate import tabulate

import itertools

from ggene import seqs, draw
from ggene.seqs import bio
from ggene.seqs.bio import reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE, ALIASES, ALIASES_INV, ALIASES_REV

rc = reverse_complement
mutation_specs = list(ALIASES.keys())

from ggene.genomemanager import GenomeManager

# class SymmetricDict(dict):
    
#     def __setitem__(self, v1, v2):
#         if v1 in self:
#             v1s = self.__getitem__(v1)
#         else:
#             v1s = ''
#         if v2 in self:
#             v2s = self.__getitem__(v2)
#         else:
#             v2s = ''
        
#         if not v1 in v2s:
#             v2s = v2s + v1
#         if not v2 in v1s:
#             v1s = v1s + v2
        
#         dict.__setitem__(self, v1, v1s)
#         dict.__setitem__(self, v2, v2s)
    
#     def __missing__(self, k):
#         return k

# SymDict = SymmetricDict

class MutationMap(dict):

    _cache = {}
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._spec = ""
    
    def mutate(self, base, choice = None):
        out_bases = self[base]
        if choice and len(out_bases) > 1:
            return out_bases[choice]
        elif not choice:
            return random.choice(out_bases)
        else:
            return out_bases
    
    def mutate_each(self, base):
        out_bases = self[base]
        return list(out_bases)
        
    def mutate_position(self, seq, pos):
        seqlist = list(seq)
        seqlist[pos] = self.mutate(seq[pos])
        return "".join(seqlist)
    
    def mutate_positions(self, seq, positions):
        seqlist = list(seq)
        for pos in positions:
            seqlist[pos] = self.mutate(seq[pos])
        return "".join(seqlist)
    
    def get_all_mutation_products(self, seq):
        
        all_seqs = []
        
        for p in range(len(seq)):
            
            new_b = self.mutate_each(seq[p])
            
            seqlist = list(seq)
            for b in new_b:
                seqlist[p] = b
                all_seqs.append("".join(seqlist))
        
        return list(set(all_seqs))
    
    def get_all_mutations(self, seq):
        
        all_seqs = {}
        
        for p in range(len(seq)):
            seqlist = list(seq)
            new_b = self.mutate_each(seq[p])
            # print(f"seq[{p}] = {seq[p]} -> {new_b}")
            for b in new_b:
                if b == seq[p]:
                    continue
                seqlist[p] = b
                mut = bio.get_alias(seq[p], b)
                all_seqs[(mut, p)] = "".join(seqlist)
        
        return all_seqs
    
    def __setitem__(self, v1, v2):
        if v1 in self:
            v1s = self.__getitem__(v1)
        else:
            v1s = ''
        if v2 in self:
            v2s = self.__getitem__(v2)
        else:
            v2s = ''
        
        if not v1 in v2s:
            v2s = bio.get_alias(v2s, v1, exclude = v2)
        if not v2 in v1s:
            v1s = bio.get_alias(v1s, v2, exclude = v1)
        
        dict.__setitem__(self, v1, v1s)
        dict.__setitem__(self, v2, v2s)
    
    def __getitem__(self, v1):
        vout = dict.__getitem__(self, v1)
        return ALIASES.get(vout, vout)
    
    def __missing__(self, k):
        return k
    
    def __repr__(self):
        
        parts = []
        for v1, v1s in self.items():
            parts.append(f"{v1}→{v1s}")
        
        return f"MutationMap({self._spec}: {", ".join(parts)})"
    
    @classmethod
    def from_spec(cls, mutation_spec):
        
        if mutation_spec in cls._cache:
            return cls._cache.get(mutation_spec)
        
        sd = cls()
        sd._spec = mutation_spec
        
        for s in mutation_spec:
            bout = ALIASES.get(s, s)
            
            for bi in bout:
                for bj in bout:
                    if bi!=bj:
                        sd[bi] = bj
        
        cls._cache[mutation_spec] = sd
        return sd

MM = MutationMap
    
def load_genome():
    return GenomeManager()

def get_transcript_cds(gm, gene_name, chr, gene = None, txid=""):
    
    if not gene:
        gene = gm.assemble_gene(gene_name, chr)
        start = gene.start
        end = gene.end
    else:
        start=  int(gene.get("start",-1))
        end = int(gene.get("end", -1))
        
    features = ["gene","transcript","CDS"]
    feat_dicts = [d for d in gm.gene_map.fetch(chr,start,end, features =features)]
    cdss = [f for f in feat_dicts if f.get("feature") == "CDS"]
    if not txid:
        txid = cdss[0].get("info",{}).get("transcript_id")
    
    tx_feats = [f for f in feat_dicts if f.get("feature") == "transcript" and f.get("info",{}).get("transcript_id") == txid]
    if len(tx_feats) > 0:
        tx_feat = tx_feats[0]
    else:
        tx_feat = {}
    
    tx_cdss = [cds for cds in cdss if cds.get("info",{}).get("transcript_id") == txid]
    
    full_seq = []
    
    for cds in tx_cdss:
        info = cds.get("info",{})
        frame = int(info.get("frame", 0))
        cds_seq = gm.get_sequence(chr, cds.get("start"), cds.get("end"), frame=frame)
        full_seq.append(cds_seq)
        
    if tx_feat.get("strand", "+") == "-":
        full_seq = [reverse_complement(f) for f in reverse(full_seq)]
    
    return full_seq, tx_feat

def get_alt_codon_highlight(codon, marker = "*", high = True, width = 5):
    
    frmstr = "{pre}{codon}{post}{delimit}"
    
    if isinstance(marker, int):
        pre_norm = "\x1b[38;5;240m"
        pre_hi = f"\x1b[38;5;{marker}m"
        post = "\x1b[0m"
        delimit = " "*(width - len(codon))
    elif isinstance(marker, str):
        pre_norm =" "
        pre_hi = marker
        post = ""
        delimit = " "*(width - len(codon) - 1)
    
    if high:
        pre = pre_hi
    else:
        pre = pre_norm
        
    return frmstr.format(pre=pre, codon = codon, post = post, delimit = delimit)
    
def display_alt_codons(cds_seq, use_color = False, suppress = False, max_print = 5, fname = ""):
    
    # max_print = 10
    codons_per_aa = 6
    codon_width = 3
    col_width = 4
    
    if use_color:
        marker_hi = 40
    else:
        marker_hi = "*"
    
    chunksz = codon_width*(212//col_width)
    nchunks = len(cds_seq)//chunksz + 1
    
    lognstates = 0
    aacodes = []
    
    outlines = []
    
    for n in range(nchunks):
        cds_subseq = cds_seq[n*chunksz:(n+1)*chunksz]
        aaline = []
        lines = [[] for i in range(codons_per_aa)]
        for i in range(0, len(cds_subseq), 3):
            
            codon = cds_subseq[i:i+3]
            aa = CODON_TABLE.get(codon,".")
            codons = get_aa_codons(aa)
            lognstates += np.log2(len(codons)) if codons else 1
            if len(codons) < codons_per_aa:
                codons = codons + ["   "]*(codons_per_aa-len(codons))
            
            aaind = codons.index(codon) if not aa == "." else " "
            aacode = aa+str(aaind)
            aacodes.append(aacode)
            aahead = get_alt_codon_highlight(aacode, high=False, width=5)
            aaline.append(aahead)
            
            for cd, line in zip(codons, lines):
                hi = cd==codon
                cd_hi = get_alt_codon_highlight(cd, marker = marker_hi, high = hi, width = 5)
                line.append(cd_hi)
        
        if n < max_print:
            outlines.append("".join(aaline))
            outlines.extend(["".join(l) for l in lines])
            outlines.append("\n")
        
    if nchunks - max_print > 0:
        outlines.append(f" ... and {nchunks-max_print} more rows..")
    outlines.append(f"this protein sequence has 2^{lognstates:0.0f} possible encoding DNA sequences")
    
    if not suppress:
        for line in outlines:
            print(line)
        print()
    
    if fname:
        fdir = Path("./data/outputs")
        fpath = fdir / fname
        
        with open(fpath, "w+") as f:
            for line in outlines:
                f.write(draw.scrub_ansi(line) + "\n")
        
    return aacodes, outlines

def display_alt_codon_entropy(aacodedict, sampledict = {}, labels = False, suppress = False, fname = ""):
    
    bg, fg = draw.get_color_scheme("test")
    
    lines =[]
    
    aacode_norm = {aa:{cd:f / sum(cdmap.values()) for cd, f in cdmap.items()} for aa, cdmap in aacodedict.items()}
    aacode_entropy = {aa:sum([-f*np.log(f) for f in cdmap.values()]) for aa, cdmap in aacode_norm.items()}
    
    lines.append("{:<10}{:<10}{:<12}{:<10}{:<10}{:<10}".format(*["AA", "Entropy", "Max Entropy", "ncodons","nsamps", "dist"]))
    
    for aa, ent in aacode_entropy.items():
        ncodons = len(aacodedict[aa])
        nsamps = sampledict.get(aa, 1)
        max_entropy = np.log(ncodons)
        dist = draw.scalar_plot_distribution(aacodedict[aa], bg_color = bg, fg_color = fg-12, labels = labels)
        
        lines.append("{:<10}{:<10.3f}{:<12.3f}{:<10}{:<10}{:<10}".format(*[aa, ent, max_entropy, ncodons, nsamps, dist[0]]))
        if labels:
            lines.append("{:<52}{:<10}".format(*["", "".join(dist[1])]))
    
    if not suppress:
        for line in lines:
            print(line)
    
    if fname:
        fdir = Path("./data/outputs")
        fpath = fdir / fname
        
        with open(fpath, "w+") as f:
            for line in lines:
                f.write(draw.scrub_ansi(line) + "\n")
    
    return aacode_entropy

def make_mutation_adj_map(muts = ['N']):
    
    mut_to_map = {}
    cod_to_aa = {}
    aa_to_adj = {}
    
    for mut in muts:
        for aa in ORDER_AA:
            
            cod, adj = get_adjacent_codons(aa, mutations = mut)
            
            cod_to_aa.update(cod)
            aa_to_adj[aa] = list(adj.keys())
        mut_to_map[mut] = (cod_to_aa, aa_to_adj)
        
    return mut_to_map

def disp_seq(seq, max_disp = 256):
    seq_len = len(seq)
    ndisp = seq_len // max_disp
    for n in range(ndisp):
        print(n, seq[n*max_disp:(n+1)*max_disp])
    print(n+1, seq[(n+1)*max_disp:])
    print()

def test_mute_map():
    
    tests = 10
    nspecs_max = 5
    nbases = 15
    npos = 10
    
    for n in range(tests):
        
        nspecs = random.randint(1, nspecs_max)
        mspecs = "".join(random.sample(mutation_specs, k = nspecs))
        mm = MM.from_spec(mspecs)
        print(repr(mm))
        
        seq = "".join(random.choices(ORDER[:4], k = nbases))
        
        mpos = [random.randint(0, nbases-1) for p in range(npos)]
        newseq = mm.mutate_positions(seq, mpos)
        print(f"{mspecs}: {seq} → {newseq}")
        print()
        
    pass

def test_codon_stability(mutation = "N", aas = [], suppress = False, index_func = None):
    
    mm = MM.from_spec(mutation)
    # print(repr(mm))
    
    fg = 53
    
    if not index_func:
        index_func = lambda aa, cd: get_aa_codons(aa).index(cd)
    
    cds = list(CODON_TABLE.keys())
    if not aas:
        aas = list(ORDER_AA)
    
    header = ["Codon", "Same", "Other", "Unique","Total"]
    tabs = {}
    all_dists = {}
    
    for aa in aas:
        tab = []
        dist = {}
        dist_unique = {}
        cds = get_aa_codons(aa)
        cds = sorted(cds, key = lambda cd:index_func(aa, cd))
        for cd in cds:
            num_same = num_total = 0
            unique = set()
            cd_muted = mm.get_all_mutations(cd)
            for (m, i), new_cd in cd_muted.items():
                new_aa = CODON_TABLE.get(new_cd)
                num_total += 1
                if aa == new_aa:
                    num_same += 1
                unique.add(new_aa)
            num_unique = len(unique)
            row = [cd, num_same, num_total - num_same, num_unique, num_total]
            tab.append(row)
            dist[cd] = num_same
            dist_unique[cd] = num_unique
        
        bars = draw.scalar_plot_distribution(dist, key_order = cds, edges = True, maxval = num_total,fg_color = fg)[0]
        bars_unique = draw.scalar_plot_distribution(dist_unique, key_order = cds, edges = True, maxval = num_total,fg_color = fg)[0]
        if not suppress:
            print(f"AA = {aa} {bars} {bars_unique}")
            print(tabulate(tab, headers = header, tablefmt = "simple"))
            print()
        tabs[aa] = tab
        all_dists[aa] = (bars, bars_unique)
    return tabs, all_dists

def test_all_codon_stability(fname = "", aas = [], index_func = None):
    test_mutations = ["N","YR","SW","KM","B","D","H","V"]
    
    if not aas:
        aas = ORDER_AA
    
    headers = ["AA", "Data"] + [f"{mut} Dist." for mut in test_mutations]
    rows_same = []
    rows_uni = []
    for i in range(len(aas)):
        rows_same.append([aas[i], "Same"])
        rows_uni.append([aas[i], "Unique"])
    
    for tm in test_mutations:
        
        tabs, bars = test_codon_stability(tm, aas = aas, suppress = True, index_func = index_func)
        
        for i,aa in enumerate(bars):
            rows_same[i].append(bars[aa][0])
            rows_uni[i].append(bars[aa][1])
    
    rows_all = rows_same + [" "] + rows_uni
    tab = tabulate(rows_all, headers = headers, tablefmt = "rounded_grid", stralign = "center")
    # tab_uni = tabulate(rows_uni, tablefmt = "rounded_grid", stralign = "center")
    if fname:
        fdir = Path("./data/outputs")
        fpath= fdir / fname
        with open(fpath,"w+") as f:
            for line in tab:
                f.write(line)
    
    print(tab)
    # print(tab_uni)
    print()

def permute(seq, n):
    return seq[n:] + seq[:n]

def test_codon_transformations():
    
    txs = [
        rc,
        reverse,
        complement,
        lambda cdn: permute(cdn, 1),
        lambda cdn: permute(cdn, 2),
        lambda cdn: rc(permute(cdn, 1)),
        lambda cdn: rc(permute(cdn, 2)),
        lambda cdn: reverse(permute(cdn, 1)),
        lambda cdn: reverse(permute(cdn, 2)),
    ]
    
    hdr = ["Codon", "RC", "Rev","Comp","P1","P2","RCP1","RCP2*","RP1", "RP2*"]
    for aa in ORDER_AA:
        tab = []
        cdns = get_aa_codons(aa)
        print(f"AA: {aa}")
        for cdn in cdns:
            row = [cdn]
            for tx in txs:
                tx_cdn = tx(cdn)
                txaa = CODON_TABLE[tx_cdn]
                row.append(txaa)
            tab.append(row)
        print(tabulate(tab, headers = hdr, tablefmt = "simple"))
        print()
    
def test_rcp2():
    
    txs = [
        lambda cdn: rc(permute(cdn, 2)),
    ]
    
    hdr = ["AA", "Codon", "RCP2", "New AA"]
    for aa in ORDER_AA:
        tab = []
        cdns = get_aa_codons(aa)
        for cdn in cdns:
            row = [aa, cdn]
            for tx in txs:
                tx_cdn = tx(cdn)
                txaa = CODON_TABLE[tx_cdn]
                row.append(tx_cdn)
                row.append(txaa)
            tab.append(row)
        print(tabulate(tab, headers = hdr, tablefmt = "simple"))
        print()

def random_walk_mutations(mutation_spec = "N", start = "ATG", posses = [0,1,2], nsteps = 10):
    
    mm = MM.from_spec(mutation_spec)
    start_aa = CODON_TABLE.get(start)
    path = [(start, start_aa)]
    hist_cdn = {cdn:0 for cdn in CODON_TABLE}
    hist_aa = {aa:0 for aa in ORDER_AA}
    seq = start
    for n in range(nsteps):
        pos = random.choice(posses)
        res = mm.mutate_position(seq, pos)
        aa_res = CODON_TABLE.get(res)
        path.append((res, aa_res))
        
        hist_cdn[res] += 1
        hist_aa[aa_res] += 1
        seq = res
    return path, hist_cdn, hist_aa

def plot_transition_hist(all_stops):
    
    aa_hist = {(aa, ab):0 for aa, ab in itertools.product(ORDER_AA, ORDER_AA)}
    
    lastaa = ""
    for cdn, aa in all_stops:
        
        if not lastaa:
            lastaa = aa
            continue
        
        aa_hist[(lastaa, aa)] += 1
        
        lastaa = aa
    
    unused = []
    for k, v in aa_hist.items():
        if v == 0:
            unused.append(k)
            
    for k in unused:
        del aa_hist[k]
    
    bars = draw.scalar_plot_distribution(aa_hist, bit_depth=16)
    # for r in bars:
    #     print(r)
    # print()
    
    return aa_hist, unused
    
def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    # question: how unstable is any codon related to any other? are there sub-circuits? does it depend on allowed mutation mode? i.e. only S<->W, or Y<->R
    
    # def ord(aa, cdn, offset, sign, vocab = "ATGC"):
    #     return sum([vocab.index(s)*2**((offset + sign*i) % 3) for i, s in enumerate(cdn)])
        
    # test_all_codon_stability(index_func = ord)
    
    # ind_f = lambda aa, cdn: ord(aa, cdn, 0, 1)
    
    # ind_f = lambda aa, cdn: ord(aa, cdn, 0, 1, vocab = "AGCT")
    # test_all_codon_stability(aas = ["I","L","R","S","T"],index_func = ind_f)
    
    # start = "ATC"
    start = "ATG"
    
    path, hist_cdn, hist_aa = random_walk_mutations(mutation_spec="N", start = start, nsteps = 200)
    
    all_stops = list(set(path))
    print(f"{len(all_stops)} codons hit")
    
    bars_cdn = draw.scalar_plot_distribution(hist_cdn, bit_depth = 16)
    print("Codons:")    
    for r in bars_cdn:
        print(r)
    print()
    bars_aa = draw.scalar_plot_distribution(hist_aa, bit_depth = 16)
    print("Amino acids:")
    for r in bars_aa:
        print(r)
    print()
    
    min_cdn = [k for k,v in hist_cdn.items() if v == min(hist_cdn.values())][0]
    max_cdn = [k for k,v in hist_cdn.items() if v == max(hist_cdn.values())][0]
    min_aa = [k for k,v in hist_aa.items() if v == min(hist_aa.values())][0]
    max_aa = [k for k,v in hist_aa.items() if v == max(hist_aa.values())][0]
    
    print(f"AA min {min_aa} with {hist_aa[min_aa]}, max {max_aa} with {hist_aa[max_aa]}")
    print(f"Codon min {min_cdn} with {hist_cdn[min_cdn]}, max {max_cdn} with {hist_cdn[max_cdn]}")
    
    tx_hist, unused = plot_transition_hist(path)
    
    min_aatx = [k for k,v in tx_hist.items() if v == min(tx_hist.values())][0]
    max_aatx = [k for k,v in tx_hist.items() if v == max(tx_hist.values())][0]
    
    print(f"{len(tx_hist)} amino acid transitions did occur, {len(unused)} never occurred")
    print(f"most traveled amino acid transition: {max_aatx[0]}=>{max_aatx[1]}, with {tx_hist[max_aatx]} hits")
    
    key_order = ORDER_AA
    thr = {}
    for k0 in key_order:
        if not k0 in thr:
            thr[k0] = {}
        for k1 in key_order:
            d0 = thr[k0]
            if not k1 in d0:
                d0[k1] = tx_hist.get((k0, k1), 0)
        
    plot2d = draw.scalar_plot_2d_distribution(thr, key_order = key_order, min_color = 16, max_color = 117)
    for line in plot2d:
        print(line)
    
    # test_rcp2()
    """
    generates some pairs:
    A ⟵⟶ A
    C ⟵⟶ Q
    D ⟵⟶ S (TCR)
    E ⟵⟶ S (GAR)
    F ⟵⟶ K 
    G ⟵⟶ P 
    I ⟵⟶ I (ATW)
    N ⟵⟶ L (TTR)
    L ⟵⟶ S
    L ⟵⟶ R
    R ⟵⟶ R
    T ⟵⟶ V
    Y ⟵⟶ * (TAR)
    
    outliers:
    H CAT ⟵⟶ * TGA
    H CAC ⟵⟶ W TGG
    M ATG ⟵⟶ I ATC
    
    maybe ind(cdn) + ind(rc(perm(cdn,2))) = something
    
    """
    
    
    # mmyr = MM.from_spec("N")
    # print(repr(mmyr))
    
    # cds = list(CODON_TABLE.keys())
    
    # for cd in cds:
    #     aa = CODON_TABLE.get(cd)
    #     cd_muted = mmyr.get_all_mutations(cd)
    #     for (m, i), new_cd in cd_muted.items():
    #         new_aa = CODON_TABLE.get(new_cd)
    #         print(f"{m}, {i}: {aa} -> {new_aa}")
    #     print()
    
    # seq = "ACGCTACGCTA"
    
    # all_prods = mmyr.get_all_mutation_products(seq)
    # print(all_prods)
    
    # all_mutes = mmyr.get_all_mutations(seq)
    # for (m, i), new_seq in all_mutes.items():
    #     print(f"{m}, {i}: {seq} -> {new_seq}")
    # print()
    
    # test_mute_map()
    
    # cdns = ["AAA","ATG","GTC","CGA"]
    # muts = ["R","S","K","YR"]
    # pos = [1, 2, 0, 1]
    
    # for cdn, p, m in zip(cdns, pos, muts):
    #     # newcdn = mutate_codon(cdn, p, m)
    #     mm = MM.from_spec(m)
    #     newcdn = mm.mutate_position(cdn, p)
    #     print(f"{cdn} --{m}{p}-> {newcdn}")
    








if __name__=="__main__":
    main()



















"""
codon indexer desired properties:
num_codons = 4^3 = 64
num_aas = 20
num_slots = 6

0 <= ind(cd) < num_slots    

ind(cd) = 6 - ind(rc(cd))

no redundancies, i.e.
ind(cd_aa_i) != ind(cd_aa_j)


CODON_TABLE = {
    TTT: F, TTC: F, TTA: L, TTG: L,
    TCT: S, TCC: S, TCA: S, TCG: S,
    TAT: Y, TAC: Y, TAA: *, TAG: *,
    TGT: C, TGC: C, TGA: *, TGG: W,
    
    CTT: L, CTC: L, CTA: L, CTG: L,
    CCT: P, CCC: P, CCA: P, CCG: P,
    CAT: H, CAC: H, CAA: Q, CAG: Q,
    CGT: R, CGC: R, CGA: R, CGG: R,
    
    ATT: I, ATC: I, ATA: I, ATG: M,
    ACT: T, ACC: T, ACA: T, ACG: T,
    AAT: N, AAC: N, AAA: K, AAG: K,
    AGT: S, AGC: S, AGA: R, AGG: R,
    
    GTT: V, GTC: V, GTA: V, GTG: V,
    GCT: A, GCC: A, GCA: A, GCG: A,
    GAT: D, GAC: D, GAA: E, GAG: E,
    GGT: G, GGC: G, GGA: G, GGG: G
}

Order: 
A<T<G<C
R<Y

?
A < G
^   ^
T < C

6f:
L   CTN|TTR
S   TCN|AGY
R   CGN|AGR

4f:
T   ACN
V   GTN
G   GGN
A   GCN
P   CCN

3f:
I   AT[ATC] = ATH
*   TAR|TGA

2f:
K   AAR
N   AAY
Y   TAY
F   TTY
C   TGY
E   GAR
D   GAY
Q   CAR
H   CAY

1f:
M   ATG
W   TGG

zero bit first pairs:
AC, TC, GT, GG, GC, CT, CG, CC
(CT CG CC) | (GT GG GC) | (AC TC)
[GC][TGC] | [AT]C

one bit first pairs:
AA, AG, TA, TT, TG, GA, CA
A[A|G] | T[A|T] | [G|C]A

two bit first pairs:
AT, TG

sixfold codons of form [F|B](C[FBCD]|D[FBCD])
one bit for F|B
one bit for C..|D..
two bits for [FBCD]

L:
[C|T]T[ATGC]
- CT[ATGC]  CTN     N:2b
- TT[AG]    TTR     R:1b
=> (CTN|TTR)

TTA < TTG < CTA < CTT < CTG < CTC

TTA < CTA < CTT < CTG < CTC < TTG


S:
[A|T](G[ATGC]|C[ATGC])
- TC[ATGC]  TCN     N:2b
- AG[TC]    AGY     Y:1b
=> (TCN|AGY)

R:
[C|A]G[ATGC]
- CG[ATGC]  CGN     N:2b
- AG[AG]    AGR     R:1b
=> (CGN|AGR)

pairish:

S, R: ((TC|CG)N)|(AG[R|Y])

between row:
P, G: [GG|CC]N


T, C*W ew

within row:
D, E: GA[Y|R]
H, Q: CA[Y|R]


# get_cds_stats(gm)
all chromes + filter
in units of bp:
Mean: 161.669
SD: 276.470
Min: 1
Max: 2e+04

in units of amino acids:
Mean: 53.890
SD: 92.157
Min: 0.3
Max: 7e+03


"""