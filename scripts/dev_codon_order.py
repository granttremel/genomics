
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

from ggene.genomemanager import GenomeManager

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


def get_codon_frequency_stats(gm, chrs = ["12"], max_num_genes = None):
    
    aa_codon_map = {aa:get_aa_codons(aa) for aa in ORDER_AA}
    codon_code_map = {}
    for aa, cds in aa_codon_map.items():
        codon_code_map.update({cd: aa_codon_map.get(aa).index(cd) for cd in cds})
    
    if chrs:
        pass
        # chrs = [chr]
    else:
        chrs = gm.gene_map.chromes
        
    freqs = {aa:{i:0 for i in range(len(cds))} for aa, cds in aa_codon_map.items()}
    
    txs = []
    
    ng = 0
    for chr in chrs:
        
        genes = gm.gene_map.fetch(chr, start=0, features = ("gene",))
        print(f"analyzing chrom {chr}")
        
        for g in genes:
            ginfo = g.get("info",{})
            gene_name = ginfo.get("gene_name")
            cdss, tx = get_transcript_cds(gm, gene_name, chr, gene=g)
            txinfo = tx.get("info",{})
            txid = txinfo.get("transcript_id")
            txs.append(tx)
            cdsseq = "".join(cdss)
            max_cd = 3*(len(cdsseq)//3)
            
            for i in range(0, max_cd, 3):
                cd = cdsseq[i:i+3]
                aa = CODON_TABLE[cd]
                freqs[aa][codon_code_map[cd]] += 1
            ng +=1
            
            if ng % 100 == 0:
                print(f"analyzed gene {ng} {gene_name}, assembled transcript {txid} with number cdss {len(cdss)} and {len(cdsseq)} bps")
                
            if max_num_genes and ng >= max_num_genes:
                break
        if max_num_genes and ng >= max_num_genes:
            break
    
    freqs_total = {aa:sum(fs.values()) for aa, fs in freqs.items()}
    freqs_norm = {aa:{i:f/freqs_total.get(aa) for i,f in fs.items()} for aa, fs in freqs.items()}
    
    return freqs_norm, freqs_total, txs

def display_codon_data(aas = [], fields = {}, disp_fields = [], width = 8, index_func = None, color_func = None):
    
    if not aas:
        aas = ORDER_AA
    if not index_func:
        index_func = lambda aa, cd: get_aa_codons(aa).index(cd)
    if not color_func:
        color_func = lambda aa, cd, f, fval: False
    
    num_codons_max = 6
    
    rst ="\x1b[0m"
    col_hi = "\x1b[38;5;40m"
    
    fkeys = list(fields.keys())
    if not disp_fields:
        disp_fields = fkeys
    ncols = 1 + len(disp_fields)
    
    all_cds = ["".join([ORDER[q] for q in ijk]) for ijk in itertools.product(range(4), range(4), range(4))]
    cd_aas = [CODON_TABLE.get(cd) for cd in all_cds]
    
    header = []
    rows = [[] for n in range(num_codons_max)]
    
    for aa in ORDER_AA:
        if not aa in aas:
            continue
        
        header.extend([aa] + disp_fields)
        aacds = get_aa_codons(aa)
        for i in range(num_codons_max):
            if i>=len(aacds):
                rows[i].extend([" "]*ncols)
            else:
                cd = aacds[i]
                cdind = index_func(aa, cd)
                
                cdstr = cd
                if color_func(aa, cd, {}, None):
                    cdstr = f"{col_hi}{cd}{rst}"
                row = [cdstr]
                
                cdfields = {}
                for fk in fkeys:
                    fval = fields.get(fk)(aa, cd, cdfields)
                    fstr = str(fval)
                    if color_func(aa, cd, cdfields, fval):
                        fstr = f"{col_hi}{fval}{rst}"
                    cdfields[fk] = fval
                    if fk in disp_fields:
                        row.append(fstr)
                
                rows[cdind].extend(row)
        
    max_row = 256
    max_row_items = max_row // width
    max_cols = ncols * ((max_row_items) // ncols)
    num_groups = len(header)//max_cols + 1
    
    for ncol in range(num_groups):
        subtab = []
        
        subhdr = ["Ind"] + header[ncol*max_cols:ncol*max_cols+max_cols]
    
        for i in range(num_codons_max):
            r = rows[i]
            subrow = [str(i)] + r[ncol*max_cols:(ncol+1)*max_cols]
            subtab.append(subrow)
        
        print(tabulate(subtab, headers = subhdr, tablefmt = "simple"))
        print()
    
    return all_cds, cd_aas

def expand_codons():
    
    all_cds = ["".join([ORDER[q] for q in ijk]) for ijk in itertools.product(range(4), range(4), range(4))]
    cd_aas = [CODON_TABLE.get(cd) for cd in all_cds]
    
    full_codons = []
    
    for aa in ORDER_AA:
        print(aa)
        for cd in get_aa_codons(aa):
            cdfull = cd + rc(cd)
            rccdfull = rc(cdfull)
            cdfull2 = rc(cd) + cd
            rccdfull2 = rc(cdfull2)
            if cdfull == cdfull2:
                print(cdfull, cdfull2,"symm")
            if cdfull == rc(cdfull2):
                print(cdfull, cdfull2,"antisymm")
            else:
                print(cdfull, cdfull2)
            
            # ord_full = ord2(aa, cdfull, {})
            
            # ofm6 = ord_full%6
            # offp6 = (ord_full//6)%6
            # print(" ", cdfull, ofm6, offp6, (ofm6-offp6)%6)
            # print(" ", cdfull, (ofm6+offp6)%6, (ofm6-offp6)%6)
            
            full_codons.append(cdfull)
            full_codons.append(cdfull2)
            full_codons.append(rccdfull)
            full_codons.append(rccdfull2)
            
        print()
    max_full = 4**6 # 64^2
    fc_dd = list(set(full_codons))
    print(f"{len(full_codons)} full codons out of {max_full}, {len(fc_dd)} after deduplication")
    
    pass

def ord(a, c,f):
    return sum([ORDER.index(b)*2**(len(c)-i-1) for i, b in enumerate(c)])

def ord2(a, c,f):
    return sum([ORDER.index(b)*2**i for i, b in enumerate(c)])

def ord3(a, c, f):
    b1, b2 ,b3 = list(c)
    b12 = b1+b2
    inds= []
    if b12 in ["AC", "TC","GT","GG","GC","CT","CG","CC"]: # totally defined by final base
        inds.append(b3 in ALIASES.get("R")) # final base is pur or pyr?
        inds.append(b3 in ALIASES.get("W")) # final base is AT or GC?
    elif b12 in ["AA", "AG","TA","TG","GA","CA"]: # partially defined by final base
        inds.append(b3)    
    elif b12 in ["AT","TG"]:
        inds.append(b3)
    
    
def main():
    
    bg, fg = draw.get_color_scheme("test")
    

    gm = load_genome()
    genemap = gm.gene_map

    prote_vc = 'ABCDEF'
    
    # gene_name = "ERC1"
    # chr = 12
    gene_name = "AGRN"
    chr = 1
    
    if False:
        ran = 50000
        gene_feats = [f for f in genemap.fetch(1, 1055000-ran, 1055100+ran, features =("gene",))]
        gene = [g for g in gene_feats if g.get("info").get("gene_name")==gene_name][0]
        seq, txs = get_transcript_cds(gm, gene_name, chr, gene = gene)
        seq = "".join(seq)
        
        aacodes, lines = display_alt_codons(seq, use_color = True, suppress = False)
        _, _ = display_alt_codons(seq, suppress = True, fname = "codon_choices_AGRN.txt")
        
        aacodedict = {}
        aa_reord = {aa: [] for aa in ORDER_AA}
        
        for c in aacodes:
            aa = c[0]
            cd = int(c[1])
            codon = get_aa_codons(aa)[cd]
            if not aa in aacodedict:
                aacodedict[aa] = dict()
            if not cd in aacodedict[aa]:
                aacodedict[aa][cd] = 0
            aa_reord[aa].append(codon)
            aacodedict[aa][cd] += 1
        
        display_alt_codon_entropy(aacodedict, suppress = True)
    
    # freqs_norm, freqs_ct, txs = get_codon_frequency_stats(gm, chrs = ["12"])
    # display_alt_codon_entropy(freqs_norm, freqs_ct, fname = "codon_entropy_chr12.txt")
    
    # codons = "".join("".join(get_aa_codons(aa)) for aa in ORDER_AA)
    # display_alt_codons(codons, use_color = True)
    # cdgcs = {aa:[cd.count("G") + cd.count("C") for cd in get_aa_codons(aa)] for aa in ORDER_AA}
    # for aa, cdgcs in cdgcs.items():
    #     print(aa, cdgcs)
    
    # orddict = {ORDER[i]:2**(i-4) for i in range(4)}
    

        
        
        # ind = ORDER.index(b3)
        
        
        
    
    fields = {
        "GC":lambda a,c,f: c.count("G") + c.count("C"),
        
        "Ord": ord,
        "Ord//6": lambda a,c,f:f.get("Ord")//6,
        "Ord%6": lambda a,c,f:f.get("Ord")%6,
        
        "Ord2": ord2,
        "Ord2%6": lambda a,c,f:f.get("Ord2")%6,
        
        # "Ord3":ord3,
        # "Ord4":lambda a,c,f: ord3(a,c,f) + ord3(a, reverse_complement(c), f),
        
        "RC": lambda a, c,f: reverse_complement(c),
        "RC_AA": lambda a,c,f: CODON_TABLE.get(f.get("RC")),
        "RCOrd": lambda a,c,f: ord(a, f.get("RC"), f),
        "RCOrd%6": lambda a,c,f:f.get("RCOrd")%6,
        
        "Rev_AA": lambda a,c,f: CODON_TABLE.get(reverse(c)),
        "P1_AA": lambda a,c,f: CODON_TABLE.get(c[1:] + c[:1]),
        "P2_AA": lambda a,c,f: CODON_TABLE.get(c[2:] + c[:2]),
        
        "SumOrds": lambda a,c,f: f.get("Ord") + f.get("RCOrd"),
        "SumOrds%6": lambda a,c,f:f.get("SumOrds")%6,
        "DiffOrds": lambda a,c,f: f.get("Ord") - f.get("RCOrd"),
        "DiffOrds%6": lambda a,c,f:f.get("DiffOrds")%6,
        
    }
    
    disp_fields = [
        # "GC",
        # "Ord",
        # "Ord//6",
        "Ord%6",
        # "Ord2%6",
        # "Ord2",
        # "RC",
        # "RC_AA",
        # "RCOrd",
        "RCOrd%6",
        # "Rev_AA",
        # "P1_AA",
        # "P2_AA",
        # "SumOrds",
        "SumOrds%6",
        # "DiffOrds",
        "DiffOrds%6",
        ]
    
    color_func= lambda aa, cd, f, fv: fv in get_aa_codons(aa)
    
    # display_codon_data(
    #     aas = "LSR",
    #     fields = fields, 
    #     disp_fields = disp_fields, 
    #     color_func = color_func
    #     )
    
    expand_codons()
    
    """

    
    
    
    """
    
    # get_cds_stats(gm)
    """
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





"""