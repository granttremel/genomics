
from math import e
import numpy as np
import random
from matplotlib import pyplot as plt
import os
from pathlib import Path
import pandas as pd

import itertools
from tabulate import tabulate

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE
from ggene.seqs.heal import Healer


import reprlib

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def get_cds(gm, gene_name, chr, cds_id):
    
    gene = gm.assemble_gene(gene_name, chr)
    cdss = gene.get_feature("CDS", pf = lambda i,f : f.sfid==cds_id)
    print(f"{len(cdss)} cds identified")
    # print("".join(c.sfid for c in cdss))
    
    cds = cdss[0]
    return gm.get_sequence(chr, cds.start, cds.end)

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
        rst = "\x1b[0m"
        colstr = "\x1b[38;5;240m"
        col_hi = f"\x1b[38;5;{marker}m"
        pre_norm = colstr
        pre_hi = col_hi
        post = rst
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
    

def display_alt_codons(cds_seq, use_color = False, suppress = False, cdmap = {}, fname = ""):
    
    max_print = 10
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

def random_walk_mutation_map(c2a, a2adj, start, nsteps):
    
    path = []
    cod_path = [start]
    curr = c2a.get(start)
    curr_cod = start
    for n in range(nsteps):
        
        adj = a2adj.get(curr)
        choice = random.choice(adj)
        next = c2a.get(choice)
        
        print(f"mutated {curr} -> {next}, {curr_cod} -> {choice}")
        
        path.append(next)
        cod_path.append(choice)
        curr = next
        curr_cod = choice
        
    return path, cod_path

def disp_seq(seq, max_disp = 256):
    seq_len = len(seq)
    ndisp = seq_len // max_disp
    for n in range(ndisp):
        print(n, seq[n*max_disp:(n+1)*max_disp])
    print(n+1, seq[(n+1)*max_disp:])
    print()

def write_seq(seq, fname, max_line = 256):
    import os
    fpath = os.path.join("./data/outputs", fname)
    seq_len = len(seq)
    ndisp = seq_len // max_line
    with open(fpath, "w") as f:
        for n in range(ndisp):
            f.write(" ".join((str(n), seq[n*max_line:(n+1)*max_line])) + "\n")
        f.write(" ".join((str(n+1), seq[(n+1)*max_line:])))

def check_written_seq(fname):
    
    import os
    fpath = os.path.join("./data/outputs", fname)
    
    with open(fpath, "r") as f:
        for line in f:
            print(line)
            print(repr(line))
    
    

def color_protein_metasequence(mseq, vocab, colors):
    if not len(colors) == len(vocab):
        return
    coldict = {v:c for v,c in zip(vocab, colors)}
    print(coldict)
    
    reset = '\x1b[0m'
    bold = '\x1b[1m'
    
    max_disp = 128
    # seq_len = len(mseq)
    # ndisp = seq_len // max_disp
    
    colseq = []
    initcol = coldict.get(vocab[0])
    lastcol = initcol
    n = 0
    for s in mseq:
        
        scol = coldict.get(s)
        if not scol == lastcol:
            new_col = bold + scol
            colseq.append(new_col)
            lastcol = scol
            
        colseq.append(s)
        n +=1
        if n%max_disp == 0:
            colseq.append(reset)
            print("".join(colseq))
            colseq = [bold+lastcol]        
    colseq.append(reset)
    print("".join(colseq))
    
    return colseq

def display_all_codons():
    
    aas = bio.ORDER_AA
    codons = {aa:get_aa_codons(aa) for aa in aas}
    
    header = [f" {aa} " for aa in aas]
    lines = [["   " * len(aas)] for a in range(6)]
    
    for aa in aas:
        [f" {aa} " for aa in aas]
        for k, v in codons.items():
            
            vind = bio.get_codon_index(v)
            
            pass
    
    
    pass

def test_color_scales(num_tests):
    reset = "\x1b[0m"
    print(reset)
    
    max_num = 0
    
    col_format = "\x1b[38;5;{c}m"
    # startcols = [231,230,195,225,229,194,159,224,189,219]
    # endcols = [232 - sc + 16 for sc in startcols]
    startcols = [[random.randint(0, 40-1) for i in range(3)] for n in range(num_tests)]
    endcols = [[random.randint(0, 40-1) for i in range(3)] for n in range(num_tests)]
    
    for n in range(num_tests):
        startcol = random.choice(startcols)
        endcol = random.choice(endcols)
        num_steps = random.randint(10, 15)
        
        # cscale = draw.Colors().get_color_scale(startcol, endcol, num_steps)
        cscale = draw.Colors().get_color_scale_16b(startcol, endcol, num_steps)
        
        for i in range(len(cscale)):
            c1 = cscale[i]
            c1code = col_format.format(c = c1)
            print(f"{c1code}this is color {i} with int {c1}!{reset}")
            
        # print(f"desired number: {num_steps}, actual number: {len(cscale)}")
        # print()
        
        if len(cscale) > max_num:
            max_num = len(cscale)
    print(f"max number of colors observed: {max_num}")
    
def render_stats_table(data, label1, data2=None, label2=""):
    outlines = []
    marg = ""*10
    outlines.append("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format(*[marg,"Mean", "SD", "CV", "Min", "Max"]))
    outlines.append(f"{label1:<10}{np.mean(data):<10.2g}{np.std(data):<10.2g}{np.std(data)/np.mean(data):<10.2g}{min(data):<10.2g}{max(data):<10.2g}")
    if data2 is not None:
        marg2 = label2
        outlines.append(f"{marg2:<10}{np.mean(data2):<10.2g}{np.std(data2):<10.2g}{np.std(data2)/np.mean(data2):<10.2g}{min(data2):<10.2g}{max(data2):<10.2g}")
    
    for l in outlines:
        print(l)
    print()
    
def get_cds_stats(gm):
    
    all_cds=[]
    for chr in gm.gene_map.chromes:
        cds_inds = [(d.get("start"), d.get("end")) for d in gm.gene_map.fetch(chr, 0, gm.gene_map.max_indices[chr], features=  ["CDS"])]
        chunksz = int(np.sqrt(len(cds_inds)))
        nchunks = len(cds_inds)//chunksz
        
        for n in range(nchunks):
            
            cds_lens= [b-a for a,b in cds_inds if b-a > 0]
            seen_inds = set()
            nf = 0
            for (st, en), l in zip(cds_inds, cds_lens):
                if any([not (st > sen or sst > en) for sst, sen in seen_inds]):
                    continue
                else:
                    seen_inds.add((st, en))
                    all_cds.append(l)
                nf +=1
                if nf >= (n+1)*chunksz:
                    break
            cds_inds = cds_inds[nf:]
            cds_lens = cds_lens[nf:]
                
    
    meanlen = np.mean(all_cds)
    sdlen = np.std(all_cds)
    minlen = np.min(all_cds)
    maxlen = np.max(all_cds)
    n = len(all_cds)
    
    print("in units of bp:")
    print(f"Mean: {meanlen:0.3f}")
    print(f"SD: {sdlen:0.3f}")
    print(f"Min: {minlen:0.0g}")
    print(f"Max: {maxlen:0.0g}")
    print()
    print("in units of amino acids:")
    print(f"Mean: {meanlen/3:0.3f}")
    print(f"SD: {sdlen/3:0.3f}")
    print(f"Min: {minlen/3:0.0g}")
    print(f"Max: {maxlen/3:0.0g}")
    all_cds = np.array(all_cds)
    
    lo, hi = np.percentile(all_cds, [1, 99])
    
    f, ax = plt.subplots()
    ax.hist(all_cds, bins = 100, range = (lo, hi))
    ax.set_title("Length of CDS regions over all chromosomes")
    txt1 = f"n={len(all_cds):2.1e}\nseq len > 0\nno overlap"
    ax.text( 0.8, 0.88,txt1, transform = ax.transAxes, ha = 'right', va = "center")
    f.savefig("./data/outputs/cds_hist_filt.png")
        
    f, ax = plt.subplots()
    ax.hist(all_cds/3, bins = 100, range = (lo/3, hi/3))
    ax.set_title("Length of CDS regions over all chromosomes (in AA)")
    ax.text( 0.8, 0.88,txt1, transform = ax.transAxes, ha = 'right', va = "center")
    f.savefig("./data/outputs/cds_aa_hist_filt.png")
    return all_cds

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
    
    chrdata = {}
    
    txs = []
    
    ng = 0
    for chr in chrs:
        
        genes = gm.gene_map.fetch(chr, start=0, features = ("gene",))
        print(f"analyzing chrom {chr}")
        genedata = {}
        
        for g in genes:
            freqs = {aa:{cd:0 for cd in cds} for aa, cds in aa_codon_map.items()}
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
                freqs[aa][cd] += 1
            
            genedata[gene_name] = freqs
            
            ng +=1
            
            if ng % 100 == 0:
                print(f"analyzed gene {ng} {gene_name}, assembled transcript {txid} with number cdss {len(cdss)} and {len(cdsseq)} bps")
            
            if max_num_genes and ng >= max_num_genes:
                break
        
        chrdata[chr] = genedata
        if max_num_genes and ng >= max_num_genes:
            break
    
    return chrdata

def format_value(val, width = None):
    
    if isinstance(val, str):
        fcode = ""
    elif isinstance(val, int):
        fcode = "0g"
    elif isinstance(val, float):
        if int(val) == val:
            val = int(val)
            fcode = "0g"
        elif (val > 0 and abs(val) < 1e-5) or abs(val) > 1e5:
            fcode = "0.3e"
        else:
            fcode= "0.3f"
    else:
        fcode = ""
    
    if width is not None:
        wdthcode = "<" + str(width)
    else:
        wdthcode = ""
    
    fstr = fcode
    outstr = format(format(val, fstr), wdthcode)
    # if len(outstr) > width:
    #     outstr = outstr[:width-3] + ".. "
    return outstr
    
    
def write_codon_frequency_stats(chrdata, fname):
    
    fdir = Path("./data/outputs/codon_frequencies")
    fpath = fdir / fname
    fpath_summ = fpath.with_stem(fpath.stem + "_summary")
    
    rawlines_chr = {}
    summlines = []
    
    fmtr = lambda s: format_value(s, width = 16)
    
    header = ["Codon", "AA", "Chr", "Gene","Frequency", "Normalized Frequency"]
    header_summ = ["Codon", "AA", "Chr","Mean","SD","Min","Max","n"]
    
    chrs = list(chrdata.keys())
    print(f"chrs: {", ".join(chrs)}")
    
    codons = list(itertools.chain.from_iterable([get_aa_codons(aa) for aa in ORDER_AA]))
    ncdns = len(codons)
    print(f"codons: {", ".join(codons)}")
    
    nchrs = len(chrs)
    num_genes = [len(gd) for chr, gd in chrdata.items()]
    max_num_genes = max(num_genes)
    num_all_genes = sum(num_genes)
    all_genes = itertools.chain.from_iterable([gd.keys() for chr, gd in chrdata.items()])
    print(f"{num_all_genes} genes: {", ".join([g for g in all_genes if g][:20])} ...")
    
    alldata = np.zeros((ncdns, nchrs, max_num_genes))
    
    for codon in codons:
        cdnind = codons.index(codon)
        aa = CODON_TABLE[codon]
        npercdn = 0
        
        for chr in chrs:
            nchr = chrs.index(chr)
            genedata = chrdata[chr]
            num_genes = len(genedata)
            nperchr = 0
            
            rawlines = rawlines_chr.get(chr, [])
            # rawlines.append([codon, aa]+["-------"]*4)
            
            for ig, gene in enumerate(genedata):
                aadata = genedata[gene]
                
                cdndata = aadata[aa] # cdn -> freqs
                ns = len(cdndata)
                norm = sum(cdndata.values())
                f = cdndata[codon]
                if f == 0:
                    continue
                
                fnorm = f / norm if norm > 0 else -1
                
                alldata[cdnind, nchr, ig] = f
                
                rawlines.append([codon, aa, chr, gene, f, fnorm])
                nperchr += ns
            
            rawlines_chr[chr] = rawlines
            
            alldata_slc = alldata[cdnind, nchr]
            mean = np.mean(alldata_slc[:num_genes])
            sd= np.std(alldata_slc[:num_genes])
            minv = np.min(alldata_slc[:num_genes])
            maxv = np.max(alldata_slc[:num_genes])
            
            summlines.append([codon, aa, chr, mean, sd, minv, maxv, nperchr])
            npercdn += nperchr
        
        alldata_slc = alldata[cdnind]
        mean = np.mean(alldata_slc[alldata_slc > 0])
        sd= np.std(alldata_slc[alldata_slc > 0])
        minv = np.min(alldata_slc[alldata_slc > 0])
        maxv = np.max(alldata_slc[alldata_slc > 0])

        summlines.append([codon, aa, "------",  mean, sd, minv, maxv, npercdn])
    
    for chr in chrs:
        rawlines = rawlines_chr[chr]
        
        rawdf = pd.DataFrame(rawlines, columns = header, dtype= "object")
        rawdf_nostop = rawdf[rawdf.loc[:,"AA"] != "*"].sort_values(["AA","Gene"])
        rawdf_stop = rawdf[rawdf.loc[:,"AA"] == "*"].sort_values(["Gene"])
        rawdf_full = pd.concat((rawdf_nostop, rawdf_stop), axis = 0)
        
        fpath_chr = fpath.with_stem(fpath.stem + f"_chr{chr}")
        with open(fpath_chr, "w+") as f:
            f.write("".join(map(fmtr,header))+ "\n")
            # for rowvals in rawlines:
            for i,row in rawdf_full.iterrows():
                
                line = "".join(map(fmtr, row))
                f.write(line + "\n")
    
    with open(fpath_summ, "w+") as f:
        f.write("".join(map(fmtr,header_summ)) + "\n")
        
        for rowvals in summlines:
            line = "".join(map(fmtr, rowvals))
            f.write(line + "\n")
    
    return alldata

def parse(value):
    res = value
    try:
        res= float(value)
        if res == int(res):
            res = int(res)
    except:
        try:
            res = int(value)
        except:
            pass
    return res

def load_codon_frequency_data(chr, gene = ""):
    
    fdir = Path("./data/outputs/codon_frequencies")
    fname = f"codon_freqs_chr{chr}.txt"
    fpath = fdir / fname
    if not fpath.exists():
        return []
    
    data = []
    
    with open(fpath, "r") as f:
        gene_ind = 3
        for i,line in enumerate(f):
            linevals = [v.strip() for v in line.split("  ") if v.strip()]
            if i == 0:
                header = linevals
                gene_ind= header.index("Gene")
            if gene and linevals[gene_ind] != gene:
                continue
            linevals = [parse(v) for v in linevals]
            data.append(linevals)
    return data
    

def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    
    # aa ='Q'
    # cod, adj = get_adjacent_codons(aa)
    # print(cod)
    # print(adj)
    
    # # print("only A<->G or C<->T")
    
    # aa ='N'
    # cod, adj = get_adjacent_codons(aa, mutations = "N")
    # print(cod)
    # print(adj)
    
    # mutation_map = make_mutation_adj_map(muts = ["N","YR","SW","KM","B","D","H","V"])
    # mutation_map = make_mutation_adj_map(muts = ["YR"])
    
    # print(mutation_map)
    # c2a, a2adj = mutation_map.get("YR")
    
    # start = "AAA"
    # nsteps = 20
    
    # path, cpath = random_walk_mutation_map(c2a, a2adj, start, nsteps)
    

    gm = load_genome()
    genemap = gm.gene_map

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
    
    # prote_vc = draw._protein_metavocab
    prote_vc = 'ABCDEF'
    # gene_name = "ERC1"
    # chr = 12
    # gene_name = "AGRN"
    # chr = 1
    gene_name = "ZNF804A"
    chr = 2

    ran = 50000
    gene_feats = [f for f in genemap.fetch(chr, start=0, features =("gene",))]
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
    
    display_alt_codon_entropy(aacodedict)
    
    # chrs = []
    # chrdata = get_codon_frequency_stats(gm, chrs = chrs)
    # alldata = write_codon_frequency_stats(chrdata, "codon_freqs.txt")
    
    # print([parse(v) for v in ["1","1.1","hello","Gene"]])
    
    # cdn_freq_Y = load_codon_frequency_data("Y")
    # print(tabulate(cdn_freq_Y))
    
    cdn_freq_znf = load_codon_frequency_data(2, gene = "ZNF804A")
    print(tabulate(cdn_freq_znf))
    
    
    
    # freqs_norm, freqs_ct, txs = get_codon_frequency_stats(gm, chrs = ["12"])
    chrdata = get_codon_frequency_stats(gm, chrs = ["12"])
    display_alt_codon_entropy(freqs_norm, freqs_ct, fname = "codon_entropy_chr12.txt")
    
    return
    codons = "".join("".join(get_aa_codons(aa)) for aa in ORDER_AA)
    
    display_alt_codons(codons, use_color = True)
    
    # cdgcs = {aa:{cd:cd.count("G") + cd.count("C")} for cd, aa in CODON_TABLE.items()}
    cdgcs = {aa:[cd.count("G") + cd.count("C") for cd in get_aa_codons(aa)] for aa in ORDER_AA}
    
    for aa, cdgcs in cdgcs.items():
        # aa = CODON_TABLE.get(cd)
        print(aa, cdgcs)
    
    
    # freqs_norm, freqs_ct, txs = get_codon_frequency_stats(gm, chrs = [])
    # display_alt_codon_entropy(freqs_norm, freqs_ct,  fname = "codon_entropy.txt")
    
    # for aa, fdict in freqs.items():
    #     dist = draw.scalar_plot_distribution(fdict, bg_color = bg, fg_color = fg-12)[0]
    #     print(aa, dist)
        # for i, f in fdict.items():
        
    # print()
    # print(f"AGRN codon metalanguage with {len(aacode)} symbols:")
    # colseq = color_protein_metasequence(aacode, prote_vc, [f"\x1b[38;5;{random.randint(30,220)}m" for c in range(6)])
    # colseq = color_protein_metasequence(aacode, prote_vc)
    # print("".join(aacode))
    # disp_seq("".join(aacode))
    # write_seq("".join(aacode), "codon_seq.txt")
    # check_written_seq("codon_seq.txt")
    
    
    
    

if __name__=="__main__":
    main()

