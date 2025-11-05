

import numpy as np
import random
from matplotlib import pyplot as plt

import itertools
from tabulate import tabulate

import os
from pathlib import Path
import pandas as pd

# from scipy.stats import wasserstein
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE

from ggene.genomemanager import GenomeManager

def get_random_sequence(seq_len, minval = -100, maxval = 100, var = 1, bias = 0, nint = 1):
    
    rand_seq_d = [(2*random.random()-1)*var+bias for i in range(seq_len)]
    
    for n in range(nint):
        rand_seq = [sum(rand_seq_d[:i]) for i in range(seq_len)]
        rand_seq_d = rand_seq
    
    _max= max(rand_seq)
    _min = min(rand_seq)
    _ran = _max - _min
    
    tgt_ran = maxval - minval
    
    rand_seq_sc = [tgt_ran*(r-_min)/_ran+minval for r in rand_seq]
    return rand_seq_sc

def test_ruler():
    ruler_params = [
        (0, 500, 50, 5, 5, 0, ".0f"),
        (0, 500, 250, 10, 5, 0, ".0f"),
        (0, 1, 250, 10, 5, 0, ".1f"),
        (0, 1, 250, 10, 5, 2, ".1f"),
        (0, 10, 249, 11, 6, 4, ".1f"),
        (0, 14, 215, 10, 7, 3, ".1f"),
        (-0.48, 0.51, 243, 3, 11, 3, ".2f"),
    ]
    
    for prms in ruler_params:
        rseq = get_random_sequence(prms[2])
        sctxt = draw.scalar_to_text_nb(rseq, bit_depth = 8)
        testruler = draw.make_ruler(*prms)
        print(sctxt[0])
        print(testruler)
        print()

def load_genome():
    return GenomeManager()

def get_chr_gc(gm, chr, chunksz = 10e6, start = 1e6, allow_aliases = False):
    
    def seq_gc(seq, feats):
        gc = (seq.count("G") + seq.count("C"))
        if allow_aliases:
            gc += seq.count("R") / 2 + seq.count("Y") / 2 + seq.count("N")/2
        gc /= len(seq)
        return gc
    
    gcs, starts = get_chromosomal_quantity(gm, chr, seq_gc, chunksz = chunksz, start = start)
    
    return gcs, starts

def get_chromosomal_quantity(gm, chr, seq_fn, chunksz = 10e6, start = 1e6):
    
    chrmax = gm.gene_map.max_indices.get(str(chr))
    num_chunks = int((chrmax - start)//chunksz)
    step = chunksz
    starts = []
    qts = []
    
    chrstr = str(chr)
    
    for n in range(num_chunks):
        end = start + step
        seq = gm.get_sequence(chrstr, start, end)
        ufeats = gm.get_all_annotations(chrstr, start, end)
        feats = [uf.to_dict() for uf in ufeats]
        if len(seq) < 1:
            start = end
            continue
        # gc = (seq.count("G") + seq.count("C"))
        # if allow_aliases:
        #     gc += seq.count("R") / 2 + seq.count("Y") / 2 + seq.count("N")/2
        # gc /= len(seq)
        qt = seq_fn(seq, feats )
        starts.append(start)
        qts.append(qt)
        start = end
    return qts, starts

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
    
    return "".join(full_seq), tx_feat

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
            linevals = [v.strip() for v in line.split(" ") if v.strip()]
            if i == 0:
                header = linevals
                gene_ind= header.index("Gene")
                continue
            if gene and linevals[gene_ind] != gene:
                continue
            linevals = [parse(v) for v in linevals]
            data.append(linevals)
    return header, data

def make_codon_vecs(data, genes = None):
    
    exclude = ["M","W","*"]
    cd_exclude = itertools.chain.from_iterable([bio.get_aa_codons(aa) for aa in exclude])
    cdi_exclude= sorted([bio.get_codon_index(cd) for cd in cd_exclude], reverse = True)
    ncodons = len(CODON_TABLE)
    codons = [bio.index_to_codon(icdn) for icdn in range(ncodons) if icdn not in cdi_exclude]
    cd_vecs = {}
    
    for row in data:
        try:
            cd, aa, chr, gene, f, fn = row
        except:
            cd, aa, chr, genef, fn = row
            gene, f = genef[:-1], int(genef[-1])
            continue
        
        if genes and not gene in genes:
            continue
        if not gene in cd_vecs:
            cd_vecs[gene] = [0]*ncodons
        
        cdind = bio.get_codon_index(cd)
        cd_vecs[gene][cdind] = fn
    
    arr = np.zeros((len(cd_vecs), ncodons - len(cdi_exclude)))
    genes = []
    for i, g in enumerate(cd_vecs):
        vec = cd_vecs[g]
        for iexc in cdi_exclude:
            del vec[iexc]
        arr[i, :] = vec
        genes.append(g)
        
    return arr, genes, codons

def do_pca(cdvecs, genes, codons, standardize = False):
    
    if standardize:
        mean = np.mean(cdvecs, axis = 0)
        sd = np.std(cdvecs, axis = 0)
        cdvecs = (cdvecs - mean) / sd
    
    pca = PCA()
    pca.fit(cdvecs)
    
    return pca, pca.transform(cdvecs)[:, :20]

def cluster_kmeans_range(pca_vecs, min_k = 2, max_k = 20):
    
    lbls = []
    inertias = []
    sils = []
    for k in range(min_k, max_k+1):
        kmeans = KMeans(n_clusters=k, random_state = 167)
        lbl = kmeans.fit_predict(pca_vecs)
        lbls.append(lbl)
        inertias.append(kmeans.inertia_)
        sils.append(silhouette_score(pca_vecs, lbl))
    
    return lbls, inertias, sils

def cluster_kmeans(pca_vecs, k):
    
    kmeans = KMeans(n_clusters=k, random_state = 167)
    lbl = kmeans.fit_predict(pca_vecs)
    inertia = kmeans.inertia_
    sil = silhouette_score(pca_vecs, lbl)
    
    return kmeans, lbl, inertia, sil

def get_label_centroids(data, lbls, num_labels):
    
    cents = [np.zeros(data[0].shape) for lbl in range(num_labels)]
    nums = [0 for lbl in range(num_labels)]
    
    for d, lbl in zip(data, lbls):
        
        cents[lbl] += d
        nums[lbl] += 1
    
    for lbl in range(num_labels):
        cents[lbl]/=nums[lbl]
        
    vars = [np.zeros(data[0].shape) for lbl in range(num_labels)]
    for d, lbl in zip(data, lbls):
        vars[lbl] += np.power((d - cents[lbl]),2)
    
    for lbl in range(num_labels):
        vars[lbl]/=nums[lbl]
    
    return np.array(cents), np.array(vars)

def tabulate_labels(centroids, k, codons):
    tab = []
    for aa in ORDER_AA:
        if aa in "MW*":
            continue
        for cdn in bio.get_aa_codons(aa):
            row = [aa, cdn]
            for i,cent in enumerate(centroids):
                row.append(cent[codons.index(cdn)])
            tab.append(row)
        tab.append([" "])
    print(tabulate(tab, headers = range(1, k+1)))
    return tab

def show_labeled_codons(centroids, codons, vars = None):
        
    do_flip = True
        
    dists = []
    for aa in ORDER_AA:
        if aa in "MW*":
            continue
        for cdn in bio.get_aa_codons(aa):
            icdn = codons.index(cdn)
            dist = draw.scalar_to_text_nb(centroids[:,icdn], minval = 0, bit_depth = 8)
            dists.append([aa, cdn, dist[0]])
            if vars is not None:
                vdist = draw.scalar_to_text_nb(vars[:,icdn], minval = 0, bit_depth = 8, fg_color = 65, flip = do_flip)
                dists.append([" ", " ", vdist[0]])
        dists.append([" "])
    
    headers = ["AA","Codon","Label"]
    print(tabulate(dists, headers = headers))
    return dists, headers

def show_label_vectors(centroids, k, codons, vars = None, samples = None, sample_labels = None):
    # sample_labels is tuple of kmeans_label, real_label?
    
    do_flip = True
    
    do_vars = True
    if vars is None:
        do_vars = False
        vars = centroids
    
    dists = []
    for lbl in range(k):
        centvec = []
        varvec = []
        for aa in ORDER_AA:
            if aa in "MW*":
                continue
            cdns = bio.get_aa_codons(aa)
            scale = len(cdns) / 2
            for cdn in cdns:
                icdn = codons.index(cdn)
                centvec.append(scale*centroids[lbl, icdn])
                varvec.append(vars[lbl, icdn])
            centvec.append(0)
            varvec.append(0)
        maxval = 1.0
        cdist = draw.scalar_to_text_nb(centvec, minval = 0, maxval = maxval, bit_depth = 8)
        dists.append([lbl, "Cent.", cdist[0]])
        if do_vars:
            vdist = draw.scalar_to_text_nb(varvec, minval = 0, bit_depth = 8, flip = do_flip, fg_color = 65)
            dists.append([" ", "Var", vdist[0]])
    
    if samples is not None:
        dists.append([" "])
        if sample_labels is None:
            sample_labels = [("?", f"Sample{n}") for n in range(len(samples))]
        for n in range(len(samples)):
            kmlbl, infolbl = sample_labels[n]
            sampvec = []
            for aa in ORDER_AA:
                if aa in "MW*":
                    continue
                cdns = bio.get_aa_codons(aa)
                for cdn in cdns:
                    icdn = codons.index(cdn)
                    sampvec.append(samples[n, icdn])
                sampvec.append(0)
            maxval = 1.0
            sdist = draw.scalar_to_text_nb(sampvec, minval = 0, maxval = maxval, bit_depth = 8)
            dists.append([kmlbl, infolbl, sdist[0]])
            # dists.append([" "])
    
    aahdr = []
    for aa in ORDER_AA:
        if aa in "MW*":
            continue
        ncds = len(bio.get_aa_codons(aa))
        aahdr.append(aa)
        aahdr.append(" "*(ncds))
    aahdr = "".join(aahdr)
    
    headers = ["Label","Data", aahdr]
    print(tabulate(dists, headers = headers))
    return dists, headers

def group_genes(genes, labels, k):
    
    grouped = [[] for i in range(k)]
    
    for g, lbl in zip(genes, labels):
        grouped[lbl].append(g)
    
    return grouped

def get_scaled_dists(centroids, vars, k):
    
    res = {}
    res2 = {}
    
    for i in range(k):
        ci =centroids[i]
        vi = vars[i]
        dvi2 = np.dot(vi, vi)
        for j in range(i+1, k):
            cj = centroids[j]
            vj = vars[j]
            dvj2 = np.dot(vj, vj)
            
            delta = cj - ci
            dist2_unscaled = np.dot(delta, delta)
            res2[(i,j)] = np.sqrt(dist2_unscaled)
            
            vipdelta = np.dot(delta, vi) / dist2_unscaled
            vjpdelta = np.dot(delta, vj) / dist2_unscaled
            vsum = vipdelta - vjpdelta
            vsumnorm = np.sqrt(np.dot(vsum, vsum))
            
            res[(i,j)] = np.sqrt(dist2_unscaled / vsumnorm)
    
    return res, res2

def label_genes(chr, pca, kmeans, genes = []):
    
    all_cddata = []
    samples = []
    sample_labels = []
    
    if genes:
        for gene in genes:
            hdr,cddata = load_codon_frequency_data(chr, gene = gene)
            all_cddata.append(cddata)
        all_cddata = np.vstack(all_cddata)
    else:
        hdr,cddata = load_codon_frequency_data(chr)
        all_cddata = cddata

    cdvecs, genes, codons = make_codon_vecs(all_cddata)
    cd_pca = pca.transform(cdvecs)[:,:20]
    cdvecs_lbl = kmeans.predict(cd_pca)
    samples = cdvecs
    for lbl, gene in zip(cdvecs_lbl, genes):
        sample_labels.append((lbl, gene))
    return samples, sample_labels

def gather_gene_data(gm, gene_names, chr, genes = []):
    
    if not genes:
        gene_feats = [f for f in gm.gene_map.fetch(chr, start=0, features =("gene",))]
        genes = [g for g in gene_feats if g.get("info").get("gene_name") in gene_names]
    
    headers = ["Gene Name", "CG", "Start", "Length"]
    gene_data = []
    gns = []
    for gene in genes:
        if not gene:
            continue
        ginfo = gene.get("info")
        gn = ginfo.get("gene_name")
        if not gn:
            continue
        seq, txfeat = get_transcript_cds(gm, gn, chr, gene = gene)
        gc = (seq.count("G") + seq.count("C"))/len(seq)
        start = gene.get("start")
        length = gene.get("end") - start
        gene_data.append([gc, start, length])
        gns.append(ginfo.get("gene_name"))
    return np.array(gene_data, dtype = 'float32'), headers, gns

def what_is_this(gm, chr, pos):
    feats = gm.get_all_annotations(chr, pos, pos+1)
    print(f"{len(feats)} features located at {chr}:{pos}")
    for ufeat in feats:
        feat = ufeat.to_dict()
        finf = feat.get("attributes",{})
        print(f"{feat.get("type")}: {finf.get("name","?")}")
        for k in feat.keys():
            if k in ["type", "attributes"]:
                continue
            print(f"  {k}: {feat.get(k)}")
        print("  attributes:")
        for k in finf:
            print(f"    {k}: {finf[k]}")
    return feats

def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    gm = load_genome()
    
    # what_is_this(gm, 13, 102760340)
    
    # return
    
    test_chr = "X"
    
    hdr,cddata_chr11 = load_codon_frequency_data(test_chr)
    
    cdvecs, genes, codons = make_codon_vecs(cddata_chr11)
    
    pca, pcares = do_pca(cdvecs, genes, codons, standardize = False)
    
    print(pcares.shape)
    
    min_k = 2
    max_k = 20
    klbls, ins, sils = cluster_kmeans_range(pcares, min_k = min_k, max_k = max_k)
    
    ins_txt = draw.scalar_to_text_nb(ins, minval = 0, add_range = True)
    ins_txt, _ = draw.add_ruler(ins_txt, min_k, max_k, num_labels = 4, ticks = 0, minor_ticks = 100, fstr= ".0f")
    sils_txt = draw.scalar_to_text_nb(sils, minval = 0, add_range = True)
    sils_txt, _ = draw.add_ruler(sils_txt, min_k, max_k, num_labels = 4, ticks = 0, minor_ticks = 100, fstr = ".0f")
    
    print("inertias:")
    for r in ins_txt:
        print(r)
    # print(str(min_k) + "╵"*(max_k- min_k - 1) + str(max_k))
    print()
    
    print("silhouettes:")
    for r in sils_txt:
        print(r)
    # print(str(min_k) + "╵"*(max_k- min_k - 1) + str(max_k))
    print()
    
    
    k = 7
    kmeans, klbls, inrt, sil = cluster_kmeans(pcares, k)
    """
    notes on clusters for chr 11:
    5 leans G/C, and GC in particular. loves CGG, CGC, and AGC. uniquely loves AGC Ser. along with 2, loves Leu CTG and CTC. however, 5 takes Ala GCC over GCG!
    
    6 doesn't like GC and evades it in A, Q, R, S, and T codons. uniquely uniform in C, D, E, F, G codons
    
    3 and 4 are common in their C, D codons (preferring TGT, GAC over TGC, GAT), but differ in their N codon (3 prefers AAC, 4 is uniform) and in S, 4 leans G/C
    
    """
    
    
    

    
    cents, vars = get_label_centroids(cdvecs, klbls, k)
    
    num_show = 30
    genes_grp = group_genes(genes, klbls, k)
    print("Grouped genes:")
    for i, row in enumerate(genes_grp):
        print(i, ", ".join(row[:num_show]))
    print()
    
    # _ = tabulate_labels(cents, k, codons)
    
    # rows, hdr = show_labeled_codons(cents, codons, vars = vars)
    
    rows2, hdr2 = show_label_vectors(cents, k, codons, vars = vars)
    
    cpgdata = []
    for aa in ORDER_AA:
        if aa in "MW*":
            continue
        cdns = bio.get_aa_codons(aa)
        scale = len(cdns) / 2
        for cdn in cdns:
            cpgdata.append(int("CG" in cdn))
        cpgdata.append(0)
        
    cpgres = draw.scalar_to_text_8b(cpgdata)
    for r in cpgres:
        print(" "*17 + r)
    for i in range(3):
        row = []
        for aa in ORDER_AA:
            if aa in "MW*":
                continue
            cdns = bio.get_aa_codons(aa)
            row.extend([cdn[i] for cdn in cdns])
            row.append(" ")
        print(" "*17 + "".join(row))
    print()
    
    # for icdn in range(64):
    #     print(icdn, bio.index_to_codon(icdn), CODON_TABLE.get(bio.index_to_codon(icdn)))
    
    
    test_chr = 13
    
    samples, sample_labels = label_genes(test_chr, pca, kmeans)
    cats = [lbl for lbl, gn in sample_labels if gn != 'None']
    gns = [gn for _, gn in sample_labels if gn != 'None']
    
    gene_dicts = gm.gene_map.fetch(test_chr, start = 0, features = ("gene",))
    gdict = {g.get("info",{}).get("gene_name","?"): g for g in gene_dicts}
    genes = [gdict.get(gn) for gn in gns if gn in gdict]
    
    genedata, hdrs, gene_names = gather_gene_data(gm, [], test_chr, genes = genes)
    genedata_lbl = [[gn] + list(gd) for gn, gd in zip(gene_names, genedata)]
    
    print(tabulate(genedata_lbl[:20], headers = hdrs))
    
    max_disp = 256
    chunksz = int(32e3)
    start = int(16e6)
    ran = gm.gene_map.max_indices[str(test_chr)] - start
    num_chunks = ran//chunksz
    
    
    def seq_cpg(seq, feats):
        return seq.count("CG")
    
    def seq_genes(seq, feats):
        return sum([1 for f in feats if f.get("type","") == "gene"])
    
    def seq_exons(seq, feats):
        return sum([1 for f in feats if f.get("type","") == "exon"])
    
    def seq_len_cds(seq, feats):
        if not feats:
            return 0
        return sum([f.get("end") - f.get("start") for f in feats if f.get("type","") == "CDS"])/len(feats)
    
    def seq_feats(seq, feats):
        return sum([1 for f in feats])
    
    def seq_max_run(seq, feats):
        # scale = round(np.sqrt(len(seq)))
        scale = 256
        
        maxrun = 0
        for strt in range(0, len(seq), scale):
            subseq = seq[strt:strt+scale]
            runs, _, _ = seqs.process.correlate_longest_subseq(subseq, reverse_complement(subseq), scale = 32)
            maxrun = max(maxrun, max(runs))
        return maxrun
    
    # seq_fn = seq_cpg
    # qt_name = "CpG"
    # seq_fn = seq_genes
    # qt_name = "ngenes"
    # seq_fn = seq_len_cds
    # qt_name = "len(cds)"
    seq_fn = seq_max_run
    qt_name = "Run"
    
    # gcs, starts = get_chr_gc(gm, test_chr, chunksz = chunksz, start = start, allow_aliases = True)
    qts, starts = get_chromosomal_quantity(gm, test_chr, seq_fn, chunksz = chunksz, start = start)
    
    print(f"{qt_name} range: {min(qts):0.2f}-{max(qts):0.2f}")
    
    # minval = 0.28
    # maxval = 0.65
    minval = 0
    maxval = max(qts)
    num_disp_chunks = int(num_chunks // max_disp)
    
    print("")
    print()
    
    print(f"{qt_name} for chr{test_chr}")
    for i in range(num_disp_chunks):
        seq_ran = ran // (num_disp_chunks)
        seq_start = start + i*seq_ran
        
        gcs_crop = qts[i*max_disp:(i+1)*max_disp]
        if len(gcs_crop) == 0:
            break
        
        gc_bars = draw.scalar_to_text_nb(gcs_crop, minval = minval, maxval = maxval, bit_depth = 24, add_range = True)
        gc_bars, dists = draw.add_ruler(gc_bars, seq_start, seq_start + seq_ran, num_labels = 10, ticks = 2, minor_ticks = 5, genomic = True)
        print(f"Section {i+1} | {seq_start/1e6:.1f}M - {(seq_start+seq_ran)/1e6:.1f}M (major = {dists[1]/1e3:2.0f}k, minor = {dists[2]/1e3:2.0f}k)")
        for r in gc_bars:
            print(r)
        print()
    
    return
    
    centr_pos = 17.7e6
    bands = [
        (17.e6, 18.9e6),
        (54.7e6, 59e6),
        (68.1e6, 72.8e6),
        (78.5e6, 87.1e6),
        (89.4e6, 94.4e6),
        (101.1e6, 104.2e6),
        (106.4e6, 109.6e6),
    ]
    
    f,ax = plt.subplots(1,1)
    ax.scatter(cats, genedata[:,0], marker = ',', s = 1)
    ax.set_title(f"{qt_name} vs Codon label")
    f.savefig(f"./data/outputs/cg_vs_codonlbl_chr{test_chr}.png")
    
    f2,ax2 = plt.subplots(1,1)
    ax2.scatter(cats, genedata[:,1], marker = ',', s = 1)
    xlim = ax2.get_xlim()
    ax2.plot([xlim[0], xlim[1]], [centr_pos, centr_pos], 'r--', linewidth = 1)
    for bstart, bend in bands:
        ax2.plot([xlim[0], xlim[1]], [bstart, bstart], 'k--', linewidth = 1)
        ax2.plot([xlim[0], xlim[1]], [bend, bend], 'k--', linewidth = 1)
    
    ax2.set_title("Start vs Codon label")
    f2.savefig(f"./data/outputs/start_vs_codonlbl_chr{test_chr}.png")
    
    f3,ax3 = plt.subplots(1,1)
    ax3.scatter(cats, genedata[:,2], marker = ',', s = 1)
    ax3.set_title("Length vs Codon label")
    ax3.set_yscale("log")
    f3.savefig(f"./data/outputs/length_vs_codonlbl_chr{test_chr}.png")
    
    f3,ax3 = plt.subplots(1,1)
    ax3.scatter(genedata[:,0], genedata[:,1], marker = ',', s = 1)
    xlim = ax3.get_xlim()
    ax3.plot([xlim[0], xlim[1]], [centr_pos, centr_pos], 'r--', linewidth = 1)
    ax3.set_title(f"{qt_name}  vs Start")
    f3.savefig(f"./data/outputs/cg_vs_start_chr{test_chr}.png")
    
    
    # genedata = [
    #     (1, "AGRN"),
    #     (2, "ZNF804A"),
    #     (1, "HES4"),
    #     (16, "IRX3"),
    #     (2, "CREB1"),
    #     (11, "BDNF"),
    #     (17, "SLC6A4"),
    #     (17, "BRCA1")
    # ]
    
    # samples = []
    # sample_labels = []
    # for chr, gene in genedata:
    #     hdr,cddata_agrn = load_codon_frequency_data(chr, gene = gene)
    
    #     agrnvec, genes, codons = make_codon_vecs(cddata_agrn)
    #     agrn_pca = pca.transform(agrnvec)[:,:20]
    #     agrn_lbl = kmeans.predict(agrn_pca)[0]
    #     samples.append(agrnvec)
    #     sample_labels.append((agrn_lbl, gene))
    # samples = np.vstack(samples)
    
    # rows2, hdr2 = show_label_vectors(cents, k, codons, vars = vars, samples = samples, sample_labels = sample_labels)
    
    # dists, dists_unscaled = get_scaled_dists(cents, vars, k)
    # for i, j in dists:
    #     print(i, j, format((dists[(i,j)]), "0.5f"), format((dists_unscaled[(i,j)]), "0.5f"))
    
    # fdir = Path("./data/outputs/codon_frequencies")
    # fp1 = fdir / "labeled_codons_chr11_k7.txt"
    # fp2 = fdir / "codon_vectors_k7_samples.txt"
    
    # with open(fp1, "w+") as f:
    #     f.write("".join([format(v, "<8") for v in hdr]) + "\n")
    #     for line in rows:
    #         f.write(draw.scrub_ansi("".join([format(v, "<8") for v in line])) + "\n")
    
    # with open(fp2, "w+") as f:
    #     f.write("".join([format(v, "<8") for v in hdr2]) + "\n")
    #     for line in rows2:
    #         f.write(draw.scrub_ansi("".join([format(v, "<8") for v in line])) + "\n")
    
    
    
    
if __name__=="__main__":
    main()

