

import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

from tabulate import tabulate

import numpy as np

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import bio, process, find, align, heal
from ggene.seqs.bio import ALIASES, ALIASES_REV

@dataclass
class IntronData:
    gene:Dict[str,Any]
    gene_index: int
    exon1:Dict[str,Any]
    exon_index: int
    exon2:Dict[str,Any]
    intron_seq:str
    
    @property
    def gene_name(self):
        return self.gene.get("info",{}).get("gene_name","?")

    @property
    def seq(self):
        return self.intron_seq
    
    @property
    def start(self):
        return self.exon1.get("end", -1)
    
    @property
    def end(self):
        return self.exon1.get("start", -1)

    @property
    def chr(self):
        return self.gene.get("chrom", "")

    def __repr__(self):
        return f"IntronData({self.gene_name}, bps={len(self.intron_seq)})"
    
    def __str__(self):
        return f"intron {self.exon_index} at {self.chr}:{self.start}-{self.end} on gene {self.gene_name} (gene index {self.gene_index}) with {len(self.intron_seq)} bases, "

    def to_tuple(self):
        return self.gene_name, self.gene.get("chrom", ""), self.exon_index

@dataclass
class SimilarityResult:
    exon_data:IntronData
    similarity:float
    cds_aa_seq:str

base_score = {    
    ("A","A"):1.0,
    ("G","G"):1.0,
    ("T","T"):1.0,
    ("C","C"):1.0,
    
    # G-U, common
    ("T","G"):0.75,
    ("G","T"):0.75,
    
    # G-A, common
    ("A","G"):0.5,
    ("G","A"):0.5,
    
    ("A","A"):0.0,
    ("T","T"):0.0,
    ("G","G"):0.0,
    ("C","C"):0.0,
    
    ("A","C"):0.0,
    ("T","C"):0.0,
    ("C","A"):0.0,
    ("C","A"):0.0,
}


def load_genome():
    return GenomeManager()

def get_intron_seq(gm, nchr, exon1, exon2):
    
    intron_seq = gm.get_sequence(nchr, exon1.get("end"), exon2.get("start"))
    
    if exon1.get("strand") == "-":
        intron_seq = bio.reverse_complement(intron_seq)
    
    return intron_seq

def get_feature_by_index(gm, nchr, feature_ind, start = 0, end = None, feature_type = 'gene'):
    
    features = (feature_type,)
    
    nf = 0
    feat = None
    for g in gm.gene_map.fetch(nchr, start, end=end, features = features):
        if nf < feature_ind:
            nf +=1
        else:
            feat = g
            break
    return feat, nf

def get_gene_by_index(gm, nchr, gene_ind):
    
    gene, ngs = get_feature_by_index(gm, nchr, gene_ind, feature_type = 'gene')
    if gene is None:
        gene_ind = gene_ind % ngs
        gene, ngs = get_feature_by_index(gm, nchr, gene_ind, feature_type = 'gene')
    
    return gene, ngs

def get_gene_by_name(gm, nchr, gene_name):
        
    features = ('gene',)
    
    nf = 0
    feat = None
    for g in gm.gene_map.fetch(nchr, 0, end=None, features = features):
        if not g.get("info",{}).get("gene_name","") == gene_name:
            nf+=1
            continue
        else:
            break
        
    return feat, nf

def get_intron(gm:GenomeManager, nchr, gene_index, exon_index):
    
    gene, ng = get_gene_by_index(gm, nchr, gene_index)
    exon1, nintr = get_feature_by_index(gm, nchr, exon_index, feature_type = "exon")
    exon2, nintr = get_feature_by_index(gm, nchr, exon_index + 1, feature_type = "exon")
    
    intron_seq = get_intron_seq(gm, nchr, exon1, exon2)
    intrdat = IntronData(gene, ng, exon1, nintr, exon2, intron_seq)
    
    # if len(intrdat.seq) < 1:
    #     print(f"empyt intron {intrdat}")
    #     return None
    
    return intrdat

def get_gene_introns(gm:GenomeManager, nchr, gene_index = None, gene_name = "", max_num = None, min_len = 256, max_len = 2048):
    
    if not gene_index:
        gene, ngs = get_gene_by_name(gm, nchr, gene_name)
    else:
        gene, ngs = get_feature_by_index(gm, nchr, gene_index, feature_type = 'gene')
    
    gene_index = ngs
    
    introns = []
    
    exind = 0
    done = False
    while not done:
        
        intr = get_intron(gm, nchr, gene_index, exind)
        
        if intr in introns:
            done = True
            break
        elif introns and intr.start < introns[-1].end:
            done = True
            break
        elif max_num and len(introns) > max_num:
            done = True
            break
        
        if not intr or len(intr.seq) < min_len or len(intr.seq) > max_len:
            exind += 1
        else:
            print(f"added {intr}")
            introns.append(intr)
        
        exind += 1
    
    return introns

def get_random_intron(gm:GenomeManager, min_len = None, max_len = None):
    
    done = False
    while not done:
        nchr = random.choice(list(gm.gene_map.max_indices.keys()))
        num_gene = random.randint(0, 500)

        gene, ng = get_gene_by_index(gm, nchr, num_gene)
        nexon1= random.randint(0, 10)
        exon1, nintr = get_feature_by_index(gm, nchr, nexon1, feature_type = "exon")
        exon2, nintr = get_feature_by_index(gm, nchr, nexon1 + 1, feature_type = "exon")
        intron_seq = get_intron_seq(gm, nchr, exon1, exon2)
        
        if min_len and len(intron_seq) < min_len:
            continue
        elif max_len and len(intron_seq) > max_len:
            continue
        break
        
    intrdat = IntronData(gene, ng, exon1, nintr, exon2, intron_seq)
    
    return intrdat

def score_sequences_corrs(seqa, seqb, topk = 5):
    min_len = min(len(seqa), len(seqb))
    corrs, rccorrs = process.correlate(seqa, seqb)
    topk_corrs = sorted(corrs)[:topk]
    topk_rccorrs = sorted(rccorrs)[:topk]
    return sum(topk_corrs)/topk, sum(topk_rccorrs)/topk

def score_sequences_runs(seqa, seqb, topk = 5, max_err = 16):
    min_len = min(len(seqa), len(seqb))
    runs, _, shifts, _ = process.correlate_longest_subseq_err(seqa, seqb, max_err)
    topk_runs = sorted([(r,s) for r,s in zip(runs, shifts) if s > 3], key = lambda k:-k[0])[:topk]    
    rcruns, _, rcshifts, _ = process.correlate_longest_subseq_err(seqa, bio.reverse_complement(seqb), max_err)
    topk_rcruns = sorted([(r,s) for r,s in zip(rcruns, rcshifts) if s > 3], key = lambda k:-k[0])[:topk]    
    return sum([r for r, s in topk_runs])/min_len/topk, sum([r for r, s in topk_rcruns])/min_len/topk

def score_sequences_algn(seqa, seqb):
    min_len = min(len(seqa), len(seqb))
    score = align.score_sequences(seqa, seqb)
    rcscore = align.score_sequences(seqa, bio.reverse_complement(seqb))
    return score/min_len, rcscore/min_len

def downsample_scores(scores):
    
    scores = np.array(scores)
    nr, nc = scores.shape
    
    ds_score = 0.5*scores[1:nr-1, 1:nc-1]
    ds_score += 0.25*scores[2:nr, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 2:nc]
    ds_score += 0.25*scores[0:nr-2, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 0:nc-2]
    
    return ds_score

def compare_introns(seqa, seqb, chunksz = 128, score_modes = ["alignment", "runs", "corrs"], resample = True):
    
    if len(seqb) > len(seqa):
        seqa, seqb = seqb, seqa
    
    q = 1
    if resample:
        q = 2
        chunksz = chunksz//2
    
    nchksa = len(seqa) // chunksz
    nchksb = len(seqb) // chunksz
    
    score_funcs = []
    if "alignment" in score_modes:
        score_funcs.append(score_sequences_algn)
    if "runs" in score_modes:
        score_funcs.append(score_sequences_runs)
    if "corrs" in score_modes:
        score_funcs.append(score_sequences_corrs)
    
    scores = [[[] for n in range(nchksb)] for sf in score_funcs]
    rcscores = [[[] for n in range(nchksb)] for sf in score_funcs]
    row_lbls = [f"b{n}" for n in range(nchksb)]
    
    for ia in range(nchksa):
        ssa = seqa[chunksz*ia:chunksz*(ia+q)]
        
        for ib in range(nchksb):
            
            ssb = seqb[chunksz*ib:chunksz*(ib+q)]
            if not ssa or not ssb:
                continue
                
            for i, sf in enumerate(score_funcs):
                score, rcscore = sf(ssa, ssb)
                scores[i][ib].append(score)
                rcscores[i][ib].append(rcscore)
    
    all_hms = [[] for r in range(len(scores[0])+3)]
    all_rchms = [[] for r in range(len(rcscores[0])+3)]
    
    for nsm,sm in enumerate(score_modes):
        
        if nsm > 0:
            row_lbls = []
        
        if resample:
            scores_rs = downsample_scores(scores[nsm])
            rcscores_rs = downsample_scores(rcscores[nsm])
        else:
            scores_rs = scores[nsm]
            rcscores_rs = rcscores[nsm]
        
        hm = draw.heatmap(scores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra")
        
        for nr in range(len(hm)):
            all_hms[nr].append(hm[nr])
        
        rchm = draw.heatmap(rcscores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra")
        for nr in range(len(rchm)):
            all_rchms[nr].append(rchm[nr])
    
    frm = "{:^32}       " * len(score_modes)
    for row in all_hms:
        if row:
            print(frm.format(*row))
    for row in all_rchms:
        if row:
            print(frm.format(*row))
    
    return scores

def plot_run_hist(runs, inds, shifts, max_overlap = 0.9):
    
    run_data = []
    
    for r,i,s in zip(runs, inds, shifts):
        
        if r < 1:
            continue
        
        overlap = (r - abs(s)) / r
        if overlap > max_overlap:
            continue
        
        run_data.append(r)
    
    rhist, bins = np.histogram(run_data, bins = max(run_data), range = (0, max(run_data)), density = True)
    ScalarPlot(rhist, add_range = True, minval = 0).show()

def test_err_run(seqa, seqb):
    
    mean_err_ratios = []
    mean_quals = []
    maxmean_quals = []
    
    last_mmq = -1
    
    for max_err in range(0, 20):
        
        # print(f"max err = {max_err}")
        
        runs, inds, shifts, errs = process.correlate_longest_subseq_err(seqa, seqb, fill = 0, max_err = max_err)
        
        quals = [r-e for r, e in zip(runs, errs)]
        err_ratios = [(r-e)/r if r>0 else 0 for r, e in zip(runs, errs)]
        
        mean_err_ratios.append(np.mean(err_ratios))
        mean_quals.append(np.mean(quals))
        mmq = max(quals) - np.mean(quals)
        maxmean_quals.append(mmq)
        
        # plot_run_hist(runs, inds, shifts)
        
        # if max_err in [10]:
        #     print(f"max_err: {max_err}, mmq: {mmq:0.3f}, last mmq: {last_mmq:0.3f}")
            
        #     buffer = 0
        #     topk = 1
        #     tops, topdatas = process.extract_max_runs(seqa, seqb, runs, inds, shifts, topk, buffer = buffer)
        #     draw.highlight_run(seqa, seqb, tops, topdatas, suppress = False)
            
        #     diffs = align.align_sequences(*tops[0])
        #     net_offset = sum([abs(a-b) for a,b in diffs])
        #     diffs2 = align.align_sequences(tops[0][1], tops[0][0])
        #     net_offset2 = sum([abs(a-b) for a,b in diffs])
            
        #     seqafill, seqbfill = align.fill_aligned(*tops[0], diffs)
        #     print("aligned (?)")
        #     print(seqafill)
        #     print(seqbfill)
        #     seqafill, seqbfill = align.fill_aligned(tops[0][1],tops[0][0], diffs2)
        #     print("aligned (?)")
        #     print(seqafill)
        #     print(seqbfill)
        #     print()
        
        last_mmq = mmq
        
    print("mean err ratios")
    ScalarPlot(mean_err_ratios, add_range = True, minval = 0).show()
    print("mean quals")
    ScalarPlot(mean_quals, add_range = True, minval = 0).show()
    print("max to mean quals")
    ScalarPlot(maxmean_quals, add_range = True, minval = 0).show()

def get_intron_run_metrics(intron_data, topk = 1, max_err = 16):
    
    buffer = 0
    
    seqa = intron_data.intron_seq
    seqb = seqa
    runs, inds, shifts, errs = process.correlate_longest_subseq_err(seqa, seqb, fill = 0, max_err = max_err)
    
    print("forward runs")
    ScalarPlot(runs, add_range = True, minval = 0).show_chunks()
    tops, topdatas = process.extract_max_runs(seqa, seqb, runs, inds, shifts, topk, buffer = buffer)
    draw.highlight_run(seqa, seqb, tops, topdatas, suppress = False)
    
    print("run hist")
    frhist, bins = np.histogram(runs, bins = max(runs), range = (0, max(runs)), density = True)
    ScalarPlot(frhist, add_range = True, minval = 0).show()
    
    cons = bio.merge(*tops[0])
    cons_entr = bio.consensus_entropy(cons)
    print(f"consensus: {cons}, entropy = {cons_entr:0.3f}")
    
    ##### reverse comp #####
    rc_seqb = bio.reverse_complement(seqb)
    rcruns, rcinds, rcshifts, rcerrs = process.correlate_longest_subseq_err(seqa, rc_seqb, max_err = max_err)
    
    print("rc runs")
    ScalarPlot(rcruns, add_range = True, minval = 0).show_chunks()
    rctops, rctopdatas = process.extract_max_runs(seqa, rc_seqb, rcruns, rcinds, rcshifts, topk, buffer = buffer)
    draw.highlight_run(seqa, rc_seqb, rctops, rctopdatas, suppress = False)
    
    print("rc run hist")
    rcrhist, bins = np.histogram(rcruns, bins = max(rcruns), range = (0, max(rcruns)), density = True)
    ScalarPlot(rcrhist, add_range = True, minval = 0).show()
    
    cons = bio.merge(*rctops[0])
    cons_entr = bio.consensus_entropy(cons)
    print(f"consensus: {cons}, entropy = {cons_entr:0.3f}")

def get_intron_corr_metrics(intron1_data, intron2_data, topk = 1, scale = 64):
    
    seqa = intron1_data.intron_seq
    seqb = intron2_data.intron_seq
    
    # map = {"A":"R","G":"R","C":"Y","T":"Y"}
    # comp_func = lambda a, b: map.get(a) == map.get(b)
    # score_func = lambda a, b, sc: 1.0
    comp_func = None
    score_func = None
    
    corr, rccorr = process.correlate(seqa, seqb, fill = 0, comparison_func = comp_func, score_func = score_func, scale = scale)
    
    num_bins = max(16, int(np.sqrt(len(seqa))))
    hist_max = 0.5
    
    cplot = ScalarPlot(corr, add_range = True, minval = 0)
    rcplot = ScalarPlot(rccorr, add_range = True, minval = 0)
    ScalarPlot.show_paired(cplot, rcplot, chunksz = 256)
    
    print("###### forward ######")
    chist, bins = np.histogram(corr, bins = num_bins, range = (0, hist_max), density = True)
    ScalarPlot(chist, add_range = True, xmin = bins[0], xmax = bins[-1], ruler = True, num_labels = 2, ticks = 0).show()
    
    mspec = {}
    corr_seqs = process.extract_top_correlated(seqa, seqb, corr, topk, scale=scale)
    for sa, sb, scorr in corr_seqs:
        mspec = heal.get_mutation_spectrum(sa, sb, mutation_spec = mspec)

        print(f"correlation {scorr:0.3f}")
        sah, sbh = draw.highlight_matching(sa, sb, suppress = True)
        print()
        
        print("aligned:")
        algn = align.align_sequences(sa, sb)[0]
        algn.print()
    
    print("mutation spectrum:")
    ScalarPlot(mspec, labels = True, bit_depth = 16, add_range = True, space = 1).show()
    print()
    
    print("###### reverse ######")
    rcmspec = {}
    rcchist, bins = np.histogram(rccorr, bins = num_bins, range = (0, hist_max), density = True)
    ScalarPlot(rcchist, add_range = True, xmin = bins[0], xmax = bins[-1], ruler = True, num_labels = 2, ticks = 0).show()
    
    rcseqb = bio.reverse_complement(seqb)
    rccorr_seqs = process.extract_top_correlated(seqa, rcseqb, rccorr, topk, scale=scale)
    for sa, sb, srccorr in rccorr_seqs:
        rcmspec = heal.get_mutation_spectrum(sa, sb, mutation_spec = rcmspec, b_is_rc = True)
        
        print(f"correlation {srccorr:0.3f}")
        sah, sbh = draw.highlight_matching(sa, sb)
        print()
    
    print("RC mutation spectrum:")
    ScalarPlot(rcmspec, labels = True, bit_depth = 16, add_range = True, space = 1).show()
    
    
def get_intron_alignment(intron1_data, intron2_data):
    
    algns= align.align_sequences(intron1_data.intron_seq, intron2_data.intron_seq)
    algn = algns[0]
    algn.print()

def find_intron(gm, min_len = 0, max_len = 0):
    
    while True:
        rintr = get_random_intron(gm)
        seq_len = len(rintr.intron_seq)
        if seq_len < min_len:
            continue
        elif max_len and seq_len > max_len:
            continue
        
        print(rintr)
        
        if len(rintr.intron_seq) < 253:
            print(rintr.intron_seq)
        else:
            print(rintr.intron_seq[:253], "...")
            print("..." + rintr.intron_seq[len(rintr.intron_seq) - 253:])
        
        res = input()
        if 'n' in res.lower():
            pass
        else:
            break
        
    return rintr

def measure_random_introns(gm, min_len = 0, max_len =0):
    
    keeps = []
    while True:
        
        intr1 = get_random_intron(gm, min_len, max_len)
        print(str(intr1))
        
        scores = compare_introns(intr1.seq, intr1.seq, chunksz = 64)
        
        res = input()
        if "y" in res.lower():
            keeps.append(intr1)
        elif 'n' in res.lower():
            break
    
    return keeps

def compare_random_introns(gm, min_len = 0, max_len = 0):
    
    intr1 = get_random_intron(gm, min_len, max_len)
    
    keeps = []
    while True:
        
        intr2 = get_random_intron(gm, min_len, max_len)
        
        print(f"comparing introns:")
        print(str(intr1))
        print(str(intr2))
    
        scores = compare_introns(intr1.seq, intr2.seq, chunksz = 64)
    
        # get_intron_corr_metrics(intr1, intr2)
        # get_intron_alignment(intr1, intr2)
        
        res = input()
        if "y" in res.lower():
            keeps.append((intr1, intr2))
        elif 'n' in res.lower():
            break
            
        intr1 = intr2
    
    return keeps

def compare_gene_introns(gm, nchr, gene_index=None, gene_name = "", max_num = 8):
    
    introns = get_gene_introns(gm, nchr, gene_index=gene_index, gene_name = gene_name, max_num = max_num)
    nintrons = len(introns)
    
    for n in range(nintrons):
        intr1 = introns[n]
        print(f"self similarity of intron {intr1}")
        scores = compare_introns(intr1.seq, intr1.seq, chunksz = 64)
        input()
        
        for nn in range(n+1, nintrons):
            intr2 = introns[nn]
            print("comparing introns:")
            print(str(intr1))
            print(str(intr2))
            scores = compare_introns(intr1.seq, intr2.seq, chunksz = 64)
            
            input()
    

def main():
    
    gm = load_genome()
    
    
    # nchr = 12
    # gene_index = 176
    # exon_index = 7
    
    # nchr = 14
    # gene_index = 288
    # exon_index = 9
    
    # nchr = 11
    # gene_index = 11
    # exon_index = 0
    
    # nchr = 'Y'
    # gene_index = 9
    # exon_index = 8
    
    # intr = get_intron(gm, nchr, gene_index, exon_index)
    # print("alignment scores:")
    # scores = compare_introns(intr.seq, intr.seq, chunksz = 64)
    # print("runs:")
    # scores = compare_introns(intr.seq, intr.seq, chunksz = 64, score_mode = "runs")
    
    # get_intron_run_metrics(intr)
    # get_intron_corr_metrics(intr, intr, scale = None)
    
    # keeps = measure_random_introns(gm, min_len = 512, max_len = 2048)
    # for a in keeps:
    #     print(str(a))
    
    # keeps = compare_random_introns(gm, min_len = 500, max_len = 2000)
    # for a, b in keeps:
    #     print("intron pair:")
    #     print(str(a))
    #     print(a.to_tuple())
    #     print(str(b))
    #     print(b.to_tuple())
    #     print()
    
    compare_gene_introns(gm, 20, max_num = 8, 
                        gene_index = 1200,
                        # gene_name = "KISS1"
                    )
    
    
    interesting_pairs = [
        
        # contain ACACATGTAT, ATACATGTGT motif that may form stem loop
        [("ART1",11, 5),
        ("SEZ6",17, 2)],
        
        # CAGCGGTG, plus tons of RC homology
        [("NXF2","X", 9),
        ("NALF1",13, 2)],
        
    ]
    
    
    """
    impressively similar:
    intron 9 on gene MAOA X:43654907-43746817 (gene index 177) with 1480 bases
    intron 9 on gene BCLAF3 X:19912860-19991061 (gene index 96) with 1480 bases
    
    big strip:
    intron 7 on gene RASD1 17:17494437-17496395 (gene index 297) with 843 bases
    intron 2 on gene NOP10 15:34339159-34343180 (gene index 73) with 654 bases
    
    its like fish lines!!:
    intron 101 at X:346775-346720 on gene CLDN34 (gene index 37) with 460 bases, 
    
    not much goin on?
    self similarity of intron intron 95 at 19:371322-371205 on gene STK11 (gene index 39) with 1310 bases, 
    
    little stars
    intron 109 at 19:408401-405445 on gene STK11 (gene index 39) with 606 bases, 
    
    bigg strip
    intron 10 at 20:157593-157454 on gene ? (gene index 39) with 1182 bases, 
    intron 35 at 20:278429-278330 on gene ? (gene index 39) with 744 bases, 
    compare:
    intron 32 at 20:277383-277226 on gene ? (gene index 39) with 948 bases, 
    intron 35 at 20:278429-278330 on gene ? (gene index 39) with 744 bases, 
    lump!
    self similarity of intron intron 35 at 20:278429-278330 on gene ? (gene index 39) with 744 bases, 

    scaly
    self similarity of intron intron 32 at 20:277383-277226 on gene ? (gene index 39) with 948 bases, 
    
    """

if __name__=="__main__":
    main()
