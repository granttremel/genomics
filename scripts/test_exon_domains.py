
import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

from tabulate import tabulate

import numpy as np

from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.translate import Ribosome
from ggene.seqs import bio, process, find

@dataclass
class ExonData:
    gene:Dict[str,Any]
    exon:Dict[str,Any]
    cds:Dict[str,Any]
    exon_seq:str
    cds_seq:str
    cds_aa_seq:str
    
    @property
    def gene_name(self):
        return self.gene.get("info",{}).get("gene_name","?")

    def __repr__(self):
        return f"ExonData({self.gene_name}, {self.cds_aa_seq})"

@dataclass
class SimilarityResult:
    exon_data:ExonData
    similarity:float
    cds_aa_seq:str

def load_genome():
    return GenomeManager()

def get_cds_exon(gm, cds):
    
    nchr = cds.get("chrom")
    
    cds_start = cds.get("start")
    cds_end = cds.get("end")
    cds_exon_num = cds.get("info",{}).get("exon_number")
    
    exon = None
    for ex_feat in gm.gene_map.fetch(nchr, cds.get("start"), cds.get("end"), features = ('exon',)):
        exon_start = ex_feat.get("start")
        exon_end = ex_feat.get("end")
        exon_num = ex_feat.get("info",{}).get("exon_number")
        if exon_start <= cds_start and exon_end >= cds_end and exon_num == cds_exon_num:
            exon = ex_feat
            break
    
    if exon is None:
        print(f"failed to identify exon for cds {cds}")
    
    exon_seq = gm.get_sequence(nchr, exon.get("start"), exon.get("end"))
    cds_seq = gm.get_sequence(nchr, cds.get("start"), cds.get("end"))
    
    exon_seq = bio.to_dna(exon_seq)
    cds_seq = bio.to_dna(cds_seq)
    
    if exon.get("strand") == "-":
        exon_seq = bio.reverse_complement(exon_seq)
        cds_seq = bio.reverse_complement(cds_seq)
    
    cds_aa_seq, res = translate_cds(cds_seq, int(cds.get("frame", 0)))
    
    return exon, exon_seq, cds_seq, cds_aa_seq, res

def translate_cds(cds_seq, frame = 0):
    
    cd_table = bio.get_codon_table(rna = False)
    
    aas = []
    i = 0
    for i in range(frame, len(cds_seq), 3):
        
        cdn = cds_seq[i:i+3]
        if len(cdn) < 3:
            break
        aa = cd_table.get(cdn, "X")
        if aa == "X":
            print(f"failed to translate codon {cdn} at position {i} of cds")
        aas.append(aa)
    
    resid = cds_seq[i:]
    if len(resid) == 3:
        resid = ""
    
    return "".join(aas), resid

def compare_cds_seqs_simple(aa_seq_a, aa_seq_b, topk = 5, err_tol = 2):
    
    # runs, inds, shifts = process.correlate_longest_subseq(aa_seq_a, aa_seq_b, err_tol = err_tol)
    
    aa_cls_rev = bio.AA_CLASSES_REV
    
    def score_func(a, b, sc):
        if a==b:
            return 1
        else:
            a_cls = aa_cls_rev.get(a)
            b_cls = aa_cls_rev.get(b)
            if len(set(a_cls).intersection(b_cls)) > 0:
                return 0.5
        return 0.0
    
    sims, _ = process.correlate(aa_seq_a, aa_seq_b, score_func = score_func)
    
    # score = sum([r for r in runs if r >= run_thresh]) / max(len(aa_seq_a),len(aa_seq_b))
    score = sum(sorted(sims, reverse = True)[:topk])
    return score
    

def compare_cds_seqs(aa_seq_a, aa_seq_b, run_thresh_f = 0.1, err_tol = 5):
    
    if len(aa_seq_a) > len(aa_seq_b):
        aa_seq_a, aa_seq_b = aa_seq_b, aa_seq_a
    
    len_a = len(aa_seq_a)
    len_b = len(aa_seq_b)
    rev_a = bio.reverse(aa_seq_a)
    rev_b = bio.reverse(aa_seq_b)
    
    min_seq_len = len(aa_seq_a)
    run_thresh = int(run_thresh_f * min_seq_len)
    
    runs, inds, shifts = process.correlate_longest_subseq(aa_seq_a, aa_seq_b)
    
    score = 0
    if max(runs) < run_thresh:
        print(f"early out due to low similarity: max run {max(runs)}")
        return score
    
    for rlen, ri, rshift in zip(runs, inds, shifts):
        
        if rlen < run_thresh:
            continue
        
        subseq_a_fwd = aa_seq_a[ri:]
        subseq_b_fwd = aa_seq_b[ri + rshift:]
        match, err = find.count_matching(subseq_a_fwd, subseq_b_fwd, err_tol = err_tol)
        score += match / len_a
        
        subseq_a_rev = rev_a[len_a - ri:]
        subseq_b_rev = rev_b[len_b - rshift - ri:]
        match, err = find.count_matching(subseq_a_rev, subseq_b_rev, err_tol = err_tol)
        score += match / len_a
        
    return score

def scan_cds_similarity(gm, ref_aa_seq, chr, max_genes = 10, max_cds = None, start = 0, min_len = 10):
    
    sim_data = []
    
    ngs = 0
    for gene in gm.gene_map.fetch(chr, start = start, end = None, features = ('gene',)):
        
        ncds = 0
        last_end = -1
        
        gene_start = gene.get("start")
        gene_end = gene.get("end")
        gene_name = gene.get("info",{}).get("gene_name")
        if not gene_name:
            continue
        
        for test_cds in gm.gene_map.fetch(chr, start = gene_start, end = gene_end, features = ('CDS',)):
            
            # exon, exon_seq, test_aa_seq, res = get_cds_exon(gm, test_cds)
            
            cds_start = test_cds.get("start")
            cds_end = test_cds.get("end")
            cds_ex_num = test_cds.get("info",{}).get("exon_number")
            if cds_ex_num is None:
                continue
            
            if last_end > cds_start:
                continue
                
            last_end = cds_end
                
            test_seq = gm.get_sequence(chr, cds_start, cds_end)
            
            if test_cds.get("strand") == '-':
                test_seq = bio.reverse_complement(test_seq)
            
            test_aa_seq, res = translate_cds(test_seq, int(test_cds.get("frame", 0)))
            if len(test_aa_seq) < min_len:
                continue
            
            test_exdat = ExonData(gene, {}, test_cds, "", test_seq, test_aa_seq)
            
            # score = compare_cds_seqs(ref_aa_seq, test_aa_seq)
            score = compare_cds_seqs_simple(ref_aa_seq, test_aa_seq, err_tol = 3)
            
            # if score > 0:
            #     score_func = lambda a, b, sc: int(a==b)
            #     corr, _ = process.correlate(ref_aa_seq, test_aa_seq, score_func = score_func)
            #     ScalarPlot(corr, add_range = True).show()
        
            sim_data.append(SimilarityResult(test_exdat, score, test_aa_seq))
            ncds += 1
            
            if max_cds and ncds >= max_cds:
                break
        
        ngs += 1
        
        if max_genes and ngs >= max_genes:
            break
    sim_data = list(sorted(sim_data, key = lambda s:-s.similarity))
    return sim_data

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

def get_cds_by_index(gm, nchr, gene, cds_ind):
    
    start = gene.get("start", 0)
    end = gene.get("end")
    
    cds, ncd = get_feature_by_index(gm, nchr, cds_ind, start=start, end=end, feature_type = 'CDS')
    if cds is None:
        cds_ind = cds_ind % ncd
        cds, ncd = get_feature_by_index(gm, nchr, cds_ind, start=start, end=end, feature_type = 'CDS')
    
    return cds, ncd

def get_exon(gm:GenomeManager, nchr, gene_index, cds_index):
    
    gene, ng = get_gene_by_index(gm, nchr, gene_index)
    cds, ncd = get_cds_by_index(gm, nchr, gene, cds_index)
    exon, exon_seq, cds_seq, cds_aa_seq, res = get_cds_exon(gm, cds)
    exdat = ExonData(gene, exon, cds, exon_seq, cds_seq, cds_aa_seq)
    
    print(f"selected cds {ncd} from exon {exon.get('info',{}).get('exon_number', "?")} on gene {gene.get('info',{}).get('gene_name','?')} (gene index {gene_index}) with {len(exon_seq)} bases and {len(cds_aa_seq)} amino acids")
    
    return exdat

def get_random_exon(gm:GenomeManager):
    
    nchr = random.choice(list(gm.gene_map.max_indices.keys()))
    
    num_gene = random.randint(0, 500)
    
    gene, ng = get_gene_by_index(gm, nchr, num_gene)
    
    ncd = random.randint(0, 10)
    cds, ncd = get_cds_by_index(gm, nchr, gene, ncd)
    
    exon, exon_seq, cds_seq, cds_aa_seq, res = get_cds_exon(gm, cds)
    
    exdat = ExonData(gene, exon, cds, exon_seq, cds_seq, cds_aa_seq)
    
    print(f"selected cds {ncd} from exon {exon.get('info',{}).get('exon_number', "?")} on gene {gene.get('info',{}).get('gene_name','?')} (gene index {num_gene}) with {len(exon_seq)} bases and {len(cds_aa_seq)} amino acids")
    
    return exdat
    
def aa_seq_montecarlo(num_tests = 10, len_seqs = 128, amino_acid = True):
    
    if amino_acid:
        ref_seq = bio.get_random_aa_sequence(len_seqs)
    else:
        ref_seq = bio.get_random_sequence(len_seqs)
    
    res = []
    
    for nt in range(num_tests):
        if amino_acid:
            test_seq = bio.get_random_aa_sequence(len_seqs)
        else:
            test_seq = bio.get_random_sequence(len_seqs)
        
        data, inds, shifts = process.correlate_longest_subseq(ref_seq, test_seq)
        # data, rcdata = process.correlate_v2(ref_seq, test_seq)
        
        max_run = max(data)
        res.append(max_run)
    
    mean = np.mean(res)
    sd = np.std(res)
    vmin = min(res)
    vmax = max(res)
    
    print(f"monte carlo results with seq len {len_seqs} over {num_tests} tests:")
    print(f"mean max run: {mean:0.3f} ({mean/len_seqs:0.3f} per BP)")
    print(f"sd max run: {sd:0.3f} ({sd/np.sqrt(len_seqs):0.3f} per root BP)")
    print(f"min max run: {vmin}")
    print(f"max max run: {vmax}")
    
    sp = ScalarPlot(data, add_range = True)
    sp.show()
    
    print(repr(sp), repr(sp.ruler))
    
    return res, mean, sd

def run_aa_seq_montecarlo():
    means = []
    sds = []
    lens = [32, 64, 128, 256, 512]
    nt_fact = 1000
    # try_cds_similarity()
    for _len in lens:
        # num_tests = max(250, 250*256 // _len)
        num_tests = nt_fact*256 // _len
        _, mean, sd = aa_seq_montecarlo(num_tests = num_tests, len_seqs = _len, amino_acid = False)
        means.append(mean)
        sds.append(sd)
    
    from matplotlib import pyplot as plt
    
    f, ax = plt.subplots()
    ax.scatter(lens, means)
    ax.set_title("mean max run length")
    f.savefig("./data/outputs/mean_max_run.png")
    
    f, ax = plt.subplots()
    ax.scatter(lens, sds)
    ax.set_title("sd max run length")
    f.savefig("./data/outputs/sd_max_run.png")

def try_cds_similarity():
    gm = load_genome()
    
    # exdat = get_random_exon(gm)
    exdat = get_exon(gm, "2", 231, 3)
    
    # print(exdat.exon_seq)
    print(exdat.cds_aa_seq)
    
    test_aa_seq1 = exdat.cds_aa_seq
    test_aa_seq2 = exdat.cds_aa_seq[:20] + "A" + exdat.cds_aa_seq[20:]
    test_aa_seq3 = exdat.cds_aa_seq[:20] + "AAAAAAAAAAAA" + exdat.cds_aa_seq[20:120] + "BBBBBBBBBBB" + exdat.cds_aa_seq[120:]
    test_aa_seq4 = bio.get_random_aa_sequence(len(exdat.cds_aa_seq))
    
    test_seqs = [test_aa_seq1, test_aa_seq2, test_aa_seq3, test_aa_seq4]
    
    for i in range(len(test_seqs)):
        ts = test_seqs[i]
        runs, inds, shifts = process.correlate_longest_subseq(exdat.cds_aa_seq, ts)
        sc1 = compare_cds_seqs_simple(exdat.cds_aa_seq, test_aa_seq1)
        print(f"similarity to test seq {i}:")
        print(f"score: {sc1}")
        print(f"runs: {runs}")
        # print(f"shifts: {shifts}")
        print()
    
    # sims = scan_cds_similarity(gm, exdat.cds_aa_seq, "2", max_genes=10, max_cds = 3, start = 10e6)
    # for sim in sims:
    #     print(sim)
    
    
    pass

def display_cds_similarity(ref_seq, test_seq):
    
    runs, inds, shifts = process.correlate_longest_subseq(ref_seq, test_seq, err_tol = 3)
    
    ScalarPlot(runs, minval = 0, add_range = True).show()
    
    maxrun = max(runs)
    argmax = runs.index(maxrun)
    maxind = inds[argmax]
    maxshift = shifts[maxind]
    
    ref_pre = " "*abs(maxshift)
    test_pre = ""
    if maxshift > 0:
        test_pre, ref_pre = ref_pre, test_pre
    
    center_diff = len(ref_seq)//2 - len(test_seq) // 2
    
    match_start_ref = maxind
    match_end_ref = maxind + maxrun
    
    match_start_test = maxind - maxshift
    match_end_test = maxind + maxrun - maxshift
    
    print(ref_pre, ref_seq[:match_start_ref], ref_seq[match_start_ref:match_end_ref], ref_seq[match_end_ref:])
    print(test_pre, test_seq[:match_start_test], test_seq[match_start_test:match_end_test], test_seq[match_end_test:])

def test_corr_subseq():
    
    seq1 = "VQKTQPRSALKLYQDLSLLHANELLLNRGWFCHLRNDSHYVVYTRELDGIDRIFIVVLNFGESTLLNLHNMISGLPAKMRIRLSTNSADKGSKVDTSGIFLDKGEGLIFEHNTKNLLHRQTAFRDRCFVSNRACYSSVLNILYTSC"
    seq2 = "                             VRVWRYLKGKDLVARESLLDGGNKVVISGFGDPLICDNQVSTGDTRIFFVNPAPPYLWPAHKNELMLNSSLMRITLRNLEEVEFCVE                              "
    # seq1 = "VQKTQPRSALKLYQDLSL LHANELLLN RGWFCHLRNDSHYVVYTRELDGIDRIFIVVLNFGESTLLNLHNMISGLPAKMRIRLSTNSADKGSKVDTSGIFLDKGEGLIFEHNTKNLLHRQTAFRDRCFVSNRACYSSVLNILYTSC"
    # seq2 = "VRVWRYLKGKDLVARESLLDGGNKVVISGFGDPLICDNQVSTGDTRIFFVNPAPPYLWP AHKNELMLN SSLMRITLRNLEEVEFCVE"
    
    
    match1 = "LHANELLLN"
    match2 = "AHKNELMLN"
    
    m1ind = seq1.index(match1)
    m1ind_c = m1ind - len(seq1)//2
    
    m2ind = seq2.index(match2)
    m2ind_c = m2ind - len(seq2)//2
    
    print(m1ind, m1ind_c)
    print(m2ind, m2ind_c)
    
    print(seq2.count(" "))
    
    display_cds_similarity(seq1, seq2)
    
    """
    VQKTQPRSALKLYQDLSL LHANELLLN RGWFCHLRNDSHYVVYTRE
    TGDTRIFFVNPAPPYLWP AHKNELMLN SSLMRITLRNLEEVEFCVE
    
    .                             .
    .                              .
    
    """
    
    # runs, inds, shifts = process.correlate_longest_subseq(seq1, seq2, err_tol = 3)
    # ScalarPlot(runs, add_range = True).show()
    
    
    
    pass

def get_similarity_matrix(aa_seqs):
    
    num_seqs = len(aa_seqs)
    arr = np.zeros((num_seqs, num_seqs))
    
    for i in range(len(aa_seqs)):
        for j in range(i, num_seqs):
            sim = compare_cds_seqs_simple(aa_seqs[i], aa_seqs[j])
            arr[i, j] = sim
    
    return arr

def main():
    
    ref_chr = "1"
    # gene_index = 1740
    gene_index = 1739
    cds_index = 0
    
    gm = load_genome()
    exdat = get_exon(gm, ref_chr, gene_index, cds_index)
    exdat2 = get_exon(gm, ref_chr, gene_index, cds_index+1)
    
    full_exon_seq = exdat2.exon_seq + exdat.exon_seq
    
    # exdat = get_random_exon(gm)
    
    # i = 0
    # for g in gm.gene_map.fetch(ref_chr, 0, end=None, features = ("gene",)):
    #     # gi,_ = get_gene_by_index(gm, ref_chr, i)
        
    #     gn = g.get("info",{}).get("gene_name")
    #     i += 1
    #     print(i,gn)
        
    #     if gn == "KISS1":
    #         break
    #     pass
    
    # return
    
    print(f"cds index {cds_index}")
    print(exdat.gene)
    print(exdat.exon)
    print(exdat.cds)
    
    full_exon_seq = exdat2.exon_seq + exdat.exon_seq
    full_cds_seq = exdat2.cds_seq + exdat.cds_seq
    full_cds_aa, _ = translate_cds(full_cds_seq, int(exdat2.cds.get("info",{}).get("frame", 0)))
    
    print("exon:", full_exon_seq)
    print("cds:", full_cds_seq)
    print("cds:", full_cds_aa)
    
    cds_index += 1
    
    test_chr = "2"
    
    res = scan_cds_similarity(gm, exdat.cds_aa_seq, test_chr, max_genes = None)
    sims = [r.similarity for r in res]
    mean = np.mean(sims)
    sd = np.std(sims)
    maxsim = max(sims)
    num_zeros = sum([1 for s in sims if s ==0])
    pcts = np.percentile(sims, [1, 99.5])
    hist, bins = np.histogram(sims, bins = int(np.sqrt(len(sims))), range = (pcts[0], max(sims)))
    
    print(f"reference sequence: {exdat.cds_aa_seq}")
    print(f"similarity mean = {mean:0.3f}, sd = {sd:0.3f}, max = {maxsim:0.3f}, number zero = {num_zeros}, n = {len(res)}")
    ScalarPlot(hist, add_range = True, minval = 0, ruler = True, xmin = bins[0], xmax = bins[-1], ticks = 0, minor_ticks = 0).show()
    print()
    for r in res[:10]:
        if r.similarity == 0:
            pass
        else:
            print(r)
    
    print()
    
    gene_names = [exdat.gene_name] + [r.exon_data.gene_name for r in res[:20]]
    sim_mat = get_similarity_matrix([exdat.cds_aa_seq] + [r.exon_data.cds_aa_seq for r in res[:20]])
    print(tabulate(sim_mat, headers = gene_names, floatfmt = "0.2f"))
    
    display_cds_similarity(exdat.cds_aa_seq, res[0].exon_data.cds_aa_seq)



if __name__=="__main__":
    main()