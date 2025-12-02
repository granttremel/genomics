

import random
from typing import List, Dict, Union, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab, compare
from ggene.seqs.bio import reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE

from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.motifs.motif import MotifDetector

import random
from ggene.genomemanager import GenomeManager

@dataclass
class GeneResult:
    gene:Dict[str,Any]
    introns: List[Tuple[Dict[str, Any]]]
    intron_seqs: List[str]
    cds: List[Dict[str,Any]]
    cds_seqs: List[str]
    aa_seq: str

def load_genome():
    return GenomeManager()

def get_gene_data(gm:GenomeManager, chr, gene_name, start, end = None):
    
    for f in gm.gene_map.fetch(chr, start, end=end, features = ('gene',)):
        
        gn = f.get("info",{}).get("gene_name","")
        if gn == gene_name:
            break
    
    gene = f
    
    cds = []
    cds_seqs = []
    
    last_cds = None
    for f in gm.gene_map.fetch(chr, gene.get("start"), end = gene.get("end"), features = ('CDS',)):
        
        if last_cds:
            if last_cds.get("end") > f.get("start"):
                continue
        
        cds.append(f)
        cds_seq = gm.get_sequence(chr, f.get("start"), f.get("end"), frame = int(f.get("frame", 0)))
        cds_seqs.append(cds_seq)
        
        last_cds = f
    
    full_cds = "".join(cds_seqs)
    aa_seq, r = translate_cds_seq(full_cds)
    print(f"remaining: {r}")
    
    introns = []
    intron_seqs = []
    
    last_exon = None
    for f in gm.gene_map.fetch(chr, gene.get("start"), end = gene.get("end"), features = ('exon',)):
        
        if not last_exon:
            last_exon = f
            continue
        
        if last_exon.get("end") > f.get("start"):
            continue
            
        introns.append((last_exon, f))
        
        intron_seq = gm.get_sequence(chr, last_exon.get("end"), f.get("start"))
        intron_seqs.append(intron_seq)
        
        last_exon = f
        
    res = GeneResult(gene, introns, intron_seqs, cds, cds_seqs, aa_seq)
    return res
        
def translate_cds_seq(cds_seq):
    
    aa_seq = []
    i=0
    for i in range(0, len(cds_seq), 3):
        cdn = cds_seq[i:i+3]
        aa = CODON_TABLE.get(cdn,"X")
        aa_seq.append(aa)
    
    rem = ""
    if len(cds_seq) - i < 3:
        rem = cds_seq[-i-1:]
    
    return "".join(aa_seq), rem

def get_random_seq(seq_len = 2048):
    return [vocab.VOCAB[random.randint(0, 3)] for i in range(seq_len)]    

def compare_random_seqs(seq_len = 2048, chunksz = 64, score_modes = ["runs"], scale = None, shift = 1):
    rand_seq = get_random_seq(seq_len = seq_len)
    compare.compare_sequences(rand_seq, rand_seq, chunksz=chunksz, score_modes = score_modes, resample = True, scale = scale, shift_step = shift)

def get_comparison_stats(num_samples, seq_len, auto = False, cycle = False):
    
    nruns = []
    nrcruns =[]
    
    rs1 = get_random_seq(seq_len)
    for n in range(num_samples):
        
        rs2 = get_random_seq(seq_len)
        
        if auto:
            rs1 = rs2
        elif cycle:
            pass
        else:
            rs1 = get_random_seq(seq_len)
            
        
        # data, rcdata = compare.score_sequences_runs(rs1, rs2)
        data, rcdata = compare.score_sequences_corrs(rs1, rs2)
        
        
        nruns.append(data)
        nrcruns.append(rcdata)
        
        rs2 = rs1
    
    mean = np.mean(nruns)
    sd = np.std(nruns)
    rcmean = np.mean(nrcruns)
    rcsd = np.std(nrcruns)
    
    print(f"runs over {num_samples} samples, {seq_len} sequence length: mean={mean:0.3f}, sd={sd:0.3f}, min={min(nruns):0.3f}, max={max(nruns):0.3f}")
    print(f"runs over {num_samples} samples, {seq_len} sequence length: mean={rcmean:0.3f}, sd={rcsd:0.3f}, min={min(nrcruns):0.3f}, max={max(nrcruns):0.3f}")

    return mean, sd, rcmean, rcsd

def analyze_gene(gm, chr, gene_name, start = 0, gene_chunks = 128, score_modes = ["runs"], scale = None, shift = 1):
    
    generes = get_gene_data(gm, chr, gene_name, start=start - 1000)
    
    gene_start = generes.gene.get("start")
    gene_end = generes.gene.get("end")
    
    intron_lens = [len(intr) for intr in generes.intron_seqs]
    
    print(f"gene {gene_name} has {len(generes.intron_seqs)} introns (mean {np.mean(intron_lens):0.1f}, {min(intron_lens)}-{max(intron_lens)}), {len(generes.cds_seqs)} cds sequences encoding protein with {len(generes.aa_seq)} residues")
    
    print(f"final intron spans {generes.introns[-1][0].get("end")}-{generes.introns[-1][0].get("start")}")
    
    max_disp = 256
    num_chunks = 3*max_disp
    chunksz = (gene_end - gene_start)//gene_chunks
    start = gene_start - chunksz*(num_chunks - gene_chunks)//2
    length = num_chunks * chunksz
    
    gm.display_chromosomal_quantity(chr, "cds_len", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    gm.display_chromosomal_quantities(chr, ["cg","cpg"], chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    
    gene_seq = gm.get_sequence(chr, gene_start, gene_end)
    
    intron_start = generes.introns[-1][0].get("end") - start
    intron_end = generes.introns[-1][1].get("start") - start
    
    compare.compare_sequences(generes.intron_seqs[-1], generes.intron_seqs[-1], chunksz=chunksz, score_modes = score_modes, resample = True, scale = scale, shift_step = shift,
                              add_ruler = True, xmin = 0, xmax = 62, num_labels = 31)
    pass

def find_long_vars(gm, testchr, topk = 10, context_len = 256):
    
    vars = []
    
    for v in gm.annotations.stream_by_types(feature_types =['variant'], chrom = testchr, start = 0):
        
        vlen = len(v.attributes.get("alt")[0]) - len(v.attributes.get("ref"))
        
        if abs(vlen) < 10:
            continue
        
        in_gene = False
        gene_name = ""
        in_exon = False
        in_cds = False
        
        for ff in gm.annotations.stream_by_types(feature_types = ["gene","exon","cds"], chrom = testchr, start = v.start, end = v.start+1):
            
            if ff.feature_type == "gene":
                in_gene = True
                gene_name = ff.name if ff.name else ff.id
            elif ff.feature_type == "exon":
                in_exon = True
            elif ff.feature_type == "cds":
                in_cds = True
            
        v.attributes["in_gene"] = in_gene
        v.attributes["gene_name"] = gene_name
        v.attributes["in_exon"] = in_exon
        v.attributes["in_cds"] = in_cds
        
        upstream = gm.get_sequence(testchr, v.start - context_len, v.start)
        downstream = gm.get_sequence(testchr, v.end, v.end + context_len)
        v.attributes["upstream"] = upstream
        v.attributes["downstream"] = downstream
        
        vars.append(v)
    
    print(f"found {len(vars)} variants")
    
    vars = sorted(vars, key= lambda k:len(k.attributes.get("alt")[0]) - len(k.attributes.get("ref")))
    ins = vars[:topk]
    dels = vars[len(vars)-topk-1:]
    return ins, dels

def print_vars(gm):
    
    chrs = range(23)
    
    all_ins = []
    all_dels = []
    
    for chr in chrs:
        ins, dels = find_long_vars(gm, chr)
        print(f"variants on chromosome {chr}")
        all_ins.extend(ins)
        all_dels.extend(dels)
    
    all_ins = list(sorted(all_ins, key = lambda k: len(k.attributes.get("alt")[0]) - len(k.attributes.get("ref"))))[:10]
    all_dels = list(sorted(all_dels, key = lambda k: len(k.attributes.get("alt")[0]) - len(k.attributes.get("ref"))))[-10:]
    
    for v in all_ins + all_dels:
        vlen = len(v.attributes.get("alt")[0]) - len(v.attributes.get("ref"))
        vtype = "ins" if vlen > 0 else "del"
        
        host_feat_str = []            
        if v.attributes.get("in_gene"):
            host_feat_str.append(f"in gene {v.attributes.get("gene_name")}")
        if v.attributes.get("in_exon"):
            host_feat_str.append("in exon")
        if v.attributes.get("in_cds"):
            host_feat_str.append("in cds")
        # print()
        print(f"Variant ({vtype}) at {v.chrom}:{v.start}-{v.end} with length {abs(vlen)}, qual {v.attributes.get("qual"):0.3f}, zygosity {v.attributes.get("zygosity")}, {", ".join(host_feat_str)}")
            
        print(f"ref: {v.attributes.get("ref")}")
        print(f"alt: {v.attributes.get("alt")[0]}")
        print(f"upstream: {v.attributes.get("upstream")}")
        print(f"downstream: {v.attributes.get("downstream")}")
        
        print()

    



def get_random_stats():
    num_samples = 1024
    seq_len = 512
    
    allo_stats = []
    rcallo_stats = []
    auto_stats = []
    rcauto_stats = []
    
    num_seq_samps = 12
    
    metric = "corrs"
    min_seq_len = 8
    max_seq_len = 64
    yticks = np.linspace(0.01, 0.17, 8)
    sdyticks = np.linspace(0.01, 0.17, 8)
    
    # metric = "runs"
    # min_seq_len = 64
    # max_seq_len = 128
    # yticks = np.linspace(24, 32, 8)
    # sdyticks = np.linspace(0.95, 1.5, 8)
    
    
    odd = False
    seq_lens = 2*np.exp(np.linspace(np.log(min_seq_len/2), np.log(max_seq_len/2), num_seq_samps)).astype(int) + int(odd)
    
    for seq_len in seq_lens:
        
        m, sd, rcm, rcsd = get_comparison_stats(num_samples, seq_len = int(seq_len))
        allo_stats.append((m, sd))
        rcallo_stats.append((rcm, rcsd))
        
        m, sd, rcm, rcsd = get_comparison_stats(num_samples, seq_len = int(seq_len), auto = True)
        auto_stats.append((m, sd))
        rcauto_stats.append((rcm, rcsd))
    
    from matplotlib import pyplot as plt
    
    log_seq_lens = np.log(seq_lens)
    
    
    f, ax = plt.subplots()
    ax.plot(log_seq_lens, [np.log(m) for m, sd in allo_stats])
    ax.plot(log_seq_lens, [np.log(m) for m, sd in rcallo_stats])
    ax.plot(log_seq_lens, [np.log(m) for m, sd in auto_stats])
    ax.plot(log_seq_lens, [np.log(m) for m, sd in rcauto_stats])
    ax.set_xticks(log_seq_lens, labels = [format(int(sl), "") for sl in seq_lens])
    ax.set_yticks(np.log(yticks), labels = [format(yt, "0.2f") for yt in yticks])
    ax.legend(["allo", "rc allo", "auto", "rc auto"])
    ax.set_xlabel("sequence length")
    ax.set_ylabel("log mean autocorrelation response")
    f.savefig(f"./data/random_seq_{metric}_means.png")
    
    f, ax = plt.subplots()
    ax.plot(log_seq_lens, [np.log(sd) for m, sd in allo_stats])
    ax.plot(log_seq_lens, [np.log(sd) for m, sd in rcallo_stats])
    ax.plot(log_seq_lens, [np.log(sd) for m, sd in auto_stats])
    ax.plot(log_seq_lens, [np.log(sd) for m, sd in rcauto_stats])
    ax.set_xticks(log_seq_lens, labels = [format(int(sl), "") for sl in seq_lens])
    ax.set_yticks(np.log(sdyticks), labels = [format(yt, "0.2f") for yt in sdyticks])
    ax.legend(["allo", "rc allo", "auto", "rc auto"])
    ax.set_xlabel("sequence length")
    ax.set_ylabel("log std autocorrelation response")
    
    f.savefig(f"./data/random_seq_{metric}_sds.png")
    
    pass

def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    gm = load_genome()
    
    test_chr = 2
    gene_start = 184598470
    gene_end = 184939492
    gene_chunks = 128
    gene_name = "ZNF804A"
    
    num_steps = 16
    # chunksz = 4*len(gene_seq) // max_disp
    chunksz = 64
    scale = None
    shift = 1
    # score_modes = ["corrs"]
    score_modes = ["runs"]
    
    # analyze_gene(gm, test_chr, gene_name, start = gene_start)
    
    # compare_random_seqs(seq_len = 2047, chunksz = 128)
    
    # chrs = [7]
    
    from ggene.unified_stream import RepeatMaskerStream, BEDStream
    from ggene.sequence_stream import FASTAStream
    
    # rpts = BEDStream("./data/repeatmasker/repeats.sorted.bed")
    rpts = BEDStream("./data/repeatmasker/repeats.sorted.bed.gz")
    # rpts = BEDStream("./data/repeatmasker/repeats.bed")
    print(rpts)
    print(rpts.tabix)
    print(rpts.tabix.contigs)
    
    for rpt in rpts.tabix.fetch("chr2", start = 25160800):
        print(rpt)
        input()
    
    # for rpt in rpts.stream("1", 204.1e6):
    #     print(rpt)
    #     input()
    
    # gm.annotations.add_repeatmasker("./data/repeatmasker/hg38.fa.out.gz")
    # rpt_annos = gm.annotations.streams.get("repeats")
    
    # for rpt in rpt_annos.stream(2, 1e6):
    #     print(rpt)
    #     input()
    
    # print(gm.annotations.streams)
    

    

if __name__=="__main__":
    main()


