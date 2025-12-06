
from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np
import random

from Bio import Align

from ggene import DATA_DIR
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align
from ggene.seqs.vocab import VOCAB

@dataclass
class Repeat:
    repeat:Dict[str, Any]
    motif:str
    chr:str
    start:int
    end:int
    phase:int
    # context:int
    is_rc:bool
    seq:str = ""
    parent_gene:str = ""
    
    @property
    def seq_len(self):
        return len(self.seq)
    
    @property
    def motif_len(self):
        return len(self.motif)
    
    def __eq__(self, other):
        
        if len(self.motif) != len(other.motif):
            return False
            
        ext = self.motif * 2
        if other.motif in ext:
            return True
        elif bio.reverse_complement(other.motif) in ext:
            return True
        return False
    
    def combine(self, other):
        return Repeat(self.motif, self.instances+other.instances, self.count+other.count, self.length+other.length)
        
    def __hash__(self):
        return hash(self.motif)

def load_genome():
    gm = GenomeManager()
    rpts = unified_stream.BEDStream(DATA_DIR / "repeatmasker/repeats.sorted.bed.gz", feature_type = "repeats")
    gm.annotations.add_source("repeats", rpts)
    return gm, rpts

def extract_repeats(gm, chr, repeat_name = "", repeat_type = "", motif = "", start = None, end = None):
    
    rpts = []
    for rpt in gm.annotations.streams.get("repeats").stream(chr, start=start, end=end):
        
        if repeat_name and rpt.name != repeat_name:
            continue
        if repeat_type and rpt.attributes.get("type") != repeat_type:
            continue
        if motif and rpt.attributes.get("motif") != motif:
            continue
        
        seq = gm.annotations.get_sequence(chr, rpt.start, rpt.end)
        rpt.attributes["seq"] = seq
        rpts.append(rpt)
    
    return rpts

def get_motif_stats(gm, chr, min_length = 3, combine_rc = True, topk = 5, rank_by_len = False, rank_by_insts = False):
    
    hist_bins = 256
    motif_count = {}
    motif_len = {}
    motif_insts = {}
    rtn_groups = {}
    
    rpts = gm.annotations.streams.get("repeats")
    for rpt in rpts.stream(chr, start=0):
        
        motif = rpt.attributes.get("motif")
        if not motif:
            continue
        
        mtflen = len(motif)
        rptlen = rpt.end - rpt.start
        if mtflen < min_length:
            continue
        
        mkey, ind, _ = bio.get_least_seq_index(motif, do_rc = combine_rc)
        
        if not mkey in motif_insts:
            motif_count[mkey] = 0
            motif_len[mkey] = 0
            motif_insts[mkey] = []
            rtn_groups[mkey] = set()
        
        motif_count[mkey] += rptlen / len(motif)
        motif_len[mkey] += rptlen
        motif_insts[mkey].append(rpt.start)
        rtn_groups[mkey].add(motif)
    
    print(f"found {len(motif_count)} motifs on chr{chr}")
    if rank_by_len:
        top = sorted(list(motif_len.keys()), key = lambda k:-motif_len[k])
        print("top by length:")
    elif rank_by_insts:
        top = sorted(list(motif_insts.keys()), key = lambda k:-len(motif_insts[k]))
        print("top by instances:")
    else:
        top = sorted(list(motif_count.keys()), key = lambda k:-motif_count[k])
        print("top by count:")
    
    for i in range(topk):
        mtf = top[i]
        print(f"{mtf}: count = {motif_count[mtf]:0.1f}, len = {motif_len[mtf]:0.1f}, instances = {len(motif_insts[mtf])}")
        print(f"rotation group: {",".join(rtn_groups[mtf])}")
        disthist, bins = np.histogram(motif_insts[mtf], bins = hist_bins)
        ScalarPlot(disthist, add_range = True, xmin = bins[0], xmax = bins[-1], ruler = True, num_labels = 11, ticks = 0, genomic = True).show()
        show_motif_symms(mtf)
        print()
    
    return motif_count, motif_len, motif_insts

def pool_motif_stats(mtf_cts_a, mtf_len_a, mtf_insts_a, mtf_cts_b, mtf_len_b, mtf_insts_b):
    return {mtf:mtf_cts_a.get(mtf,0) + mtf_cts_b.get(mtf,0) for mtf in mtf_cts_a}, {mtf:mtf_len_a.get(mtf,0) + mtf_len_b.get(mtf,0) for mtf in mtf_len_a},{mtf:mtf_insts_a.get(mtf,[]) + mtf_insts_b.get(mtf,[]) for mtf in mtf_insts_a}

def get_motif_asymm(mtf):
    
    wgts = [n if b in 'GC' else 0 for n, b in enumerate(mtf)]
    cent = sum(wgts) / len(mtf)
    var = sum([(n-cent)**2 for n in wgts]) / len(mtf)*len(mtf)
    asymm = np.sqrt(var)/cent if cent else 0
    return cent, asymm

def show_motif_symms(mtf):
    
    asymms = []
    rcasymms = []
    for i in range(len(mtf)):
        rmtf = mtf[i:] + mtf[:i]
        rcrmtf = bio.reverse_complement(rmtf)
        m, sd = get_motif_asymm(rmtf)
        asymms.append(sd)
        rcm, rcsd = get_motif_asymm(rcrmtf)
        rcasymms.append(rcsd)
        print(i, rmtf, rcrmtf)
        
    print("asymmetry")
    ScalarPlot(asymms, add_range = True).show()
    print("RC asymmetry")
    ScalarPlot(rcasymms, add_range = True).show()

# def place_motif(mtf_tree:Dict[str,List[Any]], mtf, depth = 0):
def place_motif(mtf_tree:List[Dict[str, Any]], mtf, depth = 0):
    max_depth = -1
    max_mtf = None
    for subtree in mtf_tree:
        for submotif in subtree:
            if mtf in 2*submotif:
                new_tree, new_depth = place_motif(subtree[submotif].copy(), mtf, depth = depth + 1)
                
                if new_depth > max_depth:
                    max_depth = new_depth
                    max_mtf = submotif
        
    if max_mtf and max_depth > depth:
        
        
        pass
    else:
        mtf_tree.append({mtf:[]})
        max_depth = depth
        
    return mtf_tree, max_depth

def get_motif_hamdist(mtfa, mtfb):
    
    
    
    
    pass

def view_repeats(gm, chr, start, motifs_only = False):
    rpts = gm.annotations.streams.get("repeats")
    for rpt in rpts.stream(chr, start=start):
        if motifs_only and not rpt.attributes.get("motif"):
            continue
        print(rpt)
        input()

optimal_ordering_labels = [
    "None", # 0000
    "Fwd", # 0001
    "RC",  # 0010
    "Fwd+RC", # 0011
    "Comp", # 0100
    "Fwd+Comp", # 0101
    "Comp+RC", # 0110
    "Fwd+RC+Comp", # 0111
    "Rev", # 1000
    "Fwd+Rev", # 0001
    "RC+Rev",  # 0010
    "Fwd+RC+Rev", # 0011
    "Comp+Rev", # 0100
    "Fwd+Comp+Rev", # 0101
    "Comp+RC+Rev", # 0110
    "Fwd+RC+Comp+Rev", # 0111
]

def test_rcind_stats(seq_len, num_trials, random_seqs = True):
    
    stats = [0 for n in range(2**4)]
    log_opt_inds = []
    dyad_inds = set()
    ndyads = 0
    
    for n in range(num_trials):
        
        if random_seqs:
            test_seq = "".join([VOCAB[random.randint(0, 3)] for n in range(seq_len)])
        else:
            test_seq = bio.index_to_seq(n, seq_len = seq_len)
        
        opt_seq, opt_ind, v = bio.get_least_seq_index(test_seq, do_rc = True, do_rev = False, do_comp = False)
        
        l2v = np.log2(v)
        if l2v != int(l2v):
            lbl = optimal_ordering_labels[v]
            ndyads += 1
            dyad_inds.add(opt_ind)
            print(test_seq, opt_seq, lbl)
        
        # if opt_ind > 1e16:
        #     opt_ind = float(opt_ind)
        
        # log_opt = np.log(float(opt_ind)) if opt_ind > 0 else -1
        log_opt = opt_ind
        log_opt_inds.append(log_opt)
        stats[v] += 1
    
    hist, bins = np.histogram(log_opt_inds, bins = int(np.sqrt(num_trials)))
    
    ScalarPlot(hist, add_range = True).show()
    print(f"{bins[0]:0.0e} - {bins[-1]:0.0e}")
    print()
    
    return stats, ndyads, len(dyad_inds)

def gather_repeat_seqs(gm:GenomeManager, rpt_proto, chrs = [], context = 128, length_range = 0):
    
    rpt_proto, ind, arg_min = bio.get_least_seq_index(rpt_proto)
    rpt_test = 2*rpt_proto
    rpt_len = len(rpt_proto)
    
    rpts = gm.annotations.streams.get("repeats")
    
    out_rpts = []
    
    if not chrs:
        chrs = [str(i) for i in range(24)] + ["X","Y"]
    
    for chr in chrs:
        testchr = str(chr)
        
        for rpt in rpts.stream(testchr, start = None):
            new_seq = ""
            mtf = rpt.attributes.get("motif","")
            if abs(rpt_len - len(mtf)) <= length_range:
                rmtf, ind, mtf_arg_min = bio.get_least_seq_index(mtf)
                is_rc = False
                if mtf in rpt_test:
                    new_seq = gm.annotations.get_sequence(testchr, rpt.start - context, rpt.end + context)
                elif bio.reverse_complement(mtf) in rpt_test:
                    new_seq = bio.reverse_complement(gm.annotations.get_sequence(testchr, rpt.start - context, rpt.end + context))
                    is_rc = True
                else:
                    continue
                # phase = len(mtf) - mtf_arg_min
                phase = get_phase(rpt_proto, new_seq[context:len(new_seq) - context])
                new_rpt = Repeat(rpt, mtf, rpt.chrom, rpt.start, rpt.end, phase, is_rc, new_seq)
                out_rpts.append(new_rpt)
    
    for rpt in out_rpts:
        for g in gm.annotations.stream_by_types(["gene"], rpt.chr, start = rpt.start, end = rpt.end):
            gn = g.get("info",{}).get("gene_name","")
            if gn:
                rpt.parent_gene_name = gn
                break
    
    return rpt_proto, out_rpts

def get_phase(rpt_proto, rpt_seq):
    
    mlen = len(rpt_proto)
    nrpts = len(rpt_seq)//len(rpt_proto)
    prpt = (rpt_proto * (nrpts + 1))[:len(rpt_seq)]
    corr, rcorr = process.correlate(prpt, rpt_seq, scale = 2*mlen)
    maxshift = (corr.index(corr==max(corr))) % len(rpt_proto)
    # ScalarPlot(corr).show()
    # print(f"maxshfit: {maxshift}")
    return maxshift

def analyze_repeat(motif, rpt):
    
    # motif, _, _ = bio.get_seq_index_abs(rpt.motif)
    
    if not rpt.seq:
        return
    
    print(rpt.repeat)
    
    mlen = rpt.motif_len
    num_rpts = 3
    # test_motif, _, _ = bio.get_least_seq_index(rpt.motif)
    # test_motif = rpt.motif
    test_motif = motif
    dbl_rpt = test_motif * num_rpts
    
    aligner = Align.PairwiseAligner()
    
    scores = []
    algns = []
    for ib in range(mlen * (num_rpts - 1), rpt.seq_len - mlen * (num_rpts - 1)):
        
        subseq = rpt.seq[ib - mlen:ib + mlen]
        score = aligner.score(dbl_rpt, subseq)
        scores.append(score)
        if score == 2*rpt.motif_len:
            algn = aligner.align(dbl_rpt, subseq)[0]
            algns.append(algn)
    
    ScalarPlot(scores, add_range = True).show_chunks()
    # ScalarPlot(scores, add_range = True).show(fmt = " "*mlen*(num_rpts - 1) + "{}")
    # print(rpt.seq)

def analyze_repeat_comp_scales(test_rpt, context):
    
    comps = ["GC","AG","AC"]
    cols = [30, 90, 70]
    seq_lbls = ["upstream","repeat","downstream"]
    
    scale = 4
    step = 1
    keep = None
    shift_step = 1
    
    seqs = [test_rpt.seq[:context], test_rpt.seq[context:test_rpt.seq_len - context], test_rpt.seq[-context-1:]]
    
    for step in [2, 3, 4, 5, 6, 7]:
        print(f"************ step = {step} ************")
        for ns, seq in enumerate(seqs):
            print(seq_lbls[ns])
            for nc,comp in enumerate(comps):
                cgdata = process.correlate_composition(seq,comp_type = comp, scale = scale, step = step, keep = keep, shift_step = shift_step, fill=None)[0]
                ScalarPlot(cgdata, add_range = True, maxval = 1, minval = 0, fg_color = cols[nc]).show(plot_label = comp)
                # print(f"sd = {np.std(cgdata):0.3f}")
        input()

def display_motifs(gm, chr, rpt_context = 128, max_disp = 256, rpt_skip = 10):
    
    mc, ml, mi = get_motif_stats(gm, chr, min_length = 5, combine_rc = True, topk = 1)
    
    mc = sorted(mc.items(), key = lambda k: -k[1])
    
    for mtf, cnt in mc:
        test_motif, rpts = gather_repeat_seqs(gm, mtf, chrs = [chr], context = rpt_context)
        rpts = list(sorted(rpts, key = lambda rpt:-rpt.end + rpt.start))
        nrpts = len(rpts)//rpt_skip
        data = np.zeros((2, 3, max_disp))
        
        print(f"test motif: {test_motif} with {nrpts} repeat instances")
        
        nr = 0
        for test_rpt in rpts[::rpt_skip]:
            
            print(test_rpt.repeat, f"parent gene = {test_rpt.parent_gene}" if test_rpt.parent_gene else "no parent gene :(", f"phase {test_rpt.phase}")
            print(test_rpt.seq[test_rpt.motif_len*2:test_rpt.motif_len*10], "rc" if test_rpt.is_rc else "")
            print(test_rpt.seq[test_rpt.motif_len*2 - test_rpt.phase:test_rpt.motif_len*10 - test_rpt.phase], "rc" if test_rpt.is_rc else "")
            # phase_seq = test_rpt.seq[rpt_context-test_rpt.phase:rpt_context+test_rpt.motif_len*3-test_rpt.phase]
            # print(test_rpt.seq[rpt_context+test_rpt.phase:rpt_context+test_rpt.motif_len+test_rpt.phase], "(+ phase offset)")
            # print(test_rpt.seq[rpt_context-test_rpt.phase:rpt_context+test_rpt.motif_len-test_rpt.phase], "(- phase offset)")
            
            # analyze_repeat(mtf, test_rpt)
            # gm.display_chromosomal_quantities(testchr, ["gc", "ag"], chunksz=chunksz, start = start, length = length, resample = True, maxval = 1)
            
            analyze_repeat_comp_scales(test_rpt, rpt_context)
            
        #     qts, _ = gm.get_chromosomal_quantities(
        #         testchr, 
        #         ["gc", "ag", "ac"], 
        #         chunksz=chunksz, 
        #         start = test_rpt.start - context, 
        #         length = context, 
        #         resample = True, 
        #         do_rc = test_rpt.is_rc
        #     )
        #     data[0,:,:] += np.array(qts)
            
        #     qts, _ = gm.get_chromosomal_quantities(
        #         testchr, 
        #         ["gc", "ag", "ac"], 
        #         chunksz=chunksz, 
        #         start = test_rpt.end, 
        #         length = context, 
        #         resample = True, 
        #         do_rc = test_rpt.is_rc
        #     )
        #     data[1,:,:] += np.array(qts)
        
        # data = data / nrpts
        
        # for i in range(3):
        #     for s in range(2):
        #         print(f"{ilbls[i]}, {slbls[s]}")
        #         ScalarPlot(data[s, i, :], add_range = True, minval = 0, maxval = 1).show()
        #         print()
        # input()


def main():
    
    gm, rpts = load_genome()
    
    # view_repeats(gm, "1", 10e6, motifs_only = True)
    
    # test_motif = "TCAGCC"
    # rpts = extract_repeats(gm, "1", motif = test_motif)
    # print(f"found {len(rpts)} repeats")
    
    # for r in rpts:
    #     print(r)
    
    topk = 3
    testchr = "19"
    
    max_disp = 256
    chunksz = 2
    
    # mtf = "TGTCTC"
    # mtf = "AAAGGG"
    
    # show_motif_symms(mtf)
    # show_motif_symms(2*mtf)
    # show_motif_symms(3*mtf)
    
    view_repeats(gm, "1", 10e6)
    
    return
    
    motif_cts, motif_lens, motif_insts = {}, {}, {}
    for chr in list(range(1, 24)) + ["X","Y"]:
        mc, ml, mi = get_motif_stats(gm, str(chr), min_length = 5, combine_rc = True, topk = 10, rank_by_insts = True)
        motif_cts, motif_lens, motif_insts = pool_motif_stats(mc, ml, mi, motif_cts, motif_lens, motif_insts)
        input()
    
    
    
    # mc = sorted(mc.items(), key = lambda k: -k[1])
    
    # res = input("repeat:")
    # res = "TGTCTC"
    # res = "AGAGAC"
    
    # test_motif = res
    # test_motif, _, _ = bio.get_least_seq_index(test_motif)

    
    # for test_rpt in rpts:
    #     analyze_repeat(test_motif, test_rpt)
    #     for chunksz in [8, 16, 32, 64]:
    #         # chunksz = 16
    #         length = chunksz * max_disp
    #         start = test_rpt.start + test_rpt.seq_len // 2 - length // 2
    #         gm.display_chromosomal_quantity(testchr, "gc", chunksz=chunksz, start = start, length = length, resample = True, maxval = 1)
            
    #         num_rpt_chunks = test_rpt.seq_len // chunksz
    #         num_border_chunks = (max_disp - num_rpt_chunks)//2
            
    #         print("X"*num_border_chunks + '-'*num_rpt_chunks + "X"*num_border_chunks)
    #         input()
    
    # test_rpt = rpts[irpt + rpt_ind]
    
    # print(test_rpt.repeat)
    
    # mlen = test_rpt.motif_len
    # num_rpts = 3
    # test_motif, _, _ = bio.get_least_seq_index(test_rpt.motif)
    # dbl_rpt = test_motif * num_rpts
    
    # aligner = Align.PairwiseAligner()
    
    # scores = []
    # algns = []
    # for ib in range(mlen * (num_rpts - 1), test_rpt.seq_len - mlen * (num_rpts - 1)):
        
    #     subseq = test_rpt.seq[ib - mlen:ib + mlen]
    #     score = aligner.score(dbl_rpt, subseq)
    #     scores.append(score)
    #     if score == 2*test_rpt.motif_len:
    #     # if True:
    #         algn = aligner.align(dbl_rpt, subseq)[0]
    #         algns.append(algn)
    
    # ScalarPlot(scores, add_range = True).show()
    # print(" "*test_rpt.motif_len + test_rpt.seq)
    
    
    # na = 0
    # for i, a in enumerate(algns):
    #     print(f"score: {a.score} ({i} / {test_rpt.seq_len})")
    #     print(str(a))
    #     input()
    
    # bio.ORDER = "GATC"
    
    # seq_len = 5
    # num_seqs = 4**seq_len
    # for i in range(num_seqs):
    #     seq = bio.index_to_seq(i, seq_len = seq_len)
    #     rcseq = bio.reverse_complement(seq)
    #     rcind = bio.get_seq_index(seq)
    #     is_pal = (i + rcind) == num_seqs-1
    #     is_dyad = (bio.get_least_seq_index(seq) == bio.get_least_seq_index(rcseq))
        
        
    #     print(i, seq, "pal" if is_pal else "", "dyad" if is_dyad else "")
    
    
    
    # nts = 512
    # seq_lens = [4, 6, 8, 12]
    # # nts = 1
    # # nts = 1024
    # # seq_lens = [32]
    # for seq_len in seq_lens:
    
    # for i in range(3, 8):
    #     seq_len = i
    #     nts = 4**i
    #     stats, ndyads, nunique = test_rcind_stats(seq_len, nts)
            
    #     print(f"sequence length {seq_len}, with {nts} trials: {ndyads} dyads, {nunique} unique")
    #     for i, v in enumerate(stats):
    #         if not v:
    #             continue
    #         lbl = optimal_ordering_labels[i]
    #         print(f"{lbl}: {v} ({v/nts:0.2%})")
    #     print()
        # print(f"seq len {seq_len}:")
        # print(f"  Fwd: {num_fwd / nts:0.2%}, Comp: {num_rev / nts:0.2%}, Rev: {num_Comp / nts:0.2%}, RC: {num_rc / nts:0.2%}, tie: {num_tie / nts:0.2%}")
    
    # testseq = "AAGCTGATCGA"
    # testseq = "GGGAGG"
    # testseq = "GAGGGGGGGGAG"
    
    # seq_len = 8
    
    # while True:
    #     test_seq = "".join([VOCAB[random.randint(0, 3)] for n in range(seq_len)])
    #     init_ind = bio.get_seq_index_abs(test_seq)
    #     opt_seq, opt_ind = bio.get_seq_least_index(test_seq)
    #     # print(f"{test_seq}, {init_ind} -> {opt_seq}, {opt_ind}")
    #     # input()
        
    pass

if __name__=="__main__":
    main()

