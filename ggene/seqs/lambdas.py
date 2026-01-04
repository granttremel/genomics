

import re
import regex
import functools

from ggene import seqs
from ggene.seqs.bio import reverse_complement
from ggene.seqs.find import consensus_to_re
from ggene.motifs.motif import MotifDetector

md:MotifDetector = MotifDetector()
md.setup_default_motifs()

def seq_gc(seq, feats):
    return (seq.count("G") + seq.count("C")) / len(seq)

def seq_at(seq, feats):
    return (seq.count("A") + seq.count("T")) / len(seq)

def seq_ag(seq, feats):
    """
    AG ratio
    """
    return (seq.count("A") + seq.count("G")) / len(seq)

def seq_ac(seq, feats):
    """
    AC ratio (amino vs keto)
    """
    return (seq.count("A") + seq.count("C")) / len(seq)

def seq_cpg(seq, feats):
    return seq.count("CG")

def seq_cpg_to_gc(seq, feats):
    return seq.count("CG") / max(1, seq.count("C") + seq.count("G"))

def seq_polya(seq, feats):
    return _seq_polyn(seq, feats, "A", do_rc = True)

def seq_polyt(seq, feats):
    return _seq_polyn(seq, feats, "T", do_rc = True)

def seq_polyg(seq, feats):
    return _seq_polyn(seq, feats, "G", do_rc = True)

def seq_polyc(seq, feats):
    return _seq_polyn(seq, feats, "C", do_rc = True)

def seq_polyy(seq, feats):
    return _seq_polyn(seq, feats, "Y", do_rc = False)

def seq_polyr(seq, feats):
    return _seq_polyn(seq, feats, "R", do_rc = False)

def seq_polys(seq, feats):
    return _seq_polyn(seq, feats, "S", do_rc = False)

def seq_polyw(seq, feats):
    return _seq_polyn(seq, feats, "W", do_rc = False)

@functools.lru_cache(maxsize = 32)
def _get_polyn_pattern(b, do_rc):
    
    bs = [b]
    if do_rc:
        bs.append(reverse_complement(b))
    
    ptrn = "|".join([consensus_to_re("(%s{5,})" % bb) for bb in bs])
    
    return re.compile(ptrn)

def _seq_polyn(seq, feats, b, do_rc = False):
    
    ptrn = _get_polyn_pattern(b, do_rc)
    
    ms = re.finditer(ptrn, seq)
    
    max_run = 0
    for m in ms:
        st, en = m.span()
        max_run = max(max_run, en-st)
    return max_run

def seq_feats(seq, feats):
    return sum(1 for f in feats)

def seq_genes(seq, feats):
    return sum(1 for f in feats if f.feature_type == "gene")

def seq_exons(seq, feats):
    return sum(1 for f in feats if f.feature_type == "exon")

def seq_motifs(seq, feats):
    return sum(1 for f in feats if f.feature_type == "motif")

def seq_pseudo(seq, feats):
    return sum(1 for f in feats if f.feature_type == "pseudogene")

def seq_lncrna(seq, feats):
    return sum(1 for f in feats if f.feature_type == "lncRNA")

def seq_ncrna(seq, feats):
    return sum(1 for f in feats if f.feature_type == "ncRNA")

# ?
def seq_nongenes(seq, feats):
    return sum(1 for f in feats if f.feature_type in ["lncRNA","ncRNA","pseudogene"])

def seq_simple_rpts(seq, feats):
    return sum(1 for f in feats if f.feature_type == "repeat" and f.attributes.get("type","") == "Simple_repeat")

def seq_cds_len(seq, feats):
    if not feats:
        return 0
    return sum([f.get("end") - f.get("start") for f in feats if f.get("type","") == "CDS"])/len(feats)

def seq_cds_pct(seq, feats):
    return seq_cds_len(seq, feats) / len(seq)

def _seq_biotype(seq, feats, biotype):
    
    min_pos = 1e9
    max_pos = 0
    out = 0
    for f in feats:
        if f.attributes.get("gene_biotype","") == biotype:
            out += f.end - f.start
        
        min_pos = min(min_pos, f.start)
        max_pos = max(max_pos, f.end)
    
    return out / (max_pos - min_pos)

def seq_pc_bt(seq, feats):
    """
    protein coding biotype
    """
    return _seq_biotype(seq, feats, "protein_coding")

def seq_lncrna_bt(seq, feats):
    """
    long noncoding RNA biotype
    """    
    return _seq_biotype(seq, feats, "lncRNA")


def seq_max_run(seq, feats):
    # scale = round(np.sqrt(len(seq)))
    scale = 256
    
    maxrun = 0
    for strt in range(0, len(seq), scale):
        subseq = seq[strt:strt+scale].replace("N","")
        runs, _, _ = seqs.process.correlate_longest_subseq(subseq, reverse_complement(subseq), scale = 32)
        maxrun = max(maxrun, max(runs))
    return maxrun

def _seq_motif(seq, feats, motif_name):
    ptrn = md.motifs.get(motif_name)
    if not ptrn:
        return -1
    else:
        return ptrn.count_instances(seq)

def seq_hammerhead_st1(seq, feats):
    return _seq_motif(seq, feats, "hammerhead_stem1")

def _seq_te(seq, feats, te_name):
    return sum([1 for f in feats if f.get("feature_type") == "dfam_hit" and te_name in f.get("name")])

def seq_alu(seq, feats):
    return _seq_te(seq, feats, "Alu")

def seq_aluJ(seq, feats):
    return _seq_te(seq, feats, "AluJ")

def seq_aluY(seq, feats):
    return _seq_te(seq, feats, "AluY")

def seq_aluS(seq, feats):
    return _seq_te(seq, feats, "AluS")

def seq_line1(seq, feats):
    return _seq_te(seq, feats, "L1")

def seq_L1HS(seq, feats):
    return _seq_te(seq, feats, "L1HS")

def seq_L1PA(seq, feats):
    return _seq_te(seq, feats, "L1PA")

def seq_L1P(seq, feats):
    return _seq_te(seq, feats, "L1P")

def seq_L1M(seq, feats):
    return _seq_te(seq, feats, "L1M")

def seq_line2(seq, feats):
    return _seq_te(seq, feats, "L2")

def seq_LTR(seq, feats):
    return _seq_te(seq, feats, "LTR")

def seq_MIR(seq, feats):
    return _seq_te(seq, feats, "MIR")

def seq_SVA(seq, feats):
    return _seq_te(seq, feats, "SVA")

def seq_ALR(seq, feats):
    return _seq_te(seq, feats, "ALR")

def seq_MER(seq, feats):
    return _seq_te(seq, feats, "MER")

def _seq_repeat_motif(seq, feats, motif):
    dblmotif = motif*2
    rcdblmotif = reverse_complement(dblmotif)
    
    res = 0
    
    for f in feats:
        if f.feature_type != "repeat":
            continue
        
        test_mtf = f.attributes.get("motif","")
        if not test_mtf:
            continue
        
        if test_mtf in dblmotif or test_mtf in rcdblmotif:
            res += (f.end - f.start)/len(test_mtf)
    
    return res

def seq_ttctt(seq, feats):
    return _seq_repeat_motif(seq, feats, "TTCTT")

def seq_tataat(seq, feats):
    return _seq_repeat_motif(seq, feats, "TATAAT")

def seq_tatatc(seq, feats):
    return _seq_repeat_motif(seq, feats, "TATATC")

def seq_ttga(seq, feats):
    return _seq_repeat_motif(seq, feats, "TTGA")

def seq_aattaag(seq, feats):
    return _seq_repeat_motif(seq, feats, "AATTAAG")

def seq_atgt(seq, feats):
    return _seq_repeat_motif(seq, feats, "ATGT")

def seq_gaggat(seq, feats):
    return _seq_repeat_motif(seq, feats, "GAGGAT")

def seq_ccctgc(seq, feats):
    return _seq_repeat_motif(seq, feats, "CCCTGC")

def _seq_pattern(seq, feats, ptrn):
    return len(re.findall(ptrn, seq))

def _seq_fuzzy_pattern(seq, feats, ptrn, max_err):
    
    ptrn_str = "(%s){e<=%s}" % (str(ptrn), str(max_err))
    eptrn = regex.compile(ptrn_str)
    return len(regex.findall(eptrn, seq))
    
def _seq_repeats(seq, feats, rptlen):
    seqq = seq.replace("N","")
    runs, _, _ = seqs.process.correlate_longest_subseq(seqq, seqq, scale = rptlen+1)
    return sum([1 for r in runs if r==rptlen])

def seq_tetra_repeats(seq, feats):
    rptlen = 4
    return _seq_repeats(seq, feats, rptlen)

def seq_penta_repeats(seq, feats):
    rptlen = 5
    return _seq_repeats(seq, feats, rptlen)

def seq_hexa_repeats(seq, feats):
    rptlen = 6
    return _seq_repeats(seq, feats, rptlen)

def _seq_repeat(seq, feats, repeat):
    ptrn = re.compile(consensus_to_re(repeat))
    ms = re.findall(ptrn, seq)
    return len(ms)

def _seq_repeat_run(seq, feats, repeat, gap = 0):
    ptrn = re.compile(consensus_to_re(repeat))
    ms = re.finditer(ptrn, seq)
    
    maxr = 0
    r = 0
    last_end = 0
    for m in ms:
        st, en = m.span()
        if st - last_end <= gap:
            r += 1
        elif st < last_end:
            continue
        else:
            r = 0
        maxr = max(maxr, r)
        last_end = en
    
    return maxr

def needs_features(seq_spec):
    
    if seq_spec in lambda_map:
        seq_spec = lambda_map.get(seq_spec)
    
    if seq_spec in [seq_genes, seq_pc_bt]:
        return ["gene"]
    elif seq_spec == seq_pseudo:
        return ["pseudogene"]
    elif seq_spec == "nongenes":
        return ["pseudogene","ncRNA","lncRNA"]
    elif seq_spec in [seq_lncrna_bt, seq_ncrna]:
        return ["ncRNA", "lncRNA"]
    elif seq_spec in [seq_exons]:
        return ["exon"]
    elif seq_spec in [seq_motifs]:
        return ["motif"]
    elif seq_spec in [seq_cds_len, seq_cds_pct]:
        return ["CDS"]
    elif seq_spec in [seq_feats]:
        return ["all"]
    elif seq_spec in [_seq_te, seq_alu, seq_aluJ, seq_aluY, seq_aluS, seq_line1, seq_L1M, seq_L1P, seq_L1HS, seq_L1PA]:
        return ["dfam_hit"]
    elif seq_spec in [_seq_repeat_motif, seq_ttctt, seq_ttga, seq_tataat, seq_tatatc, seq_aattaag, seq_atgt, seq_gaggat, seq_ccctgc]:
        return ["repeat","simple_repeat"]
    elif not seq_spec in lambda_map.values():
        return ["all"]
    else:
        return []

lambda_map = {
    "gc":seq_gc,
    "at":seq_at,
    "ag":seq_ag,
    "ac":seq_ac,
    "cpg":seq_cpg,
    "cpg_to_gc":seq_cpg_to_gc,
    
    "polya":seq_polya,
    "polyt":seq_polyt,
    "polyg":seq_polyg,
    "polyc":seq_polyc,
    "polyy":seq_polyy,
    "polyr":seq_polyr,
    "polys":seq_polys,
    "polyw":seq_polyw,
    
    "features":seq_feats,
    "genes":seq_genes,
    "exons":seq_exons,
    "motifs":seq_motifs,
    "pseudo":seq_pseudo,
    "lncRNA":seq_lncrna,
    "ncRNA":seq_ncrna,
    "nongenes":seq_nongenes,
    "simple_repeats":seq_simple_rpts,
    
    "Alu":seq_alu,
    "AluJ":seq_aluJ,
    "AluY":seq_aluY,
    "AluS":seq_aluS,
    "Line1":seq_line1,
    "L1HS":seq_L1HS,
    "L1PA":seq_L1PA,
    "L1P":seq_L1P,
    "L1M":seq_L1M,
    "L2":seq_line2,
    "LTR":seq_LTR,
    "MIR":seq_MIR,
    "SVA":seq_SVA,
    
    "TTCTT":seq_ttctt,
    "TATAAT":seq_tataat,
    "TATATC":seq_tatatc,
    "TTGA":seq_ttga,
    "AATTAAG":seq_aattaag,
    "ATGT":seq_atgt,
    "GAGGAT":seq_gaggat,
    "CCCTGC":seq_ccctgc,
    
    "cds_len":seq_cds_len,
    "cds_pct":seq_cds_pct,
    "features":seq_feats,
    "protein_coding":seq_pc_bt,
    "lncrna":seq_lncrna_bt,
    
    "max_run":seq_max_run,
    
    "hammerhead_st1": seq_hammerhead_st1,
    
    "tetra_repeats":seq_tetra_repeats,
    "penta_repeats":seq_penta_repeats,
    "hexa_repeats":seq_hexa_repeats,
}
