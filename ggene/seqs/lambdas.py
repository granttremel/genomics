

import re
import regex

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

def _seq_polyn(seq, feats, b, do_rc = False):
    
    bs = [b]
    if do_rc:
        bs.append(reverse_complement(b))
    
    ptrn = "|".join([consensus_to_re("(%s{5,})" % bb) for bb in bs])
    ptrn = re.compile(ptrn)
    
    ms = re.finditer(ptrn, seq)
    
    max_run = 0
    for m in ms:
        st, en = m.span()
        max_run = max(max_run, en-st)
    return max_run


def seq_feats(seq, feats):
    return len(feats)

def seq_genes(seq, feats):
    return sum([1 for f in feats if f.get("type","") == "gene"])

def seq_exons(seq, feats):
    return sum([1 for f in feats if f.get("type","") == "exon"])

def seq_motifs(seq, feats):
    return sum([1 for f in feats if f.get("type","") == "motif"])

def seq_cds_len(seq, feats):
    if not feats:
        return 0
    return sum([f.get("end") - f.get("start") for f in feats if f.get("type","") == "CDS"])/len(feats)

def seq_cds_pct(seq, feats):
    return seq_cds_len(seq, feats) / len(seq)

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
    
    if seq_spec in [seq_genes, seq_exons, seq_motifs, seq_cds_len, seq_cds_pct, seq_feats]:
        return True
    else:
        return False

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
    
    "genes":seq_genes,
    "exons":seq_exons,
    "motifs":seq_motifs,
    
    "cds_len":seq_cds_len,
    "cds_pct":seq_cds_pct,
    "features":seq_feats,
    "max_run":seq_max_run,
    
    "hammerhead_st1": seq_hammerhead_st1,
    
    "tetra_repeats":seq_tetra_repeats,
    "penta_repeats":seq_penta_repeats,
    "hexa_repeats":seq_hexa_repeats,
}
