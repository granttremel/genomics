

from ggene import seqs
from ggene.seqs.bio import reverse_complement
from ggene.motifs.motif import MotifDetector

md:MotifDetector = MotifDetector()
md.setup_default_motifs()

def seq_cg(seq, feats):
    return (seq.count("C") + seq.count("G")) / len(seq)

def seq_cpg(seq, feats):
    return seq.count("CG")

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

def _seq_pattern(seq, feats, motif_name):
    ptrn = md.motifs.get(motif_name)
    if not ptrn:
        return -1
    else:
        return ptrn.count_instances(seq)

def seq_hammerhead_st1(seq, feats):
    return _seq_pattern(seq, feats, "hammerhead_stem1")

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

lambda_map = {
    "cg":seq_cg,
    "cpg":seq_cpg,
    
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
