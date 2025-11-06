

from ggene import seqs
from ggene.seqs.bio import reverse_complement
    
# def _seq_cpg(seq, feats, allow_aliases):
    
#     pass

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

def _seq_repeats(seq, feats, rptlen):
    runs, _, _ = seqs.process.correlate_longest_subseq(seq, seq, scale = rptlen+1)
    return sum([1 for r in runs if r==rptlen])

def seq_penta_repeats(seq, feats):
    rptlen = 5
    return _seq_repeats(seq, feats, rptlen)
    # scale = 5
    # runs, _, _ = seqs.process.correlate_longest_subseq(seq, reverse_complement(seq), scale = scale+1)
    # return sum([1 for r in runs if r==scale])

lambda_map = {
    "cpg":seq_cpg,
    "genes":seq_genes,
    "exons":seq_exons,
    "len_cds":seq_len_cds,
    "features":seq_feats,
    "max_run":seq_max_run,
    "penta_repeats":seq_penta_repeats
}
