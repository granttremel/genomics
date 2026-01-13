
import numpy as np

from ggene.seqs import bio, find, process, align, heal, compare
from ggene.seqs.bio import reverse_complement

class miRNA:
    
    def __init__(self, feat, seq = ""):
        
        self.name = feat.name
        self.chrom = feat.chrom
        self.start = feat.start
        self.end = feat.end
        self.strand = feat.strand
        
        self.seq = seq
        
        self._feature = feat
    
    def infer_guide(self):
        
        
        
        pass
    

def get_gc_bias(seq):
    
    seq_len = len(seq)
    summ = 0
    mean = 0
    var = 0
    for i in range(seq_len):
        bi = seq[i]
        if bi in 'GC':
            summ += 1
            mean += i
            var += i*i
    
    mean = mean / summ
    var = var - mean*mean
    sd = np.sqrt(abs(var))
    return mean/seq_len, sd/seq_len

def infer_guide(seq, rcseq, gcm, rcgcm):
    
    pref_fwd = False
    pref_rc = False
    
    delta = gcm - rcgcm
    if seq[0] in 'AT' or rcseq[0] in 'AT':
        if seq[0] in 'AT' and gcm >= 0.5:
            pref_fwd = True
        
        if rcseq[0] in 'AT' and rcgcm >= 0.5:
            pref_rc = True
        
    else:
        if gcm >= 0.5:
            pref_fwd = True
        elif rcgcm >= 0.5:
            pref_rc = True
    
    if pref_fwd or not pref_rc:
        pref_fwd = delta > 0
        pref_rc = delta <=0
    
    return pref_fwd, pref_rc

def get_mirna_sequence(gm, mirna_feat):
    
    seq = gm.get_sequence(mirna_feat.chrom, mirna_feat.start, mirna_feat.end)
    rcseq = bio.reverse_complement(seq)
    
    if mirna_feat.strand == '-':
        seq, rcseq = rcseq, seq
    
    return seq, rcseq

def analyze_mirna(seq, allow_swap = True, did_swap = False, file = None):
    
    rcseq = reverse_complement(seq)
    
    gcm, gcsd = get_gc_bias(seq)
    rcgcm, rcgcsd = get_gc_bias(rcseq)
    gc_delta = gcm - rcgcm
    
    pref_fwd, pref_rc = infer_guide(seq, rcseq, gcm, rcgcm)
    
    if allow_swap and pref_rc and not pref_fwd:
        return analyze_mirna(rcseq, allow_swap = False, did_swap = True)
    
    print(f"Fwd has GC bias {gcm:0.3f}, RC has GC bias {rcgcm:0.3f}, delta {gc_delta:0.3f}, preferred orientations: fwd {pref_fwd}, rc {pref_rc}", file=file)
    algn = align.align_sequences(seq, rcseq)[0]
    algn.print(chunksz = 256, emph_indels = False, color_subs = True, file=file)
    
    blen = len(algn.target_algn)- algn.gaps
    cnt = round(gcm * blen)
    print(" "*(9+cnt - 1) + "^" , file=file)
    
    align.print_aligned_rna(algn.target_algn, algn.query_algn)
    
    ms = algn.get_mutation_spectrum()
    hm = heal.plot_mutation_spectrum(ms, do_heatmap = True, relative = False, colorbar = False)
    
    
    bs = algn.get_gap_spectrum(separate = False)
    scp = heal.plot_gap_spectrum(bs)
    
    return ms, bs, gcm, rcgcm, did_swap

def display_mirna(feat = {}, feat_name = "", gm = None):
    
    if not feat:
        
        strm = gm.annoations.streams.get("genes")
        for f in strm.stream():
            
            if f.name == feat_name:
                feat = f
                break
    
    seq = feat.attributes.get("seq", "")
    
    if not seq:
        if gm:
            # seq = gm.annotations.get_sequence(feat.chrom, feat.start, feat.start + feat.length//2)
            seq = gm.annotations.get_sequence(feat.chrom, feat.start, feat.end)
        else:
            print(f"miRNA feature {feat} has no bound sequence")
            return
    # else:
    #     seq = seq[:len(seq)//2]
    
    analyze_mirna(seq, allow_swap = True)

