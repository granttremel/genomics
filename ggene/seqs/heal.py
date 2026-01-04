
import numpy as np

from ggene import seqs
from ggene import draw
from ggene.draw import Heatmap, ScalarPlot
from ggene.seqs import bio, find, process, manipulate
from .bio import reverse_complement
from .vocab import VOCAB
import itertools


class Healer:
    
    def __init__(self, gm, chr, start, seq_len):
        
        self.gm = gm
        self.chr = chr
        self.start= start
        self.pos = start
        self.base_seq = gm.get_sequence(chr, start, start+seq_len)
        self.seqs = [self.base_seq]
        self.rcbase_seq = reverse_complement(self.base_seq)
        self.seq_len = len(self.base_seq)
        self.motifs = []
        self.changes = []
        
        self.runs, self.run_inds, self.run_shifts = process.correlate_longest_subseq(self.base_seq, self.base_seq)
        self.rcruns, self.rcrun_inds, self.rcrun_shifts = process.correlate_longest_subseq(self.base_seq, self.rcbase_seq)
        self.corr, self.rcorr = process.correlate(self.base_seq, self.base_seq)
        self.corr8, self.rcorr8 = process.correlate(self.base_seq, self.base_seq, scale = 8)
        self.corr16, self.rcorr16 = process.correlate(self.base_seq, self.base_seq, scale = 16)
        self.corr32, self.rcorr32 = process.correlate(self.base_seq, self.base_seq, scale = 32)
        
        self.templates = {}
        self.rctemplates = {}
    
    def advance(self, seq_len = None):
        if not seq_len:
            seq_len = self.seq_len
        self.pos += seq_len
        new_seq = self.gm.get_sequence(self.chr, self.pos, self.pos+seq_len)
        self.seqs.append(new_seq)
        # extract templates..?
    
    def broaden(self, length):
        seq = self.gm.get_sequence(self.chr, self.pos - length, self.pos + self.seq_len + length)
        return seq
    
    def locate_repeats(self, zscore = 2, span = 8, do_rc = False, do_low = True, do_high = True):
        
        if do_rc:
            data = self.rcruns
            inds = self.rcrun_inds
            shifts = self.rcrun_shifts
            seq = self.rcbase_seq
        else:
            data = self.runs
            inds = self.run_inds
            shifts = self.run_shifts
            seq = self.base_seq
        
        cmean = np.mean(data)
        csd = np.std(data)
        
        maxthr = cmean + zscore*csd
        
        highs = []
        
        for i,(v, j, sh) in enumerate(zip(data, inds, shifts)):
            if i == self.seq_len//2:
                continue
            
            if v >= maxthr and do_high:
                rpseq, ineg, ipos = self.expand_repeat(j, v, sh)
                rplen = ipos - ineg
                highs.append((ineg, rplen, v, rpseq))
            
        return highs
    
    def expand_repeat(self, index, length, shift):
        
        done = False
        i = index+length
        while not done:
            b1 = self.base_seq[i]
            if i+shift >= self.seq_len:
                b2 = self.base_seq[i-shift]
            else:
                b2 = self.base_seq[i+shift]
            if b1 == b2:
                pass
            else:
                done = True
                break
            
            i += 1
            if i >= self.seq_len:
                break
        ipos = min(i, self.seq_len)
        
        done = False
        i = index
        while not done:
            
            b1 = self.base_seq[i]
            if i+shift >= self.seq_len:
                b2 = self.base_seq[i-shift]
            else:
                b2 = self.base_seq[i+shift]
            if b1 == b2:
                pass
            else:
                done = True
                break
            
            i -= 1
            if i <= 0:
                break
        ineg = max(i,0)
        return self.base_seq[ineg:ipos], ineg, ipos
    
    def add_motif(self, motif):
        self.motifs.append(motif)
    
    def identify_motif(self, motif, err_tol = 1):
        locs = []
        mods = []
        mlen = len(motif)
        for i in range(self.seq_len - mlen):
            test_seq = self.base_seq[i:i+mlen]
            err = find.compare_sequences(test_seq, motif, err_tol = err_tol)
            if err <= err_tol:
                locs.append((i, err))
                if err > 0:
                    mods.append(test_seq)
        
        return locs, mods
            
def identify_mutation(ba, bb):
    
    m = bio.ALIASES_REV.get(ba + bb)
    if not m:
        m = bio.ALIASES_REV.get(bb + ba)
    
    if not m:
        return "e"
    else:
        return m

def identify_mutation_type(ba, bb):
    m = identify_mutation(ba, bb)
    
    if m in "YR":
        return "I"
    elif m in "SW":
        return "C"
    elif m in "KM":
        return "V"
    else:
        return " "

def get_mutation_spectrum(seqa, seqb, mutation_spec = {},  b_is_rc = False, allowed_mutations = 'RYSWKM'):
    
    cons = bio.merge(seqa, seqb)
    rccons = bio.merge(seqa, bio.complement(seqb))
    
    ms = {k:mutation_spec.get(k, 0) for k in allowed_mutations}
    
    for i,b in enumerate(cons):
        
        if b in VOCAB:
            continue
        
        real_b = b
        # if b_is_rc:
        #     real_b = rccons[i]
        # else:
        #     real_b = b
        
        if b in ms:
            ms[real_b] += 1
        
        for k in allowed_mutations:
            v = bio.ALIASES.get(k)
            if b in v and not k == 'N':
                ms[real_b] += 1
    
    return ms

def convert_mutation_spectrum(mutation_spec, aliases):
    
    alias_map = bio.get_alias_map(aliases)
    new_ms = {k+kk:0 for k,kk in itertools.product(alias_map.values(), alias_map.values()) if k!=kk}
    
    for k, v in mutation_spec.items():
        b, bb  = alias_map.get(k[0], "N"), alias_map.get(k[1], "N")
        if b==bb:
            continue
        new_ms[b+bb] += v
    
    return new_ms

def make_mutation_spec_relative(mutation_spec, mutation_rate = None):
    
    mtn_spc_rel = {}
    
    num_mtns = len(mutation_spec)
    
    if not mutation_rate:
        mutation_rate = sum(mutation_spec.values())
    
    for mt, r in mutation_spec.items():
        
        rr = r / mutation_rate - 1/num_mtns 
        mtn_spc_rel[mt] = rr
    
    return mtn_spc_rel

def plot_mutation_spectrum(mutation_spec, do_heatmap = False, aliases = "", relative = True, suppress = False, **kwargs):
    """
    for dist, first line is "to", second line is "from"
    
    """
    
    if aliases:
        mutation_spec = convert_mutation_spectrum(mutation_spec, aliases)
        kwargs["ab"] = aliases
    
    if relative:
        mutation_spec = make_mutation_spec_relative(mutation_spec)
        kwargs['center'] = None
        kwargs['minval'] = None
        kwargs['symmetric_color'] = True
    
    if do_heatmap:
        plt = _mutation_heatmap(mutation_spec, suppress=suppress, **kwargs)
    else:
        plt = _mutation_dist(mutation_spec, suppress=suppress, **kwargs)
    
    return plt

def _mutation_heatmap(mutation_spec, suppress = False, **kwargs):
    
    hm_data = []
    
    kwargs["center"] = kwargs.pop("center", None)
    kwargs["minval"] = kwargs.pop("minval", 0)
    kwargs["col_width"] = kwargs.pop("col_width", 2)
    kwargs["col_space"] = kwargs.pop("col_space", 0)
    kwargs["row_height"] = kwargs.pop("row_height", 1)
    kwargs["row_space"] = kwargs.pop("row_space", 0)
    kwargs["symmetric_color"] = kwargs.pop("symmetric_color", False)
    kwargs["colorbar"] = kwargs.get("colorbar", True)
    
    ab = kwargs.get("ab", VOCAB)
    
    for b in ab:
        hm_row = []
        for bb in ab:
            hm_row.append(mutation_spec.get(b + bb))
        hm_data.append(hm_row)
    
    hm = draw.Heatmap(hm_data, row_labels = ab, col_labels = ab, **kwargs)
    hm.show(suppress = suppress)
    print()
    
    return hm
    

def _mutation_dist(mutation_spec, suppress = False, **kwargs):
    
    data = []
    lbl1s = []
    lbl2s = []
    
    minval = kwargs.pop("minval", 0)
    ab = kwargs.get("ab", VOCAB)
    
    for b in ab:
        for bb in ab:
            
            if b==bb:
                continue
            
            data.append(mutation_spec[b+bb])
            lbl2s.append(bb)
        
        data.append(minval)
        lbl1s.append(b+" "*(len(ab)-1))
        lbl2s.append(" ")
    
    scd = ScalarPlot(data, minval=minval, add_range = True, **kwargs)
    scd.show(suppress = suppress)
    
    if not suppress:
        print("".join(lbl2s))
        print("".join(lbl1s))
        print()
    
    return scd

def get_variant_data(vars, subs = {}, ins = [], dels = []):
    
    for var in vars:
        ref = var.ref
        try:
            alt0, alt1 = var.alt
        except:
            print(f"failed to get alts from variant {var}")
            print(f"var raw line: {var.attributes.get("_raw_line","")}")
            continue
        
        if len(ref) == len(alt0):
            sub = ref + alt0
            if not sub in subs:
                subs[sub] = 0
            subs[sub] += 1
        else:
            if len(ref) > len(alt0):
                dels.append(ref[1:])
            else:
                ins.append(alt0[1:])
    
    return subs, ins, dels

def show_variant_data(vars = [], subs = {}, ins = [], dels = [], suppress = False):
    
    if not subs:
        subs, ins, dels = get_variant_data(vars)
    
    # plot_mutation_spectrum(subs, do_heatmap = True, relative = False, color_scheme = 'ember_glow', add_middle = False)
    plt = plot_mutation_spectrum(subs, do_heatmap = True, relative = True, color_scheme = "moss", suppress = suppress, colorbar = False)
    
    mean_ins= np.mean([len(i) for i in ins])
    mean_del= np.mean([len(d) for d in dels])
    
    ins_comps = {b:0 for b in "ATGC"}
    del_comps = {b:0 for b in "ATGC"}
   
    for i in ins:
        for b in i:
            ins_comps[b] += 1
    
    for d in dels:
        for b in d:
            del_comps[b] += 1
    
    sum_ins = sum(ins_comps.values())
    sum_dels = sum(del_comps.values())
    ins_comps = {b:v/sum_ins for b, v in ins_comps.items()}
    del_comps = {b:v/sum_dels for b, v in del_comps.items()}
    
    print(f"Insertions: mean length {mean_ins:0.1f}")
    sc_ins = ScalarPlot(ins_comps, mode = "distribution", labels = list('ATGC'), key_order = 'ATGC', space = 1, add_range = True)
    sc_ins.show(suppress=suppress)
    
    print(f"Deletions: mean length {mean_del:0.1f}")
    sc_dels = ScalarPlot(del_comps, mode = "distribution", labels = list('ATGC'), key_order = 'ATGC', space = 1, add_range = True)
    sc_dels.show(suppress=suppress)
    
    return plt, sc_ins, sc_dels
