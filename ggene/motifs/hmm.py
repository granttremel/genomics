
from typing import Any, Optional, List, Dict
from dataclasses import dataclass

import pyhmmer
from pyhmmer.plan7 import HMM, Pipeline, Domain, TopHits, Hit, Builder, Background
from pyhmmer.easel import DigitalSequence, DigitalSequenceBlock, MSA, DigitalMSA

from ggene.config import DATA_DIR
from ggene import draw
from .motif import BaseMotif

import h5py



hmm_path = DATA_DIR / "dfam" / "human_repeats.hmm"

AB = pyhmmer.easel.Alphabet.dna()
BG = Background(AB, uniform = True)

@dataclass
class SearchResult:
    target_name:str
    query_name:str
    score:float
    norm_score:float
    
    target_to:int
    target_from:int
    target_length:int
    
    
    query_to:int
    query_from:int
    query_length:int
    
    target_cov: float = None
    query_cov: float = None
    
    _target_seq:str = ""
    _query_seq:str = ""
    _id_seq:str = ""
    
    _best_domain:Any=None
    _alignment:Any = None
    
    @classmethod
    def from_hit(cls, hit:Hit, target_seq:str = ""):
        
        algn = hit.best_domain.alignment
        
        idseq = hit.best_domain.alignment.identity_sequence
        idlen = len(idseq.replace("+","").replace(" ",""))
        norm_score = hit.score / hit.best_domain.alignment.target_length
        
        target_cov = idlen/algn.target_length
        query_cov = idlen/algn.hmm_length
        
        return SearchResult(
            hit.name.decode(),
            hit.best_domain.alignment.hmm_name.decode(),
            hit.score,
            norm_score,
            
            algn.target_to,
            algn.target_from,
            algn.target_length,
            algn.hmm_to,
            algn.hmm_from,
            algn.hmm_length,
            
            target_cov = target_cov,
            query_cov = query_cov,
            
            _target_seq = target_seq,
            _id_seq = idseq,
            _alignment = algn
        )
    
class HMMMotif(BaseMotif):
    
    def __init__(self, name, hmm, scoring_function, allow_rc = True, motif_class = ""):
        super().__init__(name)
        
        self.hmm = hmm
        self.score_func = scoring_function
        self.allow_rc = allow_rc
        self.motif_class = motif_class
    
    @property
    def length(self):
        return self.hmm.M
    
    
    
    def __call__(self, seq):
        return 0
    
    def score(self, seq, return_positions=False):
        return None
    
    def count_instances(self, seq):
        return 0
    
    def find_instances(self, seq, threshold=None):
        
        insts = []
        
        res = search_hmm(self.hmm, [seq])
        
        for hit in res:
            
            if threshold and hit.score < threshold:
                continue
            algn = hit.best_domain.alignment
            insts.append((algn.target_from, algn.target_to, hit.score))
        
        return insts

def get_hmmmotifs(names = [], classes = [], families = []):
    
    motifs = []
    
    hmms = load_hmms()
    
    for hmm in hmms:
        
        name = hmm.name.decode("utf-8")
        cn, fn = get_repeat_class(name)
        if classes and cn not in classes:
            continue
        
        if families and not fn in families:
            continue
        
        if names and name not in names:
            continue
        
        new_motif = HMMMotif(name, hmm, None)
        motifs.append(new_motif)
    
    return motifs


def load_hmms(hmm_path = hmm_path, name_filt = None):
    
    hmms = []
    
    name_filt = bytes(name_filt, encoding="utf-8") if name_filt else ""
    
    with pyhmmer.plan7.HMMFile(hmm_path) as hmf:
        for hh in hmf:
            if name_filt and not name_filt in hh.name:
                continue
            
            hmms.append(hh)
    return hmms

hmm_cls = ["Alu","LTR","MER", "L1", "L2", "L3", "L4", "L5", "CR1", "Charlie", "Tigger", "Arthur", "Zaphod", "Ricksha", "Kanga", "MamGyp", "MST","UCON", "HERV", "Eulor", "MLT", "tRNA", "Eutr", "EUTRE", "MARE", "MADE", "ERVL", "U", "HAL", "HSAT", "GSAT", "THE1", "HUERS", "X", "hAT", "SVA", "PRIMA", "MamTip", "MamRep", "MIR", "PABL", "HY", "Eut", "COMP", "BSR"]

def load_hmm_classes(classes = []):
    
    out_hmms = {}
    
    if not classes:
        classes = hmm_cls
    
    hmms = load_hmms()
    
    for hmm in hmms:
        name = hmm.name.decode("utf-8")
        
        cn, fn = get_repeat_class(name)
        
        if not cn in classes:
            continue
        
        if not cn in out_hmms:
            out_hmms[cn] = {}
        
        cnd = out_hmms[cn]
        
        if not fn in cnd:
            cnd[fn] = []
        
        cnd[fn].append(hmm)
    
    return out_hmms

def get_repeat_class(name):
    
    for cls_name in sorted(hmm_cls, key = lambda k:-len(k)):
        
        if name.startswith(cls_name):
            
            fam_name = name.removeprefix(cls_name)
            
            return cls_name, fam_name
    
    return "other", ""

def load_alus():
    return load_hmms(name_filt = "Alu")

def load_line1s():
    return load_hmms(name_filt = "L1")
    
def check_truncation(repeat, hmm):
    
    rpt_len = repeat.end - repeat.start
    rpt_hmm_start = repeat.attributes.get("start_hmm", 0)
    rpt_hmm_end = repeat.attributes.get("end_hmm", 0)
    
    hmm_len = len(hmm.consensus)
    
    if repeat.strand == '-':
        rpt_hmm_start, rpt_hmm_end = rpt_hmm_end, rpt_hmm_start
    
    tr_5p = rpt_hmm_start / hmm_len
    tr_3p = 1 - rpt_hmm_end / hmm_len
    
    # if repeat.strand == '-':
    #     tr_5p, tr_3p = tr_3p, tr_5p
    
    repeat.attributes["motif_len"] = hmm_len
    repeat.attributes["motif_start"] = repeat.start - rpt_hmm_start
    repeat.attributes["motif_end"] = repeat.end + rpt_hmm_end
    repeat.attributes["5p_truncation"] = tr_5p
    repeat.attributes["3p_truncation"] = tr_3p
    
    return repeat

def search_hmm(hmm, seqs):
    
    dseqs = []
    for seq in seqs:
        
        seq_enc = AB.encode(seq)
        
        dseq = pyhmmer.easel.DigitalSequence(AB, sequence = seq_enc)
        dseqs.append(dseq)
    
    
    dsb = pyhmmer.easel.DigitalSequenceBlock(AB, iterable = dseqs)
    ppl = pyhmmer.plan7.Pipeline(hmm.alphabet)
    res = ppl.search_hmm(hmm, dsb)
    
    return res

def print_repeats(rpts, start_offset = 0):
    
    for rpt in rpts:
        
        parts = []
        
        parts.append(f"repeat_len={rpt.end - rpt.start}")
        
        motif_len = rpt.get("motif_len", 0)
        if motif_len:
            parts.append(f"motif_len={motif_len}")
        
        for k in ["bits","e-value","bias","kimura_div", "start_hmm", "end_hmm"]:
            v = rpt.attributes.get(k, None)
            if v:
                parts.append(f"{k}={v}")
        
        for k in ["5p_truncation", "3p_truncation"]:
            
            v = rpt.attributes.get(k, 0)
            if v > 0.1:
                parts.append(f"{k}={v:0.0%}")
        
        start_hmm = rpt.attributes.get("start_hmm",0)
        end_hmm = rpt.attributes.get("end_hmm", 0)
        
        inst_pos_str = f"{rpt.chr}:{rpt.start}-{rpt.end}"
        # inst_pos_off_str = f"{rpt.start-start_offset}-{rpt.end-start_offset}"
        # hmm_pos_str = f"{start_hmm}-{end_hmm}/{motif_len}"
    
        print(f"{rpt.name} at {inst_pos_str}, ({", ".join(parts)})")

def search_hmms(seqs:Dict[str,str], hmms, norm_score_thresh = 0.5):
    
    scores = {}
    
    dseqs = []
    names = []
    for name, seq in seqs.items():
        dseq = DigitalSequence(AB, sequence = AB.encode(seq), name = name.encode())
        dseqs.append(dseq)
        names.append(name)
        
    db = DigitalSequenceBlock(AB, dseqs)
    
    ppl = Pipeline(AB)
    
    for hmm_name, hmm in hmms.items():
        
        res = ppl.search_hmm(hmm, db)
        
        if not res:
            continue
        
        max_res = max(res, key = lambda r:r.score)
        if max_res:
            target_seq = seqs.get(max_res.name.decode(), "")
            sres = SearchResult.from_hit(max_res, target_seq = target_seq)
            
            if sres.norm_score < norm_score_thresh:
                continue
            
            if not hmm.name in scores:
                scores[hmm_name] = []
            scores[hmm_name].append(sres)
        
    return scores
    

def make_msa(seqs:Dict[str,str]):
    
    dseqs = []
    for name, seq in seqs.items():
        
        dseq = DigitalSequence(AB, name = name.endoce("utf-8"), sequence = AB.encode(seq))
        dseqs.append(dseq)
    
    dsb = DigitalSequenceBlock(AB, dseqs)
    
    msa = DigitalMSA(AB, sequences = dsb)
    
    return msa

def build_hmm(msa=None, seqs:Dict[str,str] = {}):
    
    if not msa:
        msa = make_msa(seqs)
    
    builder = Builder(AB)
    hmm, prof, optprof = builder.build_msa(msa, BG)
    
    return hmm, prof, optprof
    
    


