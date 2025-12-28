

from typing import Dict, List, Tuple, Optional, Any
from tabulate import tabulate
from pathlib import Path
import numpy as np

from ggene.seqs import bio
from ggene.motifs.motif import MatchResult

class Jaspar:
    
    bases = "ATGC"
    
    def __init__(self, name, pwm, id = "", path = ""):
        self.name = name
        self.pwm = pwm # (base_index (ATGC 0123), position)
        self.id = id
        self.num_pos = len(self.pwm[0])
        self.path = path
        
        # self.base_map = [bases.index(b) for b in "ATGC"]
        
    @classmethod
    def from_file(cls, file_path):
        
        bases = []
        all_wgts = None
        
        header = None
        
        i=0
        with open(file_path, 'r') as f:
            for line in f:
                
                if not header:
                    header = line
                    continue
                
                ls = line.split()
                base=ls.pop(0).strip()
                _ = ls.pop(0)
                _ = ls.pop(-1)
                
                bases.append(base)
                weights = np.array([int(v.strip()) for v in ls])
                if all_wgts is None:
                    all_wgts = np.zeros((4, len(weights)))
                all_wgts[i] = weights
                i += 1
        
        base_map = [bases.index(b) for b in list(cls.bases)]
        pwm = all_wgts[base_map]
        wsum = pwm.sum(axis = 0)
        pwm_norm = pwm/wsum[None, :]
        
        _id, name = header.split()
        _id = _id.strip().lstrip(">")
        name = name.strip()
        
        return cls(name, pwm_norm, id=_id, path = file_path)
    
    @classmethod
    def from_path(cls, file_path, max_motifs = None):
        
        file_path = Path(file_path).absolute()
        
        jlib = JasparLibrary()
        
        nj = 0
        for jasp in file_path.iterdir():
            nj += 1
            
            if jasp.suffix == ".jaspar":
                
                try:
                    jsp = cls.from_file(str(jasp))
                except Exception as e:
                    print(f"error encountered with jasper at {str(jasp)}: {str(e)}")
                    continue
                
                if not jsp:
                    continue
                
                jlib.add_jaspar(jsp)
            
            if max_motifs and jlib.num_jaspars >= max_motifs:
                break
        
        return jlib
    
    def score(self, seq):
        if len(seq)>self.num_pos:
            seq = seq[:self.num_pos]
        
        basemap = {b:i for i,b in enumerate("ATGC")}
        
        score = 0
        
        for i in range(self.num_pos):
            
            bi = seq[i]
            bind = basemap.get(bi, 0)
            v = self.pwm[bind, i]
            score += v
        
        return score / self.num_pos
        
    def scan(self, seq):
        
        seq_len = len(seq)
        scores = []
        
        for i in range(0, seq_len - self.num_pos):
            ssq = seq[i:i+self.num_pos]
            sc = self.score(ssq)
            scores.append(sc)
        
        return np.array(scores)
    
    @property
    def entropy(self):
        
        entr = 0
        
        for ib in range(self.num_pos):
            
            p_base = self.pwm[:, ib]
            p_base = p_base[p_base > 0]
            
            par_entr = np.sum(-p_base * np.log(p_base))
            entr += par_entr
        
        return entr/self.num_pos

    def to_consensus(self, min_p_ratio = 2):
        
        cons = []
        
        for ib in range(self.num_pos):
            
            p_bases = self.pwm[:, ib]
            
            p_bases = {b:p for b, p in zip(self.bases, p_bases)}
            nb = 4
            
            done = False
            while not done and nb > 1:
                min_b, min_p = min(p_bases.items(), key = lambda k:k[1])
               
                if min_p < 1/nb/min_p_ratio:
                    nb-=1
                    del p_bases[min_b]
                else:
                    done = True 
                
            bs = p_bases.keys()
            b_cons = bio.get_alias(*bs)
            
            cons.append(b_cons)
        
        return "".join(cons)

    def print(self):
        
        tab = []
        tab.append(["Name",self.name])
        tab.append(["ID",self.id])
        tab.append(["Num pos.", self.num_pos])
        tab.append(["Cons", self.to_consensus()])
        
        print(tabulate(tab))
        
        tab = []
        headers=["B"] + [f"P{i}" for i in range(self.num_pos)]
        
        for bi in range(len(self.bases)):
            
            b = self.bases[bi]
            tab.append([b] + list(self.pwm[bi, :]))
            
        print(tabulate(tab, headers = headers, floatfmt = "0.3f"))

class JasparLibrary:
    
    bases = "ATGC"
    
    def __init__(self):
        self.jaspars:Dict[str, Jaspar] = {}
        self.jaspar_ids:List[str] = []
        self.big_pwm = None
        
    @property
    def num_jaspars(self):
        return len(self.jaspars)
    
    def __len__(self):
        return self.num_jaspars
        
    def add_jaspar(self, jaspar):
        self.jaspars[jaspar.id] = jaspar
        self.jaspar_ids.append(jaspar.id)
        self.big_pwm = None
    
    def score_many(self, seq):
        
        scores = []
        
        seq_oh = self.seq_to_onehot(seq)
        
        num_pwms = len(self.jaspars)
        max_npos = max([jsp.num_pos for jsp in self.jaspars.values()])
        
        pwms = np.zeros((num_pwms, 4, max_npos))
        for i, jsp in enumerate(self.jaspars.values()):
            npos = jsp.num_pos
            pwms[i, :, :npos] = jsp.pwm
        
        subseq = seq_oh[:, :max_npos]
        
        for n in range(len(self.jaspars)):
            pwm = pwms[n]
            prod = pwm * subseq
            score = np.sum(prod)
            scores.append(score)
            
        return np.array(scores)

    def scan_many(self, seq):
        
        seq_len = len(seq)
        
        big_pwm, pwm_lens = self.assemble_pwms()
        pwm_len = big_pwm.shape[2]
        scores = np.zeros((big_pwm.shape[0], seq_len))
        
        len_mask = np.ones((big_pwm.shape[0]))
        
        seq_oh = self.seq_to_onehot(seq, len_fill = pwm_len)
        
        for p in range(len(seq)):
            
            len_mask = (pwm_lens <= seq_len - p)
            subseq = seq_oh[:, p:p+pwm_len]
            subseq_exp = subseq[None,:]
            
            prod = big_pwm * subseq_exp
            score = np.sum(prod, axis = (1,2))
            pscore = len_mask * score
            scores[:, p] = pscore
            
        scores /= pwm_lens[:,None]
    
        return scores
    
    def find_instances(self, seq, min_score = 0.85) -> Dict[str,List[MatchResult]]:
        
        all_scores:np.ndarray = self.scan_many(seq)
        
        insts=  {}
        
        # for (_id, name), scores in all_scores.items():
        for i in range(all_scores.shape[0]):
            _id = self.jaspar_ids[i]
            jspr = self.jaspars[_id]
            name = jspr.name
            mtf_len = jspr.num_pos
            
            scores = all_scores[i]
            
            jis = np.argwhere(scores > min_score)
            if not jis.size > 0:
                continue
            matches = []
            for ji in jis:
                start = ji.item()
                end = start + mtf_len
                score = scores[ji].item()
                subseq = seq[start:end]
                
                mres = MatchResult("", start, end, score, name, seq=subseq)
                matches.append(mres)
            
            insts[_id] = matches
        
        return insts
        
    @classmethod
    def seq_to_onehot(cls, seq, len_fill = 0):
        
        oh = np.zeros((4, len(seq) + len_fill))
        
        for i, b in enumerate(seq):
            if b in cls.bases:
                bi = cls.bases.index(b)
                oh[bi, i] = 1
        
        return oh
    
    def get_assembled_pwm(self):
        
        if self.big_pwm is None:
            self.big_pwm = self.assemble_pwms()
        return self.big_pwm
    
    def assemble_pwms(self):
        
        num_pwms = len(self.jaspars)
        max_len = max([jsp.num_pos for jsp in self.jaspars.values()])
        
        big_pwm = np.zeros((num_pwms, len(self.bases), max_len))
        pwm_lens = np.zeros((num_pwms,))
        
        for i,jsp_id in enumerate(self.jaspar_ids):
            
            jsp = self.jaspars[jsp_id]
            
            big_pwm[i, :, :jsp.num_pos] = jsp.pwm
            pwm_lens[i] = jsp.num_pos
        
        return big_pwm, pwm_lens
    
    def print(self, chunksz = 256):
        
        tab = []
        
        for jsp_id in self.jaspar_ids:
            jsp = self.jaspars[jsp_id]
            tab.append([jsp.id, jsp.name, jsp.num_pos, jsp.entropy, jsp.to_consensus(), jsp.path])
        
        print(tabulate(tab, headers = ["ID", "Name", "Num Pos.", "Entropy", "Consensus", "Path"]))