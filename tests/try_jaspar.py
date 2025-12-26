

import numpy as np
from tabulate import tabulate

from ggene.config import other_paths, DATA_DIR
from ggene.database.genome_manager import GenomeManager
from ggene.database import annotations
from ggene.seqs import gen
from ggene import draw
from ggene.draw.scalar_plot import ScalarPlot

class Jaspar:
    
    bases = "ATGC"
    
    def __init__(self, name, pwm, id = ""):
        self.name = name
        self.pwm = pwm # ATGC
        self.id = id
        self.num_pos = len(self.pwm[0])
        
        # self.base_map = [bases.index(b) for b in "ATGC"]
        
    @classmethod
    def from_file(cls, file_path):
        
        # pwm = {}
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
        wsum = pwm.sum(axis = 1)
        pwm_norm = pwm/wsum[:,None]
        
        _id, name = header.split()
        _id = _id.strip().lstrip(">")
        name = name.strip()
        
        return cls(name, pwm_norm, id=_id)
    
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
        
        return score
        
    def scan(self, seq):
        
        seq_len = len(seq)
        scores = []
        
        for i in range(0, seq_len - self.num_pos):
            ssq = seq[i:i+self.num_pos]
            sc = self.score(ssq)
            scores.append(sc)
        
        return scores
    
    @classmethod
    def score_many(cls, jsprs, seq):
        
        scores = []
        
        seq_oh = cls.seq_to_onehot(seq)
        
        num_pwms = len(jsprs)
        max_npos = max([jsp.num_pos for jsp in jsprs.values()])
        
        pwms = np.zeros((num_pwms, 4, max_npos))
        for i, jsp in enumerate(jsprs.values()):
            npos = jsp.num_pos
            pwms[i, :, :npos] = jsp.pwm
        
        # subseq = seq_oh[:, :max_npos][None, :]
        subseq = seq_oh[:, :max_npos]
        
        for n in range(len(jsprs)):
            pwm = pwms[i]
            prod = pwm * subseq
            score = np.sum(prod)
            scores.append(score)
        # prod = pwms * subseq
        # scores = np.sum(prod, axis = (1,2))
            
        return scores

    @classmethod
    def scan_many(cls, jsprs, seq):
        
        seq_oh = cls.seq_to_onehot(seq)
        
        
        
        pass
    
    @classmethod
    def seq_to_onehot(cls, seq):
        
        oh = np.zeros((4, len(seq)))
        
        for i, b in enumerate(seq):
            bi = cls.bases.index(b)
            oh[bi, i] = 1
        
        return oh

def load_genome():
    gm = GenomeManager()
    return gm

def load_jaspar(filepath):
    jsp = None
    try:
        jsp = Jaspar.from_file(filepath)
    except Exception as e:
        print(str(e))
    
    return jsp

def load_jaspars(jaspar_path, n = None):
    
    jsps = {}
    
    nj = 0
    for jasp in jaspar_path.iterdir():
        nj += 1
        
        if jasp.suffix == ".jaspar":
            
            jsp = load_jaspar(str(jasp))
            
            if not jsp:
                continue
            
            jsps[jsp.id]=jsp
        
        if n and len(jsps) >= n:
            break
    
    print(f"checked {nj} jaspar files")
    
    return jsps

def test_jaspars(jsps):
    for jsp in jsps:
        print(jsp.name)
        print(tabulate(jsp.pwm))
        
        seq = gen.get_random_sequence(5*jsp.num_pos)
        scores = jsp.scan(seq)
        
        ScalarPlot(scores, add_range = True).show()
        
        # for i in range(5):
        #     seq = gen.get_random_sequence(jsp.num_pos)
        #     score = jsp.score(seq)
        #     print(score, seq)
        print()

def main():
    
    gm =load_genome()
    
    jaspar_path = DATA_DIR / "jaspar"
    
    jsps = load_jaspars(jaspar_path, n = 25)
    
    seq = gen.get_random_sequence(256)
    
    res = Jaspar.score_many(jsps, seq)
    print(res)
    # print(len(jsps))
    
    # for n, jsp in jsps.items():
    #     print(jsp.name, jsp.id, jsp.num_pos)




if __name__=="__main__":
    main()

