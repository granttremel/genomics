
import random
import numpy as np

from ggene.seqs import bio
from ggene.seqs.bio import reverse_complement, reverse, complement
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE, ALIASES, ALIASES_FULL

class Base:
    
    _bases = ORDER[:4]
    
    def __init__(self, value = None, fdict = {}):
        
        self.pdict = {}
        if fdict:
            self.num_obs = int(sum(fdict.values()))
            self.pdict = {k:v/self.num_obs for k,v in fdict.items() if v}
        elif value:
            self.num_obs = 1
            self.pdict[value] = 1.0
        else:
            raise ValueError
        
        self.index = 0
        self.alias = bio.get_minimal_alias2(*self.pdict.keys())
        self.parent = None
    
    @property
    def RC(self):
        fdict = {complement(b):p*self.num_obs for b, p in self.pdict.items()}
        return Base(fdict = fdict)
    
    @property
    def mode(self):
        maxp = 0
        maxb = ""
        for b, p in self.pdict.items():
            if p > maxp:
                maxp = p
                maxb = b
        return maxb
    
    @property
    def entropy(self):
        ent = 0
        for k, p in self.pdict.items():
            if p > 0:
                ent += -p*np.log(p)
        return ent
    
    def joint_entropy(self, other_base):
        
        ent = 0
        for ba in self._bases:
            for bb in self._bases:
                pa = self.pdict.get(ba, 0)
                pb = other_base.pdict.get(bb, 0)
                if pa and pb:
                    fact = 1 if ba == bb else 0.5
                    ent += -fact*pa*pb*np.log(fact*pa*pb) # yeah idk
    
        return ent
    
    def mutual_information(self, other_base):
        return self.entropy + other_base.entropy - self.joint_entropy(other_base)
    
    def KL_div(self, other_base):
        kldiv = 0
        for b in self._bases:
            pa = self.pdict.get(b,0)
            pb = other_base.pdict.get(b,0)
            
            if (pa and not pb) or (pb and not pa):
                kldiv +=1
                continue
            elif not pa and not pb:
                continue
            
            kldiv += pa * np.log(pa / pb)
        return kldiv
    
    def KL_div_symm(self, other_base):
        return self.KL_div(other_base) + other_base.KL_div(self)
    
    def lambda_div(self, other_base, lamb):
        
        PQ_pdict = {k:lamb*self.pdict.get(k,0) + (1-lamb)*other_base.pdict.get(k,0) for k in self._bases}
        PQ_num_obs = int(lamb*self.num_obs + (1-lamb)*other_base.num_obs)
        lambPQ = Base.from_probs(PQ_pdict, PQ_num_obs)
        
        return lamb*self.KL_div(lambPQ) + (1-lamb)*other_base.KL_div(lambPQ)
        
    def JS_div(self, other_base):
        return self.lambda_div(other_base, 0.5)
        
    def set_index(self, ind):
        self.index = ind
    
    def set_parent(self, seq):
        self.parent = seq
    
    def update(self, **fdict):
        no = self.num_obs + sum(fdict.values())
        old_fdict = {k:v*self.num_obs for k,v  in self.pdict.items()}
        new_fdict = {k: fdict.get(k, 0) + old_fdict.get(k,0) for k in self._bases}
        
        self.num_obs = no
        self.pdict = {k:v/self.num_obs for k, v in new_fdict.items() if v}
        self.alias = bio.get_minimal_alias(*self.pdict.keys())
    
    def combine(self, other_base, in_place = False):
        if in_place:
            other_fdict = {b:other_base.pdict.get(b, 0) * other_base.num_obs for b in self._bases}
            self.update(**other_fdict)
        else:
            fdict = {b:self.pdict.get(b,0)*self.num_obs + other_base.pdict.get(b,0)*other_base.num_obs for b in self._bases}
            return Base(fdict = fdict)
    
    def trim(self, p_min, keep_obs = False):
        
        delobs = 0
        delbs = []
        for b, p in self.pdict.items():
            if p < p_min:
                delobs += p*self.num_obs
                delbs.append(b)
        
        new_obs = self.num_obs if keep_obs else (self.num_obs - delobs)
        robs = float(new_obs) / self.num_obs
        
        self.num_obs = int(new_obs)
        
        for b in delbs:
            del self.pdict[b]
        
        for k in self.pdict:
            self.pdict[k] = robs*self.pdict[k]
    
    def copy(self):
        return Base.from_probs(self.pdict, self.num_obs)
    
    def __repr__(self):
        parts = []
        for k, v in self.pdict.items():
            parts.append(f"{k}:{v:0.2f}")
        parts.append(f"n={self.num_obs}")
        return f"Base({",".join(parts)})"

    def crepr(self):
        parts = []
        for k, v in self.pdict.items():
            if v == 1.0:
                return k
            else:
                parts.append(f"{k}:{v:0.2f}")
        return f"[{",".join(parts)}]"

    def __str__(self):
        return self.alias
    
    def __len__(self):
        return 1
    
    def __eq__(self, other):
        return self.alias == other.alias
    
    def to_consensus(self):
        return self.alias
    
    def sample(self):
        rn = random.random()
        pmarg = 0
        for b, p in self.pdict.items():
            pmarg += p
            if rn < pmarg:
                return b
        return ""
    
    @classmethod
    def from_probs(cls, pdict, num_obs):
        return cls(fdict = {k:num_obs*p for k,p in pdict.items()})



class Sequence:
    
    _seqs = {}
    
    def __init__(self, bases=[]):
        self.bases = [b if isinstance(b, Base) else Base(value = b) for b in bases]
        self.refresh()
        
    def substitute(self, pos, base):
        self.bases[pos]= base
        
    @property
    def RC(self):
        bases = [b.RC for b in reversed(self.bases)]
        return Sequence(bases=bases)
    
    @property
    def consensus(self):
        return "".join(b.to_consensus() for b in self.bases)
    
    @property
    def num_obs(self):
        return max([b.num_obs for b in self.bases])
    
    def trim(self, min_p, keep_obs = False):
        for b in self.bases:
            b.trim(min_p, keep_obs = keep_obs)
    
    def sample(self):
        return "".join(b.sample() for b in self.bases)
    
    def append(self, base):
        if not isinstance(base, Base):
            base = Base(base)
        self.bases.append(base)
        base.set_index(len(self))
        base.set_parent(self)
    
    def merge(self, other_seq, in_place = False):
        if in_place:
            for i in range(len(self)):
                self.substitute(i, self.bases[i].combine(other_seq.bases[i], in_place = False))
        else:
            return Sequence(bases = [ba.combine(bb) for ba, bb in zip(self.bases, other_seq.bases)])
    
    def hamming(self, other_seq, mode = "alias"):
        """
        modes:
          alias
          mode
          kl_div
          js_div
        """
        
        if mode == "alias":
            return self._hamming_alias(self.bases, other_seq.bases)
        elif mode == "mode":
            return self._hamming_mode(self.bases, other_seq.bases)
        elif mode == "kl_div":
            return self._hamming_kl_div(self.bases, other_seq.bases)
        elif mode == "js_div":
            return self._hamming_js_div(self.bases, other_seq.bases)
        else:
            return -1
    
    def compare(self, other_seq, allow_alias = False, err_tol = 0):
        
        print(f"comparing {repr(self)} with {repr(other_seq)}")
        
        err = 0
        for ba, bb in zip(self.bases, other_seq.bases):
            
            if ba == bb:
                pass
            elif allow_alias:
                print(f"comparing {repr(ba)} with {repr(bb)}")
                baali = ba.alias
                bbali = bb.alias
                print(f"ba alias {baali}, bb alias {bbali}") 
                
                if ba.alias in ALIASES_FULL.get(bb.alias):
                    print(f"equal")
                    pass
                elif bb.alias in ALIASES_FULL.get(ba.alias):
                    print(f"equal")
                    pass
                else:
                    err += 1
            else:
                err += 1
            
            if err > err_tol:
                return err
        return err
    
    def fit_consensus(self, consensus):
        err = self.compare(consensus, allow_alias = True, err_tol = 0)
        return err==0
    
    @classmethod
    def _hamming_mode(cls, bases_a, bases_b):
        dist = 0
        for ba, bb in zip(bases_a, bases_b):
            if ba.mode != bb.mode:
                dist += 1
        return dist
    
    @classmethod
    def _hamming_kl_div(cls, bases_a, bases_b):
        dist = 0
        for ba, bb in zip(bases_a, bases_b):
            dist += ba.KL_div(bb)
        return dist
    
    @classmethod
    def _hamming_js_div(cls, bases_a, bases_b):
        dist = 0
        for ba, bb in zip(bases_a, bases_b):
            dist += ba.JS_div(bb)
        return dist
    
    @classmethod
    def _hamming_alias(cls, bases_a, bases_b):
        dist = 0
        for ba, bb in zip(bases_a, bases_b):
            balis = bio.ALIASES.get(ba.alias)
            bbalis = bio.ALIASES.get(bb.alias)
            ints = set(balis).symmetric_difference(bbalis)
            dist += len(ints) / 2
            
        return dist
    
    def refresh(self):
        for i, b in enumerate(self.bases):
            b.set_index(i)
            b.set_parent(self)
    
    def copy(self):
        return Sequence(bases = [b.copy() for b in self.bases])
    
    def slice(self, start, end, keep_ref = False):
        new_seq = Sequence()
        new_seq.bases = self.bases[start:end]
        if not keep_ref:
            new_seq = new_seq.copy()
            new_seq.refresh()
        return new_seq
    
    def isolate(self):
        hsh = hash(self)
        if hsh in self._seqs:
            return self._seqs.get(hsh)
        else:
            self._seqs[hsh] = self
            return self
    
    def __len__(self):
        return len(self.bases)

    def __getitem__(self, ind):
        bs = self.bases[ind]
        if len(bs) > 1:
            return Sequence(bases = bs)
        else:
            return bs

    def __getitems__(self, indexer):
        return Sequence(self.bases[indexer])
    
    def __str__(self):
        return "".join(str(b) for b in self.bases)
    
    def __repr__(self):
        max_bases =256
        bases = self.bases[:max_bases]
        parts = []
        for b in bases:
            parts.append(repr(b))
        return f"Sequence({",".join(parts)})"
    
    def crepr(self):
        bstrs = [b.crepr() for b in self.bases]
        return f"Seq({"".join(bstrs)})"
    
    def __eq__(self, other):
        return self.consensus == other.consensus
    
    def __hash__(self):
        return hash(self.consensus)
    
    @classmethod
    def from_string(cls, seq_str):
        bases = []
        for b in seq_str:
            base = Base(b)
            bases.append(base)
        return Sequence(bases=bases)
    

