
from typing import List, Dict, Any

from .base import Base, Sequence

class NGramTree:
    
    def __init__(self, n, ns, la = 0):
        
        self.n = n
        self.ns = ns
        self.la = la
        
        self.vocab:List[Sequence] = []
        self.ngs:Dict[Sequence,int] = {}
        self.txs:Dict[Sequence,List[Sequence]] = {}
        self.base_seq:str = ""
    
    def from_sequence(self, seq):
        
        if isinstance(seq, str):
            seq = Sequence.from_string(seq)
            
        self.base_seq = seq.consensus
        
        self._get_ngrams(seq)
        self._get_ngram_transitions(seq)
    
    def _get_ngrams(self, seq):
        
        self.ngs = {}
        
        for i in range(len(seq) - self.n):
            ng = seq[i:i+self.n]
            if not ng in self.ngs:
                self.vocab.append(ng.isolate())
                self.ngs[ng.isolate()] = 0
            self.ngs[ng] += 1
            
        return self.ngs
    
    def _get_ngram_transitions(self, seq):
        
        if not self.la:
            self.la = self.n
    
        self.txs = {}
        
        for i in range(len(seq)-self.n-self.ns-self.la):
            
            sref = seq[i:i+self.n]
            ssq = seq[i+self.la:i+self.la+self.ns]
            
            if not sref in self.txs:
                self.txs[sref.isolate()] = []
            else:
                sref.isolate().merge(sref, in_place = True)
            
            if ssq in self.txs[sref]:
                for ngs in self.txs[sref]:
                    if ssq == ngs:
                        ngs.merge(ssq,in_place = True)
                        break
            else:
                self.txs[sref].append(ssq)
            
        return self.txs
    
    def generalize_ngram(self, cons):
        
        if isinstance(cons, str):
            cons = Sequence.from_string(cons)
        
        # gather subsequents from consensus
    
        for ref in self.txs.keys():
            
            sqs = self.txs.get(ref)
            
            delvals = []
            f_agg = 0
            
            for sq in sqs:
                if sq.fit_consensus(cons):
                    
                    delvals.append(sq)
                    cons.merge(sq, in_place = True)
                    print(f"merged {sq} into {cons}")
                    # f_agg += sqs.get(sq, 0)
                
            for dv in delvals:
                sqs.remove(dv)
            sqs.append(cons)
    
        
    def print(self):
        
        ngtx = self.txs
        for sr in ngtx:
            if sr.num_obs < 10:
                continue
            print(f"{sr.crepr()}: {sr.num_obs}")
            for ssq in ngtx[sr]:
                print(f"  {ssq.crepr()}: {ssq.num_obs}")
            print()
        pass
    
    def to_array(self):
        
        
        
        pass
    
    