

from abc import ABC, abstractmethod
import numpy as np
import re

class MotifDetector:

    def __init__(self):
        
        #function from string to float on [0,1]
        # self.match_criterion = match_criterion
        self.motif={}
        self.thresholds={}
    
    def add_motif(self, motif, threshold):
        self.motif[motif.name]=motif
        self.thresholds[motif.name]=threshold
    
    def qualify(self, seq, motifs=[]):
        
        if not motifs:
            motifs = self.motif.keys()
        
        results={}
        
        for pn in motifs:
            res = self.motif[pn](seq)
            if res > self.thresholds[pn]:
                results[pn]=res
                
        return results
    

class BaseMotif(ABC):
    
    def __init__(self, name):
        self.name=name
        
    @abstractmethod
    def evaluate(self, seq):
        pass
    
    @abstractmethod
    def __call__(self, seq):
        return self.evaluate(seq)
    
    @abstractmethod
    def scan(self, seq, window, stride=1):
        pass
    
    @abstractmethod
    def find_max(self, seq):
        pass

class MetricMotif(BaseMotif):
    
    def __init__(self, name, func):
        super().__init__(name)
        self.func=func
        self.summary_stat=np.mean
        
    def __call__(self, seq):
        return self.func(seq)
    
    def scan(self, seq, window, stride=1):
        if window<1:
            return self(seq)
        outinds=list()
        outvs = list()
        for i in range(window//2,len(seq)-window//2,stride):
            subseq=seq[i:i+window]
            v = self.func(subseq)
            outinds.append((i,i+window))
            outvs.append(v)
        return outinds,outvs
    
    def find_max(self, seq, stride=1):
        inds,vs = self.scan(seq,stride=stride)
        imax = np.argmax(vs)[0]
        vmax = vs[imax]
        return imax, vmax
    
class PatternMotif(BaseMotif):
    
    def __init__(self, name, pattern, score=None):
        super().__init__(name)
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        self.pattern=pattern
        if not score:
            self.score = lambda res:1
        else:
            self.score=score
    
    def __call__(self, seq):
        res = re.search(self.pattern, seq)
        if res:
            return self.score(res.groups()[0])
        return 0
    
    def scan(self, seq, window, stride=1):
        matchiter = re.finditer(self.pattern,seq)
        outinds=list()
        outvs=list()
        for res in matchiter:
            outinds.append(res.start)
            outvs.append(self.score(res.groups()[0]))
        return outinds, outvs
    
    def find_max(self, seq, window, stride=1):
        inds,vs = self.scan(seq,window,stride=stride)
        imax = np.argmax(vs)[0]
        vmax = vs[imax]
        return imax, vmax
    

class Motif:
    def __init__(self, name, start, end, score):
        self.type="motif"
        self.sfid=name
        self.start=start
        self.end=end
        self.score=score
        
        self.parent_gene=None
        self.parents=[]
        
    def set_parent_gene(self, parent_gene):
        self.parent_gene=parent_gene
        
    def set_parent(self, feature):
        if not feature in self.parents:
            self.parents.append(feature)


gcfunc = lambda seq: (seq.count('G') + seq.count('C'))/len(seq)
gcpattern = MetricMotif('gc_content',gcfunc, length=None)


splice_donor_ptrn = 'GG(G)U[AG]AGU'
splice_donor_motif = PatternMotif("splice_donor", splice_donor_ptrn)


motif_detector = MotifDetector()
motif_detector.add_motif(gcpattern, 0.5)
motif_detector.add_motif(splice_donor_motif, 0)

