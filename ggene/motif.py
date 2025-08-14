
from ggene import aliases_dna, aliases_rna, to_rna, reverse_complement, is_rna, is_dna
from abc import ABC, abstractmethod
import numpy as np
import re

def string_to_re(pattern_string):
    
    if is_rna(pattern_string):
        aliases = aliases_rna
    elif is_dna(pattern_string):
        aliases = aliases_dna
    else:
        aliases = aliases_dna
    
    pattern_list = list(pattern_string)
    
    for i in range(len(pattern_list)):
        if pattern_list[i] in aliases:
            restr = '|'.join(aliases[pattern_list[i]])
            pattern_list[i] = f'[{restr}]'
    
    return ''.join(pattern_list)

def repeat_to_re(repeat_str, subs = {}, max_gap = 0):
    
    repeat_list = list(repeat_str)
    
    for s,v in subs.items():
        pass
    
    repeat_list.append('*')
    if max_gap>0:
        return '[]*'+'.'
    
    pass

class MotifDetector:

    def __init__(self):
        self.motifs={}
    
    def add_motif(self, motif):
        self.motifs[motif.name]=motif
    
    def identify(self, seq, motifs=[]):
        if not motifs:
            motifs=self.motifs.keys()
            
        for m in motifs:
            pass
    
    def score(self, seq, motifs=[]):
        pass

class BaseMotif(ABC):
    
    def __init__(self, name):
        self.name=name
        
    @abstractmethod
    def score(self, seq, return_positions=False):
        
        pass
    
    @abstractmethod
    def find_instances(self, seq, threshold=None):
        """Return list of (start, end, score) tuples"""
        pass
    
    def to_features(self, seq, window_size=None):
        """Convert to ML-ready features"""
        if window_size:
            return self._windowed_features(seq, window_size)
        return self._sequence_features(seq)

class MetricMotif(BaseMotif):
    
    def __init__(self, name, id_function, scoring_function):
        super().__init__(name)
        self.identify_func = id_function
        self.score_func=scoring_function
    
    def score(self, seq, return_positions=False):
        
        return self.score_func(seq)
    
    def find_instances(self, seq, threshold=None):
        
        
        
        pass
    
    def __call__(self, seq):
        return self.func(seq)

class RepeatMotif(BaseMotif):
    
    def __init__(self, name, repeat, scoring_function):
        super().__init__(name)
        self.repeat = repeat
        self.repeat_len= len(self.repeat)
        self.score_func = scoring_function
    
    def score(self, seq):
        
        match_ct={}
        score=0
        
        rptre = re.compile(self.repeat)
        match_ct[0] = len(re.findall(rptre, seq))
        score += match_ct[0]
        
        rpt = list(self.repeat)
        for i in range(self.repeat_len):
            _rpt = rpt.copy()
            _rpt[i] = '.'
            
            rptre = re.compile(''.join(_rpt))
            
            match_ct[i] = len(re.findall(rptre, seq))
            score += match_ct[i]*0.5
        score /= len(seq)
        return score
    
    def to_re(self):
        
        restr_consec = self.repeat + '+'
        restr_nonconsec = self.repeat + '+.*'
        
        pass
    
    def find_instances(self, seq):
        
        matches = []
        rptre = re.compile(self.repeat)
        for m in re.finditer(rptre, seq):
            start = m.start
            end = m.end
            substr = seq[start:end]
            score = self.score(m.groups()[0])
            
            pass
        
        
        pass
    
    

class PatternMotif(BaseMotif):
    
    def __init__(self, name, pattern, scoring_function):
        super().__init__(name)
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        self.pattern=pattern
        self.score_func = scoring_function
    
    def __call__(self, seq):
        res = re.search(self.pattern, seq)
        if res:
            return self.score(res.groups()[0])
        return 0
    
    def score(self, seq, return_positions=False):
        res = re.search(self.pattern, seq)
        return self.score_func(res)
    
    def find_instances(self, seq, threshold=None):
        
        
        pass
    
    

class Motif:
    
    _counter = {}
    
    def __init__(self, name, chrom, start, end, score):
        self.type="motif"
        self.name = name
        ind = self._counter.get(name,0)
        self._counter[name] = ind + 1
        self.sfid=f'{name}-{ind}'
        self.chrom=chrom
        self.start=start
        self.end=end
        self.relative_start = 0
        self.relative_end = self.end - self.start
        
        self.score=score
        
        self.parent_gene=None
        self.parents=[]
        
    def set_parent_gene(self, parent_gene):
        self.parent_gene=parent_gene
        
    def set_parent(self, feature):
        if not feature in self.parents:
            self.parents.append(feature)
            
    def refer_to(self, reference):
        self.start_relative = self.start - reference.start
        self.end_relative = self.end - reference.start
    def __repr__(self) -> str:
        """String representation of the feature."""
        
        parts=[f"id={self.sfid}"]
        parts.append(f"{self.chrom}:{self.start}-{self.end}")
        parts.append(f"score={self.score}")
            
        content = ','.join(parts)
        return f'{type(self).__name__}({content})'

gcfunc = lambda seq: (seq.count('G') + seq.count('C'))/len(seq)
gcpattern = MetricMotif('gc_content',gcfunc, length=None)

repeat1='CA'
repeat2 = 'CCCAT'
repeat3 = ''

splice_donor_ptrn = 'GGGU[AG]AGU'
splice_donor_motif = PatternMotif("splice_donor", splice_donor_ptrn)

splice_branch_ptrn = '[CU]U[AG]AC'
splice_branch_motif = PatternMotif("splice_branch", splice_branch_ptrn)

splice_acceptor_ptrn = '[CU][AUCG]CAGG'
splice_acceptor_motif = PatternMotif("splice_acceptor", splice_acceptor_ptrn)

motif_detector = MotifDetector()
motif_detector.add_motif(gcpattern, 0.5)
motif_detector.add_motif(splice_donor_motif, 0)
motif_detector.add_motif(splice_donor_motif, 0)

