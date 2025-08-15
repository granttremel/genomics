
from motifs.motif import BaseMotif

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

gcfunc = lambda seq: (seq.count('G') + seq.count('C'))/len(seq)
gcpattern = MetricMotif('gc_content',gcfunc, length=None)
