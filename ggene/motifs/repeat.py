
import re

from .motif import BaseMotif

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
    
    def find_instances(self, seq, threshold=0.5):
        """Find tandem repeats and microsatellites"""
        instances = []

        # Find exact repeats
        pattern = f"({self.repeat}){{2,}}"  # At least 2 repeats
        for match in re.finditer(pattern, seq):
            score = len(match.group()) / self.repeat_len  # Number of repeats
            if score >= threshold:
                instances.append((match.start(), match.end(), score))

        # Find approximate repeats (with mismatches)
        # Use dynamic programming or suffix arrays for efficiency

        return instances