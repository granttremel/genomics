
import re
from ggene import is_rna, is_dna, aliases_rna, aliases_dna
from .motif import BaseMotif

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
        """Find all instances of the pattern in a sequence.
        
        Args:
            seq: DNA/RNA sequence to search
            threshold: Minimum score threshold (optional)
            
        Returns:
            List of (start, end, score) tuples
        """
        instances = []
        
        # Handle both compiled patterns and string patterns
        if hasattr(self.pattern, 'finditer'):
            # It's a compiled pattern
            matches = self.pattern.finditer(seq)
        else:
            # It's a string pattern
            matches = re.finditer(self.pattern, seq, re.IGNORECASE)
        
        # Find all matches
        for match in matches:
            start = match.start()
            end = match.end()
            matched_seq = match.group()
            
            # Calculate score for this match
            score = self.score_func(matched_seq) if self.score_func else 1.0
            
            # Apply threshold if provided
            if threshold is None or score >= threshold:
                instances.append((start, end, score))
        
        return instances
    

splice_donor_ptrn = 'GGGU[AG]AGU'
splice_donor_motif = PatternMotif("splice_donor", splice_donor_ptrn, lambda x: 1.0)

splice_branch_ptrn = '[CU]U[AG]AC'
splice_branch_motif = PatternMotif("splice_branch", splice_branch_ptrn, lambda x: 1.0)

splice_acceptor_ptrn = '[CU][AUCG]CAGG'
splice_acceptor_motif = PatternMotif("splice_acceptor", splice_acceptor_ptrn, lambda x: 1.0)

# motif_detector = MotifDetector()
# motif_detector.add_motif(gcpattern, 0.5)
# motif_detector.add_motif(splice_donor_motif, 0)
# motif_detector.add_motif(splice_donor_motif, 0)

