
import re
from ggene.seqs.find import consensus_to_re, reverse_complement_re
from .motif import BaseMotif

class PatternMotif(BaseMotif):
    
    def __init__(self, name, pattern, scoring_function, allow_rc = True, motif_class = ""):
        super().__init__(name)
        
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        
        self.pattern=pattern
        self.score_func = scoring_function
        self.allow_rc = allow_rc
        self.motif_class = motif_class
        
    def __call__(self, seq):
        res = re.search(self.pattern, seq)
        if res:
            return self.score(res.groups()[0])
        return 0
    
    def score(self, seq, return_positions=False):
        res = re.search(self.pattern, seq)
        return self.score_func(res)
    
    def count_instances(self, seq):
        return len(self.pattern.findall(seq))
    
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
