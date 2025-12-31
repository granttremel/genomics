"""Regex pattern-based motif detection."""

import re
from typing import List, Dict, Optional, Callable
import numpy as np

from ggene.motifs.motif import BaseMotif, MatchResult, MotifLibrary
from ggene.seqs.find import consensus_to_re

class PatternMotif(BaseMotif):
    """Regex pattern-based motif.

    Uses compiled regex patterns for efficient sequence scanning.
    Supports IUPAC ambiguity codes through consensus_to_re conversion.
    """

    def __init__(self, name: str, pattern: str,
                 scoring_function: Optional[Callable] = None,
                 motif_id: str = "",
                 allow_rc: bool = True,
                 motif_class: str = ""):
        super().__init__(name, motif_id=motif_id, allow_rc=allow_rc, motif_class=motif_class)

        # Compile pattern if string
        if isinstance(pattern, str):
            self._pattern_str = pattern
            self.pattern = re.compile(pattern, re.IGNORECASE)
        else:
            self._pattern_str = pattern.pattern
            self.pattern = pattern

        self.score_func = scoring_function or (lambda x: 1.0)

    @property
    def length(self) -> int:
        """Estimate pattern length (may vary for patterns with quantifiers)."""
        # Remove regex special chars and estimate
        clean = re.sub(r'[^ATGCUNRYWSMKHBVD]', '', self._pattern_str.upper())
        return max(len(clean), 1)

    def scan(self, seq: str) -> np.ndarray:
        """Score each position (1.0 if pattern matches at position, else 0)."""
        scores = np.zeros(len(seq))
        for match in self.pattern.finditer(seq):
            scores[match.start()] = self.score_func(match.group())
        return scores

    def find_instances(self, seq: str, threshold: float = 0.0,
                      chrom: str = "", offset: int = 0) -> List[MatchResult]:
        """Find all pattern matches."""
        results = []

        for match in self.pattern.finditer(seq):
            matched_seq = match.group()
            score = self.score_func(matched_seq)

            if score >= threshold:
                results.append(MatchResult(
                    start=offset + match.start() if offset else match.start(),
                    end=offset + match.end() if offset else match.end(),
                    score=score,
                    name=self.name,
                    motif_id=self.id,
                    strand="+",
                    seq=matched_seq,
                    chrom=chrom,
                ))

        return results

    def count_instances(self, seq: str) -> int:
        """Count pattern matches."""
        return len(self.pattern.findall(seq))


class PatternLibrary(MotifLibrary):
    """Collection of pattern motifs with batch scanning."""

    def __init__(self):
        self._patterns: Dict[str, PatternMotif] = {}
        self._pattern_ids: List[str] = []

    @property
    def motif_ids(self) -> List[str]:
        return self._pattern_ids

    def add(self, motif: BaseMotif) -> None:
        """Add a pattern motif."""
        if not isinstance(motif, PatternMotif):
            raise TypeError(f"Expected PatternMotif, got {type(motif)}")
        self._patterns[motif.id] = motif
        if motif.id not in self._pattern_ids:
            self._pattern_ids.append(motif.id)

    def add_pattern(self, name: str, pattern: str,
                   scoring_function: Optional[Callable] = None,
                   motif_class: str = "") -> PatternMotif:
        """Convenience method to add a pattern by name/regex."""
        motif = PatternMotif(name, pattern, scoring_function, motif_class=motif_class)
        self.add(motif)
        return motif

    def get_motif(self, motif_id: str) -> PatternMotif:
        return self._patterns[motif_id]

    def scan_all(self, seq: str) -> Dict[str, np.ndarray]:
        """Scan with all patterns."""
        return {pid: self._patterns[pid].scan(seq) for pid in self._pattern_ids}

    @classmethod
    def from_dict(cls, patterns: Dict[str, str],
                  motif_class: str = "") -> 'PatternLibrary':
        """Create library from {name: pattern} dict."""
        lib = cls()
        for name, pattern in patterns.items():
            lib.add_pattern(name, pattern, motif_class=motif_class)
        return lib

    def load_defaults(self, class_names = []):
        
        if not class_names:
            class_names = pattern_classes.keys()
        
        for clsn in class_names:
            
            ptrn_names = pattern_classes.get(clsn, [])
            
            for ptrn_name in ptrn_names:
                pre_ptrn = default_patterns.get(ptrn_name)
                ptrn = consensus_to_re(pre_ptrn)
                
                ptrn_mtf = PatternMotif(ptrn_name, ptrn, allow_rc = False, motif_class = clsn)
                
                self.add(ptrn_mtf)
    

default_patterns = {
    
    "splice_donor": "YAGGTRAGT", # ACTYACCTR
    "splice_branch": "CTCAY",
    "splice_acceptor": "Y{5,}YAG",
    
    "TATA_box":"TATAYAY",
    "CAAT_box":"GGCCAATCT",
    "E_box":"CAGCTG|CACGTG",
    
    "AU_rich":"ATTTA",
    # "Shine-Dalgarno":"AGGAGGT", # prokaryotes and archaea ... 
    "Kozak":"ACCATGG",
    "PolyA":"A{7,}|T{7,}",
    
    "hammerhead":"YYRRGCCGTTACCTRCAGCTGATGAGCTCCAARAAGAGCGAAACCNRNYAGGTCCTGYAGTAYTGGCYNR",
    "hammerhead_stem1":"RGCCGN{45,65}TGGCY", # nominal 55
    "hammerhead_stem2":"CTRCAGN{30,42}CTGYAG", # nominal 36
    "hammerhead_big_loop":"GCTGATGAGN{12,20}CGAAAAN{9,14}TCC",
    
    "SRP_Alu_stemloop_3":"GCCGGGCGCG(N{6,12})CGTGCCT",
    "SRP_Alu_stemloop_4":"GTCCC(N{8,16})GGGAG",
    "SRP_Alu_stem_3":"CGGGCGCG",
    "SRP_Alu_stem_5a":"GGCTGAGGCTGGA",
    "SRP_stem_5d":"CTGGGCTGTAGTGNGCTATGC",
    "SRP_stemloop_5d":"CTGGGCTGTAGTGCGCTATGC(N{134,174})GAATAGCCACTGCACTCCAG", # nominal 154, total 195
    "SRP_S_stemloop_5e":"CGGGTGTCCGCACTA(N{97,137})TAGTGGGATCGCGCCTG", # nominal 117, total 149
    "SRP_S_stemloop_5f":"GTTCGGCAT(N{80,110})GTGCTGATC", # nominal 95, total 114
    "SRP_S_stemloop_6":"CAATATGGTGACCTCCCG(N{3,6})CGGGGGACCACCAGGTTG", # nominal 4, total 40
    "SRP_S_stemloop_8a":"AGGGGTG(N{26,38})AACTCCC", # nominal 32, total 48
    "SRP_S_stemloop_8b":"GGCCCAGGTCG(N{3,6})CGGAGCAGGTC", # nominal 4, total 26
    
    "telomerase_pseudoknot": "CGACTGTAAAAAAN{10,20}GCGGGCGACTTTCAGTCGCTCTTTTTGTCGCGCGC",
    
    "telomerase_pseudoknot_hTR44-184": "GGCTAGGC(N{3,7})GCTTTT(N{3,7})CCGCGCGCTG(N{4,12})GCTGACTTTCAGCGGGCGGAAAAGC(N{2,6})GCCTGCCGCC(N{22,34})AAATGTCAGCT",
    "telomerase_pseudoknot_hTR44-184-1": "GGCTAGGC(N{3,7})GCTTTT",
    "telomerase_pseudoknot_hTR44-184-2": "GCTTTT(N{3,7})CCGCGCGCTG",
    "telomerase_pseudoknot_hTR44-184-3": "CCGCGCGCTG(N{4,12})GCTGACTTT",
    "telomerase_pseudoknot_hTR44-184-4": "GCTGACTTTCAGCGGGCGGAAAAGC",
    "telomerase_pseudoknot_hTR44-184-5": "CAGCGGGCGGAAAAGC",
    "telomerase_pseudoknot_hTR44-184-6": "AAAAGC(N{2,6})GCCTGCCGCC",
    "telomerase_pseudoknot_hTR44-184-7": "GCCTGCCGCC(N{22,34})AAATGTCAGCT",
    
    "telomerase_pseudoknot_hTR44-184-loop5": "TTCCACCGTTCATTCTAGAGCAAACAAA",
    
    "telomerase_pseudoknot_hTR171-184": "AAAAAATGTCAGCT",
    
    "msat2": "CACGTG",
    "msat_telomerase": "TTAGGG",
    "msat_telomerase-2": "CTAAC",
    
    "rand1":"YYYRYYYRYYYRYYYRRRYYYYRRYRRYRRYY",
    # "misc1": "TATATATATGGGAGA",
    # "misc1": "TATATATATGTGGGAAA",
    # "misc1": "TATATATACGGGAAA",
    # "misc1": "TATATATATGGGAGA",
}

pattern_classes = {
    "splice":["splice_donor","splice_acceptor"],
    
    "promoter":[
        "TATA_box",
        "CAAT_box",
        "E_box", 
        "PolyA",
        "AU_rich",
        # "Shine-Dalgarno",
        "Kozak"
    ],
    
    "hammerhead":["hammerhead","hammerhead_stem1","hammerhead_stem2","hammerhead_big_loop"],
    
    "SRP":[
        "SRP_Alu_stemloop_3",
        "SRP_Alu_stemloop_4",
        "SRP_Alu_stem_3",
        "SRP_Alu_stem_5a",
        "SRP_stem_5d",
        "SRP_stemloop_5d",
        "SRP_S_stemloop_5e",
        "SRP_S_stemloop_5f",
        "SRP_S_stemloop_6",
        "SRP_S_stemloop_8a",
        "SRP_S_stemloop_8b",
    ],
    
    "SRP_Alu":[
        "SRP_Alu_stemloop_3",
        "SRP_Alu_stemloop_4",
        "SRP_Alu_stem_3",
        "SRP_Alu_stem_5a",
    ],
    "SRP_S":[
        "SRP_S_stemloop_5e",
        "SRP_S_stemloop_5f",
        "SRP_S_stemloop_6",
        "SRP_S_stemloop_8a",
        "SRP_S_stemloop_8b",
    ],
    
    "pseudoknot":[
        "telomerase_pseudoknot",
        "telomerase_pseudoknot_hTR44-184",
        "telomerase_pseudoknot_hTR44-184-1",
        "telomerase_pseudoknot_hTR44-184-2",
        "telomerase_pseudoknot_hTR44-184-3",
        "telomerase_pseudoknot_hTR44-184-4",
        "telomerase_pseudoknot_hTR44-184-5",
        "telomerase_pseudoknot_hTR44-184-6",
        "telomerase_pseudoknot_hTR44-184-7",
        
        "telomerase_pseudoknot_hTR44-184-loop5",
        
        "telomerase_pseudoknot_hTR171-184",
    ],
    
    "msat":[
        # "msat_1",
        "msat2",
        "msat_telomerase",
        "msat_telomerase-2",
    ],
    
    "misc":[
        "rand1",
        
        
    ],
}


