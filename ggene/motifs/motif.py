
from ggene.seqs.bio import to_rna, reverse_complement, is_rna, is_dna
from abc import ABC, abstractmethod
import numpy as np
import re

from .motifio import MotifIO
from ..seqs.find import consensus_to_re

def calculate_significance(self, seq, motif_score):
      """Calculate p-value for motif occurrence"""
      # Shuffle sequence to create null distribution
      # Or use pre-computed background model
      pass

class MotifDetector:

    def __init__(self):
        self.motifs={}
        self.io = MotifIO()
    
    def add_motif(self, motif):
        self.motifs[motif.name]=motif
    
    def identify(self, seq, rcseq = None):
        """
        returns dict: motif_name -> list[(motif_start, motif_end, score, is_rc)]
        """
        
        seq_len = len(seq)
        if not rcseq:
            rcseq = reverse_complement(seq)
        
        all_insts = {}
        
        for motif_name, motif in self.motifs.items():
            
            if not motif_name in all_insts:
                all_insts[motif_name] = []
            
            for motif_start, motif_end, score in motif.find_instances(seq):
                all_insts[motif_name].append((motif_start, motif_end, score, False))
            
            if motif.allow_rc:
                for rcmotif_start, rcmotif_end, score in motif.find_instances(rcseq):
                    all_insts[motif_name].append((seq_len - rcmotif_end, seq_len - rcmotif_start,  score, True))
        
        return all_insts
    
    def score(self, seq, motifs=[]):
        pass

    def setup_default_motifs(self, class_names = []):
        
        if not class_names:
            class_names = motif_classes.keys()
        
        from .pattern import PatternMotif
        
        for k1 in class_names:
            cls_motifs = motif_classes[k1]
            for k2 in cls_motifs:
                motif = default_motifs.get(k2)
                if motif:
                    motif_re = consensus_to_re(motif)
                    try:
                        pm = PatternMotif(k2, motif_re, lambda a:1.0, allow_rc = True)
                    except Exception as e:
                        print("failed to initialize pattern:",k1, k2, motif, motif_re, f"with error {str(e)}")
                        continue
                    self.add_motif(pm)
                

class IndexedMotifDetector:
    """Pre-index genome for fast motif searches"""

    def __init__(self):
        self.suffix_array = None
        self.bwt_index = None

    def build_index(self, genome_sequence):
        """Build suffix array or BWT for O(m) searches"""
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
        return []
    
    def to_features(self, seq, window_size=None):
        """Convert to ML-ready features"""
        if window_size:
            return self._windowed_features(seq, window_size)
        return self._sequence_features(seq)

default_motifs = {
    
    # "splice_donor": "GTRAGT",
    "splice_donor": "YAGGTRAGT",
    "splice_branch": "CTCAY",
    "splice_acceptor": "Y{5,}YAG",
    
    "TATA_box":"TATAYAY",
    "CAAT_box":"GGCCAATCT",
    "E_box":"CAGCTG|CACGTG",
    
    "AU_rich":"ATTTA",
    "Shine-Dalgarno":"AGGAGGT",
    "Kozak":"ACCATGG",
    
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
    
    # "telomerase_pseudoknot_hTR44-184": "GGCTAGGC[GCCGT]GCTTTT[GCTCC]CCGCGCGCTG[TTTTTCTC]GCTGACTTTCAGCGGGCGGAAAAGC[CTCG]GCCTGCCGCC[TTCCACCGTTCATTCTAGAGCAAACAAA]AAATGTCAGCT",
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
    
    # "msat_1": "AAGCATGGACCATTCCTTCAGGATGGAGCCTGCCAAGTGTGGACCATTCCTTCAGGATGCAGCCAGGTAAGCATCAGCCATTCCTTCATAATGCGGCCAGGTAAGCAT",
    "msat_telomerase": "TTAGGG",
    "msat_telomerase-2": "CTAAC",
    
    "misc1": "TATATATATGGGAGA",
    "misc1": "TATATATATGTGGGAAA",
    "misc1": "TATATATACGGGAAA",
    "misc1": "TATATATATGGGAGA",
}

motif_classes = {
    "splice":["splice_donor","splice_acceptor"],
    
    "promoter":["TATA_box","CAAT_box","E_box", "AU_rich","Shine-Dalgarno","Kozak"],
    
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
        "msat_telomerase",
        "msat_telomerase-2",
    ],
    
    "misc":[],
}


