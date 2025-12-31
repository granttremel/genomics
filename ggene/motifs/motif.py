"""Base classes for motif detection.

This module defines the core interfaces for motif detection:
- MatchResult: Standardized result for motif matches
- BaseMotif: Interface for individual motif definitions
- MotifLibrary: Interface for collections with batch operations

Concrete implementations:
- pattern.py: PatternMotif, PatternLibrary (regex-based)
- pwm.py: PWMMotif, PWMLibrary (position weight matrices)
- hmm.py: HMMMotif, HMMLibrary (hidden Markov models)

Integration with UGenomeAnnotations:
- database/motifs.py: MotifStream classes that use these libraries
"""

from typing import List, Dict, Tuple, Optional, Any, Iterator, TYPE_CHECKING
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np

if TYPE_CHECKING:
    from ggene.database.ufeature import UFeature


@dataclass
class MatchResult:
    """Standardized result for a motif match.

    Attributes:
        start: Start position in scanned sequence (0-based)
        end: End position in scanned sequence (0-based, exclusive)
        score: Match score (interpretation depends on motif type)
        name: Motif name
        motif_id: Unique motif identifier
        strand: Strand of match ('+' or '-')
        seq: Matched sequence
        chrom: Chromosome (set when converting to genomic coordinates)
    """
    start: int
    end: int
    score: float
    name: str
    motif_id: str = ""
    strand: str = "+"
    seq: str = ""
    chrom: str = ""

    def to_genomic(self, chrom: str, offset: int) -> 'MatchResult':
        """Convert to genomic coordinates.

        Args:
            chrom: Chromosome name
            offset: Genomic position of seq[0] (1-based)
        """
        return MatchResult(
            start=offset + self.start,
            end=offset + self.end,
            score=self.score,
            name=self.name,
            motif_id=self.motif_id,
            strand=self.strand,
            seq=self.seq,
            chrom=chrom,
        )

    def to_ufeature(self, feature_type: str = "motif",
                    source: str = "Motif") -> 'UFeature':
        """Convert to UFeature for streaming integration."""
        from ggene.database.ufeature import UFeature
        return UFeature({
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end - 1,  # UFeature uses inclusive end
            'score': self.score,
            'name': self.name,
            'id': self.motif_id,
            'strand': self.strand,
            'feature_type': feature_type,
            'source': source,
            'matched_seq': self.seq,
        })


class BaseMotif(ABC):
    """Base class for individual motif definitions.

    Subclasses implement different motif types:
    - PatternMotif: Regex patterns
    - PWMMotif: Position weight matrices
    - HMMMotif: Hidden Markov models

    Attributes:
        name: Human-readable motif name
        id: Unique identifier
        allow_rc: Whether to also scan reverse complement
        motif_class: Category/class of motif (e.g., "promoter", "splice")
    """
    name: str
    id: str
    allow_rc: bool = True
    motif_class: str = ""

    def __init__(self, name: str, motif_id: str = "",
                 allow_rc: bool = True, motif_class: str = ""):
        self.name = name
        self.id = motif_id or name
        self.allow_rc = allow_rc
        self.motif_class = motif_class

    @property
    @abstractmethod
    def length(self) -> int:
        """Expected/typical motif length (for overlap calculation in chunked scanning)."""
        pass

    @abstractmethod
    def scan(self, seq: str) -> np.ndarray:
        """Score each position in sequence.

        Args:
            seq: DNA/RNA sequence to scan

        Returns:
            Array of scores, length = len(seq) - self.length + 1
        """
        pass

    def find_instances(self, seq: str, threshold: float = 0.0,
                      chrom: str = "", offset: int = 0) -> List[MatchResult]:
        """Find all matches above threshold.

        Args:
            seq: Sequence to scan
            threshold: Minimum score threshold
            chrom: Chromosome for genomic coordinates (optional)
            offset: Genomic offset for coordinates (optional)

        Returns:
            List of MatchResult objects
        """
        scores = self.scan(seq)
        results = []

        for i, score in enumerate(scores):
            if score >= threshold:
                end = i + self.length
                result = MatchResult(
                    start=offset + i if offset else i,
                    end=offset + end if offset else end,
                    score=float(score),
                    name=self.name,
                    motif_id=self.id,
                    strand="+",
                    seq=seq[i:end],
                    chrom=chrom,
                )
                results.append(result)

        return results

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name!r}, len={self.length})"


class MotifLibrary(ABC):
    """Base class for motif collections with batch operations.

    Libraries provide efficient batch scanning across many motifs.
    Subclasses can optimize using vectorization (numpy), FFT convolution, etc.

    Usage:
        library = PWMLibrary()
        library.load("/path/to/jaspar/")

        # Scan sequence with all motifs
        for result in library.find_all_instances(seq, threshold=0.9):
            print(result.name, result.start, result.score)

        # Integrate with UGenomeAnnotations
        stream = library.to_stream()
        annotations.add_motifs(stream, "pwm")
    """

    @property
    @abstractmethod
    def motif_ids(self) -> List[str]:
        """List of motif IDs in this library."""
        pass

    @property
    def num_motifs(self) -> int:
        """Number of motifs in library."""
        return len(self.motif_ids)

    @property
    def max_length(self) -> int:
        """Maximum motif length (for chunked scanning overlap)."""
        return max((self.get_motif(mid).length for mid in self.motif_ids), default=0)

    @abstractmethod
    def add(self, motif: BaseMotif) -> None:
        """Add a motif to the library."""
        pass

    @abstractmethod
    def get_motif(self, motif_id: str) -> BaseMotif:
        """Get a motif by ID."""
        pass

    @abstractmethod
    def scan_all(self, seq: str) -> Dict[str, np.ndarray]:
        """Scan sequence with all motifs.

        Args:
            seq: Sequence to scan

        Returns:
            Dict mapping motif_id -> score array
        """
        pass

    def find_all_instances(self, seq: str, threshold: float = 0.0,
                          chrom: str = "", offset: int = 0) -> Iterator[MatchResult]:
        """Find instances of all motifs above threshold.

        Args:
            seq: Sequence to scan
            threshold: Minimum score threshold
            chrom: Chromosome for genomic coordinates
            offset: Genomic offset for coordinates

        Yields:
            MatchResult objects sorted by position
        """
        all_scores = self.scan_all(seq)
        results = []

        for motif_id, scores in all_scores.items():
            motif = self.get_motif(motif_id)
            for i, score in enumerate(scores):
                if score >= threshold:
                    end = i + motif.length
                    results.append(MatchResult(
                        start=offset + i if offset else i,
                        end=offset + end if offset else end,
                        score=float(score),
                        name=motif.name,
                        motif_id=motif_id,
                        strand="+",
                        seq=seq[i:end],
                        chrom=chrom,
                    ))

        # Sort by position for heapq.merge compatibility
        results.sort(key=lambda r: (r.start, r.end))
        yield from results

    def to_stream(self, feature_type:str, source_name:str, min_score: float = 0.0):
        """Create a MotifStream for UGenomeAnnotations integration.

        Args:
            min_score: Minimum score threshold for matches

        Returns:
            A MotifStream instance bound to this library
        """
        from ggene.database.motifs import LibraryMotifStream
        return LibraryMotifStream(self, min_score=min_score, feature_type = feature_type, source_name = source_name)

    def __len__(self) -> int:
        return self.num_motifs

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.num_motifs} motifs)"


# =============================================================================
# Legacy code below - kept for backwards compatibility
# =============================================================================

# Legacy imports
try:
    from ggene.seqs.bio import reverse_complement
    from .motifio import MotifIO
    from ..seqs.find import consensus_to_re
    _LEGACY_IMPORTS_OK = True
except ImportError:
    _LEGACY_IMPORTS_OK = False


class MotifDetector:
    """Legacy motif detector - use MotifLibrary subclasses instead."""

    def __init__(self):
        self.motifs = {}
        if _LEGACY_IMPORTS_OK:
            self.io = MotifIO()
    
    def add_motif(self, motif):
        self.motifs[motif.name]=motif
    
    def identify(self, seq, rcseq = None):
        """
        returns dict: motif_name -> list[(motif_start, motif_end, score, is_rc)]
        """
        
        if not rcseq and _LEGACY_IMPORTS_OK:
            rcseq = reverse_complement(seq)
        
        all_insts = {}
        
        for motif_name, motif in self.motifs.items():
            
            if not motif_name in all_insts:
                all_insts[motif_name] = []
            
            for res in motif.find_instances(seq):
                all_insts[motif_name].append(res)
            
            if motif.allow_rc:
                for res in motif.find_instances(rcseq):
                    all_insts[motif_name].append(res)
        
        return all_insts
    
    def score(self, seq, motifs=[]):
        pass

    def setup_default_motifs(self, pattern_classes = [], hmm_classes = []):
        
        if pattern_classes:
            self.setup_default_patterns(pattern_classes)
        
        if hmm_classes:
            self.setup_default_hmms(hmm_classes)
        
    
    def setup_default_hmms(self, hmm_classes):
        
        from .hmm import get_hmmmotifs
        
        hmm_mtfs = get_hmmmotifs(classes = hmm_classes)
        
        for hmm in hmm_mtfs:
            self.add_motif(hmm)
            
    
    def setup_default_patterns(self, class_names=[]):
        if not class_names:
            class_names = motif_classes.keys()
        
        from .pattern import PatternMotif
        for k1 in class_names:
            cls_motifs = motif_classes[k1]
            for k2 in cls_motifs:
                motif = default_motifs.get(k2)
                if motif and _LEGACY_IMPORTS_OK:
                    motif_re = consensus_to_re(motif)
                    try:
                        pm = PatternMotif(k2, motif_re, scoring_function=lambda _: 1.0,
                                         allow_rc=True, motif_class=k1)
                    except Exception as e:
                        print("failed to initialize pattern:", k1, k2, motif, motif_re,
                              f"with error {str(e)}")
                        continue
                    self.add_motif(pm)
            
        
    def get_motif_class(self, motif_name = "", motif = ""):
        
        if motif:
            for name, _motif in default_motifs.items():
                if motif == _motif:
                    motif_name = name
                    break
        
        for mtf_cls, mtfs in motif_classes.items():
            if motif_name in mtfs:
                return mtf_cls

        return ""

class IndexedMotifDetector:
    """Pre-index genome for fast motif searches"""

    def __init__(self):
        self.suffix_array = None
        self.bwt_index = None

    def build_index(self, genome_sequence):
        """Build suffix array or BWT for O(m) searches"""
        pass


default_motifs = {
    
    # "splice_donor": "GTRAGT",
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
    "msat2": "CACGTG",
    "msat_telomerase": "TTAGGG",
    "msat_telomerase-2": "CTAAC",
    
    "rand1":"YYYRYYYRYYYRYYYRRRYYYYRRYRRYRRYY",
    # "misc1": "TATATATATGGGAGA",
    # "misc1": "TATATATATGTGGGAAA",
    # "misc1": "TATATATACGGGAAA",
    # "misc1": "TATATATATGGGAGA",
}

motif_classes = {
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


