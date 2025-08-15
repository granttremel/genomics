"""Advanced regex patterns for genomics.

These patterns go beyond simple matching to capture complex biological patterns.
"""

import re
import regex
from typing import List, Tuple, Optional

class GenomicRegexPatterns:
    """Collection of useful regex patterns for genomic analysis."""
    
    # === BASIC PATTERNS ===
    
    # DNA/RNA bases with IUPAC ambiguity codes
    IUPAC_DNA = {
        'R': '[AG]',   # puRine
        'Y': '[CT]',   # pYrimidine  
        'S': '[GC]',   # Strong (3 H bonds)
        'W': '[AT]',   # Weak (2 H bonds)
        'K': '[GT]',   # Keto
        'M': '[AC]',   # aMino
        'B': '[CGT]',  # not A
        'D': '[AGT]',  # not C
        'H': '[ACT]',  # not G
        'V': '[ACG]',  # not T
        'N': '[ACGT]', # aNy
    }
    
    # === REPEAT PATTERNS ===
    
    @staticmethod
    def tandem_repeat(unit: str, min_repeats: int = 2, max_repeats: int = None) -> str:
        """Create pattern for tandem repeats.
        
        Examples:
            'CAG' repeats (Huntington's): (CAG){10,}
            Microsatellites: (CA){5,20}
        """
        if max_repeats:
            return f"({unit}){{{min_repeats},{max_repeats}}}"
        return f"({unit}){{{min_repeats},}}"
    
    @staticmethod
    def fuzzy_repeat(unit: str, mismatches: int = 1) -> str:
        """Allow mismatches in repeat units.
        
        Example: CA repeat with 1 mismatch allowed
        """
        if len(unit) == 2 and mismatches == 1:
            # For dinucleotide with 1 mismatch
            return f"({unit}|{unit[0]}.|.{unit[1]})*"
        
        # General case - build alternation
        patterns = [unit]
        for i in range(len(unit)):
            pattern = list(unit)
            pattern[i] = '.'
            patterns.append(''.join(pattern))
        
        return f"({'|'.join(patterns)})*"
    
    @staticmethod
    def interrupted_repeat(unit: str, max_gap: int = 3) -> str:
        """Repeats with interruptions.
        
        Example: CAG repeats with up to 3bp interruptions
        """
        return f"(({unit}){{2,}}(.{{0,{max_gap}}}({unit}){{2,}})*)"
    
    # === STRUCTURAL PATTERNS ===
    
    @staticmethod
    def g_quadruplex(min_g_run: int = 3, max_loop: int = 7) -> str:
        """G-quadruplex forming sequences.
        
        Pattern: G3+N1-7G3+N1-7G3+N1-7G3+
        """
        g_run = f"G{{{min_g_run},}}"
        loop = f".{{1,{max_loop}}}"
        return f"({g_run}{loop}){{3}}{g_run}"
    
    @staticmethod
    def inverted_repeat(min_stem: int = 6, max_loop: int = 20) -> str:
        """Palindromic sequences that can form hairpins.
        
        Uses backreferences to match reverse complement
        """
        # This is complex - we need to capture and reverse complement
        # Simplified version that finds potential regions
        return f"([ACGT]{{{min_stem}}}).{{0,{max_loop}}}(?=(.*)).*"
    
    @staticmethod
    def z_dna_region() -> str:
        """Regions likely to form Z-DNA (alternating purines/pyrimidines)."""
        return r"([AG][CT]){5,}"
    
    # === REGULATORY PATTERNS ===
    
    @staticmethod
    def tata_box() -> str:
        """TATA box and variants."""
        return r"TATA[AT]A[AT][AG]"
    
    @staticmethod
    def caat_box() -> str:
        """CAAT box promoter element."""
        return r"GG[CT]CAATCT"
    
    @staticmethod
    def gc_box() -> str:
        """GC box (Sp1 binding site)."""
        return r"GGGCGG|CCGCCC"
    
    @staticmethod
    def e_box() -> str:
        """E-box (bHLH binding site)."""
        return r"CA[ACGT]{2}TG"
    
    @staticmethod
    def kozak_sequence() -> str:
        """Kozak consensus for translation initiation."""
        return r"[AG]..ATGG"
    
    # === SPLICE SITES ===
    
    @staticmethod
    def splice_donor() -> str:
        """5' splice site (donor)."""
        return r"[AC]AG\|GT[AG]AGT"  # | marks exon/intron boundary
    
    @staticmethod
    def splice_acceptor() -> str:
        """3' splice site (acceptor)."""
        return r"[CT]{10,}[ATGC]CAG\|G"
    
    @staticmethod
    def branch_point() -> str:
        """Branch point for splicing."""
        return r"[CT]T[AG]A[CT]"
    
    # === EPIGENETIC PATTERNS ===
    
    @staticmethod
    def cpg_island(min_length: int = 200) -> str:
        """CpG islands (unmethylated regions)."""
        # Simplified - real CpG island detection needs statistics
        return f"(CG.*?){{10,}}"
    
    @staticmethod
    def methylation_site() -> str:
        """CpG dinucleotides (methylation targets)."""
        return r"CG"
    
    # === POLYMORPHISM PATTERNS ===
    
    @staticmethod
    def variable_tandem_repeat(units: List[str]) -> str:
        """VNTRs with variable units.
        
        Example: ['CAG', 'CAA'] for polyglutamine
        """
        unit_pattern = '|'.join(units)
        return f"({unit_pattern}){{3,}}"
    
    @staticmethod
    def snp_context(ref_base: str, alt_base: str, context: int = 2) -> str:
        """SNP with surrounding context.
        
        Captures context around a known SNP
        """
        return f"([ACGT]{{{context}}})[{ref_base}{alt_base}]([ACGT]{{{context}}})"
    
    # === RESTRICTION SITES ===
    
    RESTRICTION_ENZYMES = {
        'EcoRI': r'GAATTC',
        'BamHI': r'GGATCC',
        'HindIII': r'AAGCTT',
        'PstI': r'CTGCAG',
        'SmaI': r'CCCGGG',
        'XbaI': r'TCTAGA',
        'NotI': r'GCGGCCGC',
        'SfiI': r'GGCC.{5}GGCC',  # Has variable middle region
    }
    
    # === COMPLEX PATTERNS ===
    
    @staticmethod
    def microsatellite_instability() -> str:
        """Detect microsatellite regions prone to instability."""
        patterns = [
            r"(A){8,}",           # Mononucleotide
            r"([AC]T){5,}",       # Dinucleotide
            r"(CAG|CTG){4,}",     # Trinucleotide (disease-associated)
            r"([ACGT]{4}){3,}",   # Tetranucleotide
        ]
        return '|'.join(patterns)
    
    @staticmethod
    def transposon_terminal_repeat() -> str:
        """Terminal inverted repeats of transposons."""
        # Simplified - real TIRs need palindrome matching
        return r"([ACGT]{10,50}).*\1"
    
    @staticmethod
    def polya_signal() -> str:
        """Polyadenylation signals."""
        return r"AATAAA|ATTAAA|AGTAAA"
    
    @staticmethod
    def shine_dalgarno() -> str:
        """Ribosome binding site in prokaryotes."""
        return r"AGGAGG.{5,9}ATG"
    
    # === UTILITY FUNCTIONS ===
    
    @staticmethod
    def find_all_overlapping(pattern: str, sequence: str) -> List[Tuple[int, int, str]]:
        """Find all matches including overlapping ones.
        
        Standard re.findall misses overlaps.
        """
        matches = []
        regex = re.compile(pattern)
        
        pos = 0
        while pos < len(sequence):
            match = regex.search(sequence, pos)
            if not match:
                break
            matches.append((match.start(), match.end(), match.group()))
            pos = match.start() + 1  # Move just 1 position to find overlaps
        
        return matches
    
    @staticmethod
    def find_imperfect_matches(pattern: str, sequence: str, max_mismatches: int = 2) -> List:
        """Find approximate matches allowing mismatches.
        
        Uses regex module (install with: pip install regex)
        """
        try:
            fuzzy_pattern = f"({pattern}){{e<={max_mismatches}}}"
            return regex.findall(fuzzy_pattern, sequence)
        except ImportError:
            print("Install 'regex' module for fuzzy matching: pip install regex")
            return []
    
    @staticmethod
    def expand_iupac(pattern: str) -> str:
        """Expand IUPAC codes in pattern to regex."""
        expanded = pattern
        for code, expansion in GenomicRegexPatterns.IUPAC_DNA.items():
            expanded = expanded.replace(code, expansion)
        return expanded


# === EXAMPLES ===

def demonstrate_patterns():
    """Show example usage of genomic regex patterns."""
    
    patterns = GenomicRegexPatterns()
    
    # Example sequence with various features
    test_seq = "GAATTCGGGGGGATGCAGCAGCAGCAGCAGCAGCAGAATAAATATAAAAAAAAAAACGTAGT"
    
    print("Testing genomic regex patterns:\n")
    
    # Find restriction sites
    print("EcoRI sites:")
    for site, pattern in patterns.RESTRICTION_ENZYMES.items():
        matches = re.finditer(pattern, test_seq)
        for match in matches:
            print(f"  {site} at position {match.start()}: {match.group()}")
    
    # Find repeats
    print("\nCAG repeats (Huntington's disease):")
    cag_pattern = patterns.tandem_repeat("CAG", min_repeats=3)
    matches = re.finditer(cag_pattern, test_seq)
    for match in matches:
        repeat_count = len(match.group()) // 3
        print(f"  Found {repeat_count} CAG repeats at position {match.start()}")
    
    # Find polyA signal
    print("\nPolyA signals:")
    polya_pattern = patterns.polya_signal()
    matches = re.finditer(polya_pattern, test_seq)
    for match in matches:
        print(f"  PolyA signal at position {match.start()}: {match.group()}")
    
    # Find G-rich regions (potential G-quadruplex)
    print("\nG-rich regions:")
    g_pattern = patterns.g_quadruplex(min_g_run=4, max_loop=3)
    matches = re.finditer(g_pattern, test_seq)
    for match in matches:
        print(f"  Potential G-quadruplex at position {match.start()}")
    
    # Find homopolymer runs
    print("\nHomopolymer runs:")
    for base in "ACGT":
        pattern = f"{base}{{5,}}"
        matches = re.finditer(pattern, test_seq)
        for match in matches:
            print(f"  {base} run of length {len(match.group())} at position {match.start()}")

if __name__ == "__main__":
    demonstrate_patterns()