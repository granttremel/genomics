#!/usr/bin/env python3
"""
te_parser.py - Lightweight TE subfamily name parser

This module parses TE family names into hierarchical components without
requiring access to the Dfam H5 database. Useful for working directly
with RepeatMasker output or Dfam hits files.

Usage:
    from te_parser import TEName, parse_te_name, build_te_index
    
    # Parse a single name
    te = parse_te_name("AluYk11")
    print(te.family, te.subfamily, te.variant)  # Alu Y k11
    
    # Build an index from a hits file
    index = build_te_index_from_hits("hg38_dfam.nrph.bed.gz")
"""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import gzip


@dataclass
class TEName:
    """Parsed transposable element name"""
    full_name: str           # Original full name
    family: str              # Major family (Alu, L1, MIR, etc.)
    subfamily: str           # First-level subfamily (Y, S, J for Alu; HS, PA2 for L1)
    variant: str             # Further subdivision (a5, k11, etc.)
    classification: str      # Repeat class if known (SINE/Alu, LINE/L1, etc.)
    suffix: str              # Any trailing info (_3end, _5end, etc.)
    
    @property
    def hierarchy_tuple(self) -> Tuple[str, ...]:
        """Return tuple for sorting/grouping: (family, subfamily, variant)"""
        return (self.family, self.subfamily, self.variant)
    
    @property
    def short_name(self) -> str:
        """Return family + subfamily (e.g., 'AluY', 'L1HS')"""
        if self.subfamily:
            return f"{self.family}{self.subfamily}"
        return self.family
    
    def matches_filter(self, family: str = None, subfamily: str = None) -> bool:
        """Check if this TE matches the given filter criteria"""
        if family and self.family != family:
            return False
        if subfamily and self.subfamily != subfamily:
            return False
        return True


# Classification lookup for common families
FAMILY_CLASSIFICATIONS = {
    'Alu': 'SINE/Alu',
    'MIR': 'SINE/MIR',
    'L1': 'LINE/L1',
    'L2': 'LINE/L2',
    'L3': 'LINE/L3',
    'CR1': 'LINE/CR1',
    'SVA': 'Retroposon/SVA',
    'HERV': 'LTR/ERV1',
    'LTR': 'LTR',
    'THE': 'LTR/ERV1',
    'MLT': 'LTR/ERVL-MaLR',
    'MST': 'LTR/ERVL-MaLR',
    'Charlie': 'DNA/hAT-Charlie',
    'Tigger': 'DNA/TcMar-Tigger',
    'MER': 'DNA',
    'MADE': 'DNA/TcMar-Mariner',
    'Mariner': 'DNA/TcMar-Mariner',
}


def parse_te_name(name: str) -> TEName:
    """
    Parse a TE name into its hierarchical components.
    
    Examples:
        parse_te_name("AluYk11")     -> TEName(family='Alu', subfamily='Y', variant='k11')
        parse_te_name("L1PA2")       -> TEName(family='L1', subfamily='PA2', variant='')
        parse_te_name("L1MC4a_3end") -> TEName(family='L1', subfamily='MC4a', variant='', suffix='_3end')
        parse_te_name("MIR3")        -> TEName(family='MIR', subfamily='3', variant='')
    """
    original = name
    suffix = ''
    
    # Extract common suffixes
    suffix_match = re.search(r'(_\d*end|_int|_I|_LTR|_orf\d*)$', name)
    if suffix_match:
        suffix = suffix_match.group(0)
        name = name[:suffix_match.start()]
    
    # Alu family: complex subfamily structure
    # AluY, AluYa5, AluYk11, AluSx, AluSp, AluJo, AluJb
    alu_match = re.match(r'^Alu([YSJC])([a-z]*)(\d*)$', name)
    if alu_match:
        subfamily = alu_match.group(1)
        variant = (alu_match.group(2) or '') + (alu_match.group(3) or '')
        return TEName(original, 'Alu', subfamily, variant, 'SINE/Alu', suffix)
    
    # L1 family: L1HS, L1PA2, L1P1, L1M1, L1MC4a, L1ME, L1MEg
    l1_match = re.match(r'^L1(HS|PA\d+|P[AB]?\d*|M[A-Z]*\d*[a-z]*|ME[a-z]?\d*)$', name)
    if l1_match:
        subfamily = l1_match.group(1)
        return TEName(original, 'L1', subfamily, '', 'LINE/L1', suffix)
    
    # L2, L3 families
    l2_match = re.match(r'^(L[23])([a-z]?)$', name)
    if l2_match:
        return TEName(original, l2_match.group(1), l2_match.group(2), '', f'LINE/{l2_match.group(1)}', suffix)
    
    # MIR family: MIR, MIRb, MIRc, MIR3
    mir_match = re.match(r'^MIR(\d?)([a-z]?)$', name)
    if mir_match:
        subfamily = mir_match.group(1) or mir_match.group(2) or ''
        variant = mir_match.group(2) if mir_match.group(1) else ''
        return TEName(original, 'MIR', subfamily, variant, 'SINE/MIR', suffix)
    
    # SVA family: SVA_A, SVA_B, etc.
    sva_match = re.match(r'^SVA_?([A-F]?)$', name)
    if sva_match:
        return TEName(original, 'SVA', sva_match.group(1), '', 'Retroposon/SVA', suffix)
    
    # HERV families: HERVK, HERV-K, HERVH, etc.
    herv_match = re.match(r'^HERV[-_]?([A-Z0-9]+)(.*)$', name)
    if herv_match:
        return TEName(original, 'HERV', herv_match.group(1), herv_match.group(2), 'LTR/ERV', suffix)
    
    # LTR elements: LTR7, LTR7b, LTR12C
    ltr_match = re.match(r'^(LTR|THE|MLT|MST)(\d+)([A-Za-z]*)$', name)
    if ltr_match:
        family = ltr_match.group(1)
        subfamily = ltr_match.group(2)
        variant = ltr_match.group(3)
        classification = FAMILY_CLASSIFICATIONS.get(family, 'LTR')
        return TEName(original, family, subfamily, variant, classification, suffix)
    
    # DNA transposons: Charlie1, Charlie15a, Tigger1, Tigger16a
    dna_match = re.match(r'^(Charlie|Tigger|MADE|Mariner)(\d+)([a-z]?)$', name)
    if dna_match:
        family = dna_match.group(1)
        subfamily = dna_match.group(2)
        variant = dna_match.group(3)
        classification = FAMILY_CLASSIFICATIONS.get(family, 'DNA')
        return TEName(original, family, subfamily, variant, classification, suffix)
    
    # MER elements: MER1A, MER5B, MER20
    mer_match = re.match(r'^MER(\d+)([A-Z]?)$', name)
    if mer_match:
        return TEName(original, 'MER', mer_match.group(1), mer_match.group(2), 'DNA', suffix)
    
    # Generic fallback: try to split name into alpha prefix and rest
    generic_match = re.match(r'^([A-Za-z]+)(\d*)([A-Za-z]*)$', name)
    if generic_match:
        family = generic_match.group(1)
        subfamily = generic_match.group(2) or ''
        variant = generic_match.group(3) or ''
        classification = FAMILY_CLASSIFICATIONS.get(family, '')
        return TEName(original, family, subfamily, variant, classification, suffix)
    
    # No pattern matched
    return TEName(original, name, '', '', '', suffix)


def build_te_index(names: List[str]) -> Dict[str, Dict[str, List[TEName]]]:
    """
    Build a hierarchical index from a list of TE names.
    
    Returns:
        Dict mapping family -> subfamily -> list of TEName objects
    
    Example:
        index = build_te_index(["AluY", "AluYa5", "AluYk11", "AluSx"])
        # index['Alu']['Y'] = [TEName(AluY), TEName(AluYa5), TEName(AluYk11)]
        # index['Alu']['S'] = [TEName(AluSx)]
    """
    index: Dict[str, Dict[str, List[TEName]]] = {}
    
    for name in names:
        te = parse_te_name(name)
        
        if te.family not in index:
            index[te.family] = {}
        
        subfamily_key = te.subfamily or '_root'
        if subfamily_key not in index[te.family]:
            index[te.family][subfamily_key] = []
        
        index[te.family][subfamily_key].append(te)
    
    return index


def build_te_index_from_hits(hits_file: str, name_col: int = 3) -> Dict[str, Dict[str, List[TEName]]]:
    """
    Build TE index from a Dfam/RepeatMasker hits file.
    
    Args:
        hits_file: Path to hits file (can be gzipped)
        name_col: 0-based column index for family name (default 3 for BED format)
    
    Returns:
        Hierarchical index of TE names
    """
    names = set()
    
    opener = gzip.open if hits_file.endswith('.gz') else open
    with opener(hits_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) > name_col:
                names.add(parts[name_col])
    
    return build_te_index(list(names))


def summarize_index(index: Dict[str, Dict[str, List[TEName]]]) -> None:
    """Print a summary of the TE index"""
    print(f"{'Family':<15} {'Subfamilies':<12} {'Total':<8} Classification")
    print("-" * 60)
    
    for family in sorted(index.keys()):
        subfamilies = index[family]
        total = sum(len(v) for v in subfamilies.values())
        
        # Get classification from first entry
        first_te = next(iter(subfamilies.values()))[0]
        classification = first_te.classification or 'Unknown'
        
        print(f"{family:<15} {len(subfamilies):<12} {total:<8} {classification}")


def get_subfamily_tree(index: Dict[str, Dict[str, List[TEName]]], family: str) -> None:
    """Print a tree view of subfamilies for a given family"""
    if family not in index:
        print(f"Family '{family}' not found")
        return
    
    print(f"{family}")
    subfamilies = index[family]
    
    for i, (subfamily, members) in enumerate(sorted(subfamilies.items())):
        is_last = i == len(subfamilies) - 1
        prefix = "└── " if is_last else "├── "
        child_prefix = "    " if is_last else "│   "
        
        print(f"{prefix}{subfamily or '(base)'} ({len(members)} members)")
        
        # Show first few variants
        variants = sorted(set(m.variant for m in members if m.variant))
        if variants:
            var_str = ", ".join(variants[:10])
            if len(variants) > 10:
                var_str += f", ... +{len(variants)-10} more"
            print(f"{child_prefix}variants: {var_str}")


# ============================================================================
# Command-line interface
# ============================================================================

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("""
TE Name Parser - Parse transposable element family names

Usage:
    python te_parser.py <command> [args]

Commands:
    parse <name>       - Parse a single TE name
    demo               - Show parsing examples
    index <file>       - Build index from hits file (BED/TSV)
    tree <file> <fam>  - Show subfamily tree for a family

Examples:
    python te_parser.py parse AluYk11
    python te_parser.py demo
    python te_parser.py index hg38_dfam.nrph.bed.gz
    python te_parser.py tree hg38_dfam.nrph.bed.gz Alu
        """)
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "parse":
        if len(sys.argv) < 3:
            print("Usage: te_parser.py parse <name>")
            sys.exit(1)
        te = parse_te_name(sys.argv[2])
        print(f"Full name:       {te.full_name}")
        print(f"Family:          {te.family}")
        print(f"Subfamily:       {te.subfamily or '(none)'}")
        print(f"Variant:         {te.variant or '(none)'}")
        print(f"Suffix:          {te.suffix or '(none)'}")
        print(f"Classification:  {te.classification or '(unknown)'}")
        print(f"Hierarchy tuple: {te.hierarchy_tuple}")
    
    elif command == "demo":
        examples = [
            "AluY", "AluYa5", "AluYk11", "AluSx", "AluJo", "AluJb",
            "L1HS", "L1PA2", "L1PA3", "L1P1", "L1M1", "L1MC4a", "L1ME", "L1MEg",
            "L1MC4a_3end", "L1PA2_5end",
            "MIR", "MIRb", "MIRc", "MIR3",
            "SVA_A", "SVA_B", "SVA_F",
            "HERVK", "HERV-K", "HERVH",
            "LTR7", "LTR7b", "LTR12C",
            "THE1A", "THE1B", "MLT1A", "MST",
            "Charlie1", "Charlie15a", "Tigger1", "Tigger16a",
            "MER1A", "MER5B", "MER20",
            "L2", "L2a", "L2b", "L3", "L3b",
        ]
        
        print(f"{'Name':<20} {'Family':<10} {'Subfam':<8} {'Var':<6} {'Classification'}")
        print("-" * 70)
        for name in examples:
            te = parse_te_name(name)
            print(f"{name:<20} {te.family:<10} {te.subfamily:<8} {te.variant:<6} {te.classification}")
    
    elif command == "index":
        if len(sys.argv) < 3:
            print("Usage: te_parser.py index <hits_file>")
            sys.exit(1)
        
        print(f"Building index from {sys.argv[2]}...")
        index = build_te_index_from_hits(sys.argv[2])
        print()
        summarize_index(index)
    
    elif command == "tree":
        if len(sys.argv) < 4:
            print("Usage: te_parser.py tree <hits_file> <family>")
            sys.exit(1)
        
        index = build_te_index_from_hits(sys.argv[2])
        get_subfamily_tree(index, sys.argv[3])
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
