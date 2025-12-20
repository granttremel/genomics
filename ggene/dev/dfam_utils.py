#!/usr/bin/env python3
"""
dfam_utils.py - Utilities for working with Dfam/FamDB data

This module provides:
1. H5 file querying for family metadata
2. Subfamily name parsing for common TE families (Alu, L1, etc.)
3. Classification hierarchy building

Usage:
    from dfam_utils import DfamDB, parse_subfamily_name
    
    db = DfamDB("/path/to/dfam38-1_full.0.h5")
    info = db.get_family("AluY")
    
    # Parse subfamily names into structured tuples
    parsed = parse_subfamily_name("AluYk11")  
    # -> {'family': 'Alu', 'subfamily': 'Y', 'subsubfamily': 'k11', 'full': 'AluYk11'}
"""

import h5py
import json
import re
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class TEFamily:
    """Represents a transposable element family with metadata"""
    accession: str
    name: str
    classification: str = ""
    repeat_type: str = ""
    repeat_subtype: str = ""
    length: int = 0
    description: str = ""
    clades: List[int] = field(default_factory=list)
    
    # Parsed subfamily info
    family_group: str = ""      # e.g., "Alu", "L1", "MIR"
    subfamily: str = ""          # e.g., "Y", "S", "J" for Alu
    subsubfamily: str = ""       # e.g., "k11", "a2" 
    
    def to_dict(self) -> Dict:
        return asdict(self)


class DfamDB:
    """
    Interface to FamDB H5 files for querying repeat family information.
    
    Example:
        db = DfamDB("/path/to/dfam38-1_full.0.h5")
        
        # Get single family
        alu = db.get_family("AluY")
        
        # Get all Alu families
        alus = db.get_families_by_prefix("Alu")
        
        # Build classification lookup
        lookup = db.build_classification_lookup()
    """
    
    def __init__(self, h5_path: str):
        self.h5_path = Path(h5_path)
        if not self.h5_path.exists():
            raise FileNotFoundError(f"FamDB file not found: {h5_path}")
        
        # Cache for family data
        self._name_to_acc: Dict[str, str] = {}
        self._families_cache: Dict[str, TEFamily] = {}
        self._index_built = False
    
    def _build_name_index(self):
        """Build index mapping family names to accessions"""
        if self._index_built:
            return
        
        with h5py.File(self.h5_path, 'r') as f:
            if 'Families' not in f:
                raise ValueError("No Families group in H5 file")
            
            families = f['Families']
            for acc in families.keys():
                fam = families[acc]
                if 'name' in fam.attrs:
                    name = fam.attrs['name']
                    if isinstance(name, bytes):
                        name = name.decode('utf-8')
                    self._name_to_acc[name] = acc
        
        self._index_built = True
        print(f"Indexed {len(self._name_to_acc)} families")
    
    def _extract_family_info(self, acc: str, fam_group: h5py.Group) -> TEFamily:
        """Extract TEFamily from H5 group"""
        info = TEFamily(accession=acc, name="")
        
        # Extract attributes
        attr_map = {
            'name': 'name',
            'classification': 'classification', 
            'repeat_type': 'repeat_type',
            'repeat_subtype': 'repeat_subtype',
            'description': 'description',
            'length': 'length',
        }
        
        for h5_attr, field_name in attr_map.items():
            if h5_attr in fam_group.attrs:
                val = fam_group.attrs[h5_attr]
                if isinstance(val, bytes):
                    val = val.decode('utf-8', errors='replace')
                elif hasattr(val, 'item'):  # numpy scalar
                    val = val.item()
                setattr(info, field_name, val)
        
        # Handle clades specially (may be array)
        if 'clades' in fam_group.attrs:
            clades = fam_group.attrs['clades']
            if hasattr(clades, 'tolist'):
                info.clades = clades.tolist()
            else:
                info.clades = [clades] if clades else []
        
        # Parse subfamily structure
        parsed = parse_subfamily_name(info.name)
        info.family_group = parsed.get('family', '')
        info.subfamily = parsed.get('subfamily', '')
        info.subsubfamily = parsed.get('subsubfamily', '')
        
        return info
    
    def get_family(self, name_or_acc: str) -> Optional[TEFamily]:
        """Get family info by name or accession"""
        self._build_name_index()
        
        # Check cache first
        if name_or_acc in self._families_cache:
            return self._families_cache[name_or_acc]
        
        # Determine accession
        if name_or_acc.startswith('DF') or name_or_acc.startswith('DR'):
            acc = name_or_acc
        else:
            acc = self._name_to_acc.get(name_or_acc)
            if acc is None:
                return None
        
        with h5py.File(self.h5_path, 'r') as f:
            families = f['Families']
            if acc not in families:
                return None
            
            info = self._extract_family_info(acc, families[acc])
            self._families_cache[name_or_acc] = info
            self._families_cache[acc] = info
            return info
    
    def get_families_by_prefix(self, prefix: str) -> List[TEFamily]:
        """Get all families whose names start with a prefix (e.g., 'Alu')"""
        self._build_name_index()
        
        matching = [name for name in self._name_to_acc if name.startswith(prefix)]
        return [self.get_family(name) for name in matching]
    
    def get_families_by_classification(self, classification: str) -> List[TEFamily]:
        """Get all families matching a classification string (e.g., 'SINE/Alu')"""
        results = []
        
        with h5py.File(self.h5_path, 'r') as f:
            families = f['Families']
            for acc in families.keys():
                fam = families[acc]
                if 'classification' in fam.attrs:
                    cls = fam.attrs['classification']
                    if isinstance(cls, bytes):
                        cls = cls.decode('utf-8')
                    if classification in cls:
                        results.append(self._extract_family_info(acc, fam))
        
        return results
    
    def get_all_families(self) -> List[TEFamily]:
        """Get all families (may be slow for large databases)"""
        results = []
        
        with h5py.File(self.h5_path, 'r') as f:
            families = f['Families']
            total = len(families.keys())
            
            for i, acc in enumerate(families.keys()):
                if i % 500 == 0:
                    print(f"Loading families: {i}/{total}")
                results.append(self._extract_family_info(acc, families[acc]))
        
        return results
    
    def build_classification_lookup(self) -> Dict[str, Dict]:
        """
        Build a hierarchical lookup structure for all families.
        
        Returns dict like:
        {
            'Alu': {
                'classification': 'SINE/Alu',
                'subfamilies': {
                    'Y': ['AluY', 'AluYa5', 'AluYb8', 'AluYk11', ...],
                    'S': ['AluSx', 'AluSp', ...],
                    'J': ['AluJo', 'AluJb', ...],
                }
            },
            'L1': {
                'classification': 'LINE/L1',
                'subfamilies': {...}
            }
        }
        """
        lookup = {}
        
        for fam in self.get_all_families():
            family_group = fam.family_group or fam.name
            subfamily = fam.subfamily or '_root'
            
            if family_group not in lookup:
                lookup[family_group] = {
                    'classification': fam.classification,
                    'subfamilies': {}
                }
            
            if subfamily not in lookup[family_group]['subfamilies']:
                lookup[family_group]['subfamilies'][subfamily] = []
            
            lookup[family_group]['subfamilies'][subfamily].append(fam.name)
        
        return lookup
    
    def list_family_names(self) -> List[str]:
        """Get list of all family names"""
        self._build_name_index()
        return list(self._name_to_acc.keys())


# ============================================================================
# Subfamily name parsing
# ============================================================================

# Patterns for common TE families
SUBFAMILY_PATTERNS = {
    # Alu family: AluY, AluSx, AluJo, AluYa5, AluYk11, etc.
    'Alu': re.compile(r'^Alu([YSJC]?)([a-z]*)(\d*)(.*)$'),
    
    # L1 family: L1HS, L1PA2, L1M1, L1MC4a, L1MEg, etc.
    'L1': re.compile(r'^L1(HS|PA\d+|P\d*|M[A-Z]*\d*[a-z]*|ME[a-z]*\d*)(.*)$'),
    
    # MIR family: MIR, MIRb, MIRc, MIR1_Amn, etc.  
    'MIR': re.compile(r'^MIR(\d*)([a-z]?)(_.*)?$'),
    
    # SVA family: SVA_A, SVA_B, SVA_C, etc.
    'SVA': re.compile(r'^SVA_?([A-F]?)(.*)$'),
    
    # HERV families: HERV-K, HERVH, etc.
    'HERV': re.compile(r'^HERV[-_]?([A-Z0-9]+)(.*)$'),
    
    # LTR elements: LTR7, LTR12, THE1, etc.
    'LTR': re.compile(r'^(LTR|THE|MLT|MST)(\d+)([A-Za-z]*)(.*)$'),
    
    # DNA transposons: Charlie, Tigger, etc.
    'Charlie': re.compile(r'^Charlie(\d+)([a-z]?)(.*)$'),
    'Tigger': re.compile(r'^Tigger(\d+)([a-z]?)(.*)$'),
    'MER': re.compile(r'^MER(\d+)([A-Za-z]?)(.*)$'),
}


def parse_subfamily_name(name: str) -> Dict[str, str]:
    """
    Parse a TE family name into its component parts.
    
    Examples:
        parse_subfamily_name("AluYk11") 
        -> {'family': 'Alu', 'subfamily': 'Y', 'subsubfamily': 'k11', 'full': 'AluYk11'}
        
        parse_subfamily_name("L1PA2")
        -> {'family': 'L1', 'subfamily': 'PA2', 'subsubfamily': '', 'full': 'L1PA2'}
        
        parse_subfamily_name("MIRb")
        -> {'family': 'MIR', 'subfamily': 'b', 'subsubfamily': '', 'full': 'MIRb'}
    
    Args:
        name: The TE family name (e.g., "AluYa5", "L1HS", "MIR3")
    
    Returns:
        Dict with keys: family, subfamily, subsubfamily, full, extra
    """
    result = {
        'family': '',
        'subfamily': '',
        'subsubfamily': '',
        'extra': '',
        'full': name
    }
    
    # Try each pattern
    for family_prefix, pattern in SUBFAMILY_PATTERNS.items():
        if name.startswith(family_prefix) or (family_prefix == 'LTR' and 
            any(name.startswith(p) for p in ['LTR', 'THE', 'MLT', 'MST'])):
            match = pattern.match(name)
            if match:
                groups = match.groups()
                result['family'] = family_prefix
                
                if family_prefix == 'Alu':
                    # Alu special handling: AluY, AluYa5, AluYk11
                    result['subfamily'] = groups[0] if groups[0] else ''
                    result['subsubfamily'] = (groups[1] or '') + (groups[2] or '')
                    result['extra'] = groups[3] if len(groups) > 3 else ''
                    
                elif family_prefix == 'L1':
                    # L1 special handling
                    result['subfamily'] = groups[0] if groups[0] else ''
                    result['extra'] = groups[1] if len(groups) > 1 else ''
                    
                elif family_prefix in ['MIR', 'SVA', 'HERV']:
                    result['subfamily'] = groups[0] if groups[0] else ''
                    result['subsubfamily'] = groups[1] if len(groups) > 1 and groups[1] else ''
                    
                elif family_prefix == 'LTR':
                    # LTR/THE/MLT/MST handling
                    result['family'] = groups[0]
                    result['subfamily'] = groups[1] if groups[1] else ''
                    result['subsubfamily'] = groups[2] if len(groups) > 2 and groups[2] else ''
                    
                else:
                    # Generic handling for Charlie, Tigger, MER
                    result['subfamily'] = groups[0] if groups[0] else ''
                    result['subsubfamily'] = groups[1] if len(groups) > 1 and groups[1] else ''
                
                return result
    
    # No pattern matched - try generic parsing
    # Look for common patterns: Name + Number + Letter
    generic = re.match(r'^([A-Za-z]+)(\d*)([A-Za-z]*)(_.*)?$', name)
    if generic:
        result['family'] = generic.group(1)
        result['subfamily'] = generic.group(2) or ''
        result['subsubfamily'] = generic.group(3) or ''
        result['extra'] = generic.group(4) or ''
    else:
        result['family'] = name
    
    return result


def build_subfamily_hierarchy(names: List[str]) -> Dict:
    """
    Build a hierarchical structure from a list of TE family names.
    
    Args:
        names: List of TE names like ["AluY", "AluYa5", "AluSx", "L1HS", ...]
    
    Returns:
        Nested dict structure organizing by family -> subfamily -> members
    """
    hierarchy = {}
    
    for name in names:
        parsed = parse_subfamily_name(name)
        family = parsed['family'] or 'Unknown'
        subfamily = parsed['subfamily'] or '_root'
        
        if family not in hierarchy:
            hierarchy[family] = {}
        
        if subfamily not in hierarchy[family]:
            hierarchy[family][subfamily] = []
        
        hierarchy[family][subfamily].append({
            'name': name,
            'subsubfamily': parsed['subsubfamily'],
            'extra': parsed['extra']
        })
    
    return hierarchy


# ============================================================================
# Command-line interface
# ============================================================================

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("""
Dfam Utilities - Work with FamDB H5 files

Usage:
    python dfam_utils.py <h5_file> <command> [args]

Commands:
    list              - List all family names
    info <name>       - Get info for a family
    search <prefix>   - Find families starting with prefix
    parse <name>      - Parse a family name into components
    hierarchy <file>  - Build hierarchy from list of names (one per line)
    export <outfile>  - Export all families to JSON

Examples:
    python dfam_utils.py dfam38-1_full.0.h5 info AluYa5
    python dfam_utils.py dfam38-1_full.0.h5 search Alu
    python dfam_utils.py - parse AluYk11
        """)
        sys.exit(1)
    
    h5_file = sys.argv[1]
    command = sys.argv[2] if len(sys.argv) > 2 else "list"
    
    # Parse command doesn't need H5 file
    if command == "parse":
        name = sys.argv[3] if len(sys.argv) > 3 else ""
        if name:
            result = parse_subfamily_name(name)
            print(json.dumps(result, indent=2))
        else:
            # Demo with some examples
            examples = ["AluY", "AluYa5", "AluYk11", "AluSx", "AluJo",
                       "L1HS", "L1PA2", "L1M1", "L1MC4a",
                       "MIR", "MIRb", "MIR3", 
                       "SVA_A", "HERVK", "LTR7", "THE1A",
                       "Charlie15a", "Tigger16a", "MER5B"]
            for ex in examples:
                result = parse_subfamily_name(ex)
                print(f"{ex:20s} -> family={result['family']:8s} sub={result['subfamily']:5s} subsub={result['subsubfamily']}")
        sys.exit(0)
    
    # Other commands need the H5 file
    if h5_file == "-":
        print("H5 file required for this command")
        sys.exit(1)
    
    db = DfamDB(h5_file)
    
    if command == "list":
        names = db.list_family_names()
        for name in sorted(names):
            print(name)
    
    elif command == "info":
        name = sys.argv[3] if len(sys.argv) > 3 else ""
        if not name:
            print("Usage: dfam_utils.py <h5> info <family_name>")
            sys.exit(1)
        info = db.get_family(name)
        if info:
            print(json.dumps(info.to_dict(), indent=2))
        else:
            print(f"Family not found: {name}")
    
    elif command == "search":
        prefix = sys.argv[3] if len(sys.argv) > 3 else ""
        if not prefix:
            print("Usage: dfam_utils.py <h5> search <prefix>")
            sys.exit(1)
        families = db.get_families_by_prefix(prefix)
        print(f"Found {len(families)} families starting with '{prefix}':")
        for fam in sorted(families, key=lambda x: x.name):
            print(f"  {fam.name:25s} {fam.classification}")
    
    elif command == "export":
        outfile = sys.argv[3] if len(sys.argv) > 3 else "families.json"
        families = db.get_all_families()
        with open(outfile, 'w') as f:
            json.dump([fam.to_dict() for fam in families], f, indent=2)
        print(f"Exported {len(families)} families to {outfile}")
    
    elif command == "hierarchy":
        # Build hierarchy from names in H5 or from a file
        if len(sys.argv) > 3:
            with open(sys.argv[3]) as f:
                names = [line.strip() for line in f if line.strip()]
        else:
            names = db.list_family_names()
        
        hier = build_subfamily_hierarchy(names)
        print(json.dumps(hier, indent=2))
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
