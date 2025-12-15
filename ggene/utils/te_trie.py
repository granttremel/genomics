#!/usr/bin/env python3
"""
te_trie.py - Discover TE family/subfamily structure from names alone

This module builds a prefix tree (trie) from TE names and compresses it
to reveal natural family/subfamily boundaries without any prior knowledge
of TE nomenclature.

The algorithm:
1. Build a character-level trie from all names
2. Collapse linear chains (single-child nodes) into single nodes
3. Identify branching points as family/subfamily boundaries
4. Apply minimum support threshold to filter noise

Usage:
    from te_trie import TETrieIndex
    
    names = ["AluY", "AluYa5", "AluYk11", "AluSx", "AluJo", "L1HS", "L1PA2"]
    index = TETrieIndex(names, min_support=2)
    
    # Get the discovered hierarchy
    hierarchy = index.get_hierarchy()
    
    # Parse a name according to discovered structure
    parts = index.parse_name("AluYa5")  # -> ['Alu', 'Y', 'a5']
"""

from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field
from collections import defaultdict
import json


@dataclass
class TrieNode:
    """Node in the character trie"""
    char: str = ""
    children: Dict[str, 'TrieNode'] = field(default_factory=dict)
    is_terminal: bool = False  # True if a full name ends here
    count: int = 0             # Number of names passing through this node
    names: Set[str] = field(default_factory=set)  # Names terminating here
    
    def is_branch_point(self) -> bool:
        """True if this node has multiple children"""
        return len(self.children) > 1
    
    def is_linear(self) -> bool:
        """True if this node has exactly one child"""
        return len(self.children) == 1
    
    def single_child(self) -> Optional['TrieNode']:
        """Return the single child if linear, else None"""
        if len(self.children) == 1:
            return next(iter(self.children.values()))
        return None


@dataclass 
class CompressedNode:
    """Node in the compressed trie where linear chains are collapsed"""
    label: str                  # The collapsed string (e.g., "Alu" instead of "A"->"l"->"u")
    children: Dict[str, 'CompressedNode'] = field(default_factory=dict)
    is_terminal: bool = False
    count: int = 0
    terminal_names: Set[str] = field(default_factory=set)
    
    def to_dict(self) -> dict:
        """Convert to nested dict for JSON serialization"""
        result = {
            'label': self.label,
            'count': self.count,
            'is_terminal': self.is_terminal,
        }
        if self.terminal_names:
            result['names'] = sorted(self.terminal_names)
        if self.children:
            result['children'] = {k: v.to_dict() for k, v in sorted(self.children.items())}
        return result


class TETrieIndex:
    """
    Discovers TE family/subfamily structure from a list of names.
    
    Example:
        names = ["AluY", "AluYa5", "AluYk11", "AluSx", "AluJo"]
        index = TETrieIndex(names, min_support=2)
        
        # Discovered structure:
        # Alu (root, 5 names)
        #   ├── Y (3 names)
        #   │   ├── a5 (terminal)
        #   │   └── k11 (terminal)
        #   ├── Sx (terminal)
        #   └── Jo (terminal)
    """
    
    def __init__(self, names: List[str], min_support: int = 2, 
                 branch_threshold: float = 0.3):
        """
        Build the index from a list of TE names.
        
        Args:
            names: List of TE family names
            min_support: Minimum count for a branch to be considered significant.
                         Branches with fewer than this many descendants are 
                         collapsed into their parent.
            branch_threshold: Minimum fraction of parent's count that a branch
                             must have to be considered a real split (default 0.3).
                             This prevents rare prefixes from creating spurious families.
        """
        self.names = list(set(names))  # Deduplicate
        self.min_support = min_support
        self.branch_threshold = branch_threshold
        
        # Build character-level trie
        self._root = TrieNode()
        for name in self.names:
            self._insert(name)
        
        # Compress the trie with smart boundary detection
        self._compressed_root = self._compress_smart(self._root, "")
        
        # Promote children of stub nodes (like "L" -> "L1", "LTR")
        self._promote_children()
        
        # Apply minimum support filtering
        if min_support > 1:
            self._filter_by_support(self._compressed_root)
        
        # Build lookup for fast parsing
        self._build_parse_lookup()
    
    def _insert(self, name: str) -> None:
        """Insert a name into the character trie"""
        node = self._root
        node.count += 1
        
        for char in name:
            if char not in node.children:
                node.children[char] = TrieNode(char=char)
            node = node.children[char]
            node.count += 1
        
        node.is_terminal = True
        node.names.add(name)
    
    def _compress(self, node: TrieNode, accumulated: str) -> CompressedNode:
        """
        Compress the trie by collapsing linear chains.
        
        A linear chain like A -> l -> u becomes a single node "Alu".
        Branching points are preserved as separate nodes.
        """
        # Accumulate characters while we're in a linear chain
        current = node
        label = accumulated
        
        while current.is_linear() and not current.is_terminal:
            child = current.single_child()
            label += child.char
            current = child
        
        # Now current is either: terminal, branch point, or leaf
        compressed = CompressedNode(
            label=label,
            is_terminal=current.is_terminal,
            count=current.count,
            terminal_names=current.names.copy()
        )
        
        # Recursively compress children
        for char, child in current.children.items():
            child_compressed = self._compress(child, char)
            compressed.children[child_compressed.label] = child_compressed
        
        return compressed
    
    def _compress_smart(self, node: TrieNode, accumulated: str) -> CompressedNode:
        """
        Smart compression that identifies natural family boundaries.
        
        Uses multiple heuristics:
        1. Branching balance - real splits have multiple significant branches  
        2. Terminal density - names ending here signal a natural grouping
        3. Character transitions - some transitions indicate word boundaries
        """
        current = node
        label = accumulated
        
        def is_boundary_char_transition(prev_char: str, next_char: str) -> bool:
            """Detect character transitions that suggest word boundaries"""
            if not prev_char or not next_char:
                return False
            # Number after letters: split at Alu|1 but keep L|1 together if L1 is the family
            if prev_char.isalpha() and next_char.isdigit():
                return True
            # Underscore is always a boundary
            if next_char == '_' or prev_char == '_':
                return True
            return False
        
        # Accumulate through the trie, stopping at natural boundaries
        while True:
            # Always stop at terminals - a name ends here
            if current.is_terminal:
                break
            
            # No children = leaf  
            if len(current.children) == 0:
                break
            
            # Single child - usually continue, but check for char transitions
            if len(current.children) == 1:
                child = list(current.children.values())[0]
                
                # Check if this is a natural boundary based on character type
                prev_char = label[-1] if label else ''
                next_char = child.char
                
                # If this looks like a word boundary AND the child has
                # multiple branches, stop here
                if is_boundary_char_transition(prev_char, next_char):
                    if len(child.children) > 1 or child.is_terminal:
                        break
                
                label += child.char
                current = child
                continue
            
            # Multiple children - is this a real branch point?
            children = list(current.children.values())
            counts = sorted([c.count for c in children], reverse=True)
            total = sum(counts)
            
            # Calculate branching metrics
            top_fraction = counts[0] / total if total > 0 else 1
            significant_branches = sum(1 for c in counts 
                                       if c >= max(self.min_support, 2))
            
            # Real branch point conditions:
            if significant_branches >= 2 and top_fraction < 0.85:
                break
            elif significant_branches >= 3:
                break
            else:
                # Follow the dominant branch
                dominant = max(children, key=lambda c: c.count)
                label += dominant.char
                current = dominant
        
        compressed = CompressedNode(
            label=label,
            is_terminal=current.is_terminal,
            count=current.count,
            terminal_names=current.names.copy()
        )
        
        # Recursively compress children
        for char, child in current.children.items():
            child_compressed = self._compress_smart(child, char)
            compressed.children[child_compressed.label] = child_compressed
        
        return compressed
    
    def _should_split_to_root(self, node: CompressedNode) -> bool:
        """
        Check if children of this node should be promoted to root level.
        
        This handles cases like "L" grouping "L1" and "LTR" when they
        should be separate top-level families.
        """
        if not node.children:
            return False
        
        # If node label is very short (1-2 chars) and children form
        # distinct "words", promote them
        if len(node.label) <= 2 and len(node.children) >= 2:
            # Check if children look like complete family names
            # (have their own subfamilies or significant count)
            distinct_families = sum(
                1 for c in node.children.values()
                if len(c.children) >= 2 or c.count >= 4
            )
            return distinct_families >= 2
        
        return False
    
    def _promote_children(self) -> None:
        """
        Post-process to promote children of stub nodes to root level.
        Applies recursively until no more promotions are needed.
        """
        changed = True
        while changed:
            changed = False
            promotions = []
            
            for label, child in list(self._compressed_root.children.items()):
                if self._should_split_to_root(child):
                    promotions.append(label)
            
            for label in promotions:
                parent = self._compressed_root.children[label]
                del self._compressed_root.children[label]
                changed = True
                
                for child_label, grandchild in parent.children.items():
                    new_label = parent.label + child_label
                    grandchild.label = new_label
                    self._compressed_root.children[new_label] = grandchild
                
                # Also handle any terminals at the parent level
                if parent.terminal_names:
                    # These become their own leaf node
                    leaf = CompressedNode(
                        label=parent.label,
                        is_terminal=True,
                        count=len(parent.terminal_names),
                        terminal_names=parent.terminal_names
                    )
                    self._compressed_root.children[parent.label] = leaf
    
    def _filter_by_support(self, node: CompressedNode) -> None:
        """
        Filter branches with insufficient support.
        
        Branches with count < min_support are merged back into parent
        by appending their label to form longer terminal names.
        """
        to_remove = []
        to_add_terminals = []
        
        for label, child in node.children.items():
            if child.count < self.min_support and not child.children:
                # Low-support leaf - merge into parent as terminal
                to_remove.append(label)
                for name in child.terminal_names:
                    to_add_terminals.append(name)
            else:
                # Recurse into children with sufficient support
                self._filter_by_support(child)
        
        for label in to_remove:
            del node.children[label]
        
        for name in to_add_terminals:
            node.terminal_names.add(name)
            node.is_terminal = True
    
    def _build_parse_lookup(self) -> None:
        """Build a lookup structure for parsing names"""
        self._prefixes: List[Tuple[str, int]] = []  # (prefix, depth)
        self._collect_prefixes(self._compressed_root, "", 0)
        # Sort by length descending for greedy matching
        self._prefixes.sort(key=lambda x: len(x[0]), reverse=True)
    
    def _collect_prefixes(self, node: CompressedNode, path: str, depth: int) -> None:
        """Collect all prefixes from the compressed trie"""
        current_path = path + node.label
        if current_path:
            self._prefixes.append((current_path, depth))
        
        for child in node.children.values():
            self._collect_prefixes(child, current_path, depth + 1)
    
    def parse_name(self, name: str) -> List[str]:
        """
        Parse a TE name into its hierarchical components based on discovered structure.
        
        Args:
            name: A TE family name (e.g., "AluYa5")
        
        Returns:
            List of components (e.g., ["Alu", "Y", "a5"])
        """
        parts = []
        remaining = name
        node = self._compressed_root
        
        while remaining and node:
            # Try to match children
            matched = False
            for child_label, child in sorted(node.children.items(), 
                                              key=lambda x: len(x[0]), reverse=True):
                if remaining.startswith(child_label):
                    parts.append(child_label)
                    remaining = remaining[len(child_label):]
                    node = child
                    matched = True
                    break
            
            if not matched:
                # No child matches - rest is the final component
                if remaining:
                    parts.append(remaining)
                break
        
        return parts if parts else [name]
    
    def get_hierarchy(self) -> dict:
        """
        Get the discovered hierarchy as a nested dictionary.
        
        Returns:
            Dict with structure like:
            {
                'Alu': {
                    'count': 100,
                    'children': {
                        'Y': {'count': 50, 'children': {...}},
                        'S': {'count': 30, 'children': {...}},
                        'J': {'count': 20, 'children': {...}},
                    }
                },
                'L1': {...}
            }
        """
        return self._compressed_root.to_dict()
    
    def get_families(self) -> List[str]:
        """Get list of top-level family names"""
        return [child.label for child in self._compressed_root.children.values()]
    
    def get_subfamilies(self, family: str) -> List[str]:
        """Get subfamilies for a given family"""
        for child in self._compressed_root.children.values():
            if child.label == family:
                return [c.label for c in child.children.values()]
        return []
    
    def print_tree(self, max_depth: int = 4) -> None:
        """Print the discovered hierarchy as a tree"""
        self._print_node(self._compressed_root, "", True, 0, max_depth)
    
    def _print_node(self, node: CompressedNode, prefix: str, is_last: bool, 
                    depth: int, max_depth: int) -> None:
        """Helper for tree printing"""
        if depth > max_depth:
            return
        
        connector = "└── " if is_last else "├── "
        if depth == 0:
            print(f"(root) [{node.count} names]")
        else:
            terminal_mark = " *" if node.is_terminal else ""
            print(f"{prefix}{connector}{node.label} [{node.count}]{terminal_mark}")
        
        if depth < max_depth:
            child_prefix = prefix + ("    " if is_last else "│   ")
            children = list(node.children.values())
            for i, child in enumerate(sorted(children, key=lambda x: -x.count)):
                self._print_node(child, child_prefix, i == len(children) - 1, 
                               depth + 1, max_depth)
    
    def to_flat_index(self) -> Dict[str, Dict[str, List[str]]]:
        """
        Convert to flat family -> subfamily -> [names] index.
        
        This produces the same structure as te_parser.build_te_index()
        but discovered automatically from the names.
        """
        index = {}
        
        def collect(node: CompressedNode, family: str, subfamily: str):
            if node.is_terminal:
                for name in node.terminal_names:
                    if family not in index:
                        index[family] = {}
                    sub_key = subfamily or '_root'
                    if sub_key not in index[family]:
                        index[family][sub_key] = []
                    index[family][sub_key].append(name)
            
            for child in node.children.values():
                if not family:
                    # Top level - this becomes the family
                    collect(child, child.label, "")
                elif not subfamily:
                    # Second level - this becomes subfamily
                    collect(child, family, child.label)
                else:
                    # Deeper levels - keep same family/subfamily
                    collect(child, family, subfamily)
        
        collect(self._compressed_root, "", "")
        return index


def build_te_index_from_names(names: List[str], min_support: int = 2) -> TETrieIndex:
    """
    Convenience function to build index from names.
    
    Args:
        names: List of TE family names
        min_support: Minimum support for a branch (default 2)
    
    Returns:
        TETrieIndex with discovered hierarchy
    """
    return TETrieIndex(names, min_support)


def load_names_from_file(filename: str, column: int = 0, delimiter: str = '\t') -> List[str]:
    """
    Load TE names from a file.
    
    Args:
        filename: Path to file (can be gzipped)
        column: Column index containing names (0-based)
        delimiter: Field delimiter
    
    Returns:
        List of unique names
    """
    import gzip
    
    names = set()
    opener = gzip.open if filename.endswith('.gz') else open
    
    with opener(filename, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split(delimiter)
            if len(parts) > column:
                names.add(parts[column])
    
    return list(names)


# ============================================================================
# Command-line interface
# ============================================================================

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("""
TE Trie Index - Discover family/subfamily structure from names

Usage:
    python te_trie.py <command> [args]

Commands:
    demo                     - Run demo with example names
    file <path> [col] [min]  - Build index from file
    parse <index_json> <name> - Parse a name using saved index

Examples:
    python te_trie.py demo
    python te_trie.py file names.txt 0 2
    python te_trie.py file hits.bed.gz 3 3
        """)
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "demo":
        # Demo with realistic TE names
        names = [
            # Alu family
            "AluY", "AluYa5", "AluYa8", "AluYb8", "AluYb9", 
            "AluYc", "AluYc3", "AluYk11", "AluYk12", "AluYk13",
            "AluSx", "AluSx1", "AluSx3", "AluSp", "AluSq", "AluSq2", "AluSz", "AluSz6",
            "AluSc", "AluSc5", "AluSc8", "AluSg", "AluSg4", "AluSg7",
            "AluJo", "AluJb", "AluJr", "AluJr4",
            # L1 family
            "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7",
            "L1P1", "L1P2", "L1P3", "L1P4",
            "L1M1", "L1M2", "L1M3", "L1M4", "L1M5",
            "L1MA1", "L1MA2", "L1MA3", "L1MA4",
            "L1MB1", "L1MB2", "L1MB3", "L1MB4",
            "L1MC1", "L1MC2", "L1MC3", "L1MC4", "L1MC4a", "L1MC5",
            "L1MD1", "L1MD2", "L1MD3",
            "L1ME1", "L1ME2", "L1ME3", "L1MEa", "L1MEb", "L1MEc", "L1MEd", "L1MEg",
            # MIR family
            "MIR", "MIRb", "MIRc", "MIR3",
            # SVA family
            "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F",
            # LTR elements
            "LTR7", "LTR7A", "LTR7B", "LTR7C", "LTR7Y",
            "LTR12", "LTR12B", "LTR12C", "LTR12D", "LTR12E",
            "THE1A", "THE1B", "THE1C", "THE1D",
            "MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1D",
            # DNA transposons
            "Charlie1", "Charlie2", "Charlie3", "Charlie15a", "Charlie16a",
            "Tigger1", "Tigger2", "Tigger3", "Tigger16a", "Tigger17a",
            "MER1A", "MER1B", "MER2", "MER5A", "MER5B", "MER20",
        ]
        
        print(f"Building index from {len(names)} names...\n")
        index = TETrieIndex(names, min_support=2)
        
        print("=== Discovered Hierarchy ===\n")
        index.print_tree(max_depth=3)
        
        print("\n=== Top-level Families ===")
        for fam in index.get_families():
            subs = index.get_subfamilies(fam)
            print(f"  {fam}: {len(subs)} subfamilies - {subs[:5]}{'...' if len(subs) > 5 else ''}")
        
        print("\n=== Parsing Examples ===")
        test_names = ["AluYa5", "AluYk11", "AluSx1", "L1PA2", "L1MC4a", "L1MEg", "SVA_A", "THE1B"]
        for name in test_names:
            parts = index.parse_name(name)
            print(f"  {name:15s} -> {parts}")
        
        print("\n=== Flat Index (family -> subfamily -> names) ===")
        flat = index.to_flat_index()
        for family in sorted(flat.keys())[:3]:
            print(f"\n  {family}:")
            for subfam, members in sorted(flat[family].items())[:4]:
                print(f"    {subfam}: {members[:3]}{'...' if len(members) > 3 else ''}")
    
    elif command == "file":
        if len(sys.argv) < 3:
            print("Usage: te_trie.py file <path> [column] [min_support]")
            sys.exit(1)
        
        filepath = sys.argv[2]
        column = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        min_support = int(sys.argv[4]) if len(sys.argv) > 4 else 2
        
        print(f"Loading names from {filepath} (column {column})...")
        names = load_names_from_file(filepath, column)
        print(f"Found {len(names)} unique names")
        
        print(f"\nBuilding index (min_support={min_support})...")
        index = TETrieIndex(names, min_support)
        
        print("\n=== Discovered Hierarchy ===\n")
        index.print_tree(max_depth=3)
        
        print("\n=== Top-level Families ===")
        for fam in index.get_families():
            subs = index.get_subfamilies(fam)
            print(f"  {fam}: {len(subs)} subfamilies")
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
