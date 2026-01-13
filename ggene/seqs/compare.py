
from pathlib import Path
import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

import subprocess
import numpy as np

from Bio.Phylo import NewickIO

from tabulate import tabulate

from scipy import cluster

# from ggene.database.genome_manager import GenomeManager
from ggene import draw
from ggene.draw import Heatmap, ScalarPlot
from ggene.seqs import bio, process, align
from ggene.seqs.bio import ALIASES, ALIASES_REV, get_aliases, VOCAB

@dataclass
class ClusterHit:
    query_label:str
    index:int
    length:int
    pct_id:float
    strand:str
    compressed_alignment:str
    cluster_label:str
    target_label:str
    seq:str = ""
    is_centroid:bool = False

    @classmethod
    def from_line(cls, line, seq = "", curr_cluster = ""):
        
        parts = line.split("\t")
        rt = parts[0].strip()
        
        if rt not in ["S","H"]:
            # print("wrong entry type")
            return None
        
        query_label = parts[8].strip()
        index = int(parts[1].strip())
        length = int(parts[2].strip())
        pct_id = float(parts[3].strip()) if not "*" in parts[3] else 100.0
        strand = parts[4].strip()
        cmpa = parts[7].strip() if not "*" in parts[5] else ""
        target_label = int(parts[9].strip()) if not "*" in parts[7] else ""
        is_centroid = rt=="S"
        if is_centroid:
            curr_cluster = query_label
        
        return cls(query_label, index, length, pct_id, strand, cmpa, curr_cluster, target_label, seq = seq, is_centroid = is_centroid)

@dataclass
class Cluster:
    query_label:str
    index:int
    num_hits:int
    
    @property
    def is_centroid(self):
        return False
    
    @classmethod
    def from_line(cls, line):
        
        parts = line.split("\t")
        rt = parts[0].strip()
        if rt != "C":
            # print("wrong entry type")
            return None
        
        query_label = parts[8].strip()
        index = int(parts[1].strip())
        num_hits = int(parts[2].strip())
        
        return cls(query_label, index, num_hits)
        

class ClusterResult:
    
    def __init__(self):
        self.hits = []
        self.clusters = []
        self.clustered_hits = {}
        
    def add_hit(self, hit):
        self.hits.append(hit)
        if hit.is_centroid:
            self.clustered_hits[hit.query_label] = []
        if not hit.cluster_label in self.clustered_hits:
            print(f"can't add hit {hit} to clustered hits..")
        else:
            self.clustered_hits[hit.cluster_label].append(hit)
    
    def add_cluster(self, cluster):
        self.clusters.append(cluster)

    def print_clusters(self, show_sequences: bool = False, max_seq_len: int = 60):
        """Print clustering results in a tree-like format.

        Note: This is NOT a phylogenetic tree - clusters are independent groups
        based on similarity threshold, not evolutionary relationships.
        """
        print(f"\n{'='*80}")
        print(f"Clustering Results: {len(self.clustered_hits)} clusters, {len(self.hits)} total sequences")
        print(f"{'='*80}\n")

        # Sort clusters by size (largest first)
        sorted_clusters = sorted(
            self.clustered_hits.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )

        for cluster_idx, (centroid_label, hits) in enumerate(sorted_clusters, 1):
            # Find the centroid hit
            centroid = next((h for h in hits if h.is_centroid), None)
            if not centroid:
                continue

            cluster_size = len(hits)
            print(f"Cluster {cluster_idx}: {cluster_size} sequences (id ≥ threshold)")
            print(f"├─ Centroid: {centroid.query_label}")

            if show_sequences and centroid.seq:
                print(f"│  └─ Seq: {centroid.seq[:max_seq_len]}{'...' if len(centroid.seq) > max_seq_len else ''}")

            # Print hits (non-centroid members)
            non_centroid_hits = [h for h in hits if not h.is_centroid]

            for hit_idx, hit in enumerate(non_centroid_hits):
                is_last = (hit_idx == len(non_centroid_hits) - 1)
                connector = "└─" if is_last else "├─"

                print(f"│  {connector} Hit: {hit.query_label} ({hit.pct_id:.1f}% id to centroid)")

                if show_sequences and hit.seq:
                    extension = "   " if is_last else "│  "
                    print(f"│  {extension}└─ Seq: {hit.seq[:max_seq_len]}{'...' if len(hit.seq) > max_seq_len else ''}")

            print()

    def print_summary(self):
        """Print summary statistics about the clustering."""
        cluster_sizes = [len(hits) for hits in self.clustered_hits.values()]

        print(f"\nClustering Summary:")
        print(f"  Total sequences: {len(self.hits)}")
        print(f"  Total clusters: {len(self.clustered_hits)}")
        print(f"  Largest cluster: {max(cluster_sizes) if cluster_sizes else 0} sequences")
        print(f"  Smallest cluster: {min(cluster_sizes) if cluster_sizes else 0} sequences")
        print(f"  Average cluster size: {np.mean(cluster_sizes) if cluster_sizes else 0:.1f}")
        print(f"  Singletons (clusters of 1): {sum(1 for s in cluster_sizes if s == 1)}")
        print()

    def get_cluster_by_label(self, label: str):
        """Get all sequences in a cluster by centroid label."""
        return self.clustered_hits.get(label, [])

    def add_from_line(self, line, seq = "", curr_cluster = ""):
        
        hit = ClusterHit.from_line(line, seq=seq, curr_cluster=curr_cluster)
        if hit:
            self.add_hit(hit)
            return hit
        
        clstr = Cluster.from_line(line)
        if clstr:
            self.add_cluster(clstr)
            return clstr
        
    @classmethod
    def from_uc(cls, uc_path):
        
        comms = []
        
        cres = cls()
        curr_cluster = ""
        with open(uc_path) as f:
            for line in f:
                
                if line.startswith("#"):
                    comms.append(line)
                    continue
                
                rec = cres.add_from_line(line, curr_cluster = curr_cluster)
                
                if rec.is_centroid:
                    curr_cluster = rec.cluster_label
                
        
        return cres, comms

class PhyloLeaf:

    def __init__(self, name, branch_length, data = None, parent = None):
        self.name = name
        self.branch_length = branch_length
        self.data = data
        self.parent = parent

    @property
    def is_leaf(self):
        return True

    @property
    def is_branch(self):
        return False

    @property
    def id(self):
        return self.name

    @property
    def num_leaves(self):
        return 1

    @property
    def depth(self):
        return self.parent.depth + 1 if self.parent else 0

    def set_parent(self, parent):
        self.parent = parent

    def bind_data(self, data:Dict[str,Any]):
        if self.name in data:
            self.data = data[self.name]

    def branch(self, other, branch_length = 0.0, confidence = 0.0):
        return PhyloBranch(branch_length, confidence, self, other)

    def draw(self, prefix="", is_last=True, formatter=None, show_length=True):
        """Draw this leaf with tree characters.

        Args:
            prefix: String prefix for indentation
            is_last: Whether this is the last child
            formatter: Optional function to format self.data, e.g. lambda d: d[:50]
            show_length: Whether to show branch lengths
        """
        connector = "└─ " if is_last else "├─ "

        # Format the data
        if formatter and self.data:
            data_str = formatter(self.data)
        elif self.data:
            data_str = str(self.data)
            if len(data_str) > 60:
                data_str = data_str[:57] + "..."
        else:
            data_str = ""

        # Build the line
        branch_str = f" [{self.branch_length:.4f}]" if show_length else ""
        line = f"{prefix}{connector}{self.name}{branch_str}"
        if data_str:
            line += f": {data_str}"

        print(line)

    def __iter__(self):
        """Iterate over this node (yields self)."""
        yield self

    def __len__(self):
        return 1

    @classmethod
    def from_clade(cls, clade, data:Dict[str,Any] = {}):

        if not clade.is_terminal():
            return PhyloBranch.from_clade(clade, data=data)
        else:
            leafdata = data.get(clade.name)
            leaf = cls(clade.name, clade.branch_length, data = leafdata)
            return leaf

    def print(self, tabs = 0, only_data = False):
        tabstr = "  " * tabs
        if only_data:
            print(tabstr, repr(self.data))
        else:
            print(tabstr, repr(self))

    def __repr__(self):

        parts = []
        parts.append(f"branch_length={self.branch_length:0.1%}")
        if self.data:
            datastr = repr(self.data)
            if len(datastr) > 64:
                datastr = datastr[:62] + ".."
            parts.append(f"data={datastr}")
        parts.append(f"parent={self.parent.id[-6:]}")
        parts.append(f"d={self.depth}")

        return f"PhyloLeaf({self.name}, {", ".join(parts)})"
    
class PhyloBranch:
    
    def __init__(self, branch_length, confidence, *limbs, parent = None):
        self.branch_length = branch_length
        self.confidence = confidence
        self.limbs = []
        self.parent = parent
        for lmb in limbs:
            self.add_limb(lmb)
    
    def add_limb(self, limb):
        self.limbs.append(limb)
        limb.set_parent(self)
        
    def set_parent(self, parent):
        self.parent = parent
    
    def bind_data(self, data:Dict[str,Any]):
        for limb in self.limbs:
            limb.bind_data(data)
        
    @property
    def id(self):
        return hex(id(self))[-4:]

    @property
    def is_leaf(self):
        return False
    
    @property
    def is_branch(self):
        return True
    
    @property
    def depth(self):
        return self.parent.depth + 1 if self.parent else 0
    
    @property
    def num_leaves(self):
        return sum(lmb.num_leaves for lmb in self.limbs)
    
    def __getitem__(self, ind):
        """Get node by index or name (recursive search)."""
        if isinstance(ind, int):
            return self.limbs[ind]
        elif isinstance(ind, str):
            # Check if this branch's ID matches
            if self.id == ind:
                return self

            # Recursively search all limbs
            for limb in self.limbs:
                if limb.is_leaf:
                    if limb.name == ind:
                        return limb
                else:
                    # Recursively search this branch
                    try:
                        result = limb[ind]
                        if result is not None:
                            return result
                    except (KeyError, IndexError):
                        continue

            # Not found in any limb
            raise KeyError(f"Node '{ind}' not found in tree")
        else:
            raise TypeError(f"Index must be int or str, not {type(ind)}")

    def __iter__(self):
        """Iterate over all nodes in this branch (depth-first)."""
        # Yield self first
        yield self
        # Then recursively yield from all limbs
        for limb in self.limbs:
            yield from limb

    def __len__(self):
        return len(self.limbs) + sum(len(limb) for limb in self.limbs)

    @classmethod
    def from_clade(cls, clade, data:Dict[str, Any]={}):
        
        if clade.is_terminal():
            return PhyloLeaf.from_clade(clade, data=data)
        else:
            limbs = []
            for subclade in clade.clades:
                sb = cls.from_clade(subclade, data=data)
                limbs.append(sb)
            
            branch = cls(clade.branch_length, clade.confidence, *limbs)
            return branch
    
    def iterleaves(self):
        """Iterate over all leaf nodes in this branch."""
        for lmb in self.limbs:
            if lmb.is_leaf:
                yield lmb
            else:
                yield from lmb.iterleaves()

    def draw(self, prefix="", is_last=True, formatter=None, show_length=True, show_confidence=True):
        """Draw this branch with tree characters.

        Args:
            prefix: String prefix for indentation
            is_last: Whether this is the last child
            formatter: Optional function to format leaf data
            show_length: Whether to show branch lengths
            show_confidence: Whether to show confidence values
        """
        connector = "└─ " if is_last else "├─ "
        extension = "   " if is_last else "│  "

        # Build branch info
        info_parts = []
        if show_length and self.branch_length:
            info_parts.append(f"len={self.branch_length:.4f}")
        if show_confidence and self.confidence:
            info_parts.append(f"conf={self.confidence:.2f}")
        info_parts.append(f"{self.num_leaves} leaves")

        info_str = ", ".join(info_parts) if info_parts else ""
        line = f"{prefix}{connector}Branch {self.id[-4:]}"
        if info_str:
            line += f" ({info_str})"

        print(line)

        # Draw children
        for i, limb in enumerate(self.limbs):
            is_last_limb = (i == len(self.limbs) - 1)
            limb.draw(
                prefix=prefix + extension,
                is_last=is_last_limb,
                formatter=formatter,
                show_length=show_length
            )

    def print(self, tabs = 0, only_data = False):
        tabstr = "  " * tabs
        if not only_data:
            print(tabstr, repr(self))
        for lmb in self.limbs:
            lmb.print(tabs=tabs+1)

    def __repr__(self):
        
        parts = []
        if self.branch_length:
            parts.append(f"branch_length={self.branch_length:0.1%}")
        if self.confidence:
            parts.append(f"confidence={self.confidence:0.2f}")
        parts.append(f"num_limbs={len(self.limbs)}")
        parts.append(f"parent={self.parent.id}")
        
        return f"PhyloBranch({self.id}, {", ".join(parts)})"

class PhyloTree:

    def __init__(self, name, root = None):
        self.name = name
        self.root:Optional[PhyloBranch] = None
        self._node_index = None  # Lazy-built index of name -> node
        if root:
            self.set_root(root=root)

    def set_root(self, root):
        self.root = root
        self.root.set_parent(self)
        self._node_index = None  # Invalidate index

    def bind_data(self, data:Dict[str,Any]):
        self.root.bind_data(data)

    def _build_index(self):
        """Build name->node index for fast lookup."""
        self._node_index = {}
        for node in self:
            if node.is_leaf:
                self._node_index[node.name] = node
            else:
                self._node_index[node.id] = node

    @property
    def id(self):
        return hex(id(self))[-4:]

    @property
    def num_leaves(self):
        return self.root.num_leaves

    @property
    def depth(self):
        return 0

    def __getitem__(self, ind):
        """Get node by index or name (uses cached index for speed)."""
        if isinstance(ind, str):
            # Build index on first use
            if self._node_index is None:
                self._build_index()

            if ind in self._node_index:
                return self._node_index[ind]
            else:
                raise KeyError(f"Node '{ind}' not found in tree")
        else:
            return self.root[ind]

    def __iter__(self):
        """Iterate over all nodes in the tree."""
        return iter(self.root)

    def iterleaves(self):
        """Iterate over all leaf nodes."""
        return self.root.iterleaves()

    def draw(self, formatter=None, tree_width=None, consensus = ""):
        """Draw the phylogenetic tree with ASCII art.

        Args:
            formatter: Optional function to format leaf data
            tree_width: Width for tree structure (depth in characters)
        """
        print(f"PhyloTree: {self.name} ({self.num_leaves} leaves)\n")

        # Assign vertical row position to each leaf
        positions = {}
        self._assign_vertical_positions(self.root, positions)

        # Assign horizontal depth to each node and calculate max depth
        depths = {}
        max_depth = self._assign_horizontal_depths(self.root, 0, depths, max_width=999)

        # Auto-calculate tree_width if not provided (give 2-3 chars per depth level)
        if tree_width is None:
            tree_width = max_depth * 2 + 2

        # Build 2D grid for drawing
        leaves = list(self.iterleaves())
        num_rows = len(leaves)
        grid = [[' ' for _ in range(tree_width + 1)] for _ in range(num_rows)]

        # Draw tree structure into grid
        self._draw_tree_structure(self.root, grid, positions, depths, tree_width, max_depth)

        if consensus:
            print(" "*(tree_width + 23) + consensus)

        # Print each row with leaf info
        for i, leaf in enumerate(leaves):
            row = positions[id(leaf)]
            tree_str = ''.join(grid[row]).rstrip()

            # Format data
            if formatter and leaf.data:
                data_str = formatter(leaf.data)
            elif leaf.data:
                data_str = str(leaf.data)
            else:
                data_str = ""

            print(f" {tree_str} {leaf.name:>12} ({leaf.branch_length:0.1%}) {data_str}")

    def _assign_vertical_positions(self, node, positions, counter=None):
        """Assign vertical row positions to nodes (leaves only get positions)."""
        if counter is None:
            counter = [0]

        if node.is_leaf:
            positions[id(node)] = counter[0]
            counter[0] += 1
            return counter[0] - 1
        else:
            child_positions = []
            for limb in node.limbs:
                child_positions.append(self._assign_vertical_positions(limb, positions, counter))
            # Branch position is midpoint of children
            positions[id(node)] = int(sum(child_positions) / len(child_positions))
            return positions[id(node)]

    def _assign_horizontal_depths(self, node, depth, depths, max_width):
        """Assign horizontal depth (column) to each node."""
        depths[id(node)] = min(depth, max_width)
        max_d = depth

        if not node.is_leaf:
            for limb in node.limbs:
                max_d = max(max_d, self._assign_horizontal_depths(limb, depth + 1, depths, max_width))

        return max_d

    def _draw_tree_structure(self, node, grid, positions, depths, tree_width, max_depth):
        """Recursively draw tree branches into grid with proper spacing."""
        # Calculate column positions (map depth to grid column with spacing)
        chars_per_level = tree_width / max(max_depth, 1)

        def depth_to_col(d):
            """Convert tree depth to grid column."""
            return int(d * chars_per_level)

        if node.is_leaf:
            # Draw horizontal line to leaf (fill remaining width)
            row = positions[id(node)]
            start_col = depth_to_col(depths[id(node)])
            for x in range(start_col, len(grid[0])):
                if grid[row][x] == ' ':
                    grid[row][x] = '─'
        else:
            # Draw branch point and children
            node_row = positions[id(node)]
            node_col = depth_to_col(depths[id(node)])

            # Get children positions
            child_rows = [positions[id(limb)] for limb in node.limbs]
            min_row = min(child_rows)
            max_row = max(child_rows)

            # Draw vertical connector between first and last child
            for r in range(min_row, max_row + 1):
                if grid[r][node_col] == ' ':
                    grid[r][node_col] = '│'

            # Draw each child branch
            for i, limb in enumerate(node.limbs):
                child_row = positions[id(limb)]
                child_col = depth_to_col(depths[id(limb)])
                is_last_child = (i == len(node.limbs) - 1)

                # Mark branch connection at node column
                if child_row == min_row and i == 0:
                    grid[child_row][node_col] = '┌'
                elif is_last_child:
                    grid[child_row][node_col] = '└'
                else:
                    grid[child_row][node_col] = '├'

                # Draw horizontal line from node to child
                for x in range(node_col + 1, child_col):
                    if grid[child_row][x] == ' ':
                        grid[child_row][x] = '─'

                # Mark child branch point (if it's a branch node)
                if not limb.is_leaf and child_col < len(grid[0]):
                    # Place '┬' at the child's column to show it branches
                    grid[child_row][child_col] = '┬'

                # Recursively draw children
                self._draw_tree_structure(limb, grid, positions, depths, tree_width, max_depth)

    @classmethod
    def from_newick(cls, name, newick_tree:str, data:Dict[str, Any]={}):

        prs = NewickIO.Parser
        ntree = next(iter(prs.from_string(newick_tree).parse()))

        root = PhyloBranch.from_clade(ntree.root, data=data)
        tree = cls(name, root=root)

        return tree

    def print(self, only_data = False):
        print(repr(self))
        self.root.print(tabs=1, only_data = only_data)

    def __repr__(self):
        return f"PhyloTree({self.name}, num_limbs={len(self)})"

    def __len__(self):
        return len(self.root)


def count_matching_wobble(s1: str, s2: str, err_tol: Optional[int] = None, wobble_tol = None, allow_alias = False) -> int:
    """Helper to compare two sequences and count mismatches."""
    if len(s1) < len(s2):
        return len(s1)
    
    if err_tol is None:
        err_tol = len(s1)
    if wobble_tol is None:
        wobble_tol = len(s1)

    ia, iar = get_aliases(VOCAB)
    
    nerr = 0
    nwob = 0
    for a, b in zip(s1, s2):
        if a == b:
            pass
        elif allow_alias and a in ia and b in ia[a]:
            pass
        elif allow_alias and b in ia and a in ia[b]:
            pass
        elif a+b in ["TG","GT"]:
                nwob += 1
        elif a+b in ["AG","GA"]:
                nwob += 1
        elif a+b in ["AC","CA"]:
                nwob += 1
        else:
            nerr += 1
            
        if nerr > err_tol:
            return nerr
        if nwob > wobble_tol:
            return nwob
        
    return nerr

def count_matching(s1: str, s2: str, err_tol: Optional[int] = None) -> int:
    """Helper to compare two sequences and count mismatches."""
    
    if err_tol is None:
        err_tol = len(s1)
    
    match = 0
    nerr = 0
    for a, b in zip(s1, s2):
        if a == b:
            match += 1
            pass
        else:
            nerr += 1
            
        if nerr > err_tol:
            break
        
    return match, nerr


def score_sequences_corrs(seqa, seqb, topk = 5, scale = None, **kwargs):
    min_len = min(len(seqa), len(seqb))
    corrs, rccorrs = process.correlate(seqa, seqb, scale = scale, fill = 0.25)
    # corrs, rccorrs = process.correlate_fast(seqa, seqb)
    topk_corrs = sorted(corrs)[:topk]
    topk_rccorrs = sorted(rccorrs)[:topk]
    return sum(topk_corrs)/topk, sum(topk_rccorrs)/topk

def score_sequences_runs(seqa, seqb, topk = 5, max_err = 16, scale = None, **kwargs):
    runs, _, shifts, _ = process.correlate_longest_subseq_err(seqa, seqb, max_err, scale=scale)
    topk_runs = sorted([(r,s) for r,s in zip(runs, shifts) if s > 3], key = lambda k:-k[0])[:topk]    
    rcruns, _, rcshifts, _ = process.correlate_longest_subseq_err(seqa, bio.reverse_complement(seqb), max_err, scale=scale)
    topk_rcruns = sorted([(r,s) for r,s in zip(rcruns, rcshifts) if s > 3], key = lambda k:-k[0])[:topk]    
    return sum([r for r, s in topk_runs])/topk, sum([r for r, s in topk_rcruns])/topk

def score_sequences_algn(seqa, seqb, **kwargs):
    min_len = min(len(seqa), len(seqb))
    score = align.score_sequences(seqa, seqb)
    rcscore = align.score_sequences(seqa, bio.reverse_complement(seqb))
    return score/min_len, rcscore/min_len

def downsample_scores(scores):
    
    scores = np.array(scores)
    nr, nc = scores.shape
    
    ds_score = 0.5*scores[1:nr-1, 1:nc-1]
    ds_score += 0.25*scores[2:nr, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 2:nc]
    ds_score += 0.25*scores[0:nr-2, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 0:nc-2]
    
    return ds_score

def get_score_function(score_mode):
    if score_mode == "alignment":
        return score_sequences_algn
    elif score_mode == "runs":
        return score_sequences_runs
    elif score_mode == "corrs":
        return score_sequences_corrs
    else:
        return None

def calc_sequence_comparison(seqa, seqb, nchksa, nchksb, chunksz, score_func, q = 1):
    
    scores = [[] for n in range(nchksb)]
    rcscores = [[] for n in range(nchksb)]
    
    for ia in range(nchksa):
        ssa = seqa[chunksz*ia:chunksz*(ia+q)]
        
        for ib in range(nchksb):
            
            ssb = seqb[chunksz*ib:chunksz*(ib+q)]
            if not ssa or not ssb:
                continue
                
            score, rcscore = score_func(ssa, ssb)
            scores[ib].append(score)
            rcscores[ib].append(rcscore)
    
    return scores, rcscores

def compare_sequences(seqa, seqb, chunksz = 128, score_modes = ["alignment", "runs", "corrs"], resample = True, **kwargs):
    
    # if len(seqb) > len(seqa):
    #     seqa, seqb = seqb, seqa
    
    q = 1
    if resample:
        q = 2
        chunksz = chunksz//2
    
    nchksa = len(seqa) // chunksz
    nchksb = len(seqb) // chunksz
    
    score_funcs = [get_score_function(sm) for sm in score_modes]
    
    scores = [[[] for n in range(nchksb)] for sf in score_funcs]
    rcscores = [[[] for n in range(nchksb)] for sf in score_funcs]
    row_lbls = [f"b{n}" for n in range(nchksb)]
    
    file = kwargs.pop("file", None)
    ruler = kwargs.pop("add_ruler", False)
    xmin = kwargs.pop("xmin", None)
    xmax = kwargs.pop("xmax", None)
    num_labels = kwargs.pop("num_labels", 5)
    half_block = kwargs.pop("half_block", True)
    
    cs = kwargs.pop("color_scheme", "moss")
    
    for ia in range(nchksa):
        ssa = seqa[chunksz*ia:chunksz*(ia+q)]
        
        for ib in range(nchksb):
            
            ssb = seqb[chunksz*ib:chunksz*(ib+q)]
            if not ssa or not ssb:
                continue
                
            for i, sf in enumerate(score_funcs):
                score, rcscore = sf(ssa, ssb, **kwargs)
                scores[i][ib].append(score)
                rcscores[i][ib].append(rcscore)
    
    all_hms = [[] for r in range(len(scores[0])+3)]
    all_rchms = [[] for r in range(len(rcscores[0])+3)]
    
    for nsm,sm in enumerate(score_modes):
        
        if nsm > 0:
            row_lbls = []
        
        if resample:
            scores_rs = downsample_scores(scores[nsm])
            rcscores_rs = downsample_scores(rcscores[nsm])
        else:
            scores_rs = scores[nsm]
            rcscores_rs = rcscores[nsm]
        
        hm = Heatmap(scores_rs, center = None, minval = 0, row_labels = row_lbls, suppress = True, color_scheme = cs, col_space = 0, row_space = 0, add_middle = False, half_block = half_block , colorbar = True)
        
        rows = hm.get_rows()
        for nr in range(len(rows)):
            all_hms[nr].append(rows[nr])
        
        rchm = Heatmap(rcscores_rs, center = None, minval = 0, row_labels = row_lbls, suppress = True, color_scheme = cs, col_space = 0, row_space = 0, add_middle = False, half_block = half_block, colorbar = True)
        
        rcrows = rchm.get_rows()
        for nr in range(len(rcrows)):
            all_hms[nr].append(rcrows[nr])
    
    frm = "{:^32}       " * 2
    for row in all_hms:
        if row:
            print(frm.format(*row), file = file)
    # for row in all_rchms:
    #     if row:
    #         print(frm.format(*row))
    
    return scores


def walk_sequences(seqa, seqb, chunksz, num_steps, **kwargs):
    
    asublen = len(seqa)//num_steps
    bsublen = len(seqb)//num_steps
    
    print(f"a sublength: {asublen}, b sublength: {bsublen}")
    
    for n in range(num_steps):
        
        subseqa = seqa[n*asublen:(n+1)*asublen]
        subseqb = seqb[n*bsublen:(n+1)*bsublen]
        
        compare_sequences(subseqa, subseqb, chunksz=chunksz, **kwargs)


def organize_sequences(seqs, score_func=None, ref = ""):
    
    if score_func is None:
        score_func = lambda s1, s2: count_matching(s1, s2)
    
    if not ref:
        ref = seqs.pop(0)
    
    scores = [score_func(ref, seq) for seq in seqs]
    
    new_ord = sorted(enumerate(scores), key = lambda k: -k[1])
    seqs_org = [ref] + [seqs[i] for i, sc in new_ord]
    
    return seqs_org, [i for i, sc in new_ord]

def cluster_sequences(seqs:Dict[str,str], name = "clstr", id_thresh = 0.9):
    """Cluster sequences using usearch cluster_fast (greedy clustering).

    This performs FLAT clustering based on similarity threshold, NOT phylogenetic analysis.

    Algorithm:
      1. First sequence becomes first centroid (S)
      2. Each subsequent sequence compared to existing centroids
      3. If similarity ≥ threshold → assigned as Hit (H) to closest centroid
      4. If similarity < threshold → becomes new centroid (S)
      5. Result: independent clusters, ordered by size or creation

    Record types in output:
      S: Centroid/Seed - representative sequence of cluster
      H: Hit - sequence assigned to a centroid
      C: Cluster summary record
      N: No hit (if sequence doesn't match any cluster)

    Args:
        seqs: Dict of {name: sequence}
        name: Name for temp files
        id_thresh: Identity threshold (0.0-1.0)

    Returns:
        ClusterResult with clusters organized by centroid

    For PHYLOGENETIC TREES (evolutionary relationships), use build_phylogenetic_tree() instead.
    """

    temp_in = Path(f"/tmp/{name}_input.fasta")
    temp_out = Path(f"/tmp/{name}_clusters.uc")

    with open(temp_in, 'w') as f:
        for seq_name, seq in seqs.items():
            f.write(f">{seq_name}\n{seq}\n")

    subprocess.run([
        "usearch", "-cluster_fast",
        str(temp_in),
        "-id", format(id_thresh, "0.1f"),
        "-uc", str(temp_out)
    ], check=True, capture_output = True)

    clres, comms = ClusterResult.from_uc(temp_out)

    temp_in.unlink()
    temp_out.unlink()

    return clres


def build_distance_matrix(seqs: Dict[str, str], metric: str = 'hamming'):
    """Build pairwise distance matrix for phylogenetic analysis.

    Args:
        seqs: Dict of {name: sequence}
        metric: 'hamming' (substitutions only) or 'edit' (with indels)

    Returns:
        names (list), distance_matrix (2D numpy array)
    """
    names = list(seqs.keys())
    n = len(names)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            seq1, seq2 = seqs[names[i]], seqs[names[j]]

            if metric == 'hamming':
                # Simple hamming distance (requires equal length)
                min_len = min(len(seq1), len(seq2))
                dist = sum(a != b for a, b in zip(seq1[:min_len], seq2[:min_len]))
                dist += abs(len(seq1) - len(seq2))  # Add length difference
            elif metric == 'edit':
                # Edit/Levenshtein distance
                from Levenshtein import distance
                dist = distance(seq1, seq2)
            elif metric == 'align':
                pass
                
            else:
                raise ValueError(f"Unknown metric: {metric}")

            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return names, dist_matrix


def build_phylogenetic_tree(seqs: Dict[str, str], method: str = 'fasttree', **kwargs):
    """Build a proper phylogenetic tree (evolutionary relationships).

    This is different from cluster_sequences() which does flat clustering.

    Args:
        seqs: Dict of {name: sequence}
        method: 'fasttree' (fast ML), 'muscle+nj' (neighbor-joining), or 'distance' (simple NJ from distances)

    Returns:
        Tree structure (depends on method)

    Recommended workflow for phylogenetics:
      1. Align sequences with MUSCLE: align.align_sequences_muscle()
      2. Build tree with FastTree or RAxML
      3. Visualize with ete3 or BioPython
    """
    if method == 'fasttree':
        return build_phylogenetic_tree_fasttree(seqs, **kwargs)
        
    elif method == 'distance':
        # Simple distance-based tree from hamming/edit distances
        names, dist_matrix = build_distance_matrix(seqs, metric='hamming')
        print(f"Distance matrix built for {len(names)} sequences")
        print("Use scipy.cluster.hierarchy or scikit-bio for tree building from distances")
        # cluster.hierarchy.fcluster()
        
        return names, dist_matrix

    else:
        raise ValueError(f"Unknown method: {method}")

def build_phylogenetic_tree_fasttree(seqs:Dict[str,str], name = "phylo_tree", use_stdin: bool = True, **kwargs):
    """Build phylogenetic tree using FastTree (maximum likelihood).

    Args:
        seqs: Dict of sequences (will be aligned with MUSCLE first)
        name: Name for temp files
        use_stdin: If True, pipe FASTA to stdin (no temp file). If False, use temp file.

    Returns:
        Newick format tree string
    """

    # Align sequences first
    seqs_aligned = align.align_sequences_muscle(seqs, name=name)

    if use_stdin:
        # Option 1: Pipe FASTA directly to stdin (no temp files!)
        fasta_str = ""
        for seq_name, seq in seqs_aligned.items():
            fasta_str += f">{seq_name}\n{seq}\n"

        result = subprocess.run(
            ["fasttree", "-nt"],  # Read from stdin
            input=fasta_str,
            capture_output=True,
            text=True,
            check=True
        )

    else:
        # Option 2: Use temp file (more compatible with large datasets)
        temp_in = Path(f"/tmp/{name}_algn.fasta")

        with open(temp_in, 'w') as f:
            for seq_name, seq in seqs_aligned.items():
                f.write(f">{seq_name}\n{seq}\n")

        result = subprocess.run(
            ["fasttree", "-nt", str(temp_in)],
            capture_output=True,
            text=True,
            check=True
        )

        temp_in.unlink()

    # FastTree outputs Newick tree to stdout
    tree_newick = result.stdout.strip()
    
    pt = PhyloTree.from_newick("testtree", tree_newick, data = seqs_aligned)
        
    return pt
    