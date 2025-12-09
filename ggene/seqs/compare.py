

import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

import numpy as np

from tabulate import tabulate


from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.seqs import bio, process, find, align, heal
from ggene.seqs.bio import ALIASES, ALIASES_REV, get_aliases, VOCAB


def count_matching_wobble(s1: str, s2: str, err_tol: Optional[int] = None, wobble_tol = None, allow_alias = False) -> int:
    """Helper to compare two sequences and count mismatches."""
    if len(s1) < len(s2):
        return len(s1)
    
    if err_tol is None:
        err_tol = len(s1)
    if wobble_tol is None:
        # wobble_tol = len(s1)
        allow_wobble = False
    else:
        allow_wobble = True

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
        elif allow_wobble and a+b in ["TG","GT"]:
                nwob += 1
        elif allow_wobble and a+b in ["AG","GA"]:
                nwob += 1
        elif allow_wobble and a+b in ["AC","CA"]:
                nwob += 1
        else:
            nerr += 1
            
        if nerr > err_tol:
            return nerr
        if allow_wobble and nwob > wobble_tol:
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

def compare_sequences(seqa, seqb, chunksz = 128, score_modes = ["alignment", "runs", "corrs"], resample = True, **kwargs):
    
    if len(seqb) > len(seqa):
        seqa, seqb = seqb, seqa
    
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
    
    ruler = kwargs.pop("add_ruler", False)
    xmin = kwargs.pop("xmin", None)
    xmax = kwargs.pop("xmax", None)
    num_labels = kwargs.pop("num_labels", 5)
    
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
        
        hm = draw.heatmap(scores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra", col_space = 0, row_space = 0,
                          ruler = ruler, xmin = xmin, xmax = xmax, num_labels=num_labels)
        
        for nr in range(len(hm)):
            all_hms[nr].append(hm[nr])
        
        rchm = draw.heatmap(rcscores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra", col_space = 0, row_space = 0,
                          ruler = ruler, xmin = xmin, xmax = xmax, num_labels=num_labels)
        
        for nr in range(len(rchm)):
            all_rchms[nr].append(rchm[nr])
    
    frm = "{:^32}       " * len(score_modes)
    for row in all_hms:
        if row:
            print(frm.format(*row))
    for row in all_rchms:
        if row:
            print(frm.format(*row))
    
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



class PhyloLeaf:

    _print_frm = "{info:<32}{seq} ({extra})"
    
    def __init__(self, seq, score=0.0):
        self.seq = seq
        self.parent = None
        self.score = score

    @property
    def depth(self):
        return self.parent.depth + 1 if self.parent else 0

    def get_all_leaves(self):
        return [self]

    def compare(self, new_seq):
        score, consensus = Phylo.score(self.seq, new_seq)
        return score, consensus

    def insert(self, new_seq):
        score, consensus = self.compare(new_seq)
        new_leaf = PhyloLeaf(new_seq, score=score)
        new_node = PhyloNode(self, new_leaf, parent=self.parent, consensus=consensus)
        return new_node

    def recompute_consensus(self):
        return self.seq

    def print_tree(self, prefix="", is_last=True, show_length = 32):   
        connector = "└─ " if is_last else "├─ "
        extension = "   " if is_last else "│  "
        
        info = f"{prefix}{connector}Leaf: "
        extra = f"score={self.score:.2f}"
        
        row = self._print_frm.format(info=info, seq=self.seq[:show_length], extra=extra)
        print(row)

    def __hash__(self):
        return hash(self.seq)

    def __str__(self):
        return self.seq

    def __repr__(self):
        return f"PhyloLeaf({self.seq[:32]}...,score={self.score:.2f},depth={self.depth})"

class PhyloNode:

    _print_frm = "{info:<32}{seq} ({extra})"

    def __init__(self, childa, childb, parent=None, consensus="", score = 0.0):

        self.childa: Union[PhyloLeaf, 'PhyloNode'] = childa
        self.childa.parent = self

        self.childb: Union[PhyloLeaf, 'PhyloNode'] = childb
        self.childb.parent = self

        self.parent = parent
        self.consensus = consensus
        self.score = score

    @property
    def depth(self):
        return self.parent.depth + 1 if self.parent else 0

    def get_all_leaves(self):
        return self.childa.get_all_leaves() + self.childb.get_all_leaves()

    def compare(self, new_seq):
        sca, consa = Phylo.score(str(self.consensus), new_seq)
        return sca, consa

    def insert(self, new_seq):
        # Compare new sequence to both children
        sca, _ = Phylo.score(str(self.childa), new_seq)
        scb, _ = Phylo.score(str(self.childb), new_seq)

        # Insert into the child with better score
        if sca > scb:
            self.childa = self.childa.insert(new_seq)
            self.childa.parent = self
        else:
            self.childb = self.childb.insert(new_seq)
            self.childb.parent = self
        # Update this node's consensus based on all leaves beneath it
        
        return self

    def _update_consensus(self, new_seq):
        cons = bio.merge_with_gaps(self.consensus, new_seq)
        return cons

    def recompute_consensus(self):
        
        refa = self.childa.recompute_consensus()
        refb = self.childb.recompute_consensus()
        
        score, cons = Phylo.score(refa, refb, ref_is_consensus = True)
        self.consensus = cons
        
        return cons

    def print_tree(self, prefix="", is_last=True, show_length = 32):
        connector = "└─ " if is_last else "├─ "
        extension = "   " if is_last else "│  "
        
        info = f"{prefix}{connector}Node: " 
        extra = f"depth={self.depth}, {len(self.get_all_leaves())} leaves"
        row = self._print_frm.format(info=info, seq = self.consensus[:show_length], extra = extra)
        print(row)

        # Print children
        self.childa.print_tree(prefix + extension, is_last=False, show_length=show_length)
        self.childb.print_tree(prefix + extension, is_last=True, show_length=show_length)

    def __repr__(self):
        return f"PhyloNode(consensus={self.consensus[:32]}..., depth={self.depth}, leaves={len(self.get_all_leaves())})"

    def __str__(self):
        return self.consensus
    
class Phylo:

    def __init__(self, refa, refb):
        sc, cons = self.score(refa, refb)
        lrefa = PhyloLeaf(refa, score=sc)
        lrefb = PhyloLeaf(refb, score=sc)
        self.root = PhyloNode(lrefa, lrefb, parent=self, consensus=cons)

    @property
    def depth(self):
        return -1  # Root is at depth -1 so first nodes are at depth 0

    @property
    def leaves(self):
        return self.root.get_all_leaves()

    def insert(self, new_seq):
        self.root = self.root.insert(new_seq)

    @classmethod
    def score(cls, ref, seq, ref_is_consensus=False):
        if ref_is_consensus:
            algn = align.align_sequences(str(ref), str(seq))[0]
        algn = align.align_sequences(str(ref), str(seq))[0]
        return algn.score, algn.consensus

    def print(self, show_length = 32):
        print(f"Phylo Tree ({len(self.leaves)} leaves)")
        self.root.print_tree(prefix="", is_last=True, show_length = show_length)

    def recompute_consensus(self):
        _ = self.root.recompute_consensus()

    def __repr__(self):
        return f"Phylo({len(self.leaves)} leaves)"
    

