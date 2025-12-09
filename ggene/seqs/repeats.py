

from typing import Dict, List, Any, Union, Optional
from dataclasses import dataclass

import re
import regex
import numpy as np


from ggene.genomemanager import GenomeManager
from ggene.seqs import bio, align, process, compare

def hamming_distance_algn(motifa, motifb):
    if len(motifb) < len(motifa):
        motifa, motifb = motifb, motifa
    
    algn = align.align_sequences(motifa, motifb)[0]
    
    ham = algn.length - algn.identities
    
    return ham

def hamming_distance_re(motifa, motifbs, max_err = 1, min_err = 1):
    
    ptrn_str = "(%s){e<=%s}" % (str(motifa), str(max_err))
    motifa_ptrn = regex.compile(ptrn_str, regex.BESTMATCH)
    
    res = []
    for motifb in motifbs:
        # match = regex.search(motifa_ptrn, str(motifb))
        # match = regex.match(motifa_ptrn, str(motifb))
        match = regex.fullmatch(motifa_ptrn, str(motifb))
        if match and sum(match.fuzzy_counts) >= min_err:
            s, i, d = match.fuzzy_counts
            res.append((True, s, i, d))
        else:
            res.append((False, -1, -1, -1))
    return res
    

@dataclass
class Loc:
    chr:str
    start:int
    end:int
    strand:str
    
    @classmethod
    def from_feature(cls, feat):
        
        return cls(
            feat.get("chr", feat.get("chrom", "")),
            feat.get("start", 0),
            feat.get("end"),
            feat.get("strand","+")
        )
    
    def __hash__(self):
        return hash(self.chr) + hash(self.start) + hash(self.end) + hash(self.strand)

    def __eq__(self, other:Union['Loc',Any]):
        return self.chr == other.chr and self.start==other.start and self.end==other.end and self.strand==other.strand

@dataclass
class RepeatSeq:
    location:Loc
    repeat:str
    count:int
    upstream:str = ""
    downstream:str = ""
    gene:str = ""
    raw_repeat:Dict[str, Any] = None
    
    @property
    def context(self):
        return len(self.upstream)

    def __len__(self):
        return len(self.repeat)
    
    @classmethod
    def get_instance(cls, gm:GenomeManager, rpt_feat, context_sz = 128):
        
        loc = Loc.from_feature(rpt_feat)
        rpt_seq = gm.get_sequence(loc.chr, loc.start, loc.end)
        mtf_len = len(rpt_feat.get("motif"))
        count = len(rpt_seq)//mtf_len if mtf_len else -1
        upstr = gm.get_sequence(loc.chr, loc.start - context_sz, loc.start)
        downstr = gm.get_sequence(loc.chr, loc.end, loc.end + context_sz)
        
        if rpt_feat.get("strand","+") in ["C","-"]:
            
            rpt_seq = bio.reverse_complement(rpt_seq)
            upstr = bio.reverse_complement(upstr)
            downstr = bio.reverse_complement(downstr)
        
        gene = ""
        for g in gm.annotations.query_range(loc.chr, loc.start, loc.end):
            gn = g.get("info",{}).get("gene_name","")
            if gn:
                gene = gn
        
        rpt_inst = cls(loc, rpt_seq, count, upstr, downstr, gene=gene, raw_repeat = rpt_feat)
        
        return rpt_inst

class Repeat:
    
    def __init__(self, name, type, base_motif = "", motif = ""):
        
        if not motif:
            motif, do_rc = self.normalize_motif(base_motif)
        self.motif = motif
        
        self.name = name
        self.type = type
        self.base_motif = base_motif
        
        self.locs:List[Loc] = []
        self.instances:Dict[Loc, RepeatSeq] = {}
        self.neighbors:List['Repeat'] = []
        
    @property
    def motif_len(self):
        return len(self.motif)
    
    @property
    def num_instances(self):
        return len(self.instances)
    
    def add_instances(self, repeat_features: List[Dict[str,Any]], context_sz = 0):
        for rptf in repeat_features:
            self.add_instance(rptf, context_sz = context_sz)
    
    def add_instance(self, repeat_feat:Dict[str,Any], rpt_inst:RepeatSeq):
        
        if rpt_inst.location in self.locs:
            return
        
        rpt_motif, _ = self.normalize_motif(repeat_feat.get("motif"))
        
        if rpt_motif == self.motif:
            loc = Loc.from_feature(repeat_feat)
            self.instances[loc] = rpt_inst
            return rpt_inst
        else:
            return None
    
    def get_instance(self, gm:GenomeManager, rpt_feat, context_sz = 128):
        
        loc = Loc.from_feature(rpt_feat)
        rpt_seq = gm.get_sequence(loc.chr, loc.start, loc.end)
        mtf_len = len(rpt_feat.get("motif"))
        count = len(rpt_seq)//mtf_len if mtf_len else -1
        upstr = gm.get_sequence(loc.chr, loc.start - context_sz, loc.start)
        downstr = gm.get_sequence(loc.chr, loc.end, loc.end + context_sz)
        
        rpt_inst = RepeatSeq(loc, rpt_seq, count, upstr, downstr)
        
        if rpt_inst.location in self.locs:
            return
        else:
            self.add_instance(rpt_feat, rpt_inst)
        
    @classmethod
    def check_neighbor(cls, rpta:'Repeat', rptb:'Repeat'):
        
        if abs(len(rpta.motif) - len(rptb.motif)) > 1:
            return False
        
        score = rpta.compare_motif(rptb)
        
        if score >= rpta.motif_len - 1:
            return True
        else:
            return False
        
    def get_hamming_distance(self, other_rpt:'Repeat'):
        return self.hamming_distance(self.motif, other_rpt.motif)
    
    @classmethod
    def hamming_distance(cls, mtfa:str, mtfb:str):
        
        if len(mtfb) < len(mtfa):
            mtfa, mtfb = mtfb, mtfa
        
        algn = align.align_sequences(mtfa, mtfb)[0]
        gaps, ids, mms = algn.gaps, algn.identities, algn.mismatches
        
        ham = algn.length - algn.identities
        
        return ham
    
    def compare_motif(self, other_rpt:'Repeat'):
        sc1 = align.score_sequences(self.motif, 2*other_rpt.motif)
        sc2 = align.score_sequences(other_rpt.motif, 2*self.motif)
        return max(sc1, sc2)
    
    def get_phase(self, mtf_proto = "", rpt_seq = ""):
        
        if not mtf_proto:
            mtf_proto = self.motif
        
        mlen = len(mtf_proto)
        nrpts = len(rpt_seq)//len(mtf_proto)
        prpt = (mtf_proto * (nrpts + 1))[:len(rpt_seq)]
        corr, rcorr = process.correlate(prpt, rpt_seq, scale = 2*mlen)
        maxshift = (corr.index(corr==max(corr))) % len(mtf_proto)
        # ScalarPlot(corr).show()
        # print(f"maxshfit: {maxshift}")
        return maxshift
    
    @classmethod
    def from_feature(cls, rpt_feat:Dict[str,Any]):
        
        name = rpt_feat.get("name","")
        type = rpt_feat.get("type","")
        base_motif = rpt_feat.get("motif","")
        
        motif, do_rc = cls.normalize_motif(base_motif)
        
        return cls(name, type, base_motif, motif=motif)
    
    @classmethod
    def normalize_motif(cls, motif):
        
        motif_nrm, _, _, ind = bio.get_least_seq_index(motif)
        
        if motif_nrm is None:
            return None, None
        
        if (ind//2)%2 == 1:
            do_rc = True
        else:
            do_rc = False
            
        return motif_nrm, do_rc
        
    def __repr__(self):
        return f"Repeat({self.motif}, name={self.name}, type={self.type}, num_instances={len(self.instances)})"
    
class RepeatLibrary:
    
    def __init__(self, motifs_only = True):
        
        self.types:List[str] = []
        self.motifs:List[str] = []
        self.repeats:Dict[str,Repeat] = {}
        self.neighbors = []
        
        self.motifs_only = motifs_only
    
    @property
    def num_repeats(self):
        return len(self.repeats)
    
    @property
    def num_instances(self):
        return sum([rpt.num_instances for rpt in self.repeats.values()])
    
    def add_repeat(self, rpt:Repeat):
        
        if self.motifs_only and not rpt.motif:
            return
        
        if not rpt.type in self.types:
            self.types.append(rpt.type)
        
        if not rpt.motif in self.motifs:
            self.motifs.append(rpt.motif)
        
        if rpt.motif and not rpt.motif in self.repeats:
            self.repeats[rpt.motif] = rpt
        
    
    def add_instance(self, rpt_feat:Dict[str, Any], rpt_inst: RepeatSeq):
        
        rpt = Repeat.from_feature(rpt_feat)
        if rpt is None:
            return
        
        if self.motifs_only and not rpt.motif:
            return
        
        if not rpt in self:
            self.add_repeat(rpt)    
        
        lib_rpt = self.repeats.get(rpt.motif)
        rpt_inst = lib_rpt.add_instance(rpt_feat, rpt_inst)
        
        return rpt_inst
    
    def get_instance(self, gm:GenomeManager, rpt_feat, context_sz = 128):
        
        loc = Loc.from_feature(rpt_feat)
        rpt_seq = gm.get_sequence(loc.chr, loc.start, loc.end)
        mtf_len = len(rpt_feat.get("motif"))
        count = len(rpt_seq)//mtf_len if mtf_len else -1
        upstr = gm.get_sequence(loc.chr, loc.start - context_sz, loc.start)
        downstr = gm.get_sequence(loc.chr, loc.end, loc.end + context_sz)
        
        rpt_inst = RepeatSeq(loc, rpt_seq, count, upstr, downstr)
        rpt_inst = self.add_instance(rpt_feat, rpt_inst)
    
        return rpt_inst
    
    def build_neighborhood(self, hd_delta = 1, len_diff = 1):
        
        for i, mtfi in enumerate(self.motifs):
            for j, mtfj in enumerate(self.motifs[i+1:]):
                
                if abs(len(mtfi) - len(mtfj)) > len_diff:
                    continue
                
                hd = Repeat.hamming_distance(mtfi, mtfj)
            
                if hd > 0 and hd <= hd_delta:
                    self.neighbors.append((mtfi, mtfj))
        
        return self.neighbors
    
    def get_top(self, topk = 10):
        return sorted(self.repeats.values(), key=lambda k:-k.num_instances)[:topk]
    
    def __contains__(self, rpt:Repeat):
        return rpt.motif in self.repeats


def fold_repeat(motif, repeat_seq, num_repeats = None):
    
    if not num_repeats:
        num_repeats = len(repeat_seq) // len(motif)
    
    corrs = []
    
    for i in range(len(repeat_seq) - len(motif)):
        compres, _ = compare.count_matching(motif, repeat_seq[i:i+len(motif)])
        corrs.append(compres)
        
    # from ggene import scalar_plot
    # scalar_plot.ScalarPlot(corrs, minval = 0, add_range = True).show()
    
    top_corrs = sorted(enumerate(corrs), key = lambda k:-k[1])
    
    folded = []
    
    for n in range(num_repeats):
        
        ind, c = top_corrs[n]
        
        rpt_seq = repeat_seq[ind: ind+len(motif)]
        
        folded.append(rpt_seq)
  
    return folded


########################################## for tree

class Motif:
    
    def __init__(self, motif):
        
        self.motif = motif
        self.motif_len = len(motif)
        self.parents = set()
        self.peers = set()
        self.children = set()
        
        self.motif_ptrn = re.compile(f'({motif}){{e==1}}')
    
    @property
    def is_leaf(self):
        return len(self.children)==0
    
    def get_hamming_algn(self, other_motif):
        return hamming_distance_algn(self.motif, other_motif.motif)
    
    def get_hamming_re(self, other_motifs):
        return hamming_distance_re(self, other_motifs)
    
    def join(self, other_motif):
        if self.motif_len < other_motif.motif_len:
            self.add_child(other_motif)
            other_motif.add_parent(self)
        elif self.motif_len > other_motif.motif_len:
            self.add_parent(other_motif)
            other_motif.add_child(self)
        else:
            self.add_peer(other_motif)
            other_motif.add_peer(self)
    
    def add_parent(self, other_motif):
        self.parents.add(other_motif)
    
    def add_peer(self, other_motif):
        self.peers.add(other_motif)
    
    def add_child(self, other_motif):
        self.children.add(other_motif)
    
    def print_children(self, print_all = True, tab = 0):
        
        if tab > 0 and not print_all:
            return
        
        print(f"{"  "*tab} {self.motif} ({len(self.parents)} parents, {len(self.children)} children, {len(self.peers)} peers)")
        for child in self.children:
            child.print_children(print_all=print_all, tab = tab+1)
            
    
    def __hash__(self):
        return hash(self.motif)
    
    def __len__(self):
        return self.motif_len
    
    def __str__(self):
        return self.motif

class MotifNetwork:
    
    def __init__(self):
        self.motifs = []
        self.motif_lens = {}
        
    def add_to_tree(self, new_motif):
        
        self.motifs.append(new_motif)
        if not new_motif.motif_len in self.motif_lens:
            self.motif_lens[new_motif.motif_len] = []
        self.motif_lens[new_motif.motif_len].append(new_motif)
        
        parent_motifs = self.motif_lens.get(new_motif.motif_len - 1,[])
        res = new_motif.get_hamming_re(parent_motifs)
        
        for mtf, (r,s,i,d) in zip(parent_motifs, res):
            if not r:
                continue
            mtf.join(new_motif)

    
    def get_roots(self):
        roots = []
        for mtf in self.motifs:
            if len(mtf.parents) == 0:
                roots.append(mtf)
        return roots
    
    def join_peers(self):
        
        for i, mtfi in enumerate(self.motifs):
            
            peer_motifs = self.motif_lens.get(mtfi.motif_len, [])
            res = mtfi.get_hamming_re(peer_motifs)
            for mtf, (r,s,i,d) in zip(peer_motifs, res):
                if not r:
                    continue
                mtf.join(mtfi)
    
    def print(self):
        
        roots = self.get_roots()
        
        for mtf in roots:
            mtf.print_children(print_all = True, tab = 0)
            
    def print_stats(self):
        
        npars = []
        npeers = []
        nchilds = []
        
        for mtf in self.motifs:
            npars.append(len(mtf.parents))
            npeers.append(len(mtf.peers))
            nchilds.append(len(mtf.children))
        
        num_edges = (sum(npars) + sum(npeers) + sum(nchilds))//2
        
        print(f"Network with {len(self.motifs)} motifs, {num_edges} edges")
        print(f"Parents: {sum(npars)} parent nodes with {sum(npars)} connections, mean={np.mean(npars):0.3f}, sd={np.std(npars):0.3f}, min={min(npars)}, max={max(npars)}")
        print(f"Peers: {sum(npeers)} peer nodes with {sum(npeers)} connections, mean={np.mean(npeers):0.3f}, sd={np.std(npeers):0.3f}, min={min(npeers)}, max={max(npeers)}")
        print(f"Children: {sum(nchilds)} child nodes with {sum(nchilds)} connections, mean={np.mean(nchilds):0.3f}, sd={np.std(nchilds):0.3f}, min={min(nchilds)}, max={max(nchilds)}")
        