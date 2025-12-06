

from typing import Dict, List, Any, Union, Optional
from dataclasses import dataclass

from ggene.genomemanager import GenomeManager
from ggene.seqs import bio, align, process

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
    
    @property
    def context(self):
        return len(self.upstream)

    def __len__(self):
        return len(self.repeat)

class Repeat:
    
    def __init__(self, name, type, base_motif = "", motif = ""):
        
        # if not base_motif:
        #     base_motif = motif = name
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
