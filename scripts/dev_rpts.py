
from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np

from Bio import Align

from ggene import DATA_DIR
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw

from ggene.seqs import bio, process, find, align

@dataclass
class Motif:
    motif:str
    instances:List[int]
    count:float
    length:float
    
    def __eq__(self, other):
        
        if len(self.motif) != len(other.motif):
            return False
            
        ext = self.motif * 2
        if other.motif in ext:
            return True
        elif bio.reverse_complement(other.motif) in ext:
            return True
        return False
    
    def combine(self, other):
        return Motif(self.motif, self.instances+other.instances, self.count+other.count, self.length+other.length)
        
    def __hash__(self):
        return hash(self.motif)


def load_genome():
    gm = GenomeManager()
    rpts = unified_stream.BEDStream(DATA_DIR / "repeatmasker/repeats.sorted.bed.gz", feature_type = "repeats")
    gm.annotations.add_source("repeats", rpts)
    return gm, rpts

def extract_repeats(gm, chr, repeat_name = "", repeat_type = "", motif = "", start = None, end = None):
    
    rpts = []
    for rpt in gm.annotations.streams.get("repeats").stream(chr, start=start, end=end):
        
        if repeat_name and rpt.name != repeat_name:
            continue
        if repeat_type and rpt.attributes.get("type") != repeat_type:
            continue
        if motif and rpt.attributes.get("motif") != motif:
            continue
        
        seq = gm.annotations.get_sequence(chr, rpt.start, rpt.end)
        rpt.attributes["seq"] = seq
        rpts.append(rpt)
    
    return rpts

def get_motif_stats(gm, chr, min_length = 3, combine_rc = True, combine_perms = True):
    
    topk = 5
    hist_bins = 128
    motif_count = {}
    motif_len = {}
    motif_insts = {}
    perm_groups = {}
    
    rpts = gm.annotations.streams.get("repeats")
    for rpt in rpts.stream(chr, start=0):
        
        motif = rpt.attributes.get("motif")
        if not motif:
            continue
        
        mtflen = len(motif)
        rptlen = rpt.end - rpt.start
        if mtflen < min_length:
            continue
        
        if combine_rc:
            motif_rc = bio.reverse_complement(motif)
        else:
            motif_rc = ""
        
        proto_motif = motif
        if combine_perms:
            for seen_mtf in perm_groups.keys():
                dbl = seen_mtf + seen_mtf
                if motif in dbl or motif_rc in dbl:
                    proto_motif = seen_mtf
                    break
        else:
            if motif_rc in perm_groups:
                proto_motif = motif_rc
        
        if not motif in motif_insts:
            motif_count[proto_motif] = 0
            motif_len[proto_motif] = 0
            motif_insts[proto_motif] = []
            perm_groups[proto_motif] = set()
        
        motif_count[proto_motif] += rptlen / len(motif)
        motif_len[proto_motif] += rptlen
        motif_insts[proto_motif].append(rpt.start)
        perm_groups[proto_motif].add(motif)
    
    top_cts = sorted(list(motif_count.keys()), key = lambda k:-motif_count[k])
    top_lens = sorted(list(motif_len.keys()), key = lambda k:-motif_len[k])
    
    print(f"found {len(motif_count)} motifs on chr{chr}")
    
    print("top by count:")
    for i in range(topk):
        mtf = top_cts[i]
        print(f"{mtf}: count = {motif_count[mtf]:0.1f}, len = {motif_len[mtf]:0.1f}, instances = {len(motif_insts[mtf])}")
        print(f"permutation group: {",".join(perm_groups[mtf])}")
        disthist, bins = np.histogram(motif_insts[mtf], bins = hist_bins)
        sctxt = draw.scalar_to_text_nb(disthist, minval = 0)
        for r in sctxt:
            print(r)
        print()
    
    print("top by length:")
    for i in range(topk):
        mtf = top_lens[i]
        print(f"{mtf}: count = {motif_count[mtf]:0.1f}, len = {motif_len[mtf]:0.1f}, instances = {len(motif_insts[mtf])}")
        print(f"permutation group: {",".join(perm_groups[mtf])}")
        disthist, bins = np.histogram(motif_insts[mtf], bins = hist_bins)
        sctxt = draw.scalar_to_text_nb(disthist, minval = 0)
        for r in sctxt:
            print(r)
        print()
    
    # print("permutation groups:")
    # for proto, perms in perm_groups.keys():
        
    #     pass
    
    return motif_count, motif_len, motif_insts

def get_motif_length_stats(gm):
    
    motif_count = {}
    motif_len = {}
    
    rpts = gm.annotations.streams.get("repeats")
    for rpt in rpts.stream(chr, start=0):
        
        motif = rpt.attributes.get("motif")
        if not motif:
            continue
        
        
        
        pass
    
    pass

def view_repeats(gm, chr, start, motifs_only = False):
    rpts = gm.annotations.streams.get("repeats")
    for rpt in rpts.stream(chr, start=start):
        if motifs_only and not rpt.attributes.get("motif"):
            continue
        print(rpt)
        input()



def main():
    
    gm, rpts = load_genome()
    
    # view_repeats(gm, "1", 10e6, motifs_only = True)
    
    # test_motif = "TCAGCC"
    # rpts = extract_repeats(gm, "1", motif = test_motif)
    # print(f"found {len(rpts)} repeats")
    
    # for r in rpts:
    #     print(r)
    
    testchr = "1"
    
    # print("combined perms:")
    # get_motif_stats(gm, testchr, min_length = 4)
    # print()
    
    print("no combined perms:")
    get_motif_stats(gm, testchr, min_length = 5, combine_perms = False)
    
    pass

if __name__=="__main__":
    main()


