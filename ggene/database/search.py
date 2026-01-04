
from typing import TYPE_CHECKING
import re

from ggene import motifs
from ggene.seqs import bio
from ggene.seqs.find import consensus_to_re

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager



class GenomeSearch:
    
    def __init__(self, genome_manager:'GenomeManager'):
        self.gm = genome_manager

    
    def simple_sequence_search(self, sseq, chromes, start, end, batch_len = 2**16,  personal = False):
        rcsseq = bio.reverse_complement(sseq)
    
        for chrom in chromes:
            print(f"starting chrom {chrom}")
            seq = "hi"
            pos = int(start)
            while seq:
                
                if personal:
                    seq = gm.annotations.get_personal_sequence(chrom, pos, pos + batch_len)
                else:
                    seq = gm.annotations.get_sequence(chrom, pos, pos + batch_len)
                
                n_insts = seq.count(sseq)
                n_rcinsts = seq.count(rcsseq)
                n_tot = n_insts + n_rcinsts
                
                if n_tot > 0:
                    
                    ind = None
                    if n_insts > 0:
                        ind = seq.index(sseq)
                    rcind = None
                    if n_rcinsts> 0:
                        rcind = seq.index(rcsseq)
                    
                    print(pos, n_tot, f"FWD: {n_insts}, RC: {n_rcinsts}")
                    if ind:
                        print("FWD", pos+ind, seq[ind-16:ind+len(sseq)+16])
                        feats = gm.annotations.query_range(chrom, pos+ind, pos+ind+1)
                        
                        for f in feats:
                            if f.feature_type == "dfam_hit":
                                print(f)
                    if rcind:
                        print("RC ",pos+rcind, seq[rcind-16:rcind+len(sseq)+16])
                        feats = gm.annotations.query_range(chrom, pos+rcind, pos+rcind+1)
                        for f in feats:
                            if f.feature_type == "dfam_hit":
                                print(f)
                    
                    
                    if not feats:
                        print("no feats in range")
                    
                    input()
                
                pos += batch_len
                
                if end and pos > end:
                    break
            
    
    def feature_search(self, chrom, start, end, conditions):
        
        if chrom:
            chromes = [chrom]
        else:
            chromes = [chrom for chrom in self.gm.iter_chromes()]
        
        rejected = []
        
        q_any = conditions.pop("any", None)
        
        for chrom in chromes:
            print(f"searching on chromosome {chrom}")
            
            for f in self.gm.annotations.stream_all(chrom, start = start, end = end):
                
                if q_any:
                    for k, v in f.to_dict().items():
                        
                        if isinstance(v,str) and q_any in v:
                            # return f
                            rejected.append(f)
                            yield f
                
                for c, v in conditions.items():
                    
                    if getattr(f, c) != v:
                        break
                    
                    # return f
                    rejected.append(f)
                    yield f
        
        # return None
    
    def consensus_search(self, chrom, start, end, cons):
        
        ptrn = re.compile(consensus_to_re(cons))
        
        for seq in self.iter_seq(chrom, start, end, chunksz = 1024):
            
            for res in re.finditer(ptrn, seq):
                
                if res:
                    yield res
            
        
        
    
    def iter_seq(self, chrom, start, end, chunksz):
        
        while start < end:
            
            chunk_end = min(end, start+chunksz)
            seq = self.gm.get_sequence(chrom, start, chunk_end)
            
            start += chunksz
            
            yield seq
        
        
    