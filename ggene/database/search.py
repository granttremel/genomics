
from typing import TYPE_CHECKING
import re

from ggene import motifs
from ggene.seqs.find import consensus_to_re

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager



class GenomeSearch:
    
    def __init__(self, genome_manager:'GenomeManager'):
        self.gm = genome_manager

    
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
        
        
    