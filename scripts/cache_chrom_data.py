

import os
import numpy as np
from ggene import seqs, draw
import random
from ggene.seqs import bio, find, process, heal, align, vocab, lambdas
from ggene.seqs.bio import reverse_complement


from ggene.database.genome_manager import GenomeManager
from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.draw.scalar_plot import ScalarPlot
from ggene.processing.chrome_mapper import ChromeMapper
from ggene.database.cache import GenomeCache

from ggene.database.annotations import chr_lens

def load_genome():
    return GenomeManager()

def map_chromosomes(gm, display_width = 256):
    
    seq_specs = ['gc','genes']
    
    max_len = max([chr_len for chr_len in chr_lens.values()])
    
    scale = display_width / max_len
    
    chunksz = int(max_len / display_width)+1
    
    
    sg, fg = ChromeMapper.get_generators(gm, seq_specs)
    
    cm = ChromeMapper(sg, fg)
    
    for chrom in gm.iter_chromes():
        chr_len = chr_lens.get(chrom, -1)
        
        qts, _ = cm.get_chromosomal_quantities(chrom, seq_specs, chunksz, start=1e6, length = chr_len, needs_feats = ['gene'])

        print(f"Chromosome {chrom}. gc | genes")
        sc1 = ScalarPlot(qts[0], add_range = True, minval = 0.3, maxval = 0.6)
        sc2 = ScalarPlot(qts[1], add_range = True, minval = 0, xmin = 0, xmax = chr_len, ruler = True, num_labels = 5, num_ticks = 0, num_minor_ticks = 0, fg_color = 124)
        ScalarPlot.show_paired(sc1, sc2)
        
        
    
def main():
    
    gm = load_genome()

    gc = GenomeCache(gm, base_resolution = 4096, dtype = 'float32')

    seq_specs = []
    # seq_specs += ['gc']
    # seq_specs += ['at', 'ag', 'ac', 'cpg', 'cpg_to_gc', 'polya', 'polyt', 'polyg', 'polyc', 'polyy', 'polyr', 'polys', 'polyw']
    # seq_specs += ['features', 'genes', 'exons', 'motifs', 'pseudo', 'lncRNA', 'ncRNA', 'simple_repeats', ]
    # seq_specs += [  'Alu', 'AluJ', 'AluY', 'AluS', 'Line1', 'L1M', 'L1P', 'L2', 'LTR', 'MIR', 'SVA', ]
    # seq_specs += ['TTCTT']
    # seq_specs += [ 'TTCTT', 'TATAAT', 'TATATC', 'TTGA', 'AATTAAG', 'ATGT', 'GAGGAT', 'CCCTGC']
    # seq_specs += ['protein_coding', 'lncrna', 'cds_len', 'cds_pct']
    # seq_specs += ["L1HS","L1PA"]
    seq_specs = ["miRNA"]

    gc.cache_all_chromes(seq_specs, overwrite = True)
    




    pass

if __name__=="__main__":
    main()

