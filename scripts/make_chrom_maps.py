

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


def check_colors():
    
    teststr = "TACTATATACACCCACTAC"
    
    for c in range(16, 231):
        col = f'\x1b[38;5;{c}m'
        print(c, col + teststr + "\x1b[0m")

def get_random_seq(seq_len):
    return "".join("ATGC"[random.randint(0, 3)] for i in range(seq_len))
    
def find_nice_colors(seq_len = 256):
    
    # seq = list(range(-100, 0, 14)) + list(range(0, -50, 10)) + list(range(-50,25,9)) + list(range(25,100,18)) + list(range(100, 50, -19)) + [50]
    n = 0
    while True:
        seq = get_random_seq(seq_len, var = 4)
        offset = 24
        # fc = 170
        bc = 234
        # bc = random.randint(20, 230-offset)
        fc = random.randint(20,230)
        # fc = min(bc+offset, 230)
        
        if n == 0:
            fc = 236
            bc = 244
        
        res = draw.scalar_plot.scalar_to_text_nb(seq, bg_color = bc, fg_color = fc, bit_depth = 16)
            
        print(f"bg color = {bc}, fg color = {fc}")
        for r in res:
            print(r)
        print()
        n += 1
        
        if n%5 == 0:
            res=input()
            if 'n' in res.lower():
                break
    
    return bc

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
    
    # bg, fg = draw.get_color_scheme("test")

    gm = load_genome()

    gc = GenomeCache(gm)
    
    lambda_dict = {}
    
    for k in lambdas.lambda_map:
        
        needs = tuple(lambdas.needs_features(k))
        
        if not needs in lambda_dict:
            lambda_dict[needs] = []
        lambda_dict[needs].append(k)
        
    
    
        
    
    # seq_specs = lambda_dict.get(tuple(), [])
    old_seq_specs = ['gc', 'at', 'ag', 'ac', 'cpg', 'cpg_to_gc', 'polya', 'polyt', 'polyg', 'polyc', 'polyy', 'polyr', 'polys', 'polyw']
    # seq_specs = ['features', 'genes', 'exons', 'motifs', 'pseudo', 'lncRNA', 'ncRNA', 'simple_repeats', ]
    # seq_specs = [  'Alu', 'AluJ', 'AluY', 'AluS', 'Line1', 'L1M', 'L1P', 'L2', 'LTR', 'MIR', 'SVA', ]
    # seq_specs = [ 'TTCTT', 'TATAAT', 'TATATC', 'TTGA', 'AATTAAG', 'ATGT', 'GAGGAT', 'CCCTGC']
    # seq_specs = [ 'cds_len', 'cds_pct', ]
    
    # seq_specs = ['protein_coding', 'lncrna']
    
    # seq_specs += [ 'max_run', 'hammerhead_st1', 'tetra_repeats', 'penta_repeats', 'hexa_repeats']
    
    
    # gc.cache_all_chromes(seq_specs, overwrite = True)
    
    
    
    
    
    seq_specs = ['gc', 'at', 'ag', 'ac', 'cpg', 'cpg_to_gc', 'polya', 'polyt', 'polyg', 'polyc', 'polyy', 'polyr', 'polys', 'polyw']
    seq_specs += ['features', 'genes', 'exons', 'motifs', 'pseudo', 'lncRNA', 'ncRNA', 'simple_repeats', ]
    seq_specs += [  'Alu', 'AluJ', 'AluY', 'AluS', 'Line1', 'L1M', 'L1P', 'L2', 'LTR', 'MIR', 'SVA', ]
    seq_specs += [ 'TTCTT', 'TATAAT', 'TATATC', 'TTGA', 'AATTAAG', 'ATGT', 'GAGGAT', 'CCCTGC']
    seq_specs += [ 'cds_len', 'cds_pct']
    
    # to_del = ["simple_repeats", "motifs","TTCTT","TATAAT","TATATC","TTGA","AATTAG","ATGT","GAGGAT","CCCTGC"]
    # to_del = ["AATTAAG"]
    to_del = []
    
    for chrom in gm.iter_chromes():
        
        for seq_spec in seq_specs:
            
            res = gc.load_quantity_cache(chrom, seq_spec)
            
            if res is None:
                continue
            
            needs = lambdas.needs_features(seq_spec)
            
            print(f"chr{chrom} seq_spec {seq_spec} ({", ".join(needs)}) mean {np.mean(res)}, range ({np.min(res)}-{np.max(res)})")
            
    
    
    
    
    
    
    # data = gc.load_quantity_cache("Y","gc")
    
    # print(data.shape)
    
    # sc1 = ScalarPlot(data).show_chunks(chunksz = 256)
    

    # print(gm.annotations.streams.keys())

    # map_chromosomes(gm)


















    pass

if __name__=="__main__":
    main()

