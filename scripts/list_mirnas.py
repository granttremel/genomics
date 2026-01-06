import numpy as np

from ggene.config import other_paths
from ggene.database.genome_manager import GenomeManager
from ggene.database.cache import GenomeCache
from ggene.database.ufeature import UFeat

from ggene.genome import genes, ncrna

from ggene.seqs import process, bio, vocab, lambdas, align

from ggene import draw
from ggene.draw.scalar_plot import ScalarPlot
from ggene.draw.colored_sequence import ColoredSequence, StyledColor


def load_genome():
    return GenomeManager(skip_clinvar = True)

def list_mirnas(gm, chrom, limit = 128):
    
    gs = gm.annotations.streams.get("genes")
    
    feats = []
    n = 0
    for f in gs.stream(chrom, start = int(1e6), end = None):
        
        gn = f.gene_name if f.gene_name else ""
        bt = f.gene_biotype if f.gene_biotype else ""
        
        if f.gene_name and "mirna" in f.gene_name.lower():
            feats, n = process_mirna(f, feats=feats,n=n)
        elif f.gene_biotype and "mirna" in f.gene_biotype.lower():
            feats, n = process_mirna(f, feats=feats,n=n)
        else:
            pass
        
        if limit and n > limit:
            break
    
    gene_org = {}
    feats_org = {}
    num_in_exon = 0
    num_in_3p = 0
    num_in_5p = 0
    
    for f in feats:
        pg, gbt, in_exon, in_3putr, in_5putr  = show_mirna(f, gm=gm)
        print()
        if not gbt in feats_org:
            feats_org[gbt] = []
        feats_org[gbt].append(f)
        if gbt:
            num_in_exon += int(in_exon)
            num_in_3p += int(in_3putr)
            num_in_5p += int(in_5putr)
        if pg:
            pgn = pg.name if pg.name else pg.id
            if not pgn in gene_org:
                gene_org[pgn] = []
            gene_org[pgn].append(f)
            
    for gbt, fs in feats_org.items():
        print(f"biotype {gbt} contains {len(fs)} miRNA genes")
    print()
    
    for pgn, fs in gene_org.items():
        print(f"gene {pgn} contains {len(fs)} miRNA genes")
    print()
    
    print(f"MIR genes found in: exons {num_in_exon}, 3p UTR {num_in_3p}, 5p UTR {num_in_5p}")
    
    return feats

def process_mirna(f, feats = [], n = 0):
    
    if f.gene_name != f.name:
        return feats, n
    
    feats.append(f)
    n+=1
    
    return feats, n

def show_mirna(f, gm = None):
    
    in_ex = in_3putr = in_5putr = False
    seq = ""
    bt = "lone"
    
    ctx = 4
    if gm:
        pg, bt, in_ex, in_3putr, in_5putr = get_parent_gene(gm.annotations.streams.get("genes"), f)
        fseq = gm.get_sequence(f.chrom, f.start-ctx, f.end+ctx)
        seq = fseq[ctx:-ctx-1]
        rcseq = bio.reverse_complement(seq)
    else:
        pg = None
    
    print(f)
    mc = StyledColor.from_8bit(fg=0)
    mmc = StyledColor.from_8bit(fg=2)
    
    if pg:
        print(f"  parent gene: {pg.name} ({pg.feature_type}, {pg.gene_biotype})")
        print(f"    in exon: {in_ex}, in 3putr: {in_3putr}, in 5putr: {in_5putr}")
        # bt = pg.gene_biotype
    else:
        print("  no parent gene")
    
    if seq:
        algns = align.align_sequences(seq, rcseq)
        algn = next(iter(algns))
        algn.print(chunksz = 256, emph_indels = True, color_subs = True)
    
    return pg, bt, in_ex, in_3putr, in_5putr 

    

def get_parent_gene(gs, f):
    
    g = None
    bt = "lone"
    in_ex = False
    in_3putr = False
    in_5putr = False
    
    for f in gs.stream(f.chrom, f.start, f.end):
        
        # print(f)
        
        if f.feature_type == 'exon':
            in_ex = True
        elif f.feature_type == "three_prime_utr":
            in_3putr = True
        elif f.feature_type == "five_prime_utr":
            in_5putr = True
        elif f.feature_type in ['gene', 'lncRNA'] and not g:
            # print(f)
            if f.gene_biotype and f.gene_biotype == "protein_coding":
                g= f
                bt = f.gene_biotype
            elif f.feature_type == "lncRNA":
                g= f
                bt = f.feature_type
    
    return g, bt, in_ex, in_3putr, in_5putr

    

def main():
    
    gm = load_genome()
    

    fs = list_mirnas(gm, "1", limit = 4)




if __name__=="__main__":
    main()


"""
miR-376 uses ADAR to switch target
miR-142 uses it to prevent drosha cleavage


biotype lone contains 15 miRNA genes
biotype protein_coding contains 108 miRNA genes
biotype lncRNA contains 25 miRNA genes

MIR genes found in: exons 25, 3p UTR 5, 5p UTR 5

gene DNM3 contains 3 miRNA genes
"""