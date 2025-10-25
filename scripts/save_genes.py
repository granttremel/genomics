
import os

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def get_genes_range(geneman, chr, start, end):
    
    annos = geneman.get_all_annotations(chr, start, end, include_motifs = False)
    
    genes = [f for f in annos if f.feature_type == "gene"]
    
    fts = [f.feature_type for f in annos]
    print("feature types:",list(set(fts)))
    return genes

def get_genes_pos(geneman, chr, pos):
    
    genes = geneman.get_features_at_position(chr, pos, include_types = ["gene"])
    return genes

def save_gene(geneman, gene_name, chr):
    
    gene = geneman.assemble_gene(gene_name, chr)
    
    if not gene:
        print(f"unable to save gene {gene_name}")

def list_known_genes():
    
    libdir = "./data/library"
    
    genefiles = os.listdir(libdir)
    gns = [gf.removesuffix('.json') for gf in genefiles]
    return gns

def main():
    

    genes = [
        # 'C2orf80',
        # 'STMN3',
        "Metazoa_SRP"
    ]
    
    gene_idx = [
    ]
    
    gm = load_genome()
    known_genes = list_known_genes()
    
    offset = 184598366 - 185463093 # my dataset - GRCh37
    
    genes_dedup = [g for g in genes if not g in known_genes]
    print(f"after deduplication, finding {len(genes_dedup)} genes by name")
    gm.find_and_assemble_genes(genes_dedup)
    
    for gene_name, gchr, gstart, gend, refid in gene_idx:
        
        if refid== "GRCh37":
            gstart += offset
            gend += offset
        
        # idx = gm._make_index(gchr, gstart, gend)
        gpos = (gstart + gend)//2
        idx = gm._make_index(gchr, gpos, gpos+1)
        # grch37_idx = gm._make_index(gchr, gstart - offset, gend - offset)
        
        print(f"probing genes at {idx}")
        genes = get_genes_range(gm, gchr, gstart, gend)
        # genes = get_genes_pos(gm, gchr, gpos)
        print(f"found {len(genes)} genes at {idx}")
        
        for gene_feat in genes:
            _gene_name = gene_feat.name
            if _gene_name and not _gene_name in known_genes:
                save_gene(gm, _gene_name, gchr)

if __name__=="__main__":
    main()

