
import os

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def collect_vars(geneman, chr, start, end):
    raw_vars = list(geneman.annotations.streams["variants"].vcf(geneman._make_index(chr, start, end)))
    return raw_vars

def print_var(var_id, var, gene_name):
    
    atts = ["REF", "ALT", "QUAL", "gt_bases"]
    
    print(f"variant data for {var_id} on gene {gene_name}")
    
    for a in atts:
        print(a, getattr(var, a))
        
    if var.num_het > 0:
        print("heterozygous")
    elif var.num_hom_ref > 0:
        print("homozygous for reference allele")
    elif var.num_hom_alt > 0:
        print("homozygous for alternate allele")

def print_all_var(var):
    
    atts = dir(var)
    for a in atts:
        if a.startswith('_'):
            continue
        print(a, getattr(var, a))
    
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
    
    
    var_data = [
        # ("rs3131296",6, 32205216, "NOTCH4",""), # NOTCH4
        # ("rs12807809",11, 124736389, "",""),
        # ("rs9960767",18, 55487771, "TCF4",""), # TCF4
    ]
    
    offset = 184598366 - 185463093 # my dataset - GRCh37
    # print(offset) # -864727
    
    search_window = 10
    
    gm = load_genome()
    known_genes = list_known_genes()
    
    for var_id, vchr, vpos, gene_name, refid in var_data:
        
        search_start = vpos - search_window
        search_end = vpos + search_window
        
        if refid== "GRCh37":
            search_start += offset
            search_end += offset
        
        idx = gm._make_index(vchr, search_start, search_end)
        grch37_idx = gm._make_index(vchr, search_start - offset, search_end - offset)
        
        vars = collect_vars(gm, vchr, search_start, search_end)
        
        print(f"variants identified for {var_id} on gene {gene_name}: {len(vars)}")
        print(f"my index {idx}")
        print(f"GRCh37 index {grch37_idx}")
        
        if len(vars) > 0:
            for v in vars:
                print_var(var_id, v, gene_name)
        
            for v in vars:
                print_all_var(v)
        
        if gene_name and not gene_name in known_genes:
            save_gene(gm, gene_name, vchr)

if __name__=="__main__":
    main()

