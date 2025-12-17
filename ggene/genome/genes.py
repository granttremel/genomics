

from ggene import DEFAULT_LIBRARY
from ggene.genomemanager import GenomeManager

from ggene.unified_stream import UnifiedFeature

LIB_DIR = "ncache"

def list_feature_atts(gm:GenomeManager, att = "feature_type"):
    
    feat_vals = set()
    
    getter = None
    if att in UnifiedFeature.__dataclass_fields__:
        getter = lambda f:getattr(f, att)
    else:
        getter = lambda f:f.attributes.get(att,"")
    
    for chr in gm.iter_chromes():
        print(f"(chr{chr})")
        for feat in gm.annotations.stream_all(chr, start = 1e6):
            
            val = getter(feat)
            if not val:
                pass
            elif not val in feat_vals:
                feat_vals.add(val)
                print(val)
    
    return feat_vals

def get_first_gene(gm, chr, pos, upstream = False):
    
    step = 1024
    
    gene = None
    while not gene:
        
        feats = gm.get_features_at_position(chr, pos, include_types = ['gene'])
        
        if feats and any(feats):
            gene = feats[0]
        
        pos += (-1 if upstream else 1)*step
    
    return gene

def get_gene_names(gm, chr):
    gns = set()
    for g in gm.annotations.stream_by_types(['gene'], chr, 0, end=None):
        bt = g.attributes.get("gene_biotype","")
        
        if not g or not bt=="protein_coding" or not g.name:
            continue
        gns.add((g.name, chr))
    
    return gns
    
def get_all_gene_names(gm):
    
    gns = set()
    
    for chr in gm.iter_chromes():
        new_gns = get_gene_names(gm, chr)
        gns = set.union(gns, new_gns)
            
    return list(gns)

def save_gene_names(gene_names):
    fname = DEFAULT_LIBRARY / LIB_DIR / 'gene_names.txt'
    with open(fname.absolute(),"w+") as f:
        f.write("\n".join([f"{gn}, chr{chr}" for gn, chr in gene_names]))
    return gene_names

def load_gene_names(keep_chrom = False):
    
    if keep_chrom:
        return load_gene_names_chr()
    
    fname = DEFAULT_LIBRARY / LIB_DIR / 'gene_names.txt'
    all_gene_names = []
    with open(fname.absolute(),"r") as f:
        for line in f:
            gn, chr = line.strip().split(", ")
            all_gene_names.append(gn)
    return all_gene_names


def load_gene_names_chr():
    all_gene_names = {}
    fname = DEFAULT_LIBRARY / LIB_DIR / 'gene_names.txt'
    with open(fname.absolute(),"r") as f:
        for line in f:
            gn, chr = line.strip().split(", ")
            if not chr in all_gene_names:
                all_gene_names[chr] = []
            
            all_gene_names[chr].append(gn)
    
    return all_gene_names











