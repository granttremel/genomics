
import json

from ggene import DEFAULT_LIBRARY
from ggene.database.genomemanager import GenomeManager
from ggene.database.unified_stream import UFeature

LIB_DIR = "ncache/pseudogenes"

pg_types = [

    {
        "TC":'transcribed', 
        "TL":'translated'
    },    
    {
        "UP":'unprocessed',
        "P":'processed',
        "U":'unitary',
        "IG":'IG',
        "TR":'TR',
        "rRNA":'rRNA'
    },
    {
        "C":"C",
        "J":"J",
        "V":"V",
        "D":"D",
    }
]

pg_subclass = [
    "P", 
    "UP", 
    "U",
    "TC_P", 
    "TC_U", 
    "TC_UP",  
    "TL_P", 
    "IG"
    "IG_C", 
    "IG_J", 
    "IG_V", 
    "TR_V", 
    "TR_J", 
    "rRNA", 
]

def get_pseudogene_subclass(pg_biotype):

    pg_tup = []
    
    for pg_group in pg_types:
        for acro, opt in pg_group.items():
            if opt in pg_biotype:
                pg_tup.append(acro)
                break
    
    return "_".join(pg_tup)

def gather_all_pseudogenes(gm:GenomeManager):
    rna_genes = {}
    for chr in gm.iter_chromes():
        rna_genes[chr] = gather_pseudogenes(gm, chr)
    
    return rna_genes


def gather_pseudogenes(gm:GenomeManager, chrom):
    
    pg_types = []
    pgs ={}
        
    print(f"######### starting chr{chrom} ####### \n")
    
    for gene in gm.annotations.stream_by_types(['gene'], chrom, start = 2e6):
        
        gn = gene.get("gene_name", gene.get("name",""))
        bt = gene.get("attributes",{}).get("gene_biotype")
        if not "pseudo" in bt:
            continue
        
        if not bt in pg_types:
            pg_types.append(bt)
            pgs[bt] = []
        pgs[bt].append(gene)
    
    return pgs

def write_pseudogenes(pgs, ftag = ""):
    
    for pgt, genes in pgs.items():
        json_str = json.dumps({pgt:[g.to_dict() for g in genes]})
        
        for repstr in ["[",r"}},"]:
            json_str = json_str.replace(repstr, repstr + "\n")
            pass
        
        fname = DEFAULT_LIBRARY / LIB_DIR / (f"pseudogenes_{pgt}" + ("_" + ftag if ftag else "") + ".json")
        with open(fname.absolute(),"w+") as f:
            f.write(json_str)

def load_pseudogenes(pg_type, ftag = ""):
    
    fname = DEFAULT_LIBRARY  / LIB_DIR /  (f"pseudogenes_{pg_type}" + ("_" + ftag if ftag else "") + ".json")
    
    with open(fname.absolute(),'r') as f:
        pgenes = json.loads(f.read())
    
    pgenes = pgenes.get(pg_type,[])
    pg_feats = []
    for gd in pgenes:
        gd["feature_type"] = gd.pop("type", "")
        pg_feats.append(UFeature(**gd))
    
    return pg_feats

def characterize_pseudogenes(gm:GenomeManager):
    
    counts = {}
    
    for chr in gm.iter_chromes():
        for f in gm.annotations.stream_all(chr, start = 0):
            bt = f.attributes.get("gene_biotype")
            if not bt or not "pseudo" in bt:
                continue
            
            pg_tup = get_pseudogene_subclass()
            
            if not pg_tup in counts:
                counts[pg_tup] = 0
            counts[pg_tup] += 1
    
        for tup, cts in counts.items():
            print(f"{", ".join([v for v in tup if v])}: {cts}")
            
        print()
    
    return counts
