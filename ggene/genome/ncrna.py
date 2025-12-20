
import re
import json
from difflib import SequenceMatcher

from ggene import DEFAULT_LIBRARY
from ggene.database.genomemanager import GenomeManager
from ggene.database.unified_stream import UFeature

LIB_DIR = "ncache/ncrnas"

rna_types = ['rRNA','rRNA_pseudogene','lncRNA','snRNA','snoRNA','miRNA','scaRNA','scRNA','sRNA','vault_RNA','misc_RNA']

rna_subclass={
    
    'rRNA':{
        "RNA5SP{n}":re.compile(r'RNA5SP([0-9])([0-9]+)'),
        "RNA5S{n}":re.compile(r'RNA5S([0-9]+)'), 
        "RNA5-8S[NP]{n}":re.compile(r'RNA5-8S([NP])([0-9]+)'),
        '5_8S_rRNA':'5_8S_rRNA', 
        '5S_rRNA':'5S_rRNA',
    },
    
    'lncRNA':{
        "TTTY":re.compile(r'TTTY([0-9][AB]?)'),
        "gene + LP{n}":re.compile(r'([A-Z0-9]+)([PL][0-9])'),
        "idot":re.compile(r'([A-Z0-9orf\-]+)-[IDO]T([0-9]?)'),
        
        "long intergenic noncoding":re.compile(r'LINC([0-9]+)'),
        "long intergenic noncoding 2":re.compile(r'([A-Z0-9]+)LNC'),
        "long intergenic noncoding 3":re.compile(r'LNC([A-Z0-9]+)'),
        "antisense":re.compile(r'([A-Z0-9orf\-]+)-?AS([0-9]?)'),
        "micro RNA":re.compile(r'MIR([A-Z0-9]+)'),
        "U RNA":re.compile(r'RNU([A-Z0-9\-]+)'),
        "orf":re.compile(r'([A-Z0-9]+orf[A-Z0-9]+)'),
    },
    
    'snRNA':{
        'RNU':re.compile(r'RNU([12467]|5[A-F])-?([0-9]+)(P)?'),
        'RNVU':re.compile(r'RNVU1-([0-9]{1,2})'),
        'RNU + ATAC':re.compile(r'RNU([46])ATAC([0-9]{1,2}P)?'),
        'RN7SK':'RN7SK',
        'RNU6V':"RNU6V",
        'U':re.compile(r'U([12467])'),
    },
    
    'snoRNA':{
        'SNO':re.compile(r'SNOR([AD])([0-9]{1,3})([A-Z]?)'),
        'SCA':re.compile(r'SCARNA(10|12|18)(B?)'),
        'U':re.compile(r'U([38])'),
        'RNU':re.compile(r'RNU105([BC])'),
        'snoZ196':'snoZ196',
    },
    
    'miRNA':{
        'MIR':re.compile(r'MIR([0-9]+)-?([1-9])?'),
        'MIRLET7':re.compile(r'MIRLET7([A-I])([1-3])?'),
        'SNORD138':'SNORD138',
        'hsa-mir-8069-1':'hsa-mir-8069-1',
        
    },
    
    'scaRNA':{
        'SCA':re.compile(r'SCARNA([0-9]+)(B)?'),
    },
    
    'scRNA':{
        'BCYRN1':'BCYRN1'
    },
    
    'sRNA':{
    },
    
    'vault_RNA':{
        'VT':re.compile("VTRNA([0-9])-([0-9])(P)?")
    },
    
    'misc_RNA':{
        'RN7SL(P)':re.compile(rf'RN7SL([0-9]+)(P)?'),
        'RNY':re.compile(r'RNY([134])(P[1-9])?'),
        'VT':re.compile("VTRNA([0-9])-([0-9])(P)?"),
        'Y_RNA':'Y_RNA',
        '7SK':'7SK',
        'Vault':'Vault',
        'Metazoa_SRP':'Metazoa_SRP',
        'Telomerase-vert':'Telomerase-vert',
        'NRON':'NRON',
        
    }
}
rna_subclass['rRNA_pseudogene'] = rna_subclass['rRNA']


def gather_all_rna_genes(gm:GenomeManager):
    rna_genes = {}
    for chr in gm.iter_chromes():
        rna_genes[chr] = gather_rna_genes(gm, chr)
    
    return rna_genes


def gather_rna_genes(gm:GenomeManager, chrom):
    
    rna_genes = {bt:[] for bt in rna_types}
    
    print(f"######### starting chr{chrom} ####### \n")
    
    ng = 0
    nrg = 0
    for gene in gm.annotations.stream_by_types(['gene'], chrom):
        ng += 1
        gn = gene.get("gene_name", gene.get("name",""))
        bt = gene.get("attributes",{}).get("gene_biotype")
        
        if not bt in rna_types:
            continue
        
        # if not "RNA" in bt:
        #     continue
        
        sc, v = assign_rna_subclass(gene)
        
        if not bt in rna_genes:
            # rna_types.append(bt)
            rna_genes[bt] = []
        rna_genes[bt].append(gene)
        nrg += 1
    
    print(f"searchged {ng} genes on chr{chrom} and found {nrg} rna genes")
    return rna_genes


def save_rna_genes(rna_genes, ftag = ""):
    
    for rnt, genes in rna_genes.items():
        json_str = json.dumps({rnt:[g.to_dict() for g in genes]})
        
        for repstr in ["[",r"}},"]:
            json_str = json_str.replace(repstr, repstr + "\n")
            pass
        
        
        fname = DEFAULT_LIBRARY / LIB_DIR / (f"rna_genes_{rnt}" + ("_" + ftag if ftag else "") +".json")
        # fname.mkdir(exist_ok = True)
        # print(f"saving rna genes to {fname}")
        with open(fname,"w+") as f:
            f.write(json_str)

def load_rna_genes(rna_type, ftag = "", subclasses = None):
    
    fname = DEFAULT_LIBRARY  / LIB_DIR /  (f"rna_genes_{rna_type}" + ("_" + ftag if ftag else "") + ".json")
    
    with open(fname,'r') as f:
        rna_genes = json.loads(f.read())
    
    rna_genes = rna_genes.get(rna_type,[])
    rna_feats = []
    for gd in rna_genes:
        gd["feature_type"] = gd.pop("type", "")
        feat = UFeature(**gd)
        
        sc, v = assign_rna_subclass(feat)
        
        if subclasses and not sc in subclasses:
            continue
        
        rna_feats.append(feat)
    
    return rna_feats


def save_rna_subclasses(rna_subclasses, subclass_type,  ftag = ""):
    rna_cls_names = {clsn:[f.name if f.name else f.id for f in feats] for clsn, feats in rna_subclasses.items()}
    
    fname = DEFAULT_LIBRARY  / LIB_DIR /  (f"{subclass_type}_subclasses" + ("_" + ftag if ftag else "") + ".json")
    
    with open(fname, "w+") as f:
        json.dump(rna_cls_names, f, indent = 1)

def load_rna_subclasses(subclass_type, ftag = ""):
    
    fname = DEFAULT_LIBRARY  / LIB_DIR / (f"{subclass_type}_subclasses" + ("_" + ftag if ftag else "") + ".json")
    
    with open(fname, "r") as f:
        subclass_names = json.load(f)
    
    return subclass_names

def assign_rna_subclass(rna_feat, gene_names = [], do_search = False):
    
    name = rna_feat.name
    if not name:
        return "None", tuple()
    
    sc_ptrns = rna_subclass.get(rna_feat.attributes.get("gene_biotype",""),[])
    if not sc_ptrns:
        return None, None
    
    
    subclass = None
    val = None
        
    for scn, sc in sc_ptrns.items():
        if isinstance(sc, re.Pattern):
            
            m = re.match(sc, rna_feat.name)
            if m:
                subclass = scn
                val = m.groups()
                rna_feat.attributes["subclass"] = tuple(v if v else "" for v in val)
                break
        
        elif isinstance(sc, list):
            if rna_feat.name in sc:
                subclass = scn
                val = tuple()
                break
        
        elif sc in rna_feat.name:
                subclass = scn
                val = tuple()
                rna_feat.attributes["subclass"] = val
                break
    
    if not subclass:
        subclass, val = is_rna_subclass_gene(rna_feat, gene_names = gene_names, do_search = do_search)
    
    return subclass, val

def is_rna_subclass_gene(rna_feat, gene_names = [], do_search = False):
    
    if not rna_feat or not rna_feat.name:
        return None, None
    import string
    
    rna_name = rna_feat.name
    
    if rna_name in gene_names:
        return "gene name", rna_name
    
    max_ratio = 0
    max_gn = ""
    for gn in gene_names:
        if rna_name.startswith(gn):
            return "prefixed gene name", gn
        elif rna_name.startswith(gn.rstrip(string.digits)):
            print(rna_name, gn)
            return "prefixed gene + n", gn
    
    if do_search:
        for gn in gene_names:
            s = SequenceMatcher(None, rna_name, gn)
            r = s.ratio()
            if r > 0.6:
                return "near gene name", gn
            
            if r > max_ratio:
                max_ratio = r
                max_gn = gn
    
        rna_feat.attributes["nearest_gene"] = max_gn
        print(f"closest match to rna name {rna_name} is {max_gn} with ratio {max_ratio:0.3f}")
    
    return None, None

def classify_rnas(rnas, rna_classes = {}, maxct = 100, gene_names = [], do_search = False):
    
    unknowns = []
    
    ct = 0
    for rna in rnas:
        
        f, v = assign_rna_subclass(rna, gene_names = gene_names, do_search=do_search)
        
        if f is None:
            unknowns.append(rna)
            ct += 1
        else:
            if not f in rna_classes:
                rna_classes[f] = []
            rna_classes[f].append(rna)
        
        if ct >= maxct:
            break
    
    return rna_classes, unknowns