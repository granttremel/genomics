
import re
import json

from ggene.config import DEFAULT_LIBRARY
from ggene.database.genome_manager import GenomeManager
from ggene.database.annotations import UFeature

LIB_DIR = "ncache/TEs"

te_types = [


]
"""

























"""

te_subclass = {
    
    "L1":re.compile(r"L1(P[AB]?|HS|M[123A-E]?||)([A-Ga-z0-9]{0,3}|REC2|15-16)(_[35]end|_orf2)"), # perf
    
    "L2-4":re.compile(r"(L[234])(.*)"), # 
    
    "AluJ":re.compile(r"AluJ([bor]4?)"), # perf
    "AluS":re.compile(r"AluS([a-z][0-9]?)"), # perf
    "AluY":re.compile(r"AluY([a-z][0-9]?)?"), # perf
    
    "SVA":re.compile(r"SVA_?([A-Z])"), # perf
    "M[EI]R":re.compile(r"M[EI]R([0-9]*)([A-Fa-d][0-9]*)?(-int|_Amn)?"), # perf
    # "MER":re.compile(r"MER([0-9]*)([A-Za-z]*)?"),
    # "MIR":re.compile(r"MIR([abc13])?(_Amn)?"),
    "LTR":re.compile(r"LTR([0-9]*[A-Z]?)(-int)?"), # perf
    "MLT":re.compile(r"MLT([12])([A-Z])([0-5])?(-int)?"), # perf
    "MST":re.compile(r"MST([A-D])([1-2])?"), # perf
    "MA[RD]E":re.compile(r"(MA[RD]E)([0-9])?"),
    "UCON":re.compile(r"UCON_?(1?[0-9]{1,2})"), # perf
    "HUERS":re.compile(r"HUERS-P([123]b?)"), # perf
    "U":re.compile(r"U([0-9]3?)"), # perf
    
    "THE1":re.compile(r"THE1([A-D])?(-int)?"), # perf
    "HY":re.compile(r"HY([1-5])"),
    "Eutr":re.compile(r"Eutr([12]?[0-9B])"),
    "EUTREP":re.compile(r"EUTREP(1?[0-9])"),
    "ORSL":re.compile(r"ORSL(-2[ab])?"),
    "SATR":re.compile(r"SATR(1|2)"),
    "PABL":re.compile(r"PABL_(A|B)(-int)?"),
    "HERV":re.compile(r"H?ERV([A-Z0-9]*)(-Fc1|-int)?(_LTR1)?"),
    "Eut":re.compile(r"Eut(hAT|Tc1)?(-N?[012]a?)"),
    "Amn":re.compile(r"Amn(L2-1|Harb1|SINE[12])"),
    "Alr":re.compile(r"Alr([ab]?)"),
    
    "CR1":re.compile(r"CR1(-[A-Za-z0-9]+)?_?(AMi|Croc|Amn(i-1)?)"),
    "PRIMA":re.compile(r"PRIMA(4l?|X)(-int|_LTR)?"), # perf
    
    "Mam":re.compile(r"Mam(Rep|RTE|Gypsy|Gyp(LTR)?|SINE1|Tip|_R4)([0-9]*)?(-int)?"),
    "_Mam":re.compile(r"(Helitron|Chap1a?|DNA1|CR1|)([123]N[ab])?_Mam"),
    
    "FXAM":re.compile(r"F[LR]?AM(_[AC])?"),
    "[GH]SAT":re.compile(r"[GH]SAT(5|II)?"),
    
    "X-LINE/DNA":re.compile(r"X([0-9]+[A-Za-z]?)_(LINE|DNA)"),
    
    'named':re.compile(r'(Zaphod|Charlie|Tigger|Ricksha|Eulor|Kanga|Arthur)_?([A-Za-z0-9]+)?'), # DNA transposons
    'hAT':re.compile(r'hAT-([0-9A-Za-z]*)(?:-|_)(Mam|Ther)|OldhAT1'),
    'tRNA':re.compile(r'tRNA(?:-|_)(?:i|([A-Za-z]{3})(-[ATCGY]{3})?(_v|-i)?)'),
    
    'Looper':'Looper',
    'Cheshire':'Cheshire',
    'Harlequin':'Harlequin',
    'BLACKJACK':'BLACKJACK',
    'other':re.compile(r'CER|L5|FordPrefect(_a)?|5S|7SK|7SLRNA|Alr[ab]?|LOR1(a|b|-int)||BC200|BSR[ad]|MSR1|Plat_L3|HSMAR[12]|HAL1(M[8E]|b)|LSU-rRNA_(Cel|Hsa)|HAL1|LFSINE_Vert|MARNA|MST-int|REP522|SS(T1|U-rRNA_Hsa|PrimLTR79)')
    # 'other':'CER|L5|FordPrefect|FordPrefect_a|5S|7SK|7SLRNA|Alr|Alra|Alrb|LOR1(a|b|-int)|BC200|BSR[ad]|MLT1?(-int)?|MSR1|Plat_L3|HSMAR[12]|HAL1(M[8E]|b)|LSU-rRNA_(Cel|Hsa)|HAL1|LFSINE_Vert|MARNA|MST-int|REP522|SS(T1|U-rRNA_Hsa|PrimLTR79'.split("|")
    
}


def gather_all_tes(gm:GenomeManager):
    rna_genes = {}
    for chr in gm.iter_chromes():
        rna_genes[chr] = gather_tes(gm, chr)
    
    return rna_genes



def gather_tes(gm:GenomeManager, chrom):
    
    te_types = []
    tes ={}
        
    print(f"######### starting chr{chrom} ####### \n")
    
    for gene in gm.annotations.stream_by_types(['gene'], chrom, start = 2e6):
        
        gn = gene.get("gene_name", gene.get("name",""))
        bt = gene.get("attributes",{}).get("gene_biotype")
        if not "pseudo" in bt:
            continue
        
        if not bt in te_types:
            te_types.append(bt)
            tes[bt] = []
        tes[bt].append(gene)
    
    return tes

def write_tes(tes, ftag = ""):
    
    for tet, genes in tes.items():
        json_str = json.dumps({tet:[g.to_dict() for g in genes]})
        
        for repstr in ["[",r"}},"]:
            json_str = json_str.replace(repstr, repstr + "\n")
            pass
        
        fname = DEFAULT_LIBRARY / LIB_DIR / (f"tes_{tet}" + ("_" + ftag if ftag else "") + ".json")
        with open(fname.absolute(),"w+") as f:
            f.write(json_str)

def load_tes(te_type, ftag = ""):
    
    fname = DEFAULT_LIBRARY  / LIB_DIR /  (f"tes_{te_type}" + ("_" + ftag if ftag else "") + ".json")
    
    with open(fname.absolute(),'r') as f:
        tes = json.loads(f.read())
    
    tes = tes.get(te_type,[])
    te_feats = []
    for gd in tes:
        gd["feature_type"] = gd.pop("type", "")
        te_feats.append(UFeature(**gd))
    
    return te_feats

def get_te_subclass(te_biotype):

    te_tup = []
    
    for te_group in te_types:
        for acro, opt in te_group.items():
            if opt in te_biotype:
                te_tup.append(acro)
                break
    
    return "_".join(te_tup)

def examine_subclass(subclasses, subclass, show_unique_insts = False, show_insts = False):
        
    sc_types = [set() for i in range(9)]
    
    insts = []
    
    for feat in subclasses.get(subclass):
        scs = feat.attributes.get("subclass_data",[])
        if show_insts:
            print(feat.name, subclass, feat)
        
        inst = []
        for scval, scset in zip(scs, sc_types):
            inst.append(scval)
            scset.add(scval)
        insts.append((feat.name, tuple(inst)))
        
    # print()

    if show_insts:
        insts = sorted(insts)
        for inst in insts:
            print(inst)
    elif show_unique_insts:
        insts = sorted(list(set(insts)), key = lambda k:"".join(k[1]))
        for inst in insts:
            print(f"{inst[0]}: {inst[1]}")
    # print()
    
    sc_types = [sc for sc in sc_types if sc]
    
    print("subclass values:")
    for i, sct in enumerate(sc_types):
        if sct:
            scts = [str(t) for t in sct if t]
            print(i, "\t", ", ".join(sorted(scts)))
        else:
            print(i, "\t", "-")
    
    return sc_types

def assign_subclass(rpt):
    
    name = rpt.name
    if not name:
        return "None", tuple()
    
    sc_ptrns = te_subclass
    
    subclass = "None"
    val = tuple()
        
    for scn, sc in sc_ptrns.items():
        if isinstance(sc, re.Pattern):
            
            m = re.match(sc, rpt.name)
            if m:
                subclass = scn
                val = m.groups()
                rpt.attributes["subclass"] = tuple(v if v else "" for v in val)
                break
        
        elif isinstance(sc, list):
            if rpt.name in sc:
                subclass = scn
                val = tuple()
                break
        
        elif sc in rpt.name:
                subclass = scn
                val = tuple()
                rpt.attributes["subclass"] = val
                break
    
    
    val = [str(v) if val is not None else "" for v in val]
    
    return subclass, val

def characterize_tes(gm:GenomeManager):
    
    counts = {}
    
    for chr in gm.iter_chromes():
        for f in gm.annotations.stream_all(chr, start = 0):
            bt = f.attributes.get("gene_biotype")
            if not bt or not "pseudo" in bt:
                continue
            
            te_tup = get_te_subclass()
            
            if not te_tup in counts:
                counts[te_tup] = 0
            counts[te_tup] += 1
    
        for tup, cts in counts.items():
            print(f"{", ".join([v for v in tup if v])}: {cts}")
            
        print()
    
    return counts
