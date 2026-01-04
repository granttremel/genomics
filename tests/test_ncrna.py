
import numpy as np

from ggene import other_paths
from ggene.database.genome_manager import GenomeManager
from ggene.database.annotations import UFeat

from ggene.genome import genes, ncrna

from ggene.seqs import process, bio, vocab

from ggene import draw
from ggene.scalar_plot import ScalarPlot


def load_genome():
    return GenomeManager(dfam_path = other_paths.get("dfam_path"))

def initialize_rna_library(gm, chrom = ""):
    
    rna_genes_chr = {}
    if chrom:        
        rna_genes_chr[chrom] = ncrna.gather_rna_genes(gm, chrom)
    else:
        rna_genes_chr = ncrna.gather_all_rna_genes(gm)
    
    rna_genes = {}
    for chr, scdict in rna_genes_chr.items():
        for sc, feats in scdict.items():
            if not sc in rna_genes:
                rna_genes[sc] = []
            rna_genes[sc].extend(feats)

    print(f"loaded {sum([len(rnas) for rnas in rna_genes.values()])} rna genes")
    
    for st, rnas in rna_genes.items():
        
        print(f"subclass {st} has {len(rnas)} genes")
        for rna in rnas:
            if not isinstance(rna, UFeat):
                print(f"erroneous RNA: {rna} with subclass {st}")
    
    ncrna.save_rna_genes(rna_genes)
    
    for sc in rna_genes.keys():
        ncrna.save_rna_subclasses(rna_genes, sc)
        
def test_vtrna(gm:GenomeManager, context_sz = 128, pct = 90, show_scplots = False, show_heatmaps = False):
    
    rna_types = [
        "vault_RNA",
        "misc_RNA", 
        # "snoRNA"
    ]
    rna_scs = [
        "VT", 
        "Vault", 
        # "U"
    ]
    
    cts = {rnsc:0 for rnsc in rna_scs}
    feats = []
    
    max_per_sc = 5
    
    for rna_type in rna_types:
        rna_feats = ncrna.load_rna_genes(rna_type, subclasses = rna_scs)
        
        for f in rna_feats:
            
            sc, v = ncrna.assign_rna_subclass(f)
            if cts.get(sc, 0) > max_per_sc:
                continue
            
            if sc == "U" and v[0] ==3:
                continue
            
            full_seq = gm.get_sequence(f.chrom, f.start - context_sz, f.end + context_sz)
            
            f.attributes["full_seq"] = full_seq
            f.attributes["seq"] = full_seq[context_sz:len(full_seq) - context_sz]
            f.attributes["upstream"] = full_seq[:context_sz]
            f.attributes["downstream"] = full_seq[-context_sz:]
            
            feats.append(f)
            cts[sc] += 1
            
        
    print(f"collected {len(feats)} rna features")
    
    cf = lambda a,b:a==b
    # sf = lambda a,b,sc: 1/sc
    
    hmdata = {"fwd_pct":[], "rev_pct":[], "fwd_sd":[],"rev_sd":[]}
    
    for i in range(len(feats)):
        
        fi = feats[i]
        seqi = fi.attributes.get("full_seq","")
        if fi.strand == '-':
            seqi = bio.reverse_complement(seqi)
        
        for v in hmdata.values():
            v.append([])
        
        for j in range(i+1,len(feats)):
            pass
        
            fj = feats[j]
            seqj = fj.attributes.get("seq","")
            if fj.strand == '-':
                seqj = bio.reverse_complement(seqj)
            
            
            max_len = max(len(seqi), len(seqj))
            min_len = min(len(seqi), len(seqj))
            len_ratio = max_len / min_len
            seqii = process.pad_sequence(seqi, max_len)
            seqjj = process.pad_sequence(seqj, max_len)
            
            sf = lambda a,b,sc: len_ratio/sc
            # sf = lambda a,b,sc: (3 if a in "GC" else 2)*len_ratio/sc
            
            res, rcres = process.correlate(seqii, seqjj, comparison_func = cf, score_func = sf, fill = None)
            
            xlabel = draw.highlight_matching(seqii, seqjj, colors = (65+6, 53), do_both = True, suppress = True, color_bg = True)
            
            # print(xlabel[0])
            # print(xlabel[1])
            
            if show_scplots:
                print(f"Plot of gene {fi.name}, {fj.name}")
                sc1 = ScalarPlot(res, add_range = True, minval = 0)
                sc2 = ScalarPlot(rcres, add_range = True, minval = 0)
                ScalarPlot.show_paired(sc1, sc2, chunksz = 256, xlabel = xlabel, center_xlabel = True)
            
            hmdata["fwd_pct"][-1].append(np.percentile(res, [pct])[0])
            hmdata["rev_pct"][-1].append(np.percentile(rcres, [pct])[0])
            hmdata["fwd_sd"][-1].append(np.std(res))
            hmdata["rev_sd"][-1].append(np.std(rcres))
            
            # input()
    
    def map_label(feat, col_width):
        
        if "VT" in feat.name:
            sc_keys = feat.attributes.get("subclass", [])
            return [f"VT",f"{sc_keys[0]}-{sc_keys[1]}", "hii"]
        
        elif "Vault" in feat.name:
            return ["","Vlt", "!!!"]
        else:
            
            return [feat.name[:col_width], feat.name[-col_width:], "RNA"]

    if show_heatmaps:
        row_height = 1
        col_width = 3
        row_lbls = [f.name for f in feats]
        col_lbls = [map_label(f, col_width) for f in feats]
        
        print(col_lbls)
        
        minval = maxval = None
        hmkwargs = dict(row_labels = row_lbls, col_labels=col_lbls, center = None, minval = minval, maxval = maxval, col_width = col_width, row_height = row_height)
        
        for k, data in hmdata.items():
            print(f"Heatmap of {k}")
            color_scheme = "unterra" if "rev" in k else "terra"
            draw.make_heatmap(data, color_scheme = color_scheme, **hmkwargs)
        
            
        
        
    


def main():
    
    
    gm = load_genome()
    
    # test_vc = "ᖰᖱᕫᕬ"
    # test_vc = "ᘔᘖᐁᐃ" # slay!
    # test_vc = "ᖸᖹᘉᘈ" # has potential
    # test_vc = "ᑭᑯᖗᖘ"
    # test_vc = "ᘂᘃᖼᖽ"
    
    # initialize_rna_library(gm, chrom = "")
    
    test_vtrna(gm, context_sz = 64, show_scplots = True)


    # sick_chars = draw._sick_chars
    # step = 32
    # for i in range(len(sick_chars) // step):
    #     print(" ".join(sick_chars[i*step:(i+1)*step]))





    pass


if __name__=="__main__":
    main()










