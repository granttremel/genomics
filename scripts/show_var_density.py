
import os
from matplotlib import pyplot as plt
import numpy as np
plt.set_loglevel("CRITICAL")

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def count_variants(geneman, chr, start, end):
    
    annos = geneman.get_all_annotations(chr, start, end, include_motifs = False)
    genes = [f for f in annos if f.feature_type == "gene"]
    cds = [f for f in annos if f.feature_type == "CDS"]
    vars = [v for v in annos if v.feature_type == "variant"]
    
    print(f"identified {len(vars)} variants with {len(genes)} genes and {len(cds)} CDS in range")
    
    var_pos = np.array([v.start for v in vars])
    var_lens = np.array([max(1,abs(len(v.attributes.get("alt")[0]) - len(v.attributes.get("ref")))) for v in vars])
    gene_pos = np.array([[g.start, g.end] for g in genes])
    gene_names = [g.name for g in genes if g.name]
    cds_pos = np.array([[cc.start, cc.end] for cc in cds])
    return var_pos, var_lens, gene_pos, gene_names, cds_pos

def merge_feature_lines(feat_pos, min_dist):
    
    out_feat_pos = []
    
    if len(feat_pos) < 2:
        return feat_pos
    
    curr_fp = feat_pos[0]
    for i in range(1, len(feat_pos)):
        fp = feat_pos[i]
        if fp[0] - curr_fp[-1] < min_dist:
            curr_fp[-1] = fp[-1] # extend current
        else:
            out_feat_pos.append(curr_fp) # can't be extended so commit
            curr_fp = fp
    out_feat_pos.append(curr_fp)
    
    print(f"merged {len(feat_pos)} features into {len(out_feat_pos)}")
    
    return out_feat_pos

def plot_feature_lines(ax, feat_pos, feat_y, binrange, **kwargs):
    
    min_dist = (binrange[-1] - binrange[0]) / 100
    min_res = min_dist / 500
    feat_pos = merge_feature_lines(feat_pos, min_dist)
    
    for fp in feat_pos:
        fp_clip = [max(fp[0], binrange[0]), min(fp[1], binrange[-1])]
        if fp_clip[-1] - fp_clip[0] < min_res:
            fp_mid = (fp_clip[-1] + fp_clip[0]) / 2
            ax.plot(fp_mid, feat_y,'.', markersize = 1, **kwargs)
        else:
            ax.plot(fp_clip, [feat_y, feat_y], linewidth = 3,  **kwargs)
    return ax

def add_gene_names(ax, gene_names):
    gene_name_str = '\n'.join(["Genes:"] + gene_names)
    ax.text(0.04, 0.95, gene_name_str, transform = ax.transAxes, fontsize = 9.0, va = "top")
    return ax

def plot_density(var_pos, var_lens, gene_pos = None, gene_names = None, cds_pos = None, idx = "", fname_suffix = ""):
    
    f,ax = plt.subplots()
    vs, bs, _ = ax.hist(var_pos, weights = var_lens, bins = 100)
    brange = [bs[0], bs[-1]]
    
    ylim = ax.get_ylim()
    p = 0.005
    gene_y = ylim[1] * p
    cds_y = ylim[1]*(p + 0.015)
    
    if gene_pos is not None:
        plot_feature_lines(ax, gene_pos, gene_y, brange, c="y")
    
    if cds_pos is not None:
        plot_feature_lines(ax, cds_pos, cds_y, brange, c="firebrick")
    
    if gene_names is not None and len(gene_names) < 15:
        add_gene_names(ax, gene_names)
    
    ax.set_title(idx)
    
    if idx and not fname_suffix:
        fname_suffix = idx
    
    save_fig(f, fname_suffix)
    
    return f, ax

def save_fig(f, fname_suffix = ""):
    fname = f"variant_density_{fname_suffix}.png"
    f.savefig(os.path.join("./out", fname))
    pass

def plot_range(geneman, chr, start, windowsz):
    gene_range = (chr, start, start + windowsz)
    idx = geneman._make_index(*gene_range)
    
    vp, vlens, gp, gns, cdsp = count_variants(geneman, *gene_range)
    
    f,ax = plot_density(vp, vlens, gene_pos=gp, gene_names=gns, cds_pos = cdsp, idx=idx)
    return f, ax

def plot_ranges(geneman, chr, start, windowsz):
    
    res = True
    while res:
        
        print(f"current index: {chr}:{start}-{start+windowsz}")
        f, ax = plot_range(geneman, chr, start, windowsz)
        save_fig(f, fname_suffix = "wkng")
        
        rr = input("keep going?")
        if 'n' in rr.lower():
            res = False
        elif 'r' in rr.lower():
            start -= windowsz
        else:
            start += windowsz
    

def main():
    
    chr = 19
    windowsz = 125000
    start = 1000000
    
    gm = load_genome()
    
    # f, ax =plot_range(gm, chr, start, windowsz)
    plot_ranges(gm, chr, start, windowsz)
    
    

if __name__=="__main__":
    main()

