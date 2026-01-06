
from typing import Dict, List
import time
from pathlib import Path
import gzip
from tabulate import tabulate
import random

import numpy as np

from ggene.display.artists.line_artist import LineArtist, LineArtistParams
from ggene.database.annotations import chr_lens
from ggene.config import get_config, get_paths, DATA_DIR
from ggene.draw.scalar_plot import ScalarPlot
from ggene.draw.heatmap import make_heatmap, Heatmap

from ggene.seqs import process
from ggene.motifs.hmm import load_hmms, load_hmm_classes, AB

from ggene.draw.colors import Colors

from ggene.seqs import process, bio, gen

from ggene.database.genome_manager import GenomeManager
# from scripts.make_chrom_maps import get_random_seq


def load_genome():
    gm = GenomeManager()
    return gm

def random_crop_seq(seq, to_length):
    
    seq_len = len(seq)
    
    if seq_len <= to_length:
        return seq
    
    n = random.randint(0, seq_len - to_length)
    
    return seq[n:n+to_length]

def print_feat_summary(feats):
    
    # print(feats)
    
    fd = {}
    
    for f in feats:
        if not f.feature_type in fd:
            fd[f.feature_type] = []
        fd[f.feature_type].append(f)
    
    # print(fd)
    
    for ft, fs in sorted(fd.items(), key = lambda k:-max(kk.length for kk in k[1])):
        print(f"feature type {ft} with {len(fs)} members.")
        fm = sorted(fs, key = lambda f:-f.length)[0]
        print(f"  {fm}")

def extra():
    
    # smat, srcmat = process.get_combined_inter_matrix(kernel, kernel)
    # print(smat.shape)
    # shrsmat = process.shear_matrix(smat, fill = np.nan)
    # shrsrcmat = process.shear_matrix(srcmat, fill = np.nan)
    # Heatmap(shrsmat, half_block = True, colorbar = True, minval=-3, maxval=3, colors=cols).show()
    # Heatmap(shrsrcmat, half_block = True, colorbar = True, minval=-3, maxval=3, colors=cols).show()
    # input()
    # Heatmap(smat==0, half_block = True, colorbar = True, colors=cols2).show()
    # Heatmap(srcmat==0, half_block = True, colorbar = True, colors=cols2).show()
    
    # mini_smat = smat[:32,:32]
    # smat_shear = process.shear_matrix(mini_smat, axis = 0)
    # smat_shear1 = process.shear_matrix(mini_smat, axis = 1)
    
    # Heatmap(mini_smat, half_block = True, colorbar = True, minval=-3, maxval=3, colors=cols).show()
    # Heatmap(smat_shear, half_block = True, colorbar = True, minval=-3, maxval=3, colors=cols).show()
    # Heatmap(smat_shear1, half_block = True, colorbar = True, minval=-3, maxval=3, colors=cols).show()
    
    # input()
    pass


def scan_kernel(gm, kernel, chrom, start, end=None, do_random = False, cmap = "terra", seq_len = 64, hmm_str = ""):
    
    # la = LineArtist("line_artist", LineArtistParams(display_width = 256, display_height = 8, use_global_features = False))
    
    if do_random:
        kernel = gen.get_random_sequence(len(kernel))
    
    # print(f"kernel: {kernel}")
    
    seq = gm.get_sequence(chrom, start, end=end)
    
    # c1, c2 = Colors.get_color_scheme_24b(cmap)
    # cols = Colors.add_middle(c1, c2, brightness = 0.8, saturation = 0.8)
    
    cols2 = [
        [20,20,20],
        [200,200,200],
    ]
    
    cols3 = [
        [120,20,20],
        [20,20,20],
        [20,120,120],
    ]
    # pretty!
    cols5 = [
        [20,20,20],
        [120,20,20],
        [120,120,20],
        [120,20,120],
        [20,120,120],
    ]
    
    cols7 = [
        [20,120,120],
        [120,120,20],
        [120,20,120],
        # [20,20,20],
        [150,150,150],
        [20,120,20],
        [20,20,170],
        [120,20,20],
    ]
    
    all_bs = "KWRNYSM"
    
    focus_ind = -1 # 0=all, 1=RY, 2=SW, 3=MK
    bs = all_bs
    cols = cols7
    
    color_key = ", ".join(["{}{}".format(Colors.get_color(c), f"██ alias = {b} ██") for c,b in zip(cols, bs)]) + Colors.RESET
    
    num_chunks = int(len(seq)/seq_len) + 1
    chunksz = int(len(seq) / num_chunks)
    col_labels = list(kernel)
    
    if focus_ind > -1:
        ffmat= process.get_full_inter_matrix(kernel, kernel)
        fmat = ffmat[focus_ind, 0]
    else:
        ffmat = process.get_combined_inter_matrix(kernel,kernel)
        fmat = ffmat[0]
    
    hm0 = Heatmap(fmat, half_block = True, colorbar = True, colors = cols, col_labels = col_labels, row_labels=col_labels)
    hm0.show()
    # rows = hm0.get_rows()
    shrmat = process.shear_matrix(fmat, fill=np.nan, axis = 0)
    hm = Heatmap(shrmat, half_block = True, colorbar = True, colors = cols, col_labels = col_labels, row_labels=col_labels)
    hm.show()
    # rows = hm.get_rows()
    
    rows = hm0.get_rows()+hm.get_rows()
    
    print(color_key)
    
    with open("./notes/fotos/alu_auto", 'a+') as f:
        f.write(hmm_str)
        for row in rows:
            f.write(row + "\n")
    
    input()
    
    for n in range(num_chunks):
        
        sst = n*chunksz
        een = (n+1)*chunksz
        sseq = seq[sst:een]
        
        if do_random:
            print("doing random seq")
            sseq = gen.get_random_sequence(len(sseq))
        
        row_labels = list(sseq)
        
        feats = [f for f in gm.annotations.stream_all(chrom, start+sst, start+een)]
        
        # ffmat= process.get_combined_inter_matrix(sseq,kernel)
        
        if focus_ind > -1:
            ffmat= process.get_full_inter_matrix(sseq, kernel)
            fmat = ffmat[focus_ind, 0]
        else:
            ffmat = process.get_combined_inter_matrix(sseq,kernel)
            fmat = ffmat[0]
        
        hm0 = Heatmap(fmat, half_block = True, colorbar = True, colors = cols, col_labels = col_labels, row_labels=row_labels)
        rows = hm0.get_rows()
        
        for r in rows:
            time.sleep(0.01)
            print(r)
        print()
        
        shrmat = process.shear_matrix(fmat, fill=np.nan, axis = 0)
        
        hm = Heatmap(shrmat, half_block = True, colorbar = True, colors = cols, col_labels = col_labels, row_labels=row_labels)
        rows = hm.get_rows()
        # rows = hm.show(suppress = True)
        
        for r in rows:
            time.sleep(0.01)
            print(r)
        print()
        
        print(color_key)
        
        # frcmat = ffmat[focus_ind, 1]
        # rcshrmat = process.shear_matrix(frcmat, fill=np.nan, axis = 1)
        # rchm = Heatmap(rcshrmat, half_block = True, colorbar = True, colors = cols)
        # rcrows = rchm.show(suppress = False)

        # print(kernel, 'kernel')
        # print(sseq ,'seq')
        
        print_feat_summary(feats)
        
        res = input("keep?")

        if 'y' in res.lower():
            with open("./notes/fotos/blobby", 'a+') as f:
                f.write(f"chr{chrom}:{start+sst}:{start+een}")
                for row in hm.get_rows():
                    f.write(row + "\n")

        n+=1

def test_hmms(num_rounds = 10, max_len = 512, do_random = False):
    
    all_hmms = load_hmms()
    hmms = [hmm for hmm in all_hmms if "L1P" in hmm.name.decode("utf-8")]
    hmm2s = [hmm for hmm in all_hmms if "AluY" in hmm.name.decode("utf-8")]
    
    done = set()
    
    n = 0
    while n < num_rounds:
        
        hmm1 = random.choice(hmms)
        hmm2 = random.choice(hmm2s)
        
        if (hmm1.name, hmm2.name) in done:
            continue
        
    
        if do_random:
            print("doing random")
            seqa = gen.get_random_sequence(max_len)
            seqb = gen.get_random_sequence(max_len)
            
        else:
            print(f"chose hmms {hmm1.name} and {hmm2.name} with M's {hmm1.M}, {hmm2.M}")
    
            seqa = hmm1.consensus
            seqb = hmm2.consensus
            
            if len(seqa) > max_len:
                seqa = random_crop_seq(seqa, max_len)
            if len(seqb) > max_len:
                seqb = random_crop_seq(seqb, max_len)
        
        done.add((hmm1.name, hmm2.name))
        
        mat = process.get_inter_matrix(seqa, seqb)
        mat_aproj = np.mean(mat, axis = 0)
        mat_bproj = np.mean(mat, axis = 1)
        
        ran = max(np.max(np.abs(mat_aproj)), np.max(np.abs(mat_bproj)))
        center = 0
        
        hm = Heatmap(mat, half_block = True, colorbar = True, symmetric_color = True, color_scheme = "lava")
        rows = hm.show(suppress = True)
        
        btrim = Heatmap(mat_bproj[:,None], half_block = False, center=center, minval=-ran, colorbar = False, color_scheme = "lava").show(suppress = True)
        
        for row, end in zip(rows, btrim):
            print(row + "  " + end.strip().replace(" ", "▀▀"))
        
        atrim = Heatmap(mat_aproj[None,:], half_block = False, center=center, minval = -ran, colorbar = True, color_scheme = "lava").show()
        
        print(f"chose hmms {hmm1.name} and {hmm2.name} with M's {hmm1.M}, {hmm2.M}, seq lens {len(seqa)}, {len(seqb)}")
        print(hmm1.name, hmm1.description)
        print(hmm2.name, hmm2.description)
        
        res = input("keep?")
        
        if 'y' in res.lower():
            
            with open("./notes/fotos/idk", 'a+') as f:
                
                f.write(f"HMM1 {hmm1.name} {hmm1.description}\n")
                f.write(f"HMM2 {hmm2.name} {hmm2.description}\n")
                
                for row in hm.get_rows():
                    
                    f.write(row + "\n")
        
        n+=1

def get_random_gene(gm):
    
    chrom = random.choice([chrom for chrom in gm.iter_chromes()])
    
    chr_len = chr_lens.get(chrom)
    
    pos = int(random.random()*chr_len)
    
    gene = None
    for f in gm.annotations.stream_by_types(['gene'], chrom, pos, end=None):
        
        if f and f.name and f.feature_type=='gene':
            gene = f
    
    return gene

def main():

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--name","-n", type = str, default = "AluJr",)

    args = parser.parse_args()

    kern_name = args.name

    cmap = 'aqua_vitae'
    
    cfg = get_config()
    vars, genes, seqs, lib, other = get_paths()
    gm = load_genome()

    # test_hmms(num_rounds = 5, do_random = True)
    # test_hmms()
    
    hmms = load_hmms()
    
    # kern_name = "MER4"
    # kern_name = "MIR3"
    # kern_name = "ALRa"
    
    kernel = None
    for hmm in hmms:
        if kern_name in hmm.name.decode():
            # print(hmm.name, hmm.description)
            cons = hmm.consensus
            if len(cons) > 350:
                continue
            
            kernel = cons
            break
    
    # kernel = random_crop_seq(kernel, 256)
    
    # 8:145052465-145066685
    # chrom = "8"
    # start = 145052465
    # end = start + 16384
    chrom, start = gm.get_random_location()
    length = 16384
    # gene = get_random_gene(gm)
    
    # chrom=gene.chrom
    # start=  gene.start
    # end = gene.end
    
    # print(f"randomly selected gene {gene}")
    # input()
    
    scan_kernel(gm, kernel, chrom, start, start+length, do_random = False, cmap = cmap, seq_len = 2*len(kernel), hmm_str = f"{hmm.name}: {hmm.description}")
    



if __name__=="__main__":
    main()



