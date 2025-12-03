

import random
from typing import List, Dict, Union, Any, Optional
from dataclasses import dataclass

import numpy as np

from tabulate import tabulate


from ggene.genomemanager import GenomeManager
from ggene import draw
from ggene.scalar_plot import ScalarPlot
from ggene.seqs import bio, process, find, align, heal
from ggene.seqs.bio import ALIASES, ALIASES_REV

def score_sequences_corrs(seqa, seqb, topk = 5, scale = None, **kwargs):
    min_len = min(len(seqa), len(seqb))
    corrs, rccorrs = process.correlate(seqa, seqb, scale = scale, fill = 0.25)
    topk_corrs = sorted(corrs)[:topk]
    topk_rccorrs = sorted(rccorrs)[:topk]
    return sum(topk_corrs)/topk, sum(topk_rccorrs)/topk

def score_sequences_runs(seqa, seqb, topk = 5, max_err = 16, scale = None, **kwargs):
    runs, _, shifts, _ = process.correlate_longest_subseq_err(seqa, seqb, max_err, scale=scale)
    topk_runs = sorted([(r,s) for r,s in zip(runs, shifts) if s > 3], key = lambda k:-k[0])[:topk]    
    rcruns, _, rcshifts, _ = process.correlate_longest_subseq_err(seqa, bio.reverse_complement(seqb), max_err, scale=scale)
    topk_rcruns = sorted([(r,s) for r,s in zip(rcruns, rcshifts) if s > 3], key = lambda k:-k[0])[:topk]    
    return sum([r for r, s in topk_runs])/topk, sum([r for r, s in topk_rcruns])/topk

def score_sequences_algn(seqa, seqb, **kwargs):
    min_len = min(len(seqa), len(seqb))
    score = align.score_sequences(seqa, seqb)
    rcscore = align.score_sequences(seqa, bio.reverse_complement(seqb))
    return score/min_len, rcscore/min_len

def downsample_scores(scores):
    
    scores = np.array(scores)
    nr, nc = scores.shape
    
    ds_score = 0.5*scores[1:nr-1, 1:nc-1]
    ds_score += 0.25*scores[2:nr, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 2:nc]
    ds_score += 0.25*scores[0:nr-2, 1:nc-1]
    ds_score += 0.25*scores[1:nr-1, 0:nc-2]
    
    return ds_score

def get_score_function(score_mode):
    if score_mode == "alignment":
        return score_sequences_algn
    elif score_mode == "runs":
        return score_sequences_runs
    elif score_mode == "corrs":
        return score_sequences_corrs
    else:
        return None

def compare_sequences(seqa, seqb, chunksz = 128, score_modes = ["alignment", "runs", "corrs"], resample = True, **kwargs):
    
    if len(seqb) > len(seqa):
        seqa, seqb = seqb, seqa
    
    q = 1
    if resample:
        q = 2
        chunksz = chunksz//2
    
    nchksa = len(seqa) // chunksz
    nchksb = len(seqb) // chunksz
    
    score_funcs = [get_score_function(sm) for sm in score_modes]
    
    scores = [[[] for n in range(nchksb)] for sf in score_funcs]
    rcscores = [[[] for n in range(nchksb)] for sf in score_funcs]
    row_lbls = [f"b{n}" for n in range(nchksb)]
    
    ruler = kwargs.pop("add_ruler", False)
    xmin = kwargs.pop("xmin", None)
    xmax = kwargs.pop("xmax", None)
    num_labels = kwargs.pop("num_labels", 5)
    
    # print(f"comparing seqa length {len(seqa)}, {nchksa} chunks with seqb length {len(seqb)}, {nchksb} chunks")
    
    for ia in range(nchksa):
        ssa = seqa[chunksz*ia:chunksz*(ia+q)]
        
        for ib in range(nchksb):
            
            ssb = seqb[chunksz*ib:chunksz*(ib+q)]
            if not ssa or not ssb:
                continue
                
            for i, sf in enumerate(score_funcs):
                score, rcscore = sf(ssa, ssb, **kwargs)
                scores[i][ib].append(score)
                rcscores[i][ib].append(rcscore)
    
    all_hms = [[] for r in range(len(scores[0])+3)]
    all_rchms = [[] for r in range(len(rcscores[0])+3)]
    
    for nsm,sm in enumerate(score_modes):
        
        if nsm > 0:
            row_lbls = []
        
        if resample:
            scores_rs = downsample_scores(scores[nsm])
            rcscores_rs = downsample_scores(rcscores[nsm])
        else:
            scores_rs = scores[nsm]
            rcscores_rs = rcscores[nsm]
        
        hm = draw.heatmap(scores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra", col_space = 0, row_space = 0,
                          ruler = ruler, xmin = xmin, xmax = xmax, num_labels=num_labels)
        
        for nr in range(len(hm)):
            all_hms[nr].append(hm[nr])
        
        rchm = draw.heatmap(rcscores_rs, center = None, row_labels = row_lbls, suppress = True, color_scheme = "terra", col_space = 0, row_space = 0,
                          ruler = ruler, xmin = xmin, xmax = xmax, num_labels=num_labels)
        
        for nr in range(len(rchm)):
            all_rchms[nr].append(rchm[nr])
    
    frm = "{:^32}       " * len(score_modes)
    for row in all_hms:
        if row:
            print(frm.format(*row))
    for row in all_rchms:
        if row:
            print(frm.format(*row))
    
    return scores


def walk_sequences(seqa, seqb, chunksz, num_steps, **kwargs):
    
    asublen = len(seqa)//num_steps
    bsublen = len(seqb)//num_steps
    
    print(f"a sublength: {asublen}, b sublength: {bsublen}")
    
    for n in range(num_steps):
        
        subseqa = seqa[n*asublen:(n+1)*asublen]
        subseqb = seqb[n*bsublen:(n+1)*bsublen]
        
        compare_sequences(subseqa, subseqb, chunksz=chunksz, **kwargs)
    