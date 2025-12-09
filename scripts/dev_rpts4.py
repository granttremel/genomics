

from dataclasses import dataclass
from typing import Dict, List, Set, Optional, Any

import numpy as np
import random

import networkx as nx

from Bio import Align

from ggene import DATA_DIR, other_paths
from ggene.genomemanager import GenomeManager
from ggene import unified_stream, sequence_stream
from ggene import draw
from ggene.scalar_plot import ScalarPlot

from ggene.seqs import bio, process, find, align, repeats, lambdas
from ggene.seqs.vocab import VOCAB


def load_genome():
    gm = GenomeManager(**other_paths)
    rpts = gm.annotations.streams.get("repeats")
    return gm, rpts


def id_repeat(gm, rpts, motif, chrs = None, max_ham = 0):
    
    mtf_min, _, _, _ = bio.get_least_seq_index(motif)
    
    print(f'searching for motif {mtf_min} (originally {motif})')
    
    if not chrs:
        chrs = [str(i) for i in range(1, 23)] + ['X','Y']
    
    matches = []
    
    nrpts = 0
    for chr in chrs:
        print(f"searching chr {chr}")
        
        test_rpts = []
        test_mtfs = []
        
        for rpt in rpts.stream(chr, start = 1e6):
            
            if not rpt.get("type") == "Simple_repeat":
                continue
            
            rpt_mtf = rpt.get("attributes",{}).get("motif","")
            rpt_mtf_min, _, _, _ = bio.get_least_seq_index(rpt_mtf)
            
            if abs(len(rpt_mtf_min) - len(rpt_mtf)) <= max_ham:
                test_rpts.append(rpt)
                test_mtfs.append(rpt_mtf_min)
            
            nrpts += 1
        
        res = repeats.hamming_distance_re(mtf_min, test_mtfs, max_err = max_ham, min_err = 0)
        
        for tr, (r,s,i,d) in zip(test_rpts,res):
            if r:
                matches.append(tr)
                # print(f"errors: {s}, {i}, {d}")
        
        print(f"identified {len(matches)} matches after chr{chr} ({nrpts} total)")
        
    return matches

def get_all_variations(motif):
    
    vars = []
    
    ns = 0
    ni = 0
    nd = 0
    
    for i in range(len(motif)):
        
        currb = motif[i]
        
        #del
        if i > 0:
            new_mtf = motif[:i-1]+motif[i:]      
            new_mtf, _, _, _ = bio.get_least_seq_index(new_mtf)          
            if not new_mtf in vars:
                nd += 1
                vars.append(new_mtf)
        
        for b in 'ATGC':
            
            #sub
            if not b==currb and i > 0:
                new_mtf = motif[:i-1]+b+motif[i:]
                new_mtf, _, _, _ = bio.get_least_seq_index(new_mtf)
                if not new_mtf in vars:
                    ns += 1
                    vars.append(new_mtf)
            
            #ins
            new_mtf = motif[:i]+b+motif[i:]     
            new_mtf, _, _, _ = bio.get_least_seq_index(new_mtf)     
            if not new_mtf in vars:
                ni += 1
                vars.append(new_mtf)
        
    vars_min = []
    
    for v in vars:
        vmin, _, _, _ = bio.get_least_seq_index(v)
        vars_min.append(vmin)
    
    uvars = list(set(vars))
    uvars_min = list(set(vars_min))
    
    print(f"identified {len(vars)} variations of motif {motif}, {len(uvars)} unique, {len(uvars_min)} after minimizing; {ns} from subs, {ni} from ins, {nd} from dels")
    
    return vars

def get_random_sequence(seq_len = 6):
    return "".join([bio.ORDER[random.randint(0, 3)] for n in range(seq_len)])

def plot_chr_motif(gm, chr, motif, num_rows = 5):
    max_disp = 256
    start = 1e6
    length = gm.gene_map.max_indices.get(chr) - 2*start
    chunksz = length // (num_rows * max_disp)
    print(f"chunksz: {chunksz}")
    
    tm_min, _, _, _ = bio.get_least_seq_index(motif)
    
    ptrn = find.consensus_to_re(tm_min)
    
    seq_spec = lambda seq, feats: lambdas._seq_pattern(seq, feats, ptrn)
    seq_spec = lambda seq, feats: lambdas._seq_fuzzy_pattern(seq, feats, ptrn, max_err = 1)
    seq_specs = [seq_spec, "gc"]
    gm.display_chromosomal_quantities(chr, seq_specs, chunksz = chunksz, start=start, minvals = [None, 0.25], maxvals = [None, 0.75])

class RepeatNode:
    
    def __init__(self, motif):
        self.motif = motif
        self.num_samples = 1
        
        self.scs = {}
        self.ics = {}
        self.dcs = {}
        self.ccs = {}

    @property
    def num_connections(self):
        return len(self.scs) + len(self.ics) + len(self.dcs)

    def join(self, other, subs = 0, ins = 0, dels = 0, cont = False):
        
        if subs > 0:
            self.scs[other.motif] = other
            other.scs[self.motif] = self
        elif ins > 0:
            self.ics[other.motif] = other
            other.dcs[self.motif] = self
        elif dels > 0:
            self.dcs[other.motif] = other
            other.ics[self.motif] = self
        elif cont:
            self.ccs[other.motif] = other
        else:
            pass
        pass

class RepeatNetwork:
    def __init__(self):
        self.nodes = {}
        self.motifs = []
        
    def insert(self, motif):
        
        if motif in self.motifs:
            self.nodes[motif].num_samples += 1
            return
        
        new_node = RepeatNode(motif)
        
        res = repeats.hamming_distance_re(motif, self.motifs)
        
        for mtf, (r, s, i, d) in zip(self.motifs, res):
            
            cont = new_node.motif in 2*mtf
            
            if not r and not cont:
                continue
            
            node = self.nodes.get(mtf)
            new_node.join(node, s, i , d, cont)
        
        self.nodes[new_node.motif] = new_node
        self.motifs.append(new_node.motif)

    def print(self):
        
        nodes_srt = sorted(self.nodes.items(), key = lambda k:len(k[0]))
        
        orphs = []
        print(f"RepeatNetwork with {len(self.nodes)} nodes")
        for i, (mtf, node) in enumerate(nodes_srt):
            if node.num_connections < 1:
                orphs.append(node)
                continue
            
            print(f"  Node {i} with motif {node.motif}, {node.num_connections} connections, {node.num_samples} samples")
            print(f"    Joined by s to {",".join(node.scs.keys())}")
            print(f"    Joined by i to {",".join(node.ics.keys())}")
            print(f"    Joined by d to {",".join(node.dcs.keys())}")
            print(f"    Joined by c to {",".join(node.ccs.keys())}")
        
        print(f"Orphan nodes: {",".join([n.motif for n in orphs])}")
    
    def to_network(self):
        
        Gs = nx.Graph()
        Gi = nx.Graph()
        Gd = nx.Graph()
        Gc = nx.Graph()
        
        for mtf, node in self.nodes.items():
            
            Gs.add_node(mtf)
            Gi.add_node(mtf)
            Gd.add_node(mtf)
            Gc.add_node(mtf)
            
            for omtf in node.scs.keys():
                Gs.add_edge(mtf, omtf)
            for omtf in node.ics.keys():
                Gi.add_edge(mtf, omtf)
            for omtf in node.dcs.keys():
                Gd.add_edge(mtf, omtf)
            for omtf in node.ccs.keys():
                Gc.add_edge(mtf, omtf)
        
        return Gs, Gi, Gd, Gc

def build_rpt_network(gm, rpts, chr, max_rpts = None):
    
    net = RepeatNetwork()
    
    nrpt = 0
    for rpt in rpts.stream(chr, start = 1e6):
        
        if not rpt.get("type") == "Simple_repeat":
            continue
        
        mtf = rpt.get("attributes",{}).get("motif","")
        
        mtf_min, _, _, _ = bio.get_least_seq_index(mtf)
        
        net.insert(mtf_min)
        
        nrpt += 1
        if max_rpts and nrpt > max_rpts:
            break
        
        if nrpt % 500 ==0:
            print(f"analyzed {nrpt} repeats")
    
    return net

def main():
    
    gm, rpts = load_genome()
    
    # chrs = None
    chrs = ["1","3","17"]
    chr = chrs[0]
    # test_motif = "GAGCTC"
    # test_motif = "GAGGCAC"
    # test_motif = "GAGGGAG"
    # test_motif = "AGGCGG"
    test_motif = "AATGGAATCG"
    # test_motif = get_random_sequence(seq_len = 7)
    
    # tm_min, _, _, _ = bio.get_least_seq_index(test_motif)
    tm_min = test_motif
    print(f"test motif: {test_motif} ({tm_min})")
    print(bio.consensus_to_SW(tm_min), bio.consensus_to_RY(tm_min), bio.consensus_to_MK(tm_min))
    
    net = build_rpt_network(gm, rpts, chr, max_rpts = None)
    net.print()
    
    Gs, Gi, Gd, Gc = net.to_network()
    
    from matplotlib import pyplot as plt
    for n in range(4):
        graph = [Gs, Gi, Gd, Gc][n]
        lbl = ["subs","ins","dels","conts"][n]
        f, ax = plt.subplots()
        nx.draw_networkx(graph, ax=ax, node_size = 1, node_shape = ',', with_labels = False)
        f.savefig(f"./data/rpts_chr{chr}_G{lbl}.png")
    
    # res = id_repeat(gm, rpts, test_motif, chrs = chrs, max_ham = 1)
    # for r in res:
    #     print(r)
    
    # plot_chr_motif(gm, chr, test_motif)
    
    # for n in range(20):
    #     test_motif = get_random_sequence()
    #     vars = get_all_variations(test_motif)
    #     print(", ".join(vars))


if __name__=="__main__":
    main()

