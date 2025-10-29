
import random
import string
import numpy as np
import matplotlib.pyplot as plt
import re

from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, splice_donor, splice_branch, splice_acceptor
from ggene.translate import Ribosome
from ggene.motifs import dyad
from ggene.motifs import utils
from ggene.motifs.utils import Color
from ggene.motifs.dyad import Dyad, find_all_dyads, frequency_rank_dyads, reverse, reverse_complement
from ggene.genomemanager import GenomeManager

bases = "AUGC"
dna_bases = "ATGC"

def convert_vocab(seq, init_vocab, to_vocab):
    for v, vp in zip(init_vocab, to_vocab):
        seq = seq.replace(v, vp)
    return seq

def make_dyad(dyad_len, loop_len, total_seq_len):
    
    total_dyad_len = 2*dyad_len + loop_len
    
    dyad_start = random.randint(0, total_seq_len - total_dyad_len)
    
    vc = dyad.vocab
    vc_len = len(vc) - 1
    
    dyad = [vc[random.randint(0, vc_len)] for n in range(dyad_len)]
    loop = [vc[random.randint(0,vc_len)] for nl in range(loop_len)]
    full_dyad = dyad + loop + list(reverse_complement(dyad))
    s0 = [vc[random.randint(0, vc_len)] for ns in range(dyad_start)]
    s1 = [vc[random.randint(0,vc_len)] for ns in range(total_seq_len - total_dyad_len - dyad_start)]
    
    return "".join(s0 + full_dyad + s1), dyad_start

def make_random(total_seq_len):
    return "".join([dyad.vocab[random.randint(0,len(dyad.vocab) - 1)] for nl in range(total_seq_len)])

def get_hammerhead():
    
    seq = "YYRRGCCGUUACCURCAGCUGAUGAGCUCCAARAAGAGCGAAACCXRXYAGGUCCUGYAGUAYUGGCYXRXXXXXX"
    subs = {"Y":"C","R":"G","C":"C","G":"G","U":"U","A":"A"}
    out = []
    xbases = "ACGUACGUA"
    
    for s in seq:
        
        ssub = subs.get(s, s)
        if ssub == "X":
            # ssub = random.sample(bases, 1)[0]
            ssub = xbases[0]
            xbases = xbases[1:]
        
        out.append(ssub)
    return "".join(out)

def load_genome():
    return GenomeManager()

def get_linc(geneman, vocab = None):
    
    chr = 2
    start = 207662375
    end = 207679116
    
    seq = geneman.get_sequence(chr, start, end)
    if vocab:
        for b,v in zip(dna_bases, vocab):
            seq = seq.replace(b, v)
    return seq

def get_introns(gm, gene_name = "ZNF804A"):
    
    gene = gm.load_gene(gene_name)
    intrs = gene.get_feature("intron")
    
    seqs = []
    refseqs = []
    for intr in intrs:
        seq = gm.get_feature_sequence(intr, as_rna = True)
        refseq = gm.get_feature_sequence(intr, personal = False, as_rna = True)
        seqs.append(seq)
        refseqs.append(refseq)
    print(f"loaded {intr.type} {intr.sfid} with {len(seq)} BPs")
    return seqs, refseqs
    # return seq, refseq

def gather_introns(gm:GenomeManager, chr, num_introns = 25, num_regions = 5):
    
    introns = {}
    
    gtfstream = gm.annotations.streams.get("genes")
    maxind = gtfstream.max_indices.get(str(chr))
    nfeats_per_region = num_introns // num_regions
    
    while len(introns) < num_introns:
        
        start = random.randint(0, int(0.9*maxind/1000))*1000
        
        feats = gtfstream.stream(chr, start)
        
        last_ex = None
        nfeats = 0
        for i,ex in enumerate(feats):
            
            if not ex or not ex.name:
                continue
            if not ex.feature_type == "exon":
                continue
            if not ex.attributes.get("gene_biotype") == "protein_coding":
                continue
            if ex.name.startswith("LINC"):
                continue
            
            if not last_ex:
                last_ex = ex
                continue
            
            if ex.attributes.get("exon_number") == last_ex.attributes.get("exon_number"):
                continue
            
            if abs(ex.start - last_ex.end) < 500:
                continue
            
            if ex.name == last_ex.name:
                if not ex.name in introns:
                    introns[ex.name] = []
                introns[ex.name].append((last_ex.end+1, ex.start-1, last_ex, ex))
                nfeats += 1
            else:
                print(f"searching gene {ex.name} for introns")
            if len(introns) >= num_introns:
                break
            
            if nfeats > nfeats_per_region:
                print("changing regions")
                break
            
            last_ex = ex
        
        
        if len(introns) >= num_introns:
            break
        
    return introns

def gather_features(gm:GenomeManager, chr, feature_type = "intron", num_features = 25, max_tries = 5):
    
    genes = gm.gene_map.list_genes(chr)
    
    all_features = gm.annotations.stream_by_types(["intron"], chrom = chr)
    feats = {}
    
    nfeats = 0
    nlast = 0
    num_tries = 0
    while nfeats < num_features and num_tries < max_tries:
        gene_name = random.sample(genes, 1)
        gene = gm.assemble_gene(gene_name, chr)
        
        if not gene or not gene.name:
            print(f"skipping {gene}")
            continue
        
        print(f"collecting features from gene {gene_name}")
        gene_feats = gene.get_feature(feature_type)
        feats[gene_name] = gene_feats
        nfeats += len(gene_feats)
        if nfeats == nlast:
            num_tries += 1
        print(f"collected {len(gene_feats)} features")
        
    return feats

def aggregate_dyads(seqs, max_len = 4096):
    
    allsups= []
    nc = 0
    for i, seq in enumerate(seqs):
        nchunk = len(seq) // max_len
        for n in range(nchunk):
            
            if n == nchunk-1:
                chunk = seq[n*max_len:]
            else:
                chunk = seq[n*max_len:(n+1)*max_len]
            
            if not dyad.vocab == bases:
                chunk = convert_vocab(chunk, bases, dyad.vocab)
            
            ds = dyad.find_all_dyads(chunk, 3, max_stem = 4, max_loop = 8, min_loop = 3)
            sups = Dyad.build_superdyads(ds)
            allsups.extend(sups)
            nc += 1
            print(f"collected {len(sups)} dyads from chunk {n} of sequence {i}, {nc*max_len} bp total")
    
    return allsups


def aggregate_seqs(seqs, freqs, max_len = 4096, err_tol = 1):
    
    forms = {}
    top_subseqs = freqs.keys()
    for i, seq in enumerate(seqs):
        nchunk = len(seq) // max_len
        for n in range(nchunk):
            if n == nchunk-1:
                chunk = seq[n*max_len:]
            else:
                chunk = seq[n*max_len:(n+1)*max_len]
            
            subseq_forms, subseq_pos = dyad.find_subsequences_fuzzy(chunk, top_subseqs, err_tol = err_tol)
            
            forms = merge_template_dicts(forms, subseq_forms)
            
    return forms

def merge_template_dicts(temp, temp_new):
    
    for f, ts in temp_new.items():
        if f in temp:
            for tts in ts:
                temp[f].add(tts)
        else:
            temp[f] = ts
    return temp

def get_frequency(seqs, do_rc = False, min_len = 3, topk = 5):
    
    # seqs = [d.extract_loop() for d in dyads]
    freqs, _, maxseq = dyad.frequency_rank(seqs, min_len = min_len, topk = topk, do_rc = do_rc)
    top_subseqs = list(freqs.keys())
    
    for k,v in freqs.items():
        print(f"{k}: {v}")
    print(f"ones: {", ".join([k for k,v in freqs.items() if v==1])}")
    print(f"max sequence: {maxseq}")
    print()
    
    # utils.highlight_sequences(dyads[0].sequence, top_subseqs, do_rc = False)
    # if do_rc:
    #     print("with reverse complement:")
    #     utils.highlight_sequences(dyads[0].sequence, top_subseqs, do_rc = do_rc)

    return freqs

def validate_motifs(intron_seq, freqs, templates, max_len = 4096):
    
    len_seq = len(intron_seq)
    subseqs = list(freqs.keys())
    subseq = intron_seq[(len_seq-max_len)//2:(len_seq+max_len)//2]
    
    print(len(subseq), len(intron_seq))
    
    # temps = list(templates.keys())
    ndisp = 10
    for i in range(len(subseqs)//ndisp):
        ss = subseqs[ndisp*i:ndisp*(i+1)]
        if len(ss) < 1:
            print("emptyy")
            continue
        utils.highlight_sequences(subseq, ss)
    
    # idk
    

def get_intron_motifs(gm, introns, chr = 2, num_introns = 25, num_test = 1, do_loops = True, do_stems = False):
    
    
    intron_seqs = []
    for g in introns:
        for start, end, ex1, ex2 in introns[g]:
            seq = gm.get_sequence(chr, start, end)
            seq = seq.replace("T","U")
            intron_seqs.append(seq)
    
    test_intron_seq = intron_seqs.pop(-1)
    
    all_dyads = aggregate_dyads(intron_seqs, max_len = 4096)
    
    seqs = []
    if do_loops:
        seqs.extend([d.extract_loop() for d in all_dyads])
    if do_stems:
        seqs.extend([d.extract_stems()[0] for d in all_dyads])
    
    freqs = get_frequency(seqs, topk = 10)
    
    temps = aggregate_seqs(seqs, freqs)
    
    # new_introns = gather_introns(gm, chr+1, num_introns = num_test)
    
    # print(new_introns)
    
    validate_motifs(test_intron_seq, freqs, temps)
    
    summarize_motif_results(freqs, temps)
    
    return freqs, temps

def summarize_motif_results(freqs, temps):
    
    print(f"returned {len(freqs)} subsequences and {len(temps)} templates and forms")
    
    print("top subsequences and frequencies:")
    subseqs = sorted([k for k,v in freqs.items()])
    for ss in subseqs:
        print(ss, freqs[ss])
    
    print("subseqence templates:")
    print_tempforms(temps)

def print_tempforms(combined):
    
    out = []
    for pos in combined:
        pos_str = []
        if len(pos) == 1:
            pos_str.append(list(pos.keys())[0])
        else:
            for k, v in pos.items():
                pos_str.append(f"{v}{k}")
        out.append(f"[{".".join(pos_str)}]")
    print("".join(out))


def detect_alu(gm, chr = 1, start = 10e6, step = 1e5):
    
    patterns = [
        # gm.motif_detector.motifs.get("SRP_Alu"),
        # gm.motif_detector.motifs.get("SRP_Alu_stem1"),
        gm.motif_detector.motifs.get("SRP_Alu_stemloop1"),
        gm.motif_detector.motifs.get("SRP_Alu_stemloop2"),
        gm.motif_detector.motifs.get("SRP_Alu_stem2"),
        gm.motif_detector.motifs.get("SRP_Alu_stem3"),
    ]
     
    
    query_freq = 25
    
    all_insts = {p.name:[] for p in patterns}
    
    nfound = 0
    n=0
    while True:
        
        seq = gm.get_sequence(chr, start, start+step)
        rcseq = reverse_complement(seq)
        
        if len(seq) < 1:
            print("ran out of gene")
            input()
            break
        
        nnew ={}
        for pattern in patterns:
            insts = pattern.find_instances(seq)
            rcinsts = pattern.find_instances(rcseq)
            
            all_insts[pattern.name].extend(insts)
            all_insts[pattern.name].extend(rcinsts)
            nnew[pattern.name] = len(insts) + len(rcinsts)
            nfound == len(insts) + len(rcinsts)
        
        found_str = ", ".join([f"{k}: {v}" for k,v in nnew.items()])
        print(f"found {found_str} instances at {chr}:{start}-{start+step}")
        
        start+=step
        
        n += 1
        
        if n % query_freq == 0:
            res = input("keep going?")
            if "n" in res.lower():
                break
            
    
    pass


def convolve(seq1, seq2):
    
    seq_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:seq_len], seq2[seq_len]
    
    start = seq_len // 2
    ips = []
    comp_ips = []
    
    for t in range(-start, start):
        
        sslen = seq_len - abs(t)
        seq1t = seq1[max(t, 0):max(t, 0) + sslen]
        seq2t = seq2[max(-t, 0):max(-t, 0) + sslen]
        summ = 0
        csumm = 0
        
        for sa, sb in zip(seq1t, seq2t):
            
            csb = dyad.COMPLEMENT_MAP.get(sb)
            
            if sa == sb:
                summ +=1
            elif sa == csb:
                csumm += 1
            
        ips.append(summ)
        comp_ips.append(csumm)
    
    return summ, csumm

def test_convolve(gm):
    
    c = 1124574
    seq = gm.get_sequence(1, c-120, c+120)
    
    summ, csumm = convolve(seq, seq)
    
    
    
    
    
    pass

def main():
    
    gm = load_genome()
    
    
    # find_nice_colors()
    
    # bg_color = "\x1b[48;5;212m"
    # f, s = scalar_to_text([1,134,5,6,54,6,5,36,6,47,43,65,73])
    # print(f)
    # print(s)
    
    # detect_alu(gm, chr = 2, start=10e6, step = 5e5)
    
    # get_chrom_data(gm)
    # num_introns = 5
    # chr = str(2)
    # introns = gather_introns(gm, chr, num_introns = num_introns)
    
    # freqs, temps = get_intron_motifs(gm, introns, num_introns = 5)
    
    # freqs2, temps2 = get_intron_motifs(gm, introns, do_loops = False, do_stems = True)
    
    # vocab = dyad._comp_strs[:2] + dyad._comp_strs[4:6]
    # cm = {vocab[i]:vocab[i+1-2*(i%2)] for i in range(len(vocab))}
    # dyad.set_vocab(vocab, cm)
    # print(dyad.vocab, dyad.COMPLEMENT_MAP, dyad.aliases)
    
    # try_seq_norm()
    
    # seq = test_seq(gm)
    
    # compare_intron_cds(gm)
    
    
    
    # examine_cds(gm)
    
    # test_compare()
    
    # find_nice_colors(subseq, num_tries = 50)
    
    

if __name__=="__main__":
    main()
