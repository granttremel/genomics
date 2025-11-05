
import random
import string
import numpy as np
import matplotlib.pyplot as plt
import re

from ggene import CODON_TABLE_DNA, CODON_TABLE_RNA, splice_donor, splice_branch, splice_acceptor
from ggene.translate import Ribosome
from ggene.motifs import dyad
from ggene.motifs import utils
from ggene.motifs.utils import Colors
from ggene.motifs.dyad import Dyad, find_all_dyads, frequency_rank_dyads, reverse, reverse_complement
from ggene.genomemanager import GenomeManager

bases = "AUGC"
dna_bases = "ATGC"

RESET = '\033[0m'

CS = Colors.from_specs(text_spec=250, text_bright = True, effect_spec ="")
CD = Colors.from_specs(text_spec="yellow", effect_spec ="")
CL = Colors.from_specs(text_spec="cyan",effect_spec ="")
CB = Colors.from_specs(text_spec="blue",effect_spec ="")
CC = Colors.from_specs(text_spec="cyan",effect_spec ="")

def set_colors(tail=None, dyad=None, loop=None, seq=None, cseq=None, bright=False, background = False, effect = None):
    
    global CS, CD, CL, CB, CC
    
    if tail:
        CS.set(tail, bright=bright, background = background, effect = effect)
    if dyad:
        CD.set(dyad, bright=bright, background = background, effect = effect)
    if loop:
        CL.set(loop, bright=bright, background = background, effect = effect)
    if seq:
        CB.set(seq, bright=bright, background = background, effect = effect)
    if cseq:
        CC.set(cseq, bright=bright, background = background, effect = effect)

set_colors(seq = 174, cseq = 66, background = True)


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
    return seq, "".join(out)

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

def get_cds(gm, extra = 0, all = True):
    
    gene_name = "ZNF804A"
    gene = gm.load_gene(gene_name)
    cdss = gene.get_feature("CDS")
    
    fullseq = []
    fullref = []
    
    for cds in cdss:
        # print(cds.sfid)
        seq = gm.get_feature_sequence(cds, upstream = extra, downstream = extra, as_rna = True)
        refseq = gm.get_feature_sequence(cds, upstream = extra, downstream = extra, personal = False, as_rna = True)
        fullseq.append(seq)
        fullref.append(refseq)
    
    print(f"loaded {cds.type} {cds.sfid} with {len(seq)} BPs")
    
    if not all:
        return fullseq[3], fullref[3]
    
    return fullseq, fullref

def compare_results(seq, input, results):
    
    dyad_len, ds, loop_len, ds_rc = input
    print(ds, loop_len, ds+dyad_len+loop_len)
    
    print("ground truth:")
    d.print()
    print(seq)
    
    print(f"{len(results)} possible dyads identified")
    for result in results:
        dl_res, ds_res, loop_len_res, ds_rc_res = result 
        
        print("search result:")
        print(ds_res, loop_len_res, ds_rc_res)
        d = Dyad(*result, sequence=seq)
        d.print()
        print(seq)
        print(f"result: {result == input}")
    
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
    

def get_intron_motifs(gm, chr = 2, num_introns = 25, num_test = 1):
    
    # introns = gather_features(gm, chr, feature_type = "intron", num_features = num_introns)
    introns = gather_introns(gm, chr, num_introns = num_introns)
    
    intron_seqs = []
    for g in introns:
        for start, end, ex1, ex2 in introns[g]:
            seq = gm.get_sequence(chr, start, end)
            seq = seq.replace("T","U")
            intron_seqs.append(seq)
    
    all_dyads = aggregate_dyads(intron_seqs, max_len = 4096)
    
    seqs = [d.extract_loop() for d in all_dyads] + [d.extract_stems()[0] for d in all_dyads]
    
    freqs = get_frequency(seqs, topk = 10)
    
    temps = aggregate_seqs(seqs, freqs)
    
    new_introns = gather_introns(gm, chr+1, num_introns = num_test)
    
    print(new_introns)
    
    for g in new_introns:
        for start, end, ex in new_introns[g]:
            test_intron_seq = gm.get_sequence(chr+1, start, end)
            test_intron_seq = test_intron_seq.replace("T","U")
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

def compare_intron_cds(gm):
    
    intrs, intrrefs = get_introns(gm)
    cds, cdsref = get_cds(gm)
    fullcds = "".join(cds)
    
    max_len = 2048
    
    intr1 = intrs[0]
    
    print("intron:")
    # test_seq(intr1, max_len = max_len)
    
    # print("cds (ref):")
    # for seq in cdsref:
    #     test_seq(seq, max_len = max_len)
    
def examine_cds(gm):
    
    extra = 32
    cdsseqs, cdsrefs = get_cds(gm, extra = extra)
    frame = 0
    
    colors = {"donor":50, "branch":60, "acceptor":70}
    
    for i in range(len(cdsseqs)):
        cs = cdsseqs[i]
        print(f"exon {i}")
        st, en = find_splices(cs)
        utils.highlight_features(cs,("donor","branch","acceptor"), st, en, colors = colors)
    
    for i in range(len(cdsseqs)):
        cs = cdsseqs[i]
        cs_cds = cs[extra:len(cs) - extra]
        print(f"exon {i}")
        utils.highlight_sequences_in_frame(cs_cds, ["AUC","AUU","AUA"], frame, min_len = 0) # isoleucines idk
        real, prot = translate(cs_cds)
        
        print(prot)
        frame  = (frame + len(cdsseqs[i])) % 3


def find_splices(seq):
    
    dlen = 8
    blen = 5
    alen = 6
    
    spd_re = make_re(splice_donor)
    spb_re = make_re(splice_branch)
    spa_re = make_re(splice_acceptor)
    
    # dres = re.search(seq, )
    dres_iter = re.finditer(spd_re, seq)
    dres = [d.start for d in dres_iter]
    
    bres_iter = re.finditer(spb_re, seq)
    bres = [d.start() for d in bres_iter]
    
    ares_iter = re.finditer(spa_re, seq)
    ares = [d.start() for d in ares_iter]
    
    feat_starts = {}
    feat_ends = {}
    
    feat_starts, feat_ends = utils.make_start_ends("donor", dres, dlen, starts = feat_starts, ends = feat_ends)
    feat_starts, feat_ends = utils.make_start_ends("branch", bres, blen, starts = feat_starts, ends = feat_ends)
    feat_starts, feat_ends = utils.make_start_ends("acceptor", ares, alen, starts = feat_starts, ends = feat_ends)
    
    return feat_starts, feat_ends

def make_re(pattern):
    
    p = pattern.replace("*","")
    p = p.replace("R","[UC]")
    p = p.replace("Y","[AG]")
    p = p.replace("N","[AUGC]")
    return re.compile(p)

def find_start(seq):
    
    for i in range(len(seq)):
        if seq[i:i+3] == "AUG":
            return i
    return -1

def translate(seq, force_start = False, frame = 0):
    start = frame
    if not seq[frame:frame+3] == "AUG" and not force_start:
        start = find_start(seq)
        print(f"starting from {start}")
    aaseq = []
    real = []
    
    for i in range(start, len(seq)//3):
        
        subseq = seq[i:i+3]
        aa = CODON_TABLE_RNA.get(subseq)
        if aa == "*":
            real = aaseq.copy()
            # break
        
        aaseq.append(aa)
    
    return "".join(real), "".join(aaseq)

def try_seq_norm():
    
    s_base = "AUGCGU"
    tab = []
    
    for i in range(len(bases)):
        b0 = bases[i]
        row = []
        row.append(b0)
        for j in range(len(bases)):
            bl = bases[j]
            
            s = b0 + s_base + bl
            
            sout = dyad.normalize_sequence(s)
            soutout = dyad.normalize_sequence(sout)
            char = ""
            if sout == soutout:
                char = "i"
            else:
                char = " "
            
            if sout == s:
                suff = "no rc"
                char += " "
            else:
                suff = "did rc"
                char += "X"
            row.append(char)
            
            print(b0, bl, suff)
            print(s)
            print(sout)
            print()
        tab.append(row)
    
    cols = [" "] + [" "+b for b in bases]
    print(" ".join(cols))
    for row in tab:
        print(" ".join(row))


def test_compare():
    
    s0 = 'AUGCAUGC'
    s1 = 'AUGCNNNNABDALF'
    s2 = 'RYRYNNNNACIAJ'
    s3 = 'RYCYRNNRHIII'
    
    # spp = dyad.to_pyrpur(s0)
    # print("pyrpur:", spp)    
    
    for s in [s1, s2, s3]:
        print(s)
        err = dyad.compare_sequences(s0, s)
        print(err)
        
        # spp = dyad.to_pyrpur(s)
        # print("pyrpur:", spp)    
        
def get_chrom_data(gm):
    
    gtfstream= gm.annotations.streams.get("genes")
    tb = gtfstream.tabix
    
    chrom_data = {}
    
    # for chr in tb.contigs:
        
    #     if len(chr) > 3:
    #         continue
        
    #     print(f"chromosome {chr}")
        
    #     for i in tb.fetch(chr):
            
    #         pass
        
    #     print(f"chromosome {chr}, last row: ")
    #     print(i)
    #     line = gtfstream._parse_line(i)
    #     print(line)
    #     max_index = line.get("end")
    #     step = 1000000
    #     ii = max_index + step
    #     while True:
    #         print(f"position {ii}")
    #         try:
    #             out = tb.fetch(chr, start = ii)
    #             outt = next(out)
    #             print(outt)
    #         except Exception as e:
    #             print(str(e))
            
    #         res = input()
    #         if 'n' in res.lower():
    #             break
            
    #         ii += step
        
    #     chrom_data[chr] = max_index
    
    # print(chrom_data)
    
    introns = []
    
    lastex = None
    print("starting intron search")
    for line in gtfstream.tabix.fetch("2"):
        
        feat = gtfstream._parse_line(line)
        if feat and feat.name:
            if feat.feature_type == "exon":
                
                if not lastex:
                    lastex = feat
                    continue
                
                if not feat.attributes.get("gene_biotype") == "protein_coding":
                    continue
                
                if lastex.name == feat.name:
                    # print(f"intron! at gene {feat.name} after exon {feat.attributes.get("exon_number")}")
                    introns.append((feat.start, lastex.end, feat))
                
                lastex = feat
    
    print("finished intron search")
    
    for start, end, ex in introns:
        pass
        # print(ex.name, start, end)
    
    pass



def main():
    
    gm = load_genome()
    
    hh, _ = get_hammerhead()
    
    hh = hh.replace("U","T")
    hh = hh.replace("Y","[TC]")
    hh = hh.replace("R","[AG]")
    # hh = hh.replace("X","[ATGC]")
    hh_nox = hh.rstrip("X")
    nxs = len(hh) - len(hh_nox)
    
    hh = hh_nox.replace("X","[ATGC]")
    # hh += f"[ATGC]"+"{" + str(nxs) + "}"
    
    print(hh)
    
    # get_chrom_data(gm)
    
    # freqs, temps = get_intron_motifs(gm, num_introns = 5)
    
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
