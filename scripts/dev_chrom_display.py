


from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE

from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.motifs.motif import MotifDetector

import random
from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def test_re_conversion():
    
    testseqs = motif.default_motifs
    
    
    for tsname, testseq in testseqs.items():
    
        testre = find.consensus_to_re(testseq)
        # rctestseq = find.reverse_complement_consensus(testseq)
        # rctestre = find.reverse_complement_re2(testre)
        
        print(f"{tsname}: {testseq} -> {testre}")
        # print(f"reverse complement: {testseq} -> {rctestseq}")
        print()

def test_rcmotif(motif_detector):
    
    
    seq = "TACTATATACACCCACTAC"
    full_seq = seq + reverse_complement(seq)
    print(full_seq)
    
    res = motif_detector.identify(full_seq)
    
    for k, motifs in res.items():
        print(k, motifs)
    
    pass

def check_colors():
    
    teststr = "TACTATATACACCCACTAC"
    
    for c in range(16, 231):
        col = f'\x1b[38;5;{c}m'
        print(c, col + teststr + "\x1b[0m")

def get_test_seq(ind):
    
    if ind == 0:
        return "TGCTTCGGACATTAGAGTCGGCTGGGAATGAAGATGCCTGCAGGCCCCAGAGCACAGGCAGGGCTGAAGTGCGGGCGGCGAGTGGGCTGGGAATGAAGACACCGGAGGGCCCCAGAGTGGAGGCAGGGCTGAGGCACAGGTGGCGAGTGGGCCGGGAATGAAGACGCCGGCAGGCCCCAGAGTGGAGGCAGGGCTGAGGTGCAGGTGGCAGGTGCATCGGCTGGTGCCGTGGTGACGGCG"

def get_random_seq(seq_len):
    return "".join("ATGC"[random.randint(0, 3)] for i in range(seq_len))
    

def test_correlate():
    
    test_seq = get_test_seq(0)
    
    seq_len = len(test_seq)
    
    cres, rcres = process.correlate(test_seq, test_seq)
    
    rmax = max(cres)
    rmaxind = cres.index(rmax) - len(test_seq)//2
    
    print("Forward:")
    
    print(f"max: {rmax:0.3f} at {rmaxind}")
    
    sctxt = draw.scalar_to_text_nb(cres, minval = 0, add_range = True)
    sctxt, _ = draw.add_ruler(sctxt, xmin = 0, xmax = len(test_seq), num_labels = 9, fstr = ".0f")

    for r in sctxt:
        print(r)
    print()
    
    t = rmaxind
    space = " "*abs(t)
    if t < 0:
        seqa = space + test_seq
        seqb = test_seq + space
    else:
        seqa = test_seq + space
        seqb = space + test_seq
    ha, hb = draw.highlight_matching(seqa, seqb, suppress = True)
    
    print(ha.strip())
    print(hb.strip())
    print()
    
    
    rcrmax = max(rcres)
    rcrmaxind = rcres.index(rcrmax) - len(test_seq)//2
    
    if rcrmax < 0.4:
        return
    
    print("Reverse complement:")
    print(f"rc max: {rcrmax:0.3f} at {rcrmaxind}")
    
    rcsctxt = draw.scalar_to_text_nb(rcres, minval = 0, fg_color=  53-12, add_range = True)
    rcsctxt, _ = draw.add_ruler(rcsctxt, xmin = 0, xmax = len(test_seq), num_labels = 9, fstr = ".0f")
    
    for r in rcsctxt:
        print(r)
    print()
    
    t = rcrmaxind
    
    space = " "*abs(t)
    if t < 0:
        seqa = space + test_seq
        seqb = test_seq + space
    else:
        seqa = test_seq + space
        seqb = space + test_seq
    ha, hb = draw.highlight_matching(seqa, seqb, suppress = True)
    
    print(ha)
    print(hb)
    
    print()

def test_dyads():
    
    test_seq = get_test_seq(0)
    
    seeds, ds = dyad.find_all_dyads_build(test_seq, 5, min_loop = 3)
    
    print(f"identified {len(ds)} dyads, max stem {max([d.stem_length for d in ds])}, max loop {max(d.loop_length for d in ds)}")
    
    # for d in seeds:
    #     print(d.to_tuple())
    
    for d in seeds:
        if d.loop_length < 6:
            continue
        if d.stem_start < 40 or d.end_position > 100:
            continue
        d.print(total_len = len(test_seq))
    
    # stems = list(set([d.extract_stems()[0] for d in ds] + [d.extract_stems()[1] for d in ds]))
    
    # draw.highlight_sequences(test_seq, stems, suppress = False)
    
    # hd = draw.highlight_dyads(test_seq, ds)
    # print(hd)

def get_random_sequence(seq_len, minval = -100, maxval = 100, var = 1, bias = 0, nint = 1):
    
    rand_seq_d = [(2*random.random()-1)*var+bias for i in range(seq_len)]
    
    for n in range(nint):
        rand_seq = [sum(rand_seq_d[:i]) for i in range(seq_len)]
        rand_seq_d = rand_seq
    
    _max= max(rand_seq)
    _min = min(rand_seq)
    _ran = _max - _min
    
    tgt_ran = maxval - minval
    
    rand_seq_sc = [tgt_ran*(r-_min)/_ran+minval for r in rand_seq]
    return rand_seq_sc

def find_nice_colors(seq_len = 256):
    
    # seq = list(range(-100, 0, 14)) + list(range(0, -50, 10)) + list(range(-50,25,9)) + list(range(25,100,18)) + list(range(100, 50, -19)) + [50]
    n = 0
    while True:
        seq = get_random_sequence(seq_len, var = 4)
        offset = 24
        # fc = 170
        bc = 234
        # bc = random.randint(20, 230-offset)
        fc = random.randint(20,230)
        # fc = min(bc+offset, 230)
        
        if n == 0:
            fc = 236
            bc = 244
        
        res = draw.scalar_to_text_nb(seq, bg_color = bc, fg_color = fc, bit_depth = 16)
            
        print(f"bg color = {bc}, fg color = {fc}")
        for r in res:
            print(r)
        print()
        n += 1
        
        if n%5 == 0:
            res=input()
            if 'n' in res.lower():
                break
    
    return bc

# def random_seq():
#     "AAGCATGGACCATTCCTTCAGGATGGAGCCTGCCAAGTGTGGACCATTCCTTCAGGATGCAGCCAGGTAAGCATCAGCCATTCCTTCATAATGCGGCCAGGTAAGCAT"
#     "CCTTCAGGATGGAGCCTGCCAAGTGTGGACCATTCCTTCAGGATGCAGCCAGGTAAGCATCAGCCATTCCTTCATAATGCAGCCAGGTAAGCAT CAGCCATTCCTT"
#     "GATTGGTGCGTTTACAAACCTTGAGCTAGACACAA GGTGCTGATTGGTGCGTTTACCAACCTTGAGCTAGACACAGGGTGCTGATT"
#         "GGTGCATTTACAATCCTTTAGCTAGACATAAAAGTTCTCCAAGTCCCCACCAGATTAGCTAGATACAGAGTGCTTATTTGGTGCTTCCATGATCCCCGAGCTAGATACAGAGTGATGATTGGTGTA"
#     pass

def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    # test_re_conversion()
    
    # mdet = MotifDetector()
    # mdet.setup_default_motifs()
    
    # test_rcmotif(mdet)
    # check_colors()
    
    # test_correlate()
    
    # test_dyads()
    
    # find_nice_colors()
    
    # seq_len = 8*256
    
    # seq = get_random_seq(256)
    # print("autocorrelation")
    
    # data, rcdata = process.correlate(seq, seq)
    
    # dbars = draw.scalar_to_text_16b(data, minval = 0, fg_color = 131)
    # rcdbars = draw.scalar_to_text_16b(rcdata, minval = 0, fg_color = 91, flip=True)
    
    # for r in dbars:
    #     print(r)
    # for r in rcdbars:
    #     print(r)
    # print()
    
    # import numpy as np
    # print(f"mean = {np.mean(data):0.3f}, rcmean = {np.mean(rcdata):0.3f}")
    # print()
        
    # seqa = get_random_seq(seq_len)
    # seqb = get_random_seq(seq_len)
    
    # print("crosscorrelation")
    
    # data, rcdata = process.correlate(seqa, seqb)
    
    # dbars = draw.scalar_to_text_16b(data[:256], minval = 0, fg_color = 131)
    # rcdbars = draw.scalar_to_text_16b(rcdata[:256], minval = 0, fg_color = 91, flip=True)
    
    # for r in dbars:
    #     print(r)
    # for r in rcdbars:
    #     print(r)
    # print()
    
    # print(f"mean = {np.mean(data):0.3f}, rcmean = {np.mean(rcdata):0.3f}")
    
    # print()
    # return
    
    gm = load_genome()
    gm.motif_detector.motifs.clear()
    gm.motif_detector.setup_default_motifs(class_names = ["msat"])
    
    test_chr = 2
    
    gene_start = 184598470
    gene_end = 184939492
    gene_chunks = 128
    
    max_disp = 256
    num_chunks = 3*max_disp
    # chunksz = int(64e3)
    # start = int(16e6)
    chunksz = (gene_end - gene_start)//gene_chunks
    start = gene_start - chunksz*(num_chunks - gene_chunks)//2
    length = num_chunks * chunksz
    
    # zoom
    # start = int(86e6)
    # chunksz = int(8e3)
    # length = 16e6
    
    # ran = gm.gene_map.max_indices[str(test_chr)] - start
    # num_chunks = ran//chunksz
    
    # gm.display_chromosomal_quantity(test_chr, "cpg", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    # gm.display_chromosomal_quantity(test_chr, "cds_pct", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    
    # gm.display_chromosomal_quantity(test_chr, "motifs", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    
    # gm.display_chromosomal_quantity(test_chr, "genes", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    # gm.display_chromosomal_quantity(test_chr, "penta_repeats", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    # gm.display_chromosomal_quantity(test_chr, "hammerhead_st1", chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    
    gm.display_chromosomal_quantities(test_chr, ["cg","cpg"], chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    # gm.display_chromosomal_quantities(test_chr, ["genes","penta_repeats"], chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    # gm.display_chromosomal_quantities(test_chr, ["genes","quad_repeats","penta_repeats","hexa_repeats"], chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    
    # gm.display_chromosomal_quantities(test_chr, ["genes","hammerhead_st1"], chunksz = chunksz, start = start, max_disp = max_disp, length = length)
    pass    


if __name__=="__main__":
    main()

