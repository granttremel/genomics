

import re
import random

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align, vocab
from ggene.seqs.bio import is_consensus, reverse_complement, reverse, complement, get_aa_codons, get_adjacent_codons, hamming
from ggene.seqs.bio import ORDER, ORDER_AA, CODON_TABLE, ALIASES, ALIASES_FULL

from ggene.motifs import motif
from ggene.motifs import dyad
from ggene.motifs.motif import MotifDetector

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()

def extract_loops(gm, chrs, start, end=None, chunksz = 4e6, ptrn_spec = "SRP_S_stemloop_6", check_dyads = False):
    
    ptrn = gm.motif_detector.motifs.get(ptrn_spec).pattern
    print(f"pattern: {ptrn}")
    ptrn = re.compile(ptrn)
    
    
    if isinstance(chrs,str):
        chrs = [chrs]
    
    all_loops = []
    
    for chr in chrs:
        
        if not end:
            en = gm.gene_map.max_indices.get(str(chr))
        else:
            en = end
        st = start
        
        print(f"starting chromosome {chr}")
        
        chr_loops = []
        
        while st < en:
            
            print(f"starting {chr}:{st}-{st+chunksz}")
            seq = gm.get_sequence(chr,st, st + chunksz)
            rcseq = reverse_complement(seq)
            strands = ["+","-"]
            
            new_loops = []
            for i, seqq in enumerate([seq, rcseq]):
                
                matches = ptrn.finditer(seqq)
                
                for m in matches:
                    if m.groups():
                        loops = list(m.groups())
                    else:
                        loops = None
                    if i==1:
                        loop_start = len(seq) - m.start(0) + 1 + st
                    else:
                        loop_start = m.start(0) + st
                    new_loops.append((chr, loop_start, strands[i], loops))
                
            print(f"found {len(new_loops)} loops in {chr}:{st}-{st+chunksz}")
            print()
            for n in range(min(5, len(new_loops))):
                _chr, _st, strand, loop = new_loops[n]
                print(loop, f"(start={_st}, {strand})")
                if check_dyads:
                    check_loop_dyads(_chr, _st, loop, strand)
            
            chr_loops.extend(new_loops)
            
            st += chunksz
        print(f"found {len(chr_loops)} on chromosome {chr}")
        
        all_loops.extend(chr_loops)
    print(f"found {len(all_loops)} loops in total!")
    
    return all_loops
    
def check_loop_dyads(chr, st, loop, strand):
            
    ds = dyad.find_all_dyads(loop, min_stem = 4, min_loop = 3, err_tol = 1)
    # rawds, ds = dyad.find_all_dyads_build(loop, min_stem = 4, min_loop = 3)
    
    # ds = [d for d in ds if d.end_position < len(loop)-1 and d.stem_start > 0]
    
    maxd = None
    maxlen = 0
    for d in ds:
        # if d.total_length > maxlen:
        if d.stem_length > maxlen:
            maxd = d
            maxlen = d.stem_length
    
    if maxd:
        maxd.print(total_len = len(loop))
        print()

def test_pattern(gm):
    
    ptrn = gm.motif_detector.motifs.get("SRP_S_stemloop_6").pattern
    ptrn = re.compile(ptrn)
    # test_seq = "CAATATGGTGACCTCCCGCGGGGGACCACCAGGTTG"
    test_seq =  "CAATATGGTGACCTCCCGCCATCGGGGGACCACCAGGTTG"
    test_seq2 = "CAATATGGTGACCTCCCGTGGGACGGGGGACCACCAGGTTG"
    test_seq3 = "CAATATGGTGACCTCCCGACTCGGGGGACCACCAGGTTG"
    
    for seq in [test_seq, test_seq2, test_seq3]:
        print(f"test seq: {seq}")
        
        rcseq = reverse_complement(seq)
        
        match= ptrn.search(seq)
        if match:
            print(match.groups())
        else:
            print("no match found (fwd)")
        
        match= ptrn.search(rcseq)
        if match:
            print(match.groups())
        else:
            print("no match found (rc)")
        

def main():
    
    bg, fg = draw.get_color_scheme("test")

    gm = load_genome()
    gm.motif_detector.setup_default_motifs()

    print(list(gm.motif_detector.motifs.keys()))

    # chr = 7
    chrs = list(range(1, 22+1)) + ["X","Y"]
    # chrs = [1, 2, 3, 5, "X"]
    # chrs = ["X"]
    # chrs = [1]
    start = 1.6e6
    # end = 8*start
    end = None
    
    ptrn_name = "SRP_S_stemloop_6"
    # ptrn_name = "SRP_Alu_stemloop_3"
    # ptrn_name = "SRP_Alu_stemloop_4"
    # ptrn_name = "SRP_S_stemloop_8a" # 15 on chr 1
    # ptrn_name = "SRP_S_stemloop_5f"
    # ptrn_name = "SRP_S_stemloop_5e"
    
    # ptrn_name = "telomerase_pseudoknot_hTR44-184"
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-1" # 9 on chr1, 93 total
    ptrn_name = "telomerase_pseudoknot_hTR44-184-2" # 0 on chr1, 1 total (on chr3:129600000-193600000)
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-3" # 0 on chr1, 1 total (on chr3:129600000-193600000)
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-4" # 0 on chr1, 1 total (on chr3::129600000-193600000)
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-5" # 0 on chr1, 1 total (on chr3::129600000-193600000)
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-6" # 0 on chr1, 6 total
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-7" # 0 on chr1, 1 total (on chr3::129600000-193600000)
    
    # ptrn_name = "telomerase_pseudoknot_hTR44-184-loop5" # 1 total (on chr3::129600000-193600000) 
    
    # ptrn_name = "telomerase_pseudoknot_hTR171-184" # 9 on chr1, 140 total
    
    chunksz = 64e6

    loops = extract_loops(gm, chrs, start, end=end, chunksz = chunksz, ptrn_spec = ptrn_name)
    
    # num_check = 15
    
    # for chr, st, loop, strand in random.sample(loops, k=num_check):
        
    #     ds = dyad.find_all_dyads(loop, min_stem = 4, min_loop = 3)
    #     print(f"found {len(ds)} dyads on loop at {chr}:{st} ({strand})")
        
    #     ds = [d for d in ds if d.end_position < len(loop) and d.stem_start >= 0]
    #     maxd = None
    #     maxlen = 0
    #     for d in ds:
    #         if d.total_length > maxlen:
    #             maxd = d
    #             maxlen = d.total_length
        
    #     if maxd:
    #         maxd.print(total_len = len(loop))
        
    
    # test_pattern(gm)
    

        
if __name__=="__main__":
    main()
