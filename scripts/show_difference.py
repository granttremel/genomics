



from ggene.database.genome_manager import GenomeManager

from ggene.seqs import align, bio
from ggene.motifs import hmm as ghmm



def load_genome():
    return GenomeManager(load_jaspars = True)


def print_aligned_seqs(seqs_algn, order = [], consensus = ""):
    
    if not order:
        order = list(seqs_algn.keys())
    
    # if not consensus:
    #     consensus = bio.get_consensus(seqs_algn)
    
    for i, k in enumerate(order):
        seq = seqs_algn[k]
        a, b = align.get_colored_aligned_seqs(seq, consensus.upper(), emph_subs = True, color_subs = True)
        if i==0:
            print("{:48}{}".format("cons", b))
        print("{:<48}{}".format(k, a))

def compare_repeats(rname1, rsuff1, rname2, rsuff2, hmms, width = 256):
    
    try:
        fhmm1, merged_cons1 = merge_family(rname1, rsuff1, hmms)
        fhmm2, merged_cons2 = merge_family(rname2, rsuff2, hmms)
        # rhmm1 = hmms.get(rname1)
        # rhmm2 = hmms.get(rname2)
        # r1_cons = rhmm1.consensus.upper()
        # r2_cons = rhmm2.consensus.upper()
    except Exception as e:
        print(f"hmms invalid. options are: {", ".join(hmms.keys())}")
        print(str(e))
        return None
    
    
    algn = align.align_sequences(merged_cons1, merged_cons2)[0]
    algn.print(chunksz = width, emph_match = True, emph_subs = True, emph_indels = True, color_subs = True)
    
    return algn

def merge_family(family_prefix, suffix, hmms):
    cons = {}
    for name, hmm in hmms.items():
        
        if name.startswith(family_prefix):
            
            if suffix and not name.endswith(suffix):
                continue
            
            cons[name]=hmm.consensus.upper()
    
    print(f"collected {len(cons)} conss for family {family_prefix}")
    
    # for k, v in cons.items():
    #     print(v, k)
            
    seqs_aligned = align.align_sequences_muscle(cons)
    
    # for k, v in seqs_aligned.items():
    #     print(v, k)
    
    msa = ghmm.make_msa("merged", seqs_aligned)
    fhmm, _, _ = ghmm.build_hmm(msa)
    
    # print(seqs_aligned)
    
    merged_cons = fhmm.consensus.upper()
    # merged_cons = bio.get_consensus(seqs_aligned)
    
    print_aligned_seqs(seqs_aligned, consensus = merged_cons)
    
    return fhmm, merged_cons

    

def main():
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", type = str, default="AluJr", help = "Beginning characters prefix of first repeat family you want to compare")
    parser.add_argument("-b", type = str, default="AluSc8", help = "Beginning characters prefix of second repeat family you want to compare")
    parser.add_argument("-sfa", type = str, default="", help = "Suffixes of first repeat family you want to compare")
    parser.add_argument("-sfb", type = str, default="", help = "Beginning characters prefix of second repeat family you want to compare")
    parser.add_argument("--width", "-w", type = int, default=256, help = "display width")
    
    args = parser.parse_args()
    
    rn1 = args.a.strip()
    rn2 = args.b.strip()
    rs1 = args.sfa.strip()
    rs2 = args.sfb.strip()
    # rn1, rn2 = args.repeat.split()

    hmms = ghmm.load_hmms()
    hmms = {hmm.name.decode():hmm for hmm in hmms}

    algn =  compare_repeats(rn1, rs1, rn2, rs2,  hmms, width = args.width)




if __name__=="__main__":
    main()
    
