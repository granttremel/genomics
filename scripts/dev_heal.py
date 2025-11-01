
from ggene import seqs, draw
from ggene.seqs import convolve, convolve_longest_subseq, heal, reverse, reverse_complement, COMPLEMENT_MAP

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()


def get_test_seq():
    return "AATGAATGAATGACTGAATGAATG"



def main():
    
    bg, fg = draw.get_color_scheme("dusty")
    
    seq = get_test_seq()
    subseq = seq[:4]
    
    cv, rccv = seqs.convolve(seq, subseq)
    res = draw.scalar_to_text_nb(cv)
    for r in res:
        print(r)
        
    found = seqs.find_subsequence(seq, subseq)
    comp = seqs.compare_sequences(seq, subseq)
    
        
    
    
    pass


if __name__=="__main__":
    main()



