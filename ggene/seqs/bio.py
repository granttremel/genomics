
from ggene.seqs import vocab as gvc
from ggene.seqs.vocab import VOCAB, VOCAB_DNA, VOCAB_RNA

# VOCAB_DNA = "ATGC"
# VOCAB_RNA = "AUGC"

COMPLEMENT_MAP = {'A':'T','T':'A','C':'G','G':'C','U':'A','N':'N'}
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

START_CODONS = ['ATG', 'CTG', 'GTG', 'TTG']  # Standard start codon (Methionine) = AUG and Rare alternative starts = CUG, GUG, UUG
# ALTERNATIVE_START_CODONS = ['CTG', 'GTG', 'TTG']  # 
STOP_CODONS = ['TAA', 'TAG', 'TGA']

ALIASES = {
    'R': 'AG',   # puRine
    'Y': 'CT',   # pYrimidine  
    'S': 'GC',   # Strong (3 H bonds)
    'W': 'AT',   # Weak (2 H bonds)
    'K': 'GT',   # Keto
    'M': 'AC',   # aMino
    'B': 'CGT',  # not A
    'D': 'AGT',  # not C
    'H': 'ACT',  # not G
    'V': 'ACG',  # not T
    'N': 'ACGT', # aNy
}

ALIASES_REV = {
    'A':'RWMDHVN',
    'T':'YWKBDHN',
    'G':'RSKBDVN',
    'C':'YSMBHVN',
}

def reverse_complement(seq):
    return "".join(_reverse(_complement(seq)))

def _complement(seq):
    return [COMPLEMENT_MAP.get(s,s) for s in seq]

def complement(seq):
    return "".join(_complement(seq))

def _reverse(seq):
    return [s for s in reversed(seq)]

def reverse(seq):
    return "".join(_reverse(seq))

def to_rna(seq):
    return seq.replace("T","U")

def to_dna(seq):
    return seq.replace("U","T")

def is_rna(seq):
    return "U" in seq and not "T" in seq

def is_dna(seq):
    return "T" in seq and not "U" in seq

def _convert_vocab(obj, dna = False, rna = False, from_vocab = None, do_keys = True, do_values=True):
    vc = gvc.get_vocab(dna = dna, rna = rna)
    if isinstance(obj, str):
        return gvc._convert_seq_vocab(obj, vc, from_vocab = from_vocab)
    elif isinstance(obj, list):
        return gvc._convert_list_vocab(obj, vc, from_vocab = from_vocab)
    elif isinstance(obj, str):
        return gvc._convert_dict_vocab(obj, vc, from_vocab = from_vocab, do_keys = do_keys, do_values = do_values)
    else:
        return obj

def get_complement_map(vocab = None, rna = False):
    if not vocab:
        if rna:
            vocab = vocab
        else:
            return COMPLEMENT_MAP
        
    return gvc._make_complement_map(vocab)
    
def get_codon_table(vocab=None, rna=False):
    return gvc._convert_vocab(CODON_TABLE, vocab=vocab, rna=rna, do_values = False)

def get_start_stop_codons(vocab=None, rna=False):
    if not vocab and not rna:
        return START_CODONS, STOP_CODONS
    return _convert_vocab(START_CODONS, vocab, rna=rna), _convert_vocab(STOP_CODONS, vocab, rna=rna)

def get_aliases(vocab=None, rna=False):
    
    if not vocab and not rna:
        return START_CODONS, STOP_CODONS
    return _convert_vocab(ALIASES, vocab, rna=rna, do_keys = False), _convert_vocab(ALIASES_REV, vocab, rna=rna, do_values = False)

def get_aa_codons(aa_str, vocab=None, rna=False):
    codons = get_codon_table(vocab=vocab, rna=rna)
    return {a for a,v in codons.items() if v==aa_str.upper()}
    

# simple motifs ig
#idk
splice_donor = 'GG*GURAGU'
splice_branch = 'YURAC'
splice_acceptor = 'YNCAG*G'

