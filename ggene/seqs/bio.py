
import itertools

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

# CODON_TABLE_REV = {v:k for k, v in CODON_TABLE.items()}

ORDER_AA = "ACDEFGHIKLMNPQRSTVWY*"
_aa_chron = "[GA][VD]PS[EL]TRNKQICHFMYW" # trifonov
_aa_chron2 = "[VEG]SIKT[RD]PNFQYMHWC" # 

START_CODONS = ['ATG', 'CTG', 'GTG', 'TTG']  # Standard start codon (Methionine) = AUG and Rare alternative starts = CUG, GUG, UUG
# ALTERNATIVE_START_CODONS = ['CTG', 'GTG', 'TTG']  # 
STOP_CODONS = ['TAA', 'TAG', 'TGA']

ORDER = 'ATGCWRMKYSDHVBN'
# ORDER = "AGTC"
# ORDER = "GATC"

ALIASES = {
    # 'A': 'A',
    # 'T': 'T',
    # 'G': 'G',
    # 'C': 'C',
    'R': 'AG',   # puRine
    'Y': 'TC',   # pYrimidine  
    'S': 'GC',   # Strong (3 H bonds)
    'W': 'AT',   # Weak (2 H bonds)
    'K': 'TG',   # Keto
    'M': 'AC',   # aMino
    'B': 'TGC',  # not A
    'D': 'ATG',  # not C
    'H': 'ATC',  # not G
    'V': 'AGC',  # not T
    'N': 'ATGC', # aNy
}

ALIASES_REV = {v:k for k,v in ALIASES.items()}

ALIASES_INV = {
    'A':'RWMDHVN',
    'T':'YWKBDHN',
    'G':'RSKBDVN',
    'C':'YSMBHVN',
    'R':'DVN',
    'Y':'BHN',
    'S':'BVN',
    'W':'DHN',
    'K':'BDN',
    'M':'HVN',
    'B': 'N',
    'D': 'N',
    'H': 'N',
    'V': 'N', 
    'N': 'N',
}

def get_alias(*bases, exclude = ""):
    
    bases_unalias = [ALIASES.get(b, b) for b in itertools.chain.from_iterable(bases)]
    base_key = tuple(b for b in ORDER[:4] if b in bases_unalias and not b in exclude)
    base_hash = hash(base_key)
    
    for k in ALIASES_REV:
        if hash(tuple(k)) == base_hash:
            return ALIASES_REV[k]
    return "".join(base_key)

def reverse_complement(seq, rna = None):
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

def order_sequences(seqs):
    
    nbs = len(VOCAB)
    base_idx = {ORDER[i]:i for i in range(nbs)}
    
    seq_inds = []
    for seq in seqs:
        seq_len = len(seq)
        seq_ind = sum([base_idx.get(s)*2**(seq_len-i) for i,s in enumerate(seq)])
        seq_inds.append((seq_ind, seq))
    seq_inds_srt = sorted(seq_inds, key=lambda k:k[0])
    return [k[1] for k in seq_inds_srt]

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
    return _convert_vocab(CODON_TABLE, from_vocab=vocab, rna=rna, do_values = False)

def get_start_stop_codons(vocab=None, rna=False):
    if not vocab and not rna:
        return START_CODONS, STOP_CODONS
    return _convert_vocab(START_CODONS, vocab, rna=rna), _convert_vocab(STOP_CODONS, vocab, rna=rna)

def get_aliases(vocab=None, rna=False):
    
    if not vocab and not rna:
        return START_CODONS, STOP_CODONS
    return _convert_vocab(ALIASES, vocab, rna=rna, do_keys = False), _convert_vocab(ALIASES_INV, vocab, rna=rna, do_values = False)

def get_aa_codons(aa_str, vocab=None, rna=False):
    codons = get_codon_table(vocab=vocab, rna=rna)
    cods_ord = order_sequences([a for a,v in codons.items() if v==aa_str.upper()])
    return cods_ord

def get_codon_index(codon):
    nbs = 4
    cdnlen = len(codon)
    base_idx = {ORDER[i]:i for i in range(nbs)}
    return sum([base_idx.get(s)*nbs**(cdnlen-i-1) for i,s in enumerate(codon)])

def index_to_codon(ind):
    nbs = 4
    inds = [(ind//nbs**k)%nbs for k in range(2, -1, -1)]
    return "".join(ORDER[i] for i in inds)

def get_adjacent_codons(aa_str, mutations = "N"):
    
    codons = get_aa_codons(aa_str)
    
    mutes = {}
    for mstr in mutations:
        alistr = ALIASES.get(mstr)
        mutes.update({a:alistr.replace(a,"") for a in alistr})
    
    adj_cods = []
    
    for cd in codons:
        cdl = list(cd)
        for i in range(len(cdl)):
            c = cd[i]
            possible = mutes.get(c)
            if not possible:
                continue
            
            for p in possible:
                cdl[i] = p
                cdlstr = "".join(cdl)
                if cdlstr in codons:
                    pass
                elif cdlstr in adj_cods:
                    pass
                else:
                    adj_cods.append(cdlstr)
            cdl = list(cd)
    coddict = {c:CODON_TABLE.get(c,".") for c in codons}
    adjdict = {c:CODON_TABLE.get(c,".") for c in adj_cods}
    return coddict, adjdict

def get_minimal_alias(b1, b2):
    b1 = ALIASES.get(b1)
    b2 = ALIASES.get(b2)
    bset = set(b1+b2)
    alias_key = ''.join([a for a in ORDER if a in bset])
    return ALIASES_REV.get(alias_key, "N")

def get_minimal_alias_map():
    
    aliases_rev = {v:k for k,v in ALIASES.items()}
    map = {}
    
    for a, avals in ALIASES.items():
        for b, bvals in ALIASES.items():
            if a==b:
                continue
            if (a,b) in map or (b, a) in map:
                continue
            ab = set(list(avals+bvals))
            abb = ''.join([ali for ali in ORDER if ali in ab])
            map[(a,b)] = aliases_rev[abb]
    return map

# simple motifs ig
#idk
splice_donor = 'GG*GURAGU'
splice_branch = 'YURAC'
splice_acceptor = 'YNCAG*G'

