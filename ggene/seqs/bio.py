
import random
import itertools
import numpy as np

from ggene.seqs import vocab as gvc
from ggene.seqs.vocab import VOCAB, VOCAB_DNA, VOCAB_RNA

# VOCAB_DNA = "ATGC"
# VOCAB_RNA = "AUGC"
OTHER = "Ψ"

# reflection about S-W axis
COMPLEMENT_MAP = {
    'A':'T','T':'A',
    'C':'G','G':'C',
    
    'R':'Y','Y':'R',
    'S':'S','W':'W',
    'K':'M','M':'K',
    
    'B':'V','V':'B',
    'D':'H','H':'D',
    
    'N':'N'}

# reflection about R-Y axis
TRANSITION_MAP = {
    "A":"G","G":"A",
    "C":"T","T":"C",
    "R":"R","Y":"Y",
    "S":"W","W":"S",
    "K":"M","M":"K",
    "N":"N"}

# excluding the transversions already described by complement
# reflection about K-M axis
TRANSVERSION_MAP = {
    "A":"C","C":"A",
    "G":"T","T":"G",
    "R":"Y","Y":"R",
    "S":"W","W":"S",
    "K":"K","M":"M",
    "N":"N"}

# ???
_wobble = {
    ('A','G'):'+',
    ('T','G'):'⋅'
}

START_CODONS = ['ATG', 'CTG', 'GTG', 'TTG']  # Standard start codon (Methionine) = AUG and Rare alternative starts = CUG, GUG, UUG
STOP_CODONS = ['TAA', 'TAG', 'TGA']

ORDER = 'ATGCWRMKYSDHVBN-'

ALIASES = {
    'A':'A',
    'T':'T',
    'G':'G',
    'C':'C',
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

ALIASES_FULL = {
    'A':'A',
    'T':'T',
    'G':'G',
    'C':'C',
    'R': 'AGR',   # puRine
    'Y': 'TCY',   # pYrimidine  
    'S': 'GCS',   # Strong (3 H bonds)
    'W': 'ATW',   # Weak (2 H bonds)
    'K': 'TGK',   # Keto
    'M': 'ACM',   # aMino
    'B': 'TGCYSKB',  # not A
    'D': 'ATGRWKD',  # not C
    'H': 'ATCYWMH',  # not G
    'V': 'AGCRSMV',  # not T
    'N': 'ATGCRYSWKMBDHVN', # aNy
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

ORDER_AA = "ACDEFGHIKLMNPQRSTVWY*"
_aa_chron = "[GA][VD]PS[EL]TRNKQICHFMYW" # trifonov
_aa_chron2 = "[VEG]SIKT[RD]PNFQYMHWC" # 

AA_CLASSES = {
    "basic":"HKR",
    "acidic":"DE", 
    "polar":"NQST", 
    "nonpolar":"AFILMVWY",
    "aromatic":"FWY",
    "phosphoryl":"HSTY", # H kind weird idk
    # "special":"CGP",
}

AA_CLASSES_REV = {
    "A":["nonpolar"], "C":[], "D":["acidic"], "E":["acidic"],
    "F":["nonpolar","aromatic"], "G":[],
    "H":["basic", "phosphoryl"], "I":["nonpolar"],
    "K":["basic"], "L":["nonpolar"], "M":["nonpolar"],
    "N":["polar"], "P":[], "Q":["polar"], "R":["basic"],
    "S":["polar","phosphoryl"], "T":["polar","phosphoryl"],
    "V":["nonpolar"], "W":["nonpolar","aromatic"],
    "Y":["nonpolar","aromatic","phosphoryl"],
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

def _transvert(seq):
    return [TRANSVERSION_MAP.get(s,s) for s in seq] 

def transvert(seq):
    return "".join(_transvert(seq))

def _transition(seq):
    return [TRANSITION_MAP.get(s,s) for s in seq] 

def transition(seq):
    return "".join(_transition(seq))

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

def is_seq(seq):
    voc_comb = set(VOCAB_DNA + VOCAB_RNA)
    for b in seq:
        if not b in voc_comb:
            return False
    return True
    
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
        return ALIASES, ALIASES_INV
    return _convert_vocab(ALIASES, vocab, rna=rna, do_keys = False), _convert_vocab(ALIASES_INV, vocab, rna=rna, do_values = False)

def get_random_sequence(length):
    seq = [random.choice(VOCAB) for i in range(length)]
    return "".join(seq)

def get_random_aa_sequence(length):
    codons = ["".join([random.choice(VOCAB) for ii in range(3)]) for i in range(length)]
    seq = [CODON_TABLE.get(cd) for cd in codons]
    return "".join(seq)

def get_aa_codons(aa_str, vocab=None, rna=False):
    codons = get_codon_table(vocab=vocab, rna=rna)
    cods_ord = order_sequences([a for a,v in codons.items() if v==aa_str.upper()])
    return cods_ord

def get_seq_index(seq):
    nbs = 4
    seq_len = len(seq)
    base_idx = {ORDER[i]:i for i in range(nbs)}
    return sum([base_idx.get(s)*nbs**(seq_len-i-1) for i,s in enumerate(seq)])

def index_to_seq(ind, seq_len = 4):
    nbs = 4
    if ind == 0:
        return ORDER[0]*seq_len
    inds = [(ind//nbs**k)%nbs for k in range(seq_len-1, -1, -1)]
    return "".join(ORDER[i] for i in inds)

def get_seq_index_abs(seq):
    seq = ORDER[0] + seq
    seq_len = len(seq)
    fullind = get_seq_index(seq)
    return fullind + sum([4**n for n in range(1, seq_len-1)])

def get_seq_least_index_ugh(seq):
    
    for b in VOCAB:
        sym = b
        nsym = seq.count(sym)
        if not nsym:
            pass
        else:
            break
    
    if sym == VOCAB[-1]:
        return seq, get_seq_index_abs(seq)
    
    if nsym == 1:
        opt_ind = seq.index(sym)
        opt_seq = seq[opt_ind + 1:] + seq[:opt_ind+1]
        return opt_seq, get_seq_index_abs(seq)
    
    seq_len = len(seq)
    init_seq = seq
    ind = 0
    arg_min = -1
    min_seq_ind = np.inf
    done = False
    while not done:
        symind = seq_len - seq.index(sym) - 1
        seq = seq[symind:] + seq[:symind]
        ind += symind
        
        seq_ind = get_seq_index_abs(seq)
        if seq_ind < min_seq_ind:
            arg_min = ind
            min_seq_ind = seq_ind
        
        if ind >= len(seq):
            done = True
            break
        
    return init_seq[arg_min:] + init_seq[:arg_min], min_seq_ind

def get_least_seq_index(seq, do_rc = True, do_rev = False, do_comp = False):
    
    if len(seq) < 1:
        return None, None, None, None
    
    seq_len = len(seq)
    arg_min = -1
    min_seq_ind = np.inf
    
    for i in range(len(seq)):
        test_seq = seq[i:] + seq[:i]
        test_seq_ind = get_seq_index(test_seq)
        if test_seq_ind < min_seq_ind:
            arg_min = i
            min_seq_ind = test_seq_ind
    
    cseq = rcseq = revseq = ""
    cind = rcind = revind = np.inf
    
    if do_rc:
        rcseq, rcind, rc_arg_min, _ = get_least_seq_index(reverse_complement(seq), do_rc = False, do_rev = False, do_comp = False)
    if do_comp:
        cseq, cind, c_arg_min, _ = get_least_seq_index(complement(seq), do_rc = False, do_rev = False, do_comp = False)
    if do_rev:
        revseq, revind, rev_arg_min, _ = get_least_seq_index(reverse(seq), do_rc = False, do_rev = False, do_comp = False)
    
    min_ind = min(test_seq_ind, rcind, cind, revind)
    
    code = [i == min_ind for i in [test_seq_ind, rcind, cind, revind]]
    code_int = sum([i*2**n for n,i in enumerate(code)])
    
    if rcind == min_ind:
        return rcseq, rcind, rc_arg_min, code_int
    elif cind == min_ind:
        return cseq, cind, c_arg_min, code_int
    elif revind == min_ind:
        return revseq, revind, rev_arg_min, code_int
    else:
        return seq[arg_min:] + seq[:arg_min], test_seq_ind, arg_min, code_int

def index_to_seq_abs(ind):
    seq_len = 0
    ip = ind+1
    while ip > 0:
        seq_len += 1
        ip -= 4**seq_len
    ipp = ip + 4**seq_len - 1
    seq = index_to_seq(ipp, seq_len = seq_len+1)
    if seq:
        return seq[1:]
    
    else:
        return ""
    
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

def convert_to_aliases(seq, aliases, default = "N"):
    """
    set default to "" to preserve bases with no associated alias
    
    """
    alias_map = get_alias_map(aliases)
    
    outseq = []
    
    for b in seq:
        _def = default if default else b
        newb = alias_map.get(b, _def)
        outseq.append(newb)
    
    return "".join(outseq)

def get_alias_map(aliases):
    
    alias_map = {}
    
    for b in ORDER:
        
        alis = ALIASES_INV.get(b,"")
        
        for a in alis:
            if a in aliases:
                alias_map[b] = a
                break
        
    return alias_map

def get_minimal_alias(*bs, default = "N"):
    bset = set()
    for b in bs:
        ba = ALIASES.get(b)
        bset.add(ba)
    alias_key = ''.join([a for a in ORDER if a in bset])
    return ALIASES_REV.get(alias_key, default)

def get_minimal_alias2(*bs):
    bset = set(ORDER)
    for b in bs:
        bset = bset.intersection(ALIASES_FULL.get(b))
    
    for b in ORDER:
        if b in bset:
            return b
    return ""    
    

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

def hamming(seqa, seqb, mutation_vals = {}):
    
    if not mutation_vals:
        mutation_vals = {a:1 for a in ALIASES}
    
    merged = merge(seqa, seqb)
    
    hdist = 0
    for mk in mutation_vals:
        hdist += mutation_vals.get(mk) * merged.count(mk)
    return hdist


##########3 consensuses ################

def consensus_entropy(consensus:str):
    
    if len(consensus) < 1:
        return 1.0
    
    consensus_strp = consensus.lstrip("N").rstrip("N")
    
    entr = 0
    for b in consensus_strp:
        alis = ALIASES.get(b)
        p = 1/len(alis)
        entr += -p*np.log(p)
    return entr / len(consensus)

def consensus_likelihood(consensus):
    return np.exp(consensus_logprob(consensus))

def consensus_logprob(consensus):
    
    logp = 0
    for b in consensus:
        alis = ALIASES.get(b)
        logp -= np.log(4-len(alis))
    return logp

def consensus_oddsratio(consensus, samps, obs):
    
    p_obs = obs / samps
    
    p_cons = np.exp(consensus_logprob(consensus))
    
    # mean_rand = samps * p_cons
    # var_rand = samps * p_cons * (1-p_cons)
    
    return p_obs / p_cons

def consensus_to_minimal(cns):
    
    smin = []
    
    for b in cns:
        smin.append(ALIASES.get(b,b)[0])
    
    return "".join(smin)

def merge(seqa, seqb, default = 'N'):
    outseq = []
    for sa, sb in zip(seqa, seqb):
        if sa==sb:
            outseq.append(sa)
        else:
            # mt = ALIASES_REV.get(sa+sb, ALIASES_REV.get(sb + sa, default))
            mt = get_minimal_alias(sa, sb, default=default)
            outseq.append(mt)
    return "".join(outseq)

def merge_all(seqs, default = 'N'):
    
    cons = merge(seqs.pop(0), seqs.pop(0), default = default)
    
    for seq in seqs:
        cons = merge(cons, seq, default = default)
    
    return cons

def merge_with_gaps(seqa, seqb, default = 'N'):
    
    outseq = []
    for sa, sb in zip(seqa, seqb):
        if sa==sb:
            outseq.append(sa)
        elif sa == "-" or sb == "-":
            outseq.append("-")
        else:
            mt = get_minimal_alias(sa, sb, default=default)
            outseq.append(mt)
    
    return "".join(outseq)

def get_consensus(aligned_seqs, min_p_ratio = 2):
    
    num_seqs = len(aligned_seqs)
    names = list(aligned_seqs.keys())
    seq_len = len(aligned_seqs[names[0]])
    
    cons = []
    
    for p in range(seq_len):
        
        p_bases = {k:0 for k in 'ATGC-'}
        
        for nm in names:
            base = aligned_seqs[nm][p]
            p_bases[base] += 1/num_seqs
        
        n_bases = 5

        while n_bases > 1:
            min_base, min_p = min(p_bases.items(), key=lambda x: x[1])
            threshold = 1 / n_bases / min_p_ratio

            if min_p < threshold:
                del p_bases[min_base]
                n_bases -= 1
            else:
                break
        
        remaining = list(p_bases.keys())
        
        if remaining == ['-']:
            b_ali = '-'
        else:
            b_ali = get_alias(*remaining)
        cons.append(b_ali)

    return "".join(cons)
        
def is_consensus(seq, consensus):
    
    for ba, bc in zip(seq, consensus):
        if ba == bc:
            pass
        elif ba in ALIASES.get(bc, bc):
            pass
        else:
            return False
    return True

def consensus_to_alias(cons, aliases):
    
    aliases = aliases
    
    outseq = []
    
    ali_dict = {}
    
    for ali in aliases:
        refs = ALIASES_FULL.get(ali)
        for r in refs:
            if r not in ali_dict:
                ali_dict[r] = ali
    
    for b in cons:
        outseq.append(ali_dict.get(b, 'N'))
        
    return "".join(outseq)

def consensus_to_SW(cons):
    return consensus_to_alias(cons, "SW")

def consensus_to_RY(cons):
    return consensus_to_alias(cons, "RY")

def consensus_to_MK(cons):
    return consensus_to_alias(cons, "MK")

def sample_consensus(cons):
    cs = [random.choice(ALIASES.get(b, "X")) for b in cons]
    return "".join(cs)

# simple motifs ig
#idk
splice_donor = 'GG*GURAGU'
splice_branch = 'YURAC'
splice_acceptor = 'YNCAG*G'

