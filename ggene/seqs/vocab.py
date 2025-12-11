
import string
import itertools
import random

VOCAB_DNA = "ATGC"
VOCAB_RNA = "AUGC"

VOCAB = VOCAB_DNA

"""
a good vocab would:
be asymmetric about horizontal axis, so complement is like visually rotating it. also must have 
two orienations of the character
visually distinct without being too busy or dense?
counterpunctal forms
may need two rotated versions for upper sequence and lower sequence
e.g.:
V and ^ (if they were more symmetric)
T and U ironically
A and W ?

⊓ and ⊔
Ո and Ս
Ψ ψ Ϙ Ͳ
Ш Ϣ
⋿
∈ and ∋
≪ and ≫
⋒ and ⋓
╥ and ╨
╤ ╧

ATGC = ┬╤╦╥
TACG = ╧┴╨╩

or ╪ ╫
or ╤ ╥ ╦ ╧ ╨ ╩
"""

_addl_syms = "ЂЄЉЊЋБГ"

_good_syms = [
    (9632, 9727, 1), # geo
    (9484, 9544, 4), # box, fine
    (8704, 8959, 1), # math
    (8960, 9210, 1), # technical
    (11360, 11391, 1), # latin extended c
    (9487, 9547, 4), # box, bold
    (9552, 9580, 1), # box, double
    (9728, 9983, 1), # misc symbols
    (9985, 10175, 1), # dingbats
]

def get_symbol(i):
    
    if i < len(string.ascii_uppercase):    
        return string.ascii_uppercase[i]
    else:
        i-=len(string.ascii_uppercase)
        
        for start, stop, step in _good_syms:
            
            n = (stop+1-start)//step
            if i > n:
                i-= n
                continue
            else:
                return chr(start + step * i)
    
    return ""

def screen_symbols(start, end):
    
    syms = []
    for i in range(start, end):
        
        print(i, chr(i))
        res = input("\nkeep?\n")
        
        if 'y' in res.lower():
            syms.append(i)

    print(f"new_symbols=[{",".join(syms)}]")
    return syms


def extend_vocab(vocab, k_tot):
    
    vocab = list(vocab)
    
    i = 0
    # for c in string.ascii_uppercase:
    while len(vocab) < k_tot:
        c = get_symbol(i)
        if c in vocab:
            continue
        
        vocab.append(c)
        i += 1
    
    return "".join(vocab)

def set_vocab(new_vocab):
    global VOCAB
    new_vocab_unique = list(set(new_vocab))
    if len(new_vocab)%2==0 and len(new_vocab) == len(new_vocab_unique):
        VOCAB = new_vocab
    return VOCAB

def set_rna():
    return set_vocab(VOCAB_RNA)
def set_dna():
    return set_vocab(VOCAB_DNA)

def get_vocab(dna = False, rna = False):
    if dna:
        return VOCAB_DNA
    elif rna:
        return VOCAB_RNA
    else:
        return VOCAB

def reset_vocab():
    global VOCAB
    VOCAB = VOCAB_DNA

def infer_vocab(seq, *vocabs):
    if len(vocabs) < 1:
        vocabs = [VOCAB_DNA, VOCAB_RNA]
    uniques = vocabs.copy()
    all_uniques = []
    for i in range(len(uniques)):
        for j in range(j, len(uniques)):
            vci = set(uniques[i])
            vcj = set(uniques[j])
            uniques[i] = vci.difference(vcj)
            uniques[j] = vcj.difference(vci)
            all_uniques.extend(uniques[i])
            all_uniques.extend(uniques[j])
    
    for s in seq:
        if s in all_uniques:
            for vc, uni in zip(vocabs, uniques):
                if s in uni:
                    return vc
    return None
        
def get_random_sequence(length):
    seq = [random.choice(VOCAB) for i in range(length)]
    return "".join(seq)

def _convert_seq_vocab(seq, to_vocab, from_vocab = None):
    
    if not from_vocab:
        from_vocab = VOCAB
    
    for to_v, from_v in zip(to_vocab, from_vocab):
        seq = seq.replace(from_v, to_v)
    
    return seq

def _convert_list_vocab(vc_list, to_vocab, from_vocab = None):
    if not from_vocab:
        from_vocab = VOCAB
    new_list = []
    for s in vc_list:
        for to_v, from_v in zip(to_vocab, from_vocab):
            news = s.replace(from_v, to_v)
            new_list.append(news)
    return new_list

def _convert_dict_vocab(vc_dict, to_vocab, from_vocab = None, do_keys = True, do_values = True):
    if not do_keys and not do_values:
        return vc_dict
    if not from_vocab:
        from_vocab = VOCAB
        
    new_dict = {}
    
    for k,v in vc_dict.items():
        newk = _convert_seq_vocab(k, to_vocab, from_vocab = from_vocab)
        if do_values:
            newv = _convert_seq_vocab(v, to_vocab, from_vocab = from_vocab)
        else:
            newv = v
        new_dict[newk] = newv
    
    return new_dict

def _make_complement_map(vocab):
    cm = {}
    for i in range(0, len(vocab), 2):
        b1 = vocab[i]
        b2 = vocab[i+1]
        cm[b1]=b2
        cm[b2]=b1
        
