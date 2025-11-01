

import itertools

VOCAB_DNA = "ATGC"
VOCAB_RNA = "AUGC"

VOCAB = VOCAB_DNA

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