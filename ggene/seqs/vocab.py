
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
Ω ℧
⊓ and ⊔ or ∏ ∐
 	⊻ 	⊼
Ո and Ս ⋀ ⋁ ⋂ ⋃ ⋎ ⋏
Ψ ψ Ϙ Ͳ ⋔
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

okay how about
⊑ and ⊏, ⊆ and ⊂


quads
_   0000    ----    -
▖   0001    ---A    A
▗   0010    --T-    T
▘   0100    -G--    G
▝   1000    C---    C

▄   0011    --TA    W
▌   0101    -G-A    R
▐   1010    C-T-    Y
▀   1100    CG--    S
▞   1001    C--A    M
▚   0110    -GT-    K

▙   0111    -GTA    D
▛   1101    CG-A    V
▜   1110    CGT-    B
▟   1011    C-TA    H

█   1111    CGTA    N
    



wow! these are unified canadian aboriginal syllabics:
ᐁ ᐂ ᐃ ᐄ 
ᑌ ᑍ ᑎ 
ᕫ ᕬ
ᗊ ᗋ ᗄ ᗅ
ᗯ ᗰ 
ᗨ ᗩ
ᗵ ᗶ 
ᗻ ᗼ 
ᗐ ᗑ
ᗜ ᗝ
ᙀ ᙁ
ᕰ ᕱ 
ᗢ ᗣ
(lots of these are wider than they appear)
ᙎ ᙏ ᙔ ᙕ 
ᗀ ᗃ ᗁ ᗂ ᖸ ᖹ ᖺ ᖻ ᘈ ᘉ ᖼ ᖽ ᖾ ᖿ ᖊ ᖋ ᖌ ᖍ ᖠ ᖢ ᖤ ᘂ ᘃ
ᕓ ᕕ 
ᗖ ᗗ
ᑫ ᑭ ᑯ ᑲ ᖰ ᖱ ᖗ ᖘ ᖙ ᖚ ᘤ ᘧ 
ᙈ ᙉ 
ᙡ ᙢ ᙧ ᙨ 
ᘎ ᘏ ᘐ ᘓ ᘔ ᘕ ᘖ ᘜ ᘝ


ᐁ ᐂ ᐃ ᐄ ᐌ ᐍ ᐎ ᐏ ᐐ ᐑ 
ᐫ ᐬ ᐯ ᐰ ᐱ ᐲ ᐺ ᐻ ᐼ ᐽ ᐾ ᐿ
ᑌ ᑍ ᑎ ᑏ ᑐ ᑕ ᑗ ᑘ ᑙ ᑚ ᑛ ᑜ ᑝ ᑞ ᑟ 
ᑠ ᑡ ᑢ ᑣ ᑤ ᑥ ᑧ ᑨ ᑩ ᑪ ᑫ ᑬ ᑭ ᑮ ᑯ ᑰ ᑱ ᑲ ᑳ ᑴ ᑵ ᑶ ᑷ ᑸ ᑹ ᑺ ᑻ ᑼ ᑽ ᑾ ᑿ
ᒀ ᒁ ᒂ ᒅ ᒆ ᒇ ᒈ ᒉ ᒊ ᒋ ᒌ ᒍ ᒎ ᒏ ᒐ ᒑ ᒒ ᒓ ᒔ ᒕ ᒖ ᒗ ᒘ ᒙ ᒚ ᒛ ᒜ ᒝ ᒞ  ᒠ 
ᒣ ᒤ ᒥ ᒦ ᒧ ᒨ ᒩ ᒪ ᒫ ᒬ ᒭ ᒮ ᒯ ᒰ ᒱ ᒲ ᒳ ᒴ ᒵ ᒶ ᒷ ᒸ ᒹ ᒺ 
ᓀ ᓁ ᓂ ᓃ ᓄ ᓅ ᓆ ᓇ ᓈ ᓉ ᓊ ᓋ ᓌ ᓍ ᓎ ᓏ ᓓ ᓔ ᓕ ᓖ ᓗ ᓘ ᓙ ᓚ ᓛ ᓜ ᓝ ᓞ ᓟ
ᓠ ᓡ ᓢ ᓣ ᓤ ᓥ ᓦ ᓧ ᓨ ᓩ ᓬ ᓭ ᓮ ᓯ ᓰ ᓱ ᓲ ᓳ ᓴ ᓵ ᓶ ᓷ ᓸ ᓹ ᓺ ᓻ ᓼ ᓽ ᓾ ᓿ
ᔀ ᔁ ᔂ ᔃ ᔄ ᔐ ᔑ ᔒ ᔓ ᔔ ᔕ ᔖ ᔗ ᔘ ᔙ ᔚ ᔛ ᔜ ᔝ ᔞ ᔟ
ᔠ ᔡ ᔢ ᔣ ᔤ ᔦ ᔧ ᔨ ᔩ ᔪ ᔫ ᔬ ᔭ ᔮ ᔯ ᔰ ᔱ ᔲ ᔳ ᔴ ᔵ ᔶ ᔷ ᔸ ᔹ ᔺ ᔻ ᔼ ᔽ
ᕂ ᕃ ᕄ ᕅ ᕆ ᕇ ᕈ ᕉ ᕊ ᕋ ᕌ ᕍ ᕎ ᕏ ᕒ ᕓ ᕔ ᕕ ᕖ ᕗ ᕘ ᕙ ᕚ ᕛ ᕜ ᕞ ᕟ
ᕠ ᕡ ᕢ ᕣ ᕤ ᕥ ᕦ ᕧ ᕨ ᕩ ᕫ ᕬ ᕭ ᕮ ᕯ ᕰ ᕱ ᕲ ᕳ ᕴ ᕵ ᕶ ᕷ ᕸ ᕹ ᕺ ᕼ ᕾ ᕿ
ᖀ ᖁ ᖂ ᖃ ᖄ ᖅ ᖆ ᖇ ᖈ ᖉ ᖊ ᖋ ᖌ ᖍ ᖎ ᖏ ᖐ ᖑ ᖒ ᖓ ᖔ ᖕ ᖖ ᖗ ᖘ ᖙ ᖚ ᖛ ᖜ ᖝ ᖞ
ᖠ ᖡ ᖢ ᖣ ᖤ ᖥ ᖦ ᖧ ᖨ ᖩ ᖪ ᖫ ᖬ ᖭ ᖯ ᖰ ᖱ ᖲ ᖳ ᖴ ᖵ ᖶ ᖷ ᖸ ᖹ ᖺ ᖻ ᖼ ᖽ ᖾ ᖿ
ᗀ ᗁ ᗂ ᗃ ᗄ ᗅ ᗆ ᗇ ᗈ ᗉ ᗊ ᗋ ᗌ ᗍ ᗎ ᗏ ᗐ ᗑ ᗒ ᗓ ᗔ ᗕ ᗖ ᗗ ᗘ ᗙ ᗚ ᗛ ᗜ ᗝ ᗞ ᗟ
ᗠ ᗡ ᗢ ᗣ ᗤ ᗥ ᗦ ᗧ ᗨ ᗩ ᗪ ᗫ ᗬ ᗭ ᗯ ᗰ ᗱ ᗲ ᗳ ᗴ ᗵ ᗶ ᗷ ᗸ ᗹ ᗺ ᗻ ᗼ ᗽ ᗾ ᗿ
ᘀ ᘂ ᘃ ᘄ ᘅ ᘆ ᘇ ᘈ ᘉ ᘊ ᘋ ᘌ ᘍ ᘎ ᘏ ᘐ ᘑ ᘒ ᘓ ᘔ ᘕ ᘖ ᘗ ᘘ ᘙ ᘚ ᘛ ᘜ ᘝ ᘞ ᘟ
ᘠ ᘡ ᘢ ᘣ ᘤ ᘥ ᘦ ᘧ ᘨ ᘩ ᘪ ᘫ ᘬ ᘭ ᘮ ᘯ ᘰ ᘱ ᘲ ᘳ ᘴ ᘵ ᘶ ᘷ ᘸ ᘹ ᘺ ᘻ ᘼ ᘽ ᘾ ᘿ
ᙀ ᙁ ᙂ ᙃ ᙄ ᙅ ᙈ ᙉ ᙊ ᙋ ᙌ ᙍ ᙎ ᙏ ᙐ ᙑ ᙒ ᙓ ᙔ ᙕ ᙖ ᙗ ᙘ ᙙ ᙚ ᙛ ᙜ ᙝ ᙞ ᙟ
ᙠ ᙡ ᙢ ᙣ ᙤ ᙥ ᙦ ᙧ ᙨ ᙩ ᙪ ᙫ ᙬ ᙭
ᢰ ᢱ ᢲ ᢳ ᢴ ᢵ ᢶ ᢷ ᢸ ᢹ ᢺ ᢻ ᢼ ᢽ ᢾ ᢿ ᣀ ᣁ ᣂ ᣃ ᣄ ᣅ ᣆ ᣇ ᣈ ᣉ ᣊ ᣋ ᣌ ᣍ ᣎ ᣏ
ᣐ ᣑ ᣒ ᣓ ᣠ ᣡ ᣢ ᣣ ᣤ ᣥ ᣦ ᣧ ᣨ ᣩ ᣪ ᣫ ᣬ ᣭ ᣮ ᣯ
theyre like perfect, used for writing Inuktitut, Carrier, Cree + dialects, Ojibwe, Blackfoot, Canadian Athabascan. additions for more Cree dialects, Ojibwe, Dene in extended.

"""

_latin_rotated = {
    "c":"ᴒ",
    "":"",
    "":"",
    "I":"ꟷ",
    "":"",
    "m":"ᴟ",
    "n":"ᴝ",
    "o":"ᴑ",
    "T":"Ꟶ",
    
}

_addl_syms = "ЂЄЉЊЋБГ"

quad_vc_base = "-ATGCWRYSAKDVBHN"
quad_vc_block = "_▖▗▘▝▄▌▐▀▞▚▙▛▜▟█"

quad_vc = {bs:bk for bs, bk in zip(quad_vc_base, quad_vc_block)}

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
        
