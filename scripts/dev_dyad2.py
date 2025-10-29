
import random
import string
import numpy as np
import matplotlib.pyplot as plt

from ggene.motifs import dyad
from ggene.motifs.dyad import Dyad, find_all_dyads, reverse_complement
from ggene.genomemanager import GenomeManager

bases = "AUGC"
dna_bases = "ATGC"

color_names = ["black","red","green","yellow","blue","magenta","cyan","white"]
color_ints = range(len(color_names))

class Color:
    def __init__(self, text, background, effect):
        self.text = self.get_color(text, background = False)
        self.background = self.get_color(background, background = True)
        self.effect = self.get_effect(effect)
    
    def set(self, spec, bright = False, background = False, effect = False):
        if background:
            self.set_background(spec, bright=bright)
        else:
            self.set_text(spec, bright=bright)
        
        if effect:
            self.set_effect(spec)
    
    def as_tuple(self):
        return (self.text, self.background, self.effect)
    
    def set_text(self, text_spec, bright=False):
        self.text = self.get_color(text_spec, bright=bright, background = False)
        
    def set_background(self, background_spec, bright=False):
        self.text = self.get_color(background_spec, bright=bright, background = True)
        
    def set_effect(self, effect_spec):
        self.effect = self.get_effect(effect_spec)
        
    @staticmethod
    def get_color(color_spec, bright=False, background = False):
        
        cc = 0
        if isinstance(color_spec, str) and color_spec:
            cc = color_ints[color_names.index(color_spec)]
        elif isinstance(color_spec, int):
            cc = color_spec
        
        if bright:
            cc += 8
        bgc = "38;5;"
        if background:
            bgc = "48;5;"
        
        return f'\x1b[{bgc}{cc}m'

    @staticmethod
    def get_effect(effect_spec):
        if effect_spec == "bold":
            return "\x1b[1m"
        elif effect_spec == "dim":
            return "\x1b[2m"
        elif effect_spec == "underline":
            return "\x1b[4m"
        elif effect_spec == "blink":
            return "\x1b[5m"
        elif effect_spec == "reverse":
            return "\x1b[7m"
        return ""
    
    def __str__(self):
        return str(self.effect) + str(self.background) + str(self.text)
    
    @classmethod
    def from_specs(cls, text_spec = "", text_bright = False, bg_spec = "", bg_bright = False, effect_spec = ""):
        out = cls("","","")
        out.set_text(text_spec, bright = text_bright)
        out.set_background(bg_spec, bright=bg_bright)
        out.set_effect(effect_spec)
        return out

RESET = '\033[0m'

CS = Color.from_specs(text_spec=250, text_bright = True, effect_spec ="")
CD = Color.from_specs(text_spec="yellow", effect_spec ="")
CL = Color.from_specs(text_spec="cyan",effect_spec ="")
CB = Color.from_specs(text_spec="blue",effect_spec ="")
CC = Color.from_specs(text_spec="cyan",effect_spec ="")

def set_colors(tail=None, dyad=None, loop=None, seq=None, cseq=None, bright=False, background = False, effect = None):
    
    global CS, CD, CL, CB, CC
    
    if tail:
        CS.set(tail, bright=bright, background = background, effect = effect)
        
    if dyad:
        CD.set(dyad, bright=bright, background = background, effect = effect)
        
    if loop:
        CL.set(loop, bright=bright, background = background, effect = effect)
        
    if seq:
        CB.set(seq, bright=bright, background = background, effect = effect)
        
    if cseq:
        CC.set(cseq, bright=bright, background = background, effect = effect)

set_colors(seq = 174, cseq = 66, background = True)

def find_nice_colors(seq, num_tries = 10, bright = False, background = True, effect = None):
    
    for n in range(num_tries):
        cb_col = random.randint(0, 255)
        cc_col = random.randint(0, 255)
        
        print(f"CB = {cb_col}, CC = {cc_col}")
        
        set_colors(seq = cb_col, cseq = cc_col, bright=bright, background=background, effect = effect)
        
        cseq = dyad.complement(seq)
        rcseq = dyad.reverse_complement(seq)
        rseq = dyad.reverse(seq)
        
        print_seq(seq)
        print_seq(cseq)
        print_seq(rcseq)
        print_seq(rseq)
        print()
        
        if n%5 == 4:
            input()
    

def print_seq(seq, color_complement = True):
    
    if not color_complement:
        print(seq)
        return
    
    noncomp = dyad.vocab[::2]
    comp = dyad.vocab[1::2]
    
    iscomp = -1
    colored_seq = []
    for s in seq:
        scomp = s in comp
        sym = s
        
        if scomp == iscomp:
            pass
        else:
            colored_seq.append(RESET)
            colored_seq.append(str(CC) if iscomp else str(CB))
            iscomp = scomp
        
        colored_seq.append(sym)
    
    colored_seq.pop(0)
    colored_seq.append(RESET)
    printseq = "".join(colored_seq)
    print(printseq)
    
def print_dyad(seq, dyad_len, dyad_start, loop_len, dyad_rc_start, total_len = -1, pre = "", post = ""):
    
    loop_start = dyad_start + dyad_len
    dyad_rc_start = loop_start + loop_len
    tail_start = dyad_rc_start + dyad_len
    len_seq = len(seq)
    dyad_seq_len = 2*dyad_len + loop_len
    
    parts = []
    
    head_start = 0
    tail_end = len_seq
    if total_len:
        crop = total_len - 2*dyad_len - loop_len
        head_start = dyad_start - crop//2
        tail_end = tail_start + (crop+1)//2
    
    if head_start < 0:
        tail_end += abs(head_start)
        head_start = 0
    if tail_end > len_seq:
        head_start -= (tail_end - len_seq)
        tail_end = len_seq
    
    cs = Color.get_color(0, bright=True)
    cd = Color.get_color("yellow")
    cl = Color.get_color("cyan")
    
    parts.append(pre)
    parts.append(f"{cs}{seq[head_start:dyad_start]}{RESET}")
    parts.append(f"{cd}{seq[dyad_start:loop_start]}{RESET}")
    parts.append(f"{cl}{seq[loop_start:dyad_rc_start]}{RESET}")
    parts.append(f"{cd}{seq[dyad_rc_start:tail_start]}{RESET}")
    parts.append(f"{cs}{seq[tail_start:tail_end]}{RESET}")
    parts.append(post)
    print("".join(parts))


def convert_vocab(seq, init_vocab, to_vocab):
    
    for v, vp in zip(init_vocab, to_vocab):
        seq = seq.replace(v, vp)
    
    return seq

def make_dyad(dyad_len, loop_len, total_seq_len):
    
    total_dyad_len = 2*dyad_len + loop_len
    
    dyad_start = random.randint(0, total_seq_len - total_dyad_len)
    
    vc = dyad.vocab
    vc_len = len(vc) - 1
    
    dyad = [vc[random.randint(0, vc_len)] for n in range(dyad_len)]
    loop = [vc[random.randint(0,vc_len)] for nl in range(loop_len)]
    full_dyad = dyad + loop + list(reverse_complement(dyad))
    s0 = [vc[random.randint(0, vc_len)] for ns in range(dyad_start)]
    s1 = [vc[random.randint(0,vc_len)] for ns in range(total_seq_len - total_dyad_len - dyad_start)]
    
    return "".join(s0 + full_dyad + s1), dyad_start

def make_random(total_seq_len):
    return "".join([dyad.vocab[random.randint(0,len(dyad.vocab) - 1)] for nl in range(total_seq_len)])

def get_hammerhead():
    
    seq = "YYRRGCCGUUACCURCAGCUGAUGAGCUCCAARAAGAGCGAAACCXRXYAGGUCCUGYAGUAYUGGCYXRXXXXXX"
    subs = {"Y":"C","R":"G","C":"C","G":"G","U":"U","A":"A"}
    out = []
    xbases = "ACGUACGUA"
    
    for s in seq:
        
        ssub = subs.get(s, s)
        if ssub == "X":
            # ssub = random.sample(bases, 1)[0]
            ssub = xbases[0]
            xbases = xbases[1:]
        
        out.append(ssub)
    return "".join(out)

def load_genome():
    return GenomeManager()

def get_linc_seq(geneman, vocab = None):
    
    chr = 2
    start = 207662375
    end = 207679116
    
    seq = geneman.get_sequence(chr, start, end)
    if vocab:
        for b,v in zip(dna_bases, vocab):
            seq = seq.replace(b, v)
    return seq



def print_nested_dyad(seq, dyad_map):
    
    newseq = [s for s in seq]
    
    done = False
    while not done:
        for n in range(len(newseq)):
            if newseq[n] in dyad_map:
                newseq.pop(n)
                newseq.insert(n)
    pass

def compare_results(seq, input, results):
    
    dyad_len, ds, loop_len, ds_rc = input
    print(ds, loop_len, ds+dyad_len+loop_len)
    
    print("ground truth:")
    print_dyad(seq, *input)
    
    print(f"{len(results)} possible dyads identified")
    for result in results:
        dl_res, ds_res, loop_len_res, ds_rc_res = result 
        
        print("search result:")
        print(ds_res, loop_len_res, ds_rc_res)
        print_dyad(seq, *result)
        print(f"result: {result == input}")
    

def test_find_all_dyads(seq_len, seq = None):
    
    min_stem = 3
    
    if not seq:
        seq = make_random(seq_len)
    
    alldyads = find_all_dyads(seq, min_stem, max_loop = 8)
    print(f"found {len(alldyads)} dyads")
    
    superdyads = dyad.build_superdyads(alldyads)
    if len(superdyads) > 0:
        max_len = max([s[0] for s in superdyads])
    else:
        max_len=0
    
    print(f"built {len(superdyads)} superdyads, max len {max_len}")
    
    for d in superdyads:
        print_dyad(seq, *d)
    
    return seq, alldyads, superdyads

def generate_alt_vocab(vocab_len):
    
    vocab_len = min(vocab_len, 26)
    if vocab_len %2 == 1:
        vocab_len -= 1
    
    vocab = "".join(random.sample(string.ascii_uppercase, vocab_len))
    
    cm = {}
    for vi in range(vocab_len//2):
        vv = vi + vocab_len // 2
        cm[vocab[vi]] = vocab[vv]
        cm[vocab[vv]] = vocab[vi]
        # okay.. what if complement wasn't cyclic :0
        # or subgroups of different order
    
    return vocab, cm

def generate_acyclic_vocab(vocab_len):
    
    vocab_len = min(vocab_len, 26)
    if vocab_len %2 == 1:
        vocab_len -= 1
    
    vocab = "".join(random.sample(string.ascii_uppercase, vocab_len))
    vc = vocab.copy()
    
    cm = {}
    v0 = vc.pop(0)
    v = v0
    while len(vc) > 0:
        vi = random.randint(0, len(vc) - 1)
        cm[v] = vc.pop(vi)
        v = cm[v]
    cm[v] = v0
    
    return vocab, cm

def test_compare_dyad_seqs(seq, dyads):
    
    dyad_seqs = dyad.extract_dyad_seqs(seq, dyads)
    
    same = []
    similar = []
    
    for i in range(len(dyad_seqs)):
        for j in range(i+1, len(dyad_seqs)):
            
            err = dyad.compare_seqs(dyad_seqs[i], dyad_seqs[j])
            
            if err < 1:
                same.append((i, j))
            elif err < 2:
                similar.append((i,j))
    
    return same, similar

def test_substitute(seq, dyads):
    
    for d in dyads:
        
        newseq, ds, newsym = dyad.try_substitute(seq, d)
        
        print(f"{len(ds)} new dyads found after substitution")
        
        if len(ds) < 1:
            try:
                dyad.clear_symbol(newsym)
                print(f"cleared symbol {newsym} and its complement")
            except:
                print(dyad.vocab, dyad.COMPLEMENT_MAP, dyad.aliases)
                raise

def test_linc(gm):
    
    vocab = bases
    linc_seq = get_linc_seq(gm, vocab)
    seq_len = len(linc_seq)
    
    dyads, supers = dyad.find_all_dyads_build(linc_seq, 3, max_loop = 8, min_loop = 3)
    for sup in supers[:10]:
        print_dyad(linc_seq, *sup, total_len = 30)
        
    # new_seq, shift = dyad.encode_dyad(linc_seq, supers[0])
    # print_dyad(linc_seq, *supers[0], total_len = 50)
    # print_dyad(new_seq, *supers[0], total_len = 50)
    
    same, similar = test_compare_dyad_seqs(linc_seq, supers)
    
    print(f"{len(same)} identical seqs, {len(similar)} similar seqs")
    
    test_substitute(linc_seq, supers)
    
    # dyad_lens = {}
    # for d in supers:
    #     dlen = d[0]
    #     n = dyad_lens.get(dlen, 0)
    #     dyad_lens[dlen] = n + 1
        
    # print(dyad_lens)
    
    # freqs, max_seq = dyad.frequency_rank_dyads(linc_seq, supers)
    # print(freqs)
    # print(max_seq, freqs[max_seq])
    
    # dyad_seqs = dyad.extract_dyad_seqs(linc_seq, supers)
    # dyad_seq_dd = list(set(dyad_seqs))
    # print(f"{len(supers)} original dyad seqs, {len(dyad_seq_dd)} after deduplication")

def load_intron(gm):
    
    gene_name = "ZNF804A"
    chr = 2
    pos = 184913701
    
    gene = gm.load_gene(gene_name)
    intrs = gene.get_feature("intron", sfid="intron-2")
    
    if len(intrs) > 0:
        intr = intrs[0]
    
    seq = gm.get_feature_sequence(intr, as_rna = True)
    refseq = gm.get_feature_sequence(intr, personal = False, as_rna = True)
    
    return seq, refseq
    
def test_intron(gm):
    
    fullseq, fullrefseq = load_intron(gm)
    if not dyad.vocab == bases:
        fullseq = convert_vocab(fullseq, bases, dyad.vocab)
        fullrefseq = convert_vocab(fullrefseq, bases, dyad.vocab)
    
    seq = fullseq[:1024]
    refseq = fullrefseq[:1024]
    
    cseq = "".join(dyad.complement(seq))
    rcseq = "".join(dyad.reverse_complement(seq))
    
    # allds = dyad.find_all_dyads(seq, 3, 4, 8, 3)
    # ds, sups = dyad.find_all_dyads_build(seq, 3, 8, 3)
    ds = dyad.find_all_dyads(seq, 3, max_stem = 4, max_loop = 8, min_loop = 3)
    sups = Dyad.build_superdyads(ds)
    sups2 = Dyad.build_superdyads2(ds)
    
    print(f"found {len(ds)} dyads, {len(sups)}")
    
    
    dmaxs = [d.maximize() for d in sups]
    
    for d in dmaxs:
        d.print(total_len = 64)
    
    
    return
    
    allds = list(sorted(allds, key = lambda a:a[1]))
    maxd = max([d[0] for d in allds])
    meand = np.mean([d[0] for d in allds])
    
    print(f"{len(allds)} dyads identified, mean {meand:0.1f} max {maxd}")
    
    dseqs = dyad.extract_dyad_seqs(seq, allds)
    
    tot = max([len(s) for s in dseqs]) + 8
    for d, _seq in zip(allds, dseqs):
        
        post = ""
        dexp, err = dyad.expand_dyad(seq, d)
        if err == 0:
            post = " (can expand)"
        
        print_dyad(seq, *d, total_len = tot, post = post)
    
    return seq

def test_hammerhead():
    hh = get_hammerhead()
    print(hh)
    
    dyad_map = {}
    
    dyads, supers = dyad.find_all_dyads_build(hh, 3, max_loop = 8, min_loop = 3)
    for sup in supers:
        print_dyad(hh, *sup, total_len = len(hh))
    
    da = supers[1]
    
    newseq, newds, newsym = dyad.try_substitute(hh, da)
    dyad_map[newsym] = da
    
    print(f"{len(newds)} dyads found after substitution")
    # if len(newds) > 0:
        
    print(newseq)
    
    dyads, supers = dyad.find_all_dyads_build(newseq, 3, max_loop = 14, min_loop = 3)
    for sup in supers:
        print_dyad(newseq, *sup, total_len = len(hh))
    
    db = supers[2]
    
    newseq, newds, newsym = dyad.try_substitute(newseq, db)
    dyad_map[newsym] = db
    
    print(f"{len(newds)} dyads found after substitution")
    print(newseq)
    
    dyads, supers = dyad.find_all_dyads_build(newseq, 3, max_loop = 20, min_loop = 3)
    for sup in supers:
        print_dyad(newseq, *sup, total_len = len(hh))
    
    dc = supers[0]
    
    newseq, newds, newsym = dyad.try_substitute(newseq, dc)
    dyad_map[newsym] = dc
    
    print(f"{len(newds)} dyads found after substitution")
    print(newseq)
    
    dyads, supers = dyad.find_all_dyads_build(newseq, 3, max_loop = 18, min_loop = 3)
    for sup in supers:
        print_dyad(newseq, *sup, total_len = len(hh))
    
    d0 = supers[1]
    print_dyad(newseq, *d0, total_len = 20)
    dcon, err = dyad.contract_dyad(newseq, d0)
    print(f"error after contraction: {err}")
    print_dyad(newseq, *dcon, total_len = 20)
    
    newseq, newds, newsym = dyad.try_substitute(newseq, dcon)
    dyad_map[newsym] = dcon
    
    print(newseq)
    
    
    
def try_expand_contract():
    pass
    # hp1 = supers[1]
    # hp2 = supers[-1]
    # d0 = supers[1]
    
    # print_dyad(hh, *d0, total_len = 20)
    
    # dcon, err = dyad.contract_dyad(hh, d0)
    # print(f"error after contraction: {err}")
    # print_dyad(hh, *dcon, total_len = 20)
    
    # dexp, err = dyad.expand_dyad(hh, d0)
    # print(f"error after expansion: {err}")
    # print_dyad(hh, *dexp, total_len = 20)
    
    

def main():
    
    gm = load_genome()
    
    # vocab = dyad._comp_strs[:2] + dyad._comp_strs[4:6]
    # cm = {vocab[i]:vocab[i+1-2*(i%2)] for i in range(len(vocab))}
    # dyad.set_vocab(vocab, cm)
    print(dyad.vocab, dyad.COMPLEMENT_MAP, dyad.aliases)
    
    seq = test_intron(gm)
    
    # find_nice_colors(subseq, num_tries = 50)
    
    

if __name__=="__main__":
    main()
