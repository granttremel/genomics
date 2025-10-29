
import random
import string
import numpy as np
import matplotlib.pyplot as plt

from ggene.motifs import dyad
from ggene.motifs.dyad import search_dyad, seek_dyad, seek_dyad_rev, find_all_dyads, reverse_complement

bases = "AUCG"

RESET = '\033[0m'
CS = '\033[90m' # gray
CD = '\033[93m' # Yellow
CL = '\033[96m' # cyan

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
    
def _print_dyad(seq, dyad_len, dyad_start, loop_len, dyad_rc_start):
    
    loop_start = dyad_start + dyad_len
    dyad_rc_start = loop_start + loop_len
    tail_start = dyad_rc_start + dyad_len
    len_seq = len(seq)
    
    seq_dyads = " "*dyad_start + seq[dyad_start:loop_start] + " "*loop_len + seq[dyad_rc_start:tail_start] + " "*(len_seq - tail_start)
    seq_loop = " "*(dyad_start + dyad_len) + seq[loop_start:dyad_rc_start] + " "*(len_seq - dyad_rc_start)
    seq_else = seq[:dyad_start] + " "*(tail_start - dyad_start) + seq[tail_start:]
    
    print(seq)
    print(seq_dyads)
    print(seq_loop)
    print(seq_else)

def print_dyad(seq, dyad_len, dyad_start, loop_len, dyad_rc_start):
    
    loop_start = dyad_start + dyad_len
    dyad_rc_start = loop_start + loop_len
    tail_start = dyad_rc_start + dyad_len
    len_seq = len(seq)
    
    parts = []
    parts.append(f"{CS}{seq[:dyad_start]}{RESET}")
    parts.append(f"{CD}{seq[dyad_start:loop_start]}{RESET}")
    parts.append(f"{CL}{seq[loop_start:dyad_rc_start]}{RESET}")
    parts.append(f"{CD}{seq[dyad_rc_start:tail_start]}{RESET}")
    parts.append(f"{CS}{seq[tail_start:]}{RESET}")
    print("".join(parts))

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
    

def test_make_dyad():
    
    dyad_len = 10
    loop_len = 4
    total_seq_len = 30
    
    d, ds = make_dyad(dyad_len, loop_len, total_seq_len)
    
    print_dyad(d, dyad_len, ds, loop_len)
    

def test_search_dyad():
    d, input, limits = setup_test()
    dyad_len, _, _, _ = input
    _, _, max_loop, min_loop = limits
    
    res = search_dyad(d, dyad_len, max_loop, min_loop)
    result = (dyad_len, res[0], res[1], res[2])
    
    compare_results(d, input, [result])
    
def test_search_dyad_inout():
    d, input, limits = setup_test()
    dyad_len, _, _, _ = input
    _, _, max_loop, min_loop = limits
    
    res = search_dyad(d, dyad_len, max_loop, min_loop)
    result = (dyad_len, res[0], res[1], res[2])
    
    compare_results(d, input, [result])

def test_seek_dyad():
    d, input, limits = setup_test()
    
    allres = seek_dyad(d, *limits)
    
    results = [(v[0], k[0], v[1], k[1]) for k,v in allres.items()]
    
    compare_results(d, input, results)
    
def test_seek_dyad_rev():

    d, input, limits = setup_test()
    
    res = seek_dyad_rev(d, *limits)

    compare_results(d, input, [res])

def test_find_all_dyads(seq_len, seq = None):
    
    min_dyad = 3
    
    # seq = "GUUAAAAUGGGCAAUUUAAG"
    if not seq:
        seq = make_random(seq_len)
    
    # print(f"test sequence: {seq}")
    
    alldyads = find_all_dyads(seq, min_dyad)
    print(f"found {len(alldyads)} dyads")
    
    superdyads = dyad.build_superdyads(alldyads)
    max_len = max([s[0] for s in superdyads])
    print(f"built {len(superdyads)} superdyads, max len {max_len}")
    
    for d in superdyads:
        print_dyad(seq, *d)
    
    return seq, alldyads, superdyads

def test_find_all_build(seq_len, seq = None):
    
    min_dyad = 3
    if not seq:
        seq = make_random(seq_len)
    
    dyads, superdyads = dyad.find_all_dyads_build(seq, min_dyad)
    
    max_len = max([s[0] for s in superdyads])
    
    print(f"found {len(dyads)} naive dyads and {len(superdyads)} superdyads, max len {max_len}")
    
    for d in superdyads:
        print_dyad(seq, *d)
        
    return seq, dyad, superdyads

def test_find_all_chunk(seq_len, chunksz, seq = None):
    
    min_dyad = 3
    if not seq:
        seq = make_random(seq_len)
    
    superdyads = dyad.find_all_dyads_chunk(seq, min_dyad, chunksz)
    
    max_len = max([s[0] for s in superdyads])
    
    print(f"found {len(superdyads)} superdyads, max len {max_len}")
    
    for d in superdyads:
        print_dyad(seq, *d)
        
    return seq, dyad, superdyads

def compare_find_alls(seq_len):
    
    seq = make_random(seq_len)
    
    _, alld, superd = test_find_all_dyads(seq_len, seq = seq)
    
    _, alldb, superdb = test_find_all_build(seq_len, seq = seq)
    
    pass


def test_dyad_stats(seq_len = 30, num_tests = 10):
    
    ninit = []
    init_lens = []
    nsuper = []
    super_lens = []
    diff = []
    
    max_init = max_super = (-1,-1,-1,-1)
    max_init_seq = max_super_seq = ""
    
    for n in range(num_tests):
        
        seq, alls, supers = test_find_all_dyads(seq_len)
        
        ninit.append(len(alls))
        init_lens.extend([d[0] for d in alls])
        
        nsuper.append(len(supers))
        super_lens.extend([d[0] for d in supers])
        
        diff.append(len(alls) - len(supers))
        
        for a in alls:
            if a[0] > max_init[0]:
                max_init = a
                max_init_seq = seq
        for s in supers:
            if s[0] > max_super[0]:
                max_super = s
                max_super_seq = seq
        
    ninitmean = np.mean(ninit)
    ninitsd = np.std(ninit)
    initlenmean = np.mean(init_lens)
    initlensd = np.std(init_lens)
    nsupermean = np.mean(nsuper)
    nsupersd = np.std(nsuper)
    superlenmean = np.mean(super_lens)
    superlensd = np.std(super_lens)
    diffmean = np.mean(diff)
    diffsd = np.std(diff)
    
    initlenmax = np.max(init_lens)
    superlenmax = np.max(super_lens)
    
    print(f"Results with seq_len {seq_len} and num_tests {num_tests}")
    print(f"Naive dyads:")
    print(f"Number: {ninitmean:0.1f} ({ninitsd:0.3f})")
    print(f"Length: {initlenmean:0.1f} ({initlensd:0.3f})")
    print(f"Max length: {initlenmax}")
    print_dyad(max_init_seq, *max_init)
    
    print(f"Unique superdyads:")
    print(f"Number: {nsupermean:0.1f} ({nsupersd:0.3f})")
    print(f"Length: {superlenmean:0.1f} ({superlensd:0.3f})")
    print(f"Max length: {superlenmax}")
    print_dyad(max_super_seq, *max_super)
    
    print(f"Change in #: {diffmean:0.1f} ({diffsd:0.3f})")



def test_unique_superdyads(seq, dyads):
    
    print(f"{len(dyads)} dyads to start")
    for d in dyads:
        print_dyad(seq, *d)
    
    supers = dyad.build_superdyads(dyads)
    
    print(f"{len(supers)} unique superdyads found")
    for d in supers:
        print_dyad(seq, *d)
    return supers
    

def test_check_subdyads(seq, dyads):
    
    for i in range(len(dyads)):
        da = dyads[i]
        print(f"dyad {i}")
        print_dyad(seq, *da)
        # print()
        for j in range(i+1, len(dyads)):
            db = dyads[j]
            res = dyad.check_subdyad(da, db)
            res_fast = dyad.check_subdyad_fast(da, db)
            
            print(f"dyad {j}")
            print_dyad(seq, *db)
            print(f"{j} is subdyad of {i}: {res} ({res_fast})")
            
            res_mut, dc = dyad.check_mutual_subdyad(da, db)
            if res_mut:
                print("superdyad:")
                print_dyad(seq, *dc)
            else:
                print("no superdyad..")
    

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

def test_alt_vocab():
    
    vc_sz = 12
    
    vc, cm = generate_alt_vocab(vc_sz)
    print(vc, cm)
    
    vc, cm = generate_acyclic_vocab(vc_sz)
    print(vc, cm)

def test_set_vocab():
    
    print(dyad.vocab, dyad.COMPLEMENT_MAP)
    new_vc, new_cm = generate_alt_vocab(4)
    
    dyad.set_vocab(new_vc, new_cm)
    print(dyad.vocab, dyad.COMPLEMENT_MAP)
    
    dyad.reset_vocab()
    print(dyad.vocab, dyad.COMPLEMENT_MAP)
    
def setup_test():
    max_dyad = 10
    min_dyad = 3
    max_loop = 8
    min_loop = 4
    total_seq_len = 30
    
    dyad_len = random.randint(min_dyad, max_dyad)
    loop_len = random.randint(min_loop, max_loop)
    d, ds = make_dyad(dyad_len, loop_len, total_seq_len)
    return d, (dyad_len, ds, loop_len, dyad_len + ds + loop_len), (max_dyad, min_dyad, max_loop, min_loop)

def main():
    
    # test_search_dyad()
    # test_seek_dyad()
    # test_seek_dyad_rev()
    # test_find_all_dyads(30)
    # test_find_all_build(30)
    
    # vc, cm = generate_alt_vocab(4)
    # print(f"vocab: {vc}")
    # dyad.set_vocab(vc, cm)
    
    seq_len = 200
    seq = make_random(seq_len)
    test_find_all_build(seq_len, seq = seq)
    test_find_all_chunk(seq_len, seq_len // 5, seq = seq)
    # compare_find_alls(200)
    # test_dyad_stats(seq_len=50, num_tests = 30)
    
    
    # test_alt_vocab()
    # test_set_vocab()
    
    pass

if __name__=="__main__":
    main()
