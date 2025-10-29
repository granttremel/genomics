

import numpy as np
import string

_comp_strs = "AĀEĒGḠIĪOŌUŪYȲ"
vocab = "AUGC"
_vc = vocab

COMPLEMENT_MAP = {'A':'U','C':'G','G':'C','U':'A','N':'N'}
_cm = COMPLEMENT_MAP

aliases = {}

def set_vocab(new_vocab, new_cm):
    global vocab, COMPLEMENT_MAP
    
    vocab = new_vocab
    COMPLEMENT_MAP = new_cm

def reset_vocab():
    global vocab, COMPLEMENT_MAP
    vocab = _vc
    COMPLEMENT_MAP = _cm

def get_next_symbols(num_symbols):
    
    global vocab, COMPLEMENT_MAP
    if num_symbols %2 == 1:
        num_symbols -= 1
    
    syms = string.ascii_uppercase
    out =[]
    s = ''
    for sym in syms:
        if sym in vocab:
            continue
        out.append(sym)
        if s:
            COMPLEMENT_MAP[s] = sym
            COMPLEMENT_MAP[sym] = s
            
        else:
            s = sym
        
        if len(out) == num_symbols:
            break
    vocab = vocab + "".join(out)
    
    return out

def encode_symbols(seq):
    global aliases
    
    nsyms = 2*len(seq)
    seq = seq
    cs = "".join(reverse_complement(seq))
    sym, csym = get_next_symbols(2)
    
    aliases[sym] = seq
    aliases[seq] = sym
    aliases[csym] = cs
    aliases[cs] = csym
    return sym, seq, csym, cs

def encode_dyad(seq, dyad):
    
    dyad_seq = extract_dyad_seqs(seq, [dyad])[0]
    
    sym, subseq, csym, csubseq = encode_symbols(dyad_seq)
    
    newseq, shift = substitute_symbols(seq, sym, subseq, csym, csubseq)
    
    return newseq, sym

def substitute_symbols(seq, sym, subseq, csym, csubseq):
    
    len_seq = len(seq)
    len_ss = len(subseq)
    out = []
    shift = 0
    i = 0
    while i < len_seq:
        testseq = seq[i:i+len_ss]
        if testseq == subseq:
            out.append(sym)
            shift += len_ss - 1
            i += len_ss - 1
            # print(f"subbed symbol {sym} at {i}")
        elif testseq == csubseq:
            out.append(csym)
            shift += len_ss - 1
            i += len_ss - 1
            # print(f"subbed comp symbol {csym} at {i}")
        else:
            out.append(seq[i])
        i+=1
            
    return "".join(out), shift

def clear_symbol(sym):
    global vocab, COMPLEMENT_MAP, aliases
    
    csym = COMPLEMENT_MAP[sym]
    
    vocab = vocab.replace(sym, "")
    vocab = vocab.replace(csym, "")
    del COMPLEMENT_MAP[sym]
    del COMPLEMENT_MAP[csym]
    s = aliases.get(sym)
    cs = aliases.get(csym)
    del aliases[sym]
    del aliases[csym]
    try:
        del aliases[s]
    except:
        pass
    try:
        del aliases[cs]
    except:
        pass

def to_pyrpur(seq):
    
    
    
    pass

def shift_dyad(dyad, shift):
    return (dyad[0], dyad[1]+shift, dyad[2], dyad[3]+shift)

def compare_seqs(s1, s2, err_tol = None):
    
    if not len(s1) == len(s2):
        return len(s1)
    
    if err_tol is None:
        err_tol = len(s1)
    
    nerr = 0
    for a, b in zip(s1, s2):
        if not a == b:
            nerr += 1
        if nerr > err_tol:
            return nerr
    return nerr

def test_dyad(seq, dyad, err_tol = None):
    
    s1 = seq[dyad[1]:dyad[1] + dyad[0]]
    cs1 = reverse_complement(s1)
    s2 = seq[dyad[3]:dyad[3]+dyad[0]]
    
    err = compare_seqs(cs1, s2, err_tol = err_tol)
    return err

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

def extract_dyad_seqs(seq, dyads):
    return [seq[d[1]:2*d[0]+d[1] + d[2]] for d in dyads]

def extract_dyad_seq(seq, dyad):
    return seq[dyad[1]:2*dyad[0]+dyad[1] + dyad[2]]

def frequency_rank_dyads(seq, dyads):
    
    dseqs = extract_dyad_seqs(seq, dyads)
    
    max_n = 0
    max_seq = ""
    
    freqs = {}
    for ds in dseqs:
        n = freqs.get(ds, 0) + 1
        freqs[ds] = n
        if n > max_n:
            max_n = n
            max_seq = ds
    
    return freqs, max_seq

def search_dyad(seq, dyad_len, max_loop, min_loop, err_tol = 0):
    """
    dumb algorithm to find palindromes/dyads given their length
    a dyad should be unique given a sequence.. right?
    """
    
    dyad_start = dyad_rc_start = loop_len = -1
    n_loop_iters = max_loop - min_loop
    n_dyad_iters = len(seq) - 2*dyad_len - min_loop
    
    done = False
    allres = []
    
    for nd in range(n_dyad_iters):
        
        test_seq = seq[nd:nd + dyad_len]
        test_seq_rc = "".join(reverse_complement(test_seq))
        
        for nl in range(min_loop, max_loop+1):
        
            comp_seq = seq[nd + dyad_len + nl:nd + 2*dyad_len + nl]
            err = compare_seqs(test_seq_rc, comp_seq, err_tol = 0)
            
            # if comp_seq == test_seq_rc:
            if err < 1:
                done = True
                
            if done:
                dyad_start = nd
                dyad_rc_start = nd + dyad_len + nl
                loop_len = dyad_rc_start - dyad_start - dyad_len
                break
        
        if done:
            break
    
    return dyad_len, dyad_start, loop_len, dyad_rc_start

def search_dyad_inout(seq, dyad_len, max_loop, min_loop):
    """
    dumb algorithm to find palindromes/dyads given their length
    a dyad should be unique given a sequence.. right?
    """
    
    dyad_start = dyad_rc_start = loop_len = -1
    n_loop_iters = max_loop - min_loop
    n_dyad_iters = len(seq) - 2*dyad_len - min_loop
    
    done = False
    allres = []
    
    for nl in range(min_loop, max_loop+1):
        
        for nd in range(n_dyad_iters):
            test_seq = seq[nd:nd + dyad_len]
            test_seq_rc = "".join(reverse_complement(test_seq))
            comp_seq = seq[nd + dyad_len + nl:nd + 2*dyad_len + nl]
            
            if comp_seq == test_seq_rc:
                dyad_start = nd
                dyad_rc_start = nd + dyad_len + nl
                loop_len = dyad_rc_start - dyad_start - dyad_len
                allres.append((dyad_start, loop_len, dyad_rc_start))
                done = True
        
        if done:
            break
    
    return allres
    
def seek_dyad(seq, max_dyad, min_dyad, max_loop, min_loop):
    """
    really dumb algorithm to find palindromes/dyads
    """
    allres= {}
    
    for dl in range(min_dyad, max_dyad+1):
        res = search_dyad(seq, dl, max_loop, min_loop)
        if res[0] > -1:
            dyad_len, dyad_start, loop_len, dyad_rc_start = res
            k = (dyad_start, dyad_rc_start)
            if k in allres:
                _dl, _loop_len = allres[k]
                if dl > _dl:
                    allres[k] = (dl, loop_len)
            else:
                allres[k] = (dl, loop_len)
    
    return allres
    

def seek_dyad_rev(seq, max_dyad, min_dyad, max_loop, min_loop):
    """
    start big and early exit?
    """
    for dl in range(min_dyad, max_dyad+1, -1):
        res = search_dyad(seq, dl, max_loop, min_loop)
        if res[0] > -1:
            return dl, *res
    return -1, -1, -1, -1


def find_all_dyads(seq, min_dyad, max_dyad = -1, max_loop = -1, min_loop = -1, err_tol = None):
    
    len_seq = len(seq)
    if max_dyad < 0:
        max_dyad = len_seq // 2
    
    if max_loop < 0:
        max_loop = len_seq - 2*min_dyad
    if min_loop < 0:
        min_loop = 3
    
    ncomps = 0
    
    allres = []
    
    for dl in range(min_dyad, max_dyad):
        
        for ds in range(0, len_seq - 2*dl):
            
            test_seq = seq[ds:ds+dl]
            test_seq_rc = "".join(reverse_complement(test_seq))
            
            for dds in range(ds + dl + min_loop, ds + dl + max_loop + 1):
                
                comp_seq = seq[dds:dds+dl]
                err = compare_seqs(test_seq_rc, comp_seq, err_tol = err_tol)
                if err < 1:
                    loop_len = dds - dl - ds
                    allres.append((dl, ds, loop_len, dds))
                ncomps += 1
    
    print(f"tested {ncomps} sequences")
    
    return allres

def find_all_dyads_window(seq, windowsz, min_dyad, max_dyad = -1, max_loop = -1, min_loop = -1, merge = False, err_tol = None):
    
    nchunks = len(seq) // windowsz
    
    allds = []
    
    for n in range(nchunks):
        start = n*windowsz
        subseq = seq[start:start+windowsz]
        subds = find_all_dyads(subseq, min_dyad, max_dyad=max_dyad, max_loop=max_loop, min_loop = min_loop, err_tol = err_tol)
        
        subds = [shift_dyad(d, start) for d in subds]
        
        allds.extend(subds)
        
        if merge:
            allds = build_superdyads(allds)
        
    return allds

def find_all_dyads_expand(seq, min_dyad, min_loop = 3, max_loop = 8):
    
    len_seq = len(seq)
    max_dyad = min_dyad + 1
    
    if max_loop < 0:
        max_loop = len_seq - 2*min_dyad
    if min_loop < 0:
        min_loop = 3
    
    min_semiloop = (min_loop+1)//2
    max_semiloop = max_loop // 2
    
    ncomps = 0
    
    allres = []
    
    for dc in range(min_dyad, len_seq - min_loop):
        _max_dyad = min(max_dyad, dc - max_semiloop)
        
        for dl in range(min_dyad, _max_dyad):
            
            for dsl in range(min_semiloop, max_semiloop):
                d = (dl, dc - dsl - dl, 2*dsl, dc + dsl)
                test_seq = extract_dyad_seqs()
        # for dl in range(min_dyad, max_dyad):
            
        #     for ds in range(0, len_seq - 2*dl):
                
        #         test_seq = seq[ds:ds+dl]
        #         test_seq_rc = "".join(reverse_complement(test_seq))
                
        #         for dds in range(ds + dl + min_loop, ds + dl + max_loop + 1):
                    
        #             comp_seq = seq[dds:dds+dl]
        #             err = compare_seqs(test_seq_rc, comp_seq)
        #             if err < 1:
        #                 loop_len = dds - dl - ds
                        
        #                 d_new = (dl, ds, loop_len, dds)
                        
        #                 d_new = maximize_dyad(seq, d_new) # produces maximal dyad given a center position
                        
        #                 allres.append(d_new)
        #             ncomps += 1
    
    print(f"tested {ncomps} sequences")
    
    return allres
    
    pass

def maximize_dyads(seq, dyads, nmax = None):
    
    dyads = list(sorted(dyads, key = lambda d:d[0])) # sort by len
    out = []
    for d in dyads:
        dmax = maximize_dyad(seq, d, nmax = nmax)
        out.append(dmax)
    
    out = list(sorted(out, key = lambda d:d[1]))
    return out

def maximize_dyad(seq, d, nmax = None):
    
    dinit = d
    
    if not nmax:
        nmax_ctr = d[2] // 2
        nmax_exp = min(d[1], len(seq) - d[3])
    else:
        nmax_ctr = nmax
        nmax_exp = nmax
    
    n = 0
    while n < nmax_ctr:
        dctr, err = contract_dyad(seq, dinit)
        if err == 0:
            dinit = dctr
        else:
            break
        n+=1
    
    if n > 0:
        print(f"contracted dyad by {n}")
    
    n=0
    while n < nmax_exp:
        dexp, err = expand_dyad(seq, dinit)
        if err == 0:
            dinit = dexp
        else:
            break
        n+=1
    if n > 0:
        print(f"expanded dyad by {n}")
        
    return dinit

def contract_dyad(seq, dyad):
    nd = (dyad[0] + 1, dyad[1], dyad[2] - 2, dyad[3] - 1)
    err = test_dyad(seq, nd)
    return nd, err

def expand_dyad(seq, dyad):
    nd = (dyad[0] + 1, dyad[1] - 1, dyad[2], dyad[3])
    err = test_dyad(seq, nd)
    return nd, err

def try_substitute(seq, dyad, search_window = 20):
    
    dyad_seq = seq[dyad[1]:dyad[0]+dyad[1]]
    newseq, newsym = encode_dyad(seq, dyad)
    
    subseq = newseq[dyad[1] - search_window: dyad[1] + search_window]
    
    ds = find_all_dyads(subseq, 3, 4, 8, 4)
    
    return newseq, ds, newsym


def find_all_dyads_build(seq, min_dyad, max_loop = -1, min_loop = -1, err_tol = None):
    
    max_dyad = min_dyad + 1
    
    dyads = find_all_dyads(seq, min_dyad, max_dyad = max_dyad, max_loop = max_loop, min_loop = min_loop, err_tol = err_tol)
    
    supers = build_superdyads(dyads)
    
    return dyads, supers

def find_all_dyads_chunk(seq, min_dyad, chunksz):
    
    seq_len = len(seq)
    nchunks = seq_len // chunksz
    all_supers = []
    
    for nc in range(2*nchunks - 1):
        start = nc * chunksz // 2
        subseq = seq[start:start+chunksz]
        
        ds, ss = find_all_dyads_build(subseq, min_dyad)
        
        print(f"chunk {nc}, {len(ds)} naive, {len(ss)} superdyads")
        
        all_supers.extend(ss)
        all_supers = build_superdyads(all_supers)
        print(f"chunk {nc}, {len(ss)} combined superdyads")
    
    return all_supers

def remove_subdyads(dyads):
    
    dyads = list(sorted(dyads, key = lambda a:-a[0]))
    out = [dyads.pop(0)]
    
    while len(dyads) > 0:
        slay = True
        db = dyads.pop(0)
        for da in out:
            if check_subdyad(da, db):
                slay = False
                break
            else:
                pass
        if slay:
            out.append(db)
    
    return out

def build_superdyads(dyads):
    
    if len(dyads) < 2:
        return dyads
    
    ncomps = 0
    
    dyads = list(sorted(dyads, key = lambda a:-a[0]))
    out = [dyads.pop(0)]
    while len(dyads) > 0:
        slay = True        
        db = dyads.pop(0)
        for i in range(len(out)):
            da = out[i]
            res, dc = check_mutual_subdyad(da, db)
            ncomps += 1
            if res:
                slay = False
                break
            else:
                pass
        if slay:
            out.append(db) # dyad b is mutually exclusive dyad
        else:
            out[i] = dc # dyad a and b have mutual super, overwrite a with c
    print(f"compared {ncomps} dyads")
    
    return out

def check_subdyad(da, db):
    """
    True if db is a subdyad of da (or vice versa)
    """
    if db[0] > da[0]:
        da,db = db, da
    elif db[0] == da[0]:
        return da == db
    
    dla, dsa, lla, ddsa = da
    dlb, dsb, llb, ddsb = db
    res = True
    res = res and dsb >= dsa
    res = res and dlb + dsb <= dla + dsa
    res = res and ddsb >= ddsa
    res = res and dlb + ddsb <= dla + ddsa 
    
    return res

def check_subdyad_fast(da, db):
    """
    True if db is a subdyad of da (or vice versa)
    """
    if db[0] > da[0]:
        da,db = db, da
    elif db[0] == da[0]:
        return da == db
    
    dla, dsa, lla, ddsa = da
    dlb, dsb, llb, ddsb = db
    ca = dla + dsa + lla / 2
    cb = dlb + dsb + llb / 2
    res = True
    res = res and dsb >= dsa
    res = res and dlb + dsb <= dla + dsa
    res = res and ca == cb
    return res

def check_mutual_subdyad(da, db):
    if db[0] > da[0]:
        da,db = db, da
    
    dla, dsa, lla, ddsa = da
    dlb, dsb, llb, ddsb = db
    
    ca = dla + dsa + lla / 2
    cb = dlb + dsb + llb / 2
    
    if ca == cb:
        dsc = min(dsa, dsb)
        llc = min(lla, llb)
        ddsc = min(ddsa, ddsb)
        dlc = ddsc - llc - dsc
        return True, (dlc, dsc, llc, ddsc)
    else:
        return False, (-1, -1, -1, -1)

