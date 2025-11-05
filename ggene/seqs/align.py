
import random

from ggene.seqs import process
from ggene.seqs import bio

def initialize_align(seqa, seqb, length = 3):
    
    prba = seqa[:length]
    prbb = seqb[:length]
    abstart =  len(seqa)
    bastart = len(seqb)
    
    for i in range(len(seqb) - length):
        bsub = seqb[i:i+length]
        if prba == bsub:
            abstart = i
            break
    
    for i in range(len(seqa) - length):
        asub = seqa[i:i+length]
        if prbb == asub:
            bastart = i
            break
    
    return abstart, bastart

# def overlap_sequences(seqa, seqb):
#     # like align but no insertion/deletion
    
    
    
    
#     pass

def align_sequences(seqa, seqb):
    # brute force
    
    alen = len(seqa)
    blen = len(seqb)
    i = j = 0
    diffs = []
    
    iab, iba = initialize_align(seqa, seqb, length=3)
    if iab < blen and iab < iba:
        j = iab
        diffs.append((0, iab))
    elif iba < alen and iba < iab:
        i = iba
        diffs.append((iba, 0))
    
    algn = True
    while True:
        
        ba = seqa[i]
        bb = seqb[j]
        
        print(ba, bb)
        
        if ba == bb:
            if not algn:
                algn = True
                diffs.append((i, j))
            i+=1
            j+=1
        else:
            # how to decide to incr i or j?
            # or subs?
            j += 1
            # if random.random() > 0.5:
            #     j+=1
            #     print("inc j")
            # else:
            #     i+=1
            #     print("inc i")
            algn = False
        
        if i >= alen or j >= blen:
            break
    diffs.append((alen, blen))
    
    return diffs

def align_sequences2(seqa, seqb):
        
    alen = len(seqa)
    blen = len(seqb)
    diffs = []
    
    cvreslen, cvresinds  = process.correlate_longest_subseq(seqa, seqb)
    maxlen = max(cvreslen)
    maxind = cvresinds[cvreslen.index(maxlen)]
    ilo, jlo = maxind
    ihi, jhi = maxind + maxlen
    diffs.append((ilo, jlo))
    diffs.append((ihi, jhi))
    
    reshi = align_sequences(seqa[ilo:], seqb[jlo:])
    rseqa = bio.reverse(seqb)
    rseqb = bio.reverse(seqa)
    
    reslo = align_sequences(seqa[ilo:], seqb[jlo:])
    reshi_rev = align_sequences(rseqa[alen-ihi:], rseqb[blen-jhi:])
    reslo = [(a+ilo, b+jlo) for a, b in reslo]
    reshi = [(alen - a + ilo, blen - b+ jlo) for a,b in reshi_rev]
    
    diffs.extend(reslo)
    diffs.extend(reshi)
    return list(sorted(diffs, key = lambda a:a[0]))
    
def align_sequences3(seqa, seqb, scale = 8):
        
    alen = len(seqa)
    blen = len(seqb)
    diffs = []
    
    init_scale = min(alen//8, blen//8)
    
    cvreslen, cvresinds, sh  = process.correlate_longest_subseq(seqa, seqb, scale = init_scale)
    maxlen = max(cvreslen)
    maxind = cvresinds[cvreslen.index(maxlen)]
    maxsh = sh[cvreslen.index(maxlen)]
    diffs.append((maxind, maxind+maxsh))
    diffs.append((maxind+maxlen, maxind+maxlen+maxsh))
    
    print("first pass")
    print(maxlen, maxind, maxsh)
    print(seqa[maxind:maxind+maxlen])
    print(seqa[maxind - 10:maxind+maxlen + 10])
    print(seqb[maxind+maxsh - 10:maxind+maxlen+maxsh + 10])
    
    # forward
    new_scale = init_scale//4
    off = maxind+maxlen
    fseqa = seqa[off:]
    fseqb = seqb[off+maxsh:]
    lsa, isa, ssa = process.correlate_longest_subseq(fseqa[:new_scale], fseqb, scale = new_scale)
    lsb, isb, ssb = process.correlate_longest_subseq(fseqa, fseqa[:new_scale], scale = new_scale)
    
    print("second pass")
    print(fseqa[:new_scale])
    print(fseqb[:new_scale])
    
    mla = max(lsa)
    mia = isa[lsa.index(mla)]
    msa = ssa[lsa.index(mla)]
    mlb = max(lsb)
    mib = isb[lsb.index(mlb)]
    msb = ssb[lsb.index(mlb)]   
    
    if mla > mlb:
        ni = mia
        ns = msa
        nl = mla
    elif mlb > mla:
        ns = msb
        ni = mib - ns
        nl = mlb
    else:
        nl = mla
        if msa < msb:
            ni = mia
            ns = msa
        else:
            ns = msb
            ni = mib - ns
    
    diffs.append((ni, ni+ns))
    diffs.append((ni + nl + off, ni+ns + nl + off))
    
    # reverse
    
    
    return diffs

def fill_aligned(seqa, seqb, diffs, fill = '-'):
    afill = []
    bfill = []
    i=j=0
    for da, db in diffs:
        afill.append(seqa[i:da])
        bfill.append(seqb[j:db])
        if da == len(seqa) or db == len(seqb):
            break
        diff = da - db
        if diff < 0:
            afill.append(fill*abs(diff))
        elif diff > 0:
            bfill.append(fill*abs(diff))
        else:
            print("weird")
        i, j = da, db
    
    return "".join(afill), "".join(bfill)

def print_aligned(seqa, seqb, diffs):
    print(f"sequence a: {seqa}")
    print(f"sequence b: {seqb}")
    print(f"{len(diffs)} differences")
    
    i = j = 0
    for a, b in diffs:
        print(seqa[i:a])
        print(seqb[j:b])
        i = a
        j = b
    print(seqa[i:])
    print(seqb[j:])
    
    

