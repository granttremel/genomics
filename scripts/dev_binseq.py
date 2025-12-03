
import random
from tabulate import tabulate

from ggene.seqs import binary, bio
from ggene.seqs.binary import BinBase, BinSeq

def show_iteration(n_max = 128):
    
    seqs = []
    
    for n in range(n_max):
        seq = bio.index_to_seq_abs(n)
        binseq = binary.to_binary(seq)
        print("{:<8}{:<8}{:<8}".format(str(n), seq, str(binseq)))
        for b in binseq:
            print(bin(b))
        seqs.append((n, seq, binseq))
    print()
    
    seqs_srt = sorted(seqs, key = lambda k: k[2][0])
    print("sorted by binary")
    for n, seq, binseq in seqs_srt:
        print("{:<8}{:<8}{:<8}".format(str(n), seq, str(binseq)))
        for b in binseq:
            print(bin(b))
    print()

# def test_binbase():
    
#     bbnull = BinBase.NULL
#     bbc = BinBase.C
#     bbt = BinBase.T
#     bbg = BinBase.G
#     bba = BinBase.A
#     bbn = BinBase.N
    
#     print("A or C:", bba or bbc)
#     print("G or null:", bbg or bbnull)
#     print("G and null:", bbg and bbnull)
#     print("T and N:", bbt and bbn)
    
def test_onehot():
    ohdict =  binary.get_onehot_fwd()
    
    for k, v in ohdict.items():
        print(k, v)
        
    ohdict_rev =  binary.get_onehot_rev()
    for k, v in ohdict_rev.items():
        print(k, v)

def test_binbase():
    
    bases = list(binary.ONE_HOT.keys())
    
    and_table = []
    or_table = []
    
    for b in bases:
        bb = BinBase(b)
        arow = [str(bb)]
        orow = [str(bb)]
        for b2 in bases:
            bb2 = BinBase(b2)
            arow.append(str(bb.logical_and(bb2)))
            orow.append(str(bb.logical_or(bb2)))
        and_table.append(arow)
        or_table.append(orow)
    
    print("AND")
    print(tabulate(and_table, headers=bases))
    print()
    
    print("OR")
    print(tabulate(or_table, headers=bases))
    print()

def test_binseq():
    
    rand_seq = random.choices(list(binary.base_inds.keys()), k = 16)
    rand_seq2 = random.choices(list(binary.base_inds.keys()), k = 16)
    
    binseq = BinSeq(rand_seq)
    binseq2 = BinSeq(rand_seq2)
    
    print(binseq)
    print(binseq2)
    print(binseq.logical_and(binseq2))
    print(binseq.logical_or(binseq2))
    
def test_process():
    # num_rounds = 32
    num_rounds = 5
    seq_len = 16
    step = 3
    
    bs1 = BinSeq(random.choices(list(binary.base_inds.keys()), k = seq_len))
    bs2_str = random.choices(list(binary.base_inds.keys()), k = (seq_len + num_rounds * step)//2)
    bs2 = BinSeq(bs2_str + bs2_str[::-1])
    bs2= bs2.subseq(stop = seq_len//2).to_dyad()
    print(f"starting bs1: {bs1}")
    print(f"starting bs2: {bs2}")
    
    bs1_series = []
    bs2_series = []
    
    subseq = bs2
    # print(f"subseq: {subseq}")
    for i in range(num_rounds):
        bs2_series.append(subseq)
        bs1 = (bs1 ^ subseq).reverse_complement()
        # bs1 = bs1 ^ bs3.reverse_complement()
        bs1_series.append(bs1)
        
    print("bs1 series")
    for bseq in bs1_series:
        print(bseq, bseq.multiplicity)
    print()

def test_binary():
    
    test_bs = BinSeq("ATGCGGACATGGGGGG")
    print(test_bs)
    print(test_bs.reverse())
    print(test_bs & test_bs.reverse())
    print(test_bs.logical_and(test_bs.reverse()))
    print(test_bs | test_bs.reverse())
    print(test_bs.logical_or(test_bs.reverse()))
    print(test_bs ^ test_bs.reverse())
    print(test_bs.reverse() ^ test_bs)
    
    print()
    print(test_bs)
    print(test_bs.reverse())
    print(test_bs.complement())
    print(test_bs.reverse_complement())

def test_binary_ops():
    
    
    bs1 = BinSeq(random.choices(list(binary.base_inds.keys()), k = 16))
    bs2 = BinSeq(random.choices(list(binary.base_inds.keys()), k = 8))
    
    bs2_double = bs2.double()
    bs2_concat = bs2.concat(bs2_double)
    bs2_rpt3 = bs2.repeat(3)
    bs2_rtrunc = bs2.rtrunc(4)
    
    print("bs1",bs1)
    print("bs2",bs2)
    print("bs2 double",bs2_double)
    print("bs2 concat",bs2_concat)
    print("bs2 repeat",bs2_rpt3)
    print("bs2 rtrunc",bs2_rtrunc)
    
    consts = binary.get_consts(nbits = 16, max_n = 5)
    for (ord, par), const in consts.items():
        print(f"const({ord},{par})", const)
    
    print(bs1)
    print(bs1.rotate(1))
    print(bs1.rotate(4))
    print(bs1.rotate(-4))


def bit_twiddle_b16(bseq):
    
    nbits = len(bseq)
    cs = binary.get_consts(nbits)

    # Swap adjacent bits
    bseq = ((bseq & cs[(1,1)]) << 1) | ((bseq & cs[(1,0)]) >> 1)
    print(bseq)
    # Swap pairs of bits
    bseq = ((bseq & cs[(2,1)]) << 2) | ((bseq & cs[(2,0)]) >> 2)
    print(bseq)
    # Swap nibbles
    bseq = ((bseq & cs[(3,1)]) << 4) | ((bseq & cs[(3,0)]) >> 4)
    print(bseq)
    # Swap bytes (for 16-bit and larger)
    bseq = ((bseq & cs[(4,1)]) << 8) | ((bseq & cs[(4,0)]) >> 8)
    print(bseq)

    return bseq.rtrunc(nbits)

def try_bit_twiddle():
    
    
    seq_len = 16
    
    bs1 = BinSeq(random.choices(list(binary.base_inds.keys()), k = seq_len))
    
    bs1_res = bit_twiddle_b16(bs1)
    
    print(bs1)
    print(bs1_res)
    
    
    
    
    
    
    
    pass

def find_tern_consts():
    
    consts = binary.get_consts(16 * 4)
    
    for n in range(1, 5):
        c10 = consts[(n,0)]
        print(c10)
        for p in range(4):
            
            c10p = c10.rotate(n*p)
            c10_str = "".join(["1" if str(c)=="N" else "0" for c in c10p])
            
            newc = []
            
            for i in range(len(c10_str)//4):
                c10_sub = list(c10_str[i:i+4])
                b = BinBase.from_onehot(map(int,c10_sub))
                newc.append(str(b))
            
            newseq = BinSeq("".join(newc))
            print(newseq, p, n)
    
    
    pass

def main():
    
    # binary.set_binary_map(4)
    
    # blen, bwidth = binary.get_binary_data()
    # print(f"number bin sequences: {blen}; binary width: {bwidth}")
    
    # show_iteration()
    
    # test_onehot()
    # test_binbase()
    # test_binseq()
    # test_process()
    # test_binary_ops()
    find_tern_consts()
    # try_bit_twiddle()
    
    # test_binary()
    
    
    pass

if __name__=="__main__":
    main()
