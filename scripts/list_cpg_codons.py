
import itertools

from ggene.seqs import bio, vocab
from ggene.seqs.bio import CODON_TABLE
from ggene import draw


def count_cpg_ncodons(n = 3):
    
    codons = list(CODON_TABLE.keys())
    aas = []
    for cdn in codons:
        aa = CODON_TABLE[cdn]
        if aa not in aas:
            aas.append(aa)
    
    num_cdns = len(codons)
    
    combos = itertools.product(*[range(num_cdns) for n in range(n)])
    
    cg_freqs = {cdn:0 for cdn in codons}
    # polyaa_freqs = {(aa1, aa2): 0 for aa1, aa2 in itertools.product(*[aas for n in range(n)])}
    
    for cdns in combos:
        seq = []
        for icdn in cdns:
            seq.append(codons[icdn])
            
        full_seq = "".join(seq)
        ncgs= full_seq.count("CG")
        for icdn in cdns:
            cdn = codons[icdn]
            cg_freqs[cdn] += ncgs
        
        # polyaa_freqs[tuple(CODON_TABLE[cdn] for cdn in cdns)] = ncgs
    
    baseline = min(cg_freqs.values())
    baseres = min([f - baseline for f in cg_freqs.values() if f > baseline])
    
    cg_freqs = {k: (v - baseline) / baseres for k, v in cg_freqs.items()}
    
    print_frm = "{:} {:<6} {:<6.1f}"
    
    print(f"{n}-codon cpg frequency")
    aa_freqs = {k:0 for k in CODON_TABLE.values()}
    for f, n in cg_freqs.items():
        if n < 1:
            continue
        aa_freqs[CODON_TABLE[f]] += n
        print(print_frm.format(f, CODON_TABLE[f], n, n))
    
    print("aa cpg frequency")
    for aa, n in aa_freqs.items():
        if n < 1:
            continue
        rel = (n - baseline) / baseres
        print(f"{aa}: {n}, {rel:0.1f}")
    
    hmdata = []
    row_lbls = []
    col_lbls = []
    val_lbls = []
    
    for a, b in itertools.product(range(4), range(4)):
        dn = vocab.VOCAB[a] + vocab.VOCAB[b]
        row_lbls.append(dn)
        
        row  = []
        lbl_row = []
        for c in range(4):
            mn = vocab.VOCAB[c]
            freq = cg_freqs[dn+mn]
            row.append(freq)
            
            if a==0 and b==0:
                col_lbls.append(mn)
        
            lbl_row.append(CODON_TABLE[dn+mn])
        
        hmdata.append(row)
        val_lbls.append(lbl_row)
    
    draw.heatmap(
        hmdata, 
        row_labels = row_lbls, 
        col_labels = col_lbls, 
        value_labels = val_lbls,
        show_values = True,
        center = None, 
        
        color_scheme = "energy", 
        symmetric_color = False, 
        
        col_width = 4,
        row_height = 2,
    
    )
            
        
    


def list_cpg_codons():

    cdns ={}
    cdn_freqs = {k:0 for k in CODON_TABLE}
    
    print("single codon CG's")
    for k in CODON_TABLE:
        if "CG" in k:
            cdn_freqs[k] += 1
            cdns[k] = CODON_TABLE[k]
            
    print("double codon CG's:")
    for k in CODON_TABLE:
        for kk in CODON_TABLE:
            if k[-1] == "C" and kk[0] == "G":
                cdn_freqs[k] += 1
                cdn_freqs[kk] += 1
                cdns[k + kk] = CODON_TABLE[k] + CODON_TABLE[kk]
    
    for k, v in cdns.items():
        colstr, key = draw.highlight_sequence(k, "CG", colors = {"CG":136}, suppress = True, show_key = False)
        print(f"{"".join(colstr)}: {v}")
    
    print("codon cpg frequency")
    aa_freqs = {k:0 for k in CODON_TABLE.values()}
    for f, n in cdn_freqs.items():
        if n < 1:
            continue
        aa_freqs[CODON_TABLE[f]] += n
        print(f"{f}, {CODON_TABLE[f]}: {n}")
    
    print("aa cpg frequency")
    for aa, n in aa_freqs.items():
        if n < 1:
            continue
        print(f"{aa}: {n}")
    
    pass


def main():
    

    count_cpg_ncodons(n=2)

    
    pass





if __name__=="__main__":
    main()

