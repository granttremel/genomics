
from collections import defaultdict, Counter
import numpy as np

def count_kmers(seq, k):
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    return Counter(kmers)

def find_kmers(seq, max_k = 15):

    mykmers = {}
    for k in range(3,max_k):
        top_kmers = count_kmers(seq, k).most_common(3)
        for s,k in top_kmers:
            mykmers[s] = k
        
    return mykmers

def scrub_kmers(kmers):
    
    ks = list(sorted(kmers.keys(), key = lambda k:len(k)))
    
    cont_score = {}
    eff_score = {}
    
    for i in range(len(ks)):
        k = ks[i]
        v = kmers[k]
        for j in range(i+1, len(ks)):
            kk = ks[j]
            vv = kmers[kk]
            
            #check if permute
            if len(k) == len(kk):
                kkp = kk
                prm = False
                for n in range(len(k)):
                    kkp = permute(kkp)
                    if k == kkp:
                        prm = True
                        break
                if prm:
                    pass
            
            
            #check if contains
            
        
    
    
    pass

def permute(seq):
    
    return seq[1:] + seq[0]

def find_universal_motifs(microsats, min_frequency=0.1):
    """Find motifs that appear across multiple microsatellites"""
    
    # Collect all k-mers across all microsatellites
    kmer_positions = defaultdict(list)  # kmer -> [(microsat_idx, position), ...]
    
    for k in range(3, 15):  # Try different motif lengths
        for idx, ms in enumerate(microsats):
            seq = ms['seq']
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_positions[kmer].append((idx, i))
    
    # Find motifs present in multiple microsatellites
    universal_motifs = {}
    for kmer, positions in kmer_positions.items():
        unique_microsats = len(set(pos[0] for pos in positions))
        if unique_microsats >= len(microsats) * min_frequency:
            universal_motifs[kmer] = {
                'prevalence': unique_microsats / len(microsats),
                'avg_count': len(positions) / unique_microsats,
                'positions': positions
            }
    
    return universal_motifs


def build_symbol_alphabet(kmers):
    """Assign symbols to meaningful motifs"""
    
    # Assign symbols to top motifs (avoiding overlaps)

    symbols = {}
    used_chars = set()
    
    for motif in kmers:  # Top 50 motifs
        # Assign single char or digraph symbol
        if len(symbols) < 26:
            symbol = chr(65 + len(symbols))  # A-Z
        else:
            symbol = f"{chr(97 + (len(symbols)-26)//26)}{chr(97 + (len(symbols)-26)%26)}"
        
        symbols[motif] = symbol
    
    return symbols

def encode_with_alphabet(seq, symbols, allow_overlap = False):
    
    motifs = sorted(symbols.items(), key=lambda x: len(x[0]), reverse=True)
    
    encoded = []
    i = 0
    while i < len(seq):
        matched = False
        
        for motif, symbol in motifs:
            if seq[i:i+len(motif)] == motif:
                encoded.append(symbol)
                i += len(motif)
                matched = True
                break
        
        if not matched:
            encoded.append('.')  # Unknown/variable position
            i += 1
    
    return ''.join(encoded)