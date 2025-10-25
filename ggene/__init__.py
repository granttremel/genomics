
"""Genomics Gene Analysis Library.

A comprehensive library for analyzing genomic data including VCF files,
GTF annotations, and gene structures for genome browser applications.
"""

import gzip
import logging
import os
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple, Union, Iterator, Any, Callable
import numpy as np
import pysam
import re
import pyranges as pr


__all__ = ['logger', 'DATA_DIR','DEFAULT_VCF_PATH', 'DEFAULT_GTF_PATH', 'COMPLEMENT_MAP', 'Feature', 'Gene','shorten_variant']

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = "./data"
DEFAULT_VCF_PATH = "./data/gt_vcf.gz"
DEFAULT_GTF_PATH = "./data/GRCh38_sorted.gtf.gz"
DEFAULT_FASTA_PATH = "./data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
DEFAULT_LIBRARY = "./data/library"

INFO_TAGS = (
    'AF1', 'DP', 'DP4', 'FQ', 'MQ', 'PC2', 'PCHI1', 
    'PV4', 'QCHI2', 'RP', 'CLR', 'UGT', 'VDB', 'RPB', 'HWE'
)

COMPLEMENT_MAP = {'A':'T','T':'A','C':'G','G':'C','U':'A','N':'N'}
COMPLEMENT_MAP_RNA = {'A':'U','T':'A','C':'G','G':'C','U':'A','N':'N'}
CODON_TABLE_DNA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def to_rna(seq):
    return seq.replace("T","U")
def to_dna(seq):
    return seq.replace("U","T")

def is_rna(seq):
    
    return "U" in seq and not "T" in seq

def is_dna(seq):
    return "T" in seq and not "U" in seq

CODON_TABLE_RNA = {to_rna(k):v for k,v in CODON_TABLE_DNA.items()}

START_CODONS = ['ATG']  # Standard start codon (Methionine)
ALTERNATIVE_START_CODONS = ['CTG', 'GTG', 'TTG']  # Rare alternative starts
STOP_CODONS = ['TAA', 'TAG', 'TGA']

aliases_dna = {
    'Y':'CT',
    'R':'AG',
    'N':'ATCG'
    }

aliases_rna = {
    'Y':'CU',
    'R':'AG',
    'N':'AUCG'
    }

aliases_rev = {
    'C':'YN',
    'T':'YN',
    'A':'RN',
    'G':'RN',
    'U':'YN',
    }


splice_donor = 'GG*GURAGU'
splice_branch = 'YURAC'
splice_acceptor = 'YNCAG*G'


def reverse_complement(seq, rna = True):
    
    return "".join(reversed(complement(seq, rna=rna)))

def complement(seq, rna = True):
    
    if rna:
        newseq = "".join([COMPLEMENT_MAP_RNA.get(s,s) for s in seq])
    else:
        newseq = "".join([COMPLEMENT_MAP.get(s,s) for s in  seq])
    
    return newseq

def seq_to_ptrn(seq, rna=True):
    
    if rna:
        ref=aliases_rna
    else:
        ref=aliases_dna
    
    re_comps = []
    for i in range(len(seq)):
        
        b = seq[i]
        
        if i+1 < len(seq):
            if seq[i+1]=='*':
                re_comps.append('(')

        if b in 'ATCGU':
            re_comps.append(b)
        elif b in 'YRN':
            restr = '|'.join(ref[b])
            re_comps.append(f'[{restr}]')
            
        if i > 0:
            if seq[i-1]=='*':
                re_comps.append(')')
        
    re_ptrn = ''.join(re_comps)
    return re_ptrn

def make_re(seq, rna=True):
    
    re_ptrn = seq_to_ptrn(seq, rna=rna)
            
    return re.compile(re_ptrn)

def make_splice_re(rna=True):
    
    donor_ptrn = seq_to_ptrn(splice_donor)
    branch_ptrn = seq_to_ptrn(splice_branch)
    acc_ptrn = seq_to_ptrn(splice_acceptor)
    
    bracc_ptrn = ''.join((branch_ptrn,'.{20,50}',acc_ptrn))
    
    return re.compile(donor_ptrn), re.compile(bracc_ptrn)


FEATURE_TYPES = [
    'gene', 'CDS', 'transcript', 'five_prime_utr', 'start_codon',
    'exon', 'stop_codon', 'three_prime_utr', 'variant'
]



def variant_dict_to_array(variant_map: Dict[int, Dict[int, int]]) -> np.ndarray:
    """Convert variant dictionary to numpy array.
    
    Args:
        variant_map: Nested dictionary with variant counts
        
    Returns:
        2D numpy array with variant data
    """
    if not variant_map:
        return np.array([])
        
    max_ref = max(variant_map.keys())
    max_alt = max(max(alt_dict.keys()) for alt_dict in variant_map.values())
                
    var_array = np.zeros((max_ref + 1, max_alt + 1))
    
    for ref_len, alt_dict in variant_map.items():
        for alt_len, count in alt_dict.items():
            var_array[ref_len, alt_len] = count
    
    return var_array

def shorten_variant(var: Any, keys: Optional[List[str]] = None) -> Dict[str, Any]:
    """Extract key variant information into a dictionary.
    
    Args:
        var: Variant object from cyvcf2
        keys: Additional attributes to extract
        
    Returns:
        Dictionary with variant information
    """
    if keys is None:
        keys = []
        
    out_dict = {
        'pos': var.POS,
        'ref': var.REF,
        'alt': var.ALT[0] if var.ALT else None,
    }
    
    for key in keys:
        if hasattr(var, key):
            out_dict[key] = getattr(var, key)
    
    return out_dict

               
class Motif:

    _counter = {}

    def __init__(self, name, chrom, start, end, score):
        self.type="motif"
        self.name = name
        ind = self._counter.get(name,0)
        self._counter[name] = ind + 1
        self.sfid=f'{name}-{ind}'
        self.chrom=chrom
        self.start=start
        self.end=end
        self.relative_start = 0
        self.relative_end = self.end - self.start
        
        self.score=score
        
        self.parent_gene=None
        self.parents=[]
        
    def set_parent_gene(self, parent_gene):
        self.parent_gene=parent_gene
        
    def set_parent(self, feature):
        if not feature in self.parents:
            self.parents.append(feature)
            
    def refer_to(self, reference):
        self.start_relative = self.start - reference.start
        self.end_relative = self.end - reference.start
    def __repr__(self) -> str:
        """String representation of the feature."""
        
        parts=[f"id={self.sfid}"]
        parts.append(f"{self.chrom}:{self.start}-{self.end}")
        parts.append(f"score={self.score}")
            
        content = ','.join(parts)
        return f'{type(self).__name__}({content})'

# End of module

