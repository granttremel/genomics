
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
import yaml
from pathlib import Path

from .genome.features import Feature, Gene

__all__ = ['logger', 'DATA_DIR','DEFAULT_VCF_PATH', 'DEFAULT_GTF_PATH', 'COMPLEMENT_MAP', 'Feature', 'Gene','shorten_variant']

# Configure logging
logger = logging.getLogger(__name__)

CONFIG_PATH = "./local.yaml" if os.path.exists("./local.yaml") else "./default.yaml"

cfg = {}
with open(CONFIG_PATH) as f:
    cfg = yaml.safe_load(f)

# Constants
DATA_DIR = Path(cfg.get("data_dir"), )
DEFAULT_VCF_PATH = DATA_DIR / "gt_vcf.gz"
DEFAULT_GTF_PATH = DATA_DIR / "GRCh38_sorted.gtf.gz"
DEFAULT_FASTA_PATH = DATA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
DEFAULT_LIBRARY = DATA_DIR / "library"

def get_paths():
    return str(DEFAULT_VCF_PATH), str(DEFAULT_GTF_PATH), str(DEFAULT_FASTA_PATH), str(DEFAULT_LIBRARY)

other_paths = {
    "repeatmasker_path":  DATA_DIR / "repeatmasker" / "repeats.sorted.bed.gz",
    "dfam_path": DATA_DIR / "dfam/hg38_dfam.nrph.bed.gz"
}

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

START_CODONS = ['ATG']  # Standard start codon (Methionine) = AUG
ALTERNATIVE_START_CODONS = ['CTG', 'GTG', 'TTG']  # Rare alternative starts = CUG, GUG, UUG
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


FEATURE_TYPES = [
    'gene', 'CDS', 'transcript', 'five_prime_utr', 'start_codon',
    'exon', 'stop_codon', 'three_prime_utr', 'variant'
]

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

# End of module

