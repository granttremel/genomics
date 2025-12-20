
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

from ggene.database.genome_manager import GenomeManager
from ggene.genome.features import Feature, Gene, shorten_variant
from ggene.seqs.bio import COMPLEMENT_MAP, CODON_TABLE, START_CODONS, STOP_CODONS, to_rna, to_dna, is_rna, is_dna, reverse_complement, complement
from ggene.config import get_config, get_paths

__all__ = ['logger', 
        #    'cfg', 'DATA_DIR', 'get_paths', 'other_paths', 
        'get_config','get_paths',
           'GenomeManager',  'Feature', 'Gene','shorten_variant', 'COMPLEMENT_MAP', 'CODON_TABLE','START_CODONS','STOP_CODONS','to_rna', 'to_dna', 'is_rna', 'is_dna', 'reverse_complement', 'complement']

# Configure logging
logger = logging.getLogger(__name__)

# CONFIG_PATH = "./local.yaml" if os.path.exists("./local.yaml") else "./default.yaml"

# cfg = {}
# with open(CONFIG_PATH) as f:
#     cfg = yaml.safe_load(f)

# Constants
# DATA_DIR = Path(cfg.get("data_dir"), )
# DEFAULT_VCF_PATH = DATA_DIR / cfg.get("vcf_name","")
# DEFAULT_GTF_PATH = DATA_DIR / "GRCh38_sorted.gtf.gz"
# DEFAULT_FASTA_PATH = DATA_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
# DEFAULT_LIBRARY = DATA_DIR / "library"

# def get_paths():
#     return str(DEFAULT_VCF_PATH), str(DEFAULT_GTF_PATH), str(DEFAULT_FASTA_PATH), str(DEFAULT_LIBRARY)

# other_paths = {
#     "repeatmasker_path":  DATA_DIR / "repeatmasker" / "repeats.sorted.bed.gz",
    # "dfam_path": DATA_DIR / "dfam/hg38_dfam.nrph.bed.gz"
# }

# End of module

