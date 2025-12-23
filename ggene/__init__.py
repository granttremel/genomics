
"""Genomics Gene Analysis Library.

A comprehensive library for analyzing genomic data including VCF files,
GTF annotations, and gene structures for genome browser applications.
"""

import logging

from ggene.database.genome_manager import GenomeManager
from ggene.genome.features import Feature, Gene, shorten_variant
from ggene.seqs.bio import COMPLEMENT_MAP, CODON_TABLE, START_CODONS, STOP_CODONS, to_rna, to_dna, is_rna, is_dna, reverse_complement, complement
from ggene.config import get_config, get_paths

__all__ = ['logger', 
            'get_config','get_paths',
            'GenomeManager',  'Feature', 'Gene', 'shorten_variant',
            'COMPLEMENT_MAP', 'CODON_TABLE','START_CODONS','STOP_CODONS','to_rna', 'to_dna', 'is_rna', 'is_dna', 'reverse_complement', 'complement']

# Configure logging
logger = logging.getLogger(__name__)


# End of module

