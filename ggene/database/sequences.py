"""Sequence streaming with integrated variant application.

This module provides efficient streaming of genomic sequences from FASTA files
with optional variant integration from VCF files.
"""

import pysam
from pathlib import Path
from typing import Optional, Iterator, List, Tuple, Dict, Any, Union, TYPE_CHECKING
import logging
from dataclasses import dataclass
from abc import ABC, abstractmethod

from ggene.database.annotations import UFeature
from ggene.seqs import vocab

if TYPE_CHECKING:
    from ggene.database.annotations import UFeature

logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)


@dataclass
class SequenceWindow:
    """Represents a window of genomic sequence."""
    chrom: str
    start: int  # 1-based
    end: int    # 1-based, inclusive
    reference: str
    personal: str  # With variants applied
    variants: List[Any] = None
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1
    
    @property
    def has_variants(self) -> bool:
        return self.reference != self.personal


class BaseSequenceStream(ABC):
    """Abstract base class for sequence streaming."""
    
    @abstractmethod
    def get_sequence(self, ref: str, start: int = None, end: int = None) -> str:
        """Get reference sequence for a region."""
        pass
    
    @abstractmethod
    def get_personal_sequence(self, chrom: str, start: int, end: int,
                            sample: Optional[str] = None) -> str:
        """Get sequence with variants applied."""
        pass


class FASTAStream(BaseSequenceStream):
    """Stream sequences from indexed FASTA files."""
    
    def __init__(self, filepath: str):
        """Initialize FASTA stream.
        
        Args:
            filepath: Path to FASTA file (should be indexed with .fai)
        """
        self.filepath = Path(filepath)
        self.fasta = None
        
        # Check if file exists
        if not self.filepath.exists():
            raise FileNotFoundError(f"FASTA file not found: {filepath}")
        
        # Open FASTA file
        try:
            self.fasta = pysam.FastaFile(str(self.filepath))
            logger.info(f"Opened FASTA file: {filepath}")
            
            # Check if indexed
            fai_file = Path(str(self.filepath) + '.fai')
            if not fai_file.exists():
                logger.warning(f"No index found for {filepath}. Creating index...")
                # pysam will create the index automatically
                
        except Exception as e:
            logger.error(f"Failed to open FASTA file: {e}")
            raise
    
    def get_sequence(self, ref: str, start: int=None, end: int=None) -> str:
        """Get reference sequence for a region.
        
        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            
        Returns:
            DNA sequence string
        """
        if not self.fasta:
            return ""
        
        try:
            # pysam uses 0-based half-open intervals
            seq = self.fasta.fetch(str(ref), start - 1, end)
            seq = seq.upper()
            if not vocab.VOCAB == vocab.VOCAB_DNA:
                seq = vocab._convert_seq_vocab(seq, vocab.VOCAB, from_vocab = vocab.VOCAB_DNA)
            return seq
        except Exception as e:
            logger.debug(f"Failed to fetch sequence {ref}:{start}-{end}: {e}")
            return ""
    
    def get_personal_sequence(self, chrom: str, start: int, end: int,
                            sample: Optional[str] = None) -> str:
        """Get reference sequence (no variants in base FASTA stream)."""
        return self.get_sequence(chrom, start, end)
    
    def get_chromosomes(self) -> List[str]:
        """Get list of available chromosomes."""
        if self.fasta:
            return list(self.fasta.references)
        return []
    
    def get_chromosome_length(self, chrom: str) -> int:
        """Get length of a chromosome."""
        if self.fasta and chrom in self.fasta.references:
            return self.fasta.get_reference_length(chrom)
        return 0
    
    def close(self):
        """Close the FASTA file."""
        if self.fasta:
            self.fasta.close()
            self.fasta = None


class SequenceStreamWithVariants(BaseSequenceStream):
    """Sequence stream that integrates variants from VCF."""
    
    def __init__(self, fasta_stream: FASTAStream, vcf_path: Optional[str] = None,
                 min_qual: float = 5.0):
        """Initialize sequence stream with variant integration.
        
        Args:
            fasta_stream: FASTAStream object
            vcf_path: Path to VCF file with variants
            min_qual: Minimum quality threshold for variants
        """
        self.fasta_stream = fasta_stream
        self.vcf = None
        self.min_qual = min_qual
        
        # Load VCF if provided
        if vcf_path:
            try:
                from cyvcf2 import VCF
                self.vcf = VCF(vcf_path)
                
                # Set index if available
                index_files = [
                    Path(vcf_path + '.tbi'),
                    Path(vcf_path + '.csi')
                ]
                for idx in index_files:
                    if idx.exists():
                        try:
                            self.vcf.set_index(str(idx))
                            logger.info(f"Using indexed VCF: {vcf_path}")
                        except:
                            pass
                        break
                else:
                    logger.warning(f"No index for VCF {vcf_path}, variant queries will be slow")
                    
            except ImportError:
                logger.error("cyvcf2 not installed. Cannot apply variants.")
            except Exception as e:
                logger.error(f"Failed to load VCF: {e}")
    
    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        """Get reference sequence."""
        return self.fasta_stream.get_sequence(chrom, start, end)
    
    def get_variants_in_region(self, chrom: str, start: int, end: int) -> List[Any]:
        """Get variants in a region."""
        if not self.vcf:
            return []
        
        variants = []
        try:
            # Try with chr prefix first, then without
            chrom_str = str(chrom)
            if not chrom_str.startswith('chr'):
                region = f"chr{chrom}:{start}-{end}"
            else:
                region = f"{chrom}:{start}-{end}"
                
            try:
                for var in self.vcf(region):
                    if var.QUAL and var.QUAL >= self.min_qual:
                        variants.append(var)
            except:
                # Try alternate format if first fails
                if not chrom_str.startswith('chr'):
                    region = f"{chrom}:{start}-{end}"
                else:
                    region = f"{chrom.lstrip('chr')}:{start}-{end}"
                    
                for var in self.vcf(region):
                    if var.QUAL and var.QUAL >= self.min_qual:
                        variants.append(var)
        except Exception as e:
            logger.debug(f"Failed to fetch variants for {chrom}:{start}-{end}: {e}")
        
        return variants
    
    def apply_variants_to_sequence(self, seq: str, variants: List[Any], 
                                  seq_start: int,  substitution = "") -> str:
        """Apply variants to a sequence.
        
        Args:
            seq: Reference sequence
            variants: List of cyvcf2 Variant objects
            seq_start: Genomic start position of sequence (1-based)
            
        Returns:
            Sequence with variants applied
        """
        if not variants:
            return seq
        seq_list = list(seq)
        deltas = self.variants_to_deltas(variants)
        return self._apply_deltas_to_sequence(seq_list, deltas, seq_start, substitution=substitution)
    
    def apply_variants_to_data(self, data:List[Any], variants:List[Any], seq_start:int, substitution = ""):
        deltas = self.variants_to_deltas(variants)
        data_del = self._apply_deltas_to_sequence(data, deltas, seq_start, substitution=substitution)
        return data_del
    
    def variants_to_deltas(self, variants:List[Any]):
        if not variants:
            return []
        if isinstance(variants[0], UFeature):
            return [(v.start, len(v.attributes.get("ref")), len(v.attributes.get("alt")[0]), v.attributes.get("alt")[0]) for v in variants]
        else:
            return [(v.POS, len(v.REF), len(v.ALT[0]), v.ALT[0]) for v in variants]
    
    def _apply_deltas_to_sequence(self, seq:List[Any], deltas:List[Tuple[int,int,int,Any]], seq_start:int, substitution = ""):
        """
        delta of form (start_index, reference_length, alternate_length, substitution)
        sub overrides substitution
        """
        
        if not deltas:
            return seq
        
        # Sort variants by position
        deltas.sort(key=lambda v: v[0])
        
        # Track position offset due to indels
        offset = 0
        
        for pos, ref_len, alt_len, sub in deltas:
            
            if substitution:
                sub = substitution * alt_len
            # Skip if no alt allele
            if not sub:
                continue
            
            # Calculate position in sequence (0-based)
            # seq_pos = d.POS - seq_start + offset
            seq_pos = pos - seq_start + offset
            
            # Validate position
            if seq_pos < 0 or seq_pos >= len(seq):
                continue
            
            # Apply variant
            try:
                if ref_len == alt_len:
                    # SNP - simple substitution
                    for i in range(ref_len):
                        if seq_pos + i < len(seq):
                            seq[seq_pos + i] = sub[i]
                            
                elif ref_len > alt_len:
                    # Deletion
                    del_length = ref_len - alt_len
                    # Replace ref with alt
                    for i in range(alt_len):
                        if seq_pos + i < len(seq):
                            seq[seq_pos + i] = sub[i]
                    # Delete extra bases
                    for _ in range(del_length):
                        if seq_pos + alt_len < len(seq):
                            del seq[seq_pos + alt_len]
                    offset -= del_length
                    
                else:
                    # Insertion
                    ins_length = alt_len - ref_len
                    # Replace ref bases
                    for i in range(ref_len):
                        if seq_pos + i < len(seq):
                            seq[seq_pos + i] = sub[i]
                    # Insert additional bases
                    for i in range(ref_len, alt_len):
                        seq.insert(seq_pos + i, sub[i])
                    offset += ins_length
                    
            except Exception as e:
                logger.debug(f"Failed to apply variant at {pos}: {e}")
        
        logger.debug(seq)
        return seq
    
    def get_personal_sequence(self, chrom: str, start: int, end: int,
                            sample: Optional[str] = None) -> str:
        """Get sequence with variants applied.
        
        Args:
            chrom: Chromosome
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            sample: Sample name (optional, for multi-sample VCFs)
            
        Returns:
            Personalized sequence with variants
        """
        # Get reference sequence
        ref_seq = self.get_sequence(chrom, start, end)
        if not ref_seq:
            return ""
        
        # Get and apply variants
        variants = self.get_variants_in_region(chrom, start, end)
        logger.debug(variants)
        if variants:
            return self.apply_variants_to_sequence(ref_seq, variants, start)
        
        return ref_seq
    
    def get_sequence_with_alignment(self, chrom: str, start: int, end: int,
                                   sample: Optional[str] = None) -> Tuple[str, str, List[Tuple[int, int]]]:
        """Get aligned reference and personal sequences with gaps.
        
        This is similar to GenomeIterator._apply_variants_to_window but 
        returns aligned sequences with gaps for visualization.
        
        Returns:
            Tuple of (aligned_ref, aligned_pers, variant_list)
            where variant_list contains tuples of (position, delta)
            and delta = len(ALT) - len(REF)
        """
        fill = '-'
        ref_seq = self.get_sequence(chrom, start, end)
        if not ref_seq:
            return "", "", []
        
        variants = self.get_variants_in_region(chrom, start, end)
        if not variants:
            # No variants, return identical sequences with empty variant list
            return ref_seq, ref_seq, []
        
        # Build aligned sequences with gaps and track variants
        aligned_ref = []
        aligned_pers = []
        variant_list = []  # List of (position, delta) tuples
        
        ref_idx = 0
        aligned_idx = 0
        
        # Sort variants by position
        variants.sort(key=lambda v: v.POS)
        
        for pos in range(start, end + 1):
            # Check for variant at this position
            variant_at_pos = None
            for var in variants:
                if var.POS == pos:
                    variant_at_pos = var
                    break
            
            if variant_at_pos and variant_at_pos.ALT:
                ref_allele = variant_at_pos.REF
                alt_allele = variant_at_pos.ALT[0]
                
                # Calculate delta for this variant
                delta = len(alt_allele) - len(ref_allele)
                if delta != 0:
                    variant_list.append((pos, delta))
                
                if len(ref_allele) == len(alt_allele):
                    # SNP - no gaps needed
                    for i in range(len(ref_allele)):
                        if ref_idx + i < len(ref_seq):
                            aligned_ref.append(ref_seq[ref_idx + i])
                            aligned_pers.append(alt_allele[i] if i < len(alt_allele) else 'N')
                            aligned_idx += 1
                    ref_idx += len(ref_allele)
                    
                elif len(ref_allele) > len(alt_allele):
                    # Deletion - add gaps to personal
                    for i in range(len(alt_allele)):
                        if ref_idx < len(ref_seq):
                            aligned_ref.append(ref_seq[ref_idx])
                            aligned_pers.append(alt_allele[i])
                            aligned_idx += 1
                            ref_idx += 1
                    
                    # Add gaps for deleted bases
                    for i in range(len(ref_allele) - len(alt_allele)):
                        if ref_idx < len(ref_seq):
                            aligned_ref.append(ref_seq[ref_idx])
                            aligned_pers.append(fill)
                            aligned_idx += 1
                            ref_idx += 1
                            
                else:
                    # Insertion - add gaps to reference
                    for i in range(len(ref_allele)):
                        if ref_idx < len(ref_seq):
                            aligned_ref.append(ref_seq[ref_idx])
                            aligned_pers.append(alt_allele[i])
                            aligned_idx += 1
                            ref_idx += 1
                    
                    # Add gaps for inserted bases
                    for i in range(len(ref_allele), len(alt_allele)):
                        aligned_ref.append(fill)
                        aligned_pers.append(alt_allele[i])
                        aligned_idx += 1
                        
            else:
                # No variant at this position
                if ref_idx < len(ref_seq):
                    base = ref_seq[ref_idx]
                    aligned_ref.append(base)
                    aligned_pers.append(base)
                    aligned_idx += 1
                    ref_idx += 1
        
        return ''.join(aligned_ref), ''.join(aligned_pers), variant_list
    
    def stream_windows(self, chrom: str, start: int, end: int, 
                      window_size: int = 100, 
                      stride: Optional[int] = None) -> Iterator[SequenceWindow]:
        """Stream sequence windows with variants.
        
        Args:
            chrom: Chromosome
            start: Start position
            end: End position  
            window_size: Size of each window
            stride: Step size (defaults to window_size)
            
        Yields:
            SequenceWindow objects
        """
        if stride is None:
            stride = window_size
        
        pos = start
        while pos <= end:
            window_end = min(pos + window_size - 1, end)
            
            ref_seq = self.get_sequence(chrom, pos, window_end)
            pers_seq = self.get_personal_sequence(chrom, pos, window_end)
            variants = self.get_variants_in_region(chrom, pos, window_end)
            
            yield SequenceWindow(
                chrom=chrom,
                start=pos,
                end=window_end,
                reference=ref_seq,
                personal=pers_seq,
                variants=variants
            )
            
            pos += stride
            
            # Break if we've covered the region
            if window_end >= end:
                break