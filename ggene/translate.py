#!/usr/bin/env python3
"""
Ribosome class for translating transcripts to proteins with variant awareness.
Handles proper start/stop codon detection and CDS traversal.
"""

from typing import Dict, List, Optional, Tuple, Callable, Union, TYPE_CHECKING, Any
from dataclasses import dataclass
import logging

from ggene import CODON_TABLE_DNA,CODON_TABLE_RNA,START_CODONS,ALTERNATIVE_START_CODONS,STOP_CODONS,to_rna, to_dna

if TYPE_CHECKING:
    from .features import Feature

logger = logging.getLogger(__name__)

# Standard genetic code
# CODON_TABLE = {
#     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
#     'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
#     'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
#     'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
#     'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
#     'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#     'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
#     'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#     'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#     'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
#     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#     'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#     'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
# }

# Special codons


# RNA versions
# RNA_CODON_TABLE = {k.replace('T', 'U'): v for k, v in CODON_TABLE.items()}
# RNA_START_CODONS = [c.replace('T', 'U') for c in START_CODONS]
# RNA_STOP_CODONS = [c.replace('T', 'U') for c in STOP_CODONS]


@dataclass
class TranslationResult:
    """Result of transcript translation."""
    protein_sequence: str
    start_position: int
    stop_position: int
    cds_segments: List[Tuple[int, int]]
    warnings: List[str]
    mutations: List[Dict[str, Any]]
    frame: int
    strand: str


class Ribosome:
    """
    Translates transcripts to proteins with variant awareness.
    Handles CDS traversal, start/stop codon detection, and frame shifts.
    """
    
    def __init__(self, allow_alternative_starts: bool = False,
                 verify_start_stop: bool = True,
                 handle_frameshifts: bool = True):
        """Initialize the Ribosome.
        
        Args:
            allow_alternative_starts: Whether to accept alternative start codons
            verify_start_stop: Whether to verify start/stop codons are present
            handle_frameshifts: Whether to handle frameshift mutations
        """
        self.allow_alternative_starts = allow_alternative_starts
        self.verify_start_stop = verify_start_stop
        self.handle_frameshifts = handle_frameshifts
        
        # Determine which start codons to use
        self.valid_start_codons = START_CODONS.copy()
        if allow_alternative_starts:
            self.valid_start_codons.extend(ALTERNATIVE_START_CODONS)
    
    def translate_transcript(self, transcript: 'Feature', 
                           sequence_generator: Callable,
                           personal: bool = True, verify_start_stop=True) -> TranslationResult:
        """Translate a transcript to protein sequence.
        
        Args:
            transcript: Transcript feature object
            sequence_generator: Function to get sequences (from GenomeManager)
            personal: Whether to use personal (True) or reference (False) sequence
            
        Returns:
            TranslationResult with protein sequence and metadata
        """
        warnings = []
        mutations = []
        self._vss = self.verify_start_stop
        self.verify_start_stop = verify_start_stop
        
        # Get all CDS segments for this transcript
        cds_segments = self._get_cds_segments(transcript)
        if not cds_segments:
            logger.warning(f"No CDS segments found for transcript")
            return TranslationResult(
                protein_sequence="",
                start_position=0,
                stop_position=0,
                cds_segments=[],
                warnings=["No CDS segments found"],
                mutations=[],
                frame=0,
                strand=transcript.strand if hasattr(transcript, 'strand') else '+'
            )
        
        # Sort CDS segments by position
        cds_segments.sort(key=lambda x: x.start)
        
        # Get strand information
        strand = transcript.strand if hasattr(transcript, 'strand') else '+'
        
        # Find the translation start
        start_pos, start_frame, start_warnings = self._find_translation_start(
            cds_segments, sequence_generator, strand, personal
        )
        warnings.extend(start_warnings)
        
        if start_pos is None:
            return TranslationResult(
                protein_sequence="",
                start_position=0,
                stop_position=0,
                cds_segments=[(c.start, c.end) for c in cds_segments],
                warnings=warnings,
                mutations=[],
                frame=0,
                strand=strand
            )
        
        # Build the complete CDS sequence starting from the start codon
        cds_sequence, actual_segments = self._build_cds_sequence(
            cds_segments, sequence_generator, start_pos, strand, personal
        )
        
        # Translate the CDS sequence
        protein_seq, stop_pos, translation_warnings = self._translate_cds(
            cds_sequence, start_pos, strand
        )
        warnings.extend(translation_warnings)
        
        # Check for mutations if personal sequence
        if personal:
            ref_result = self.translate_transcript(
                transcript, sequence_generator, personal=False
            )
            mutations = self._identify_mutations(protein_seq, ref_result.protein_sequence)
        
        
        self.verify_start_stop = self._vss
        
        return TranslationResult(
            protein_sequence=protein_seq,
            start_position=start_pos,
            stop_position=stop_pos if stop_pos else start_pos + len(cds_sequence),
            cds_segments=actual_segments,
            warnings=warnings,
            mutations=mutations,
            frame=start_frame,
            strand=strand
        )
    
    def _get_cds_segments(self, transcript: 'Feature') -> List['Feature']:
        """Extract CDS segments from a transcript.
        
        Args:
            transcript: Transcript feature
            
        Returns:
            List of CDS features sorted by position
        """
        cds_segments = []
        
        # Check direct subfeatures
        if hasattr(transcript, 'subfeatures'):
            for subfeature in transcript.subfeatures:
                if hasattr(subfeature, 'type') and subfeature.type.upper() == 'CDS':
                    cds_segments.append(subfeature)
        
        # Check for CDS in a different structure
        if hasattr(transcript, 'cds'):
            cds_segments.extend(transcript.cds)
        
        return cds_segments
    
    def _find_translation_start(self, cds_segments: List['Feature'],
                               sequence_generator: Callable,
                               strand: str,
                               personal: bool) -> Tuple[Optional[int], int, List[str]]:
        """Find the translation start position.
        
        Args:
            cds_segments: List of CDS features
            sequence_generator: Function to get sequences
            strand: Strand direction ('+' or '-')
            personal: Whether to use personal sequence
            
        Returns:
            Tuple of (start_position, frame, warnings)
        """
        warnings = []
        
        # Check if there's an annotated start_codon
        first_cds = cds_segments[0]
        
        # Look for start_codon feature
        if hasattr(first_cds, 'parent') and hasattr(first_cds.parent, 'subfeatures'):
            for feature in first_cds.parent.subfeatures:
                if hasattr(feature, 'type') and feature.type == 'start_codon':
                    return feature.start, 0, warnings
        
        # Otherwise, search for start codon in first CDS
        first_seq = sequence_generator(first_cds, personal=personal)
        if not first_seq:
            warnings.append("Could not get sequence for first CDS")
            return None, 0, warnings
        
        # If negative strand, work with reverse complement
        if strand == '-':
            first_seq = self._reverse_complement(first_seq)
        
        # Search for start codon
        for i in range(0, min(len(first_seq) - 2, 100), 3):  # Check first 100 bases
            codon = first_seq[i:i+3].upper()
            if codon in self.valid_start_codons:
                actual_pos = first_cds.start + i if strand == '+' else first_cds.end - i - 2
                return actual_pos, i % 3, warnings
        
        # If no start found but not verifying, use CDS start
        if not self.verify_start_stop:
            warnings.append("No start codon found, using CDS start")
            return first_cds.start, 0, warnings
        
        warnings.append("No valid start codon found")
        return None, 0, warnings
    
    def _build_cds_sequence(self, cds_segments: List['Feature'],
                          sequence_generator: Callable,
                          start_pos: int,
                          strand: str,
                          personal: bool) -> Tuple[str, List[Tuple[int, int]]]:
        """Build the complete CDS sequence from segments.
        
        Args:
            cds_segments: List of CDS features
            sequence_generator: Function to get sequences
            start_pos: Translation start position
            strand: Strand direction
            personal: Whether to use personal sequence
            
        Returns:
            Tuple of (complete_cds_sequence, segment_positions)
        """
        cds_parts = []
        segment_positions = []
        started = False
        
        for cds in cds_segments:
            # Check if we've reached the start position
            if not started:
                if strand == '+' and cds.end < start_pos:
                    continue
                elif strand == '-' and cds.start > start_pos:
                    continue
                else:
                    started = True
            
            # Get sequence for this CDS segment
            seq = sequence_generator(cds, personal=personal)
            if not seq:
                logger.warning(f"Could not get sequence for CDS segment {cds.start}-{cds.end}")
                continue
            
            # Trim if this contains the start position
            if strand == '+' and cds.start <= start_pos <= cds.end:
                trim_start = start_pos - cds.start
                seq = seq[trim_start:]
                segment_positions.append((start_pos, cds.end))
            elif strand == '-' and cds.start <= start_pos <= cds.end:
                trim_end = start_pos - cds.start + 1
                seq = seq[:trim_end]
                segment_positions.append((cds.start, start_pos))
            else:
                segment_positions.append((cds.start, cds.end))
            
            cds_parts.append(seq)
        
        # Concatenate all CDS parts
        complete_cds = ''.join(cds_parts)
        
        # Reverse complement if negative strand
        if strand == '-':
            complete_cds = self._reverse_complement(complete_cds)
        
        return complete_cds, segment_positions
    
    def _translate_cds(self, cds_sequence: str, 
                      start_pos: int,
                      strand: str) -> Tuple[str, Optional[int], List[str]]:
        """Translate CDS sequence to protein.
        
        Args:
            cds_sequence: Complete CDS DNA sequence
            start_pos: Start position in genome
            strand: Strand direction
            
        Returns:
            Tuple of (protein_sequence, stop_position, warnings)
        """
        warnings = []
        protein_parts = []
        stop_pos = None
        
        # Translate codon by codon
        for i in range(0, len(cds_sequence) - 2, 3):
            codon = cds_sequence[i:i+3].upper()
            
            # Handle unknown bases
            if 'N' in codon:
                protein_parts.append('X')
                continue
            
            # Check for stop codon
            if codon in STOP_CODONS:
                stop_pos = start_pos + i + 3 if strand == '+' else start_pos - i - 3
                if self.verify_start_stop:
                    break  # Stop at first stop codon
                else:
                    protein_parts.append('*')
            else:
                # Translate codon
                aa = CODON_TABLE_DNA.get(codon, 'X')
                protein_parts.append(aa)
        
        # Check if we found a stop codon
        if stop_pos is None and self.verify_start_stop:
            warnings.append("No stop codon found")
        
        return ''.join(protein_parts), stop_pos, warnings
    
    def _identify_mutations(self, personal_seq: str, reference_seq: str) -> List[Dict[str, Any]]:
        """Identify mutations between personal and reference sequences.
        
        Args:
            personal_seq: Personal protein sequence
            reference_seq: Reference protein sequence
            
        Returns:
            List of mutation descriptions
        """
        mutations = []
        
        min_len = min(len(personal_seq), len(reference_seq))
        
        for i in range(min_len):
            if personal_seq[i] != reference_seq[i]:
                mutation = {
                    'position': i + 1,  # 1-based
                    'reference': reference_seq[i],
                    'alternate': personal_seq[i],
                    'type': 'missense'
                }
                
                # Classify mutation type
                if personal_seq[i] == '*':
                    mutation['type'] = 'nonsense'
                elif reference_seq[i] == '*':
                    mutation['type'] = 'readthrough'
                
                mutations.append(mutation)
        
        # Check for length differences
        if len(personal_seq) != len(reference_seq):
            mutations.append({
                'position': min_len,
                'type': 'length_change',
                'reference_length': len(reference_seq),
                'alternate_length': len(personal_seq)
            })
        
        return mutations
    
    def _reverse_complement(self, dna_seq: str) -> str:
        """Get reverse complement of DNA sequence.
        
        Args:
            dna_seq: DNA sequence
            
        Returns:
            Reverse complement
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base.upper(), 'N') for base in reversed(dna_seq))
    
    def format_protein(self, protein_seq: str, line_length: int = 60) -> str:
        """Format protein sequence for display.
        
        Args:
            protein_seq: Protein sequence
            line_length: Characters per line
            
        Returns:
            Formatted protein sequence
        """
        lines = []
        for i in range(0, len(protein_seq), line_length):
            lines.append(protein_seq[i:i+line_length])
        return '\n'.join(lines)
    
    def get_mutation_summary(self, mutations: List[Dict[str, Any]]) -> str:
        """Generate a summary of mutations.
        
        Args:
            mutations: List of mutation dictionaries
            
        Returns:
            Human-readable mutation summary
        """
        if not mutations:
            return "No mutations detected"
        
        summary_parts = []
        
        for mut in mutations:
            if mut['type'] == 'missense':
                summary_parts.append(f"p.{mut['reference']}{mut['position']}{mut['alternate']}")
            elif mut['type'] == 'nonsense':
                summary_parts.append(f"p.{mut['reference']}{mut['position']}*")
            elif mut['type'] == 'readthrough':
                summary_parts.append(f"p.*{mut['position']}{mut['alternate']}")
            elif mut['type'] == 'length_change':
                diff = mut['alternate_length'] - mut['reference_length']
                if diff > 0:
                    summary_parts.append(f"p.{mut['position']}_ext{diff}")
                else:
                    summary_parts.append(f"p.{mut['position']}_del{-diff}")
        
        return ', '.join(summary_parts)