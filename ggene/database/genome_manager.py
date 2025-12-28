


from cyvcf2 import VCF
from cyvcf2 import Variant
import random
import numpy as np
import pysam
from typing import Dict, List, Optional, Tuple, Union, Iterator, Any, Callable
import os
import traceback
import sys

import json
import logging

logger = logging.getLogger(__name__)
logger.setLevel("CRITICAL")
# logger.setLevel(logging.DEBUG)

from ggene import draw
from ggene import dev
from ggene.config import get_config, get_paths
# from ggene.database.gene_map import GeneMap
from ggene.genome.features import Gene, Feature, shorten_variant
from ggene.genome.translate import Ribosome
from ggene.database.genome_iterator import UGenomeIterator
from ggene.browser.genome_browser_v2 import InteractiveGenomeBrowser
from ggene.database.annotations import UGenomeAnnotations
from ggene.database.ufeature import UFeature
# from ggene.database.annotations import UGenomeAnnotations, GTFStream, VCFStream, UFeature
from ggene.motifs import MotifDetector, PatternMotif, RepeatMotif
from ggene.seqs.lambdas import needs_features
from ggene.seqs.bio import reverse_complement, to_rna, to_dna, is_rna, is_dna, COMPLEMENT_MAP
from ggene.seqs.lambdas import lambda_map

class GenomeManager:
    """Main class for managing genomic data including VCF and GTF files."""
    
    chr_frm = "{chrom}"
    # chr_frm = "chr{chrom}"
    
    def __init__(self,
                 vcf_path: str = None,
                 gtf_path: str = None,
                 ref_path: str = None,
                 my_path: str = '',
                 library_path: str = None, 
                 **kwargs) -> None:
        """Initialize GenomeManager.

        Args:
            vcf_path: Path to VCF file
            gtf_path: Path to GTF file
            fasta_path: Path to FASTA file for genomic sequences
        """
        # Lazy import to avoid circular dependency
        # from ggene import get_paths
        DEFAULT_VCF_PATH, DEFAULT_GTF_PATH, DEFAULT_FASTA_PATH, DEFAULT_LIBRARY, other_paths = get_paths()

        # Use defaults if not provided
        if vcf_path is None:
            vcf_path = DEFAULT_VCF_PATH
        if gtf_path is None:
            gtf_path = DEFAULT_GTF_PATH
        if ref_path is None:
            ref_path = DEFAULT_FASTA_PATH
        if library_path is None:
            library_path = DEFAULT_LIBRARY

        try:
            self.vcf = VCF(vcf_path)
            
        except Exception as e:
            logger.error(f"Failed to load VCF (variants), skipping: {e}")
        
        try:
            self.library_path = library_path
            
            # Keep old GeneMap for backward compatibility
            # self.gene_map = GeneMap(gtf_path=gtf_path)
            
            # print(f"loading sequence from {DEFAULT_FASTA_PATH} (exists = {os.path.exists(DEFAULT_FASTA_PATH)})")
            
            # Initialize new unified annotation system with sequence streaming
            self.annotations = UGenomeAnnotations(
                fasta_path=ref_path if ref_path and os.path.exists(ref_path) else None,
                vcf_path=vcf_path if vcf_path and os.path.exists(vcf_path) else None
            )
            
            if gtf_path:
                self.annotations.add_genes(gtf_path)
            if vcf_path:
                self.annotations.add_variants(vcf_path)
            
            # Add GTF annotations
            # if gtf_path and os.path.exists(gtf_path):
            #     self.annotations.add_gtf(gtf_path, "genes")
            #     logger.info(f"Added GTF annotations: {gtf_path}")
            
            # Add VCF annotations separately for annotation track
            # if vcf_path and os.path.exists(vcf_path):
            #     self.annotations.add_vcf(vcf_path, "variants")
            #     logger.info(f"Added VCF annotations: {vcf_path}")
            
            # Initialize motif detector
            self.motif_detector = MotifDetector()
            self._setup_default_motifs()
            
            self.ribo = Ribosome()
            self.genes: Dict[str, 'Gene'] = {}
            
            # Initialize FASTA file if provided
            self.ref = None
            self.ref_path = ref_path
            if ref_path and os.path.exists(ref_path):
                try:
                    self.ref = pysam.FastaFile(ref_path)
                    logger.info(f"Loaded FASTA file: {ref_path}")
                except Exception as e:
                    logger.warning(f"Could not load FASTA file: {e}")
                    
            self.mine = None
            if my_path and os.path.exists(my_path):
                try:
                    self.mine = pysam.FastqFile(my_path)
                    logger.info(f"Loaded FASTQ file: {my_path}")
                except Exception as e:
                    logger.warning(f"Could not load FASTQ file: {e}")
            
            other_paths.update(kwargs)
            
            for k, v in other_paths.items():
                if k == "repeatmasker_path":
                    self.annotations.add_repeatmasker(v)
                if k == "dfam_path":
                    self.annotations.add_dfam(v)
                if k == "clinvar_path":
                    self.annotations.add_clinvar(v)
                    logger.debug(f"added clinvar from path {v}")
                
        except Exception as e:
            logger.error(f"Failed to initialize GenomeManager: {e}")
            raise
    
    def _setup_default_motifs(self):
        """Setup default motifs for detection."""
        
        # self.motif_detector.setup_default_motifs()
        # self.motif_detector.setup_default_motifs(pattern_classes = ["splice","promoter"], hmm_classes = ["Alu"])
        self.motif_detector.setup_default_motifs(pattern_classes = ["splice","promoter"], hmm_classes = [])
        
        # self.motif_detector.setup_default_motifs(class_names = ["splice","promoter"])
        # self.motif_detector.setup_default_motifs(class_names = ["hammerhead", "SRP", "pseudoknot","msat"])
        # self.motif_detector.setup_default_motifs(class_names = ["SRP_S"])

        # Add more motifs as needed
        logger.info("Initialized default motifs")
    
    # def add_annotation_source(self, name: str, filepath: str, source_type: str = "bed"):
    #     """Add a new annotation source to the unified system.
        
    #     Args:
    #         name: Name for this annotation source
    #         filepath: Path to the annotation file
    #         source_type: Type of file (bed, gff, vcf, etc.)
    #     """
    #     if source_type.lower() == "bed":
    #         from .annotations import BEDStream
    #         self.annotations.add_source(name, BEDStream(filepath))
    #     elif source_type.lower() in ["gff", "gtf"]:
    #         self.annotations.add_gtf(filepath, name)
    #     elif source_type.lower() == "vcf":
    #         self.annotations.add_vcf(filepath, name)
    #     else:
    #         raise ValueError(f"Unknown source type: {source_type}")
        
    #     logger.info(f"Added annotation source '{name}' from {filepath}")
    
    def get_all_annotations(self, chrom: str, start: int, end: int,
                           include_motifs: bool = True) -> List[UFeature]:
        """Get all annotations for a region including motifs.
        
        Args:
            chrom: Chromosome
            start: Start position (1-based)
            end: End position (1-based)
            include_motifs: Whether to scan for motifs
            
        Returns:
            List of UFeature objects
        """
        # Get annotations from databases
        annotations = self.annotations.query_range(chrom, start, end)
        
        # Scan for motifs if requested
        if include_motifs:
            seq = self.get_sequence(chrom, start, end)
            if seq:
                motif_features = self.scan_motifs(seq, chrom, start)
                annotations.extend(motif_features)
        
        return sorted(annotations)
    
    def scan_motifs(self, sequence: str, chrom: str, start_pos: int, strand = '+') -> List[UFeature]:
        """Scan a sequence for motifs.
        
        Args:
            sequence: DNA sequence to scan
            chrom: Chromosome name
            start_pos: Genomic start position of the sequence
            
        Returns:
            List of UFeature objects for found motifs
        """
        features = []
        
        all_insts = self.motif_detector.identify(sequence)
        
        for motif_name, instances in all_insts.items():
            if instances:
                # for motif_start, motif_end, score, is_rc in instances:
                for motif in instances:
                    
                    motif_start = motif.get("start", 0)
                    motif_end = motif.get("end", 0)
                    score = motif.get("score", 0)
                    is_rc = motif.get("is_rc", False)
                    mtf_cls = motif.get("class", "")
                    
                    mseq = sequence[motif_start:motif_end]
                    mtfstrand = strand
                    if is_rc:
                        mseq = reverse_complement(mseq)
                        mtfstrand = '-' if mtfstrand == '+' else '+'
                    
                    features.append(UFeature(
                        chrom=chrom,
                        start=start_pos + motif_start,
                        end=start_pos + motif_end - 1,
                        feature_type="motif",
                        source="MotifDetector",
                        score=score,
                        strand=mtfstrand,
                        name=motif_name,
                        attributes={
                            'sequence': mseq,
                            'is_rc':is_rc,
                            'motif_class':mtf_cls,
                            'caller':'GenomeManager'
                        }
                    ))
        
        return features
    
    def find_and_assemble_genes(self, genes):
        outgenes = []
        for g in genes:
            chrom, gs = self.gene_map.find_gene(g)
            if not chrom:
                logger.debug(f"failed to find chromosome for gene {g}")
                continue
            
            newgene = self.assemble_gene(g, chrom)
            outgenes.append(newgene)
            
        return outgenes
    
    def assemble_gene(self, gene_name: str, chrom: Union[str, int], 
                     min_qual: float = 5) -> Optional['Gene']:
        """Assemble a gene with its variants.
        
        Args:
            gene_name: Name of the gene
            chrom: Chromosome number or string
            min_qual: Minimum quality threshold for variants
            
        Returns:
            Gene object or None if gene not found
        """
        try:
            gene_data = self.gene_map.get_gene(gene_name, chrom)
            _gene = gene_data['gene'][0]
            
        except Exception as e:
            logger.warning(f"Could not create gene {gene_name}: {e}")
            return None
        
        try:
            idx = self._make_index(chrom, _gene['start'], _gene['end'])
            # variants = self.vcf(idx)
            variants = self.annotations.streams["variants"].vcf(idx)
            print(idx)
            filtered_variants = [v for v in variants if v.QUAL > min_qual]
            gene_data['variant'] = filtered_variants
        except Exception as e:
            logger.error(f"Failed to add variants to gene {gene_name}: {e}")
        
        new_gene = Gene(gene_data)
        self.save_gene(new_gene)
        self.genes[gene_name] = new_gene
        
        return new_gene
    
    def find_upstream_motif(self, feature : Feature, 
                           motif: str, max_distance: int = 1000) -> Optional[Dict[str, Any]]:
        """Find a motif upstream of this feature.
        
        Args:
            genome_manager: GenomeManager instance with loaded FASTA
            motif: DNA motif to search for
            max_distance: Maximum distance to search upstream
            
        Returns:
            Dictionary with position and sequence if found, None otherwise
        """
        if not hasattr(feature, 'chrom') or not hasattr(feature, 'start'):
            return None
            
        # Simple motif search without regex for basic implementation
        predicate = lambda seq: motif.upper() in seq.upper()
        
        result = self.probe_sequence_backward(
            feature.chrom, feature.start, predicate, 
            window_size=len(motif), max_distance=max_distance
        )
        
        if result:
            pos, seq = result
            return {
                'position': pos,
                'sequence': seq,
                'distance': self.start - pos - len(motif),
                'motif': motif
            }
        return None
    
    def find_downstream_motif(self, feature: Feature,
                             motif: str, max_distance: int = 1000) -> Optional[Dict[str, Any]]:
        """Find a motif downstream of this feature.
        
        Args:
            genome_manager: GenomeManager instance with loaded FASTA
            motif: DNA motif to search for
            max_distance: Maximum distance to search downstream
            
        Returns:
            Dictionary with position and sequence if found, None otherwise
        """
        if not hasattr(feature, 'chrom') or not hasattr(feature, 'end'):
            return None
            
        # Simple motif search without regex for basic implementation
        predicate = lambda seq: motif.upper() in seq.upper()
        
        result = self.probe_sequence_forward(
            feature.chrom, feature.end + 1, predicate,
            window_size=len(motif), max_distance=max_distance
        )
        
        if result:
            pos, seq = result
            return {
                'position': pos,
                'sequence': seq,
                'distance': pos - feature.end,
                'motif': motif
            }
        return None
    
    def transcribe_gene(self, gene_name: str, transcript_index: int = 0,
                       personal: bool = True, verify_start_stop=True) -> Optional[Dict[str, Any]]:
        """Transcribe a gene to protein using the Ribosome.
        
        Args:
            gene_name: Name of the gene to transcribe
            transcript_index: Which transcript to use (default: 0)
            personal: Whether to use personal (True) or reference (False) sequence
            
        Returns:
            Dictionary with protein sequence and translation details, or None if failed
        """
        try:
            gene = self.load_gene(gene_name)
            if not gene:
                logger.error(f"Could not load gene {gene_name}")
                return None
            
            # Get transcripts
            if hasattr(gene, 'transcripts'):
                transcripts = list(gene.transcripts.values())
            else:
                logger.error(f"Gene {gene_name} has no transcripts")
                return None
            
            if transcript_index >= len(transcripts):
                logger.error(f"Transcript index {transcript_index} out of range")
                return None
            
            transcript = transcripts[transcript_index]
            
            # Use Ribosome to translate
            result = self.ribo.translate_transcript(
                transcript,
                self.get_feature_sequence,
                personal=personal
            )
            
            # Return formatted result
            return {
                'gene_name': gene_name,
                'transcript_index': transcript_index,
                'result': result,
                'mutation_summary': self.ribo.get_mutation_summary(result.mutations) if result.mutations else None,
                'personal': personal
            }
            
        except Exception as e:
            logger.error(f"Failed to transcribe gene {gene_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def compare_transcribed(self, gene_name: str, transcript_index: int = 0) -> Optional[Dict[str, Any]]:
        """Compare personal and reference protein sequences for a gene.
        
        Args:
            gene_name: Name of the gene to compare
            transcript_index: Which transcript to use (default: 0)
            
        Returns:
            Dictionary with comparison results, or None if failed
        """
        # Get both personal and reference translations
        personal_result = self.transcribe_gene(gene_name, transcript_index, personal=True)
        reference_result = self.transcribe_gene(gene_name, transcript_index, personal=False)
        
        if not personal_result or not reference_result:
            return None
        
        # Compare the sequences
        personal_seq = personal_result['result'].protein_sequence
        reference_seq = reference_result['result'].protein_sequence
        
        # Find differences
        differences = []
        min_len = min(len(personal_seq), len(reference_seq))
        
        for i in range(min_len):
            if personal_seq[i] != reference_seq[i]:
                differences.append({
                    'position': i + 1,
                    'reference': reference_seq[i],
                    'personal': personal_seq[i],
                    'notation': f"p.{reference_seq[i]}{i+1}{personal_seq[i]}"
                })
        
        # Check for length differences
        length_diff = len(personal_seq) - len(reference_seq)
        
        return {
            'gene_name': gene_name,
            'transcript_index': transcript_index,
            'reference_sequence': reference_seq,
            'personal_sequence': personal_seq,
            'reference_length': len(reference_seq),
            'personal_length': len(personal_seq),
            'length_difference': length_diff,
            'amino_acid_changes': differences,
            'num_changes': len(differences),
            'percent_identity': (min_len - len(differences)) / min_len * 100 if min_len > 0 else 100,
            'mutations': personal_result.get('mutations', []),
            'mutation_summary': personal_result.get('mutation_summary', ''),
            'warnings': personal_result.get('warnings', [])
        }
    
    def _make_index(self, chromosome_number: Union[str, int], 
                   start: Optional[int] = None, end: Optional[int] = None) -> str:
        """Create chromosome index string for VCF queries.
        
        Args:
            chromosome_number: Chromosome identifier
            start: Start position (optional)
            end: End position (optional)
            
        Returns:
            Formatted index string
        """
        if start is None or end is None:
            return f'chr{chromosome_number}'
        return f'chr{chromosome_number}:{start}-{end}'
    
    def make_variant_map(self, chrom: Union[str, int], 
                        start: Optional[int] = None, 
                        end: Optional[int] = None) -> Dict[int, Dict[int, int]]:
        """Create a map of variant types by reference and alternate lengths.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (optional)
            end: End position (optional)
            
        Returns:
            Nested dictionary mapping ref_len -> alt_len -> count
        """
        var_map: Dict[int, Dict[int, int]] = {}
        
        try:
            for var in self.vcf(self._make_index(chrom, start=start, end=end)):
                n_ref = len(var.REF) - 1
                n_alt = len(var.ALT[0]) - 1 if var.ALT else 0
                
                if n_ref not in var_map:
                    var_map[n_ref] = {}
                    
                if n_alt not in var_map[n_ref]:
                    var_map[n_ref][n_alt] = 0
                    
                var_map[n_ref][n_alt] += 1
        except Exception as e:
            logger.error(f"Failed to create variant map: {e}")
            
        return var_map

    def get_quality_stats(self, chrom: Union[str, int], 
                         start: Optional[int] = None, 
                         end: Optional[int] = None) -> Tuple[Any, List[float]]:
        """Get quality statistics for variants in a region.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (optional)
            end: End position (optional)
            
        Returns:
            Tuple of (statistics object, quality values list)
        """
        try:
            quals = [v.QUAL for v in self.vcf(self._make_index(chrom, start=start, end=end))]
            stats = dev.stat.from_data(quals)
            return stats, quals
        except Exception as e:
            logger.error(f"Failed to get quality stats: {e}")
            return None, [] 

    def find_long_variants(self, chrom: Union[str, int], min_length: int, 
                          start: Optional[int] = None, 
                          end: Optional[int] = None) -> Tuple[List[Dict], List[Dict]]:
        """Find long insertion and deletion variants.
        
        Args:
            chrom: Chromosome identifier
            min_length: Minimum variant length to include
            start: Start position (optional)
            end: End position (optional)
            
        Returns:
            Tuple of (deletions list, insertions list)
        """
        long_dels = []
        long_inserts = []
        
        try:
            for var in self.vcf(self._make_index(chrom, start=start, end=end)):
                n_ref = len(var.REF) - 1
                n_alt = len(var.ALT[0]) - 1 if var.ALT else 0
                var_len = n_ref - n_alt
                
                if np.abs(var_len) > min_length:
                    features = list(self.gene_map.fetch(chrom, var.POS, var.POS + 1, features=tuple()))
                    
                    var_short = shorten_variant(var)
                    var_short['features'] = features
                    
                    if var_len > 0:
                        long_dels.append(var_short)
                    else:
                        long_inserts.append(var_short)
        except Exception as e:
            logger.error(f"Failed to find long variants: {e}")
            
        return long_dels, long_inserts
    
    @classmethod
    def get_random_location(cls, chromes = [], margin = 1e6):
        
        if not chromes:
            chromes = [chrom for chrom in cls.iter_chromes()]
        
        chrome = random.choice(chromes)
        
        from ggene.database.annotations import chr_lens
        
        max_index = chr_lens.get(chrome, 10e6)
        
        rand_pos = int(random.random()*(max_index - 2*margin) + margin)
        
        return chrome, rand_pos
    
    def get_sequence(self, chrom: Union[str, int], start: int, end: int, frame = 0) -> Optional[str]:
        """Get genomic sequence for a region.
        
        Now uses unified streaming system for sequence access.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            
        Returns:
            Sequence string or None if not available
        """
        # Use unified streaming system if available
        if self.annotations and self.annotations.sequence_stream:
            chrom_str = str(chrom)
            return self.annotations.get_sequence(chrom_str, start+frame, end+frame)
        
        # Fallback to direct FASTA access if available
        if self.ref:
            try:
                # Convert to string and ensure proper format
                chrom_str = self.chr_frm.format(chrom=chrom)
                    
                # pysam uses 0-based coordinates
                sequence = self.ref.fetch(chrom_str, start - 1, end)
                return sequence
            except Exception as e:
                logger.error(f"Failed to fetch sequence: {e}")
                return None
        
        logger.warning("No sequence source available")
        return None
    
    def get_feature_sequence(self, feature: Feature, 
                           upstream: int = 0, 
                           downstream: int = 0,
                           personal: bool = True,
                           as_rna = False) -> Optional[str]:
        """Get sequence for a feature with optional flanking regions.
        
        Now uses unified streaming system with integrated variant application.
        
        Args:
            feature: Feature object
            upstream: Bases to include upstream (5')
            downstream: Bases to include downstream (3')
            personal: If True, apply personal variants to get personalized sequence
            
        Returns:
            Sequence string or None if not available
        """
        if not hasattr(feature, 'chrom') or not hasattr(feature, 'start') or not hasattr(feature, 'end'):
            logger.warning("Feature missing position information")
            return None
            
        start = max(1, feature.start - upstream)
        end = feature.end + downstream
        chrom_str = str(feature.chrom)
        
        # Use unified streaming system if available
        if self.annotations and self.annotations.sequence_stream:
            if personal:
                seq = self.annotations.get_personal_sequence(chrom_str, start, end)
            else:
                seq = self.annotations.get_sequence(chrom_str, start, end)
        else:
            # Fallback to old method
            reference_seq = self.get_sequence(feature.chrom, start, end)
            if not reference_seq:
                return None
                
            if personal:
                seq = self._apply_variants_to_sequence(reference_seq, feature, start, end, upstream, downstream)
            else:
                seq = reference_seq
        
        if not seq:
            return None
            
        # Apply strand and RNA conversion if needed
        if feature.strand == '-':
            seq = reverse_complement(seq)
            
        if as_rna:
            seq = to_rna(seq)
        
        return seq
    
    def convolve(self, seq1, seq2):
        
        seq_len = min(len(seq1), len(seq2))
        seq1, seq2 = seq1[:seq_len], seq2[seq_len]
        
        start = seq_len // 2
        ips = []
        comp_ips = []
        
        for t in range(-start, start):
            
            sslen = seq_len - abs(t)
            seq1t = seq1[max(t, 0):max(t, 0) + sslen]
            seq2t = seq2[max(-t, 0):max(-t, 0) + sslen]
            summ = 0
            csumm = 0
            
            for sa, sb in zip(seq1t, seq2t):
                
                csb = COMPLEMENT_MAP.get(sb)
                
                if sa == sb:
                    summ +=1
                elif sa == csb:
                    csumm += 1
                
            ips.append(summ)
            comp_ips.append(csumm)
        
        return summ, csumm
    
    def probe_sequence_forward(self, chrom: Union[str, int], start: int, 
                              predicate: Callable[[str], bool], 
                              window_size: int = 1, 
                              max_distance: int = 10000,
                              step_size: int = 1) -> Optional[Tuple[int, str]]:
        """Search forward from a position for sequence matching predicate.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (1-based)
            predicate: Function that returns True when desired sequence is found
            window_size: Size of sequence window to test
            max_distance: Maximum distance to search
            step_size: Number of bases to move forward each iteration
            
        Returns:
            Tuple of (position, matching_sequence) or None if not found
        """
        # Check if we have a sequence source
        if not (self.annotations and self.annotations.sequence_stream) and not self.ref:
            logger.warning("No sequence source available")
            return None
            
        try:
            for pos in range(start, start + max_distance, step_size):
                seq = self.get_sequence(chrom, pos, pos + window_size - 1)
                if seq and predicate(seq):
                    return (pos, seq)
            return None
        except Exception as e:
            logger.error(f"Error in forward probe: {e}")
            return None
    
    def probe_sequence_backward(self, chrom: Union[str, int], start: int,
                               predicate: Callable[[str], bool],
                               window_size: int = 1,
                               max_distance: int = 10000,
                               step_size: int = 1) -> Optional[Tuple[int, str]]:
        """Search backward from a position for sequence matching predicate.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (1-based)
            predicate: Function that returns True when desired sequence is found
            window_size: Size of sequence window to test
            max_distance: Maximum distance to search
            step_size: Number of bases to move backward each iteration
            
        Returns:
            Tuple of (position, matching_sequence) or None if not found
        """
        # Check if we have a sequence source
        if not (self.annotations and self.annotations.sequence_stream) and not self.ref:
            logger.warning("No sequence source available")
            return None
            
        try:
            for pos in range(start, max(1, start - max_distance), -step_size):
                seq = self.get_sequence(chrom, pos - window_size + 1, pos)
                if seq and predicate(seq):
                    return (pos - window_size + 1, seq)
            return None
        except Exception as e:
            logger.error(f"Error in backward probe: {e}")
            return None
    
    def find_motif_around_feature(self, feature: Feature, motif: str,
                                 upstream_search: int = 1000,
                                 downstream_search: int = 1000) -> List[Dict[str, Any]]:
        """Find all occurrences of a motif around a feature.
        
        Args:
            feature: Feature to search around
            motif: DNA motif to search for (supports IUPAC codes)
            upstream_search: Distance to search upstream
            downstream_search: Distance to search downstream
            
        Returns:
            List of dictionaries with position and match information
        """
        if not hasattr(feature, 'chrom') or not hasattr(feature, 'start'):
            return []
            
        results = []
        
        # Convert IUPAC codes to regex
        iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
            'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
        }
        
        import re
        pattern = motif
        for code, bases in iupac_map.items():
            pattern = pattern.replace(code, bases)
        
        # Get sequence around feature
        seq = self.get_feature_sequence(feature, upstream_search, downstream_search)
        if not seq:
            return []
            
        # Find all matches
        start_pos = feature.start - upstream_search
        for match in re.finditer(pattern, seq, re.IGNORECASE):
            results.append({
                'position': start_pos + match.start(),
                'sequence': match.group(),
                'relative_to_feature': match.start() - upstream_search,
                'in_feature': (match.start() >= upstream_search and 
                              match.start() < upstream_search + len(feature))
            })
            
        return results
    
    def get_gc_content(self, chrom: Union[str, int], start: int, end: int) -> Optional[float]:
        """Calculate GC content for a region.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            
        Returns:
            GC content as fraction (0-1) or None if not available
        """
        seq = self.get_sequence(chrom, start, end)
        if not seq:
            return None
            
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        total = len(seq)
        
        return gc_count / total if total > 0 else 0.0
    
    def scan_gc_content(self, chrom: Union[str, int], start: int, end: int,
                       window_size: int = 100, step_size: int = 50) -> List[Tuple[int, float]]:
        """Scan GC content across a region.
        
        Args:
            chrom: Chromosome identifier
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            window_size: Size of sliding window
            step_size: Step size for sliding window
            
        Returns:
            List of (position, gc_content) tuples
        """
        results = []
        
        for pos in range(start, end - window_size + 1, step_size):
            gc = self.get_gc_content(chrom, pos, pos + window_size - 1)
            if gc is not None:
                results.append((pos, gc))
                
        return results
    
    def _apply_variants_to_sequence(self, sequence: str, feature: Feature,
                                   seq_start: int, seq_end: int,
                                   upstream: int, downstream: int) -> str:
        """Apply personal variants to a reference sequence.
        
        Args:
            sequence: Reference sequence
            feature: Feature object (may contain variants)
            seq_start: Genomic start position of the sequence
            seq_end: Genomic end position of the sequence
            upstream: Upstream bases included
            downstream: Downstream bases included
            
        Returns:
            Personalized sequence with variants applied
        """
        # Collect all variants in the region
        variants = self._collect_variants_in_region(feature, seq_start, seq_end)
        
        if not variants:
            return sequence
            
        # Sort variants by position
        variants.sort(key=lambda v: getattr(v, 'pos', getattr(v, 'start', 0)))
        
        # Convert sequence to list for easier manipulation
        seq_list = list(sequence)
        
        # Track position offset due to indels
        offset = 0
        
        for variant in variants:
            # Get variant information
            var_pos = getattr(variant, 'pos', getattr(variant, 'start', 0))
            ref_allele = getattr(variant, 'ref', '')
            alt_allele = getattr(variant, 'alt', '')
            
            if not ref_allele or alt_allele is None:
                logger.debug(f"Skipping variant at {var_pos}: missing ref/alt")
                continue
                
            # Calculate position in sequence (0-based)
            seq_pos = var_pos - seq_start + offset
            
            # Validate position
            if seq_pos < 0 or seq_pos >= len(seq_list):
                # logger.debug(f"Variant at {var_pos} outside sequence bounds")
                continue
                
            # Apply variant
            try:
                if len(ref_allele) == len(alt_allele):
                    # SNP - simple substitution
                    for i in range(len(ref_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                            
                elif len(ref_allele) > len(alt_allele):
                    # Deletion
                    del_length = len(ref_allele) - len(alt_allele)
                    # Replace ref with alt
                    for i in range(len(alt_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                    # Delete extra bases
                    for _ in range(del_length):
                        if seq_pos + len(alt_allele) < len(seq_list):
                            del seq_list[seq_pos + len(alt_allele)]
                    offset -= del_length
                    
                else:
                    # Insertion
                    ins_length = len(alt_allele) - len(ref_allele)
                    # Replace ref bases
                    for i in range(len(ref_allele)):
                        if seq_pos + i < len(seq_list):
                            seq_list[seq_pos + i] = alt_allele[i]
                    # Insert additional bases
                    for i in range(len(ref_allele), len(alt_allele)):
                        seq_list.insert(seq_pos + i, alt_allele[i])
                    offset += ins_length
                    
                logger.debug(f"Applied variant at {var_pos}: {ref_allele} -> {alt_allele}")
                
            except Exception as e:
                logger.warning(f"Failed to apply variant at {var_pos}: {e}")
                
        return ''.join(seq_list)
    
    def _collect_variants_in_region(self, feature: Feature, 
                                   start: int, end: int) -> List[Feature]:
        """Collect all variants within a genomic region.
        
        Args:
            feature: Feature that may contain variants
            start: Genomic start position
            end: Genomic end position
            
        Returns:
            List of variant features in the region
        """
        variants = []
        
        if hasattr(feature, 'parent_gene'):
            variants = feature.parent_gene.variants
        else:
            # ftype = feature.type if hasattr(feature,'type') else feature.
            
            # Check if this feature is itself a variant
            try:
                if feature.type == 'variant':
                    var_pos = getattr(feature, 'pos', getattr(feature, 'start', 0))
                    if start <= var_pos <= end:
                        variants.append(feature)
            except Exception as e:
                print(f"Exception trying to access feature {feature} type: {str(e)}")
                print(traceback.format_exc())
                raise
            # Recursively collect variants from subfeatures
            for subfeature in feature.subfeatures:
                if subfeature.type == 'variant':
                    var_pos = getattr(subfeature, 'pos', getattr(subfeature, 'start', 0))
                    if start <= var_pos <= end:
                        variants.append(subfeature)
                elif hasattr(subfeature, 'subfeatures'):
                    # Recursively check subfeatures
                    variants.extend(self._collect_variants_in_region(subfeature, start, end))
            
            # If this is a gene, also check its components
            if hasattr(feature, 'variants'):
                for variant in feature.variants:
                    var_pos = getattr(variant, 'pos', getattr(variant, 'start', 0))
                    if start <= var_pos <= end:
                        variants.append(variant)
                    
        # Remove duplicates while preserving order
        seen = set()
        unique_variants = []
        for variant in variants:
            variant_id = id(variant)
            if variant_id not in seen:
                seen.add(variant_id)
                unique_variants.append(variant)
                
        return unique_variants
    
    def get_aligned_sequences(self, chrom: Union[str, int], start: int, end: int) -> Tuple[str, str, List[Tuple[int, int]]]:
        """Get aligned reference and personal sequences with gaps for visualization.
        
        Uses unified streaming system for gap-aligned sequences needed by genome browser.
        
        Args:
            chrom: Chromosome
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            
        Returns:
            Tuple of (aligned_ref, aligned_pers, variant_list)
            where variant_list contains tuples of (position, delta)
        """
        chrom_str = str(chrom)
        
        if self.annotations and self.annotations.sequence_stream:
            return self.annotations.get_aligned_sequences(chrom_str, start, end)
        
        # Fallback to simple sequences without alignment
        ref = self.get_sequence(chrom, start, end)
        pers = ref  # No variants available
        return ref or "", pers or "", []
    
    def get_sequence_comparison(self, feature: Feature,
                               upstream: int = 0,
                               downstream: int = 0) -> Dict[str, Any]:
        """Get both reference and personal sequences with variant annotations.
        
        Args:
            feature: Feature object
            upstream: Bases to include upstream (5')
            downstream: Bases to include downstream (3')
            
        Returns:
            Dictionary with reference sequence, personal sequence, and variant info
        """
        ref_seq = self.get_feature_sequence(feature, upstream, downstream, personal=False)
        personal_seq = self.get_feature_sequence(feature, upstream, downstream, personal=True)
        
        if not ref_seq:
            return None
            
        # Collect variant information
        start = max(1, feature.start - upstream)
        end = feature.end + downstream
        variants = self._collect_variants_in_region(feature, start, end)
        
        variant_info = []
        for variant in variants:
            var_pos = getattr(variant, 'pos', getattr(variant, 'start', 0))
            variant_info.append({
                'position': var_pos,
                'ref': getattr(variant, 'ref', ''),
                'alt': getattr(variant, 'alt', ''),
                'type': getattr(variant, 'var_type', 'unknown'),
                'quality': getattr(variant, 'qual', 0),
                'relative_position': var_pos - feature.start
            })
        
        # Find differences between sequences
        differences = []
        if ref_seq and personal_seq and ref_seq != personal_seq:
            # Simple character-by-character comparison for SNPs
            min_len = min(len(ref_seq), len(personal_seq))
            for i in range(min_len):
                if ref_seq[i] != personal_seq[i]:
                    differences.append({
                        'position': i,
                        'ref_base': ref_seq[i],
                        'personal_base': personal_seq[i],
                        'genomic_position': start + i
                    })
        
        return {
            'reference_sequence': ref_seq,
            'personal_sequence': personal_seq,
            'variants': variant_info,
            'differences': differences,
            'num_variants': len(variants),
            'sequences_differ': ref_seq != personal_seq
        }
    
    def get_all_feature_sequences(self, gene_or_features: Union['Gene', List[Feature]],
                                 personal: bool = False) -> Dict[str, str]:
        """Get sequences for multiple features.
        
        Args:
            gene_or_features: Either a Gene object or list of Feature objects
            personal: If True, return personal sequences with variants
            
        Returns:
            Dictionary mapping feature IDs to sequences
        """
        sequences = {}
        
        # If it's a gene, get all its features
        if hasattr(gene_or_features, 'transcripts'):
            # It's a gene
            gene = gene_or_features
            features = []
            
            # Add exons
            for exon in gene.exons.values():
                features.append(exon)
                
            # Add CDS
            features.extend(gene.cds)
            
            # Add other subfeatures
            for subfeature in gene.subfeatures:
                if subfeature.type not in ['variant', 'intron']:
                    features.append(subfeature)
        else:
            # It's a list of features
            features = gene_or_features if isinstance(gene_or_features, list) else [gene_or_features]
        
        # Get sequences for each feature
        for feature in features:
            seq = self.get_feature_sequence(feature, personal=personal)
            if seq:
                feature_id = getattr(feature, 'sfid', f"{feature.type}_{feature.start}")
                sequences[feature_id] = seq
                
        return sequences
    
    def export_personalized_fasta(self, features: List[Feature], 
                                 output_file: str,
                                 include_reference: bool = True) -> None:
        """Export features to FASTA format with personal variants applied.
        
        Args:
            features: List of features to export
            output_file: Output FASTA file path
            include_reference: If True, also include reference sequences
        """
        try:
            with open(output_file, 'w') as f:
                for feature in features:
                    # Get feature ID
                    feature_id = getattr(feature, 'sfid', 
                                       f"{feature.type}_{feature.chrom}_{feature.start}_{feature.end}")
                    
                    # Write personal sequence
                    personal_seq = self.get_feature_sequence(feature, personal=True)
                    if personal_seq:
                        f.write(f">{feature_id}_personal\n")
                        # Write sequence in 80-character lines
                        for i in range(0, len(personal_seq), 80):
                            f.write(personal_seq[i:i+80] + "\n")
                    
                    # Write reference sequence if requested
                    if include_reference:
                        ref_seq = self.get_feature_sequence(feature, personal=False)
                        if ref_seq:
                            f.write(f">{feature_id}_reference\n")
                            for i in range(0, len(ref_seq), 80):
                                f.write(ref_seq[i:i+80] + "\n")
                                
            logger.info(f"Exported {len(features)} features to {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to export FASTA: {e}")
            raise
    
    def save_gene(self, gene: Gene,
                  include_ordered_features: bool = False) -> None:
        """Save a gene to a JSON file.
        
        Args:
            gene: Gene object to save
            output_file: Path to output JSON file
            include_ordered_features: Whether to include ordered feature lists
        """
        import json
        
        try:
            gene_dict = gene.to_dict()
            gene_dict_ser = dev.make_serializable(gene_dict)
            
            output_file = os.path.join(self.library_path, gene.name + '.json')
            # Wrap in gene name for consistency with your format
            # data = {gene.name: gene_dict}
            
            with open(output_file, 'w') as f:
                json.dump(gene_dict_ser, f, indent=3)
                
            logger.info(f"Saved gene {gene.name} to {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to save gene: {e}")
            raise
    
    def load_gene(self, gene_name: str) -> Optional[Gene]:
        """Load a gene from a JSON file.
        
        Args:
            json_file: Path to JSON file containing gene data
            
        Returns:
            Gene object or None if loading fails
        """
        
        gene_name = gene_name.upper()
        json_path = os.path.join(self.library_path, gene_name + '.json')
        
        if not os.path.exists(json_path):
            self.find_and_assemble_genes([gene_name])
        
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Extract gene data (handle both wrapped and unwrapped format)
            if len(data) == 1 and isinstance(list(data.values())[0], dict):
                # Wrapped format: {"GENE_NAME": {...}}
                gene_name = list(data.keys())[0]
                gene_dict = data[gene_name]
            else:
                # Direct format: {...}
                gene_dict = data
            
            # Reconstruct the gene
            try:
                gene = Gene.from_dict(gene_dict)
            except Exception as e:
                print(f"exception construction gene from dict: {str(e)}")
            # Add to our gene cache
            if hasattr(gene, 'name'):
                self.genes[gene.name] = gene
                
            logger.info(f"Loaded gene {getattr(gene, 'name', 'unknown')} from {json_path}")
            return gene
            
        except Exception as e:
            logger.error(f"Failed to load gene from {json_path}: {e}")
            
    
    def find_gene(self, name):
        
        avail_genes = os.listdir(self.library_path)
        for a in avail_genes:
            gn = a.rstrip('.json')
            if gn == name.upper():
                return a
        return None
    
    def load_gene_library(self,  pattern: str = "*.json") -> Dict[str, Gene]:
        """Load multiple genes from JSON files in a directory.
        
        Args:
            directory: Directory containing JSON files
            pattern: Glob pattern for JSON files
            
        Returns:
            Dictionary mapping gene names to Gene objects
        """
        import glob
        import os
        
        loaded_genes = {}
        json_files = glob.glob(os.path.join(self.library_path, pattern))
        
        for json_file in json_files:
            try:
                gene = self.load_gene(json_file)
                if gene and hasattr(gene, 'name'):
                    loaded_genes[gene.name] = gene
                    logger.debug(f"Loaded {gene.name} from {json_file}")
            except Exception as e:
                logger.warning(f"Failed to load {json_file}: {e}")
                continue
                
        logger.info(f"Loaded {len(loaded_genes)} genes from {self.library_path}")
        return loaded_genes
    
    def get_gene(self, gene_name: str, chrom: Union[str, int] = None,
                        min_qual: float = 0) -> Optional[Gene]:
        """Get a gene from cache, JSON file, or assemble it.
        
        Args:
            gene_name: Name of the gene
            chrom: Chromosome (required if assembling)
            json_file: Specific JSON file to load from
            json_dir: Directory to look for JSON files
            min_qual: Minimum quality for variants (if assembling)
            
        Returns:
            Gene object or None
        """
        # Check cache first
        if gene_name in self.genes:
            logger.debug(f"Found {gene_name} in cache")
            return self.genes[gene_name]
        
        # Try to load from JSON
        gn = self.find_gene(gene_name)
        if gn:
            gene = self.load_gene(gene_name)
            if gene:
                return gene

        # If we have chromosome info, assemble the gene
        if chrom is not None:
            logger.info(f"Assembling {gene_name} from genomic data")
            return self.assemble_gene(gene_name, chrom, min_qual)
        
        logger.warning(f"Could not find or assemble gene {gene_name}")
        return None
    
    def __getitem__(self, ind):
        self.load_gene(ind.upper())
    
    def iterate_genome(self, chrom: Union[str, int], start: int, 
                      end: Optional[int] = None, **kwargs) -> UGenomeIterator:
        """Create a genome iterator for a specific region.
        
        Args:
            chrom: Chromosome to iterate over
            start: Starting position (1-based)
            end: Ending position (1-based, inclusive)
            **kwargs: Additional arguments for GenomeIterator
            
        Returns:
            GenomeIterator instance
        """
        # return GenomeIterator(self, chrom, start, end, **kwargs)
        return UGenomeIterator(self, chrom, start, end, **kwargs)
    
    def extract_ml_features(self, features: List[Feature], 
                           upstream: int = 1000, 
                           downstream: int = 1000, **kwargs) -> List[Dict[str, Any]]:
        """Extract sequences around features for ML analysis.
        
        Args:
            features: List of features to extract sequences for
            upstream: Bases to include upstream
            downstream: Bases to include downstream
            **kwargs: Additional arguments for FeatureExtractor
            
        Returns:
            List of dictionaries with extracted sequences and metadata
        """
        return []
        # extractor = FeatureExtractor(self)
        # return extractor.extract_around_features(features, upstream, downstream, **kwargs)
    
    def get_features_at_position(self, chrom: Union[str, int], position: int,
                                exclude_types: Optional[List[str]] = None,
                                include_types: Optional[List[str]] = None) -> List[Any]:
        """Get features at a specific genomic position.
        
        Now uses unified annotation system for richer annotations.
        
        Args:
            chrom: Chromosome
            position: Genomic position (1-based)  
            exclude_types: List of feature types to exclude
            include_types: List of feature types to include
            
        Returns:
            List of features at the position
        """
        # Use unified annotations for more comprehensive results
        unified_features = self.annotations.query_point(str(chrom), position)
        
        # Convert to old Feature format for backward compatibility
        features = []
        for uf in unified_features:
            # Apply filters
            if exclude_types and uf.feature_type in exclude_types:
                continue
            if include_types and uf.feature_type not in include_types:
                continue
                
            # Convert UFeature to dict-like format expected by old code
            feature_dict = {
                'feature': uf.feature_type,
                'start': uf.start,
                'end': uf.end,
                'chrom': uf.chrom,
                'strand': uf.strand,
                'info': uf.attributes.copy()
            }
            if uf.name:
                feature_dict['info']['name'] = uf.name
                feature_dict['info']['gene_name'] = uf.name
            features.append(feature_dict)
        
        # Also get from old gene_map for complete backward compatibility
        old_features = self.gene_map.get_feature(chrom, position)
        
        if exclude_types:
            old_features = [f for f in old_features 
                       if f.get('feature') not in exclude_types]
        if include_types:
            old_features = [f for f in old_features 
                        if f.get('feature') in include_types]
        
        # Merge, avoiding duplicates
        seen = set()
        for f in old_features:
            key = (f.get('start'), f.get('end'), f.get('feature'))
            if key not in seen:
                features.append(f)
                seen.add(key)
        
        return features
    
    def browse(self, chrom: Union[str, int] = 1, position: int = 1000000,
              window_size: int = 80) -> None:
        """Launch interactive genome browser.
        
        Args:
            chrom: Starting chromosome
            position: Starting position
            window_size: Display window size in base pairs
        """
        browser = InteractiveGenomeBrowser(self)
        browser.start(chrom, position, window_size)
    
    def find_next_feature(self, chrom: Union[str, int], position: int,
                          feature_type: str = 'gene',
                          max_distance: int = 10000000) -> Optional[Dict[str, Any]]:
        """Find the next feature of a given type downstream from position.
        
        Args:
            chrom: Chromosome to search
            position: Current position (1-based)
            feature_type: Type of feature to find (e.g., 'gene', 'exon', 'CDS')
            max_distance: Maximum distance to search
            
        Returns:
            Dictionary with feature information or None if not found
        """
        try:
            # Use gene_map to fetch features downstream
            end_pos = position + max_distance
            features = list(self.gene_map.fetch(
                chrom, position + 1, end_pos, 
                features=(feature_type,) if feature_type else tuple()
            ))
            
            if features:
                # Sort by start position and get the closest
                features.sort(key=lambda f: f.get('start', float('inf')))
                next_feature = features[1]
                
                return {
                    'feature': next_feature,
                    'type': next_feature.get('feature', feature_type),
                    'start': next_feature.get('start'),
                    'end': next_feature.get('end'),
                    'name': next_feature.get('info', {}).get('gene_name', 
                            next_feature.get('info', {}).get('transcript_name', 'unknown')),
                    'distance': next_feature.get('start', 0) - position
                }
            
        except Exception as e:
            logger.warning(f"Error finding next feature: {e}")
        
        return None
    
    def find_prev_feature(self, chrom: Union[str, int], position: int,
                         feature_type: str = 'gene',
                         max_distance: int = 10000000) -> Optional[Dict[str, Any]]:
        """Find the previous feature of a given type upstream from position.
        
        Args:
            chrom: Chromosome to search
            position: Current position (1-based)
            feature_type: Type of feature to find
            max_distance: Maximum distance to search
            
        Returns:
            Dictionary with feature information or None if not found
        """
        try:
            # Use gene_map to fetch features upstream
            start_pos = max(1, position - max_distance)
            features = list(self.gene_map.fetch(
                chrom, start_pos, position - 1,
                features=(feature_type,) if feature_type else tuple()
            ))
            
            if features:
                # Sort by end position and get the closest
                features.sort(key=lambda f: f.get('end', 0), reverse=True)
                prev_feature = features[0]
                
                return {
                    'feature': prev_feature,
                    'type': prev_feature.get('feature', feature_type),
                    'start': prev_feature.get('start'),
                    'end': prev_feature.get('end'),
                    'name': prev_feature.get('info', {}).get('gene_name',
                            prev_feature.get('info', {}).get('transcript_name', 'unknown')),
                    'distance': position - prev_feature.get('end', 0)
                }
            
        except Exception as e:
            logger.warning(f"Error finding previous feature: {e}")
        
        return None
    
    def scan_for_motif(self, chrom: Union[str, int], start: int, end: int,
                      motif: str, integrate_variants: bool = True) -> List[int]:
        """Scan a genomic region for a specific motif.
        
        Args:
            chrom: Chromosome to scan
            start: Start position
            end: End position
            motif: DNA motif to search for
            integrate_variants: Whether to search in personalized sequence
            
        Returns:
            List of positions where motif was found
        """
        positions = []
        
        iterator = UGenomeIterator(
            self, chrom, start, end,
            window_size=len(motif),
            integrate_variants=integrate_variants
        )
        
        for ref_seq, personal_seq, _ in iterator:
            seq_to_check = personal_seq if integrate_variants else ref_seq
            if seq_to_check.upper() == motif.upper():
                positions.append(iterator.current_pos - iterator.stride)
        
        return positions

    ###### iterators #######
    
    @classmethod
    def iter_chromes(cls, include_mito = False):
        
        mt_chr = []
        if include_mito:
            mt_chr = ['MT']
        
        for chr in [str(i) for i in range(1, 23)] + ['X','Y'] + mt_chr:
            yield chr

