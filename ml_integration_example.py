#!/usr/bin/env python3
"""
Example integration with ML libraries (kipoi/selene) for genomic analysis.
This demonstrates how to use the genome iterator to prepare data for ML models.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional
from ggene.genomemanager import GenomeManager
from ggene.genome_iterator import FeatureExtractor

logger = logging.getLogger(__name__)


class GenomicMLPipeline:
    """Pipeline for preparing genomic data for ML models."""
    
    def __init__(self, genome_manager: GenomeManager):
        """Initialize the ML pipeline.
        
        Args:
            genome_manager: Initialized GenomeManager instance
        """
        self.gm = genome_manager
        self.extractor = FeatureExtractor(genome_manager)
    
    def prepare_regulatory_regions(self, 
                                  gene_names: List[str],
                                  promoter_upstream: int = 2000,
                                  promoter_downstream: int = 500) -> Dict[str, Any]:
        """Extract promoter regions for transcription factor binding analysis.
        
        Args:
            gene_names: List of gene names to analyze
            promoter_upstream: Bases upstream of TSS
            promoter_downstream: Bases downstream of TSS
            
        Returns:
            Dictionary with prepared sequences and metadata
        """
        promoter_sequences = []
        
        for gene_name in gene_names:
            try:
                # Find gene
                chrom, gene_info = self.gm.gene_map.find_gene(gene_name)
                if not chrom:
                    logger.warning(f"Gene {gene_name} not found")
                    continue
                
                # Get transcription start site (TSS)
                gene_data = self.gm.gene_map.get_gene(gene_name, chrom)
                if not gene_data or 'transcript' not in gene_data:
                    continue
                
                # Use first transcript's start as TSS
                transcript = gene_data['transcript'][0]
                tss = transcript['start'] if transcript['strand'] == '+' else transcript['end']
                
                # Create feature for promoter region
                from ggene.features import Feature
                promoter = Feature({
                    'feature': 'promoter',
                    'chrom': chrom,
                    'start': tss - promoter_upstream if transcript['strand'] == '+' else tss - promoter_downstream,
                    'end': tss + promoter_downstream if transcript['strand'] == '+' else tss + promoter_upstream,
                    'gene_name': gene_name,
                    'strand': transcript['strand']
                })
                
                # Extract sequence
                seq_data = self.extractor.extract_around_features(
                    [promoter], 
                    upstream=0,  # Already included in promoter definition
                    downstream=0,
                    integrate_variants=True
                )
                
                if seq_data:
                    promoter_sequences.extend(seq_data)
                    
            except Exception as e:
                logger.error(f"Error processing gene {gene_name}: {e}")
        
        # Prepare for ML
        if promoter_sequences:
            return self.extractor.prepare_for_ml(
                promoter_sequences,
                sequence_length=promoter_upstream + promoter_downstream,
                one_hot_encode=True
            )
        
        return {'sequences': np.array([]), 'metadata': []}
    
    def prepare_for_kipoi(self, sequences: np.ndarray, 
                         model_name: str = "Basset") -> np.ndarray:
        """Prepare sequences for Kipoi models.
        
        Args:
            sequences: One-hot encoded sequences (N, L, 4 or 5)
            model_name: Kipoi model name
            
        Returns:
            Formatted array for Kipoi model
        """
        # Remove the 'N' channel if present (5th channel)
        if sequences.shape[-1] == 5:
            sequences = sequences[:, :, :4]
        
        # Kipoi models often expect (N, 4, L) format
        sequences = np.transpose(sequences, (0, 2, 1))
        
        logger.info(f"Prepared {sequences.shape[0]} sequences for {model_name}")
        logger.info(f"Shape: {sequences.shape}")
        
        return sequences
    
    def prepare_for_selene(self, sequences: np.ndarray) -> np.ndarray:
        """Prepare sequences for Selene models.
        
        Args:
            sequences: One-hot encoded sequences
            
        Returns:
            Formatted array for Selene
        """
        # Selene typically expects (N, L, 4) format
        if sequences.shape[-1] == 5:
            sequences = sequences[:, :, :4]
        
        logger.info(f"Prepared {sequences.shape[0]} sequences for Selene")
        logger.info(f"Shape: {sequences.shape}")
        
        return sequences
    
    def extract_variant_context_windows(self,
                                       chrom: Union[str, int],
                                       variant_positions: List[int],
                                       window_size: int = 1000) -> Dict[str, Any]:
        """Extract sequences around variants for effect prediction.
        
        Args:
            chrom: Chromosome
            variant_positions: List of variant positions
            window_size: Total window size (variant will be centered)
            
        Returns:
            Dictionary with reference and alternate sequences
        """
        half_window = window_size // 2
        results = {
            'reference_sequences': [],
            'alternate_sequences': [],
            'variant_info': []
        }
        
        for var_pos in variant_positions:
            # Get reference sequence
            iterator_ref = self.gm.iterate_genome(
                chrom,
                var_pos - half_window,
                var_pos + half_window,
                integrate_variants=False,
                window_size=window_size
            )
            
            # Get alternate (personal) sequence
            iterator_alt = self.gm.iterate_genome(
                chrom,
                var_pos - half_window,
                var_pos + half_window,
                integrate_variants=True,
                window_size=window_size
            )
            
            try:
                ref_seq, _, _ = next(iterator_ref)
                alt_seq, _, _ = next(iterator_alt)
                
                results['reference_sequences'].append(ref_seq)
                results['alternate_sequences'].append(alt_seq)
                results['variant_info'].append({
                    'chrom': chrom,
                    'position': var_pos,
                    'window_start': var_pos - half_window,
                    'window_end': var_pos + half_window
                })
                
            except StopIteration:
                logger.warning(f"Could not extract window for variant at {var_pos}")
        
        # Convert to one-hot encoding
        ref_encoded = self._encode_sequences(results['reference_sequences'])
        alt_encoded = self._encode_sequences(results['alternate_sequences'])
        
        return {
            'reference': ref_encoded,
            'alternate': alt_encoded,
            'metadata': results['variant_info']
        }
    
    def _encode_sequences(self, sequences: List[str]) -> np.ndarray:
        """One-hot encode a list of sequences.
        
        Args:
            sequences: List of DNA sequences
            
        Returns:
            One-hot encoded array
        """
        if not sequences:
            return np.array([])
        
        max_len = max(len(s) for s in sequences)
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        
        encoded = np.zeros((len(sequences), max_len, 5))
        
        for i, seq in enumerate(sequences):
            for j, base in enumerate(seq.upper()):
                if base in base_to_idx:
                    encoded[i, j, base_to_idx[base]] = 1
                else:
                    encoded[i, j, 4] = 1  # Unknown
        
        return encoded
    
    def find_tfbs_candidates(self, 
                            chrom: Union[str, int],
                            start: int,
                            end: int,
                            tf_motifs: Dict[str, str]) -> Dict[str, List[int]]:
        """Find potential transcription factor binding sites.
        
        Args:
            chrom: Chromosome to search
            start: Start position
            end: End position
            tf_motifs: Dictionary of TF name to motif sequence
            
        Returns:
            Dictionary mapping TF names to lists of binding positions
        """
        tfbs_positions = {}
        
        for tf_name, motif in tf_motifs.items():
            positions = self.gm.scan_for_motif(
                chrom, start, end, motif, 
                integrate_variants=True
            )
            
            if positions:
                tfbs_positions[tf_name] = positions
                logger.info(f"Found {len(positions)} potential {tf_name} binding sites")
        
        return tfbs_positions


def example_tf_binding_analysis():
    """Example: Analyze transcription factor binding in promoters."""
    print("\n=== Transcription Factor Binding Analysis ===")
    
    # Initialize
    gm = GenomeManager()
    pipeline = GenomicMLPipeline(gm)
    
    # Define some common TF binding motifs
    tf_motifs = {
        'TATA': 'TATAAA',
        'CAAT': 'CCAAT',
        'GC_box': 'GGGCGG',
        'E_box': 'CACGTG',
        'AP1': 'TGACTCA'
    }
    
    # Analyze a specific gene's promoter
    gene_name = 'KISS1'  # Change to your gene of interest
    
    try:
        # Get gene information
        chrom, gene_info = gm.gene_map.find_gene(gene_name)
        if chrom:
            gene_data = gm.gene_map.get_gene(gene_name, chrom)
            if gene_data and 'transcript' in gene_data:
                transcript = gene_data['transcript'][0]
                tss = transcript['start']
                
                # Search for TF binding sites in promoter region
                tfbs = pipeline.find_tfbs_candidates(
                    chrom,
                    tss - 2000,  # 2kb upstream
                    tss + 500,   # 500bp downstream
                    tf_motifs
                )
                
                print(f"\nTF binding sites near {gene_name} TSS:")
                for tf_name, positions in tfbs.items():
                    if positions:
                        distances = [pos - tss for pos in positions]
                        print(f"  {tf_name}: {len(positions)} sites")
                        print(f"    Distances from TSS: {distances[:5]}...")
    
    except Exception as e:
        logger.error(f"Error in TF analysis: {e}")


def example_variant_effect_prediction():
    """Example: Prepare data for variant effect prediction."""
    print("\n=== Variant Effect Prediction Data Preparation ===")
    
    gm = GenomeManager()
    pipeline = GenomicMLPipeline(gm)
    
    # Find some variants to analyze
    chrom = 1
    region_start = 1000000
    region_end = 1100000
    
    # Get variant positions
    variant_positions = []
    for var in gm.vcf(gm._make_index(chrom, region_start, region_end)):
        if var.QUAL > 20:  # Quality filter
            variant_positions.append(var.POS)
        if len(variant_positions) >= 10:  # Limit for example
            break
    
    print(f"Found {len(variant_positions)} high-quality variants")
    
    if variant_positions:
        # Extract context windows
        variant_data = pipeline.extract_variant_context_windows(
            chrom,
            variant_positions[:5],  # First 5 variants
            window_size=1000
        )
        
        print(f"\nPrepared variant effect prediction data:")
        print(f"  Reference sequences shape: {variant_data['reference'].shape}")
        print(f"  Alternate sequences shape: {variant_data['alternate'].shape}")
        print(f"  Number of variants: {len(variant_data['metadata'])}")
        
        # This data can now be fed to variant effect prediction models
        # For example, with DeepSEA, ExPecto, or similar models


def example_kipoi_integration():
    """Example: Prepare data for Kipoi models."""
    print("\n=== Kipoi Model Data Preparation ===")
    
    gm = GenomeManager()
    pipeline = GenomicMLPipeline(gm)
    
    # Prepare promoter regions for multiple genes
    genes = ['KISS1', 'HTR2A', 'HTR2B']  # Your genes of interest
    
    promoter_data = pipeline.prepare_regulatory_regions(
        genes,
        promoter_upstream=1000,
        promoter_downstream=500
    )
    
    if promoter_data['sequences'].size > 0:
        # Format for Kipoi
        kipoi_data = pipeline.prepare_for_kipoi(
            promoter_data['sequences'],
            model_name="Basset"
        )
        
        print(f"\nPrepared data for Kipoi:")
        print(f"  Shape: {kipoi_data.shape}")
        print(f"  Format: (batch, channels, length)")
        
        # Now you can use this with Kipoi:
        # import kipoi
        # model = kipoi.get_model("Basset")
        # predictions = model.predict_on_batch(kipoi_data)


def main():
    """Run ML integration examples."""
    print("=" * 60)
    print("GENOMIC ML PIPELINE EXAMPLES")
    print("=" * 60)
    
    examples = [
        example_tf_binding_analysis,
        example_variant_effect_prediction,
        example_kipoi_integration
    ]
    
    for example in examples:
        try:
            example()
        except Exception as e:
            logger.error(f"Error in example: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("ML integration examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    import sys
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    main()