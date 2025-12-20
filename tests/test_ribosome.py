#!/usr/bin/env python3
"""
Test the new Ribosome class for transcript translation.
"""

import logging
from ggene.database.genomemanager import GenomeManager

# Configure logging to see what's happening
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def test_gene_translation(gm: GenomeManager, gene_name: str):
    """Test translation of a gene."""
    print(f"\n{'='*60}")
    print(f"Testing translation of {gene_name}")
    print('='*60)
    
    # First, load the gene to see what transcripts are available
    gene = gm.load_gene(gene_name)
    if not gene:
        print(f"Could not load gene {gene_name}")
        return
    
    # Check transcripts
    if hasattr(gene, 'transcripts'):
        num_transcripts = len(gene.transcripts)
        print(f"Found {num_transcripts} transcript(s) for {gene_name}")
    else:
        print(f"No transcripts found for {gene_name}")
        return
    
    # Translate the first transcript
    print(f"\nTranslating transcript 0...")
    result = gm.transcribe_gene(gene_name, transcript_index=0, personal=True)
    
    if result:
        print(f"\n--- Translation Results ---")
        print(f"Protein length: {result['protein_length']} amino acids")
        print(f"Start position: {result['start_position']}")
        print(f"Stop position: {result['stop_position']}")
        print(f"Strand: {result['strand']}")
        print(f"Frame: {result['frame']}")
        print(f"CDS segments: {len(result['cds_segments'])}")
        
        if result['warnings']:
            print(f"\nWarnings:")
            for warning in result['warnings']:
                print(f"  - {warning}")
        
        # Show first 100 amino acids
        protein_seq = result['protein_sequence']
        if protein_seq:
            print(f"\nProtein sequence (first 100 AA):")
            print(f"{protein_seq[:100]}...")
            
            # Show last 50 amino acids
            if len(protein_seq) > 100:
                print(f"\nProtein sequence (last 50 AA):")
                print(f"...{protein_seq[-50:]}")
        
        # Check for mutations
        if result['mutation_summary']:
            print(f"\nMutations detected: {result['mutation_summary']}")
    else:
        print("Translation failed")


def compare_personal_vs_reference(gm: GenomeManager, gene_name: str):
    """Compare personal and reference protein sequences."""
    print(f"\n{'='*60}")
    print(f"Comparing personal vs reference for {gene_name}")
    print('='*60)
    
    comparison = gm.compare_transcribed(gene_name, transcript_index=0)
    
    if comparison:
        print(f"\n--- Comparison Results ---")
        print(f"Reference length: {comparison['reference_length']} AA")
        print(f"Personal length: {comparison['personal_length']} AA")
        print(f"Length difference: {comparison['length_difference']}")
        print(f"Number of changes: {comparison['num_changes']}")
        print(f"Percent identity: {comparison['percent_identity']:.1f}%")
        
        if comparison['amino_acid_changes']:
            print(f"\nAmino acid changes:")
            for change in comparison['amino_acid_changes'][:10]:  # Show first 10
                print(f"  {change['notation']}")
            
            if len(comparison['amino_acid_changes']) > 10:
                print(f"  ... and {len(comparison['amino_acid_changes']) - 10} more")
        
        if comparison['mutation_summary']:
            print(f"\nMutation summary: {comparison['mutation_summary']}")
        
        # Show alignment of a region with changes
        if comparison['amino_acid_changes']:
            first_change = comparison['amino_acid_changes'][0]
            pos = first_change['position'] - 1  # 0-based
            
            # Show 20 AA around the first change
            start = max(0, pos - 10)
            end = min(pos + 10, min(len(comparison['reference_sequence']), 
                                   len(comparison['personal_sequence'])))
            
            print(f"\nAlignment around position {first_change['position']}:")
            print(f"Ref: {comparison['reference_sequence'][start:end]}")
            print(f"Per: {comparison['personal_sequence'][start:end]}")
            print(f"     {' ' * (pos - start)}^")
    else:
        print("Comparison failed")


def test_multiple_genes(gm: GenomeManager):
    """Test translation of multiple genes."""
    test_genes = [
        'KISS1',
        'HTR2A',
        'HTR2B',
        'SLC6A4',
        'FTO'
    ]
    
    print(f"\n{'='*60}")
    print(f"Testing multiple genes")
    print('='*60)
    
    results = []
    
    for gene_name in test_genes:
        print(f"\nProcessing {gene_name}...")
        
        # Try to translate
        result = gm.transcribe_gene(gene_name, transcript_index=0, personal=True)
        
        if result:
            results.append({
                'gene': gene_name,
                'length': result['protein_length'],
                'mutations': len(result.get('mutations', [])),
                'warnings': len(result.get('warnings', []))
            })
            print(f"  ✓ Success: {result['protein_length']} AA, "
                  f"{len(result.get('mutations', []))} mutations")
        else:
            print(f"  ✗ Failed")
    
    # Summary
    print(f"\n--- Summary ---")
    print(f"Successfully translated {len(results)}/{len(test_genes)} genes")
    
    if results:
        total_mutations = sum(r['mutations'] for r in results)
        print(f"Total mutations found: {total_mutations}")
        
        # Show genes with most mutations
        results.sort(key=lambda x: x['mutations'], reverse=True)
        print(f"\nGenes with most mutations:")
        for r in results[:3]:
            if r['mutations'] > 0:
                print(f"  {r['gene']}: {r['mutations']} mutations")


def main():
    """Run all tests."""
    print("Loading genome data...")
    gm = GenomeManager()
    
    # Test menu
    print("\n" + "="*60)
    print("RIBOSOME TRANSLATION TESTS")
    print("="*60)
    print("\n1. Test single gene translation (KISS1)")
    print("2. Compare personal vs reference (KISS1)")
    print("3. Test multiple genes")
    print("4. Custom gene")
    print("0. Exit")
    
    try:
        choice = input("\nSelect test (0-4): ")
        choice = int(choice)
        
        if choice == 1:
            test_gene_translation(gm, 'KISS1')
        elif choice == 2:
            compare_personal_vs_reference(gm, 'KISS1')
        elif choice == 3:
            test_multiple_genes(gm)
        elif choice == 4:
            gene_name = input("Enter gene name: ").strip().upper()
            test_gene_translation(gm, gene_name)
            compare_personal_vs_reference(gm, gene_name)
        elif choice == 0:
            print("Exiting...")
        else:
            print("Invalid choice")
            
    except KeyboardInterrupt:
        print("\nInterrupted")
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()