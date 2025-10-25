#!/usr/bin/env python3
"""Test script for motif detection in genome browser."""

import logging
import sys
from ggene.genomemanager import GenomeManager
from ggene.genome_browser import InteractiveGenomeBrowser
from ggene.genome_iterator_v2 import UnifiedGenomeIterator

# Setup logging to see debug messages
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def test_motif_detection():
    """Test that motifs are being detected and displayed."""
    
    # Initialize genome manager
    print("Initializing genome manager...")
    gm = GenomeManager()
    
    # Create iterator with motif detection
    print("\nCreating genome iterator with motif detection...")
    iterator = UnifiedGenomeIterator(
        gm,
        chrom='1',
        start=200000,
        window_size=200,
        detect_motifs=True
    )
    
    # Get a window and check for motifs
    print("\nGetting window and checking for motifs...")
    window = iterator.get_window_at(200000)
    
    print(f"Window: chr{window.chrom}:{window.start_ref}-{window.end_ref}")
    print(f"Reference sequence (first 50bp): {window.ref_seq[:50]}...")
    print(f"Alternate sequence (first 50bp): {window.alt_seq[:50]}...")
    
    if window.motifs:
        print(f"\nFound {len(window.motifs)} motifs:")
        for motif in window.motifs:
            print(f"  - {motif['name']} at {motif['start']}-{motif['end']} "
                  f"(score: {motif.get('score', 'N/A')})")
            print(f"    Sequence: {motif.get('sequence', 'N/A')}")
    else:
        print("\nNo motifs detected in this window")
    
    # Also test with a known sequence that should have motifs
    print("\n" + "="*50)
    print("Testing with known motif-containing region...")
    
    # Look for regions with genes (more likely to have regulatory motifs)
    if hasattr(gm, 'annotations'):
        genes = gm.annotations.query_range('1', 1000000, 2000000)
        if genes:
            gene = genes[0]
            print(f"\nChecking near gene at position {gene.start}...")
            
            # Create iterator near gene start (promoter region)
            iterator2 = UnifiedGenomeIterator(
                gm,
                chrom='1',
                start=max(1, gene.start - 500),
                window_size=500,
                detect_motifs=True
            )
            
            window2 = iterator2.get_window_at(max(1, gene.start - 500))
            
            if window2.motifs:
                print(f"Found {len(window2.motifs)} motifs near gene:")
                for motif in window2.motifs[:5]:  # Show first 5
                    print(f"  - {motif['name']} at {motif['start']}-{motif['end']}")
            else:
                print("No motifs found near gene")
    
    return True

def test_browser_display():
    """Test that motifs are displayed in the browser."""
    print("\n" + "="*50)
    print("Testing browser display...")
    print("Starting interactive browser. Look for motifs in the feature track.")
    print("Press 'q' to quit the browser.")
    print("="*50 + "\n")
    
    gm = GenomeManager()
    browser = InteractiveGenomeBrowser(gm)
    
    # Start browser at a position likely to have motifs
    browser.start(chrom='1', position=200000, window_size=160)

if __name__ == "__main__":
    print("Testing motif detection in genome browser")
    print("="*50)
    
    # Run detection test
    if test_motif_detection():
        print("\n✓ Motif detection test completed")
    
    # Ask if user wants to test browser
    response = input("\nDo you want to test the interactive browser? (y/n): ")
    if response.lower() == 'y':
        test_browser_display()
    
    print("\nAll tests completed!")