#!/usr/bin/env python3
"""Detailed test for motif detection - scan a larger region."""

import logging
from ggene.genomemanager import GenomeManager
from ggene.genome_iterator_v2 import UnifiedGenomeIterator

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

def scan_for_motifs():
    """Scan a region and report all motifs found."""
    
    print("Initializing genome manager...")
    gm = GenomeManager()
    
    # Create iterator with larger window to find motifs
    print("\nScanning chromosome 1 from position 200000-201000 for motifs...")
    iterator = UnifiedGenomeIterator(
        gm,
        chrom='1', 
        start=200000,
        end=201000,
        window_size=500,
        stride=500,
        detect_motifs=True
    )
    
    motifs_found = []
    positions_scanned = []
    
    # Iterate through windows
    for window in iterator:
        positions_scanned.append(f"{window.start_ref}-{window.end_ref}")
        
        if window.motifs:
            print(f"\n✓ Found {len(window.motifs)} motifs in window {window.chrom}:{window.start_ref}-{window.end_ref}")
            for motif in window.motifs:
                motifs_found.append(motif)
                print(f"  - {motif['name']:15} at position {motif['start']:10} "
                      f"(score: {motif.get('score', 1.0):.2f}) "
                      f"seq: {motif.get('sequence', 'N/A')}")
    
    print(f"\n" + "="*60)
    print(f"SUMMARY:")
    print(f"Scanned {len(positions_scanned)} windows")
    print(f"Total motifs found: {len(motifs_found)}")
    
    if motifs_found:
        # Count by type
        motif_counts = {}
        for m in motifs_found:
            name = m['name']
            motif_counts[name] = motif_counts.get(name, 0) + 1
        
        print(f"\nMotif counts by type:")
        for name, count in sorted(motif_counts.items()):
            print(f"  {name:20}: {count:3} occurrences")
        
        # Show first few examples
        print(f"\nFirst 5 motifs:")
        for m in motifs_found[:5]:
            print(f"  {m['name']} at {m['start']}: {m.get('sequence', '')}")
    
    return len(motifs_found) > 0

if __name__ == "__main__":
    print("="*60)
    print("MOTIF DETECTION DETAILED TEST")
    print("="*60)
    
    if scan_for_motifs():
        print("\n✓ SUCCESS: Motifs detected!")
    else:
        print("\n⚠ WARNING: No motifs found in scanned region")
        print("This could mean:")
        print("  1. The region doesn't contain the patterns we're looking for")
        print("  2. The patterns need adjustment for this genome")
        print("  3. Try scanning a different region")