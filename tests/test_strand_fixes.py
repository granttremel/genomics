#!/usr/bin/env python3
"""
Test the strand-aware fixes for genome browser.
"""

from ggene.database.genomemanager import GenomeManager

def main():
    print("Testing Strand-Aware Fixes")
    print("=" * 50)
    
    # Initialize GenomeManager
    print("Initializing GenomeManager...")
    gm = GenomeManager()
    
    print("\nFixes Applied:")
    print("✓ Start/stop codon highlighting now reflects for - strand")
    print("✓ Amino acid positions fixed (no more 2.5x window size)")
    print("✓ Amino acid display properly reflects for - strand")
    print("✓ CDS boundary markers are strand-aware")
    
    print("\nTest Instructions:")
    print("1. Start the browser")
    print("2. Navigate to a CDS region with variants")
    print("3. Press 'a' to toggle amino acid display")
    print("4. Press '-' to toggle between + and - strand")
    print("5. Observe that:")
    print("   - Start/stop codons move correctly")
    print("   - Amino acids align with sequence")
    print("   - CDS boundaries are in correct positions")
    
    print("\nStarting browser at a CDS region...")
    
    # Start browser in a region likely to have CDS features
    from ggene.genome_browser_v2 import InteractiveGenomeBrowser
    browser = InteractiveGenomeBrowser(gm)
    browser.start(chrom=1, position=1000000, window_size=80)

if __name__ == "__main__":
    main()