#!/usr/bin/env python3
"""Test the genome browser with motif detection display."""

from ggene.database.genome_manager import GenomeManager
from ggene.genome_browser_v2 import InteractiveGenomeBrowser

def main():
    print("Starting genome browser with motif detection...")
    print("="*60)
    print("Look for:")
    print("  1. Colored underlines showing motif locations")
    print("  2. Strand indicators (> for + strand, < for - strand)")
    print("  3. Motifs section listing detected motifs")
    print("\nMotif color scheme:")
    print("  Red underline    - Splice sites")
    print("  Magenta underline - TATA box")
    print("  Green underline   - CpG islands")
    print("  Blue underline    - PolyA signals")
    print("  Cyan underline    - Other motifs")
    print("="*60)
    print("\nPress 'q' to quit the browser\n")
    
    # Initialize and start browser
    gm = GenomeManager()
    browser = InteractiveGenomeBrowser(gm)
    
    # Start at a position known to have motifs
    browser.start(chrom='1', position=200400, window_size=160)

if __name__ == "__main__":
    main()