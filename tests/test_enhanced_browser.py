#!/usr/bin/env python3
"""
Test script for the enhanced genome browser.
Shows off new features like hotkeys, state saving, and improved compression.
"""

from ggene.genomemanager import GenomeManager

def main():
    print("Enhanced Genome Browser Test")
    print("=" * 40)
    
    # Initialize GenomeManager
    print("Initializing GenomeManager...")
    gm = GenomeManager()
    
    print("\nNew Features Added:")
    print("✓ Ctrl+→/← for fast gene jumping")
    print("✓ Ctrl+↑/↓ for save/load state")  
    print("✓ Fixed 24-line terminal height")
    print("✓ Feature type labels on left side")
    print("✓ Smart feature compression (ex1/2/3, GENE1/+2)")
    print("✓ Enhanced help system")
    
    print("\nStarting browser...")
    print("Try navigating to LINC01128 with 'g' then 'LINC01128'")
    print("Use Ctrl+→ to jump between genes quickly!")
    print("Press '?' for full help menu")
    
    # Start the enhanced browser
    from ggene.genome_browser_v2 import InteractiveGenomeBrowser
    browser = InteractiveGenomeBrowser(gm)
    browser.start(chrom=1, position=1000000, window_size=120)


if __name__ == "__main__":
    main()