#!/usr/bin/env python3
"""Quick test of genome browser variant display."""

import sys
sys.path.insert(0, '/home/gront/Documents/python/genomics-prj')

from ggene.genomemanager import GenomeManager
from ggene.genome_browser import InteractiveGenomeBrowser

def test_browser_variants():
    """Quick test to see if variants display in browser."""
    
    print("Testing Genome Browser Variant Display")
    print("=" * 60)
    
    # Initialize
    gm = GenomeManager()
    
    # Create browser
    browser = InteractiveGenomeBrowser(gm)
    
    # Set up state for testing (without starting interactive mode)
    from ggene.genome_browser import BrowserState
    browser.state = BrowserState(
        chrom=1,
        position=204195200,
        window_size=80,
        stride=20,
        show_features=True
    )
    
    # Just display the current view to test rendering
    print("\nRendering browser view at position with large insertion:")
    print("-" * 60)
    
    try:
        # This will show what the browser would display
        browser._display_current_view()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("If you see aligned sequences with gaps and variant info above,")
    print("the variant display is working correctly!")

if __name__ == "__main__":
    test_browser_variants()