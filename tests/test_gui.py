#!/usr/bin/env python3
"""
Test script for launching the PyQt genome browser.
"""
import sys
import logging
from ggene.genomemanager import GenomeManager
from ggene.genome_browser_gui import launch_genome_browser

# Set up logging
logging.basicConfig(level=logging.INFO)

def main():
    """Main function to test the genome browser GUI."""
    print("Initializing GenomeManager...")
    
    try:
        # Initialize the GenomeManager with your data files
        gm = GenomeManager()
        
        print("Launching Genome Browser GUI...")
        # Launch the GUI
        launch_genome_browser(gm)
        
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure you have the required data files (VCF, GTF, FASTA) configured")
        sys.exit(1)

if __name__ == "__main__":
    main()