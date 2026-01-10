
import argparse
import numpy as np

from ggene.database.genome_manager import GenomeManager
from ggene.config import get_paths, DATA_DIR
from ggene.browser.builds.test_browser import TestBrowser, TestBrowserState
from ggene.database.annotations import get_experiment_stream
from ggene.database.databases import find_and_load_database, find_databases, find_directory
from ggene.draw import Colors

def load_genome(**kwargs):
    gm = GenomeManager(**kwargs)
    return gm

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", "-c", default = "1", type = str,
                        help="Chromosome (1-23, X, Y, MT)")
    parser.add_argument("--position", "-p", default = int(1e6), type = int,
                        help="Window size, in bp (default 4096)")
    parser.add_argument("--window", "-w", default = 4096, type = int,
                        help="Window size, in bp (default 4096)", dest = "window_size")
    parser.add_argument("--stride", "-s", default = 256, type = int, 
                        help="Window size, in bp (default 4096)")
    parser.add_argument("--zoom", "-z", default = 0.0, type = float, 
                        help="Every additional 1 increases zoom by 4x")
    parser.add_argument("--random", "-r", default = False, action = "store_true", 
                        help="Choose random chromosome and position")
    
    parser.add_argument("--jaspars", "-j", default = False, action = "store_true", 
                        help="load jaspars?")
    parser.add_argument("--patterns", "-pt", default = False, action = "store_true", 
                        help="load patterns?")
    
    parser.add_argument("--load", "-l", default = False, action = "store_true", 
                        help="Load last state")
    parser.add_argument("--debug", "-d", default = False, action = "store_true", 
                        help="Debug mode")
                   
    display_width = 256
    args = parser.parse_args()
    
    argdict = args.__dict__
    # argdict["window_size"]
    
    a,b,c, lib_path, other_paths = get_paths()
    
    org_dir = DATA_DIR / "S_rosetta"
    gtf_path = org_dir / "GCF_000188695.1_Proterospongia_sp_ATCC50818_genomic.gtf.gz"
    ref_path = org_dir / "GCF_000188695.1_Proterospongia_sp_ATCC50818_genomic.fna.gz"
    
    gm = load_genome(gtf_path = gtf_path, ref_path = ref_path, vcf_path = "", skip_repeatmasker = True, skip_dfam = True, skip_clinvar = True, load_jaspars = False, load_patterns = False)
    
    brws = TestBrowser(gm, **argdict)
    brws.start(**argdict)     
                        
                        
                        
                    
if __name__=="__main__":
    main()

