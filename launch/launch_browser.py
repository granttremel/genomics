
import argparse
import numpy as np

from ggene.database.genome_manager import GenomeManager
from ggene.config import get_paths, DATA_DIR
from ggene.browser.builds.test_browser import TestBrowser, TestBrowserState
from ggene.browser.builds.justlines_browser import JustLinesBrowser, JustLinesBrowserState
from ggene.browser.builds.scalar_browser import ScalarBrowser, ScalarBrowserState
from ggene.browser.builds.seqs_browser import SeqsBrowser, SeqsBrowserState
from ggene.browser.builds.map_browser import MapBrowser, MapBrowserState
from ggene.browser.builds.comp_browser import CompBrowser, CompBrowserState
from ggene.database.annotations import get_experiment_stream
from ggene.database.databases import find_and_load_database, find_databases, find_directory
from ggene.draw import Colors

def get_browser(browser_type:str, gm, **kwargs):
    
    if browser_type == "test":
        return TestBrowser(gm=gm, **kwargs)
    elif browser_type == "scalar":
        return ScalarBrowser(gm=gm, **kwargs)
    elif browser_type == "seqs":
        return SeqsBrowser(gm=gm, **kwargs)
    elif browser_type == "lines":
        return JustLinesBrowser(gm=gm, **kwargs)
    elif browser_type == "comp":
        return CompBrowser(gm=gm, **kwargs)
    elif browser_type == "map":
        return MapBrowser(gm, **kwargs)
    else:
        return None

def load_genome(fast = True, **kwargs):
    a,b,c, lib_path, other_paths = get_paths()
    kwargs.update(other_paths)
    if fast:
        kwargs.update(skip_geneinfo = True, load_jaspars = False, load_patterns = False)
    
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
    parser.add_argument("--browser-type", "-b", default = "test", type = str, 
                        help="Browser type (test, scalar, lines, seqs, comp)")
    
    parser.add_argument("--fast", "-f", action = "store_true", help = "Load GenomeManager fast style? (no gene info)")
    
    display_width = 256
    
    args = parser.parse_args()
    
    argdict = args.__dict__
    gm = load_genome(fast = argdict.pop("fast", True), load_jaspars = argdict.pop("jaspars", False), load_patterns = argdict.pop("patterns", False))
    
    argdict["display_width"] = display_width
    btp = argdict.pop("browser_type","")
        
    if args.random:
        chrom, pos = GenomeManager.get_random_location()
        argdict["chrom"] = chrom
        argdict["position"] = pos
    else:
        argdict["position"] = argdict["position"] - round(argdict["window_size"]/2)
    
    if args.zoom:
        f = np.exp(np.log(4)*args.zoom)
        wsz = int(f*256)
        argdict["window_size"] = wsz
        argdict["stride"] = int(wsz/8)
    
    brws = get_browser(btp, gm, **argdict)
    brws.start(**argdict)
    
    # view = brws._rendered_view
    # with open("./data/browser/browser_views/map_alignment.txt","w") as f:
    #     for line in view:
    #         f.write(Colors.scrub_codes(line) + "\n")
    
    # for name, artist in brws.renderer.artists.items():
    #     print(name, type(artist), artist.params)
    
    # brws.save_state()
    
    # res = input("save browser?")
    # if "y" in res.lower():
    #     brws.save_history()
    

if __name__=="__main__":
    main()

