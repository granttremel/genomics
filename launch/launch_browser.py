
import argparse

from ggene.database.genome_manager import GenomeManager

from ggene.browser.builds.test_browser import TestBrowser, TestBrowserState
from ggene.browser.builds.justlines_browser import JustLinesBrowser, JustLinesBrowserState
from ggene.browser.builds.scalar_browser import ScalarBrowser, ScalarBrowserState
from ggene.browser.builds.seqs_browser import SeqsBrowser, SeqsBrowserState

def get_browser(browser_type:str, **kwargs):
    
    if browser_type == "test":
        return TestBrowser(**kwargs)
    elif browser_type == "scalar":
        return ScalarBrowser(**kwargs)
    elif browser_type == "seqs":
        return SeqsBrowser(**kwargs)
    elif browser_type == "lines":
        return JustLinesBrowser(**kwargs)
    else:
        return None

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
    parser.add_argument("--random", "-r", default = False, action = "store_true", 
                        help="Choose random chromosome and position")
    parser.add_argument("--load", "-l", default = False, action = "store_true", 
                        help="Load last state")
    parser.add_argument("--debug", "-d", default = False, action = "store_true", 
                        help="Debug mode")
    parser.add_argument("--browser-type", "-b", default = "test", type = str, 
                        help="Browser type (test, scalar, lines, seqs)")
    
    display_width = 256
    
    args = parser.parse_args()
    
    argdict = args.__dict__
    argdict["display_width"] = display_width
    btp = argdict.pop("browser_type","")
    
    if args.random:
        chrom, pos = GenomeManager.get_random_location()
        argdict["chrom"] = chrom
        argdict["position"] = pos
    
    if args.load:
        pass
    
    brws = get_browser(btp, **argdict)
    brws.start(**argdict)
    
    # brws.save_state()
    
    # res = input("save browser?")
    # if "y" in res.lower():
    #     brws.save_history()
    

if __name__=="__main__":
    main()

