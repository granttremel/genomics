
import argparse
import numpy as np

from ggene.database.genome_manager import GenomeManager
from ggene.config import get_paths, DATA_DIR
from ggene.browser.builds.test_browser import TestBrowser, TestBrowserState
from ggene.browser.builds.justlines_browser import JustLinesBrowser, JustLinesBrowserState
from ggene.browser.builds.scalar_browser import ScalarBrowser, ScalarBrowserState
from ggene.browser.builds.seqs_browser import SeqsBrowser, SeqsBrowserState
from ggene.browser.builds.exp_browser import ExpBrowser, ExpBrowserState
from ggene.database.annotations import get_experiment_stream

exp_path = DATA_DIR / "rna_seq" / "deepbase3_lncRNA_human" / "GSE43335.xls.gz"

def get_browser(browser_type:str, gm, **kwargs):
    
    if browser_type == "test":
        return TestBrowser(gm=gm, **kwargs)
    elif browser_type == "scalar":
        return ScalarBrowser(gm=gm, **kwargs)
    elif browser_type == "seqs":
        return SeqsBrowser(gm=gm, **kwargs)
    elif browser_type == "lines":
        return JustLinesBrowser(gm=gm, **kwargs)
    elif browser_type == "exp":
        exp_path = kwargs.pop("file_path", "")
        return ExpBrowser(exp_path, gm=gm, min_exp = 0.05, **kwargs)
    else:
        return None

def load_genome():
    
    a,b,c, lib_path, other_paths = get_paths()
    gm = GenomeManager(**other_paths)
    
    
    return gm

def test_exp(gm):
    
    exp_path = DATA_DIR / "rna_seq" / "deepbase3_lncRNA_human" / "GSE43335.xls.gz"
    
    strm = get_experiment_stream(exp_path, "deepbase3","exp_read", gm=gm)
    
    for f in strm.stream("2", 1e6, None):
        print(f)


def main():
    
    gm = load_genome()
    
    # test_exp(gm)
    # return
    
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
    parser.add_argument("--load", "-l", default = False, action = "store_true", 
                        help="Load last state")
    parser.add_argument("--debug", "-d", default = False, action = "store_true", 
                        help="Debug mode")
    parser.add_argument("--browser-type", "-b", default = "test", type = str, 
                        help="Browser type (test, scalar, lines, seqs, exp)")
    
    parser.add_argument("--file", "-f", default = None, type = str, dest = "file_path",
                        help = "path to experiment file to load into browser")
    
    display_width = 256
    
    args = parser.parse_args()
    
    argdict = args.__dict__
    argdict["display_width"] = display_width
    btp = argdict.pop("browser_type","")
    
    if not argdict.get("file_path",""):
        argdict["file_path"] = exp_path
    
    if args.random:
        chrom, pos = GenomeManager.get_random_location()
        argdict["chrom"] = chrom
        argdict["position"] = pos
    
    if args.zoom:
        
        f = np.exp(np.log(4)*args.zoom)
        wsz = int(f*256)
        argdict["window_size"] = wsz
        argdict["stride"] = int(wsz/8)
    
    if args.load:
        pass
    
    argdict["feature_type"] = "★exp_hit★"
    
    brws = get_browser(btp, gm, **argdict)
    brws.start(**argdict)
    
    # exp_stream = brws.exp_stream
    
    # for f in exp_stream.stream("4", start = 15220000, end = None):
    #     print(f)
    
    # brws.save_state()
    
    # res = input("save browser?")
    # if "y" in res.lower():
    #     brws.save_history()
    

if __name__=="__main__":
    main()

