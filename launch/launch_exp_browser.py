
import argparse
import numpy as np

from ggene.database.genome_manager import GenomeManager
from ggene.config import get_paths, DATA_DIR
from ggene.browser.builds.test_browser import TestBrowser, TestBrowserState
from ggene.browser.builds.justlines_browser import JustLinesBrowser, JustLinesBrowserState
from ggene.browser.builds.scalar_browser import ScalarBrowser, ScalarBrowserState
from ggene.browser.builds.seqs_browser import SeqsBrowser, SeqsBrowserState
from ggene.browser.builds.exp_browser import ExpBrowser, ExpBrowserState
from ggene.browser.builds.map_browser import MapBrowser, MapBrowserState
from ggene.database.annotations import get_experiment_stream
from ggene.database.databases import find_and_load_database, find_databases, find_directory
from ggene.draw import Colors

def get_browser(gm, **kwargs):
    exp_keys = kwargs.pop("exp_keys", [])
    return ExpBrowser(exp_keys, gm=gm, min_exp = 0.05, **kwargs)

def get_experiment_keys(exps, exp_ptrn, exp_dir=""):
    
    exp_keys = []
    
    if exp_dir:
        exp_dir = find_directory(in_dir = str(exp_dir))
    
    if exp_ptrn:
        exp_paths = find_databases(in_name = str(exp_ptrn), parent_dir = str(exp_dir))
        exp_keys = [str(fn) for fn, fpath, _ in exp_paths]
    
    if exps and not exp_keys:
        exp_keys = exps.split(',')
    
    return exp_keys

def get_experiment_filter(filter_str):
    
    if not filter_str:
        return None
    
    comps = ['<=','>=','==']
    comp_char = fa = fb = ""
    for c in comps:
        if c in filter_str:
            fa, fb = filter_str.split(c)
            comp_char = c
            break
    
    op = None
    if comp_char=='<=':
        op = '__le__'
    elif comp_char=='>=':
        op = '__ge__'
    elif comp_char=='==':
        op = '__eq__'
    
    try:
        fb = float(fb)
    except:
        print(f"failed to parse {fb}")
        return None
    
    if not op or not fa or not fb:
        return None
    
    def filter(f):
        if not f:
            return False
        
        if not hasattr(f, fa):
            return True
        
        fval = getattr(f, fa)
        if not fval:
            return True
        
        return getattr(float, op)(fval, fb)
        
    return filter
    

def load_genome(**kwargs):
    a,b,c, lib_path, other_paths = get_paths()
    kwargs.update(other_paths)
    gm = GenomeManager(**kwargs)
    return gm

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", "-c", default = "1", type = str,
                        help="Chromosome (1-23, X, Y, MT)")
    parser.add_argument("--position", "-p", default = int(1e6), type = int,
                        help="Window size, in bp (default 4096)")
    parser.add_argument("--window", "-w", default = 16000, type = int,
                        help="Window size, in bp (default 16000)", dest = "window_size")
    parser.add_argument("--stride", "-s", default = 1600, type = int, 
                        help="Window size, in bp (default 1600)")
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
    
    parser.add_argument("--exps", "-x", default = "", type = str,
                        help = "experiment database keys (comma separated)")
    parser.add_argument("--exp-ptrn", "-xp", default = "", type = str,
                        help = "experiment database pattern to discover keys (in_name)")
    parser.add_argument("--exp-dir", "-xd", default = "", type = str,
                        help = "Experiment directory to search in")
    parser.add_argument("--filter","-xf", default = "", type = str,
                        help = "Filter to refine features from experiment. form 'att>=val', 'att<=val', 'att==val'")
    
    display_width = 256
    
    args = parser.parse_args()
    
    argdict = args.__dict__
    gm = load_genome(load_jaspars = argdict.pop("jaspars", False), load_patterns = argdict.pop("patterns", False))
    
    argdict["display_width"] = display_width
    
    # exp_keys = get_experiment_keys(argdict.pop("exps",""), argdict.pop("exp_ptrn",""), argdict.pop("exp_dir",""))
    exp_keys = ()
    argdict["exp_keys"] = exp_keys
    
    print(f"enter exp_filt with {argdict.get("filter","")}")
    exp_filt = get_experiment_filter(argdict.pop("filter",""))
    argdict["exp_filter"] = exp_filt
        
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
    
    if args.load:
        pass
    
    # argdict["feature_type"] = "★exp_hit★"
    argdict["feature_type"] = "auto"
    
    brws = get_browser(gm, **argdict)
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

