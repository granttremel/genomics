

import os

import yaml
from pathlib import Path

CONFIG_PATH = "./local.yaml" if os.path.exists("./local.yaml") else "./default.yaml"

cfg = {}
with open(CONFIG_PATH) as f:
    cfg = yaml.safe_load(f)
    
# Constants
DATA_DIR = Path(cfg.get("data_dir"), )
DEFAULT_LIBRARY = DATA_DIR / cfg.get("library_dir","")
DEFAULT_VARIANTS_PATH = DATA_DIR / cfg.get("variants_filename","")
DEFAULT_GENES_PATH = DATA_DIR / cfg.get("genes_filename","")
DEFAULT_GENEINFO_PATH = DATA_DIR / cfg.get("geneinfo_filename","")
DEFAULT_SEQUENCE_PATH = DATA_DIR / cfg.get("sequence_filename","")
DEFAULT_CHRDATA_CACHE_PATH = DEFAULT_LIBRARY / "dcache"

all_paths = {k: DATA_DIR / v for k, v in cfg.items() if isinstance(v, str)}
other_paths = {k: DATA_DIR / v for k,v in cfg.get("other_paths", {}).items()}
all_paths.update(other_paths)


def get_paths():
    return str(DEFAULT_VARIANTS_PATH), str(DEFAULT_GENES_PATH), str(DEFAULT_SEQUENCE_PATH), str(DEFAULT_LIBRARY), other_paths

def get_all_paths():
    return all_paths


def get_config(default = False):
    
    if default:
        cfg_path = "./default.yaml"
    else:
        cfg_path = "/local.yaml"
        if not os.path.exists(cfg_path):
            return get_config(default = True)
        
    cfg = {}
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    
    return cfg