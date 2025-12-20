

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
DEFAULT_SEQUENCE_PATH = DATA_DIR / cfg.get("sequence_filename","")

other_paths = {
    "repeatmasker_path":  DATA_DIR / "repeatmasker" / "repeats.sorted.bed.gz",
    "dfam_path": DATA_DIR / "dfam/hg38_dfam.nrph.bed.gz"
}

def get_paths():
    return str(DEFAULT_VARIANTS_PATH), str(DEFAULT_GENES_PATH), str(DEFAULT_SEQUENCE_PATH), str(DEFAULT_LIBRARY), other_paths


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