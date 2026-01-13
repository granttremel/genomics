
from typing import TYPE_CHECKING
from pathlib import Path

from ggene.config import DATA_DIR
from ggene.database.annotations import UFeature, AnnotationStream, TabularStream, ColumnSpec, DerivedSpec
from ggene.database.sequences import FASTAStream
from ggene.database.library import prepare_fasta_file, prepare_features_file

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager

class NoncodeV6Stream(FASTAStream):
    """
    NONCODEv6
    http://www.noncode.org/download.php

    """
    
    def __init__(self, filepath):
        super().__init__(filepath)
    

class GSEStream(TabularStream):
    """
    deepbase3
    https://rna.sysu.edu.cn/deepbase3/subpages/download.php

    """
    chr_format = "chr{chrstr}"
    delimiter = '\t'
    
    columns = [
        "chrom",
        "start",
        "end",
        "id",
        ColumnSpec("feature_subtype", str),
        "name",
        "strand",
        ColumnSpec("mean_exp", float),
        ColumnSpec("total_exp", float),
        ColumnSpec("experiments_number", float),
        ColumnSpec("expressed_number", float),
    ]
    
class AllLncStream(TabularStream):
    """
    deepbase3
    https://rna.sysu.edu.cn/deepbase3/subpages/download.php

    """
    chr_format = "chr{chrstr}"
    delimiter = '\t'
    
    columns = [
        "chrom",
        "start",
        "end",
        ColumnSpec("info1", str, formatter = "|", temp=True),
        ColumnSpec("number1", int),
        "strand",
        ColumnSpec("match_start", int),
        ColumnSpec("match_end", int),
        ColumnSpec("read_types", str, formatter = lambda s:[int(ss) for ss in s.split(",")]),
        ColumnSpec("num_records", int),
        ColumnSpec("data1", str, formatter = ',', temp=True),
        ColumnSpec("data2", str, formatter = ',', temp=True),
    ]
    
    derived = [          
        DerivedSpec("name", lambda f: f.info1[0]),
        DerivedSpec("name_tag", lambda f: f.info1[0].removeprefix("hsa-lncRNA")),
        DerivedSpec("parent_feature", lambda f: (f.info1[1] if f.info1[1] != "-" else "")),
        DerivedSpec("info", lambda f: (f.info1[2] if f.info1[2] != "-" else "")),
        DerivedSpec("data", lambda f: list((int(a), int(b)) for a, b in zip(f.data1, f.data2) if a and b), verbose = True),
    ]

class PirBaseTabStream(TabularStream):
    """
    piRBase
    http://bigdata.ibp.ac.cn/piRBase/download.php
    """
    chr_format = "{chrstr}"
    feature_type = 'piRNA_hit'
    
    columns = [
        "chrom",
        "start",
        "end",
        "name",
        ColumnSpec("num", int),
        "strand"
    ]
    

class PirBaseSeqStream(FASTAStream):
    """
    piRBase
    http://bigdata.ibp.ac.cn/piRBase/download.php
    """
    
    def __init__(self, filepath):
        super().__init__(filepath)
    


def get_stream_class(filepath):
    
    if "GSE" in filepath.name:
        return GSEStream, None
    elif "allLnc" in filepath.name:
        return AllLncStream
    elif "NONCODEv6" in filepath.name:
        return NoncodeV6Stream, None
    elif "piRBase" in filepath.parent.name:
        return PirBaseTabStream, PirBaseSeqStream
    else:
        return TabularStream, FASTAStream

def load_database(filepath, source_name = "", feature_type = "", do_scrub = False):
    bed_cls, fa_cls = get_stream_class(filepath)
    
    if not source_name:
        source_name = filepath.stem.split(".")[0]
    
    if not filepath.with_suffix(".gz.tbi").exists():
        if ".xls" in filepath.suffixes:
            do_trim = True
        else:
            do_trim= False
        prepare_features_file(filepath, do_trim = do_trim, do_scrub = do_scrub)
    
    if ".fa" in filepath.suffixes:
        return fa_cls(filepath)
    else:
        return bed_cls(filepath, source_name = source_name, feature_type = feature_type)
    
def bind_database(gm:'GenomeManager', filepath, source_name = "", feature_type = ""):
    
    db = load_database(filepath, source_name=source_name, feature_type = feature_type)
    
    if isinstance(db, TabularStream):
        gm.annotations.add_source(db.source_name, db)
    # elif isinstance(db, FASTAStream):
    #     gm.annotations.
    
    return db

def get_entry(filepath, all_data = False):
    
    fn = filepath.stem.split(".")[0]
    if fn == 'dir':
        return []
    
    sfxs = filepath.suffixes
    if '.tbi' in sfxs or '.gzi' in sfxs or '.fai' in sfxs:
        return
    
    if '.gz' in sfxs:
        has_tbi = filepath.with_suffix('.gz.tbi').exists()
        return fn, filepath, has_tbi
        
    elif all_data:
        return fn, filepath, False
    else:
        return

def list_databases(parent_dir = "", all_data = False, max_items = 20, max_depth = 3, depth = 0):
    
    avoid = ["browser","etc","jaspar","idk","library","outputs","tmp","ignore"]
    
    if depth > max_depth:
        return []
    
    if not parent_dir:
        parent_dir = DATA_DIR
    else:
        parent_dir = Path(parent_dir)
    
    outfs = []

    ni = 0
    for f in parent_dir.iterdir():
        
        if f.is_dir():
            if f.stem in avoid:
                continue 
            entries = list_databases(parent_dir = f, all_data = all_data, max_items = max_items, max_depth = max_depth, depth = depth + 1)
            outfs.extend(entries)
        else:
            entry = get_entry(f)
            if entry:
                outfs.append(entry)
                ni += 1
        
        if ni >= max_items:
            break

    return outfs

def print_databases(parent_dir = "", all_data = False, max_items = 20):

    dbs = list_databases(parent_dir = parent_dir, all_data = all_data, max_items = max_items)

    # Group databases by parent directory
    grouped = {}
    for name, filepath, has_index in dbs:
        parent = filepath.parent
        if parent not in grouped:
            grouped[parent] = []
        grouped[parent].append((name, filepath, has_index))

    print(f"  {"i":2s} {"Name":24s}{"Filename":32s}{"Indexed":16s}")
    # Print grouped by directory
    ndb = 0
    for parent_path in sorted(grouped.keys()):
        print(f"{parent_path}:")
        for name, filepath, has_index in grouped[parent_path]:
            index_marker = "[indexed]" if has_index else ""
            print(f"  {ndb:2d} {name:24s}{filepath.name:32s}{index_marker:16s}")
            ndb += 1
        print()
    return dbs

def find_databases(name = "", in_name = "",has_suffix = "",  parent_dir = "", all_data =  False, max_items = 20):
    
    out_dbs = []
    dbs = list_databases(all_data = all_data, max_items = max_items)
    for fn, filepath, has_index in dbs:
        
        if parent_dir and str(filepath.parent) != str(parent_dir) and str(filepath.parent.name) != str(parent_dir):
            # print(f"{str(filepath.parent)}, {str(filepath.parent.name)} vs {parent_dir}")
            continue
        
        if name and name != fn:
            continue
        
        if in_name and in_name not in fn:
            continue
        
        if has_suffix and has_suffix not in filepath.suffixes:
            continue
        
        out_dbs.append((fn, filepath, has_index))
    
    if not out_dbs:
        print("did not identify database")
        return []
    else:
        return out_dbs


def find_database(name = "", in_name = "",has_suffix = "",  parent_dir = "", all_data =  False, max_items = 20):
    dbs = find_databases(name=name, in_name = in_name,has_suffix=has_suffix, parent_dir = parent_dir, all_data = all_data, max_items = max_items)
    return dbs[0]

def find_and_load_database(gm = None, source_name = "", feature_type = "", name = "", in_name = "", has_suffix = "", parent_dir = "", all_data = False, max_items = 20):
    db = find_database(name=name, in_name = in_name, has_suffix=has_suffix, parent_dir = parent_dir, all_data = all_data, max_items = max_items)
    if not db:
        return None
    else:
        fn, filepath, _ = db
        
        if not source_name:
            source_name = fn
        
        if gm:
            return bind_database(gm, filepath, source_name = source_name, feature_type = feature_type)
        else:
            return load_database(filepath, source_name = source_name, feature_type = feature_type)

def find_directory(in_dir, top_dir = ""):
    
    if not top_dir:
        init = DATA_DIR
    else:
        init = Path(top_dir)
    
    outdir = ""
    for ff in init.iterdir():
        
        if in_dir in ff.stem:
            outdir = ff
            break
        
        elif ff.is_dir():
            fff = find_directory(in_dir, top_dir = ff)
            if fff:
                outdir = fff
                break
        else:
            pass
    
    return outdir
