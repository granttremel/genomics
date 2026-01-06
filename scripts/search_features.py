

import argparse

from ggene.database.genome_manager import GenomeManager

def load_genome():
    gm = GenomeManager()
    return gm

def search_features(gm, chromes, position, search_str, num_features, feature_types, show_long = False, cond = None):
    
    search_str = search_str.lower()
    catt, cv = cond if cond else (None,None)
    
    for chrom in chromes:
        print(f"searching chrom {chrom}")
        
        nf = 0
        # for f in gm.annotations.stream_all(chrom, position):
        for f in gm.annotations.streams.get("genes",None).stream(chrom, start=position, end = None):
            
            if not f.feature_type in feature_types:
                continue
            
            if catt:
                if hasattr(f, catt) and not getattr(f, catt) == cv:
                    continue
            
            good = False
            
            fd = f.to_dict()
            fd.update(fd.pop("attributes",{}))
            for k, v in fd.items():
                if isinstance(v, str) and search_str in v.lower():
                    good = True
                    break
                    
            if not good:
                continue
            
            if show_long:
                print(format_feature_long(f))
            else:
                print(format_feature_short(f))
            
            nf += 1
            
            if num_features > 0 and nf > num_features:
                break
            
            if nf % 64 == 0:
                input("press enter to view more features")

def get_feature_parts(f, part_fmt = "{att}={v}"):
    
    fparts = []
    fd  = f.to_dict()
    for att, v in fd.items():
        
        # if isinstance(v, dict) or att in ["name","type"]:
        #     continue
        
        fparts.append(part_fmt.format(att=att,v=v))
    
    attparts = []
    for att, v in f.attributes.items():
        attparts.append(part_fmt.format(att=att,v=v))
    
    return fparts, attparts

def format_feature_short(f):
    
    fparts, attparts = get_feature_parts(f)
    
    return f"{f.name}, {f.feature_type}: {", ".join(fparts)}, attributes=[{", ".join(attparts)}]"

def format_feature_long(f):
    
    fparts, attparts = get_feature_parts(f, part_fmt = "{att}: {v}")
    
    out = [f"{f.name}, {f.feature_type}"]
    for part in fparts:
        out.append("  " + part)
    out.append("  attributes:")
    for part in attparts:
        out.append("    " + part)
    
    return "\n".join(out)
    

def main():
    
    parser = argparse.ArgumentParser(description='Interactive Genome Browser')
    parser.add_argument('--chrom', '-c', type=str, default='',
                       help='Chromosome to browse')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--search', '-s', type=str, default="",
                       help='Substring to search for')
    parser.add_argument('--num', '-n', type=int, default=-1,
                       help='number of features')
    parser.add_argument('--features', '-f', type=str, default="",
                       help='feature types, comma delimited (e.g. "genes,exons,repeats")')
    parser.add_argument('--condition', '-cond', type=str, default="",
                       help='condition, e.g. gene_biotype=protein_coding')
    parser.add_argument('--long', '-l', action="store_true",
                       help='flag for long format')
    
    args = parser.parse_args()
    
    features_raw = args.features
    if features_raw:
        features = [f.strip() for f in features_raw.split(',')]
    else:
        features = ["gene",'lncRNA','ncRNA',"pseudogene","exon"]
    
    cond_raw = args.condition
    if cond_raw:
        cond = tuple([c.strip() for c in cond_raw.split("=")])
    else:
        cond = None
    
    gm = load_genome()
    
    chrom, pos = args.chrom, args.position    
    if not chrom:
        chromes = [chrom for chrom in gm.iter_chromes()]
    else:
        chromes = [chrom]
    
    search_features(gm, chromes, pos, args.search, args.num, features, show_long = args.long, cond=cond)


if __name__=="__main__":
    main()