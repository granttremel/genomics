

import argparse

from ggene.database.genome_manager import GenomeManager

def load_genome():
    gm = GenomeManager()
    return gm

def browse_features(gm, chromes, position, num_features, feature_types, show_long = False, cond = None, length_range = None):
    
    catt, cv = cond if cond else (None,None)
    
    if not length_range:
        min_length, max_length = (0, 1e9)
    else:
        min_length, max_length = length_range
    
    nf = 0
    for chrom in chromes:
        for f in gm.annotations.stream_all(chrom, position):
            
            if not f.feature_type in feature_types:
                continue
            
            if f.length < min_length or f.length > max_length:
                continue
            
            if catt:
                if hasattr(f, catt):
                    v = getattr(f, catt)
                    if isinstance(v, str):
                        res = cv in v
                    else:
                        res = cv==v
                    
                    if not res:
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
    
    atts = fd.pop("attributes", {})
    for att, v in fd.items():
        
        fparts.append(part_fmt.format(att=att,v=v))
    
    attparts = []
    for att, v in atts.items():
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
    parser.add_argument('--chrom', '-c', type=str, default='1',
                       help='Chromosome to browse (default: 1)')
    parser.add_argument('--position', '-p', type=int, default=1000000,
                       help='Starting position (default: 1000000)')
    parser.add_argument('--random', '-r', action="store_true",
                       help='random position', dest = "random")
    parser.add_argument('--num', '-n', type=int, default=-1,
                       help='number of features')
    parser.add_argument('--features', '-f', type=str, default="",
                       help='feature types, comma delimited (e.g. "genes,exons,repeats")')
    parser.add_argument('--condition', '-cond', type=str, default="",
                       help='condition, e.g. gene_biotype=protein_coding')
    parser.add_argument('--verbose', '-v', action="store_true",
                       help='flag for long format')
    
    parser.add_argument('--length', '-l', type = str, default = "", help = "length range, comma separated e.g. 200,300")
    
    args = parser.parse_args()
    
    features_raw = args.features
    if features_raw:
        features = [f.strip() for f in features_raw.split(',')]
    else:
        features = ["gene","exon","variant"]
    
    cond_raw = args.condition
    if cond_raw:
        cond = tuple([c.strip() for c in cond_raw.split("=")])
    else:
        cond = None
    
    gm = load_genome()
    
    if args.random:
        chrom, pos = gm.get_random_location()
    else:
        chrom, pos = args.chrom, args.position
    
    if chrom == 'all':
        chromes = [chrom for chrom in gm.iter_chromes()]
    else:
        chromes = [chrom]
    
    length_range = None
    if args.length:
        mnstr,mxstr= args.length.strip().split(',')
        try:
            min_len = int(mnstr)
            max_len = int(mxstr)
            length_range = (min_len, max_len)
        except:
            pass
        
    browse_features(gm, chromes, pos, args.num, features, show_long = args.verbose, cond=cond, length_range = length_range)


if __name__=="__main__":
    main()