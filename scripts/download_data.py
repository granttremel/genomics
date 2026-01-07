

import argparse

from ggene.config import DATA_DIR, DEFAULT_LIBRARY
from ggene.database.download import DatabaseDownloader




def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism","-o",type=str, default="Homo sapiens", help = "organism to query")
    parser.add_argument("--source","-s", type=str, default="ncbi_genome", help = f"data source, options = {", ".join(DatabaseDownloader.cmd_formats.keys())}")
    parser.add_argument("--assembly-level", "-a", type=str, default = "", help = "filter by assembly level (e.g., 'Complete Genome', 'Chromosome')")
    parser.add_argument("--test","-t", action="store_true", help = "run in test mode")
    parser.add_argument("--dir","-d", type = str, default = "", help = "directory to create within data directory")
    
    parser.add_argument("--list","-l", action = "store_true", help = "use to list data available from data source")
    parser.add_argument("--full", "-f", action="store_true", help = "show all available fields in the assembly summary")
    
    args = parser.parse_args()
    
    dbd = DatabaseDownloader()
    
    print(args)
    
    if args.list:
        
        org= args.organism
        org = org.replace("_"," ")
        org = org.capitalize()
        print(f"searching for organism {org}")
        
        assemblies = dbd.list_available_assemblies(
            org,
            args.assembly_level if args.assembly_level else None,
            limit=100,
            print_as_found=True,
            show_all_fields=args.full
        )
    else:
        res = dbd.download_organism(args.organism, local_dir=args.dir, source=args.source, test=args.test)
        print(f"database downloaded with result {res}")

if __name__=="__main__":
    main()
