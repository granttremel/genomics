
import argparse
from pathlib import Path
import subprocess
import urllib.request
from typing import List, Dict, Optional

from ggene.config import DATA_DIR, DEFAULT_LIBRARY

refseq_ids = {
    "h_sapiens":"",
    "s_rosetta":"GCF_000188695.1",
}

organism = {
    "h_sapiens":"",
    "s_rosetta":"Proterospongia", # what's going on here
    
}

strain_ids={
    "h_sapiens":"",
    "s_rosetta":"ATCC50818",
    
}


class DatabaseDownloader:
    
    path_formats = {
        "ncbi_genome":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/695/{refseq_id}_{organism_name}_sp_{strain}/",
    }
    
    # cmds_formats = {
    #     # "ncbi_genome":"lftp -c mirror https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/695/{refseq_id}_{organism_name}_sp_{strain}/ {local_path}",
    # }
    cmd_formats = {
        
        "ncbi_genome":"lftp -c mirror {ftp_path} {local_path}",
    }
    
    
    def __init__(self, data_dir = DATA_DIR, library_dir = DEFAULT_LIBRARY):
        
        self.data_dir = Path(data_dir)
        self.library_dir = Path(library_dir)
    
    def get_ftp_path(self, source, cmd_data):
        
        if "ftp_path" in cmd_data:
            return cmd_data["ftp_path"]
        else:
            path_frm = self.path_formats[source]
            return path_frm.format(**cmd_data)
    
    def format_command(self, data_source, local_path, org_data):
        
        if not data_source in self.cmd_formats:
            print(f"data source {data_source} not found. options are: {", ".join(self.cmds.keys())}")
        
        ftp_path = self.get_ftp_path(data_source, org_data)
        
        cmd = self.cmd_formats[data_source]
        cmd = cmd.format(ftp_path = ftp_path, local_path = local_path)
        
        return cmd.split()
        
    def download_organism(self, org:str, local_dir = "", source = "ncbi_genome", test=False):
        
        org_data = self.get_organism_data(org, allow_search = True)
        
        if not org_data:
            return False
        
        if not local_dir:
            org_name = org_data.get("organism_name","")
            if not org_name:
                local_dir = org
            else:
                org_name = org_name.replace(org_data.get("strain", ""), "").strip()
                on_parts= org_name.split()
                init = on_parts[0][0].upper()
                org_name = "_".join([init] + on_parts[1:])
                local_dir = org_name
        
        full_path = self.data_dir / local_dir
        if not full_path.exists():
            full_path.mkdir(exist_ok = True)
        
        cmd = self.format_command(source, full_path.absolute(), org_data) 
        
        if test:
            print(f"(test) running command {cmd} ... ")
            s = True
        else:
            print(f"running command {cmd} ... ")
            result = subprocess.run(cmd, text=True)
            
            if result.returncode != 0:
                print(f"command failed with code {result.returncode}: {result.stderr}")
        
            s = bool(full_path.iterdir())
        
        return s
       
        
    
    def get_organism_data(self, org, allow_search = False):

        if not org in organism:
            
            if allow_search:
                
                assy_list = self.list_available_assemblies(org, limit = 5)
                
                for i, assy in enumerate(assy_list):
                    print(i, ", ".join([f"{k}: {v}" for k,v in assy.items()]))
                    
                res = input("select assembly index: ")
                
                try:
                    ind = int(res)
                    return assy_list[ind]
                except:
                    print(f"failed to parse input {res}.")
                    return {}
                
            else:
                print(f"organism {org} not found. options are: {", ".join(organism.keys())}")

        data = {
            "refseq_id":refseq_ids.get(org,""),
            "organism_name":organism.get(org, ""),
            "strain":strain_ids.get(org, "")
        }

        return data

    def list_available_assemblies(
        self,
        organism_name: Optional[str] = None,
        assembly_level: Optional[str] = None,
        limit: int = 100,
        print_as_found: bool = False,
        show_all_fields: bool = False
    ) -> List[Dict[str, str]]:
        """
        Query NCBI's assembly summary to find available genomes.
        Streams the file line-by-line for memory efficiency.

        Args:
            organism_name: Filter by organism name (case-insensitive, partial match)
            assembly_level: Filter by level (e.g., 'Complete Genome', 'Chromosome', 'Scaffold')
            limit: Maximum number of results to return
            print_as_found: Print results as they're discovered

        Returns:
            List of dictionaries with assembly metadata
        """
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

        print(f"Streaming assembly summary from NCBI...")

        try:
            with urllib.request.urlopen(url) as response:
                # Stream line by line
                headers = None
                results = []
                printed_headers = False

                for line in response:
                    line = line.decode('utf-8').strip()

                    # Find header line
                    if line.startswith('#assembly_accession'):
                        headers = line.strip("#").split('\t')

                        # Print headers if showing all fields
                        if print_as_found and show_all_fields and not printed_headers:
                            print(f"\nAvailable fields ({len(headers)} total):")
                            for i, h in enumerate(headers, 1):
                                print(f"  {i:2}. {h}")
                            print("\n" + "=" * 150)
                            printed_headers = True
                        elif print_as_found and not printed_headers:
                            print(f"\n{'RefSeq ID':<20} {'Organism':<40} {'Strain':<20} {'Level':<20}")
                            print("-" * 100)
                            printed_headers = True

                        continue

                    # Skip other comments
                    if line.startswith('#') or not line:
                        continue

                    # Need headers before processing data
                    if headers is None:
                        continue

                    # Parse TSV row
                    fields = line.split('\t')
                    if len(fields) < len(headers):
                        continue

                    row = dict(zip(headers, fields))

                    # Apply filters
                    if organism_name:
                        org_name = row.get('organism_name', '').lower()
                        if organism_name.lower() not in org_name:
                            continue

                    if assembly_level:
                        if row.get('assembly_level', '') != assembly_level:
                            continue

                    # Store result - either all fields or subset
                    if show_all_fields:
                        result = row.copy()
                        result['strain'] = self._extract_strain(row.get('infraspecific_name', ''))
                    else:
                        result = {
                            'refseq_id': row.get('assembly_accession', ''),
                            'organism_name': row.get('organism_name', ''),
                            'infraspecific_name': row.get('infraspecific_name', ''),
                            'assembly_level': row.get('assembly_level', ''),
                            'ftp_path': row.get('ftp_path', ''),
                            'strain': self._extract_strain(row.get('infraspecific_name', '')),
                        }

                    results.append(result)

                    # Print immediately if requested
                    if print_as_found:
                        if show_all_fields:
                            print(f"\nResult {len(results)}:")
                            for key, value in result.items():
                                if value:  # Only show non-empty fields
                                    print(f"  {key:<30}: {value}")
                        else:
                            print(f"{result['refseq_id']:<20} {result['organism_name']:<40} {result['strain']:<20} {result['assembly_level']:<20}")

                    # Stop when limit reached
                    if len(results) >= limit:
                        break

                if print_as_found and not show_all_fields:
                    print()

                print(f"Found {len(results)} matching assemblies")
                return results

        except Exception as e:
            print(f"Error fetching assembly data: {e}")
            return []

    def _extract_strain(self, infraspecific_name: str) -> str:
        """Extract strain ID from infraspecific_name field (e.g., 'strain=ATCC50818')"""
        if not infraspecific_name:
            return ""

        for part in infraspecific_name.split(';'):
            part = part.strip()
            if part.startswith('strain='):
                return part.split('=', 1)[1]

        return ""

    def print_assemblies(self, assemblies: List[Dict[str, str]]):
        """Pretty print assembly information"""
        if not assemblies:
            print("No assemblies found")
            return

        print(f"\n{'RefSeq ID':<20} {'Organism':<40} {'Strain':<20} {'Level':<20}")
        print("-" * 100)
        for asm in assemblies:
            print(f"{asm['refseq_id']:<20} {asm['organism_name']:<40} {asm['strain']:<20} {asm['assembly_level']:<20}")
        print()

    def get_assembly_info(self, refseq_id: str) -> Optional[Dict[str, str]]:
        """Get detailed info for a specific RefSeq assembly"""
        assemblies = self.list_available_assemblies(limit=1000000)

        for asm in assemblies:
            if asm['refseq_id'] == refseq_id:
                return asm

        return None


# def main():
    
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--organism","-o",type=str, default="Homo sapiens", help = "organism to query during list operations")
#     # parser.add_argument("--species", "-s", type=str, help="species to download")
#     parser.add_argument("--source","-s", type=str, default="ncbi_genome", help = f"data source, options = {", ".join(DatabaseDownloader.cmd_formats.keys())}")
#     parser.add_argument("--test","-t", action="store_true", help = "run in test mode")
#     parser.add_argument("--dir","-d", type = str, default = "", help = "directory to create within data directory")
    
#     parser.add_argument("--list","-l", action = "store_true", help = "use to list data available from data source")
#     parser.add_argument("--assembly-level", "-a", type=str, default = "", help = "filter by assembly level (e.g., 'Complete Genome', 'Chromosome')")
#     parser.add_argument("--full", "-f", action="store_true", help = "show all available fields in the assembly summary")
    
#     args = parser.parse_args()
    
#     dbd = DatabaseDownloader()
    
#     print(args)
    
#     if args.list:
        
#         org= args.organism
#         org = org.replace("_"," ")
#         org = org.capitalize()
#         print(f"searching for organism {org}")
        
#         assemblies = dbd.list_available_assemblies(
#             org,
#             args.assembly_level if args.assembly_level else None,
#             limit=100,
#             print_as_found=True,
#             show_all_fields=args.full
#         )
#     else:
#         res = dbd.download_organism(args.organism, local_dir=args.dir, source=args.source, test=args.test)
#         print(f"database downloaded with result {res}")

# if __name__=="__main__":
#     main()

