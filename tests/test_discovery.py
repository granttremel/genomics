#!/usr/bin/env python3
"""Quick test of assembly discovery"""

from ggene.database.download import DatabaseDownloader

dbd = DatabaseDownloader()

# Find choanoflagellate assemblies
print("=== Searching for choanoflagellate assemblies ===")
assemblies = dbd.list_available_assemblies(organism_name="choanoflagellate", limit=20)
dbd.print_assemblies(assemblies)

# Find Salpingoeca specifically
print("\n=== Searching for Salpingoeca assemblies ===")
salpingoeca = dbd.list_available_assemblies(organism_name="Salpingoeca", limit=10)
dbd.print_assemblies(salpingoeca)

# Get info for the assembly you found
print("\n=== Getting info for GCF_000188695.1 ===")
info = dbd.get_assembly_info("GCF_000188695.1")
if info:
    print(f"RefSeq ID: {info['refseq_id']}")
    print(f"Organism: {info['organism_name']}")
    print(f"Strain: {info['strain']}")
    print(f"Level: {info['assembly_level']}")
    print(f"FTP Path: {info['ftp_path']}")
