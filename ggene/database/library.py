"""
Library building utilities for genomic data.

This module provides functions for building and populating GenomeLibrary projects.
"""

import os
import numpy as np

from ggene.genome import ncrna
from ggene.config import DEFAULT_LIBRARY


def build_rna_library(gm, chromes=[], project_name="ncrna_library"):
    """
    Build an ncRNA library using the GenomeLibrary system.

    Creates indexed BED files for each RNA type with sequences,
    merged across all chromosomes into single files.

    Args:
        gm: GenomeManager instance
        chromes: List of chromosomes to process (empty = all)
        project_name: Name for the library project

    Returns:
        GenomeLibrary instance with populated datasets
    """
    from ggene.database.genome_library import GenomeLibrary

    rna_types = ncrna.rna_types
    lib = GenomeLibrary(project_name)

    gstr = gm.annotations.streams.get("genes")

    # Create a writer for each RNA type using factory pattern
    # (avoids contextmanager issues with dict-based batch writing)
    writers = {}
    writer_names = {}
    for rtp in rna_types:
        name = f"ncrna_{rtp}"
        writers[rtp] = lib.create_features_writer(
            name=name,
            columns=["chrom", "start", "end", "gene_id", "name", "strand", "seq"],
            feature_type=rtp,
            organism="homo_sapiens",
            source="ensembl",
            description=f"{rtp} genes with sequences",
            overwrite=True,
        )
        writer_names[rtp] = name

    ns = 0
    try:
        for chrom in gm.iter_chromes():
            if chromes and chrom not in chromes:
                continue

            print(f"Processing chromosome {chrom}")

            for f in gstr.stream(chrom, start=1):
                if f.gene_biotype not in rna_types:
                    continue

                fseq = gm.get_sequence(f.chrom, f.start, f.end)

                # Use gene_id as fallback if name is empty/None
                feat_name = f.name if f.name else f.gene_id
                if not feat_name:
                    feat_name = f"{f.chrom}_{f.start}_{f.end}"

                # Write directly to the appropriate writer
                writers[f.gene_biotype].write({
                    "chrom": f.chrom,
                    "start": f.start,
                    "end": f.end,
                    "gene_id": f.gene_id,
                    "name": feat_name,
                    "strand": f.strand,
                    "seq": fseq,
                })

                ns += 1
                if ns % 1000 == 0:
                    print(f"  collected {ns} sequences")

    finally:
        # Finalize all writers (closes, sorts, compresses, indexes, updates manifest)
        for rtp, writer in writers.items():
            print(f"Finalizing {writer_names[rtp]}...")
            lib.finalize_features_writer(writer_names[rtp], writer)

    print(f"\nDone! Built library with {ns} total sequences")
    print(lib.describe())

    return lib


def build_rna_fasta_library(gm, chromes=[], project_name="ncrna_sequences", overwrite = False):
    """
    Build an ncRNA FASTA library using the GenomeLibrary system.

    Creates indexed FASTA files for each RNA type.

    Args:
        gm: GenomeManager instance
        chromes: List of chromosomes to process (empty = all)
        project_name: Name for the library project

    Returns:
        GenomeLibrary instance with populated datasets
    """
    from ggene.database.genome_library import GenomeLibrary

    rna_types = ncrna.rna_types
    lib = GenomeLibrary(project_name)

    gstr = gm.annotations.streams.get("genes")

    # Create a FASTA writer for each RNA type using factory pattern
    # (avoids contextmanager issues with dict-based batch writing)
    writers = {}
    writer_names = {}
    for rtp in rna_types:
        name = f"ncrna_{rtp}"
        writers[rtp] = lib.create_fasta_writer(
            name=name,
            organism="homo_sapiens",
            source="ensembl",
            description=f"{rtp} sequences",
            overwrite=overwrite,
        )
        writer_names[rtp] = name

    ns = 0
    try:
        for chrom in gm.iter_chromes():
            if chromes and chrom not in chromes:
                continue

            print(f"Processing chromosome {chrom}")

            for f in gstr.stream(chrom, start=1):
                if f.gene_biotype not in rna_types:
                    continue

                fseq = gm.get_sequence(f.chrom, f.start, f.end)
                desc = f"{f.chrom}:{f.start}-{f.end}({f.strand}) {f.gene_id}"

                # Use gene_id as fallback if name is empty/None
                seq_id = f.name if f.name else f.gene_id
                if not seq_id:
                    seq_id = f"{f.chrom}_{f.start}_{f.end}"

                writer = writers[f.gene_biotype]
                writer.write(seq_id, fseq, description=desc)

                ns += 1
                if ns % 1000 == 0:
                    print(f"  collected {ns} sequences")

    finally:
        # Finalize all writers (closes, compresses, indexes, updates manifest)
        for rtp, writer in writers.items():
            print(f"Finalizing {writer_names[rtp]}...")
            lib.finalize_fasta_writer(writer_names[rtp], writer)

    print(f"\nDone! Built FASTA library with {ns} total sequences")
    print(lib.describe())

    return lib


# =============================================================================
# Legacy functions (kept for backwards compatibility)
# =============================================================================

def build_rna_library_legacy(gm, chromes=[]):
    """
    Original implementation - writes TSV per chromosome per RNA type.

    DEPRECATED: Use build_rna_library() instead for cleaner file organization.
    """
    rna_types = ncrna.rna_types

    fdir = DEFAULT_LIBRARY / "ncache" / "ncrnas"

    gstr = gm.annotations.streams.get("genes")

    lib = {rtp: [] for rtp in rna_types}

    ns = 0
    for chrom in gm.iter_chromes():

        if chromes and chrom not in chromes:
            continue

        print(f"starting chromosome {chrom}")

        clib = {rtp: [] for rtp in rna_types}

        for f in gstr.stream(chrom, start=1):

            if f.gene_biotype not in rna_types:
                continue

            fseq = gm.get_sequence(f.chrom, f.start, f.end)
            f.attributes["seq"] = fseq

            clib[f.gene_biotype].append(f)

            ns += 1
            if ns % 1000 == 0:
                print(f"collected {ns} sequences")

        _write_rna_library_legacy(fdir, clib, chrom)

        for rnt, rnas in clib.items():
            lib[rnt].extend(rnas)

    return lib


def _write_rna_library_legacy(fdir, lib, chrom):
    """Legacy writer - creates per-chromosome TSV files."""
    atts = ["chrom", "start", "end", "gene_id", "name", "strand"]
    header = atts + ["seq"]

    for rnt, rnas in lib.items():
        outdir = fdir / rnt / f"chr{chrom}.tsv"
        if not outdir.parent.exists():
            outdir.parent.mkdir(parents=True, exist_ok=True)

        lines = ["\t".join(header)]
        for rna in rnas:

            line = []
            for att in atts:
                v = getattr(rna, att)

                if isinstance(v, float):
                    vstr = f"{v:0.3f}"
                else:
                    vstr = str(v)

                line.append(vstr)

            line.append(rna.attributes.get("seq", ""))

            lines.append("\t".join(line))

        with open(outdir, "w") as f:
            for line in lines:
                f.write(line + "\n")

        sz = os.path.getsize(outdir)
        nbs = int(np.log10(max(1, sz)) / 3)
        szrem = int(sz / max(1, 10**nbs))

        print(f"wrote {szrem}x10^{nbs} bytes for rna type {rnt}, chrom {chrom} to {outdir}")
        