"""
Simple gene information database loader

Loads lightweight metadata about genes from public databases:
- NCBI Gene Info (descriptions, symbols)
- GO annotations (function, process, component)
- Expression data (tissue-specific)
- Protein domains (InterPro/Pfam)
"""

from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict
import gzip

from ggene.config import DATA_DIR


@dataclass
class GeneInfo:
    """Basic gene information"""
    gene_id: str
    symbol: str
    description: str
    chromosome: str = ""
    synonyms: List[str] = field(default_factory=list)
    gene_type: str = ""
    go_terms: Dict[str, List[str]] = field(default_factory=lambda: defaultdict(list))
    domains: List[str] = field(default_factory=list)
    expression: Dict[str, float] = field(default_factory=dict)

    def __str__(self):
        return f"{self.symbol}: {self.description[:80]}..."

    def summary(self) -> str:
        """Formatted summary"""
        lines = [
            f"Gene: {self.symbol} ({self.gene_id})",
            f"Description: {self.description}",
        ]

        if self.chromosome:
            lines.append(f"Location: chr{self.chromosome}")

        if self.synonyms:
            lines.append(f"Synonyms: {', '.join(self.synonyms[:5])}")

        if self.go_terms:
            lines.append("\nGO Terms:")
            for category, terms in self.go_terms.items():
                if terms:
                    lines.append(f"  {category}: {', '.join(terms[:3])}")

        if self.domains:
            lines.append(f"\nDomains: {', '.join(self.domains[:5])}")

        if self.expression:
            top_tissues = sorted(self.expression.items(), key=lambda x: -x[1])[:3]
            lines.append("\nTop Expression:")
            for tissue, tpm in top_tissues:
                lines.append(f"  {tissue}: {tpm:.1f} TPM")

        return "\n".join(lines)


class GeneInfoDB:
    """
    Lightweight gene information database

    Loads and queries gene metadata from public databases.
    """

    def __init__(self, data_dir: Optional[Path] = None):
        self.data_dir = data_dir or DATA_DIR / "gene_info"
        self.genes: Dict[str, GeneInfo] = {}
        self.symbol_to_id: Dict[str, str] = {}
        self.ensembl_to_id: Dict[str, str] = {}

    def load_ncbi_gene_info(self, filepath: Optional[Path] = None):
        """
        Load NCBI gene_info file

        Download from:
        ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

        Args:
            filepath: Path to gene_info file (auto-detects .gz)
        """
        if filepath is None:
            filepath = self.data_dir / "Homo_sapiens.gene_info.gz"
            if not filepath.exists():
                filepath = self.data_dir / "Homo_sapiens.gene_info"

        if not filepath.exists():
            print(f"Gene info file not found: {filepath}")
            print("Download from: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
            return

        opener = gzip.open if filepath.suffix == '.gz' else open

        print(f"Loading NCBI gene info from {filepath}...")
        count = 0

        with opener(filepath, 'rt') as f:
            # Skip header
            header = next(f).strip().split('\t')

            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                # Parse fields
                tax_id = fields[0]
                if tax_id != "9606":  # Human only
                    continue

                gene_id = fields[1]
                symbol = fields[2]
                synonyms = fields[4].split('|') if fields[4] != '-' else []
                dbxrefs = fields[5]
                chromosome = fields[6]
                description = fields[8]
                gene_type = fields[9]

                # Create GeneInfo
                gene_info = GeneInfo(
                    gene_id=gene_id,
                    symbol=symbol,
                    description=description,
                    chromosome=chromosome,
                    synonyms=synonyms,
                    gene_type=gene_type,
                )

                self.genes[gene_id] = gene_info
                self.symbol_to_id[symbol.upper()] = gene_id

                # Parse dbxrefs for Ensembl ID
                for xref in dbxrefs.split('|'):
                    if xref.startswith('Ensembl:'):
                        ensembl_id = xref.split(':')[1]
                        self.ensembl_to_id[ensembl_id] = gene_id

                count += 1

        print(f"Loaded {count} genes")

    def load_go_annotations(self, filepath: Optional[Path] = None, max_terms: int = 10):
        """
        Load GO annotations from GAF file

        Download from:
        http://geneontology.org/gene-associations/goa_human.gaf.gz

        Args:
            filepath: Path to GAF file
            max_terms: Max GO terms per category (to keep it lightweight)
        """
        if filepath is None:
            filepath = self.data_dir / "goa_human.gaf.gz"
            if not filepath.exists():
                filepath = self.data_dir / "goa_human.gaf"

        if not filepath.exists():
            print(f"GO annotation file not found: {filepath}")
            return

        opener = gzip.open if filepath.suffix == '.gz' else open

        print(f"Loading GO annotations from {filepath}...")
        count = 0

        with opener(filepath, 'rt') as f:
            for line in f:
                if line.startswith('!'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 14:
                    continue

                db_object_symbol = fields[2]  # Gene symbol
                go_id = fields[4]
                go_term = fields[9]  # GO term name
                aspect = fields[8]  # P, F, or C

                # Map aspect to readable name
                aspect_map = {
                    'P': 'Biological Process',
                    'F': 'Molecular Function',
                    'C': 'Cellular Component',
                }
                aspect_name = aspect_map.get(aspect, 'Other')

                # Find gene
                gene_id = self.symbol_to_id.get(db_object_symbol.upper())
                if gene_id and gene_id in self.genes:
                    gene = self.genes[gene_id]
                    if len(gene.go_terms[aspect_name]) < max_terms:
                        gene.go_terms[aspect_name].append(go_term)
                        count += 1

        print(f"Loaded {count} GO annotations")

    def load_simple_expression(self, filepath: Path, tissue_col: str = "tissue",
                              gene_col: str = "gene", tpm_col: str = "TPM"):
        """
        Load simple expression data (e.g., median TPM per tissue)

        Expects a tab-delimited file with columns for gene, tissue, and TPM.
        You'll need to prepare this from GTEx or similar.

        Args:
            filepath: Path to expression file
            tissue_col: Column name for tissue
            gene_col: Column name for gene (symbol or Ensembl ID)
            tpm_col: Column name for expression value
        """
        if not filepath.exists():
            print(f"Expression file not found: {filepath}")
            return

        print(f"Loading expression data from {filepath}...")
        count = 0

        with open(filepath) as f:
            header = next(f).strip().split('\t')
            tissue_idx = header.index(tissue_col)
            gene_idx = header.index(gene_col)
            tpm_idx = header.index(tpm_col)

            for line in f:
                fields = line.strip().split('\t')
                tissue = fields[tissue_idx]
                gene_name = fields[gene_idx]
                tpm = float(fields[tpm_idx])

                # Find gene
                gene_id = self.symbol_to_id.get(gene_name.upper())
                if not gene_id:
                    gene_id = self.ensembl_to_id.get(gene_name)

                if gene_id and gene_id in self.genes:
                    self.genes[gene_id].expression[tissue] = tpm
                    count += 1

        print(f"Loaded {count} expression values")

    def get_by_symbol(self, symbol: str) -> Optional[GeneInfo]:
        """Get gene info by symbol (case-insensitive)"""
        gene_id = self.symbol_to_id.get(symbol.upper())
        return self.genes.get(gene_id) if gene_id else None

    def get_by_id(self, gene_id: str) -> Optional[GeneInfo]:
        """Get gene info by NCBI gene ID"""
        return self.genes.get(gene_id)

    def get_by_ensembl(self, ensembl_id: str) -> Optional[GeneInfo]:
        """Get gene info by Ensembl ID"""
        gene_id = self.ensembl_to_id.get(ensembl_id)
        return self.genes.get(gene_id) if gene_id else None

    def search(self, query: str, limit: int = 10) -> List[GeneInfo]:
        """
        Search genes by symbol or description

        Args:
            query: Search string
            limit: Max results

        Returns:
            List of matching GeneInfo objects
        """
        query_upper = query.upper()
        results = []

        for gene in self.genes.values():
            # Check symbol
            if query_upper in gene.symbol.upper():
                results.append(gene)
                if len(results) >= limit:
                    break
                continue

            # Check description
            if query_upper in gene.description.upper():
                results.append(gene)
                if len(results) >= limit:
                    break
                continue

            # Check synonyms
            if any(query_upper in syn.upper() for syn in gene.synonyms):
                results.append(gene)
                if len(results) >= limit:
                    break

        return results


# Singleton instance
_global_gene_info_db: Optional[GeneInfoDB] = None


def get_gene_info_db(auto_load: bool = True) -> GeneInfoDB:
    """
    Get global gene info database instance

    Args:
        auto_load: Automatically load NCBI gene info if available

    Returns:
        GeneInfoDB instance
    """
    global _global_gene_info_db

    if _global_gene_info_db is None:
        _global_gene_info_db = GeneInfoDB()

        if auto_load:
            # Try to load NCBI gene info
            _global_gene_info_db.load_ncbi_gene_info()

    return _global_gene_info_db


def lookup_gene(query: str) -> Optional[GeneInfo]:
    """
    Quick gene lookup by symbol, ID, or Ensembl ID

    Args:
        query: Gene symbol, NCBI gene ID, or Ensembl ID

    Returns:
        GeneInfo if found, else None

    Examples:
        >>> info = lookup_gene("BRCA1")
        >>> print(info.summary())

        >>> info = lookup_gene("ENSG00000012048")
        >>> print(info.description)
    """
    db = get_gene_info_db()

    # Try symbol first
    result = db.get_by_symbol(query)
    if result:
        return result

    # Try NCBI gene ID
    result = db.get_by_id(query)
    if result:
        return result

    # Try Ensembl ID
    result = db.get_by_ensembl(query)
    return result
