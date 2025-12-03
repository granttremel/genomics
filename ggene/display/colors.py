"""
Color management for the genome browser display.

This module provides ANSI color codes and color schemes for terminal-based
visualization of genomic data.
"""


class Colors:
    """ANSI color codes for terminal display."""

    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'

    # Variant colors
    SNP = '\033[91m'        # Red for SNPs
    INSERTION = '\033[92m'   # Green for insertions
    DELETION = '\033[93m'    # Yellow for deletions

    # Feature colors
    GENE = '\033[94m'        # Blue
    TRANSCRIPT = '\033[95m'  # Magenta
    EXON = '\033[96m'        # Cyan
    CDS = '\033[93m'         # Yellow
    UTR = '\033[90m'         # Gray
    
    START_CODON = '\x1b[35m'
    STOP_CODON = '\x1b[35m'

    # Motif colors (for underlines)
    MOTIF = '\033[38;5;110m'   # Cyan for other motifs
    RCMOTIF = '\x1b[38;5;106m'
    HIGHLIGHT = '\x1b[148m' # goldish
    
    # Navigation
    POSITION = '\033[97m'    # White

    @classmethod
    def variant_color(cls, ref: str, alt: str) -> str:
        """Get color based on variant type."""
        if len(ref) == len(alt):
            return cls.SNP
        elif len(ref) > len(alt):
            return cls.DELETION
        else:
            return cls.INSERTION

    @classmethod
    def get_feature_color(cls, feature_type: str) -> str:
        """Get color for a feature type."""
        feature_type_lower = feature_type.lower()

        if 'gene' in feature_type_lower:
            return cls.GENE
        elif 'transcript' in feature_type_lower:
            return cls.TRANSCRIPT
        elif 'exon' in feature_type_lower:
            return cls.EXON
        elif 'cds' in feature_type_lower or 'coding' in feature_type_lower:
            return cls.CDS
        elif 'utr' in feature_type_lower:
            return cls.UTR
        else:
            return cls.RESET

    @classmethod
    def get_motif_color(cls, motif_type: str) -> str:
        """Get color for a motif type."""
        motif_type_lower = motif_type.lower()

        if 'splice' in motif_type_lower:
            return cls.MOTIF
        elif 'tata' in motif_type_lower:
            return cls.MOTIF
        elif 'cpg' in motif_type_lower:
            return cls.MOTIF
        elif 'polya' in motif_type_lower:
            return cls.MOTIF
        else:
            return cls.MOTIF