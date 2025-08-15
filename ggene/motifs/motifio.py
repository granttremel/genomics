
from motifs.motif import BaseMotif

class MotifIO:
    """Import/export motifs from standard formats"""

    @staticmethod
    def from_jaspar(jaspar_file):
        """Load TF binding motifs from JASPAR database"""
        pass

    @staticmethod
    def from_meme(meme_file):
        """Load motifs from MEME suite output"""
        pass

    @staticmethod
    def to_bed(motif_instances):
        """Export found motifs as BED for genome browser"""
        pass