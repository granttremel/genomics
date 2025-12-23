"""
Color management for the genome browser display.

This module provides ANSI color codes and color schemes for terminal-based
visualization of genomic data.
"""

from ggene.draw.colors import Colors



class FColors(Colors):
    """ANSI color codes for terminal display."""

    repeat_colors = {
        "AluY":141,
        "AluJ":143,
        "AluS":106,
        "Alu":142,
        
        "L1MA":53,
        "L1MB":89,
        "L1MC":126,
        "L1MD":161,
        "L1ME":162,
        "L1M":125,
        
        "L1PA":197,
        "L1PB":199,
        "L1P":198,
        
        "L1":124,
        
        "L2a":161,
        "L2b":162,
        "L2c":166,
        "L2d":167,
        "L2":160,
        
        "L3":88,
        "L4":52,
        
        "SVA_A":95,
        "SVA_B":97,
        "SVA_C":60,
        "SVA_D":132,
        "SVA":96,
        
        "MIRb":107,
        "MIRc":109,
        "MIR3":144,
        "MIR":108,
        
        
        
        "default":130
    }


    @classmethod
    def get_feature_type_color(cls, feature_type: str) -> str:
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
        elif feature_type_lower == "variant":
            return cls.INSERTION
        elif feature_type_lower == "motif":
            return cls.MOTIF
        elif feature_type_lower in ["repeat","dfam_hit"]:
            return cls(cls.repeat_colors.get("default", 0)).code
        else:
            return cls.RESET

    @classmethod
    def get_feature_color(cls, feature):
        
        if feature.feature_type in ["motif", "variant", "repeat", "dfam_hit"]:
            if feature.feature_type=="motif":
                return cls.get_motif_color(feature)
            elif feature.feature_type == "variant":
                return cls.get_variant_color(feature.attributes.get("ref",""), feature.attributes.get("alt",""))
            else:
                return cls.get_repeat_color(feature)
        elif "pseudo" in feature.attributes.get("gene_biotype",""):
            return cls.PSEUDO
        else:
            return cls.get_feature_type_color(feature.feature_type)
        
        return cls.RESET

    @classmethod
    def get_variant_color(cls, ref: str, alt: str) -> str:
        """Get color based on variant type."""
        if len(ref) == len(alt):
            return cls.SNP
        elif len(ref) > len(alt):
            return cls.DELETION
        else:
            return cls.INSERTION
    
    @classmethod
    def get_motif_color(cls, motif_feat) -> str:
        """Get color for a motif type."""
        # motif_type_lower = motif_class.lower()

        if isinstance(motif_feat, dict):
            motif_class = motif_feat.get("attributes",{}).get("","")
            motif_rc = motif_feat.get("attributes",{}).get("is_rc",False)
        else:
            motif_class = motif_feat.attributes.get("","")
            motif_rc = motif_feat.attributes.get("is_rc",False)

        if 'splice' in motif_class:
            return cls.RCMOTIF_SPL if motif_rc else cls.MOTIF_SPL
        elif 'promoter' in motif_class:
            return cls.RCMOTIF_PRO if motif_rc else cls.MOTIF_PRO
        else:
            return cls.RCMOTIF if motif_rc else cls.MOTIF
        
    @classmethod
    def get_repeat_color(cls, rpt_feat):
        
        out_col = None
        name = rpt_feat.name
        
        if not name:
            return cls.get_color(cls.repeat_colors.get("default", 0))
        
        for cn, col in cls.repeat_colors.items():
            if cn in name:
                out_col = col
                break
            
        if out_col is None:
            out_col = cls.repeat_colors.get("default", 0)
        
        return cls.get_color(out_col)
