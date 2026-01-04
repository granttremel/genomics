"""
Color management for the genome browser display.

This module provides feature color assignment using the unified Color system.
"""

from typing import Dict

from ggene.draw.colors import Colors
from ggene.draw.color import (
    Color, StyledColor, Layer,
    FeatureColorRegistry,
    RESET
)


# Module-level registry instance
_registry = FeatureColorRegistry()


class FColors(Colors):
    """Feature colors for terminal display.

    Uses the unified Color system via FeatureColorRegistry for deterministic
    feature-to-color mapping with small perceptual shifts for related varieties.

    Methods return Color objects that can be cast to ANSI strings via str().
    """

    # Legacy repeat colors for backward compatibility
    # (registry handles these now, but kept for direct access)
    repeat_colors = {
        "AluY": 141, "AluJ": 143, "AluS": 106, "Alu": 142,
        "L1MA": 53, "L1MB": 89, "L1MC": 126, "L1MD": 161, "L1ME": 162, "L1M": 125,
        "L1PA": 197, "L1PB": 199, "L1P": 198, "L1": 124,
        "L2a": 161, "L2b": 162, "L2c": 166, "L2d": 167, "L2": 160,
        "L3": 88, "L4": 52,
        "SVA_A": 95, "SVA_B": 97, "SVA_C": 60, "SVA_D": 132, "SVA": 96,
        "MIRb": 107, "MIRc": 109, "MIR3": 144, "MIR": 108,
        "default": 130
    }

    # Base category colors for reference
    SPLICE_COLOR = Color.from_hsl(200, 0.7, 0.5)
    PROMOTER_COLOR = Color.from_hsl(280, 0.6, 0.5)
    REPEAT_COLOR = Color.from_hsl(30, 0.7, 0.5)

    @classmethod
    def get_feature_type_color(cls, feature_type: str) -> Color:
        """Get color for a feature type.

        Args:
            feature_type: Feature type string (e.g., 'gene', 'exon', 'motif')

        Returns:
            Color instance
        """
        return _registry.get_color(feature_type.lower())

    @classmethod
    def get_feature_color(cls, feature) -> Color:
        """Get color for a feature object.

        Dispatches to the appropriate color method based on feature type.

        Args:
            feature: Feature object with feature_type, name, attributes

        Returns:
            Color instance
        """
        return _registry.get_feature_color(feature)

    @classmethod
    def get_styled_feature_color(cls, feature, bold: bool = False,
                                  layer: Layer = Layer.FEATURE) -> StyledColor:
        """Get styled color for a feature.

        Args:
            feature: Feature object
            bold: Whether to apply bold styling
            layer: Compositing layer

        Returns:
            StyledColor instance
        """
        color = cls.get_feature_color(feature)
        return StyledColor(fg=color, bold=bold, layer=layer)

    @classmethod
    def get_variant_color(cls, ref: str, alt: str) -> Color:
        """Get color based on variant type.

        Args:
            ref: Reference allele
            alt: Alternate allele

        Returns:
            Color instance
        """
        if len(ref) == len(alt):
            return Color.from_8bit(214)  # SNP - orange
        elif len(ref) > len(alt):
            return Color.from_8bit(202)  # Deletion - red-orange
        else:
            return Color.from_8bit(220)  # Insertion - yellow

    @classmethod
    def get_motif_color(cls, motif_feat) -> Color:
        """Get color for a motif feature.

        Uses the registry for deterministic color assignment.
        Related motifs (e.g., splice_donor, splice_acceptor) get similar
        but distinguishable colors.

        Args:
            motif_feat: Motif feature object

        Returns:
            Color instance
        """
        return _registry.get_feature_color(motif_feat)

    @classmethod
    def get_repeat_color(cls, rpt_feat) -> Color:
        """Get color for a repeat feature.

        Args:
            rpt_feat: Repeat feature object

        Returns:
            Color instance
        """
        return _registry.get_feature_color(rpt_feat)

    @classmethod
    def get_key(cls, cols: Dict[str, Color]) -> str:
        """Generate a legend/key string for named colors.

        Args:
            cols: Dict mapping names to Color instances

        Returns:
            Formatted key string with ANSI colors
        """
        key_parts = []
        for name, col in sorted(cols.items(), key=lambda k: k[0]):
            if isinstance(col, int):
                col = Color.from_8bit(col)
            key_parts.append(f"{col.to_ansi()}{name}")

        return "Key: " + ", ".join(key_parts) + RESET

    @classmethod
    def get_ansi(cls, feature) -> str:
        """Get ANSI escape string for a feature (backward compatibility).

        Args:
            feature: Feature object

        Returns:
            ANSI escape sequence string
        """
        return str(cls.get_feature_color(feature))


# Convenience functions at module level
def feature_color(feature) -> Color:
    """Get color for a feature."""
    return FColors.get_feature_color(feature)


def styled_feature_color(feature, bold: bool = False,
                         layer: Layer = Layer.FEATURE) -> StyledColor:
    """Get styled color for a feature."""
    return FColors.get_styled_feature_color(feature, bold=bold, layer=layer)
