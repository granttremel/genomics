
from .colors import Colors
from .color import (
    Color, StyledColor, Layer,
    FeatureColorRegistry, get_feature_color, get_styled_feature_color,
    RESET, BOLD, DIM, UNDERLINE
)
from .format import format_genomic
from .chars import SCALE, SCALE_H, SCALAR_PLOT, HLINES, RULER, OTHER, get_arrow

from .scalar_plot import ScalarPlot, scalar_to_text_nb, scalar_to_text_mid, scalar_plot_distribution
from .ruler import Ruler, make_ruler, add_ruler

from .heatmap import make_heatmap, Heatmap

__all__ = [
    # Colors - new unified system
    "Color", "StyledColor", "Layer",
    "FeatureColorRegistry", "get_feature_color", "get_styled_feature_color",
    "RESET", "BOLD", "DIM", "UNDERLINE",

    # Colors - legacy
    "Colors",

    "format_genomic",

    "SCALE", "SCALE_H", "SCALAR_PLOT", "HLINES", "RULER", "OTHER", "get_arrow",

    "ScalarPlot", "scalar_to_text_nb", "scalar_to_text_mid", "scalar_plot_distribution",
    "Ruler", "make_ruler", "add_ruler",
    "make_heatmap", "Heatmap"
]
