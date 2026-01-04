"""Unified color system for the genomics browser.

This module provides a single Color class that stores colors in HSL format
for easy derivation, with output converters for 8-bit, RGB, and ANSI formats.

Classes:
    Color: Immutable color with HSL storage and format converters
    StyledColor: Color with terminal styling (bold, underline, layers)
    Layer: Priority enum for color compositing
    FeatureColorRegistry: Deterministic feature-to-color mapping

Usage:
    # Create colors from any format
    c = Color.from_8bit(142)
    c = Color.from_rgb(128, 200, 50)
    c = Color.from_hsl(120, 0.7, 0.5)

    # Output to any format
    print(c.to_8bit())      # -> 142
    print(c.to_rgb())       # -> (128, 200, 50)
    print(c.to_ansi())      # -> '\\x1b[38;5;142m'
    print(str(c))           # -> '\\x1b[38;5;142m'
    print(int(c))           # -> 142

    # Derive related colors
    lighter = c.adjust_lightness(0.1)
    similar = c.shift_hue(15)  # 15 degrees

    # Feature colors
    color = get_feature_color(feature)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple, Union
from enum import IntEnum
import colorsys
import hashlib


class Layer(IntEnum):
    """Priority layers for color compositing (higher = rendered on top)."""
    BASE = 0
    FEATURE = 10
    MATCH = 20
    MOTIF = 30
    HIGHLIGHT = 40
    CURSOR = 50


# Standard 16-color palette RGB values
_STANDARD_COLORS = [
    (0, 0, 0), (128, 0, 0), (0, 128, 0), (128, 128, 0),
    (0, 0, 128), (128, 0, 128), (0, 128, 128), (192, 192, 192),
    (128, 128, 128), (255, 0, 0), (0, 255, 0), (255, 255, 0),
    (0, 0, 255), (255, 0, 255), (0, 255, 255), (255, 255, 255)
]


def _8bit_to_rgb(code: int) -> Tuple[int, int, int]:
    """Convert 256-color code to RGB tuple."""
    colint = int(code)
    if colint < 0 or colint > 255:
        raise ValueError(f"8-bit color code must be 0-255, got {colint}")

    if colint < 16:
        # Standard colors
        return _STANDARD_COLORS[colint]
    elif colint < 232:
        # 216-color cube: 6x6x6
        colint -= 16
        r = (colint // 36) * 51
        g = ((colint // 6) % 6) * 51
        b = (colint % 6) * 51
        return (r, g, b)
    else:
        # Grayscale: 24 shades
        gray = (colint - 232) * 10 + 8
        return (gray, gray, gray)


def _rgb_to_8bit(rgb: Tuple[int, int, int]) -> int:
    """Find closest 256-color code for RGB tuple."""
    r, g, b = rgb

    # Clamp values
    r = max(0, min(255, r))
    g = max(0, min(255, g))
    b = max(0, min(255, b))

    # Check grayscale first (more accurate for grays)
    if r == g == b:
        if r < 8:
            return 16  # Black in color cube
        if r > 248:
            return 231  # White in color cube
        return round((r - 8) / 10) + 232

    # 216-color cube
    ri = round(r / 51)
    gi = round(g / 51)
    bi = round(b / 51)
    return 16 + (ri * 36) + (gi * 6) + bi


@dataclass(frozen=True)
class Color:
    """Immutable color with HSL internal storage.

    Internal representation uses HSL (Hue, Saturation, Lightness) which
    allows intuitive color derivation via small parameter shifts.

    Attributes:
        h: Hue in degrees (0-360)
        s: Saturation (0-1)
        l: Lightness (0-1)
    """
    h: float = 0.0
    s: float = 0.0
    l: float = 0.5

    # --- Constructors ---

    @classmethod
    def from_hsl(cls, h: float, s: float, l: float) -> Color:
        """Create from HSL values.

        Args:
            h: Hue in degrees (0-360, wraps around)
            s: Saturation (0-1, clamped)
            l: Lightness (0-1, clamped)
        """
        return cls(
            h=h % 360,
            s=max(0.0, min(1.0, s)),
            l=max(0.0, min(1.0, l))
        )

    @classmethod
    def from_rgb(cls, r: int, g: int, b: int) -> Color:
        """Create from RGB values (0-255 each)."""
        r_norm = max(0, min(255, r)) / 255
        g_norm = max(0, min(255, g)) / 255
        b_norm = max(0, min(255, b)) / 255

        # colorsys uses HLS order (not HSL)
        h, l, s = colorsys.rgb_to_hls(r_norm, g_norm, b_norm)
        return cls(h=h * 360, s=s, l=l)

    @classmethod
    def from_8bit(cls, code: int) -> Color:
        """Create from 256-color palette code."""
        rgb = _8bit_to_rgb(code)
        return cls.from_rgb(*rgb)

    @classmethod
    def from_ansi(cls, ansi_str: str) -> Color:
        """Parse ANSI escape sequence to Color.

        Handles formats:
            - 8-bit: '\\x1b[38;5;142m' or '\\x1b[48;5;142m'
            - 24-bit: '\\x1b[38;2;128;200;50m' or '\\x1b[48;2;r;g;bm'
        """
        # Strip escape and 'm'
        s = ansi_str.strip()
        if s.startswith('\x1b[') or s.startswith('\033['):
            s = s[2:]
        if s.endswith('m'):
            s = s[:-1]

        parts = s.split(';')

        # 8-bit format: 38;5;N or 48;5;N
        if len(parts) >= 3 and parts[1] == '5':
            code = int(parts[2])
            return cls.from_8bit(code)

        # 24-bit format: 38;2;R;G;B or 48;2;R;G;B
        if len(parts) >= 5 and parts[1] == '2':
            r, g, b = int(parts[2]), int(parts[3]), int(parts[4])
            return cls.from_rgb(r, g, b)

        # Fallback: try parsing as simple color code
        if len(parts) == 1 and parts[0].isdigit():
            return cls.from_8bit(int(parts[0]))

        # Default gray if unparseable
        return cls(h=0, s=0, l=0.5)

    # --- Output Converters ---

    def to_rgb(self) -> Tuple[int, int, int]:
        """Convert to RGB tuple (0-255 each)."""
        # colorsys uses HLS order
        r, g, b = colorsys.hls_to_rgb(self.h / 360, self.l, self.s)
        return (int(r * 255), int(g * 255), int(b * 255))

    def to_8bit(self) -> int:
        """Convert to closest 256-color palette code."""
        return _rgb_to_8bit(self.to_rgb())

    def to_ansi(self, bg: bool = False, mode: str = "8bit") -> str:
        """Convert to ANSI escape sequence.

        Args:
            bg: If True, emit background color (48), else foreground (38)
            mode: "8bit" for 256-color, "24bit" for true color

        Returns:
            ANSI escape sequence string
        """
        prefix = 48 if bg else 38
        if mode == "24bit":
            r, g, b = self.to_rgb()
            return f"\x1b[{prefix};2;{r};{g};{b}m"
        else:
            return f"\x1b[{prefix};5;{self.to_8bit()}m"

    def to_hex(self) -> str:
        """Convert to hex string (e.g., '#80c832')."""
        r, g, b = self.to_rgb()
        return f"#{r:02x}{g:02x}{b:02x}"

    # --- Color Derivation ---

    def shift_hue(self, degrees: float) -> Color:
        """Return new color with hue shifted by degrees."""
        return Color(h=(self.h + degrees) % 360, s=self.s, l=self.l)

    def adjust_saturation(self, delta: float) -> Color:
        """Return new color with saturation adjusted (clamped 0-1)."""
        return Color(h=self.h, s=max(0.0, min(1.0, self.s + delta)), l=self.l)

    def adjust_lightness(self, delta: float) -> Color:
        """Return new color with lightness adjusted (clamped 0-1)."""
        return Color(h=self.h, s=self.s, l=max(0.0, min(1.0, self.l + delta)))

    def derive(self, hue_shift: float = 0, sat_delta: float = 0,
               light_delta: float = 0) -> Color:
        """Return derived color with multiple adjustments."""
        return Color(
            h=(self.h + hue_shift) % 360,
            s=max(0.0, min(1.0, self.s + sat_delta)),
            l=max(0.0, min(1.0, self.l + light_delta))
        )

    def with_lightness(self, l: float) -> Color:
        """Return new color with specific lightness."""
        return Color(h=self.h, s=self.s, l=max(0.0, min(1.0, l)))

    def with_saturation(self, s: float) -> Color:
        """Return new color with specific saturation."""
        return Color(h=self.h, s=max(0.0, min(1.0, s)), l=self.l)

    # --- Dunder Methods ---

    def __str__(self) -> str:
        """String conversion emits foreground ANSI code."""
        return self.to_ansi()

    def __int__(self) -> int:
        """Int conversion emits 8-bit code."""
        return self.to_8bit()

    def __iter__(self):
        """Iteration yields RGB tuple components."""
        return iter(self.to_rgb())

    def __repr__(self) -> str:
        return f"Color(h={self.h:.1f}, s={self.s:.2f}, l={self.l:.2f})"


@dataclass(frozen=True)
class StyledColor:
    """Color with terminal styling attributes.

    Combines foreground/background colors with text styling (bold, underline)
    and compositing priority (layer).

    Attributes:
        fg: Foreground color (optional)
        bg: Background color (optional)
        bold: Bold text
        underline: Underlined text
        dim: Dim/faint text
        layer: Compositing priority
    """
    fg: Optional[Color] = None
    bg: Optional[Color] = None
    bold: bool = False
    underline: bool = False
    dim: bool = False
    layer: Layer = Layer.BASE

    @classmethod
    def from_fg(cls, color: Union[Color, int, Tuple[int, int, int]],
                **style) -> StyledColor:
        """Create from foreground color with optional styling."""
        if isinstance(color, int):
            color = Color.from_8bit(color)
        elif isinstance(color, tuple):
            color = Color.from_rgb(*color)
        return cls(fg=color, **style)

    @classmethod
    def from_8bit(cls, fg: Optional[int] = None, bg: Optional[int] = None,
                  bold: bool = False, underline: bool = False,
                  dim: bool = False, layer: Layer = Layer.BASE) -> StyledColor:
        """Create from 8-bit color codes."""
        return cls(
            fg=Color.from_8bit(fg) if fg is not None else None,
            bg=Color.from_8bit(bg) if bg is not None else None,
            bold=bold,
            underline=underline,
            dim=dim,
            layer=layer
        )

    def to_ansi(self, mode: str = "8bit") -> str:
        """Convert to ANSI escape sequence with all styling.

        Args:
            mode: "8bit" for 256-color, "24bit" for true color

        Returns:
            ANSI escape sequence string
        """
        codes = []

        if self.bold:
            codes.append("1")
        if self.dim:
            codes.append("2")
        if self.underline:
            codes.append("4")

        if self.fg is not None:
            if mode == "24bit":
                r, g, b = self.fg.to_rgb()
                codes.append(f"38;2;{r};{g};{b}")
            else:
                codes.append(f"38;5;{self.fg.to_8bit()}")

        if self.bg is not None:
            if mode == "24bit":
                r, g, b = self.bg.to_rgb()
                codes.append(f"48;2;{r};{g};{b}")
            else:
                codes.append(f"48;5;{self.bg.to_8bit()}")

        if not codes:
            return ""
        return f"\x1b[{';'.join(codes)}m"

    def merge_over(self, other: StyledColor) -> StyledColor:
        """Merge this styled color over another (this takes precedence).

        Used for layer compositing: higher layer colors "paint over" lower.
        """
        return StyledColor(
            fg=self.fg if self.fg is not None else other.fg,
            bg=self.bg if self.bg is not None else other.bg,
            bold=self.bold or other.bold,
            underline=self.underline or other.underline,
            dim=self.dim if not self.bold else other.dim,  # dim cancelled by bold
            layer=max(self.layer, other.layer)
        )

    def with_fg(self, color: Color) -> StyledColor:
        """Return copy with different foreground."""
        return StyledColor(
            fg=color, bg=self.bg, bold=self.bold,
            underline=self.underline, dim=self.dim, layer=self.layer
        )

    def with_bg(self, color: Color) -> StyledColor:
        """Return copy with different background."""
        return StyledColor(
            fg=self.fg, bg=color, bold=self.bold,
            underline=self.underline, dim=self.dim, layer=self.layer
        )

    def with_layer(self, layer: Layer) -> StyledColor:
        """Return copy with different layer."""
        return StyledColor(
            fg=self.fg, bg=self.bg, bold=self.bold,
            underline=self.underline, dim=self.dim, layer=layer
        )

    def with_bold(self, bold: bool = True) -> StyledColor:
        """Return copy with bold styling."""
        return StyledColor(
            fg=self.fg, bg=self.bg, bold=bold,
            underline=self.underline, dim=self.dim, layer=self.layer
        )

    def __str__(self) -> str:
        return self.to_ansi()

    def __repr__(self) -> str:
        parts = []
        if self.fg:
            parts.append(f"fg={self.fg!r}")
        if self.bg:
            parts.append(f"bg={self.bg!r}")
        if self.bold:
            parts.append("bold=True")
        if self.underline:
            parts.append("underline=True")
        if self.dim:
            parts.append("dim=True")
        if self.layer != Layer.BASE:
            parts.append(f"layer={self.layer.name}")
        return f"StyledColor({', '.join(parts)})"


class FeatureColorRegistry:
    """Deterministic feature-to-color mapping with hierarchical derivation.

    Design principles:
    1. Base colors defined for feature categories (splice, promoter, repeat, etc.)
    2. Known varieties get explicit small offsets from base
    3. Unknown varieties derive via stable hash with small shifts
    4. Shifts are small (max 12deg hue, 0.08 sat, 0.06 light) for family coherence

    Usage:
        registry = FeatureColorRegistry()
        color = registry.get_color("motif", variety="TATA", motif_class="promoter")
        color = registry.get_feature_color(some_feature_object)
    """

    # Base colors for feature categories (HSL)
    CATEGORY_BASES = {
        "splice": Color.from_hsl(200, 0.7, 0.5),      # Blue family
        "promoter": Color.from_hsl(280, 0.6, 0.5),    # Purple family
        "repeat": Color.from_hsl(30, 0.7, 0.5),       # Orange family
        "gene": Color.from_hsl(220, 0.6, 0.55),       # Blue
        "transcript": Color.from_hsl(300, 0.5, 0.5),  # Magenta
        "exon": Color.from_hsl(180, 0.6, 0.5),        # Cyan
        "cds": Color.from_hsl(50, 0.7, 0.55),         # Yellow
        "utr": Color.from_hsl(0, 0, 0.5),             # Gray
        "variant": Color.from_hsl(25, 0.8, 0.55),     # Orange
        "motif": Color.from_hsl(170, 0.6, 0.5),       # Teal
        "tf_binding": Color.from_hsl(140, 0.6, 0.4),  # Green
    }

    # Known varieties with explicit small offsets: (category, hue_shift, sat_delta, light_delta)
    KNOWN_VARIETIES = {
        # Splice family (base hue ~200, blue)
        "splice_donor": ("splice", 0, 0, 0),
        "donor": ("splice", 0, 0, 0),
        "splice_acceptor": ("splice", 15, 0, 0.05),
        "acceptor": ("splice", 15, 0, 0.05),
        "splice_branch": ("splice", -10, 0.1, -0.05),
        "branch": ("splice", -10, 0.1, -0.05),

        # Promoter family (base hue ~280, purple)
        "TATA_box": ("promoter", 0, 0, 0),
        "TATA": ("promoter", 0, 0, 0),
        "CAAT_box": ("promoter", 10, 0, 0),
        "CAAT": ("promoter", 10, 0, 0),
        "E_box": ("promoter", 20, 0, 0.05),
        "Y_box": ("promoter", 30, 0, 0),
        "GC_box": ("promoter", -10, 0.1, 0),
        "Kozak": ("promoter", -20, 0, 0.05),
        "PolyA": ("promoter", 40, 0, -0.05),
        "AU_rich": ("promoter", 45, 0.05, 0),
        "Cap_site": ("promoter", -15, 0, 0.03),
        "DRE": ("promoter", 50, 0, 0),
        "TCT": ("promoter", 55, 0.05, 0),
        "tRNA_euk_1": ("promoter", -25, 0, 0),
        "tRNA_euk_2": ("promoter", -30, 0, 0.03),

        # Repeat families (base hue ~30, orange)
        "Alu": ("repeat", 0, 0, 0),
        "AluY": ("repeat", 5, 0.05, 0),
        "AluJ": ("repeat", -5, 0, 0.05),
        "AluS": ("repeat", 10, 0, 0),
        "L1": ("repeat", 30, 0, 0),
        "L1PA": ("repeat", 35, 0.05, 0),
        "L1MA": ("repeat", 25, 0, 0.05),
        "L1ME": ("repeat", 40, 0, 0),
        "L2": ("repeat", 50, 0, 0),
        "L2a": ("repeat", 50, 0, 0.03),
        "L2b": ("repeat", 52, 0, 0),
        "L2c": ("repeat", 54, 0, -0.03),
        "SVA": ("repeat", -20, 0.1, 0),
        "SVA_A": ("repeat", -18, 0.1, 0),
        "SVA_B": ("repeat", -22, 0.1, 0.03),
        "SVA_C": ("repeat", -24, 0.08, 0),
        "SVA_D": ("repeat", -26, 0.1, -0.03),
        "MIR": ("repeat", -30, 0, 0),
        "MIRb": ("repeat", -32, 0, 0.03),
        "MIRc": ("repeat", -34, 0, 0),
        "MIR3": ("repeat", -28, 0.05, 0),
    }

    # Derivation bounds for unknown varieties
    MAX_HUE_SHIFT = 12.0
    MAX_SAT_DELTA = 0.08
    MAX_LIGHT_DELTA = 0.06

    def __init__(self):
        self._cache: dict = {}

    def get_color(self, feature_type: str, variety: Optional[str] = None,
                  motif_class: Optional[str] = None) -> Color:
        """Get deterministic color for a feature.

        Args:
            feature_type: e.g., "motif", "repeat", "gene"
            variety: Specific variety name, e.g., "donor", "AluY"
            motif_class: Category hint, e.g., "splice", "promoter"

        Returns:
            Color instance
        """
        cache_key = (feature_type, variety, motif_class)
        if cache_key in self._cache:
            return self._cache[cache_key]

        color = self._compute_color(feature_type, variety, motif_class)
        self._cache[cache_key] = color
        return color

    def _compute_color(self, feature_type: str, variety: Optional[str],
                       motif_class: Optional[str]) -> Color:
        """Compute color for a feature (cache miss path)."""

        # 1. Check if this is a known variety
        if variety and variety in self.KNOWN_VARIETIES:
            category, h_shift, s_delta, l_delta = self.KNOWN_VARIETIES[variety]
            base = self.CATEGORY_BASES.get(category, self.CATEGORY_BASES.get("motif"))
            return base.derive(hue_shift=h_shift, sat_delta=s_delta, light_delta=l_delta)

        # 2. Determine base category
        category = self._determine_category(feature_type, variety, motif_class)
        base = self.CATEGORY_BASES.get(category)

        if base is None:
            base = self.CATEGORY_BASES.get(feature_type)
        if base is None:
            base = Color.from_hsl(0, 0, 0.5)  # Fallback gray

        # 3. If we have a variety name, derive from base with stable hash
        if variety:
            return self._derive_from_hash(base, variety)

        return base

    def _determine_category(self, feature_type: str, variety: Optional[str],
                           motif_class: Optional[str]) -> str:
        """Determine the color category for a feature."""
        # Use motif_class hint if provided
        if motif_class:
            mc_lower = motif_class.lower()
            if "splice" in mc_lower:
                return "splice"
            if "promoter" in mc_lower:
                return "promoter"
            if "repeat" in mc_lower:
                return "repeat"

        # Check variety name for hints
        if variety:
            v_lower = variety.lower()
            if any(x in v_lower for x in ["donor", "acceptor", "branch", "splice"]):
                return "splice"
            if any(x in v_lower for x in ["tata", "caat", "kozak", "polya", "promoter", "cap"]):
                return "promoter"
            if any(x in v_lower for x in ["alu", "line", "sine", "l1", "l2", "sva", "mir"]):
                return "repeat"

        return feature_type

    def _derive_from_hash(self, base: Color, variety: str) -> Color:
        """Derive a color from base using stable hash of variety name."""
        h = hashlib.md5(variety.encode()).digest()

        # Extract 3 bytes for our 3 dimensions
        hue_byte = h[0]
        sat_byte = h[1]
        light_byte = h[2]

        # Map to small shifts (centered around 0)
        hue_shift = (hue_byte / 255 - 0.5) * 2 * self.MAX_HUE_SHIFT
        sat_delta = (sat_byte / 255 - 0.5) * 2 * self.MAX_SAT_DELTA
        light_delta = (light_byte / 255 - 0.5) * 2 * self.MAX_LIGHT_DELTA

        return base.derive(hue_shift=hue_shift, sat_delta=sat_delta, light_delta=light_delta)

    def get_feature_color(self, feature) -> Color:
        """Get color directly from a feature object.

        Extracts feature_type, name, and motif_class from feature attributes.
        """
        feature_type = getattr(feature, 'feature_type', 'unknown')
        variety = getattr(feature, 'name', None)
        motif_class = getattr(feature, 'motif_class', None)

        # Try to get motif_class from attributes dict
        if motif_class is None and hasattr(feature, 'attributes'):
            attrs = feature.attributes
            if isinstance(attrs, dict):
                motif_class = attrs.get('motif_class')

        return self.get_color(feature_type, variety, motif_class)

    def get_styled_feature_color(self, feature, bold: bool = False,
                                  layer: Layer = Layer.FEATURE) -> StyledColor:
        """Get styled color for a feature."""
        color = self.get_feature_color(feature)
        return StyledColor(fg=color, bold=bold, layer=layer)


# Module-level singleton registry
_registry = FeatureColorRegistry()


def get_feature_color(feature) -> Color:
    """Get color for a feature using the global registry."""
    return _registry.get_feature_color(feature)


def get_styled_feature_color(feature, bold: bool = False,
                             layer: Layer = Layer.FEATURE) -> StyledColor:
    """Get styled color for a feature using the global registry."""
    return _registry.get_styled_feature_color(feature, bold=bold, layer=layer)


# Common color constants for convenience
RESET = "\x1b[0m"
BOLD = "\x1b[1m"
DIM = "\x1b[2m"
UNDERLINE = "\x1b[4m"
