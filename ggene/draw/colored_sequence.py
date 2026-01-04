"""
ColoredSequence: Object-oriented colored sequence for multi-layer highlighting.

Stores bases and colors separately, allowing multiple rounds of coloring
before final rendering to an ANSI string.
"""

from typing import List, Tuple, Dict, Optional, Union, Iterable
import regex

from .color import Color, StyledColor, Layer, RESET


# Common color presets
DEFAULT_COLORS = {
    'base': StyledColor.from_8bit(fg=240),
    'match': StyledColor.from_8bit(fg=142, bold=True),
    'match_rc': StyledColor.from_8bit(fg=168, bold=True),
    'highlight': StyledColor.from_8bit(fg=226, bold=True),
    'feature': StyledColor.from_8bit(fg=81, layer=Layer.FEATURE),
    'motif': StyledColor.from_8bit(fg=208, layer=Layer.MOTIF),
    'error': StyledColor.from_8bit(fg=196),
}


class ColoredSequence:
    """
    A sequence with per-position color information.

    Allows multiple rounds of coloring - each coloring operation can be
    applied as a layer, with higher layers taking precedence during rendering.

    Example:
        cs = ColoredSequence("ATCGATCG")
        cs.color_range(0, 4, StyledColor.from_8bit(fg=81))    # Color first 4 bases
        cs.color_matches("TCG", StyledColor.from_8bit(fg=208)) # Highlight TCG matches
        print(cs.render())                                      # Get ANSI string
    """

    def __init__(self, sequence: str, default_color: Optional[StyledColor] = None):
        """
        Initialize a ColoredSequence.

        Args:
            sequence: The base sequence string
            default_color: Default color for all positions (None = no color)
        """
        self.sequence = sequence
        self.length = len(sequence)

        # Per-position color info
        if default_color is None:
            default_color = DEFAULT_COLORS['base']
        self._colors: List[StyledColor] = [default_color for _ in range(self.length)]

        # Track named regions for legend/key generation
        self._named_regions: Dict[str, Tuple[StyledColor, List[Tuple[int, int]]]] = {}

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, idx: Union[int, slice]) -> 'ColoredSequence':
        """Get a slice of the colored sequence."""
        if isinstance(idx, int):
            new_seq = ColoredSequence(self.sequence[idx])
            new_seq._colors = [self._colors[idx]]
            return new_seq
        else:
            new_seq = ColoredSequence(self.sequence[idx])
            new_seq._colors = self._colors[idx]
            return new_seq

    def copy(self) -> 'ColoredSequence':
        """Create a deep copy."""
        new_seq = ColoredSequence(self.sequence)
        # StyledColor is immutable/frozen, so we can share references
        new_seq._colors = list(self._colors)
        new_seq._named_regions = {
            name: (spec, list(spans))
            for name, (spec, spans) in self._named_regions.items()
        }
        return new_seq

    # -------------------------------------------------------------------------
    # Core coloring methods
    # -------------------------------------------------------------------------

    def color_position(self, pos: int, color: StyledColor,
                       replace: bool = False) -> 'ColoredSequence':
        """
        Color a single position.

        Args:
            pos: Position to color
            color: Color specification
            replace: If True, replace existing color. If False, merge (color on top).
        """
        if 0 <= pos < self.length:
            if replace or self._colors[pos].layer <= color.layer:
                if replace:
                    self._colors[pos] = color
                else:
                    self._colors[pos] = color.merge_over(self._colors[pos])
        return self

    def color_range(self, start: int, end: int, color: StyledColor,
                    name: Optional[str] = None, replace: bool = False) -> 'ColoredSequence':
        """
        Color a range of positions.

        Args:
            start: Start position (inclusive)
            end: End position (exclusive)
            color: Color specification
            name: Optional name for this region (for legend)
            replace: If True, replace existing colors
        """
        start = max(0, start)
        end = min(self.length, end)

        for i in range(start, end):
            self.color_position(i, color, replace)

        if name:
            if name not in self._named_regions:
                self._named_regions[name] = (color, [])
            self._named_regions[name][1].append((start, end))

        return self

    def color_spans(self, spans: Iterable[Tuple[int, int]], color: StyledColor,
                   name: Optional[str] = None) -> 'ColoredSequence':
        """Color multiple spans with the same color."""
        for start, end in spans:
            self.color_range(start, end, color, name=name)
        return self

    def color_positions(self, positions: Iterable[int],
                        color: StyledColor) -> 'ColoredSequence':
        """Color multiple individual positions."""
        for pos in positions:
            self.color_position(pos, color)
        return self

    # -------------------------------------------------------------------------
    # Pattern-based coloring
    # -------------------------------------------------------------------------

    def color_matches(self, pattern: str, color: StyledColor,
                     name: Optional[str] = None) -> 'ColoredSequence':
        """
        Color all exact matches of a pattern.

        Args:
            pattern: Substring to find
            color: Color for matches
            name: Optional name for legend
        """
        start = 0
        pattern_upper = pattern.upper()
        seq_upper = self.sequence.upper()

        while True:
            pos = seq_upper.find(pattern_upper, start)
            if pos == -1:
                break
            self.color_range(pos, pos + len(pattern), color, name=name)
            start = pos + 1

        return self

    def color_regex(self, pattern: str, color: StyledColor,
                   name: Optional[str] = None) -> 'ColoredSequence':
        """
        Color all regex matches.

        Args:
            pattern: Regex pattern
            color: Color for matches
            name: Optional name for legend
        """
        for match in regex.finditer(pattern, self.sequence, regex.IGNORECASE):
            self.color_range(match.start(), match.end(), color, name=name)
        return self

    def color_fuzzy(self, pattern: str, color: StyledColor, max_errors: int = 1,
                   error_colors: Optional[List[StyledColor]] = None,
                   name: Optional[str] = None) -> 'ColoredSequence':
        """
        Color fuzzy matches (allowing substitutions/insertions/deletions).

        Args:
            pattern: Pattern to match
            color: Color for exact matches
            max_errors: Maximum edit distance
            error_colors: List of colors by error count (index 0 = exact, 1 = 1 error, etc.)
            name: Optional name for legend
        """
        if error_colors is None:
            # Default: color degrades with errors
            error_colors = [color]
            for i in range(1, max_errors + 1):
                # Shift color toward gray with more errors
                if color.fg:
                    dimmed = color.fg.adjust_lightness(-0.1 * i).adjust_saturation(-0.15 * i)
                    error_colors.append(StyledColor(fg=dimmed, bold=color.bold))
                else:
                    error_colors.append(color)

        pattern_str = f"({pattern}){{e<={max_errors}}}"
        compiled = regex.compile(pattern_str, regex.BESTMATCH | regex.IGNORECASE)

        for match in regex.finditer(compiled, self.sequence):
            errors = sum(match.fuzzy_counts)
            match_color = error_colors[min(errors, len(error_colors) - 1)]
            start, end = match.span()
            self.color_range(start, end, match_color, name=name)

        return self

    # -------------------------------------------------------------------------
    # Comparison coloring
    # -------------------------------------------------------------------------

    def color_matching_positions(self, other: str, match_color: StyledColor,
                                 mismatch_color: Optional[StyledColor] = None) -> 'ColoredSequence':
        """
        Color positions that match/mismatch with another sequence.

        Args:
            other: Sequence to compare against
            match_color: Color for matching positions
            mismatch_color: Color for mismatching positions (None = leave unchanged)
        """
        min_len = min(len(self), len(other))
        for i in range(min_len):
            if self.sequence[i].upper() == other[i].upper():
                self.color_position(i, match_color)
            elif mismatch_color:
                self.color_position(i, mismatch_color)
        return self

    def color_complement_matches(self, other: str, color: StyledColor) -> 'ColoredSequence':
        """Color positions where self matches complement of other (for RC detection)."""
        from ggene.seqs.bio import complement

        min_len = min(len(self), len(other))
        for i in range(min_len):
            if self.sequence[i].upper() == complement(other[i]).upper():
                self.color_position(i, color)
        return self

    # -------------------------------------------------------------------------
    # Feature-based coloring
    # -------------------------------------------------------------------------

    def color_features(self, features: Dict[str, List[Tuple[int, int]]],
                      colors: Optional[Dict[str, StyledColor]] = None,
                      auto_color: bool = True) -> 'ColoredSequence':
        """
        Color named features with their spans.

        Args:
            features: Dict mapping feature names to list of (start, end) spans
            colors: Dict mapping feature names to colors
            auto_color: If True, auto-assign colors to features without explicit colors
        """
        if colors is None:
            colors = {}

        for idx, (name, spans) in enumerate(features.items()):
            if name in colors:
                color = colors[name]
            elif auto_color:
                # Generate a consistent color based on index
                hue = (idx * 137.5) % 360  # Golden angle for good distribution
                base_color = Color.from_hsl(hue, 0.6, 0.5)
                color = StyledColor(fg=base_color, layer=Layer.FEATURE)
                colors[name] = color
            else:
                continue

            self.color_spans(spans, color, name=name)

        return self

    def color_feature(self, feature, color: Optional[StyledColor] = None) -> 'ColoredSequence':
        """
        Color a single feature object.

        Args:
            feature: Feature with start, end, name attributes
            color: Color to use (if None, auto-generates from feature)
        """
        from .color import get_styled_feature_color

        if color is None:
            color = get_styled_feature_color(feature, layer=Layer.FEATURE)

        start = getattr(feature, 'start', 0)
        end = getattr(feature, 'end', 0)
        name = getattr(feature, 'name', None)

        self.color_range(start, end, color, name=name)
        return self

    # -------------------------------------------------------------------------
    # Rendering
    # -------------------------------------------------------------------------

    def render(self, include_reset: bool = True) -> str:
        """
        Render to ANSI-colored string.

        Args:
            include_reset: If True, append reset code at end
        """
        parts = []
        last_ansi = ""

        for base, color in zip(self.sequence, self._colors):
            ansi = color.to_ansi()
            if ansi != last_ansi:
                if last_ansi:  # Reset before changing
                    parts.append(RESET)
                parts.append(ansi)
                last_ansi = ansi
            parts.append(base)

        if include_reset:
            parts.append(RESET)

        return "".join(parts)

    def render_plain(self) -> str:
        """Render without colors (just the sequence)."""
        return self.sequence

    def render_with_ruler(self, start_pos: int = 0, interval: int = 10) -> List[str]:
        """Render with a position ruler above."""
        ruler = []
        seq_line = self.render()

        # Build ruler
        for i in range(0, self.length, interval):
            pos_label = str(start_pos + i)
            ruler.append(pos_label + " " * (interval - len(pos_label)))

        return ["".join(ruler)[:self.length], seq_line]

    def get_legend(self) -> str:
        """Generate a legend/key for named regions."""
        if not self._named_regions:
            return ""

        parts = ["Key:"]
        for name, (color, spans) in self._named_regions.items():
            count = len(spans)
            ansi = color.to_ansi()
            parts.append(f" {ansi}{name}{RESET}({count})")

        return " ".join(parts)

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def reset_colors(self, default_color: Optional[StyledColor] = None) -> 'ColoredSequence':
        """Reset all colors to default."""
        if default_color is None:
            default_color = DEFAULT_COLORS['base']

        self._colors = [default_color for _ in range(self.length)]
        self._named_regions.clear()
        return self

    def get_color_at(self, pos: int) -> StyledColor:
        """Get color at a position."""
        return self._colors[pos]

    def has_color_at(self, pos: int) -> bool:
        """Check if position has non-default coloring."""
        return self._colors[pos].layer > Layer.BASE


# -----------------------------------------------------------------------------
# Factory functions for common highlighting patterns
# -----------------------------------------------------------------------------

def highlight_subsequence(seq: str, pattern: str,
                         color: Optional[StyledColor] = None) -> ColoredSequence:
    """Highlight all occurrences of a pattern in a sequence."""
    if color is None:
        # Generate color from pattern hash
        hue = hash(pattern) % 360
        base_color = Color.from_hsl(hue, 0.6, 0.5)
        color = StyledColor(fg=base_color, layer=Layer.MATCH)

    cs = ColoredSequence(seq)
    cs.color_matches(pattern, color, name=pattern)
    return cs


def highlight_comparison(seq_a: str, seq_b: str,
                        match_color: Optional[StyledColor] = None,
                        mismatch_color: Optional[StyledColor] = None) -> Tuple[ColoredSequence, ColoredSequence]:
    """
    Highlight matching/mismatching positions between two sequences.

    Returns two ColoredSequences (one for each input).
    """
    if match_color is None:
        match_color = DEFAULT_COLORS['match']
    if mismatch_color is None:
        mismatch_color = StyledColor.from_8bit(fg=240)  # Gray for mismatch

    cs_a = ColoredSequence(seq_a)
    cs_b = ColoredSequence(seq_b)

    cs_a.color_matching_positions(seq_b, match_color, mismatch_color)
    cs_b.color_matching_positions(seq_a, match_color, mismatch_color)

    return cs_a, cs_b


def highlight_features_in_sequence(seq: str,
                                  features: Dict[str, List[Tuple[int, int]]],
                                  colors: Optional[Dict[str, StyledColor]] = None) -> ColoredSequence:
    """Highlight features with named spans in a sequence."""
    cs = ColoredSequence(seq)
    cs.color_features(features, colors)
    return cs
