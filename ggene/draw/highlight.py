"""
Sequence highlighting utilities for genome browser.

This module provides functions for highlighting sequences with various patterns,
features, and comparisons. All functions use the ColoredSequence class internally
for consistent, composable highlighting.

Common Interface:
    All highlight_* functions follow this pattern:
    - Take a sequence (str) or ColoredSequence as first argument
    - Return a ColoredSequence (can be rendered with .render())
    - Support optional `name` parameter for legend generation
    - Colors can be specified as int (256-color), ColorSpec, or omitted for auto

Example:
    from ggene.draw.highlight import highlight, ColoredSequence, ColorSpec

    # Simple highlighting
    cs = highlight(seq, "TATAA", color=208)
    print(cs.render())

    # Chained highlighting
    cs = (ColoredSequence(seq)
          .color_matches("ATG", ColorSpec(fg=46, bold=True), name="start")
          .color_regex(r"TAA|TAG|TGA", ColorSpec(fg=196), name="stop")
          .color_fuzzy("TATAAA", ColorSpec(fg=226), max_errors=1))
    print(cs.render())
    print(cs.get_legend())
"""

from typing import List, Tuple, Dict, Any, Optional, Union
import random
import regex

from ggene.seqs.bio import complement
from ggene.seqs.find import find_subsequence, find_subsequences
from .colors import Colors
from .color import Color, StyledColor, Layer
from .colored_sequence import (
    ColoredSequence, DEFAULT_COLORS,
    highlight_subsequence, highlight_comparison, highlight_features_in_sequence
)

# Alias for backward compatibility
ColorSpec = StyledColor

# Re-export for convenience
__all__ = [
    'ColoredSequence', 'StyledColor', 'ColorSpec', 'Layer', 'DEFAULT_COLORS',
    'highlight', 'highlight_pattern', 'highlight_patterns',
    'highlight_fuzzy', 'highlight_features', 'highlight_spans',
    'highlight_comparison', 'highlight_dyads',
    'make_legend',
]


# ─────────────────────────────────────────────────────────────────────────────
# Unified Interface
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_colored_sequence(seq: Union[str, ColoredSequence]) -> ColoredSequence:
    """Convert string to ColoredSequence if needed."""
    if isinstance(seq, ColoredSequence):
        return seq
    return ColoredSequence(seq)


def _ensure_color_spec(color: Union[int, StyledColor, None],
                      layer: int = Layer.MATCH) -> StyledColor:
    """Convert int or None to StyledColor."""
    if color is None:
        hue = random.randint(0, 360)
        return StyledColor(fg=Color.from_hsl(hue, 0.6, 0.5), layer=layer)
    if isinstance(color, int):
        return StyledColor.from_8bit(fg=color, layer=layer)
    return color


def highlight(seq: Union[str, ColoredSequence],
              pattern: str,
              color: Union[int, ColorSpec, None] = None,
              name: Optional[str] = None) -> ColoredSequence:
    """
    Highlight all occurrences of a pattern in a sequence.

    This is the primary highlighting function for simple use cases.

    Args:
        seq: Sequence string or ColoredSequence
        pattern: Pattern to find and highlight
        color: Color as int (256-color), ColorSpec, or None for random
        name: Name for legend (defaults to pattern if not specified)

    Returns:
        ColoredSequence with highlights applied
    """
    cs = _ensure_colored_sequence(seq)
    color_spec = _ensure_color_spec(color)

    cs.color_matches(pattern, color_spec, name=name or pattern)
    return cs


def highlight_pattern(seq: Union[str, ColoredSequence],
                     pattern: str,
                     color: Union[int, ColorSpec, None] = None,
                     regex_mode: bool = False,
                     name: Optional[str] = None) -> ColoredSequence:
    """
    Highlight a pattern (exact or regex).

    Args:
        seq: Sequence string or ColoredSequence
        pattern: Pattern to highlight
        color: Color specification
        regex_mode: If True, treat pattern as regex
        name: Name for legend
    """
    cs = _ensure_colored_sequence(seq)
    color_spec = _ensure_color_spec(color)

    if regex_mode:
        cs.color_regex(pattern, color_spec, name=name)
    else:
        cs.color_matches(pattern, color_spec, name=name or pattern)

    return cs


def highlight_patterns(seq: Union[str, ColoredSequence],
                      patterns: List[str],
                      colors: Optional[Dict[str, Union[int, ColorSpec]]] = None) -> ColoredSequence:
    """
    Highlight multiple patterns in a sequence.

    Args:
        seq: Sequence string or ColoredSequence
        patterns: List of patterns to highlight
        colors: Optional dict mapping pattern -> color
    """
    cs = _ensure_colored_sequence(seq)
    if colors is None:
        colors = {}

    for pattern in patterns:
        if pattern in colors:
            color_spec = _ensure_color_spec(colors[pattern])
        else:
            color_spec = StyledColor.from_8bit(fg=random.randint(20, 230), layer=Layer.MATCH)
            colors[pattern] = color_spec

        cs.color_matches(pattern, color_spec, name=pattern)

    return cs


def highlight_fuzzy(seq: Union[str, ColoredSequence],
                   pattern: str,
                   max_errors: int = 1,
                   color: Union[int, ColorSpec, None] = None,
                   error_colors: Optional[List[Union[int, ColorSpec]]] = None,
                   name: Optional[str] = None) -> ColoredSequence:
    """
    Highlight fuzzy matches (allowing substitutions/insertions/deletions).

    Args:
        seq: Sequence string or ColoredSequence
        pattern: Pattern to match
        max_errors: Maximum edit distance
        color: Color for exact matches
        error_colors: Colors by error count [0 errors, 1 error, ...]
        name: Name for legend
    """
    cs = _ensure_colored_sequence(seq)
    color_spec = _ensure_color_spec(color)

    if error_colors:
        error_color_specs = [_ensure_color_spec(c) for c in error_colors]
    else:
        error_color_specs = None

    cs.color_fuzzy(pattern, color_spec, max_errors, error_color_specs, name=name)
    return cs


def highlight_features(seq: Union[str, ColoredSequence],
                      features: Dict[str, List[Tuple[int, int]]],
                      colors: Optional[Dict[str, Union[int, ColorSpec]]] = None) -> ColoredSequence:
    """
    Highlight named features with their spans.

    Args:
        seq: Sequence string or ColoredSequence
        features: Dict mapping feature names to list of (start, end) spans
        colors: Optional dict mapping feature names to colors
    """
    cs = _ensure_colored_sequence(seq)

    color_specs = {}
    if colors:
        for name, color in colors.items():
            color_specs[name] = _ensure_color_spec(color, layer=Layer.FEATURE)
            # print(f"validated color {color} to {color_specs[name]}")

    cs.color_features(features, color_specs, auto_color=True)
    return cs


def highlight_spans(seq: Union[str, ColoredSequence],
                   span_colors: Dict[Tuple[int, int], Union[int, ColorSpec]],
                   default_color: Union[int, ColorSpec, None] = None) -> ColoredSequence:
    """
    Highlight by explicit span -> color mapping.

    Args:
        seq: Sequence string or ColoredSequence
        span_colors: Dict mapping (start, end) tuples to colors
        default_color: Color for non-highlighted positions
    """
    if default_color:
        default_spec = _ensure_color_spec(default_color, layer=Layer.BASE)
        cs = ColoredSequence(seq if isinstance(seq, str) else seq.sequence, default_spec)
    else:
        cs = _ensure_colored_sequence(seq)

    for (start, end), color in span_colors.items():
        color_spec = _ensure_color_spec(color)
        cs.color_range(start, end, color_spec)

    return cs


def highlight_matching_bases(seq_a: Union[str, ColoredSequence],
                            seq_b: str,
                            match_color: Union[int, ColorSpec, None] = None,
                            rc_color: Union[int, ColorSpec, None] = None,
                            include_rc: bool = False) -> ColoredSequence:
    """
    Highlight positions where seq_a matches seq_b.

    Args:
        seq_a: Primary sequence to highlight
        seq_b: Comparison sequence
        match_color: Color for matching positions
        rc_color: Color for reverse complement matches
        include_rc: If True, also check for RC matches
    """
    cs = _ensure_colored_sequence(seq_a)
    match_spec = _ensure_color_spec(match_color) if match_color else DEFAULT_COLORS['match']

    cs.color_matching_positions(seq_b, match_spec)

    if include_rc:
        rc_spec = _ensure_color_spec(rc_color) if rc_color else DEFAULT_COLORS['match_rc']
        cs.color_complement_matches(seq_b[::-1], rc_spec)

    return cs


def highlight_dyads(seq: Union[str, ColoredSequence],
                   dyads: List[Any]) -> ColoredSequence:
    """
    Highlight dyad symmetry structures.

    Args:
        seq: Sequence string or ColoredSequence
        dyads: List of dyad objects with stem_start and end_position attributes
    """
    cs = _ensure_colored_sequence(seq)

    for dyad in dyads:
        color = StyledColor.from_8bit(fg=random.randint(20, 230), layer=Layer.MOTIF)
        start = dyad.stem_start
        end = dyad.end_position
        cs.color_range(start, end + 1, color)

    return cs


def make_legend(cs: ColoredSequence) -> str:
    """Generate a legend/key for a colored sequence."""
    return cs.get_legend()


# ─────────────────────────────────────────────────────────────────────────────
# Backward-compatible helper functions
# ─────────────────────────────────────────────────────────────────────────────

def make_start_ends(feature_name: str,
                   feature_positions: List[int],
                   feature_length: int,
                   starts: Optional[Dict] = None,
                   ends: Optional[Dict] = None) -> Tuple[Dict, Dict]:
    """
    Build start/end position dicts from feature positions.

    This is a helper for converting position-based features to span-based.
    """
    if starts is None:
        starts = {}
    if ends is None:
        ends = {}

    for p in feature_positions:
        if p not in starts:
            starts[p] = []
        starts[p].append(feature_name)

        end_pos = p + feature_length
        if end_pos not in ends:
            ends[end_pos] = []
        ends[end_pos].append(feature_name)

    return starts, ends


def make_spans(starts: Dict, ends: Dict) -> Dict[str, List[Tuple[int, int]]]:
    """Convert starts/ends dicts to feature -> spans dict."""
    feats = set()
    istarts = {}
    iends = {}
    spans = {}

    for s, fs in starts.items():
        for f in fs:
            if f not in istarts:
                feats.add(f)
                istarts[f] = []
            istarts[f].append(s)

    for e, fs in ends.items():
        for f in fs:
            if f not in iends:
                feats.add(f)
                iends[f] = []
            iends[f].append(e)

    for f in feats:
        st = istarts.get(f, [])
        en = iends.get(f, [])
        spans[f] = [(s, e) for s, e in zip(sorted(st), sorted(en))]

    return spans


# ─────────────────────────────────────────────────────────────────────────────
# Legacy function wrappers (for backward compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def highlight_sequence(seq: str, subseq: str,
                      color: Optional[int] = None,
                      suppress: bool = True,
                      show_key: bool = True) -> Tuple[str, str]:
    """
    [LEGACY] Highlight a single subsequence. Use highlight() instead.
    """
    cs = highlight(seq, subseq, color, name=subseq)
    result = cs.render()
    key = cs.get_legend()

    if not suppress:
        print(result)
        if show_key:
            print(key)

    return result, key


def highlight_sequences(seq: str, subseqs: List[str],
                       colors: Optional[Dict[str, int]] = None,
                       suppress: bool = True,
                       show_key: bool = True,
                       **kwargs) -> Tuple[str, str]:
    """
    [LEGACY] Highlight multiple subsequences. Use highlight_patterns() instead.
    """
    cs = highlight_patterns(seq, subseqs, colors)
    result = cs.render()
    key = cs.get_legend()

    if not suppress:
        print(result)
        if show_key:
            print(key)

    return result, key


def highlight_sequence_fuzzy(seq: str, subseq: str,
                            colors: Optional[List[int]] = None,
                            max_err: int = 1,
                            suppress: bool = True,
                            show_key: bool = True,
                            **kwargs) -> Tuple[str, str]:
    """
    [LEGACY] Highlight fuzzy matches. Use highlight_fuzzy() instead.
    """
    cs = highlight_fuzzy(seq, subseq, max_err, error_colors=colors)
    result = cs.render()
    key = cs.get_legend()

    if not suppress:
        print(result)
        if show_key:
            print(key)

    return result, key


def highlight_sequence_by_span(seq: str,
                              span_colors: Dict[Tuple[int, int], Any],
                              default_color: str = '\033[97m') -> str:
    """
    [LEGACY] Highlight by span. Use highlight_spans() instead.
    """
    # Convert old color format to new
    converted_colors = {}
    for span, color in span_colors.items():
        if isinstance(color, str):
            # Try to extract color code from ANSI string
            import re
            match = re.search(r'\[38;5;(\d+)m', color)
            if match:
                converted_colors[span] = int(match.group(1))
            else:
                converted_colors[span] = 255  # White fallback
        else:
            converted_colors[span] = color

    cs = highlight_spans(seq, converted_colors)
    return cs.render()
