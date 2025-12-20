"""
Ruler rendering class for scalar text plots.

This module provides a Ruler class that handles the complex logic of rendering
rulers with optimal label placement, including multi-row layouts to avoid label
collisions.
"""

from typing import List, Optional, Dict, Any, Callable, Union, Tuple
import numpy as np

from .format import format_genomic
from .chars import RULER, SCALE

class Ruler:
    """
    A class for creating and rendering rulers with intelligent label placement.

    This class handles:
    - Automatic or manual label/tick spacing
    - Multi-row label layout to avoid overlaps
    - Vertical drop lines connecting labels to their positions
    - Dynamic re-rendering when options change

    Example usage:
        ruler = Ruler(xmin=0, xmax=1000, num_cols=80)
        ruler.set_auto(True)
        ruler.render()
        ruler.show()

        # Add to scalar plot
        plot = ScalarPlot(data)
        plot.ruler = Ruler(xmin=0, xmax=100, num_cols=len(data))
        plot.render()
    """

    def __init__(self, xmin: float, xmax: float, num_cols: int, **kwargs):
        """
        Initialize a Ruler.

        Args:
            xmin: Minimum x value
            xmax: Maximum x value
            num_cols: Number of columns in the plot
            **kwargs: Options for rendering (see default_options)
        """
        self.xmin = xmin
        self.xmax = xmax
        self.num_cols = num_cols
        self.rows = []

        # Tick distances (populated after rendering)
        self.tick_distances = {
            'label_interval': 0,
            'tick_interval': 0,
            'minor_tick_interval': 0,
        }

        # Default options
        self.options = self._get_default_options()

        # Update with user-provided options
        self.options.update(kwargs)

        # Initial render
        self.render()

    def _get_default_options(self) -> Dict[str, Any]:
        """Get default options for ruler rendering."""
        return {
            'auto': False,
            'genomic': False,
            'num_labels': 5,
            'ticks': 5,
            'minor_ticks': 5,
            'formatter': None,  # None = auto-detect
            'multi_row_labels': False,
            'max_label_rows': 4,
            'drop_line_char': RULER.get("tick"),
            'drop_corner_char': RULER.get("rlabel_marker"),
        }

    def render(self) -> None:
        """
        Render the ruler based on current options.

        This builds all ruler rows including the main ruler line and any
        additional rows for dropped labels.
        """
        if self.options['auto']:
            params = self._auto_params()
        else:
            params = self._manual_params()

        num_labels, ticks_per_interval, minor_per_tick, formatter, use_minor_symbol = params

        # Build the main ruler structure
        label_positions, tick_positions, minor_positions = self._calculate_positions(
            num_labels, ticks_per_interval, minor_per_tick
        )

        # Calculate tick distances
        self._calculate_tick_distances(num_labels, ticks_per_interval, minor_per_tick)

        # Determine tick symbols
        tick_marker = RULER.get("minor_tick") if use_minor_symbol else RULER.get("tick")
        minor_marker = RULER.get("minor_tick")

        if self.options['multi_row_labels']:
            self.rows = self._render_multi_row(
                label_positions, tick_positions, minor_positions,
                formatter, tick_marker, minor_marker
            )
        else:
            self.rows = self._render_single_row(
                label_positions, tick_positions, minor_positions,
                formatter, tick_marker, minor_marker
            )

    def _calculate_positions(
        self,
        num_labels: int,
        ticks_per_interval: int,
        minor_per_tick: int
    ) -> Tuple[List[Tuple[int, float]], set, set]:
        """
        Calculate positions for labels, ticks, and minor ticks.

        Returns:
            (label_positions, tick_positions, minor_positions)
            where label_positions is [(col, x_value), ...]
        """
        xran = self.xmax - self.xmin

        # Calculate label positions
        label_positions = []
        for i in range(num_labels):
            x_val = self.xmin + (xran * i / (num_labels - 1))
            col_pos = int(round(self.num_cols * i / (num_labels - 1)))
            label_positions.append((col_pos, x_val))

        # Calculate tick positions (between labels)
        tick_positions = set()
        if ticks_per_interval > 0:
            for label_idx in range(num_labels - 1):
                start_col = label_positions[label_idx][0]
                end_col = label_positions[label_idx + 1][0]
                interval_cols = end_col - start_col

                # Ensure at least 2 columns per tick
                actual_ticks = min(ticks_per_interval, interval_cols // 2)

                for t in range(1, actual_ticks + 1):
                    tick_col = start_col + int(round(interval_cols * t / (actual_ticks)))
                    if tick_col not in [lp[0] for lp in label_positions]:
                        tick_positions.add(tick_col)

        # Calculate minor tick positions
        minor_positions = set()
        if minor_per_tick > 0:
            major_positions = sorted([lp[0] for lp in label_positions] + list(tick_positions))

            for idx in range(len(major_positions) - 1):
                start_col = major_positions[idx]
                end_col = major_positions[idx + 1]
                interval_cols = end_col - start_col

                # Ensure at least 2 columns per minor tick
                actual_minor = min(minor_per_tick, interval_cols // 2)

                for m in range(1, actual_minor + 1):
                    minor_col = start_col + int(round(interval_cols * m / (actual_minor)))
                    if minor_col not in major_positions:
                        minor_positions.add(minor_col)

        return label_positions, tick_positions, minor_positions

    def _calculate_tick_distances(
        self,
        num_labels: int,
        ticks_per_interval: int,
        minor_per_tick: int
    ) -> None:
        """
        Calculate the actual x-axis distances represented by each tick level.

        Note: ticks_per_interval and minor_per_tick represent the number of
        subdivisions, so the actual interval is divided by (count + 1).
        """
        xran = self.xmax - self.xmin

        # Distance between labels
        label_interval = xran / (num_labels - 1)
        self.tick_distances['label_interval'] = label_interval

        # Distance represented by each tick (ticks subdivide label intervals)
        if ticks_per_interval > 0:
            tick_interval = label_interval / (ticks_per_interval + 1)
            self.tick_distances['tick_interval'] = tick_interval
        else:
            self.tick_distances['tick_interval'] = 0

        # Distance represented by each minor tick (minor ticks subdivide tick intervals)
        if minor_per_tick > 0 and ticks_per_interval > 0:
            minor_tick_interval = tick_interval / (minor_per_tick + 1)
            self.tick_distances['minor_tick_interval'] = minor_tick_interval
        else:
            self.tick_distances['minor_tick_interval'] = 0

    def _render_single_row(
        self,
        label_positions: List[Tuple[int, float]],
        tick_positions: set,
        minor_positions: set,
        formatter: Callable,
        tick_marker: str,
        minor_marker: str
    ) -> List[str]:
        """Render ruler with all labels in a single row."""
        # Format final label to position it correctly
        
        llbl, rlbl = RULER.get("llabel_marker"), RULER.get("rlabel_marker")
        
        final_col, final_x = label_positions[-1]
        final_label_str = formatter(final_x) + rlbl
        final_label_width = len(final_label_str)
        final_label_start = self.num_cols - final_label_width


        ruler = []
        col = 0

        while col <= self.num_cols:
            # Check if at final label
            if col == final_label_start:
                ruler.append(final_label_str)
                break

            # Check for non-final label
            label_here = None
            for lp in label_positions[:-1]:
                if lp[0] == col:
                    label_here = lp
                    break

            if label_here:
                _, x_val = label_here
                label_str = llbl + formatter(x_val)
                ruler.append(label_str)
                col += len(label_str)
            elif col in tick_positions:
                ruler.append(tick_marker)
                col += 1
            elif col in minor_positions:
                ruler.append(minor_marker)
                col += 1
            else:
                ruler.append(" ")
                col += 1

        return ["".join(ruler)]

    def _render_multi_row(
        self,
        label_positions: List[Tuple[int, float]],
        tick_positions: set,
        minor_positions: set,
        formatter: Callable,
        tick_marker: str,
        minor_marker: str
    ) -> List[str]:
        """
        Render ruler with labels in multiple rows to avoid collisions.

        This uses a greedy algorithm to assign labels to rows, minimizing
        the total number of rows needed.
        """
        # Format all labels first
        formatted_labels = []
        for col, x_val in label_positions:
            label_str = formatter(x_val)
            formatted_labels.append((col, label_str))

        # Assign labels to rows
        label_rows = self._assign_label_rows(formatted_labels)

        # Build the rows
        rows = []
        max_row = max(label_rows.values())

        # Build from bottom to top (lowest row index = closest to ruler)
        for row_idx in range(max_row, -1, -1):
            row_labels = {col: label for (col, label), r in label_rows.items() if r == row_idx}

            if row_idx == 0:
                # Main ruler row with ticks and row-0 labels
                row = self._build_main_row(
                    row_labels, tick_positions, minor_positions,
                    tick_marker, minor_marker
                )
            else:
                # Drop row with vertical lines and labels
                labels_in_higher_rows = {
                    col: label for (col, label), r in label_rows.items() if r < row_idx
                }
                row = self._build_drop_row(row_labels, labels_in_higher_rows, row_idx)

            rows.append(row)

        return rows

    def _assign_label_rows(self, formatted_labels: List[Tuple[int, str]]) -> Dict[Tuple[int, str], int]:
        """
        Assign labels to rows to avoid overlaps, using minimum rows.

        Returns: {(col, label): row_index}
        """
        max_rows = self.options['max_label_rows']
        label_rows = {}

        # Sort labels by column position
        sorted_labels = sorted(formatted_labels, key=lambda x: x[0])

        for col, label in sorted_labels:
            label_width = len(label) + 1  # +1 for marker

            # Try to place in row 0 first, then row 1, etc.
            placed = False
            for row in range(max_rows):
                # Check if this label overlaps with any label already in this row
                overlaps = False
                for (other_col, other_label), other_row in label_rows.items():
                    if other_row != row:
                        continue

                    other_width = len(other_label) + 1

                    # Check for overlap
                    # This label spans [col, col + label_width)
                    # Other spans [other_col, other_col + other_width)
                    if not (col + label_width <= other_col or other_col + other_width <= col):
                        overlaps = True
                        break

                if not overlaps:
                    label_rows[(col, label)] = row
                    placed = True
                    break

            if not placed:
                # Couldn't place in any row within max_rows, use last row
                label_rows[(col, label)] = max_rows - 1

        return label_rows

    def _build_main_row(
        self,
        row_labels: Dict[int, str],
        tick_positions: set,
        minor_positions: set,
        tick_marker: str,
        minor_marker: str
    ) -> str:
        """Build the main ruler row with ticks and row-0 labels."""
        
        llbl, rlbl = RULER.get("llabel_marker"), RULER.get("rlabel_marker")
        ruler = []
        col = 0

        while col <= self.num_cols:
            # Check for label at this position
            if col in row_labels:
                label_str = llbl + row_labels[col]
                ruler.append(label_str)
                col += len(label_str)
            elif col in tick_positions:
                ruler.append(tick_marker)
                col += 1
            elif col in minor_positions:
                ruler.append(minor_marker)
                col += 1
            else:
                ruler.append(" ")
                col += 1

        return "".join(ruler)

    def _build_drop_row(
        self,
        row_labels: Dict[int, str],
        labels_in_higher_rows: Dict[int, str],
        row_idx: int
    ) -> str:
        """Build a drop row with vertical lines and labels for this row level."""
        drop_char = self.options['drop_line_char']
        corner_char = self.options['drop_corner_char']

        row = []
        col = 0

        while col <= self.num_cols:
            # Check if there's a label at this position in this row
            if col in row_labels:
                label_str = corner_char + row_labels[col]
                row.append(label_str)
                col += len(label_str)
            # Check if there's a label in a higher row that needs a drop line
            elif col in labels_in_higher_rows:
                row.append(drop_char)
                col += 1
            else:
                row.append(" ")
                col += 1

        return "".join(row)

    def _manual_params(self) -> Tuple[int, int, int, Callable, bool]:
        """Get ruler parameters from manual options."""

        num_labels = self.options['num_labels']
        ticks = self.options['ticks']
        minor_ticks = self.options['minor_ticks']
        use_minor_symbol = False

        xran = self.xmax - self.xmin
        lbl_delta = xran / (num_labels - 1)

        # Format selection
        if self.options['formatter']:
            fmtr = self.options['formatter']
            if isinstance(fmtr, str):
                frmstr = fmtr
                fmtr = lambda s: format(s, frmstr)
        elif self.options['genomic']:
            fmtr = self._get_genomic_formatter(xran)
        elif abs(xran) > 1e5:
            fmtr = lambda s: format(s, ".2e")
        elif abs(xran) < 1e-5:
            fmtr = lambda s: format(s, ".2e")
        elif abs(lbl_delta) >= 1 and all(abs(v - round(v)) < 1e-9 for v in [self.xmin, self.xmax, lbl_delta]):
            fmtr = lambda s: format(s, ".0f")
        elif abs(lbl_delta) >= 1:
            fmtr = lambda s: format(s, ".1f")
        else:
            decimals = max(1, int(-np.floor(np.log10(abs(lbl_delta)))) + 1)
            fmtr = lambda s: format(s, f".{decimals}f")

        return num_labels, ticks, minor_ticks, fmtr, use_minor_symbol

    def _auto_params(self) -> Tuple[int, int, int, Callable, bool]:
        """Get ruler parameters using automatic heuristics."""
        xran = self.xmax - self.xmin

        # Determine formatter
        if self.options['genomic']:
            fmtr = self._get_genomic_formatter(xran)
        else:
            fmtr = self._get_auto_formatter(xran)

        # Estimate label width
        test_val = self.xmax if abs(self.xmax) > abs(self.xmin) else self.xmin
        if isinstance(fmtr, str):
            test_str = format(test_val, fmtr)
            frmstr = fmtr
            fmtr = lambda s: format(s, frmstr)
        else:
            test_str = fmtr(test_val)
        avg_label_width = len(test_str) + 1

        # Determine number of labels
        max_labels_by_spacing = self.num_cols // 25
        max_labels_by_width = self.num_cols // (avg_label_width + 1)
        max_labels = max(3, min(max_labels_by_spacing, max_labels_by_width))

        nice_label_counts = [3, 4, 5, 6, 8, 10, 12, 15, 20]
        num_labels = 3
        for n in nice_label_counts:
            if n <= max_labels:
                num_labels = n
            else:
                break

        # Calculate ticks per label interval
        cols_per_label_interval = self.num_cols / (num_labels - 1)
        max_ticks_per_interval = int(cols_per_label_interval / 5)
        ticks_per_interval = max(1, max_ticks_per_interval)

        # Check if ticks are too close
        cols_per_tick = cols_per_label_interval / (ticks_per_interval + 1)
        use_minor_symbol = cols_per_tick <= 3

        # Minor ticks
        if use_minor_symbol:
            minor_per_tick = 0
        else:
            max_minor_per_tick = min(5, int(cols_per_tick) - 1)
            minor_per_tick = max(0, max_minor_per_tick)

        return num_labels, ticks_per_interval, minor_per_tick, fmtr, use_minor_symbol

    def _get_genomic_formatter(self, xran: float) -> Callable:
        """Get formatter for genomic coordinates."""
        nexp = np.log10(max(abs(self.xmax), abs(self.xmin), 1))
        if nexp > 6:
            div = 1e6
            unit = "M"
            decimals = 1 if xran/1e6 < 10 else 0
        elif nexp > 3:
            div = 1e3
            unit = "k"
            decimals = 1 if xran/1e3 < 10 else 0
        else:
            div = 1
            unit = "b"
            decimals = 0
        return lambda x: format(x/div, f".{decimals}f") + unit

    def _get_auto_formatter(self, xran: float) -> Callable:
        """Get formatter automatically based on range."""
        if abs(xran) > 1e5:
            return lambda s: format(s, ".2e")
        elif abs(xran) < 1e-3:
            return lambda s: format(s, ".2e")
        else:
            # Determine decimal places based on range
            if xran >= 100:
                decimals = 0
            elif xran >= 10:
                decimals = 1
            elif xran >= 1:
                decimals = 2
            else:
                decimals = max(0, int(-np.floor(np.log10(xran))) + 2)

            # Check if values are effectively integers
            # from ggene.draw import _get_nice_number
            nice_step = _get_nice_number(xran / 5)
            if nice_step >= 1 and all(abs(v - round(v)) < 1e-9 for v in [self.xmin, self.xmax, nice_step]):
                return lambda s: format(s, ".0f")
            else:
                return lambda s: format(s, f".{decimals}f")

    def show(self) -> None:
        """Print the rendered ruler to stdout."""
        for row in self.rows:
            print(row)

    def get_rows(self) -> List[str]:
        """Get the rendered rows without printing."""
        return self.rows.copy()

    def get_tick_distances(self) -> Dict[str, float]:
        """
        Get the x-axis distances represented by each tick level.

        Returns:
            Dict with keys: 'label_interval', 'tick_interval', 'minor_tick_interval'
        """
        return self.tick_distances.copy()

    def get_tick_legend(self, format_str: str = "g") -> str:
        """
        Generate a formatted string describing the tick distances.

        This is useful for understanding what each tick mark represents in
        terms of actual x-axis values.

        Args:
            format_str: Format specifier for the values (default: 'g' for general)

        Returns:
            A formatted string describing the tick distances

        Example:
            "Label interval: 500, Tick (╵): 100, Minor tick ('): 25"
        """
        parts = []

        label_int = self.tick_distances['label_interval']
        tick_int = self.tick_distances['tick_interval']
        minor_int = self.tick_distances['minor_tick_interval']

        if label_int > 0:
            parts.append(f"Label interval: {label_int:{format_str}}")

        if tick_int > 0:
            parts.append(f"Tick (╵): {tick_int:{format_str}}")

        if minor_int > 0:
            parts.append(f"Minor tick ('): {minor_int:{format_str}}")

        return ", ".join(parts) if parts else "No tick information available"

    def show_tick_legend(self, format_str: str = "g") -> None:
        """
        Print the tick distance legend to stdout.

        Args:
            format_str: Format specifier for the values (default: 'g' for general)
        """
        print(self.get_tick_legend(format_str))

    def calculate_position(self, num_labels: int = 0, num_ticks: int = 0, num_minor_ticks: int = 0) -> float:
        """
        Calculate the x-axis value at a given position by counting tick marks.

        This allows you to locate features by counting marks from the ruler origin.

        Args:
            num_labels: Number of label intervals to count
            num_ticks: Number of tick marks to count
            num_minor_ticks: Number of minor tick marks to count

        Returns:
            The x-axis value at the specified position

        Example:
            # Position at 2 labels + 3 ticks + 1 minor tick from the start
            x_val = ruler.calculate_position(num_labels=2, num_ticks=3, num_minor_ticks=1)
        """
        x_val = self.xmin
        x_val += num_labels * self.tick_distances['label_interval']
        x_val += num_ticks * self.tick_distances['tick_interval']
        x_val += num_minor_ticks * self.tick_distances['minor_tick_interval']
        return x_val

    def set_auto(self, auto: bool) -> 'Ruler':
        """
        Enable or disable auto mode.

        Args:
            auto: Whether to use automatic parameter selection

        Returns:
            self for method chaining
        """
        self.options['auto'] = auto
        self.render()
        return self

    def set_multi_row_labels(self, enable: bool, max_rows: int = 4) -> 'Ruler':
        """
        Enable multi-row label layout.

        Args:
            enable: Whether to use multi-row layout
            max_rows: Maximum number of label rows

        Returns:
            self for method chaining
        """
        self.options['multi_row_labels'] = enable
        self.options['max_label_rows'] = max_rows
        self.render()
        return self

    def set_range(self, xmin: float, xmax: float) -> 'Ruler':
        """
        Update the x-axis range.

        Args:
            xmin: New minimum value
            xmax: New maximum value

        Returns:
            self for method chaining
        """
        self.xmin = xmin
        self.xmax = xmax
        self.render()
        return self

    def set_option(self, **kwargs) -> 'Ruler':
        """
        Generic option setter.

        Args:
            **kwargs: Options to update

        Returns:
            self for method chaining
        """
        self.options.update(kwargs)
        self.render()
        return self

    def __repr__(self) -> str:
        mode = "auto" if self.options['auto'] else "manual"
        multi = "multi-row" if self.options['multi_row_labels'] else "single-row"
        return f"Ruler(xmin={self.xmin}, xmax={self.xmax}, cols={self.num_cols}, {mode}, {multi})"

    def __str__(self) -> str:
        return '\n'.join(self.rows)


def _get_nice_number(x, round_down=False):
    """Get a 'nice' number close to x (1, 2, 5, 10, 20, 50, 100, etc.)"""
    if x == 0:
        return 0

    exp = np.floor(np.log10(abs(x)))
    frac = abs(x) / (10 ** exp)

    if round_down:
        if frac < 1.5:
            nice_frac = 1
        elif frac < 3:
            nice_frac = 2
        elif frac < 7:
            nice_frac = 5
        else:
            nice_frac = 10
    else:
        if frac <= 1:
            nice_frac = 1
        elif frac <= 2:
            nice_frac = 2
        elif frac <= 5:
            nice_frac = 5
        else:
            nice_frac = 10

    return nice_frac * (10 ** exp) * (1 if x >= 0 else -1)


def _estimate_label_width(value, formatter):
    """Estimate the width of a formatted label."""
    if isinstance(formatter, str):
        test_str = format(value, formatter)
    else:
        test_str = formatter(value)
    # Add 1 for the tick marker
    return len(test_str) + 1


def _auto_ruler_params(xmin, xmax, num_cols, genomic=False):
    """
    Automatically determine optimal ruler parameters.

    Heuristics:
    - One label per 25 columns max, minimum 3 labels per plot
    - One tick per 5 columns max, minimum one tick per label
    - One minor tick per column max, no more than 5 minor ticks before tick/label
    - If ticks every 3 columns or fewer, use minor tick symbol and omit minor ticks

    Returns: (num_labels, ticks_per_label, minor_ticks_per_tick, formatter)
    """
    xran = xmax - xmin

    # Determine formatter first
    if genomic:
        nexp = np.log10(max(abs(xmax), abs(xmin), 1))
        if nexp > 6:
            div = 1e6
            unit = "M"
            decimals = 1 if xran/1e6 < 10 else 0
        elif nexp > 3:
            div = 1e3
            unit = "k"
            decimals = 1 if xran/1e3 < 10 else 0
        else:
            div = 1
            unit = "b"
            decimals = 0
        fmtr = lambda x: format(x/div, f".{decimals}f") + unit
    else:
        # Smart format selection
        if abs(xran) > 1e5:
            fmtr = ".2e"
        elif abs(xran) < 1e-3:
            fmtr = ".2e"
        else:
            # Determine decimal places based on range
            if xran >= 100:
                decimals = 0
            elif xran >= 10:
                decimals = 1
            elif xran >= 1:
                decimals = 2
            else:
                decimals = max(0, int(-np.floor(np.log10(xran))) + 2)

            # Check if values are effectively integers
            nice_step = _get_nice_number(xran / 5)
            if nice_step >= 1 and all(abs(v - round(v)) < 1e-9 for v in [xmin, xmax, nice_step]):
                fmtr = ".0f"
            else:
                fmtr = f".{decimals}f"

    # Estimate label width
    test_val = xmax if abs(xmax) > abs(xmin) else xmin
    avg_label_width = _estimate_label_width(test_val, fmtr)

    # Determine number of labels: one per 25 columns max, minimum 3
    max_labels_by_spacing = num_cols // 25
    max_labels_by_width = num_cols // (avg_label_width + 1)
    max_labels = max(3, min(max_labels_by_spacing, max_labels_by_width))

    # Choose a nice number of labels (prefer 3, 4, 5, 6, 8, 10, 12)
    nice_label_counts = [3, 4, 5, 6, 8, 10, 12, 15, 20]
    num_labels = 3  # default minimum
    for n in nice_label_counts:
        if n <= max_labels:
            num_labels = n
        else:
            break

    # Calculate ticks per label interval
    cols_per_label_interval = num_cols / (num_labels - 1)

    # One tick per 5 columns max, minimum one tick per label
    max_ticks_per_interval = int(cols_per_label_interval / 5)
    ticks_per_interval = max(1, max_ticks_per_interval)

    # Check if ticks would be too close (every 3 columns or fewer)
    cols_per_tick = cols_per_label_interval / (ticks_per_interval + 1)
    use_minor_symbol_for_ticks = cols_per_tick <= 3

    # Minor ticks: one per column max, no more than 5 before tick/label
    if use_minor_symbol_for_ticks:
        # Omit minor ticks when ticks are too close
        minor_per_tick = 0
    else:
        max_minor_per_tick = min(5, int(cols_per_tick) - 1)
        minor_per_tick = max(0, max_minor_per_tick)

    return num_labels, ticks_per_interval, minor_per_tick, fmtr, use_minor_symbol_for_ticks


def add_ruler(sctxt, xmin, xmax, genomic = False, auto=False, **kwargs):
    """
    Add a ruler to scalar text output.

    Args:
        sctxt: Scalar text output (list of strings)
        xmin: Minimum x value
        xmax: Maximum x value
        genomic: Use genomic formatting (kb, Mb)
        auto: Automatically determine optimal ruler parameters
        **kwargs: Manual parameters (num_labels, ticks, minor_ticks, fstr)

    Returns:
        (sctxt, dists): Updated scalar text and distance tuple
    """
    num_cols = sum(1 for s in sctxt[0] if s in SCALE)
    ran = (xmax - xmin)

    if auto:
        # Use automatic parameter selection
        num_labels, ticks_per_section, minor_per_tick, fmtr, use_minor_symbol = _auto_ruler_params(
            xmin, xmax, num_cols, genomic
        )
        ticks = ticks_per_section
        minor_ticks = minor_per_tick
    else:
        # Use manual parameters
        num_labels = kwargs.get("num_labels", 5)
        ticks = kwargs.get("ticks", 5)
        minor_ticks = kwargs.get("minor_ticks", 5)
        use_minor_symbol = False

        lbl_delta = ran / (num_labels - 1)

        # Format selection
        if kwargs.get("fstr"):
            fmtr = kwargs.get("fstr")
        elif genomic:
            fmtr = format_genomic
        elif abs(ran) > 1e5:
            fmtr = ".0e"
        elif abs(ran) < 1e-5:
            fmtr = ".0e"
        elif abs(lbl_delta) >= 1 and all(abs(v - round(v)) < 1e-9 for v in [xmin, xmax, lbl_delta]):
            fmtr = ".0f"
        elif abs(lbl_delta) >= 1:
            fmtr = ".1f"
        else:
            # For small ranges, determine appropriate decimal places
            decimals = max(1, int(-np.floor(np.log10(abs(lbl_delta)))) + 1)
            fmtr = f".{decimals}f"

    ruler, dists = make_ruler(xmin, xmax, num_cols, num_labels=num_labels,
                             ticks=ticks, minor_ticks=minor_ticks, formatter=fmtr,
                             use_minor_symbol_for_ticks=use_minor_symbol)
    sctxt.append(ruler)

    return sctxt, dists

def make_ruler(xmin, xmax, num_cols, num_labels=5, ticks=5, minor_ticks=5, formatter="0.2g", use_minor_symbol_for_ticks=False):
    """
    Create a ruler with evenly spaced labels, ticks, and minor ticks.

    Args:
        xmin: Minimum x value
        xmax: Maximum x value
        num_cols: Number of columns in the plot
        num_labels: Number of labels to place (including endpoints)
        ticks: Number of major ticks per label interval
        minor_ticks: Number of minor ticks per major tick interval
        formatter: Format string or callable for label formatting
        use_minor_symbol_for_ticks: Use minor tick symbol for ticks (when too dense)

    Returns:
        (ruler_string, (label_dist, tick_dist, minor_tick_dist))
    """
    xran = xmax - xmin
    num_labels = max(2, num_labels)

    if isinstance(formatter, str):
        
        if formatter == "genomic":
            formatter = format_genomic
        else:
            frmstr = formatter
            formatter = lambda s: format(s, frmstr)

    # Tick markers
    lbl_marker = "╰"
    tick_marker = "'" if use_minor_symbol_for_ticks else "╵"
    minor_marker = "'"

    # Calculate label positions
    label_positions = []  # (column_pos, x_value, is_last)
    for i in range(num_labels):
        x_val = xmin + (xran * i / (num_labels - 1))
        col_pos = int(round(num_cols * i / (num_labels - 1)))
        is_last = (i == num_labels - 1)
        label_positions.append((col_pos, x_val, is_last))

    # Calculate tick positions (between labels)
    tick_positions = set()
    if ticks > 0:
        for label_idx in range(num_labels - 1):
            start_col = label_positions[label_idx][0]
            end_col = label_positions[label_idx + 1][0]
            interval_cols = end_col - start_col

            # Ensure at least 2 columns per tick
            actual_ticks = min(ticks, interval_cols // 2)

            for t in range(1, actual_ticks + 1):
                tick_col = start_col + int(round(interval_cols * t / (actual_ticks + 1)))
                if tick_col not in [lp[0] for lp in label_positions]:
                    tick_positions.add(tick_col)

    # Calculate minor tick positions (between ticks and labels)
    minor_positions = set()
    if minor_ticks > 0 and not use_minor_symbol_for_ticks:
        # Get all major positions (labels + ticks)
        major_positions = sorted([lp[0] for lp in label_positions] + list(tick_positions))

        for idx in range(len(major_positions) - 1):
            start_col = major_positions[idx]
            end_col = major_positions[idx + 1]
            interval_cols = end_col - start_col

            # Ensure at least 2 columns per minor tick
            actual_minor = min(minor_ticks, interval_cols // 2)

            for m in range(1, actual_minor + 1):
                minor_col = start_col + int(round(interval_cols * m / (actual_minor + 1)))
                if minor_col not in major_positions:
                    minor_positions.add(minor_col)

    # Build the ruler string
    # First, format the final label to know its width
    final_x_val = label_positions[-1][1]
    final_label_str = formatter(final_x_val) + "╯"
    final_label_width = len(final_label_str)
    final_label_start = num_cols - final_label_width + 1

    ruler = []
    col = 0

    while col <= num_cols:
        # Check if we're at the final label position
        if col == final_label_start:
            ruler.append(final_label_str)
            break

        # Check if this column has a non-final label
        label_here = None
        for lp in label_positions[:-1]:  # Exclude last label
            if lp[0] == col:
                label_here = lp
                break

        if label_here:
            col_pos, x_val, _ = label_here
            label_str = lbl_marker + formatter(x_val)
            ruler.append(label_str)
            col += len(label_str)
        elif col in tick_positions:
            ruler.append(tick_marker)
            col += 1
        elif col in minor_positions:
            ruler.append(minor_marker)
            col += 1
        else:
            ruler.append(" ")
            col += 1

    # Calculate actual distances for return value
    label_dist = xran / (num_labels - 1)
    tick_dist = label_dist / (ticks + 1) if ticks > 0 else 0
    minor_tick_dist = tick_dist / (minor_ticks + 1) if minor_ticks > 0 and ticks > 0 else 0

    return "".join(ruler), (label_dist, tick_dist, minor_tick_dist)

