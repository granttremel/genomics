"""
Heatmap visualization for terminal display.

Provides both functional and object-oriented interfaces for creating
color-coded heatmaps in the terminal using Unicode block characters.
"""

from typing import List, Optional, Dict, Any, Union, Tuple
from dataclasses import dataclass
import numpy as np

from .colors import Colors
from .chars import SCALE, OTHER
from .ruler import make_ruler, Ruler


class Heatmap:
    """
    A class for creating and manipulating terminal heatmaps.

    Supports both standard mode (1 data row = 1+ terminal rows) and
    half-block mode (2 data rows = 1 terminal row) for doubled vertical resolution.

    Example usage:
        # Basic heatmap
        data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        hm = Heatmap(data)
        hm.show()

        # Half-block mode for compact display
        hm = Heatmap(data, half_block=True)
        hm.show()

        # With labels and colorbar
        hm = Heatmap(data, row_labels=['A', 'B', 'C'], colorbar=True)
        hm.show()
    """

    def __init__(self, data: Union[List[List[float]], np.ndarray], **kwargs):
        """
        Initialize a Heatmap.

        Args:
            data: 2D array of values (list of lists or numpy array)
            **kwargs: Rendering options (see _get_default_options)
        """
        # Convert numpy array to list
        if isinstance(data, np.ndarray):
            self.data = data.tolist()
        else:
            self.data = data

        self.rows: List[str] = []
        self.ruler: Optional[Ruler] = None

        # Default options
        self.options = self._get_default_options()
        self.options.update(kwargs)

        # Validate data
        self.nrows = len(self.data)
        self.ncols = len(self.data[0]) if self.nrows > 0 else 0

        # Calculate data range
        self._calculate_range()

        # Build color scale
        self._build_color_scale()

        # Initial render
        self.render()

    def _get_default_options(self) -> Dict[str, Any]:
        """Get default rendering options."""
        return {
            # Data range
            'minval': None,
            'maxval': None,
            'center': None,  # If set, creates symmetric range around center

            # Colors
            'colors': None,
            'color_scheme': 'terra',
            'add_middle':True,
            'brightness': 1.0,
            'contrast': 1.0,
            'hue_shift': 1.0,
            'saturation':1.0,
            'symmetric_color': True,
            'num_colors': 24,

            # Layout
            'col_width': 1,       # Characters per column
            'col_space': 0,       # Space between columns
            'row_height': 1,      # Terminal rows per data row (ignored in half_block mode)
            'row_space': 0,       # Space between rows (minimum 1 in half_block mode)

            # Half-block mode: 2 data rows per terminal row
            'half_block': False,

            # Labels
            'row_labels': None,
            'col_labels': None,

            # Value labels (only for non-half-block mode)
            'show_values': False,
            'value_fmt': '0.2f',
            'value_labels': None,  # Pre-formatted labels (2D list/array)

            # Colorbar
            'colorbar': False,
            'colorbar_width': 40,

            # Ruler
            'ruler': False,
            'xmin': None,
            'xmax': None,
        }

    def _calculate_range(self):
        """Calculate min/max values from data."""
        flat_data = [val for row in self.data for val in row if val is not None and not np.isnan(val)]

        if not flat_data:
            self._minval = 0
            self._maxval = 1
            return

        data_min = min(flat_data)
        data_max = max(flat_data)

        # Use provided values or calculate from data
        self._minval = self.options['minval'] if self.options['minval'] is not None else data_min
        self._maxval = self.options['maxval'] if self.options['maxval'] is not None else data_max

        # If center is specified, make range symmetric
        if self.options['center'] is not None:
            center = self.options['center']
            sym_range = max(abs(self._maxval - center), abs(self._minval - center))
            self._minval = center - sym_range
            self._maxval = center + sym_range

        self._range = self._maxval - self._minval
        if self._range == 0:
            self._range = 1

    def _build_color_scale(self):
        """Build the color scale for value mapping."""
        
        if self.options['colors']:
            cols = self.options['colors']
        else:
            cols = Colors.get_color_scheme_24b(
                self.options['color_scheme'], 
                brightness = self.options['brightness'],
                contrast = self.options['contrast'],
                hue_shift = self.options['hue_shift'],
                saturation = self.options['saturation'])
        
        if self.options['add_middle']:
            cols = Colors.add_middle(cols[0], cols[-1], brightness = 0.3, saturation = 0.5)
        
        colors = Colors()
        self._color_scale = colors.get_color_scale_24b(self.options['num_colors'], *cols)

    def _value_to_color(self, val: float) -> Tuple[int, int, int]:
        """Map a value to an RGB color tuple."""
        if val is None or np.isnan(val):
            return (20, 20, 20)  # Black for None values

        normalized = (val - self._minval) / self._range
        color_idx = int(normalized * (len(self._color_scale) - 1))
        color_idx = max(0, min(len(self._color_scale) - 1, color_idx))
        return self._color_scale[color_idx]

    def render(self) -> None:
        """Render the heatmap based on current options."""
        self.rows = []

        row_labels = self.options['row_labels'] or []
        col_labels = self.options['col_labels']

        max_row_label_len = max(len(str(lbl)) for lbl in row_labels) if row_labels else 0
        row_fmt = "{:<%s}{}" % str(max_row_label_len + 2)

        col_width = self.options['col_width']
        col_space = self.options['col_space']

        # Render column labels if provided
        if col_labels:
            self._render_col_labels(col_labels, row_fmt, col_width, col_space, max_row_label_len)

        # Render data rows
        if self.options['half_block']:
            self._render_half_block(row_labels, row_fmt, col_width, col_space, max_row_label_len)
        else:
            self._render_standard(row_labels, row_fmt, col_width, col_space, max_row_label_len)

        # Render ruler if requested
        if self.options['ruler']:
            self._render_ruler(max_row_label_len, col_width, col_space)

        # Render colorbar if requested
        if self.options['colorbar']:
            self._render_colorbar(max_row_label_len)

    def _render_col_labels(self, col_labels, row_fmt, col_width, col_space, max_row_label_len):
        """Render column label rows."""
        num_col_lbl_rows = len(col_labels[0]) if isinstance(col_labels[0], list) else 1

        for n in range(num_col_lbl_rows):
            header_items = ""
            for col_lbls in col_labels:
                col_lbl = col_lbls[n] if isinstance(col_lbls, list) else col_lbls
                short_lbl = str(col_lbl).center(col_width + col_space)
                header_items += short_lbl

            header_line = row_fmt.format("", header_items)  # Empty string gets padded by row_fmt
            self.rows.append(header_line)

    def _render_standard(self, row_labels, row_fmt, col_width, col_space, max_row_label_len):
        """Render in standard mode (1+ terminal rows per data row)."""
        row_height = self.options['row_height']
        row_space = self.options['row_space']
        show_values = self.options['show_values']
        value_fmt = self.options['value_fmt']
        value_labels = self.options['value_labels']

        for i, data_row in enumerate(self.data):
            row_lbl = str(row_labels[i]).ljust(max_row_label_len) if row_labels and i < len(row_labels) else ""

            # Generate value labels for this row if show_values is enabled
            value_row = None
            if show_values:
                if value_labels and i < len(value_labels):
                    value_row = value_labels[i]
                else:
                    value_row = [format(v, value_fmt) if v is not None else "" for v in data_row]

            for nr in range(row_height):
                rlbl = row_lbl if nr == 0 else ""
                vr = value_row if nr == 0 else None  # Only show values on first row

                # Use thin border for fractional row_space
                use_border = row_space > 0 and row_space < 1 and nr == row_height - 1

                line = self._make_row(data_row, col_width, col_space, row_border=use_border, value_labels=vr)
                self.rows.append(row_fmt.format(rlbl, line))

            # Add full row spaces
            if row_space >= 1:
                for _ in range(int(row_space)):
                    self.rows.append("")

    def _render_half_block(self, row_labels, row_fmt, col_width, col_space, max_row_label_len):
        """
        Render in half-block mode (2 data rows per terminal row).

        Uses upper half block (▀) with:
        - Foreground color = upper data row
        - Background color = lower data row
        """
        row_space = max(1, self.options['row_space'])  # Minimum 1 in half-block mode
        half_block = OTHER.get("upper_half", "▀")

        # Process data rows in pairs
        for i in range(0, self.nrows, 2):
            upper_row = self.data[i]
            lower_row = self.data[i + 1] if i + 1 < self.nrows else [None] * self.ncols

            # Row label from upper row
            row_lbl = ""
            if row_labels:
                if i < len(row_labels):
                    upper_lbl = str(row_labels[i])
                    lower_lbl = str(row_labels[i + 1]) if i + 1 < len(row_labels) else ""
                    row_lbl = f"{upper_lbl} {lower_lbl}".ljust(max_row_label_len)

            line = self._make_half_block_row(upper_row, lower_row, col_width, col_space, half_block)
            self.rows.append(row_fmt.format(row_lbl, line))

            # Add row spaces (minimum 1 for half-block mode)
            for _ in range(int(row_space) - 1):
                self.rows.append("")

    def _make_row(self, data_row: List[float], col_width: int, col_space: int,
                  row_border: bool = False, value_labels: Optional[List[str]] = None) -> str:
        """Create a single heatmap row string."""
        block = " " if not row_border else OTHER.get("lower_eighth", "▁")
        space = SCALE[-1]

        line = ""
        bg_code = Colors(fg=234).fg_code

        for j, val in enumerate(data_row):
            if val is None:
                line += " " * (col_width + col_space)
                continue

            color = self._value_to_color(val)
            fg_code = Colors(bg=color).bg_code

            # Use value label if provided, otherwise use block
            if value_labels and j < len(value_labels):
                val_str = value_labels[j]
                col_str = str(val_str)[:col_width].center(col_width)
            else:
                col_str = block * col_width

            space_str = space * col_space

            line += bg_code + fg_code + col_str + space_str + Colors.RESET

        return line

    def _make_half_block_row(self, upper_row: List[float], lower_row: List[float],
                              col_width: int, col_space: int, half_block: str) -> str:
        """
        Create a half-block row combining two data rows.

        Upper row data -> foreground color (top half of ▀)
        Lower row data -> background color (bottom half of ▀)
        """
        line = ""

        for j in range(len(upper_row)):
            upper_val = upper_row[j] if j < len(upper_row) else None
            lower_val = lower_row[j] if j < len(lower_row) else None

            if upper_val is None and lower_val is None:
                line += " " * (col_width + col_space)
                continue

            # Upper row = foreground (top half of ▀)
            # Lower row = background (bottom half of ▀)
            upper_color = self._value_to_color(upper_val) if upper_val is not None else (0, 0, 0)
            lower_color = self._value_to_color(lower_val) if lower_val is not None else (0, 0, 0)

            fg_code = Colors(fg=upper_color).fg_code
            bg_code = Colors(bg=lower_color).bg_code

            col_str = half_block * col_width
            space_str = " " * col_space

            line += fg_code + bg_code + col_str + Colors.RESET + space_str

        return line

    def _render_ruler(self, max_row_label_len: int, col_width: int, col_space: int):
        """Render the ruler row."""
        xmin = self.options.get('xmin', 0)
        xmax = self.options.get('xmax', self.ncols)

        num_cols = self.ncols * (col_width + col_space)
        ruler_line, _ = make_ruler(xmin, xmax, num_cols, num_labels=16, ticks=0, minor_ticks=0)

        row_fmt = "{:<%s}{}" % str(max_row_label_len + 1)
        self.rows.append(row_fmt.format("", ruler_line))

    def _render_colorbar(self, max_row_label_len: int):
        """Render the colorbar."""
        self.rows.append("")

        colorbar_len = self.options['colorbar_width']
        colorbar_line = " " * (max_row_label_len + 2)

        for i in range(colorbar_len):
            normalized = i / (colorbar_len - 1)
            val = self._minval + normalized * self._range
            color = self._value_to_color(val)
            fg = Colors.get_color(color)
            colorbar_line += fg + SCALE[-1] + Colors.RESET

        self.rows.append(colorbar_line)

        # Colorbar labels
        value_fmt = self.options['value_fmt']
        label_line = " " * (max_row_label_len + 2)
        label_line += format(self._minval, value_fmt).ljust(colorbar_len // 2)
        label_line += format(self._maxval, value_fmt).rjust(colorbar_len // 2)
        self.rows.append(label_line)

    def show(self, suppress: bool = False) -> List[str]:
        """Print the heatmap to stdout."""
        if not suppress:
            for row in self.rows:
                print(row)
            print()
        return self.rows

    def get_rows(self) -> List[str]:
        """Get rendered rows without printing."""
        return self.rows.copy()

    # === Chainable setters ===

    def set_range(self, minval: Optional[float] = None, maxval: Optional[float] = None) -> 'Heatmap':
        """Set data range for color mapping."""
        if minval is not None:
            self.options['minval'] = minval
        if maxval is not None:
            self.options['maxval'] = maxval
        self._calculate_range()
        self._build_color_scale()
        self.render()
        return self

    def set_center(self, center: float) -> 'Heatmap':
        """Set center value for symmetric color scale."""
        self.options['center'] = center
        self._calculate_range()
        self._build_color_scale()
        self.render()
        return self

    def set_color_scheme(self, scheme: str, brightness = None, contrast = None, hue = None) -> 'Heatmap':
        """Set color scheme."""
        self.options['color_scheme'] = scheme
        self._build_color_scale()
        self.render()
        return self

    def set_half_block(self, enabled: bool = True) -> 'Heatmap':
        """Enable or disable half-block mode."""
        self.options['half_block'] = enabled
        self.render()
        return self

    def set_layout(self, col_width: int = None, col_space: int = None,
                   row_height: int = None, row_space: float = None) -> 'Heatmap':
        """Set layout parameters."""
        if col_width is not None:
            self.options['col_width'] = col_width
        if col_space is not None:
            self.options['col_space'] = col_space
        if row_height is not None:
            self.options['row_height'] = row_height
        if row_space is not None:
            self.options['row_space'] = row_space
        self.render()
        return self

    def set_labels(self, row_labels: List[str] = None, col_labels: List = None) -> 'Heatmap':
        """Set row and/or column labels."""
        if row_labels is not None:
            self.options['row_labels'] = row_labels
        if col_labels is not None:
            self.options['col_labels'] = col_labels
        self.render()
        return self

    def set_value_labels(self, show_values: bool = True, value_fmt: str = None,
                        value_labels: List[List[str]] = None) -> 'Heatmap':
        """Enable/configure value labels (only works in non-half-block mode).

        Args:
            show_values: Whether to show values in cells
            value_fmt: Format string for values (e.g., '0.2f', '0.1e')
            value_labels: Pre-formatted 2D list of labels (overrides value_fmt)
        """
        if self.options['half_block']:
            print("Warning: value labels only supported in non-half-block mode")
            return self

        self.options['show_values'] = show_values
        if value_fmt is not None:
            self.options['value_fmt'] = value_fmt
        if value_labels is not None:
            self.options['value_labels'] = value_labels
        self.render()
        return self

    def enable_colorbar(self, enabled: bool = True, width: int = 40) -> 'Heatmap':
        """Enable or disable colorbar."""
        self.options['colorbar'] = enabled
        self.options['colorbar_width'] = width
        self.render()
        return self

    def enable_ruler(self, xmin: float, xmax: float) -> 'Heatmap':
        """Enable ruler with x-axis range."""
        self.options['ruler'] = True
        self.options['xmin'] = xmin
        self.options['xmax'] = xmax
        self.render()
        return self

    def set_option(self, **kwargs) -> 'Heatmap':
        """Generic option setter."""
        self.options.update(kwargs)
        self._calculate_range()
        self._build_color_scale()
        self.render()
        return self

    def __repr__(self) -> str:
        return f"Heatmap({self.nrows}x{self.ncols}, half_block={self.options['half_block']})"

    def __str__(self) -> str:
        return '\n'.join(self.rows)


# === Legacy functional interface below ===

def make_heatmap(data, row_labels=None, col_labels=None, minval=None, maxval=None, center=0, color_scheme = "terra",
            show_values=False, value_fmt="0.2f", value_labels = None, colorbar=True, suppress = False, **kwargs):

    out_rows = []

    if isinstance(data, np.ndarray):
        data = data.tolist()

    min_color, max_color = Colors.get_color_scheme_24b(color_scheme)
    mid_color, _ = Colors.get_color_scheme_24b("vscode")
    
    nrows = len(data)
    ncols = len(data[0]) if nrows > 0 else 0

    if nrows == 0 or ncols == 0:
        print("Empty data")
        return

    flat_data = [val for row in data for val in row if val is not None]

    if minval is None:
        minval = min(flat_data)
    if maxval is None:
        maxval = max(flat_data)
    
    rng = maxval - minval
    if rng == 0:
        rng = 1
    
    colors = Colors()

    num_colors = kwargs.get("num_colors", 24)
    symmetric_color = kwargs.get("symmetric_color", True)
    if center is not None:
        sym_range = max(abs(maxval - center), abs(minval - center))
        minval = center - sym_range
        maxval = center + sym_range

    else:
        if min_color is None:
            min_color = 16
        if max_color is None:
            max_color = 226
    
    if symmetric_color:
        color_scale = colors.get_color_scale_24b(num_colors//2, min_color, mid_color, max_color)
    else:
        color_scale = colors.get_color_scale_24b(num_colors, min_color, max_color)
    
    
    def value_to_color(val):
        normalized = (val - minval) / rng
        color_idx = int(normalized * (len(color_scale) - 1))
        color_idx = max(0, min(len(color_scale) - 1, color_idx))
        return color_scale[color_idx]

    if row_labels is None:
        row_labels = []
    if col_labels is None:
        col_labels = []
    
    max_row_label_len = max(len(str(lbl)) for lbl in row_labels) if row_labels else 0

    row_frm = "{:<%s}{}" % str(max_row_label_len + 1)
    
    col_width = kwargs.get("col_width", 2)
    col_space = kwargs.get("col_space", 1)
    row_height = kwargs.get("row_height", 1)
    row_space = kwargs.get("row_space", 0.5)

    if col_labels:
        
        num_col_lbl_rows = len(col_labels[0]) if isinstance(col_labels[0], list) else 1
        
        for n in range(num_col_lbl_rows):
            header_items = ""
            for i, col_lbls in enumerate(col_labels):
                col_lbl = col_lbls[n]
                short_lbl = str(col_lbl).center(col_width+col_space)
                header_items += short_lbl
            
            header_line = row_frm.format(" ", header_items)
            
            out_rows.append(header_line)
    
    for i, row in enumerate(data):
        row_lbl = str(row_labels[i]).ljust(max_row_label_len) if row_labels else ""
        line = ""
        
        value_row = None
        if show_values:
            if value_labels:
                value_row = value_labels[i]
            else:
                value_row = [value_fmt.format(v) for v in row]
        
        for nr in range(row_height):
            
            vr = value_row if nr == 0 else None
            rlbl = row_lbl if nr == 0 else ""
            
            use_border = False
            if row_space > 0 and row_space < 1 and nr == row_height - 1:
                use_border = True
            
            line = make_heatmap_row(row, value_to_color, value_labels = vr, col_width = col_width, col_space = col_space, row_border = use_border)
            out_rows.append(row_frm.format(rlbl, line))
        
        if row_space >= 1:
            for rs in range(row_space):
                out_rows.append("")
    
    if kwargs.get("ruler"):
        
        xmin = kwargs.get("xmin", 0)
        xmax = kwargs.get("xmax", len(data[0]))
        
        num_cols = len(data[0])*(col_width + col_space)
        
        ruler, _ = make_ruler(xmin, xmax, num_cols, num_labels = 16, ticks = 0, minor_ticks = 0)
        out_rows.append(row_frm.format("", ruler))
    
    if colorbar:
        out_rows.append("")
        cb_rows = make_heatmap_colorbar(max_row_label_len, value_to_color, minval, maxval, rng=rng, value_fmt = value_fmt)
        out_rows.extend(cb_rows)
        # colorbar_len = 40
        # colorbar_line = " " * (max_row_label_len + 2)

        # for i in range(colorbar_len):
        #     normalized = i / (colorbar_len - 1)
        #     val = minval + normalized * rng
        #     color_code = value_to_color(val)
        #     fg = Colors.get_color(color_code)
        #     colorbar_line += fg + SCALE[-1] + Colors.RESET

        # out_rows.append(colorbar_line)

        # # Colorbar labels
        # label_line = " " * (max_row_label_len + 2)
        # label_line += format(minval, value_fmt).ljust(colorbar_len // 2)
        # label_line += format(maxval, value_fmt).rjust(colorbar_len // 2)
        # out_rows.append(label_line)
    
    if not suppress:
        for row in out_rows:
            print(row)
        print()
    
    return out_rows

def heatrow(data, label = None, col_labels = None, **kwargs):
    
    kwargs["colorbar"] = False
    row = make_heatmap([data], row_labels = [label], col_labels = col_labels, **kwargs)
    
    return row
    

def make_heatmap_colorbar(max_row_label_len, value_to_color, minval, maxval, rng = None, value_fmt = "0.3f"):
    
    if not rng:
        rng = maxval - minval
    
    colorbar_len = 40
    colorbar_line = " " * (max_row_label_len + 2)

    colorbar = []

    for i in range(colorbar_len):
        normalized = i / (colorbar_len - 1)
        val = minval + normalized * rng
        color_code = value_to_color(val)
        fg = Colors.get_color(color_code)
        colorbar_line += fg + SCALE[-1] + Colors.RESET

    colorbar.append(colorbar_line)

    # Colorbar labels
    label_line = " " * (max_row_label_len + 2)
    label_line += format(minval, value_fmt).ljust(colorbar_len // 2)
    label_line += format(maxval, value_fmt).rjust(colorbar_len // 2)
    colorbar.append(label_line)
    
    return colorbar


def make_heatmap_row(data_row, value_to_color, value_labels = None, col_width = 1, col_space = 0, row_border = False):
    
    block = " "
    space = SCALE[-1]
    if row_border:
        block = OTHER.get("lower_eighth")
    # col_border = OTHER.get("right_eighth")
    
    line = ""
    # bg = Colors.get_color(234, background = False)
    bg = Colors(fg = 234).fg_code
    
    if col_space > 0 and col_space < 1:
        col_width -= 1
        print(col_space, col_width)
    
    for j, val in enumerate(data_row):
        if val is None:
            line += " "*(col_width + col_space)
            continue
        
        color_code = value_to_color(val)
        # fg = Colors.get_color(color_code, background = True)
        fg = Colors(bg = color_code).bg_code

        if value_labels:
            val_str = value_labels[j]
            col_str = val_str[:2].center(col_width)
        else:
            col_str = block * col_width
        
        # if col_space > 0 and col_space < 1:
        #     space_str = col_border
        # else:
        space_str = space * col_space
        
        line += bg + fg + col_str + space_str + Colors.RESET
        
    return line
    