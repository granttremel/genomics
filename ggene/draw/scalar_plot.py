"""
Object-oriented interface for scalar text plotting.

This module provides a ScalarPlot class that wraps the scalar_to_text_* functions
from draw.py, allowing for easy manipulation of plot options and re-rendering.
"""

from typing import List, Optional, Dict, Any, Union

from .colors import Colors, visible_len, visible_slice
from .chars import SCALE, SCALAR_PLOT, SCALE_H, OTHER
from .utils import quantize
from .ruler import Ruler


class ScalarPlot:
    """
    A class for creating and manipulating scalar text plots.

    This class holds scalar data and rendering options, builds output rows,
    and provides methods to adjust visualization parameters and re-render.

    Example usage:
        # Basic plot
        plot = ScalarPlot([1, 2, 3, 4, 5, 3, 2, 1])
        plot.show()

        # Adjust colors and re-render
        plot.set_color_scheme("vscode")
        plot.show()

        # Add a ruler
        plot.add_ruler(xmin=0, xmax=100)
        plot.show()

        # Distribution mode
        dist_data = {'A': 10, 'B': 20, 'C': 15}
        dist_plot = ScalarPlot(dist_data, mode='distribution', labels=True)
        dist_plot.show()
    """

    _max_line = 256

    def __init__(self, scalars: Union[List[float], Dict[str, float]], **kwargs):
        """
        Initialize a ScalarPlot.

        Args:
            scalars: Either a list of scalar values or a dict for distribution mode
            **kwargs: Options for rendering (see default_options for available options)
        """
        self.scalars = scalars
        self.rows = []
        self.ruler = None  # Will hold Ruler instance if ruler is enabled

        # Default options
        self.options = self._get_default_options()

        # Update with user-provided options
        self.options.update(kwargs)
        # self.other_kwargs = {k:v for k,v in kwargs.items() if not k in self.options}
        # print(self.other_kwargs)
        # Auto-detect mode if not specified
        if isinstance(scalars, dict) and 'mode' not in kwargs:
            self.options['mode'] = 'distribution'

        # Initial render
        self.render()

    def _get_default_options(self) -> Dict[str, Any]:
        """Get default options for all rendering modes."""
        return {
            # Common options
            'mode': 'normal',  # 'normal', 'mid', 'distribution'
            'minval': None,
            'maxval': None,
            'fg_color': 89,
            'bg_color': 234,
            'bit_depth': 24,
            'flip': False,
            'effect': None,

            # Normal mode options
            'add_range': False,
            'range_fstr': '0.2f',

            # Mid mode options
            'center': None,
            'rng': None,

            # Distribution mode options
            'key_order': [],
            'labels': False,
            'labelfmt': '',
            'space':0,

            # Ruler options
            'ruler': False,
            'xmin': None,
            'xmax': None,
            'ruler_formatter':"",
            'genomic': False,
            'auto': False,
            'num_labels': 5,
            'ticks': 5,
            'minor_ticks': None,
        }

    def render(self) -> None:
        """
        Render the plot based on current options and data.

        This method builds the output rows using the appropriate scalar_to_text_*
        function based on the mode, then creates a Ruler instance if requested.
        """
        mode = self.options['mode']

        if mode == 'normal':
            self.rows = self._render_normal()
        elif mode == 'mid':
            self.rows = self._render_mid()
        elif mode == 'distribution':
            self.rows = self._render_distribution()
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Create/update ruler if requested
        if self.options['ruler']:
            self._update_ruler()
        else:
            self.ruler = None

    def _render_normal(self) -> List[str]:
        """Render in normal mode using scalar_to_text_nb."""
        kwargs = {
            'minval': self.options['minval'],
            'maxval': self.options['maxval'],
            'fg_color': self.options['fg_color'],
            'bg_color': self.options['bg_color'],
            'bit_depth': self.options['bit_depth'],
            'flip': self.options['flip'],
            'effect': self.options['effect'],
            'add_range': self.options['add_range'],
            'range_fstr': self.options['range_fstr'],
        }
        return scalar_to_text_nb(self.scalars, **kwargs)

    def _render_mid(self) -> List[str]:
        """Render in mid mode using scalar_to_text_mid."""
        kwargs = {
            'center': self.options['center'],
            'rng': self.options['rng'],
            'fg_color': self.options['fg_color'],
            'bg_color': self.options['bg_color'],
            'effect': self.options['effect'],
        }
        
        return scalar_to_text_mid(self.scalars, **kwargs)

    def _render_distribution(self) -> List[str]:
        """Render in distribution mode using scalar_plot_distribution."""
        if not isinstance(self.scalars, dict):
            raise ValueError("Distribution mode requires dict input")

        kwargs = {
            'key_order': self.options['key_order'],
            'bit_depth': self.options['bit_depth'],
            'labels': self.options['labels'],
            'labelfmt': self.options['labelfmt'],
            'add_range': self.options['add_range'],
            'fg_color': self.options['fg_color'],
            'bg_color': self.options['bg_color'],
            'flip': self.options['flip'],
            'effect': self.options['effect'],
            'space':self.options['space']
        }
        # kwargs.update(self.other_kwargs)
        return scalar_plot_distribution(self.scalars, **kwargs)

    def _update_ruler(self) -> None:
        """Create or update the Ruler instance."""
        xmin = self.options['xmin']
        xmax = self.options['xmax']

        if xmin is None or xmax is None:
            raise ValueError("xmin and xmax must be specified for ruler")

        # Calculate number of columns from the plot
        from ggene.draw import SCALE
        num_cols = sum(1 for s in self.rows[0] if s in SCALE) if self.rows else len(self.scalars)

        ruler_kwargs = {
            'genomic': self.options['genomic'],
            'auto': self.options['auto'],
        }

        # Only pass manual parameters if not in auto mode
        if not self.options['auto']:
            ruler_kwargs['formatter'] = self.options['ruler_formatter']
            ruler_kwargs['num_labels'] = self.options['num_labels']
            ruler_kwargs['ticks'] = self.options['ticks']
            if self.options['minor_ticks'] is not None:
                ruler_kwargs['minor_ticks'] = self.options['minor_ticks']

        # Create or update the ruler
        if self.ruler is None:
            self.ruler = Ruler(xmin, xmax, num_cols, **ruler_kwargs)
        else:
            self.ruler.xmin = xmin
            self.ruler.xmax = xmax
            self.ruler.num_cols = num_cols
            self.ruler.options.update(ruler_kwargs)
            self.ruler.render()

    def show(self, plot_label = "", fmt = "") -> None:
        """Print the rendered plot to stdout."""
        if not fmt:
            fmt = "{}"
        for i, row in enumerate(self.rows):
            rlbl = plot_label if i==0 else ""
            print(f"{rlbl.ljust(len(plot_label))}{fmt.format(row)}")
        if self.ruler:
            self.ruler.show()

    def show_chunks(self, chunksz = 256):
        
        line_length = len(self.rows[0])
        num_lines = int(line_length / chunksz) + 1
        if self.options.get("maxval") is None:
            self.options["maxval"] = max(self.scalars)
        
        for n in range(num_lines):
            subdata = self.scalars[chunksz*n:chunksz*(n+1)]
            if len(subdata) < 1:
                break
            subplot = ScalarPlot(subdata, **self.options)
            subplot.show()
            print()
            
    @classmethod
    def show_paired(cls, top_plot, btm_plot, chunksz = None, xlabel = None, center_xlabel = True, suppress = False):
        
        if not chunksz:
            chunksz = len(top_plot.scalars)
        
        line_length = len(top_plot.rows[0])
        
        num_lines = (line_length // chunksz)
        if top_plot.options.get("maxval") is None:
            top_plot.options["maxval"] = max(top_plot.scalars)
        
        top_plot.options["ruler"] = False
        # top_plot.options["fg_color"] = 65
        
        btm_plot.options["flip"] = True
        
        btm_ruler = btm_plot.options.get("ruler", False)
        btm_xmin = btm_plot.options.get("xmin", -1)
        btm_xmax = btm_plot.options.get("xmax", -1)
        
        outlines = []
        
        for n in range(num_lines):
            tsubdata = top_plot.scalars[chunksz*n:chunksz*(n+1)]
            if len(tsubdata) < 1:
                break
            
            tsubplot = ScalarPlot(tsubdata, **top_plot.options)
            tsubplot.render()
            outlines.extend(tsubplot.rows)
            
            if center_xlabel and xlabel:
                if isinstance(xlabel, str):
                    xlabel = [xlabel]
                
                for _xlabel in xlabel:
                    xlbl = visible_slice(_xlabel, start = n*chunksz, stop = (n+1)*chunksz, step = 1)
                    outlines.append(xlbl)
            
            bsubdata = btm_plot.scalars[chunksz*n:chunksz*(n+1)]
            
            bsubplot = ScalarPlot(bsubdata, **btm_plot.options)
            bsubplot.render()
            outlines.extend(bsubplot.rows)
            
            if not center_xlabel and xlabel:
                if isinstance(xlabel, str):
                    xlabel = [xlabel]
                for _xlabel in xlabel:
                    xlbl = visible_slice(_xlabel, start = n*chunksz, stop = (n+1)*chunksz, step = 1)
                    outlines.append(xlbl)
            
        if not suppress:
            for line in outlines:
                print(line)
            print()
        
        return outlines
        
    def set_color(self, fg: Optional[int] = None, bg: Optional[int] = None) -> 'ScalarPlot':
        """
        Set foreground and/or background colors.

        Args:
            fg: Foreground color (0-255)
            bg: Background color (0-255)

        Returns:
            self for method chaining
        """
        if fg is not None:
            self.options['fg_color'] = fg
        if bg is not None:
            self.options['bg_color'] = bg
        self.render()
        return self

    def set_color_scheme(self, scheme_name: str) -> 'ScalarPlot':
        bg, fg = Colors.get_color_scheme(scheme_name)
        return self.set_color(fg=fg, bg=bg)

    def set_bit_depth(self, depth: int) -> 'ScalarPlot':
        self.options['bit_depth'] = depth
        self.render()
        return self

    def set_flip(self, flip: bool) -> 'ScalarPlot':
        self.options['flip'] = flip
        self.render()
        return self

    def set_range(self, minval: Optional[float] = None, maxval: Optional[float] = None) -> 'ScalarPlot':
        if minval is not None:
            self.options['minval'] = minval
        if maxval is not None:
            self.options['maxval'] = maxval
        self.render()
        return self

    def enable_range_labels(self, enable: bool = True, fstr: str = '0.2f') -> 'ScalarPlot':
        self.options['add_range'] = enable
        self.options['range_fstr'] = fstr
        self.render()
        return self

    def enable_ruler(
        self,
        xmin: float,
        xmax: float,
        genomic: bool = False,
        auto: bool = False,
        multi_row_labels: bool = False,
        max_label_rows: int = 4,
        num_labels: int = 5,
        ticks: int = 5,
        minor_ticks: Optional[int] = None
    ) -> 'ScalarPlot':
        """
        Enable and configure a ruler.

        Args:
            xmin: Minimum x value
            xmax: Maximum x value
            genomic: Use genomic formatting (kb, Mb)
            auto: Automatically determine optimal ruler parameters
            multi_row_labels: Use multiple rows to avoid label collisions
            max_label_rows: Maximum number of label rows
            num_labels: Number of labels on ruler (ignored if auto=True)
            ticks: Number of major ticks (ignored if auto=True)
            minor_ticks: Number of minor ticks (ignored if auto=True)

        Returns:
            self for method chaining
        """
        self.options['ruler'] = True
        self.options['xmin'] = xmin
        self.options['xmax'] = xmax
        self.options['genomic'] = genomic
        self.options['auto'] = auto
        if not auto:
            self.options['num_labels'] = num_labels
            self.options['ticks'] = ticks
            if minor_ticks is not None:
                self.options['minor_ticks'] = minor_ticks
        self.render()

        # Enable multi-row labels if requested
        if self.ruler and multi_row_labels:
            self.ruler.set_multi_row_labels(True, max_label_rows)

        return self

    def disable_ruler(self) -> 'ScalarPlot':
        self.options['ruler'] = False
        self.render()
        return self

    def set_mode(self, mode: str) -> 'ScalarPlot':
        """
        mode: 'normal', 'mid', or 'distribution'
        """
        if mode not in ['normal', 'mid', 'distribution']:
            raise ValueError(f"Invalid mode: {mode}. Must be 'normal', 'mid', or 'distribution'")
        self.options['mode'] = mode
        self.render()
        return self

    def set_effect(self, effect: Optional[str]) -> 'ScalarPlot':
        """
        ('border', 'bold', 'dim', etc.) or None
        """
        self.options['effect'] = effect
        self.render()
        return self

    def set_option(self, **kwargs) -> 'ScalarPlot':
        """
        Generic option setter for any option.

        Args:
            **kwargs: Options to update

        Returns:
            self for method chaining
        """
        self.options.update(kwargs)
        self.render()
        return self

    def get_rows(self) -> List[str]:
        """
        Get the rendered rows without printing (includes ruler if present).

        Returns:
            List of rendered row strings
        """
        all_rows = self.rows.copy()
        if self.ruler:
            all_rows.extend(self.ruler.get_rows())
        return all_rows

    def __repr__(self) -> str:
        mode = self.options['mode']
        if isinstance(self.scalars, dict):
            data_repr = f"dict with {len(self.scalars)} keys"
        else:
            data_repr = f"list of {len(self.scalars)} values"
        return f"ScalarPlot({data_repr}, mode='{mode}')"

    def __str__(self) -> str:
        all_rows = self.rows.copy()
        if self.ruler:
            all_rows.extend(self.ruler.get_rows())
        return '\n'.join(all_rows)


def scalar_to_text_nb(scalars, minval = None, maxval = None, fg_color = 53, bg_color = 234, bit_depth = 24, flip = False, effect = None, add_range = False, **kwargs):
    """
    ok hear me out: plot two quantities together by putting the larger "behind" the smaller using bg. e.g.:
    wait this doesn't work, two idp quantities cant share the same cell
    {bg0fg0}██{bg=fg2 fg0}█
    
    okay how about: denser labels using longer lines, like
        ▁▃▂▂█▅
       1╯││││╰6
        4╯││╰12   
         2╯╰3
    
    """
    
    if flip:
        bg, fg = Colors.get_colors_fgbg(bg_color, fg_color)
    else:
        fg, bg = Colors.get_colors_fgbg(fg_color, bg_color)
    
    add_border = False
    eff = ""
    if effect:
        if effect == "border":
            add_border = True
        else:
            eff += str(effect)
    
    
    base_bit_depth = len(SCALE) - 1
    if not bit_depth % base_bit_depth == 0:
        return ["no"]
    
    nrows = bit_depth // base_bit_depth
    ncols = len(scalars)
    nvals = base_bit_depth * nrows
    
    rows = [[fg+bg+eff] for r in range(nrows)]
    
    bit_ranges = [base_bit_depth*i for i in range(nrows)]
    
    if minval is None:
        minval = min(scalars)
    if maxval is None:
        maxval = max(scalars)
    rng = (maxval - minval)/1
    c = (minval+ maxval)/2
    if rng == 0:
        rng = 2
        c = minval
    # print([format(sc, "0.2f") for sc in scalars])
    # print(f"scalar to text with {minval}, {maxval}, {rng}, {c}")
    # print(f"data with {min(scalars)}, {max(scalars)}")
    
    for s in scalars:
        sv = int(nvals*((s - c)/rng)) + bit_depth // 2
            
        for row, bit_range in zip(rows, bit_ranges):
            if sv < bit_range:
                sym = SCALE[0]
            elif sv >= bit_range + base_bit_depth:
                sym = SCALE[-1]
            else:
                ssv = sv % base_bit_depth
                sym = SCALE[ssv]
            row.append(sym)
    
    brdc = "\x1b[38;5;6m"
    outstrs= []
    for row in rows[::-1]:
        if add_border:
            row.insert(0, brdc + OTHER.get("right_eighth",""))
            row.append(brdc + SCALE_H[1])
        row.append(Colors.RESET)
        outstrs.append("".join(row))
    
    if add_border:
        outstrs.insert(0, " " + SCALE[1]*ncols + " ")
        outstrs.append(f"{brdc} " + OTHER.get("upper_eighth","")*ncols + f" {Colors.RESET}")
    
    if add_range:
        ran_fstr = kwargs.get("range_fstr", "0.2f")
        hi, lo = SCALAR_PLOT.get("range_hi"), SCALAR_PLOT.get("range_lo")
        hi, lo = (lo, hi) if flip else (hi, lo)
        minstr = format(minval, ran_fstr)[:4]
        maxstr = format(maxval, ran_fstr)[:4]
        outstrs[0] += hi + maxstr
        if bit_depth > 8:
            outstrs[-1] += lo + minstr
        
    if flip:
        outstrs = flip_scalar_text(outstrs)
    
    return outstrs

def scalar_to_text_mid(scalars, center = None, rng = None, fg_color = 53, bg_color = 234,  effect = None):
    
    bit_depth = 16
    
    ifg, ibg = Colors.get_colors_fgbg(fg_color, bg_color)
    bg, fg = Colors.get_colors_fgbg(bg_color, fg_color)
    
    eff = ""
    if effect:
        eff += str(effect)
    
    
    base_bit_depth = len(SCALE) - 1
    if not bit_depth % base_bit_depth == 0:
        return ["no"]
    
    nrows = bit_depth // base_bit_depth
    ncols = len(scalars)
    nvals = base_bit_depth * nrows
    
    rows = [[fg+bg+eff],[ifg+ibg+eff]]
    
    bit_ranges = [base_bit_depth*i - bit_depth/2 for i in range(nrows)]
    
    if not center:
        c = 0
    else:
        c = center
    
    if not rng:
        minval = min(scalars)
        maxval  = max(scalars)
        rng = 2*max(abs(minval), abs(maxval))
    minval, maxval = c-rng/2, c+rng/2
    
    neg = False
    for s in scalars:
        sv = int(nvals*((s - c)/rng))
        if sv < 0 and not neg:
            neg = True
        elif sv >= 0 and neg:
            neg = False
        
        for row, bit_range in zip(rows, bit_ranges):
            if sv < bit_range:
                sym = SCALE[0]
            elif sv >= bit_range + base_bit_depth:
                sym = SCALE[-1]
            else:
                ssv = sv % base_bit_depth
                sym = SCALE[ssv]
            row.append(sym)
    
    outstrs= []
    for row in rows[::-1]:
        row.append(Colors.RESET)
        outstrs.append("".join(row))
        
    return outstrs


def scalar_to_text_8b(scalars, minval = None, maxval = None, fg_color = 53, bg_color = 234, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 8, flip = flip)

def scalar_to_text_16b(scalars, minval = None, maxval = None, fg_color = 53, bg_color = 234, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 16, flip = flip)


def scalar_plot_distribution(dist, key_order = [], bit_depth = 8, labels = False, labelfmt = "", add_range = False, **kwargs):
    
    if not key_order:
        key_order = sorted(dist.keys())
    scalars = [dist[ks] for ks in key_order]
    space = kwargs.get("space", 0)
    if space:
        _scalars = [0]
        for sc in scalars:
            _scalars.append(sc)
            _scalars.extend([0]*space)
        scalars = _scalars
    
    _ = kwargs.pop("minval",None)
    sctxt = scalar_to_text_nb(scalars, bit_depth = bit_depth, minval = 0, add_range = add_range, **kwargs)
    
    if labels:
        if not labelfmt:
            labelfmtr = lambda a:str(a)[0].upper()
        else:
            labelfmtr = lambda a:format(a,labelfmt)
        lbls = [] + [" "]*space
        for k in key_order:
            lbls.append(labelfmtr(k))
            lbls.extend([" "]*space)
        sctxt.append("".join(lbls))
    
    return sctxt


def flip_scalar_text(sctext):
    
    nsyms = len(SCALE)
    scale_inv = {SCALE[i]:SCALE[nsyms-i-1] for i in range(nsyms)}
    
    out = []
    for row in sctext[::-1]:
        newrow = []
        for sym in row:
            newrow.append(scale_inv.get(sym, sym))
        out.append("".join(newrow))
    return out

def hflip_scalar_text(sctext):
    
    nsyms = len(SCALE_H)
    scale_inv = {SCALE_H[i]:SCALE_H[nsyms-i-1] for i in range(nsyms)}
    
    out = []
    for row in sctext[::-1]:
        newrow = []
        for sym in row:
            newrow.append(scale_inv.get(sym, sym))
        out.append("".join(newrow))
    return out


def scalar_to_text_nbh(scalars, minval = None, maxval = None, fg_color = 7, bg_color = 212, bit_depth = 24, flip = False, effect = None):
    
    len_sc = len(scalars)
    
    if flip:
        fg = f"\x1b[48;5;{fg_color}m"
        bg = f"\x1b[38;5;{bg_color}m"
    else:
        fg = f"\x1b[38;5;{fg_color}m"
        bg = f"\x1b[48;5;{bg_color}m"
    
    eff = ""
    if effect:
        eff += str(effect)
    
    base_bit_depth = len(SCALE_H) - 1
    if not bit_depth % base_bit_depth == 0:
        return ["no"]
    
    ncols = bit_depth // base_bit_depth
    nvals = base_bit_depth * ncols
    
    rows = [[fg+bg+eff] for r in range(len_sc)]
    
    bit_ranges = [base_bit_depth*i for i in range(ncols)]
    
    if not minval:
        minval = min(scalars)
    if not maxval:
        maxval = max(scalars)
    rng = max(1, (maxval - minval)/1)
    c = (minval+ maxval)/2
    
    for s, row in zip(scalars, rows):
        sv = int(nvals*(s - c)/rng) + bit_depth // 2
        for bit_range in bit_ranges:
            if sv < bit_range:
                sym = SCALE_H[0]
            elif sv >= bit_range + base_bit_depth:
                sym = SCALE_H[-1]
            else:
                ssv = sv % base_bit_depth
                sym = SCALE_H[ssv]
            row.append(sym)
    
    outstrs= []
    for row in rows:
        row.append(Colors.RESET)
        border_row = ["▁" if r in SCALE else r for r in row]
        outstrs.append("".join(row))
        # outstrs.append("".join(border_row))
        # print(outstrs[-1])
        
    if flip:
        return hflip_scalar_text(outstrs)
    else:
        return outstrs