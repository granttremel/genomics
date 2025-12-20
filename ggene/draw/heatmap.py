

import numpy as np

from .colors import Colors
from .chars import SCALE, OTHER
from .ruler import make_ruler

def heatmap(data, row_labels=None, col_labels=None, minval=None, maxval=None, center=0, color_scheme = "terra",
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
        cs1 = colors.get_color_scale_24b(mid_color, min_color, num_colors//2)
        cs2 = colors.get_color_scale_24b(mid_color, max_color, num_colors//2)
        color_scale = cs1[::-1] + cs2
    else:
        color_scale = colors.get_color_scale_24b(min_color, max_color, num_colors)
    
    
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
        colorbar_len = 40
        colorbar_line = " " * (max_row_label_len + 2)

        for i in range(colorbar_len):
            normalized = i / (colorbar_len - 1)
            val = minval + normalized * rng
            color_code = value_to_color(val)
            fg = Colors.get_color(color_code)
            colorbar_line += fg + SCALE[-1] + Colors.RESET

        out_rows.append(colorbar_line)

        # Colorbar labels
        label_line = " " * (max_row_label_len + 2)
        label_line += format(minval, value_fmt).ljust(colorbar_len // 2)
        label_line += format(maxval, value_fmt).rjust(colorbar_len // 2)
        out_rows.append(label_line)
    
    if not suppress:
        for row in out_rows:
            print(row)
        print()
    
    return out_rows

def make_heatmap_row(data_row, value_to_color, value_labels = None, col_width = 1, col_space = 0, row_border = False):
    
    block = " "
    space = SCALE[-1]
    if row_border:
        block = OTHER.get("lower_eighth")
    # col_border = OTHER.get("right_eighth")
    
    line = ""
    bg = Colors.get_color(234, bg = False)
    
    if col_space > 0 and col_space < 1:
        col_width -= 1
        print(col_space, col_width)
    
    for j, val in enumerate(data_row):
        if val is None:
            line += " "*(col_width + col_space)
            continue
        
        color_code = value_to_color(val)
        fg = Colors.get_color(color_code, bg = True)

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
    