
from typing import List
import math
import random
import regex
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from ggene.utils import utils
from ggene.seqs.bio import reverse_complement, complement
from ggene.seqs.find import find_subsequence, find_subsequences

SCALE = " ▁▂▃▄▅▆▇█"

SCALE_H = " ▏▎▍▌▋▊▉█"
OTHER={
    "upper_half":"▀",
    "upper_eighth":"▔",
    "lower_eighth":"▁",
    "right_half":"▐",
    "left_eighth":"▏",
    "right_eighth":"▕",
    "light":"░",
    "medium":"▒",
    "dark":"▓",
    "misc":"▖▗▘▙▚▛▜▝▞▟◐◑◒◓◔◕"
}

BOX = " ╵╶└╷│┌├╴┘─┴┐┤┬┼"
_box = "⊢⊣⊤⊥"

marker = "╵"

leftright = "▗▖"
# leftright = "⌋⌊"
# leftright = "⌋⌊"
# leftright = "⌟⌞"
# leftright = ".."
# leftright = "❳❲"
# leftright = "◿◺"
# leftright = "◢◣"
# leftright= "◅▻"
# leftright= "◁▷"
# leftright = "⎦⎣"

arrows = {
        "default":[
            (0X21BD, '↽'), #          ↽-----          -----↽
            (0X21C0, '⇀'), #          ⇀-----          -----⇀
        ],
    
        "extra":[     
            (0X27F5, '⟵'), #          ⟵-----          -----⟵
            (0X27F6, '⟶'), #          ⟶-----          -----⟶
            (0X21FD, '⇽'), #          ⇽-----          -----⇽
            (0X21FE, '⇾'), #          ⇾-----          -----⇾
            (0X21E6, '⇦'), #          ⇦-----          -----⇦
            (0X21E8, '⇨'), #          ⇨-----          -----⇨
            (0X21E0, '⇠'), #          ⇠-----          -----⇠
            (0X21E2, '⇢'), #          ⇢-----          -----⇢
            (0X21D0, '⇐'), #          ⇐-----          -----⇐
            (0X21D2, '⇒'), #          ⇒-----          -----⇒
            (0X21C0, '⇀'), #          ⇀-----          -----⇀
            (0X21C1, '⇁'), #          ⇁-----          -----⇁
            (0X21C0, '⇀'), #          ⇀-----          -----⇀
            (0X21BD, '↽'), #          ↽-----          -----↽
            (0X21BC, '↼'), #          ↼-----          -----↼
            (0X2190, '←'), #          ←-----          -----←
            (0X2192, '→'), #          →-----          -----→
               
               
            
            (0X293A, '⤺'), #          ⤺-----          -----⤺
            (0X293A, '⤺'), #          ⤺-----          -----⤺
            (0X293B, '⤻'), #          ⤻-----          -----⤻
            (0X293C, '⤼'), #          ⤼-----          -----⤼
            (0X293D, '⤽'), #          ⤽-----          -----⤽
            (0X21B6, '↶'), #          ↶-----          -----↶
            (0X21B7, '↷'), #          ↷-----          -----↷
            (0X219C, '↜'), #          ↜-----          -----↜
            (0X219D, '↝'), #          ↝-----          -----↝
        ]
}

RESET = '\033[0m'

# https://www.calculators.org/math/html-arrows.php so many arrows.....
_arrows = "".join(chr(i) for i in range(0x2933, 0x2941 + 1))
_arrows2 = "".join(chr(i) for i in range(0x2962, 0x2965+1)) + str(chr(0x2970))
_arrows3 = "".join(chr(i) for i in range(0x2794, 0x27B2+1)) + str(chr(0x27BE)) # these symbols are cool in the right font
_arrows4 = "⤏⤎⤍⤌⟿⟾⟽⟹⟶⟵"

_symbols = [
    "♣","♥","⚆","⚇","⚈","⚉","⚔","⚝",
    "⚲","⚴","⛢","⚵","⚶","⚷","⚸","⚱","","","","","","","",
    "⚺","⚻","",
    "⚔","","","","","","","","","","","","","","","",
    "⍎","⍏","⍕","⍖","⍗","⍙","⍚","⍜","⍢","⍦","⍱","⍲","⎈","⎉","⎊","⎋",
    "⎌","⎍","⎏","⎐","⎑","⎒","⏀","⏁","⏂","⏃","⏄","⏅","⏚","","","",
    
    "⚙","⚛","☀","☼","☉","☸", # for proteins? idk
    "⌖","⌘","⌗","⌾","⎮","","","","","",
    "◐","◑","◒","◓","◸","◹","◺","◿","","","","","","","","",
    "◐","◑","◒","◓","◸","◹","◺","◿","","","","","","","","",
    
]

_curves = [
    "⎡","⎢","⎣","","","","","","","","","","","","","",
]

_sick_chars = [chr(i) for i in range(0x1400, 0x1680)] + [chr(i) for i in range(0x18B0, 0x18F6)]


_protein_metavocab = "".join(["☉","☼","⚙","⚛","❁","✾"])
# _protein_metavocab = "".join(["◉","◎","◴","◵","◶","◷"])

def get_arrow(name, stem_len, stem_char = '-'):
    """"
    returns L, R
    """
    left_d, right_d = arrows.get(name, arrows.get("default"))
    left_h, right_h = left_d[1], right_d[1]
    if stem_len < 1:
        return left_h, right_h
    left_f = left_h + stem_len * stem_char
    right_f = stem_len * stem_char + right_h
    
    return left_f, right_f
    
color_names = ["black","red","green","yellow","blue","magenta","cyan","white"]
color_ints = range(len(color_names))

class Colors:
    
    standard = range(8)
    high_intensity = range(8, 16)
    colors = range(16, 232)
    grays = range(232,256)
    
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'

    # Variant colors
    SNP = '\033[91m'        # Red for SNPs
    INSERTION = '\033[92m'   # Green for insertions
    DELETION = '\033[93m'    # Yellow for deletions

    # Feature colors
    GENE = '\033[94m'        # Blue
    TRANSCRIPT = '\033[95m'  # Magenta
    EXON = '\033[96m'        # Cyan
    CDS = '\033[93m'         # Yellow
    UTR = '\033[90m'         # Gray
    REPEAT = "\x1b[38;5;143m"
    
    START_CODON = '\x1b[35m'
    STOP_CODON = '\x1b[35m'

    # Motif colors (for underlines)
    MOTIF = '\033[96m'   # Cyan for other motifs
    HIGHLIGHT = '\x1b[38;5;148m' # goldish
    
    # Navigation
    POSITION = '\033[97m'    # White
    
    def __init__(self):
        
        self.text = self.get_color(45, background = False)
        self.background = self.get_color(234, background = True)
        self.effect = self.get_effect(0)
    
    def set(self, spec, bright = False, background = False, effect = False):
        if background:
            self.set_background(spec, bright=bright)
        else:
            self.set_text(spec, bright=bright)
        
        if effect:
            self.set_effect(spec)
    
    def as_tuple(self):
        return (self.text, self.background, self.effect)
    
    def set_text(self, text_spec, bright=False):
        self.text = self.get_color(text_spec, bright=bright, background = False)
        
    def set_background(self, background_spec, bright=False):
        self.text = self.get_color(background_spec, bright=bright, background = True)
        
    def set_effect(self, effect_spec):
        self.effect = self.get_effect(effect_spec)
        
    def get_color(self, color_spec, bright=False, background = False):
        
        cc = 0
        if isinstance(color_spec, str) and color_spec:
            cc = color_ints[color_names.index(color_spec)]
        elif isinstance(color_spec, int):
            cc = color_spec
        
        bgc = "38;5;"
        if background:
            bgc = "48;5;"
        
        # return f'\x1b[{bgc}{cc}m'
        return '\x1b[' + str(bgc) + str(cc) + 'm'

    def _get_color_scale(self, start, end, num_values, dx, minval, maxval):
        
        d1 = dx*dx
        d2 = dx
        if isinstance(start, int):
            start -= minval
            end -= minval
            start_vec = [start // d2, (start - (start // dx))%d2, start % d1] # all dims on [0, 6)
            end_vec = [end // d1, (end - (end // dx))%d2, end % d1]
        else:
            start_vec = [s - minval for s in start]
            end_vec = [e - minval for e in end]

        delta = [(e-s)/num_values for s,e in zip(start_vec, end_vec)]
        
        col_vecs = [[int(round(s+n*d)) for s, d in zip(start_vec, delta)] for n in range(num_values+1)]
        cols = [min(maxval, max(minval,sum([d1*cv[0] + d2*cv[1] + cv[2]]) + minval)) for cv in col_vecs]
        cols = [cols[i] for i in range(len(cols)) if i==0 or cols[i-1]!=cols[i]]
        return cols

    def get_color_scale(self, start, end, num_values):
        return self._get_color_scale(start, end, num_values, 6, self.colors[0], self.colors[-1])

    def get_color_scale_16b(self, startrgb, endrgb, num_values):
        d = 40 #?
        mincol = (2**16 - d**3)//2
        maxcol = 2**16 - mincol
        return self._get_color_scale(startrgb, endrgb, num_values, d, mincol,maxcol)
    
    def get_color_scale_24b(self, startrgb, endrgb, num_values):
        
        color_scale = []
        
        for n in range(num_values):
            cn = []
            for cs, ce in zip(startrgb, endrgb):
                col = min(255, max(0, round(n * (ce - cs) / num_values + cs)))
                cn.append(col)
            color_scale.append(cn)
        
        return color_scale

    def get_effect(self, effect_spec):
        if effect_spec == "bold":
            return "\x1b[1m"
        elif effect_spec == "dim":
            return "\x1b[2m"
        elif effect_spec == "underline":
            return "\x1b[4m"
        elif effect_spec == "blink":
            return "\x1b[5m"
        elif effect_spec == "reverse":
            return "\x1b[7m"
        return ""
    
    def __str__(self):
        return str(self.effect) + str(self.background) + str(self.text)
    
    @classmethod
    def from_specs(cls, text_spec = "", text_bright = False, bg_spec = "", bg_bright = False, effect_spec = ""):
        out = cls()
        out.set_text(text_spec, bright = text_bright)
        out.set_background(bg_spec, bright=bg_bright)
        out.set_effect(effect_spec)
        return out

RESET = '\033[0m'

CS = Colors.from_specs(text_spec=250, text_bright = True, effect_spec ="")
CD = Colors.from_specs(text_spec="yellow", effect_spec ="")
CL = Colors.from_specs(text_spec="cyan",effect_spec ="")
CB = Colors.from_specs(text_spec="blue",effect_spec ="")
CC = Colors.from_specs(text_spec="cyan",effect_spec ="")

def set_colors(tail=None, dyad=None, loop=None, seq=None, cseq=None, bright=False, background = False, effect = None):
    
    global CS, CD, CL, CB, CC
    
    if tail:
        CS.set(tail, bright=bright, background = background, effect = effect)
    if dyad:
        CD.set(dyad, bright=bright, background = background, effect = effect)
    if loop:
        CL.set(loop, bright=bright, background = background, effect = effect)
    if seq:
        CB.set(seq, bright=bright, background = background, effect = effect)
    if cseq:
        CC.set(cseq, bright=bright, background = background, effect = effect)

set_colors(seq = 174, cseq = 66, background = True)

def get_color_scheme(name):
    """
    returns bg, fg
    """
    if name == "gray":
        return 244, 236
    elif name == "blue":
        return 17, 38
    elif name == "foggy":
        return 36, 67
    elif name == "dusty":
        return 188, 138
    elif name == "ruddy":
        return 179, 131
    elif name == "icy":
        return 146, 225
    elif name == "vscode":
        return 234, 131
    elif name == "test":
        return 234, 65
    else:
        return 0,1


def get_color_scheme_24b(name):
    
    if name == "gray":
        return [63,36,97], [255,92,131]
    elif name == "coolwarm":
        return [26,158,229], [250,144,50]
    elif name == "sweet":
        return [63,36,97], [255,92,131]
    elif name == "lava":
        return [28,55,57], [196,55,57] 
    elif name == "energy":
        return [36,71,122], [245,178,37]
    elif name == "deep":
        return [20,34,78], [180,34,78]
    elif name == "terra":
        return [19,118,83], [244,143,35]
    elif name == "unterra":
        return [19,118,83], [125,25,125]
    elif name == "vscode":
        return [28,28,28], [28,28,28]

def get_ansi_color(col, bg = False):
    
    if bg:
        opt = 48
    else:
        opt = 38
    
    if isinstance(col, list):
        r,g,b = col
        code = f"\x1b[{opt};2;{r};{g};{b}m"
    else:
        code = f"\x1b[{opt};5;{col}m"
    
    return code

def get_fgbg(fg_color, bg_color):
    fg = f"\x1b[38;5;{fg_color}m"
    bg = f"\x1b[48;5;{bg_color}m"
    return fg, bg


def scalar_to_text_8b(scalars, minval = None, maxval = None, fg_color = 53, bg_color = 234, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 8, flip = flip)

def scalar_to_text_16b(scalars, minval = None, maxval = None, fg_color = 53, bg_color = 234, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 16, flip = flip)

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
        bg, fg = get_fgbg(bg_color, fg_color)
    else:
        fg, bg = get_fgbg(fg_color, bg_color)
    
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
        # return [fg+bg+eff + SCALE[-1] + RESET]

    
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
        row.append(RESET)
        outstrs.append("".join(row))
    
    if add_border:
        outstrs.insert(0, " " + SCALE[1]*ncols + " ")
        outstrs.append(f"{brdc} " + OTHER.get("upper_eighth","")*ncols + f" {RESET}")
    
    if add_range:
        # hilo = "⎴⎵"
        # hilo = "⏋⏌"
        ran_fstr = kwargs.get("range_fstr", "0.2f")
        hilo = "⌝⌟"
        hi, lo = list(hilo)
        hi, lo = (lo, hi) if flip else (hi, lo)
        minstr = format(minval, ran_fstr)
        maxstr = format(maxval, ran_fstr)
        outstrs[0] += hi + maxstr
        if bit_depth > 8:
            outstrs[-1] += lo + minstr
        
    if flip:
        outstrs = flip_scalar_text(outstrs)
    
    return outstrs

def scalar_to_text_mid(scalars, center = None, rng = None, fg_color = 53, bg_color = 234,  effect = None):
    
    bit_depth = 16
    
    ifg, ibg = get_fgbg(fg_color, bg_color)
    bg, fg = get_fgbg(bg_color, fg_color)
    
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
        row.append(RESET)
        outstrs.append("".join(row))
        
    return outstrs

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
            fmtr = utils.format_genomic
            # nexp = np.log10(max(abs(xmax), abs(xmin), 1))
            # if nexp > 6:
            #     div = 1e6
            #     unit = "M"
            # elif nexp > 3:
            #     div = 1e3
            #     unit = "k"
            # else:
            #     div = 1
            #     unit = "b"
            # fmtr = lambda x: format(x/div, ".1f") + unit
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
            formatter = utils.format_genomic
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

def heatmap(data, row_labels=None, col_labels=None, minval=None, maxval=None, center=0, color_scheme = "terra",
            show_values=False, value_fmt="0.2f", value_labels = None, colorbar=True, suppress = False, **kwargs):

    out_rows = []

    if isinstance(data, np.ndarray):
        data = data.tolist()

    min_color, max_color = get_color_scheme_24b(color_scheme)
    mid_color, _ = get_color_scheme_24b("vscode")
    
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
            fg = get_ansi_color(color_code)
            colorbar_line += fg + SCALE[-1] + RESET

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
    bg = get_ansi_color(234, bg = False)
    
    if col_space > 0 and col_space < 1:
        col_width -= 1
        print(col_space, col_width)
    
    for j, val in enumerate(data_row):
        if val is None:
            line += " "*(col_width + col_space)
            continue
        
        color_code = value_to_color(val)
        fg = get_ansi_color(color_code, bg = True)

        if value_labels:
            val_str = value_labels[j]
            col_str = val_str[:2].center(col_width)
        else:
            col_str = block * col_width
        
        # if col_space > 0 and col_space < 1:
        #     space_str = col_border
        # else:
        space_str = space * col_space
        
        line += bg + fg + col_str + space_str + RESET
        
    return line
    
def scrub_ansi(line):
    
    import re
    ansi_re = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    newline = ansi_re.sub("", line)
    return newline

def visible_len(line):
    """Get the length of a string excluding ANSI escape sequences."""
    return len(scrub_ansi(line))

def visible_slice(line, start=0, stop=None, step=1):
    """
    Slice a string based on visible character positions, preserving ANSI codes.

    This extracts visible characters from position start to stop (exclusive),
    while keeping any ANSI escape sequences that apply to those characters.
    """
    import re
    ansi_re = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")

    # Parse the string into tokens: either ANSI codes or visible characters
    tokens = []
    pos = 0
    for match in ansi_re.finditer(line):
        # Add any visible characters before this ANSI code
        if match.start() > pos:
            for char in line[pos:match.start()]:
                tokens.append(('char', char))
        # Add the ANSI code
        tokens.append(('ansi', match.group()))
        pos = match.end()
    # Add any remaining visible characters
    if pos < len(line):
        for char in line[pos:]:
            tokens.append(('char', char))

    # Count visible characters and determine slice bounds
    vlen = sum(1 for t in tokens if t[0] == 'char')
    if stop is None:
        stop = vlen

    # Extract the slice, keeping ANSI codes
    result = []
    visible_idx = 0
    last_ansi = None  # Track the last ANSI code we've seen

    for token_type, token_val in tokens:
        if token_type == 'ansi':
            # Keep track of ANSI codes - we might need them
            last_ansi = token_val
            # If we're in the slice range, include this ANSI code
            if start <= visible_idx < stop:
                result.append(token_val)
        else:  # token_type == 'char'
            if visible_idx == start and last_ansi and not result:
                # Starting a slice - include the last ANSI code to maintain color
                result.append(last_ansi)

            if start <= visible_idx < stop and (visible_idx - start) % step == 0:
                result.append(token_val)

            visible_idx += 1

            if visible_idx >= stop:
                break

    # Add RESET at the end if we have any content
    if result and RESET not in ''.join(result[-5:]):
        result.append(RESET)

    return ''.join(result)

def plot_adjacent(seqa, seqb):
    bit_depth = 16
    qseqa = quantize(seqa, bit_depth, mid=True)
    qseqb = quantize(seqb, bit_depth, mid=True)
    
    min_offset = 0
    for a, b in zip(qseqa, qseqb):
        off = a-b
        if off < min_offset:
            min_offset = off
    min_offset -= 1
    eff_bit_depth = bit_depth - min_offset
    
    
    return qseqa, [b+min_offset for b in qseqb], min_offset, eff_bit_depth

def quantize(data, bit_depth, maxval=None, minval=None, mid = False):
    if not maxval:
        maxval = max(data)
    if not minval:
        minval = min(data)
    rng = maxval-minval
    c = (minval+maxval)/2
    off = 0 if mid else 0.5
    return [int(bit_depth * (((d-c)/rng) + off)) for d in data]

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
    rng = (maxval - minval)/1
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
        row.append(RESET)
        border_row = ["▁" if r in SCALE else r for r in row]
        outstrs.append("".join(row) + "\033[1G")
        outstrs.append("".join(border_row))
        print(outstrs[-1])
        
    if flip:
        return hflip_scalar_text(outstrs)
    else:
        return outstrs

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

def clean_scalar_text(sctext):
    if isinstance(sctext, str):
        sctext = [sctext]
    elif not isinstance(sctext, list):
        sctext = list(sctext)
    
    outrows = []
    for row in sctext:
        newrow = [t for t in row if t in SCALE]
        outrows.append("".join(newrow))
    return outrows

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

def highlight_subsequences(subseqs, colors, delimit = ""):
    
    if isinstance(colors[0], int):
        colors = [f"\x1b[38;5;{c}m" for c in colors]
    
    lastcolor = colors[0]
    out = [colors[0]]
    for ss, c in zip(subseqs, colors):
        if c == lastcolor:
            pass
        else:
            out.append(RESET)
            if delimit:
                out.append(delimit)
            out.append(c)
            lastcolor = c
        out.append(ss)
    out.append(RESET)
    return "".join(out)

def make_start_ends(feature_name, feature_positions, feature_length, starts = {}, ends = {}):
    
    for p in feature_positions:
        if not p in starts:
            starts[p] = []
        starts[p].append(feature_name)
        
        end_pos = p + feature_length
        if not end_pos in ends:
            ends[end_pos] = []
        ends[end_pos].append(feature_name)
    
    return starts, ends

def make_spans(starts, ends):
    
    feats = set()
    istarts = {}
    iends = {}
    spans = {}
    
    for s, fs in starts.items():
        
        for f in fs:
            if not f in istarts:
                feats.add(f)
                istarts[f] = []
            istarts[f].append(s)
    
    for e, fs in ends.items():
        for f in fs:
            if not f in iends:
                feats.add(f)
                iends[f] = []
            iends[f].append(e)
    
    for f in feats:
        st = istarts[f]
        en = iends[f]
        spans[f] = [(s, e) for s, e in zip(sorted(st), sorted(en))]
    
    return spans

def make_key(features, colors):
    
    parts = ["key:"]
    for f in features:
        cf = colors.get(f)
        cstr = f"\x1b[38;5;{cf}m"
        fstr = "".join([cstr, f, RESET])
        parts.append(fstr)
    
    return " ".join(parts)

def highlight_matching(seqa, seqb, colors = None, do_rc = False, do_both = False, suppress = False, color_bg = False, chunksz = 256):
    
    if do_both:
        do_rc = False
    
    bg_frm = "\x1b[48;5;{}m"
    fg_frm = "\x1b[38;5;{}m"
    
    if color_bg:
        cfrm = bg_frm
    else:
        cfrm = fg_frm
    
    base_color_fg = fg_frm.format(248)
    base_color_bg = cfrm.format(234)
    base_color = Colors.RESET + base_color_fg + base_color_bg
    if not colors:
        color = cfrm.format(142)
        rcolor = Colors.MOTIF
    else:
        c, cr = colors
        color = cfrm.format(c)
        rcolor = cfrm.format(cr)
    
    if color_bg:
        color = color + Colors.BOLD
        rcolor = rcolor + Colors.BOLD
    
    seqah = [base_color]
    seqbh = [base_color]
    
    rseqb = reversed(seqb)
    
    n = 0
    for sa, sb, rsb in zip(seqa, seqb, rseqb):
        
        rcsb = complement(rsb)
        
        if sa == sb and not do_rc:
            pre = color
            post = base_color
        elif sa == rcsb and do_rc:
            pre = color
            post = base_color
        elif sa == rcsb and do_both:
            pre = rcolor
            post = base_color
        else:
            pre = ""
            post = ""
        
        seqah.append(f"{pre}{sa}{post}")
        seqbh.append(f"{pre}{sb}{post}")
        n+=1
    
    seqah.append(Colors.RESET)
    seqbh.append(Colors.RESET)
    
    seqaout = "".join(seqah)
    seqbout = "".join(seqbh)
    
    if not suppress:
        for n in range(len(seqa)//chunksz + 1):
            print(seqaout[n*chunksz:(n+1)*chunksz])
            print(seqbout[n*chunksz:(n+1)*chunksz])
    print(Colors.RESET)
    return "".join(seqah), "".join(seqbh)

def highlight_correlated(full_seq, shift, colors = None, suppress = False):
    
    if not colors:
        ca = Colors.RCMOTIF
        cab= Colors.HIGHLIGHT
        cb = Colors.MOTIF
    
    seq_len = len(full_seq)
    
    seqah = []
    seqbh = []
    
    sas = min(0, shift)
    sbs = max(0, -shift)
    
    for n in range(len(full_seq)):
        
        sf = full_seq[n]
        sa = sb = ""
        
        if sas + n < seq_len:
            sa = full_seq[sas + n]        
        if sbs + n < seq_len:
            sb = full_seq[sbs + n]
        
        if sa == sb:
            pre = ca
            post = Colors.RESET
        else:
            pre = ""
            post = ""
        
        if sa and sb:
            seqah.append(f"{cab}{sf}{post}")
        elif sa:
            seqah.append(f"{ca}{sf}{post}")
        elif sb:
            seqah.append(f"{cb}{sf}{post}")
        
    if not suppress:
        print("".join(seqah))
        print("".join(seqbh))
    return "".join(seqah+seqbh)

def highlight_features(seq, features, feature_spans = {}, feature_starts = {}, feature_ends = {}, colors = {}, show_key = True, suppress = True, break_features = False):
    """
    feature_starts: Dict:feature -> List[feature_start_pos]
    feature_ends: Dict:feature -> List[feature_end_pos]
    """
    # Build the colored string
    
    if not feature_spans:
        # feature_spans = {f:[(s, e) for s, e in zip(feature_starts[f], feature_ends[f])]for f in features}
        feature_spans = make_spans(feature_starts, feature_ends)
        # print(feature_spans)
    
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    result = []
    current_color = 0
    last_ansi = bc
    
    for s in features:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    key = make_key(features, colors)
    
    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        
        current_color = baseline_color
        
        pre_break = False
        done = False
        for f in features:
            for s, e in feature_spans[f]:
                if i < s:
                    continue
                elif i >= e:
                    continue
                else:
                    current_color = colors.get(f)
                    pre_break = i == s 
                    done = True
                if done:
                    break
            if done:
                break

        ansi = f"\x1b[38;5;{current_color}m"

        if ansi != last_ansi:
            if break_features:
                result.append(" ")
            result.append(ansi)
            last_ansi = ansi
        elif break_features and pre_break:
            result.append(" ")
        
        result.append(seq[i])        

    result.append(RESET)
    if not suppress:
        print("".join(result))
        if show_key:
            print(key)
    return result, key

def highlight_sequence(seq, subseq, color = None, suppress = True, show_key = True):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    if not color:
        color = random.randint(20, 230)
    colors = {subseq:color}

    start_pos = find_subsequence(seq, subseq)
    starts, ends = make_start_ends(subseq, start_pos, len(subseq), starts=starts, ends = ends)

    return highlight_features(seq, [subseq], feature_starts = starts, feature_ends = ends, colors=colors, suppress = suppress, show_key = show_key)

def highlight_sequence_fuzzy(seq, subseq, colors = None, max_err = 1, suppress = True, show_key = True, break_features = False):
    
    if not colors:
        colors = []
    if len(colors) <= max_err:
        colors.extend([random.randint(20, 230) for i in range(len(colors), max_err+1)])
    
    feats = [subseq]
    spans = {subseq:[]}
    clrs = {subseq:colors[0]}
    
    key_info = {}
    
    ptrn_str = "(%s){e<=%s}" % (subseq, str(max_err))
    ptrn = regex.compile(ptrn_str, regex.BESTMATCH)
    
    matches = regex.finditer(ptrn, seq)
    
    for m in matches:
        if not m:
            continue
        
        fsubseq = m.groups()[0]
        
        err = sum(m.fuzzy_counts)
        start, end = m.spans()[0]
        
        if not fsubseq in feats:
            feats.append(fsubseq)
            spans[fsubseq] = []
            clrs[fsubseq] = colors[err]
        if not err in key_info:
            key_info[err] = set()
            
        key_info[err].add(fsubseq)
        spans[fsubseq].append((start, end))
    
    out, _ = highlight_features(seq, feats, feature_spans = spans, colors = clrs, show_key = False, suppress = suppress, break_features = break_features)
    key = ", ".join([f"\x1b[38;5;{colors[err]}merr={err}:[{",".join(subseqs)}]" for err, subseqs in key_info.items()]) + "\x1b[0m"
    if show_key and not suppress:
        print(key)
    
    return out, key

def highlight_sequences(seq:str, subseqs:List[str], start_pos = None, do_rc = False, min_len = 5, colors={}, show_key = True, suppress = True):

    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    for s in subseqs:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    if not start_pos:
        start_pos = find_subsequences(seq, subseqs, do_rc = do_rc)
    for s in subseqs:
        if len(s) < min_len:
            continue

        starts, ends = make_start_ends(s, start_pos[s], len(s), starts=starts, ends = ends)

    highlight_features(seq, subseqs, feature_starts = starts, feature_ends = ends, colors = colors, show_key = show_key, suppress = suppress)

def _highlight_runs_auto(seqa, top_seqs, top_datas, **kwargs):
    all_seqs = []
    start_pos = {}
    
    for ns in range(len(top_seqs)):
        run, ind, shift = top_datas[ns]
        start_pos.update({top_seqs[ns][0]: [ind], top_seqs[ns][1]: [ind - shift]})
        all_seqs.extend(top_seqs[ns])
        
    highlight_sequences(seqa, all_seqs, start_pos = start_pos, **kwargs)

def _highlight_runs_diff(seqa, seqb, top_seqs, top_datas, **kwargs):
    
    seqs_a = []
    start_pos_a = {}
    seqs_b = []  
    start_pos_b = {}
    for ns in range(len(top_seqs)):
        run, ind, shift = top_datas[ns]
        start_pos_a.update({top_seqs[ns][0]: [ind]})
        start_pos_b.update({top_seqs[ns][1]: [ind - shift]})
        seqs_a.append(top_seqs[ns][0])
        seqs_b.append(top_seqs[ns][1])
        
    highlight_sequences(seqa, seqs_a, start_pos = start_pos_a, **kwargs)
    highlight_sequences(seqb, seqs_b, start_pos = start_pos_b, **kwargs)
    

def highlight_run(seqa, seqb, top_seqs, top_datas, suppress = False, show_key = True, **kwargs):
    
    colors = kwargs.pop("colors", {})
    
    for s in top_seqs:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    kwargs["colors"] = colors
    
    if seqa == seqb:
        _highlight_runs_auto(seqa, top_seqs, top_datas, suppress = suppress, show_key = show_key, **kwargs)
    else:
        _highlight_runs_diff(seqa, seqb, top_seqs, top_datas, suppress = suppress, show_key = show_key, **kwargs)

def highlight_sequences_in_frame(seq:str, subseqs:List[str], frame_start, min_len = 5):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there

    # Baseline color - a visible gray (color 240 is a nice medium gray)
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    colors = {}
    for s in subseqs:
        colors[s] = random.randint(20, 230)

    # Find all occurrences of each subsequence
    for s in subseqs:
        if len(s) < min_len:
            continue
        
        seq_pos = find_subsequence(seq, s, frame_start = frame_start)
        for p in seq_pos:
            if not p in starts:
                starts[p] = []
            starts[p].append(s)
        
            end_pos = p + len(s)
            if not end_pos in ends:
                ends[end_pos] = []
            ends[end_pos].append(s)

    # Build the colored string
    result = []
    active_seqs = []  # Currently active sequences
    current_color = 0
    last_ansi = bc

    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        if i in ends:
            for s in ends[i]:
                if s in active_seqs:
                    active_seqs.remove(s)

        if i in starts:
            for s in starts[i]:
                if s not in active_seqs:
                    active_seqs.append(s)

        if not active_seqs:
            current_color = baseline_color
        elif len(active_seqs) == 1:
            current_color = colors[active_seqs[0]]
        elif len(active_seqs) == 2:
            color_sum = sum(colors[s] for s in active_seqs)
            current_color = 20 + (color_sum % 211)
        else:
            current_color = 255

        ansi = f"\x1b[38;5;{current_color}m"

        if ansi != last_ansi:
            result.append(ansi)
            last_ansi = ansi

        result.append(seq[i])

    result.append(RESET)

    print("".join(result))

def highlight_sequence_by_span(seq, span_colors = {}, default_color = '\033[97m'):
    # span is (start, stop):color
    
    spans = sorted(span_colors.keys(), key=lambda k:k[0])
    ccurr = default_color
    colored_seq = []
    for i, b in enumerate(seq):
        
        in_span = False
        for st, sp in spans:
            c = span_colors[(st, sp)]
            if i>sp:
                continue
            elif i<st:
                break
            
            if i>=st and i<sp:
                in_span = True
                new_c = c
                if new_c == ccurr:
                    continue
                else:
                    colored_seq.append(new_c)
                    ccurr = new_c
        
        if not in_span:
            colored_seq.append(default_color)
        
        colored_seq.append(b)
    colored_seq.append(RESET)
    return "".join(colored_seq)

def highlight_dyads(seq, dyads):
    
    sec = {(d.stem_start, d.end_position): random.randint(20, 230) for d in sorted(dyads, key = lambda d:d.stem_start)}
    
    outseq = []
    
    for i in range(len(seq)):
        
        pre = ""
        post = ""
        for st, en in sec:
            if i < st:
                break
            if i > en:
                continue
            c = sec[(st, en)]
            pre = f"\x1b[38;5;{c}m"
            post = Colors.RESET
            
        outseq.append(f"{pre}{seq[i]}{post}")
    
    return "".join(outseq)
    
def draw_gene_structure(gene_features, gene_name, strand):
    """Draw a simple gene structure cartoon"""
    fig, ax = plt.subplots(figsize=(12, 3))
    
    # Sort features by position
    features = sorted(gene_features, key=lambda x: x['start'])
    
    # Gene baseline
    gene_start = features[0]['start']
    gene_end = features[-1]['end']
    ax.plot([gene_start, gene_end], [0.5, 0.5], 'k-', linewidth=1)
    
    # Draw features
    for feat in features:
        if feat['type'] == 'exon':
            # Determine if CDS or UTR
            height = 0.3 if 'CDS' in feat.get('extra_types', []) else 0.15
            color = 'darkblue' if 'CDS' in feat.get('extra_types', []) else 'lightblue'
            
            rect = patches.Rectangle(
                (feat['start'], 0.5 - height/2), 
                feat['end'] - feat['start'], 
                height,
                facecolor=color,
                edgecolor='black'
            )
            ax.add_patch(rect)
    
    # Add arrow for strand
    arrow_y = 0.5
    if strand == '+':
        ax.arrow(gene_end, arrow_y, 1000, 0, head_width=0.05, head_length=500, fc='red')
    else:
        ax.arrow(gene_start, arrow_y, -1000, 0, head_width=0.05, head_length=500, fc='red')
    
    # Formatting
    ax.set_ylim(0, 1)
    ax.set_xlim(gene_start - 5000, gene_end + 5000)
    ax.set_title(f"{gene_name} ({strand} strand)")
    ax.set_xlabel("Genomic Position")
    ax.get_yaxis().set_visible(False)
    
    plt.tight_layout()
    return fig

