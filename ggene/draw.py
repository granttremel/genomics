
from typing import List

import random

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from ggene.seqs.find import find_subsequence, find_subsequences

SCALE = " ▁▂▃▄▅▆▇█"

SCALE_H = " ▏▎▍▌▋▊▉█"
OTHER={
    "upper_half":"▀",
    "upper_eighth":"▔",
    "right_half":"▐",
    "right_eighth":"▕",
    "light":"░",
    "medium":"▒",
    "dark":"▓",
    "misc":"▖▗▘▙▚▛▜▝▞▟◐◑◒◓◔◕"
}

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

class Color:
    def __init__(self, text, background, effect):
        self.text = self.get_color(text, background = False)
        self.background = self.get_color(background, background = True)
        self.effect = self.get_effect(effect)
    
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
        
    @staticmethod
    def get_color(color_spec, bright=False, background = False):
        
        cc = 0
        if isinstance(color_spec, str) and color_spec:
            cc = color_ints[color_names.index(color_spec)]
        elif isinstance(color_spec, int):
            cc = color_spec
        
        # if bright:
        #     cc += 8
        bgc = "38;5;"
        if background:
            bgc = "48;5;"
        
        # return f'\x1b[{bgc}{cc}m'
        return '\x1b[' + str(bgc) + str(cc) + 'm'

    @staticmethod
    def get_effect(effect_spec):
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
        out = cls("","","")
        out.set_text(text_spec, bright = text_bright)
        out.set_background(bg_spec, bright=bg_bright)
        out.set_effect(effect_spec)
        return out

RESET = '\033[0m'

CS = Color.from_specs(text_spec=250, text_bright = True, effect_spec ="")
CD = Color.from_specs(text_spec="yellow", effect_spec ="")
CL = Color.from_specs(text_spec="cyan",effect_spec ="")
CB = Color.from_specs(text_spec="blue",effect_spec ="")
CC = Color.from_specs(text_spec="cyan",effect_spec ="")

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
def get_fgbg(fg_color, bg_color):
    fg = f"\x1b[38;5;{fg_color}m"
    bg = f"\x1b[48;5;{bg_color}m"
    return fg, bg

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
    


def scalar_to_text_8b(scalars, minval = None, maxval = None, fg_color = 212, bg_color = 7, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 8, flip = flip)

def scalar_to_text_16b(scalars, minval = None, maxval = None, fg_color = 212, bg_color = 7, flip = False):
    return scalar_to_text_nb(scalars, minval = minval, maxval = maxval, fg_color = fg_color, bg_color = bg_color, bit_depth = 16, flip = flip)

def scalar_to_text_nb(scalars, minval = None, maxval = None, fg_color = 7, bg_color = 212, bit_depth = 24, flip = False, effect = None):
    
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
        
    if flip:
        return flip_scalar_text(outstrs)
    else:
        return outstrs

def scalar_to_text_mid(scalars, center = None, rng = None, fg_color = 7, bg_color = 212,  effect = None):
    
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

def make_key(features, colors):
    
    parts = ["key:"]
    for f in features:
        cf = colors.get(f)
        cstr = f"\x1b[38;5;{cf}m"
        fstr = "".join([cstr, f, RESET])
        parts.append(fstr)
    
    return " ".join(parts)
    
def highlight_features(seq, features, feature_starts, feature_ends, colors = {}, show_key = True):
    """
    feature_starts: Dict:feature -> List[feature_start_pos]
    feature_ends: Dict:feature -> List[feature_end_pos]
    """
    # Build the colored string
    
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    result = []
    active_seqs = []  # Currently active sequences
    current_color = 0
    last_ansi = bc
    
    for s in features:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    key = make_key(features, colors)
    
    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        if i in feature_ends:
            for s in feature_ends[i]:
                if s in active_seqs:
                    active_seqs.remove(s)

        if i in feature_starts:
            for s in feature_starts[i]:
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
    if show_key:
        print(key)

def highlight_sequence(seq, subseq, colors = {}):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    if not subseq in colors:
        colors[subseq] = random.randint(20, 230)

    start_pos = find_subsequence(seq, subseq)
    starts, ends = make_start_ends(subseq, start_pos, len(subseq), starts=starts, ends = ends)

    highlight_features(seq, [subseq], starts, ends)
    pass

def highlight_sequences(seq:str, subseqs:List[str], start_pos = None, do_rc = False, min_len = 5, colors={}, show_key = True):

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

    highlight_features(seq, subseqs, starts, ends, colors = colors, show_key = show_key)

    
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

