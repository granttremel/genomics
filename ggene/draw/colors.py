
import re

from ggene.draw.chars import SCALE

class Colors:
    
    # RESET = '\033[38;5;224m\x1b[48;5;234m'
    RESET = '\x1b[0m\033[38;5;224m\x1b[48;5;234m'
    
    TEXT = "\x1b[38;5;224m"
    
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'
    HIGHLIGHT = '\x1b[38;5;148m' # goldish

    RC = '\x1b[38;5;125m'

    # Variant colors
    SNP = '\033[38;5;214m'        # Red for SNPs
    INSERTION = '\033[38;5;220m'   # Green for insertions
    DELETION = '\033[38;5;202m'    # Yellow for deletions

    # Feature colors
    GENE = '\033[94m'        # Blue
    TRANSCRIPT = '\033[95m'  # Magenta
    EXON = '\033[96m'        # Cyan
    CDS = '\033[93m'         # Yellow
    UTR = '\033[90m'         # Gray
    REPEAT = "\x1b[38;5;143m"
    PSEUDO = '\x1b[38;5;30m'
    
    START_CODON = '\x1b[35m'
    STOP_CODON = '\x1b[35m'

    # Motif colors (for underlines)
    MOTIF = '\033[96m'   # Cyan for other motifs
    
    # Motif colors (for underlines)
    MOTIF_SPL = '\033[38;5;111m'   # Cyan for other motifs
    MOTIF_PRO = '\033[38;5;112m'   # Cyan for other motifs
    RCMOTIF = '\x1b[38;5;106m'
    RCMOTIF_SPL = '\x1b[38;5;107m'
    RCMOTIF_PRO = '\x1b[38;5;108m'
    
    
    # Navigation
    POSITION = '\033[97m'    # White
    SUBTLE = '\x1b[38;5;240m'
    
    _color_frm = '\x1b[{c}m'
    _color_frm_8b = '\x1b[{b_or_f};5;{c}m'
    _color_frm_24b = '\x1b[{b_or_f};2;{c[0]};{c[1]};{c[2]}m'
    
    _fgi = 38
    _bgi = 48
    
    _24bi = 2
    _8bi = 5
    
    def __init__(self, fg = 8, bg = None, effect = None):
        
        self.fg_color = fg
        self.bg_color = bg
        if isinstance(effect, str):
            self.effect = self.get_effect_from_name(effect)
        else:
            self.effect = effect
    
    def set_color(self, fg = None, bg = None, effect = None):
        
        if fg:
            self.fg_color = fg
        if bg:
            self.bg_color = bg
        if effect:
            self.effect = effect
    
    def format_string(self, string):
        return f"{self.code}{string}{self.RESET}"
    
    @property
    def fg_code(self):
        return self.get_color(self.fg_color, background = False)
    
    @property
    def bg_code(self):
        return self.get_color(self.bg_color, background = True)
    
    @property
    def effect_code(self):
        return self.get_color(self.effect, background = False)
    
    @property
    def code(self):
        return self.effect_code + self.bg_code + self.fg_code
    
    @classmethod
    def colorize_string(cls, string, fg_color = 8, bg_color = 0, effect = 0):
        clr = cls(fg=fg_color, bg=bg_color, effect=effect)
        return clr.format_string(string)
    
    @classmethod
    def get_effect(cls, reset = False, bold = False, dim = False, underline = False):
        if reset:
            return cls.RESET
        elif bold:
            return cls.BOLD
        elif dim:
            return cls.DIM
        elif underline:
            return cls.UNDERLINE
        else:
            return cls.RESET    
    
    @classmethod
    def get_effect_from_name(cls, effect_name):
        
        if hasattr(cls, effect_name.upper()):
            return getattr(cls, effect_name.upper())
    
    @classmethod
    def get_color(cls, color_spec, background = False):
        
        if isinstance(color_spec, int):
            return cls._get_color_8b(color_spec, background = background)
        elif isinstance(color_spec, list) or isinstance(color_spec, tuple):
            if len(color_spec) == 3:
                return cls._get_color_24b(color_spec, background = background)
        elif color_spec is None:
            return ""
        else:
            return cls.RESET
    
    @classmethod
    def _get_color_8b(cls, color_spec, background = False):
        return cls._color_frm_8b.format(b_or_f = cls._bgi if background else cls._fgi ,c=color_spec)
    
    @classmethod
    def _get_color_24b(cls, color_spec, background = False):
        
        return cls._color_frm_24b.format(b_or_f = cls._bgi if background else cls._fgi, c=color_spec)
        
    @classmethod
    def get_colors_fgbg(cls, fg_color_spec, bg_color_spec):
        return cls.get_color(fg_color_spec, background = False), cls.get_color(bg_color_spec, background = True)
    
    @classmethod
    def scrub_codes(cls, line):
        
        ansi_re = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
        newline = ansi_re.sub("", line)
        return newline

    @classmethod
    def code_to_tuple(cls, code:str):
        
        code_re = re.compile(r"\x1b[([0-9]+;)?(;[0-9]+;)?([0-9]+)m")
        
        m = re.match(code_re, code)
        if m:
            tup = m.groups()
            outtup = []
            for v in tup:
                vv = ""
                if v:
                    vv = str(v)
                
                try:
                    vv = int(v)
                except:
                    vv = v
                
                outtup.append(vv)
            return tuple(outtup)
        else:
            return tuple()
    
    @classmethod
    def tuple_to_code(cls, tup):
        cstr = ";".join([str(v) for v in tup])
        return cls._color_frm.format(cstr)
    
    @classmethod
    def flip_fgbg(cls, code):
        
        a, fbgv, b = cls.code_to_tuple(code)
        
        new_fbgv = 0
        if fbgv == cls._fgi:
            new_fbgv = cls._bgi
        elif fbgv == cls._bgi:
            new_fbgv = cls._fgi
        
        return cls.tuple_to_code((a, new_fbgv, b))

    @classmethod
    def is_fg(cls, code):
        
        tup = cls.code_to_tuple(code)
        
        if tup[0] == cls._fgi:
            return True
        else:
            return False
    
    @classmethod
    def is_bg(cls, code):
        tup = cls.code_to_tuple(code)
        
        if tup[0] == cls._bgi:
            return True
        else:
            return False
    
    @classmethod
    def is_8b(cls, code):
        tup = cls.code_to_tuple(code)
        
        if len(tup) > 1 and tup[1] == cls._bgi:
            return True
        else:
            return False
    
    @classmethod
    def is_24b(cls, code):
        tup = cls.code_to_tuple(code)
        
        if len(tup) > 1 and tup[1] == cls._bgi:
            return True
        else:
            return False
    
    
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
        if effect_spec.lower() == "bold":
            return "\x1b[1m"
        elif effect_spec.lower() == "dim":
            return "\x1b[2m"
        elif effect_spec.lower() == "underline":
            return "\x1b[4m"
        elif effect_spec.lower() == "blink":
            return "\x1b[5m"
        elif effect_spec.lower() == "reverse":
            return "\x1b[7m"
        return ""
    
    def __str__(self):
        return self.code

    _color_schemes_8b = ["gray","blue","foggy","dusty","ruddy","icy","vscode","test"]

    @classmethod
    def get_color_scheme(cls, name):
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

    _color_schemes_24b = ["gray","coolwarm","sweet","lava","energy","deep","terra","unterra","vscode"]

    @classmethod
    def get_color_scheme_24b(cls, name):
        
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

    @classmethod
    def visible_len(cls, line):
        return visible_len(line)
    
    @classmethod
    def visible_slice(cls, line, start=0, stop=None, step=1):
        return visible_slice(line, start=start,stop=stop,step=step)
    
    # @classmethod
    # def show_colors(cls):
        
    #     for i in range(2**8):
            
    #         # fgcol = cls.get_color(i, background = False)
    #         # bgcol = cls.get_color(i, background = True)
    #         fgcol = "\x1b[38;5;{i}m".format(i=int(i))
    #         bgcol = "\x1b[48;5;{i}m".format(i=int(i))
            
    #         printstrs = [str(i), fgcol, "this is the color of fg", SCALE, cls.RESET, bgcol, SCALE[::-1], "this is the color of bg", cls.RESET]
    #         printstr = " ".join(printstrs)
    #         print(printstr)
    #         print(repr(printstr))


def visible_len(line):
    """Get the length of a string excluding ANSI escape sequences."""
    return len(Colors.scrub_codes(line))

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
    if result and Colors.RESET not in ''.join(result[-5:]):
        result.append(Colors.RESET)

    return ''.join(result)
