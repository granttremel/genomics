
import re

from ggene.draw.chars import SCALE

class Colors:
    
    # RESET = '\033[38;5;224m\x1b[48;5;234m'
    # RESET = '\x1b[0m\033[38;5;224m\x1b[48;5;234m'
    RESET = '\x1b[0m'
    
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
    MOTIF_TF = "\x1b[38;5;28m"
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
        
        code_re = re.compile(r"\x1b\[([0-9]+;)?([0-9]+;)?([0-9]+)m")
        
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

    @classmethod
    def get_color_scale_24b(cls, num_values, *rgbs):
        color_scale = []
        
        for i in range(1, len(rgbs)):
            startrgb = rgbs[i-1]
            endrgb = rgbs[i]
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

    _color_schemes_24b = {
        # Original schemes
        "gray": ([36, 36, 36], [220, 220, 220]),
        "coolwarm": ([26, 158, 229], [250, 144, 50]),
        "sweet": ([63, 36, 97], [255, 92, 131]),
        "peachy": ([251, 129, 124], [254, 55, 200]),
        "lava": ([28, 55, 57], [196, 55, 57]),
        "energy": ([36, 71, 122], [245, 178, 37]),
        "deep": ([20, 34, 78], [180, 34, 78]),
        "deepp": ([76, 18, 43], [14, 67, 124]),
        "bug": ([19, 118, 83], [244, 143, 35]),
        "terra": ([79,68,31],[30, 205, 40]),
        "terra_light": ([118, 55, 30], [49, 158, 83]),
        "unterra": ([19, 118, 83], [125, 25, 125]),
        "sunset": ([14, 131, 158], [234, 42, 38]),
        "aqua_vitae": ([11, 11, 160], [176, 24, 13]),
        "aqua_vitae_soft": ([33, 33, 182], [192, 58, 37]),
        "bluegreen": ([36, 17, 162], [63, 126, 7]),
        "sandy": ([79, 68, 31], [182, 171, 79]),
        "minty": ([68, 105, 46], [132, 227, 245]),
        "hott": ([121, 0, 69], [250, 11, 9]),
        "desert": ([135, 65, 38], [110, 135, 67]),
        "desert_dawn": ([119, 66, 105], [193, 130, 27]),
        "desert_dusk": ([73, 54, 43], [187, 88, 28]),
        "orangey": ([178, 91, 13], [255, 149, 25]),
        "dusky": ([30, 38, 74], [77, 53, 162]),
        "alpine": ([54, 244, 88], [108, 37, 164]),
        "vscode": ([28, 28, 28], [28, 28, 28]),

        # Ocean & Water themes
        "abyssal": ([8, 24, 58], [45, 149, 182]),
        "tideline": ([22, 78, 99], [167, 219, 216]),
        "bioluminescent": ([5, 15, 45], [72, 209, 204]),
        "coral_reef": ([255, 127, 80], [64, 224, 208]),
        "frozen_fjord": ([45, 55, 72], [176, 224, 230]),

        # Fire & Heat themes
        "ember_glow": ([45, 12, 8], [255, 140, 0]),
        "solar_flare": ([89, 13, 34], [255, 215, 0]),
        "magma_flow": ([25, 8, 8], [255, 69, 0]),
        "phoenix": ([75, 0, 30], [255, 191, 0]),
        "campfire": ([40, 20, 10], [255, 160, 50]),

        # Nature themes
        "moss": ([22, 38, 24], [144, 190, 109]),
        "autumn_canopy": ([62, 39, 35], [208, 134, 67]),
        "twilight_grove": ([28, 42, 58], [126, 188, 137]),
        "lichen": ([58, 63, 48], [168, 198, 134]),
        "fern": ([15, 42, 28], [119, 178, 85]),
        "geosmin": ([52, 61, 70], [134, 150, 167]),

        # Cosmic themes
        "nebula": ([20, 10, 38], [147, 88, 172]),
        "pulsar": ([10, 8, 28], [100, 149, 237]),
        "event_horizon": ([5, 5, 15], [88, 28, 135]),
        "aurora_borealis": ([18, 38, 58], [100, 255, 150]),
        "supernova": ([28, 18, 48], [255, 100, 120]),
        "dark_matter": ([12, 12, 18], [68, 68, 102]),

        # Warm Earth tones
        "terracotta": ([110, 47, 35], [215, 133, 95]),
        "saharan": ([92, 64, 38], [222, 184, 135]),
        "canyon_wall": ([75, 45, 35], [188, 143, 108]),
        "clay_kiln": ([58, 32, 28], [178, 102, 68]),
        "burnt_umber": ([48, 28, 18], [138, 75, 38]),

        # Cool & Muted themes
        "overcast": ([68, 72, 82], [148, 158, 172]),
        "slate_rain": ([48, 58, 68], [118, 138, 158]),
        "morning_fog": ([82, 88, 95], [178, 188, 198]),
        "glacier": ([58, 78, 98], [168, 208, 228]),
        "moonstone": ([48, 52, 68], [158, 168, 188]),

        # Vibrant & Bold themes
        "synthwave": ([28, 18, 58], [255, 20, 147]),
        "neon_jungle": ([15, 38, 28], [57, 255, 20]),
        "electric_dreams": ([18, 8, 48], [138, 43, 226]),
        "vaporwave": ([48, 28, 88], [255, 113, 206]),
        "cyberpunk": ([15, 12, 28], [0, 255, 255]),

        # Biological & Scientific themes
        "mitochondria": ([38, 58, 48], [158, 208, 88]),
        "chloroplast": ([28, 48, 28], [128, 208, 48]),
        "membrane": ([68, 48, 58], [218, 168, 188]),
        "nucleotide": ([28, 38, 68], [108, 158, 228]),
        "ribosome": ([58, 48, 38], [188, 158, 108]),

        # Gem & Mineral themes
        "amethyst": ([48, 28, 58], [153, 102, 204]),
        "malachite": ([18, 38, 28], [80, 200, 120]),
        "obsidian": ([15, 15, 18], [58, 58, 68]),
        "citrine": ([68, 48, 18], [228, 178, 48]),
        "sapphire": ([15, 28, 68], [65, 105, 225]),
        "opal": ([58, 58, 68], [188, 178, 198]),

        # Atmospheric themes
        "stratosphere": ([18, 28, 58], [88, 148, 208]),
        "smog": ([58, 55, 52], [138, 128, 118]),
        "golden_hour": ([68, 48, 38], [255, 183, 77]),
        "blue_hour": ([28, 38, 68], [98, 128, 188]),
        "thunderhead": ([38, 42, 52], [108, 118, 138]),
    }

    @classmethod
    def _adjust_brightness(cls, rgb, factor):
        """Adjust brightness of an RGB color. factor > 1 brightens, < 1 darkens."""
        return [min(255, max(0, int(c * factor))) for c in rgb]

    @classmethod
    def _adjust_contrast(cls, rgb, factor):
        """Adjust contrast relative to middle gray (128). factor > 1 increases contrast."""
        return [min(255, max(0, int(128 + (c - 128) * factor))) for c in rgb]

    @classmethod
    def _shift_hue(cls, rgb, degrees):
        """Shift hue by degrees (0-360). Converts RGB->HSV->RGB."""
        import colorsys
        r, g, b = [c / 255.0 for c in rgb]
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        h = (h + degrees / 360.0) % 1.0
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return [int(r * 255), int(g * 255), int(b * 255)]

    @classmethod
    def _adjust_saturation(cls, rgb, factor):
        """Adjust saturation. factor > 1 increases, < 1 decreases (towards gray)."""
        import colorsys
        r, g, b = [c / 255.0 for c in rgb]
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        s = min(1.0, max(0.0, s * factor))
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return [int(r * 255), int(g * 255), int(b * 255)]

    @classmethod
    def _adjust_color(cls, rgb, brightness=1.0, contrast=1.0, hue_shift=0, saturation=1.0):
        nrgb = rgb
        
        # order ?
        if hue_shift != 0:
            nrgb = cls._shift_hue(nrgb, hue_shift)
        if saturation != 1.0:
            nrgb = cls._adjust_saturation(nrgb, saturation)
        if contrast != 1.0:
            nrgb = cls._adjust_contrast(nrgb, contrast)
        if brightness != 1.0:
            nrgb = cls._adjust_brightness(nrgb, brightness)
        
        return nrgb

    @classmethod
    def add_middle(cls, startrgb, endrgb, num_mid = 1, brightness = 1.0, contrast = 1.0, hue_shift = 0, saturation = 1.0):
        
        mids = []
        n = num_mid + 1
        for nm in range(num_mid):
            cmid = [int((nm+1)*(c1 + c2)/n) for c1,c2 in zip(startrgb, endrgb)]
            
            cmid = cls._adjust_color(cmid, brightness=brightness, contrast = contrast, hue_shift = hue_shift, saturation=saturation)
            mids.append(cmid)
        
        return startrgb, *mids, endrgb

    @classmethod
    def get_color_scheme_24b(cls, name, brightness=1.0, contrast=1.0, hue_shift=0, saturation=1.0):
        """
        Get a 24-bit color scheme by name with optional adjustments.

        Args:
            name: Name of the color scheme
            brightness: Brightness multiplier (>1 brighter, <1 darker)
            contrast: Contrast multiplier (>1 more contrast, <1 less)
            hue_shift: Hue rotation in degrees (0-360)
            saturation: Saturation multiplier (>1 more saturated, <1 more muted)

        Returns:
            Tuple of (start_rgb, end_rgb) color lists
        """
        if name not in cls._color_schemes_24b:
            print(name)
            return None, None

        start, end = cls._color_schemes_24b[name]
        start, end = list(start), list(end)

        # Apply adjustments if any are non-default
        
        start = cls._adjust_color(start, brightness=brightness, contrast=contrast, hue_shift=hue_shift, saturation=saturation)
        end = cls._adjust_color(end, brightness=brightness, contrast=contrast, hue_shift=hue_shift, saturation=saturation)

        return start, end

    @classmethod
    def list_color_schemes_24b(cls):
        """Return list of available 24-bit color scheme names."""
        return list(cls._color_schemes_24b.keys())


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
