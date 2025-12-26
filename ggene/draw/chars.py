
SCALE = " ▁▂▃▄▅▆▇█"     # ⎺ ⎽ ⎯
SCALE_H = " ▏▎▍▌▋▊▉█"   # ⎸⎹

SCALAR_PLOT = {
    "range_hi":"⌝", # hilo = "⎴⎵", "⏋⏌"
    "range_lo":"⌟",
}

HLINES = {
    # "head" : "◀▶",       # Filled triangles (recommended)
    # "head": "⊲⊳",
    # "head": "⋳⋻",
    "head": "⋖⋗", # math
    # "head": "≺≻", # math
    # "head": "ᗕᗒ", # indig
    # "head": "ᗏᗌ", # indig
    # "ᓬᕒ"
    # "tail":"ᗧᗤ", # indig, almost too good
    # "tail":"⎨⎬",
    "tail":"⚟⚞",
    # "tail" : "┤├",
    
    # "body":"−−",
    "body":"――",    # horizontal bar, punctuation
    # "body" : "──",    # box
    # "body": "⋯⋯",
    "body_dash" : "╌╌",
    "body_double" : "══",
    
    "arrow": "⟵⟶",
    "tiny": "◂▸",
    "small": "ᗏᗌ",
    # "cont": "⋯⋯",
    "cont": "≪≫",
    
    # Alternative arrowheads for different styles:
    "head_chevron": "<>",      # Simple chevrons (original)
    "head_triangle": "◁▷",     # Open triangles
    "head_filled": "◀▶",       # Filled triangles
    "head_harpoon": "↽⇀",      # Half-arrow harpoons (elegant)
    "head_arrow": "←→",        # Standard arrows
    "head_heavy": "⊲⊳",        # Heavy pointing triangles
    "head_tip":"╾╼",
    
    # Vertical connectors for label tracks
    "vline": "│",
    "vline_down": "╷",         # Top connector
    "vline_up": "╵",           # Bottom connector
    "corner_tl": "╭",          # Top-left corner
    "corner_tr": "╮",          # Top-right corner
    "corner_bl": "╰",          # Bottom-left corner
    "corner_br": "╯",          # Bottom-right corner
}

MAP = {
    "pointer_top":"▼",
    "fella":"ᘝ"
    
}

RULER={
    "tick":"╵",
    "minor_tick":"'",
    "llabel_marker":"╰",
    "rlabel_marker":"╯"
}

FSYMS = {
    "gene":"◯",
    "transcript":"◉",
    "exon":"◎",
    "cds":"☉",
    
    "repeat":"▪",
    "motif":"▱",
    
    
}

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

CHARS = {
    "star1":"⋆",
    "star2":"⊹",
    "star3":"∗",
    "star4":"★",
    "star5":"⚹",
    
    "inf":"∞",
    
    "hex":"⌬",
    
    "fella":"ᘝ",
}

"""
https://cloford.com/resources/charcodes/utf-8_geometric.htm
https://cloford.com/resources/charcodes/utf-8_misc-symbols.htm
https://cloford.com/resources/charcodes/utf-8_mathematical.htm
https://cloford.com/resources/charcodes/utf-8_technical.htm

https://cloford.com/resources/charcodes/utf-8_box-drawing.htm
https://cloford.com/resources/charcodes/utf-8_block-elements.htm

https://cloford.com/resources/charcodes/utf-8_arrows.htm
https://cloford.com/resources/charcodes/utf-8_punctuation.htm

"""

can_abo = {
    "ᐃ":{
        "ᐃ":"ᐃᐄ ᐎᐏᐐᐑᐬᐂ",
        "ᐁ":"ᐁ  ᐌᐍ",
        "ᐅ":"ᐅᐆᐇᐒᐓᐔᐕᐭ ᐖ",
        "ᐊ":"ᐊ ᐗᐘᐙᐚᐮ ᐛ",
    },
    "ᐱ":{
        "ᐱ":"ᐱ",
        "ᐯ":"ᐯ",
    },
    "ᑎ":{
        "ᑎ":"ᑎ",
        "ᑌ":"ᑌ",
    },
    "diacritics":"",
    "misc":"",
}


marker = "╵"

leftright = "▗▖" # "⌋⌊", "⌟⌞" ".." "❳❲" "◿◺" "◢◣" "◅▻" "◁▷" "⎦⎣"

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
    