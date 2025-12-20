
SCALE = " ▁▂▃▄▅▆▇█"
SCALE_H = " ▏▎▍▌▋▊▉█"

SCALAR_PLOT = {
    "range_hi":"⌝", # hilo = "⎴⎵", "⏋⏌"
    "range_lo":"⌟",
}

HLINES = {
    "head" : "<>",
    "tail" : "┤├",
    "body" : "─",
}

RULER={
    "tick":"╵",
    "minor_tick":"'",
    "llabel_marker":"╰",
    "rlabel_marker":"╯"
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
    