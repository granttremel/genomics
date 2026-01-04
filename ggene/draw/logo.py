


from ggene.draw.chars import SCALE
from ggene.draw.colors import Colors




def show_consensus_logo(seq_composition, colors = {}):
    
    if not colors:
        colors = {b:c for b, c in zip("ATGC",[64, 174, 21, 155])}
    
    upper = "AT"
    lower = "GC"
    bg = "AC"
    fg = "TG"
    
    color_strs = {b:Colors.get_color(c, background = b in bg) for b,c in colors.items()}
    
    
    # can it be done with block characters .. ?    
    
