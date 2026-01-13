
from ggene.draw.color import Color, RESET
from ggene.draw.colors import Colors

SCALE = " ▁▂▃▄▅▆▇█"

def show_colors():
    
    for i in range(2**8):
        
        fgcol = Color.from_8bit(i).to_ansi()
        bgcol = Color.from_8bit(i).to_ansi(bg = True)
        
        print(i, fgcol, "this is the color of fg", SCALE, RESET, bgcol, SCALE[::-1], "this is the color of bg", RESET)
    
    print(RESET)

def main():
    show_colors()
    
if __name__=="__main__":
    main()
