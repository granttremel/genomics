
from ggene.draw.colors import Colors

SCALE = " ▁▂▃▄▅▆▇█"

def show_colors():
    
    creset = "\x1b[1m"
    print(creset)
    
    for i in range(2**8):
        
        fgcol = f"\x1b[38;5;{i}m"
        bgcol = f"\x1b[48;5;{i}m"
        
        print(i, fgcol, "this is the color of fg", SCALE, creset, bgcol, SCALE[::-1], "this is the color of bg", creset)
        # print(i, fgcol, "this is the color of fg", SCALE, bgcol, SCALE[::-1], "this is the color of bg", creset)
        # print(fgcol, "hiii hello", creset)
        # print(repr(fgcol + "hiii" + creset))
    
    print(creset)

def main():
    # Colors.show_colors()
    show_colors()
    
if __name__=="__main__":
    main()
