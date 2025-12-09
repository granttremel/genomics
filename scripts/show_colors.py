


SCALE = " ▁▂▃▄▅▆▇█"

def show_colors():
    
    creset = f"\x1b[0m"
    
    for i in range(2**8):
        
        fgcol = f"\x1b[38;5;{i}m"
        bgcol = f"\x1b[48;5;{i}m"
        
        print(i, fgcol, "this is the color of fg", SCALE, creset, bgcol, SCALE[::-1], "this is the color of bg", creset)
        

def main():
    show_colors()
    
if __name__=="__main__":
    main()
