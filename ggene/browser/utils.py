

import random
import sys
import tty
import termios


def get_user_keypress():
    
    """Get a single keypress from the user."""
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        key = sys.stdin.read(1)
        
        # Check for escape sequences (arrow keys and ctrl combinations)
        if key == '\x1b':
            next_chars = sys.stdin.read(2)
            key += next_chars
            
            prefix = suffix = ""
            # Check for Ctrl+Arrow combinations
            if next_chars == '[1':
                # Read the next character to identify Ctrl+Arrow
                extra = sys.stdin.read(2)
                if extra == ';5':
                    prefix = "ctrl_"
                elif extra == ';2':
                    prefix = "shift_"
                
                final = sys.stdin.read(1)
                if final == 'C':  # Ctrl+Right
                    suffix = 'right'
                elif final == 'D':  # Ctrl+Right
                    suffix = 'left'
                elif final == 'A':  # Ctrl+Right
                    suffix = 'up'
                elif final == 'B':  # Ctrl+Right
                    suffix = 'down'
                elif final in "c":
                    suffix = final
                
                if prefix and suffix:
                    return prefix+suffix
                    
                key += extra
        
        return key
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)



def startup_splash(display_height, display_width, splash_type = "canabo", **kwargs):
    
    print("\x1b[?1049h", end="")  # enter different terminal view sta
    
    if kwargs.get("skip", True):
        return
    
    if splash_type == "canabo":
        _splash_canabo(display_height, display_width)
    
    input("Press enter to continue...")
    
def _splash_canabo(display_height, display_width):
    cac = [chr(i) for i in range(0x1400, 0x1680)]
    
    def get_fill(n):
        return "".join([random.choice(cac) for i in range(n)])
    
    
    display_height += display_height%2
    
    msg= "genome browser"
    
    nh = (display_width//2-len(msg))//2 
    nb = (display_width//4 - 1)
    ns = display_height//4
    nf = (display_height - ns)//2
    
    for i in range(nf):
        print(get_fill(display_width))
    
    for i in range(ns):
        if i==ns//2-1:
            print(get_fill(nb) + " "*nh + msg + " "*nh + get_fill(nb))
        else:
            print(get_fill(nb) + " "*(display_width//2), get_fill(nb))
        
    for i in range(nf):
        print(get_fill(display_width))
    
