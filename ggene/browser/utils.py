

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


