
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES
from ggene.draw.colors import Colors
from ggene.draw import highlight
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs import compare, bio

@dataclass
class MultiArtistParams:
    num_spaces:int = 1
    border_char:str = ""
    border_len:int = 80
    

class MultiArtist:
    
    def __init__(self, *artists):
        
        self.artists:List[BaseArtist] = list(artists)
        
    def __call__(self, state, window):
        return self.render(state, window)
    
    def render(self, state, window):
        
        out_lines = []
        
        for artist in self.artists:
            
            if out_lines:
                spaces = self.get_spaces(self.params.num_spaces, border_char = self.params.border_char, border_len = self.params.border_len)
                out_lines.extend(spaces)
                
            new_lines = artist.render(state, window)
            out_lines.extend(new_lines)
            
        return out_lines

    def get_spaces(self, num_spaces, border_char = "", border_len = 80):
        
        if border_char:
            num_spaces = max(3, num_spaces - 1)
        
        spaces = []
        
        for i in range(num_spaces):
            
            if border_char and i == num_spaces//2:
                spaces.append(border_char * border_len)
            else:
                spaces.append("\n")
            
        return spaces
            
            
    def set_params(self, **options):
        optdict = {}
        for k,v in options.items():
            if hasattr(self.params, k):
                optdict[k] = v
        
        self.params = replace(self.params, **optdict)
        
            