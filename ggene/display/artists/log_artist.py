
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace
import logging
from datetime import datetime

from ggene.draw.chars import HLINES
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class LogArtistParams(BaseArtistParams):
    display_width:int = 256
    display_height:int = 8

class LogArtist(BaseArtist):
    
    _lmargin = 2
    _rmargin = 2
    
    def __init__(self, name, params:LogArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        # self.logger = logging.getLogger("LogArtist")
        self.log_lines = []
    
    def log(self, message, caller):
        # self.logger.log(level, message)
        log_line = self.format_log_message(message, caller)
        self.log_lines.append(log_line)
    
    def format_log_message(self, message, caller):
        ts = datetime.now().strftime("%H:%M:%S")
        return f"{ts} | {caller}: {message}"
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        txtlines = self.get_text_lines(self.params.display_height)
        
        return txtlines
    
    # def update_logs(self):
    #     # do something with logger ... ?
    #     pass
    
    def get_text_lines(self, max_height):
        
        txtlines = []
        
        for line in reversed(self.log_lines):
            
            ols = self.wrap_line(line, max_width = self.params.display_width)
            
            if len(txtlines) + len(ols) > max_height:
                break
            
            txtlines.extend(reversed(ols))
        
        return list(reversed(txtlines))
    
    def wrap_line(self, line:str, max_width, tabs = 0):
        
        if FColors.visible_len(line) < max_width:
            return [line]
        
        outlines = []
        
        lparts = list(line.split(" "))
        
        llen = 2*tabs
        split_lines = []
        while lparts:
            lp = lparts.pop(0)
            if llen + len(lp)+1 > max_width:
                outlines.append(" ".join(split_lines))
                split_lines = ["  "*(tabs + 1)]
                llen = len(split_lines[0])
                
            split_lines.append(lp)
            llen += len(lp) + 1
        
        if split_lines:
            outlines.append(" ".join(split_lines))
        
        return outlines
