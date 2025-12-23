
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES
from ggene.draw.colors import Colors
from ggene.draw import ScalarPlot
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs import process

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class ScalarArtistParams(BaseArtistParams):
    display_width:int = 256
    display_height:int = 8
    quantity:str = "correlation" # "correlation" , "longest_run", ""
    show_single:bool = False
    bit_depth:int = 24
    show_range:bool = True
    show_xlabel:bool = False
    fg1:int = 65
    fg2:int = 53
    bg:int = 242
    
    
class ScalarArtist(BaseArtist):
    
    _lmargin = 4
    _rmargin = 4
    
    def __init__(self, name, params:ScalarArtistParams, **kwargs):
        # self.params:ScalarArtistParams = params
        super().__init__(name, params, **kwargs)
        
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        out_lines = []
        
        num_chunks = self.params.display_width - (5 if self.params.show_range else 0)
        
        ref = window.ref_seq[:num_chunks]
        alt = window.alt_seq[:num_chunks]
        
        display_lines = self.get_scalar_plot(ref, alt, self.params.quantity, self.params.bit_depth, self.params.show_range, num_chunks, self.params.show_xlabel)
        out_lines.extend(display_lines)
        
        # out_lines = self.format_rows([[line] for line in out_lines])
        out_lines = self.fill_lines(out_lines, self.params.display_height)
        out_lines = self.join_rows(out_lines)
        out_lines.append(" ")
        out_lines.append(" ")
        
        return out_lines
    
    def get_scalar_plot(self, seqa, seqb, quantity, bit_depth, show_range, display_length, show_xlabel):
        
        qt, rcqt = self.get_quantity(seqa, seqb, quantity)
        
        xlabel = [seqa, seqb] if show_xlabel else [] 
        
        if self.params.show_single or not rcqt:
            plot_lines = self.format_scalar_plot_single(qt, bit_depth, show_range, display_length, xlabel)
        else:
            plot_lines = self.format_scalar_plot_double(qt, rcqt, bit_depth, show_range, display_length, xlabel)
        
        return plot_lines
        
    def get_quantity(self, seqa, seqb, quantity):
        
        if quantity == "correlation":
            return self._get_correlation(seqa, seqb)
        elif quantity == "longest_run":
            return self._get_longest_run(seqa, seqb)
        else:
            return [], []
    
    def _get_correlation(self, seqa, seqb, fill = 0, scale=None, step=1, keep=None, shift_step = 1):
        
        cf = lambda a,b: a==b
        sf = lambda a,b,sc: 1/sc
        
        d, rcd = process.correlate(seqa, seqb, comparison_func = cf, score_func = sf, fill = fill, scale=scale, step=step, keep=keep, shift_step=shift_step)
        
        # print(f"correlation has result length {len(d)}")
        
        return d, rcd
    
    def _get_longest_run(self, seqa, seqb, fill = 0, scale=None, step=1, keep=None):
        
        cf = lambda a,b: a==b
        
        d, _, _ = process.correlate_longest_subseq(seqa, seqb, comparison_func = cf, fill = fill, scale=scale, step=step, keep=keep)
        
        # print(f"run has result length {len(d)}")
        return d, []
    
    def format_scalar_plot_single(self, qt, bit_depth, show_range, display_length, show_xlabel):
        sc1 = ScalarPlot(qt, add_range = show_range, bit_depth=bit_depth, fg_color = self.params.fg1)
        sc1.render()
        return sc1.rows
    
    def format_scalar_plot_double(self, qt1, qt2, bit_depth, show_range, display_length, xlabel):
        
        sc1 = ScalarPlot(qt1, add_range = show_range, bit_depth=bit_depth, fg_color = self.params.fg1)
        sc2 = ScalarPlot(qt2, add_range = show_range, bit_depth=bit_depth, fg_color = self.params.fg2)
        plot_lines = ScalarPlot.show_paired(sc1, sc2, chunksz = display_length, xlabel = xlabel, suppress = True)
        
        return plot_lines
    
    