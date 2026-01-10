
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

import numpy as np

from ggene.draw.chars import HLINES
from ggene.draw.colors import Colors
from ggene.display.colors import FColors
from ggene.draw.heatmap import Heatmap
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs import process

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class InterArtistParams(BaseArtistParams):
    display_width:int = 256
    display_height:int = 8 # kernel is max twice this
    seq_name:str = "ref"
    color_scheme:str = "terracotta"
    focus:str = "all" # "N", "SW", "RY", "MK"
    show_rc:bool = False
    show_both:bool = False
    shear:bool = False
    show_kernel:bool = False
    show_target:bool = False
    
class InterArtist(BaseArtist):
    
    _lmargin = 4
    _rmargin = 4
    
    focuses = ["N","SW","RY","MK"]
    
    def __init__(self, name, params:InterArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        self.kernel = process.crop_sequence(kwargs.get("kernel", "ATAT"*8), 2*self.params.display_height)
        
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        out_lines = []
        seq = self.get_seq(window)
        display_lines = self.get_inter_heatmap(seq)
        out_lines.extend(display_lines)
        
        return out_lines
    
    def get_seq(self, window):
        
        if self.params.seq_name=="ref":
            seq = window.ref_seq
        elif self.params.seq_name == "alt":
            seq = window.alt_seq
        else:
            return ""
        
        if len(seq) > self.params.display_width:
            seq_hlen = len(seq)//2
            tgt_hlen = self.params.display_width//2
            seq_crop = seq[seq_hlen-tgt_hlen:seq_hlen+tgt_hlen]
        else:
            seq_crop = seq
        
        
        return seq_crop
    
    def get_inter_heatmap(self, seq):
        
        mat, rcmat = self.get_inter_matrix(seq)
        rows = []
        
        cols = self.get_colors(num_cols = int(np.max(mat) - np.min(mat))+1)
        
        row_labels = list(self.kernel) if self.params.show_kernel else []
        col_labels = list(seq) if self.params.show_target else []
        
        if not self.params.show_rc or self.params.show_both:
            
            hm_rows = self.get_heatmap_rows(mat, cols, row_labels = row_labels, col_labels = col_labels, label = "Forward")
            rows.extend(hm_rows)
            
        if self.params.show_rc or self.params.show_both:
            if rows:
                rows.append("")
            hm_rows = self.get_heatmap_rows(rcmat, cols, row_labels = row_labels, col_labels = col_labels, label = "Reverse Complement")
            rows.extend(hm_rows)
            
        return rows
    
    def get_heatmap_rows(self, mat, cols, row_labels = [], col_labels = [], label=  ""):
        
        rows = []
        if label:
            rows = [label]
        
        hm = Heatmap(mat, half_block = True, colobar = False, colors = cols, row_labels = row_labels, col_labels = col_labels)
        rows.extend(hm.get_rows())

        return rows
        
    def get_colors(self, num_cols):
        
        c1,c2 = Colors.get_color_scheme_24b(self.params.color_scheme)
        
        if num_cols > 2:
            
            _, mid, _ = Colors.add_middle(c1, c2)
            
            if np.std(mid) < 15:
                mid = [20,20,120]
            
            cols = [c1]
            hues = np.linspace(0, 360, num_cols)
            
            for h in hues[1:-1]:
                c = Colors._adjust_color(mid, hue_shift = h)
                cols.append(c)
            
            cols.append(c2)
        else:
            cols = [c1,c2]
        
        return cols
    
    def get_inter_matrix(self, seq):
        
        mat = process.get_full_inter_matrix(self.kernel, seq)
        
        if self.params.focus=="all":
            inds = np.array([0,1,2,3])[:,None,None,None]
            mat = np.sum(mat*inds, axis = 0)
            
        elif self.params.focus in self.focuses:
            find = self.focuses.index(self.params.focus)
            mat = mat[find]
        else:
            mat = mat[0]
                
        if self.params.shear:
            mat = process.shear_matrix(mat, fill = np.nan)
        
        mat = self.crop_mat(mat, self.params.display_height, self.params.display_width)
        
        return mat[0], mat[1]
    
    def crop_mat(self, mat, max_height, max_width):
        
        _, h, w = mat.shape
        
        if h > 2*max_height:
            mat = mat[:, h//2-max_height:h//2+max_height, :]
        
        if w > max_width:
            mat = mat[:, :, w-max_width//2:w+max_width//2]
        
        return mat

    @classmethod
    def display_artist(cls, seqa, seqb, **kwargs):
        
        params = InterArtistParams(**kwargs)
        
        @dataclass
        class FakeGenomeWindow:
            ref_seq:str
            alt_seq:str
        
        gw = FakeGenomeWindow(seqa, seqa)
        
        artist = InterArtist("temp_artist", params, kernel = seqb)
        
        rnd = artist.render(None, gw)
        for row in rnd:
            print(row)
