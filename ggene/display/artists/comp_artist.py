
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

import numpy as np

from ggene.draw.chars import HLINES
from ggene.draw.colors import Colors
from ggene.display.colors import FColors
from ggene.draw.heatmap import Heatmap
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs import process, compare

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class CompArtistParams(BaseArtistParams):
    display_width:int = 256
    display_height:int = 8 # kernel is max twice this
    seq_name:str = "ref"
    score_mode:str = "corrs"
    mini_chunksz:int = 32
    show_stats: bool = True
    show_minimap:bool = False
    minimap_upscale: int = 2
    
    color_scheme:str = "autumn_canopy"
    mini_color_scheme:str = "malachite"
    half_block:bool = True
    
class CompArtist(BaseArtist):
    
    _lmargin = 4
    _rmargin = 4
    
    def __init__(self, name, params:CompArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        self.active_gene_name = ""
        self.active_gene = None
        self.visited = []
        
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        if type(window).__name__ == 'DualWindow':
            window2 = window.window2
        else:
            window2 = window
        
        if self.active_gene_name and not self.active_gene:
            res = self.get_active_gene(window.features)
        
        out_lines = []
        seqa = self.get_seq(window)
        seqb = self.get_seq(window2)
        
        hmrows, rchmrows, stat, rcstat = self.get_comp_heatmap(seqa, seqb)
        if self.params.show_minimap and self.active_gene:
            mm_rows = self.make_minimap(self.active_gene, state.position, state.seqb_delta, state.window_size)
        else:
            mm_rows = []
        
        full_display_lines = self.combine_rows(hmrows, rchmrows, mm_rows, stat, rcstat)
        out_lines.extend(full_display_lines)
        
        return out_lines
    
    def set_active_gene(self, active_gene_name, active_gene = None):
        self.active_gene_name = active_gene_name
        if active_gene:
            self.visited = []
            self.active_gene = active_gene
    
    def get_active_gene(self, feats):
        ag = None
        for f in feats:
            if f.name == self.active_gene_name:
                ag = f
                break
        
        if ag:
            self.visited = []
            self.active_gene = ag
        
        return ag
        
    
    def get_seq(self, window):
        
        if self.params.seq_name=="ref":
            seq = window.ref_seq
        elif self.params.seq_name == "alt":
            seq = window.alt_seq
        else:
            return ""
        seq = seq.replace("-","")
        
        return seq
    
    def get_comp_heatmap(self, seqa, seqb):
        
        scores, rcscores = self.get_heatmap_data(seqa, seqb)
        
        hm = self.format_heatmap(scores)
        rchm = self.format_heatmap(rcscores)
        
        hmrows = hm.get_rows()
        rchmrows = rchm.get_rows()
        
        # outrows = []
        # for row, rcrow in zip(hmrows, rchmrows):
        #     vlen = Colors.visible_len(row)
        #     space = int(((self.params.display_width) - 2*vlen)/3)
        #     div = " "*space
        #     row_str = div + row + div + rcrow + div
        #     outrows.append(row_str)
        
        # if self.params.show_stats:
        #     stat_strs = [self.format_stats(data) for data in [scores, rcscores]]
        #     space = int(((self.params.display_width) - 2*len(stat_strs[0])) / 3)
        #     div = " "*space
        #     outrows.append(div + stat_strs[0] + div + stat_strs[1] + div)
        
        stat = ""
        rcstat = ""
        if self.params.show_stats:
            stat = self.format_stats(scores)
            rcstat = self.format_stats(rcscores)
        
        return hmrows, rchmrows, stat, rcstat
        
    def get_heatmap_data(self, seqa, seqb):
        nchksa = round(len(seqa) / self.params.mini_chunksz) + 1
        nchksb = round(len(seqb) / self.params.mini_chunksz) + 1
        score_func = compare.get_score_function(self.params.score_mode)
        scores, rcscores = compare.calc_sequence_comparison(seqa, seqb, nchksa, nchksb, self.params.mini_chunksz, score_func)
        scores = [row for row in scores if row]
        rcscores = [row for row in rcscores if row]
        return scores, rcscores
        
    def format_heatmap(self, scores):
        
        hm = Heatmap(scores, center = None, minval = 0, maxval = 1.0, row_labels = [], color_scheme = self.params.color_scheme, col_space = 0, row_space = 0, add_middle = False, half_block = self.params.half_block , colorbar = True)
        
        return hm

    def format_stats(self, scores):
        mean = np.mean(scores)
        sd = np.std(scores)
        min = np.min(scores)
        max = np.max(scores)
        
        return f"Mean={mean:0.3f}, sd={sd:0.3f}, min={min:0.3f}, max={max:0.3f}"

    def make_minimap(self, active_gene, pos, offset, window_sz):
    
        us = self.params.minimap_upscale
    
        start = active_gene.start
        end = active_gene.end
        
        rel_end = end - start
        rel_pos = pos - start
        
        num_blocks = int(us * rel_end / window_sz)

        data = np.zeros((num_blocks, num_blocks))
        
        for aind, bind in self.visited:
            
            # aind = int(us * a / window_sz)
            # bind = int(us * b / window_sz)
            
            if aind < num_blocks and bind < num_blocks:
                data[bind, aind] = 0.5
        
        
        aind_curr = max(0, min(num_blocks, int(us * rel_pos / window_sz)))
        bind_curr = max(0, min(num_blocks, int(us * (rel_pos + offset) / window_sz)))
        
        # pos_tup = (rel_pos, rel_pos+offset)
        pos_tup = (aind_curr, bind_curr)
        if 0 <= aind_curr < num_blocks and 0 <= bind_curr < num_blocks:
            data[bind_curr, aind_curr] = 1.0
        
        if not pos_tup in self.visited:
            self.visited.append(pos_tup)
        
        cols = ([20,20,20], [80,80,100], [170, 170, 0])
            
        hm = Heatmap(data, half_block = self.params.half_block, minval = 0, maxval = 1.0,  col_space = 0, row_space = 0, colors = cols, add_middle = False, colorbar = False)
        
        return hm.get_rows()

    def combine_rows(self, hm_rows, rchm_rows, minimap_rows, stat, rcstat):
        
        display_width = self.params.display_width
        
        if minimap_rows:
            display_width -= Colors.visible_len(minimap_rows[0])
        
        outrows = []
        for row, rcrow in zip(hm_rows, rchm_rows):
            vlen = Colors.visible_len(row) + Colors.visible_len(rcrow)
            space = int(((display_width) - vlen)/3)
            div = " "*space
            row_str = div + row + div + rcrow + div
            outrows.append(row_str)
        
        if stat or rcstat:
            vlen = len(stat) + len(rcstat)
            space = int(((display_width) - vlen)/3)
            div = " "*space
            print(stat, rcstat)
            row_str = div + stat + div + rcstat + div
            outrows.append(row_str)
        
        if len(outrows) < len(minimap_rows):
            for i in range(len(outrows), len(minimap_rows)):
                outrows.append(" "*len(outrows[0]))
        
        
        print(f"outrows {len(outrows)} mm {len(minimap_rows)}")
        
        if minimap_rows:
            for i in range(len(minimap_rows)):
                outrows[i] += minimap_rows[i]
        
        return outrows

    @classmethod
    def display_artist(cls, seqa, seqb, **kwargs):
        
        params = CompArtistParams(**kwargs)
        
        @dataclass
        class FakeGenomeWindow:
            ref_seq:str
            alt_seq:str
        
        gw = FakeGenomeWindow(seqa, seqa)
        
        artist = CompArtist("temp_artist", params, kernel = seqb)
        
        rnd = artist.render(None, gw)
        for row in rnd:
            print(row)
