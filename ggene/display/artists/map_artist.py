
from typing import List, Tuple, Dict, Any, Union, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import MAP
from ggene.draw.colors import Colors
from ggene.draw import ScalarPlot

from ggene.database.cache import GenomeCache
from ggene.processing.chrome_mapper import ChromeMapper

from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs.lambdas import needs_features, lambda_map

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class MapArtistParams(BaseArtistParams):
    display_width:int = 64
    display_height:int = 6
    scale:float = 1024 # -1 for full chromosome. larger = zoom out
    telo_buffer:int = int(1e6)
    quantity:Union[Any,str] = "gc" # "gc", "genes", ..?
    quantity2:Union[Any,str] = "genes"
    bit_depth:int = 16
    show_range:bool = True
    show_ruler:bool = True
    show_fella:bool = False
    show_paired:bool = True
    fg1:int = 52
    fg2:int = 96
    bg:int = 242
    
    
class MapArtist(BaseArtist):
    
    _lmargin = 4
    _rmargin = 4
    
    pointer = MAP.get("pointer_top","")
    fella = MAP.get("fella","")
    
    def __init__(self, name, params:MapArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        
        self.seq_gen = kwargs.get("sequence_generator", None)
        self.feat_gen = kwargs.get("feature_generator", None)
        self.current_pos = None
        self.current_chrom = ""
        self._data_cache = {}
        self.top_label = self.get_top_label()
        
        self.gnm_cache = GenomeCache()
        self.samplers = {}
            
            
    def get_top_label(self):
        
        if self.params.scale > 0:
            return f"Minimap of {self.params.quantity}/{self.params.quantity2} ({self.params.scale}x)"
        else:
            return f"Overmap of {self.params.quantity}/{self.params.quantity2}"
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        out_lines = []
        
        if not self.seq_gen or not self.feat_gen:
            return []
        
        max_length = ChromeMapper.max_indices.get(state.chrom, 10e6) - 2*self.params.telo_buffer
        
        if self.params.scale < 0:
            length = max_length
            self.current_pos = self.params.telo_buffer + length//2
            current_pos = self.update_position(state.chrom, self.current_pos, length)
        else:
            if self.params.scale < 1:
                length = int(self.params.scale * max_length)
            else:
                length = int(state.window_size * self.params.scale)
        
            if self.current_pos is None:
                self.current_pos = state.position
            
            current_pos = self.update_position(state.chrom, state.position, length)
        
        start = max(self.params.telo_buffer, current_pos - length // 2)
        length = min(max_length, length)
        
        num_chunks = self.params.display_width - (5 if self.params.show_range else 0)
        chunksz = int(length / num_chunks)
        
        qts = [self.params.quantity]
        if self.params.quantity2:
            qts.append(self.params.quantity2)
        
        marker_pos = max(0, min(num_chunks - 1, int(num_chunks * (state.position - start) / length)))
        
        map_lines = self.get_map_lines(self.seq_gen, self.feat_gen, qts, state.chrom, start, length, chunksz, marker_pos = marker_pos)
        
        map_lines = ["".join(ml) for ml in map_lines]
        
        out_lines.extend(map_lines)
        
        marker_row = self.get_marker_row(marker_pos, num_chunks)
        out_lines.insert(0, marker_row)
        
        if self.params.show_ruler:
            ruler_lines = self.get_ruler(start, start + length, display_length = num_chunks)
            out_lines.extend(ruler_lines)
        
        out_lines = self.join_rows(out_lines)
        out_lines = self.pad_rows(out_lines)
        
        return out_lines
    
    def update(self, **kwargs):
        super().update(**kwargs)
        self._data_cache = {}
    
    def get_map_lines(self, seq_gen, feat_gen, qts, chrom, start, length, chunksz, marker_pos = None):

        data = self.get_data(qts, seq_gen, feat_gen, chrom, start, length, chunksz)
        
        if len(qts) > 1:
            if self.params.show_paired:
                lines = self.format_scalar_plot_double(data[0], data[1])
            else:
                lines = []
                for n in range(len(qts)):
                    plot_lines = self.format_scalar_plot_single(data[n])
                    lines.extend(plot_lines)
                
        else:
            lines = self.format_scalar_plot_single(data[0])
        
        if marker_pos and self.params.show_fella:
            lines = self.add_fella(lines, marker_pos)
            
        return lines
    
    def get_data(self, qts, seq_gen, feat_gen, chrom, start, length, chunksz):
        
        datas = {}
        ssps = qts.copy()
        
        for seq_spec in qts:
            
            if seq_spec in self.samplers:
                
                sampler = self.samplers[seq_spec]
                
                data = sampler.sample(chrom, start, length, chunksz)
                datas[seq_spec] = data
                ssps.remove(seq_spec)
                # logger.debug(f"sampler returned data {type(data)}, {len(data)} with {type(data[0])}")
            
            else:
                
                cseq_specs = self.gnm_cache.list_cached_specs(chrom)
                
                if seq_spec in cseq_specs:
                
                    sampler = self.gnm_cache.get_sampler(seq_spec) 
                    self.samplers[seq_spec] = sampler
                    
                    data = sampler.sample(chrom, start, length, chunksz)
                    datas[seq_spec] = data
                    ssps.remove(seq_spec)
                    # logger.debug(f"sampler returned data {type(data)}, {len(data)} with {type(data[0])}")
        
        if ssps:    
            data = self.compute_data(ssps, seq_gen, feat_gen, chrom, start, length, chunksz)
            datas.update(data)
        
        return [datas.get(qt, []) for qt in qts]
    
    def compute_data(self, qts, seq_gen, feat_gen, chrom, start, length, chunksz):
        
        cm = ChromeMapper(seq_gen, feat_gen)
        
        new_data, _ = cm.get_chromosomal_quantities(chrom, qts, chunksz=chunksz, start=start, length=length)
        
        datas = {qt:data for qt, data in zip(qts, new_data)}
        return datas
        
    def get_marker_row(self, marker_index, num_chunks):

        marker = [" "] * num_chunks

        marker[marker_index] = Colors.HIGHLIGHT + Colors.BOLD + self.pointer + Colors.RESET

        return ["".join(marker)]
        
    def format_scalar_plot_single(self, qt):
        min_val = 0
        sc1 = ScalarPlot(qt, add_range = self.params.show_range, bit_depth=self.params.bit_depth, fg_color = self.params.fg1, minval = min_val)
        sc1.render()
        return sc1.get_rows()
    
    def format_scalar_plot_double(self, qt1, qt2):
        
        min_val = 0
        sc1 = ScalarPlot(qt1, add_range = self.params.show_range, bit_depth=self.params.bit_depth, fg_color = self.params.fg1, minval = min_val)
        
        min_val = 0
        max_val = None
        sc2 = ScalarPlot(qt2, add_range = self.params.show_range, bit_depth=self.params.bit_depth, fg_color = self.params.fg2, minval = min_val, maxval = max_val)
        plot_lines = ScalarPlot.show_paired(sc1, sc2, chunksz = self.params.display_width, suppress = True)
        
        return plot_lines
    
    def add_fella(self, sc_rows, marker_pos):
        
        cfella = self.fella
        for i in range(len(sc_rows)//2-1, -1, -1):
            row = sc_rows[i]
            curr_char = FColors.scrub_codes(FColors.visible_slice(row, start = marker_pos, stop = marker_pos + 1))
            
            if curr_char == " ":
                
                col = FColors(fg = 148, bg = 234).code
                prevcol = FColors(fg = self.params.fg1, bg = 234).code
                cfella = col + self.fella + prevcol
                
                new_row = FColors.visible_slice(row, start = 0, stop = marker_pos) + cfella + FColors.visible_slice(row, start = marker_pos+1, stop = FColors.visible_len(row))
                sc_rows[i] = new_row
                break
            elif curr_char == "█":
                
                col = FColors(fg = 148, bg = self.params.fg1).code
                prevcol = FColors(fg = self.params.fg1, bg = 234).code
                cfella = col + self.fella + prevcol
                
                new_row = FColors.visible_slice(row, start = 0, stop = marker_pos) + cfella + FColors.visible_slice(row, start = marker_pos+1, stop = FColors.visible_len(row))
                sc_rows[i] = new_row
                break
        
        return sc_rows


    @classmethod
    def get_generators(cls, gm:'GenomeManager', seq_specs = [], needs_feats = []):
        
        for seq_spec in seq_specs:
            feats = needs_features(seq_spec)
            needs_feats.extend(feats)
        
        needs_feats = list(set(needs_feats))
        
        seq_gen = lambda chr, start, end: gm.get_sequence(chr, start, end = end)
        feat_gen = lambda chr, start, end, feature_types: gm.annotations.stream_by_types(feature_types, chr, start, end=end)
        
        return seq_gen, feat_gen
    
    @classmethod
    def full_map(cls, qt1, qt2, seq_gen, feat_gen, **kwargs):
        
        params = MapArtistParams(quantity1 = qt1, quantity2 = qt2)
        kwargs["scale"] = -1
        params.update(**kwargs)
        
        artist = MapArtist(params, sequence_generator = seq_gen, feature_generator = feat_gen)
        
        return artist
    
    ############ no longer used but could be useful tbh ########################
    
    def cache_data(self, data):
        self._data_cache[(self.current_chrom, self.current_pos)] = data
    
    def stitch_data(self, target_position, length, num_qts, num_chunks):
        """
        Attempts to reuse cached data by stitching together overlapping chunks.

        Returns:
            stitched: Partially filled data array
            start_ind: First chunk index that needs to be fetched
            end_ind: Last chunk index (+1) that needs to be fetched
        """
        # num_chunks = self.params.display_width
        stitched = [[0] * num_chunks for n in range(num_qts)]

        # Initially assume we need to fetch everything
        need_fetch_start = 0
        need_fetch_end = num_chunks

        # The target range is [target_position - length/2, target_position + length/2]
        target_start = target_position - length // 2

        for (chrom, cached_pos), cached_data in self._data_cache.items():
            # The cached range is [cached_pos - length/2, cached_pos + length/2]
            cached_start = cached_pos - length // 2

            # Calculate the offset in base pairs
            offset = target_start - cached_start

            # Skip if no overlap
            if abs(offset) >= length:
                continue

            # Calculate chunk offset (how many chunks the data has shifted)
            chunk_offset = round(offset * num_chunks / length)

            # Determine which chunks from cached data we can reuse
            if chunk_offset > 0:
                # Target is to the right of cached, reuse right portion of cache
                copy_from = chunk_offset
                copy_to = num_chunks
                paste_at = 0
                paste_end = num_chunks - chunk_offset

                if copy_to > copy_from and paste_end > paste_at:
                    for n in range(num_qts):
                        stitched[n][paste_at:paste_end] = cached_data[n][copy_from:copy_to]
                    need_fetch_start = max(need_fetch_start, paste_end)

            elif chunk_offset < 0:
                # Target is to the left of cached, reuse left portion of cache
                copy_from = 0
                copy_to = num_chunks + chunk_offset
                paste_at = -chunk_offset
                paste_end = num_chunks

                if copy_to > copy_from and paste_end > paste_at:
                    for n in range(num_qts):
                        stitched[n][paste_at:paste_end] = cached_data[n][copy_from:copy_to]
                    need_fetch_end = min(need_fetch_end, paste_at)

            else:
                # Exact match - reuse everything
                for n in range(num_qts):
                    stitched[n] = cached_data[n][:]
                need_fetch_start = num_chunks
                need_fetch_end = 0

        return stitched, need_fetch_start, need_fetch_end
    
    def update_position(self, chrom, new_pos, length, thresh = 0.4):
        
        if abs(new_pos - self.current_pos) / length > thresh:
            self.current_pos = new_pos
            self.current_chrom = chrom
        
        return self.current_pos
    
