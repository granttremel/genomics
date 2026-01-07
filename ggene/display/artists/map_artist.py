
from typing import List, Tuple, Dict, Any, Union, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import MAP
from ggene.draw.colors import Colors
from ggene.draw import ScalarPlot

from ggene.database.cache import GenomeCache
from ggene.processing.chrome_mapper import ChromeMapper

from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs.lambdas import needs_features, lambda_map, get_feature_types, get_sources

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

# logger.setLevel("DEBUG")

@dataclass
class MapArtistParams(BaseArtistParams):
    display_width:int = 64
    display_height:int = 6
    scale:float = 1024 # -1 for full chromosome. larger = zoom out
    telo_buffer:int = int(1e6)
    quantity:Union[Any,str] = "gc" # "gc", "genes", ..?
    quantity2:Union[Any,str] = ""
    bit_depth:int = 16
    show_range:bool = True
    show_ruler:bool = True
    show_fella:bool = False
    show_paired:bool = True
    qt_label:str = ""
    qt2_label:str = ""
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
        self.top_label = kwargs.get("top_label", self.get_top_label())
        
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
        
        # logger.debug(f"num_chunks for minimap {self.name}: {num_chunks}, chunksz {chunksz:0.1f}")
        qts = [self.params.quantity]
        if self.params.quantity2:
            qts.append(self.params.quantity2)
        
        marker_pos = max(0, min(num_chunks - 1, int(num_chunks * (state.position - start) / length)))
        
        map_lines = self.get_map_lines(self.seq_gen, self.feat_gen, qts, state.chrom, start, length, chunksz, num_chunks = num_chunks, marker_pos = marker_pos)
        
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
    
    def get_map_lines(self, seq_gen, feat_gen, qts, chrom, start, length, chunksz, num_chunks = None, marker_pos = None):

        data = self.get_data(qts, seq_gen, feat_gen, chrom, start, length, chunksz, num_chunks = num_chunks)
        # logger.debug(f"data for minimap {self.name}: {len(data)}, {len(data[0])}, with length {length} so num_chunks {length/chunksz:0.1f}")
        
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
    
    def get_data(self, qts, seq_gen, feat_gen, chrom, start, length, chunksz, num_chunks = None):

        datas = {}
        ssps = qts.copy()

        for seq_spec in qts:

            if seq_spec in self.samplers:

                sampler = self.samplers[seq_spec]

                data = sampler.sample(chrom, start, length, chunksz, num_chunks = num_chunks)
                datas[seq_spec] = data
                ssps.remove(seq_spec)

            else:

                cseq_specs = self.gnm_cache.list_cached_specs(chrom)

                if seq_spec in cseq_specs:

                    sampler = self.gnm_cache.get_sampler(seq_spec)
                    self.samplers[seq_spec] = sampler

                    data = sampler.sample(chrom, start, length, chunksz, num_chunks = num_chunks)
                    datas[seq_spec] = data
                    ssps.remove(seq_spec)

        if ssps:
            data = self.compute_data(ssps, seq_gen, feat_gen, chrom, start, length, chunksz, num_chunks = num_chunks)
            datas.update(data)

        return [datas.get(qt, []) for qt in qts]
    
    def compute_data(self, qts, seq_gen, feat_gen, chrom, start, length, chunksz, num_chunks = None):

        import numpy as np

        cm = ChromeMapper(seq_gen, feat_gen)
        
        if self.current_pos in self._data_cache:
            new_data = self._data_cache[self.current_pos]
        else:
            new_data, _ = cm.get_chromosomal_quantities(chrom, qts, chunksz=chunksz, start=start, length=length)
            self.cache_data(new_data)
        
        stitched = new_data
        
        #TODO: bring back the stitching functionality
        # stitched, fetch_start, fetch_end = self.stitch_data(start, length, len(qts), num_chunks = num_chunks)
        # fetch_len = fetch_end - fetch_start
        
        # logger.debug(f"stitched returns {stitched[0]}")
        # # logger.debug(f"after stitch, fetch {fetch_start} {fetch_end} {fetch_len}")
        
        # if fetch_len > 0:
        #     # new_data, _ = cm.get_chromosomal_quantities(chrom, qts, chunksz=chunksz, start=start, length=length)
        #     fetch_start_gnm = start+int(fetch_start * chunksz)
        #     fetch_len_gnm = int(fetch_len * chunksz)
        #     logger.debug(f"fetching from {fetch_start_gnm} - {fetch_start_gnm+fetch_len_gnm} ({fetch_len_gnm})")
        #     new_data, _ = cm.get_chromosomal_quantities(chrom, qts, chunksz=chunksz, start=fetch_start_gnm, length=fetch_len_gnm)
            
        #     logger.debug(f"new data returns {new_data[0]}")
        #     # logger.debug(f"new data {len(new_data)}, {len(new_data[0])}")
            
        #     for chk in range(0, fetch_len):
        #         for i in range(len(qts)):
        #             stitched[i][chk + fetch_start] = new_data[i][chk]
            
        #     self.cache_data(stitched)
            
        #     logger.debug(f"cached stitched data at {self.current_pos}")
            
        #     # logger.debug(f"stitched data from {fetch_start} length {fetch_len} with chunk start {fetch_start}, num chunks {fetch_len}")
        #     # logger.debug(f"new_data: len {len(new_data[0])}, {sum(1 for s in new_data[0] if s!=0)} nonzero")

        # logger.debug(f"stitched: len {len(stitched[0])}, {sum(1 for s in stitched[0] if s!=0)} nonzero")
        
        # Resample to exact num_chunks if needed
        if num_chunks is not None:
            resampled_data = []
            for data in stitched:
                if len(data) != num_chunks:
                    x_old = np.linspace(0, 1, len(data))
                    x_new = np.linspace(0, 1, num_chunks)
                    data = np.interp(x_new, x_old, data)
                resampled_data.append(data)
            stitched = resampled_data

        datas = {qt: data for qt, data in zip(qts, stitched)}
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
    def get_generators(cls, gm:'GenomeManager', seq_specs = [], use_sources = True, use_features = False):
        
        seq_gen = lambda chr, start, end: gm.get_sequence(chr, start, end = end)
        if use_sources:
            
            sources = []
            for seq_spec in seq_specs:
                sources.extend(get_sources(seq_spec))
            sources = list(set(sources))
            
            # logger.debug(f"making feat gen with sources {sources}")
            
            feat_gen = lambda chr, start, end: gm.annotations.stream_by_sources(sources, chr, start, end)
        elif use_features:
            
            ftypes = []
            for seq_spec in seq_specs:
                fts = get_feature_types(seq_spec)
                ftypes.extend(fts)
            
            ftypes = list(set(ftypes))
            
            feat_gen = lambda chr, start, end: gm.annotations.stream_by_types(ftypes, chr, start, end=end)
        
        else:
            print("idk")
            feat_gen = lambda chr, start, end: gm.annotations.stream_all(chr, start, end)
        
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
        self._data_cache[self.current_pos] = data
    
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

        for cached_pos, cached_data in self._data_cache.items():
            # The cached range is [cached_pos - length/2, cached_pos + length/2]
            cached_start = cached_pos - length // 2

            # Calculate the offset in base pairs
            offset = target_start - cached_start

            # Skip if no overlap
            if abs(offset) >= length:
                continue
            
            # logger.debug(f"retreiving data from cache at {cached_pos} containing {cached_data}")

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
    
