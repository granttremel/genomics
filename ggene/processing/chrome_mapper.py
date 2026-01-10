
from typing import TYPE_CHECKING

import time

import numpy as np

from ggene.seqs.lambdas import lambda_map, needs_features, get_sources, get_feature_types
from ggene.seqs.bio import reverse_complement

from ggene.draw import add_ruler, ScalarPlot

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager

class ChromeMapper:
    max_indices = {'1': 248937043, '10': 133778498, '11': 135075908, '12': 133238549, '13': 114346637, '14': 106879812, '15': 101979093, '16': 90222678, '17': 83240391, '18': 80247514, '19': 58599303, '2': 242175634, '20': 64327972, '21': 46691226, '22': 50799123, '3': 198228376, '4': 190195978, '5': 181472430, '6': 170745977, '7': 159233377, '8': 145066516, '9': 138320835, 'MT': 16023, 'X': 156027877, 'Y': 57214397}
    def __init__(self, seq_gen, feat_gen):
        self.seq_gen = seq_gen
        self.feat_gen = feat_gen
    
    def get_chromosomal_quantities(self, chr, seq_specs, chunksz = 10e6, start = 1e6, length = None, resample = False, do_rc = False):
        
        seq_fns = []
        for seq_spec in seq_specs:
            if isinstance(seq_spec, str):
                seq_fn = lambda_map.get(seq_spec)
                seq_fns.append(seq_fn)
            elif callable(seq_spec):
                seq_fn = seq_spec
                seq_fns.append(seq_fn)
        
        if not seq_fns:
            return None, None
        
        # if not needs_feats:
        #     needs_feats = []
        #     for sfn in seq_fns:
        #         needs_feats += needs_features(sfn)
        # needs_feats = list(set([v for v in needs_feats if v is not None]))
        
        chrmax = self.max_indices.get(str(chr), 10e6)
        if not length:
            length = int(chrmax - start)
        chunksz = int(chunksz)
        num_chunks = int(length/chunksz)+1
        step = chunksz
        starts = [[] for qt in seq_fns]
        qts = [[] for qt in seq_fns]
        
        chrstr = str(chr)
        
        full_seq = self.seq_gen(chrstr, start, start + length)
        
        for n in range(num_chunks):
            end = start + step
            seq = full_seq[n*chunksz:(n+1)*chunksz]
            if do_rc:
                seq = reverse_complement(seq)
            
            ufeats =[uf for uf in self.feat_gen(chrstr, start, start+step)]
            # if needs_feats:
            #     ufeats =[uf for uf in self.feat_gen(chrstr, start, start+step)]
            # else:
            #     ufeats = []
            
            if len(seq) < 1:
                start = end
                continue
            
            for i,seq_fn in enumerate(seq_fns):
                qt = seq_fn(seq, ufeats)
                starts[i].append(start)
                qts[i].append(qt)
            start = end
            
            # if n%int(num_chunks/10)==0:
            #     print(f"n={n}/{num_chunks}")
            
        if resample:
            rsqts = [[] for i in range(len(seq_fns))]
            for n in range(1, num_chunks - 1):
                for i in range(len(seq_fns)):
                    rsqts[i].append(0.25*qts[i][n-1] + 0.5*qts[i][n] + 0.25*qts[i][n+1])
            qts = rsqts
            starts = [st[1:num_chunks-1] for st in starts]
        
        return qts, starts
    
    def get_chromosomal_quantity(self, chr, seq_spec, chunksz = 10e6, start = 1e6, length = None, resample = False, do_rc = False, needs_feats = []):
        qts, starts = self.get_chromosomal_quantities(chr, [seq_spec], chunksz = chunksz, start = start, length = length, resample = resample, do_rc = do_rc, needs_feats = needs_feats)
        return qts[0], starts[0]
    
    def display_chromosomal_quantity(self, chr, seq_spec, chunksz = 10e6, start = 1e6, max_disp = 256, length = None, resample = False, **kwargs):
        
        from ggene.genome import chrome_lengths
        if not length:
            length = chrome_lengths.get(str(chr), 10e6) - start
        
        show_hist = kwargs.get("show_hist",False)
        
        qts, starts = self.get_chromosomal_quantity(chr, seq_spec, chunksz=chunksz, start = start, length = length, resample = resample)
        
        if qts is not None:
            qt_name = str(seq_spec)
        else:
            return
        
        lines = self.get_chrom_quantity_lines(chr, qts, qt_name, chunksz = chunksz, start = start, max_disp = max_disp, length = length, suppress = False, resample = resample, **kwargs)
        
        if show_hist:
            
            data = qts
            
            minval = kwargs.get("minval", min(data))
            maxval = kwargs.get("minval", max(data))
            
            hist, bins = np.histogram(data, bins = min(max_disp, int(np.sqrt(len(data)))), range = (minval, maxval))
            
            sc1 = ScalarPlot(hist, minval = 0, add_range = True, xmin = bins[0], xmax = bins[-1], ruler = True, num_labels = 5, ticks = 0, minor_ticks = 0)
            
            print(f"histogram of {seq_spec}")
            sc1.show()
            
        return lines
    
    def display_chromosomal_quantities(self, chr, seq_specs, chunksz = 10e6, start = 1e6, max_disp = 256, length = None, resample = False, **kwargs):
        
        if not length:
            length = self.gene_map.max_indices[str(chr)] - start
        
        fg_colors = []
        
        minvals = kwargs.pop("minvals",[])
        maxvals = kwargs.pop("maxvals",[])
        
        qts, starts = self.get_chromosomal_quantities(chr, seq_specs, chunksz=chunksz, start = start, length = length, resample = resample)
        qt_names = [str(seq_spec) for seq_spec in seq_specs]
        num_qts = len(qts)
        
        all_lines=[]
        
        nice_colors = [65,65-12,173,169,60,131,55,138,91]
        qt_names_col = []
        for n in range(len(seq_specs)):
            qt_names_col.append(f"\x1b[38;5;{nice_colors[n]}m{qt_names[n]}\x1b[0m")
            
            minval = None
            if minvals:
                minval = minvals[n]
            maxval = None
            if maxvals:
                maxval = maxvals[n]
            
            do_flip = bool(n%2)
            lines = self.get_chrom_quantity_lines(chr, qts[n], qt_names[n], chunksz=chunksz, start=start, max_disp=max_disp, length=length, suppress = True, flip = do_flip, fg_color = nice_colors[n], minval=minval, maxval=maxval, **kwargs)
            all_lines.append(lines[1:])
        
        # print(qts[0][:16],qts[1][:16],qts[2][:16])
        
        qt_name_str = f" |  [{", ".join(qt_names_col)}]"
        num_chunks = length//chunksz
        num_disp_chunks = max(1, int(num_chunks // max_disp))
        lines_per_disp = 6
        
        for nsect in range(num_disp_chunks):
            for nqt in range(num_qts):
                if nqt == 0:
                    print(all_lines[nqt][nsect*lines_per_disp] + qt_name_str)
                for nl in range(1, 4):
                    print(all_lines[nqt][nsect*lines_per_disp + nl])
                if nqt == num_qts-1:
                    print(all_lines[nqt][nsect*lines_per_disp + 4])
                    print()
            
    
    def get_chrom_quantity_lines(self, chr,  qts, qt_name, chunksz = 10e6, start = 1e6, max_disp = 256, length = None, suppress = False, **kwargs):
        
        num_chunks = length//chunksz
        num_disp_chunks = max(1, int(num_chunks // max_disp))
        
        minval = kwargs.pop("minval", 0)
        maxval = kwargs.pop("maxval", max(qts))
        
        lines = []
        lines.append(f"{qt_name} for chr{chr}")
        for i in range(num_disp_chunks):
            seq_ran = int(length / num_disp_chunks)
            seq_start = start + i*seq_ran
            
            gcs_crop = qts[i*max_disp:(i+1)*max_disp]
            if len(gcs_crop) == 0:
                break
            
            gc_bars = ScalarPlot(gcs_crop, 
                                 minval = minval, maxval = maxval, bit_depth = 24, add_range = True, 
                                 xmin = seq_start, xmax = seq_start+seq_ran, num_labels = 10, ticks=2, minor_ticks = 5, genomic = True, ruler = True,
                                 **kwargs)
            
            tick_dists = gc_bars.ruler.get_tick_distances()
            major = tick_dists.get("tick_interval")
            minor = tick_dists.get("minor_tick_interval")
            
            # gc_bars = draw.scalar_to_text_nb(gcs_crop, minval = minval, maxval = maxval, bit_depth = 24, add_range = True, **kwargs)
            # gc_bars, dists = draw.add_ruler(gc_bars, seq_start, seq_start + seq_ran, num_labels = 10, ticks = 2, minor_ticks = 5, genomic = True)
            lines.append(f"Section {i+1} | {seq_start/1e6:.1f}M - {(seq_start+seq_ran)/1e6:.1f}M (major = {major/1e3:2.0f}k, minor = {minor/1e3:2.0f}k)")
            lines.extend(gc_bars.get_rows())
            lines.append("\n")
        
        if not suppress:
            for line in lines:
                print(line)
        return lines

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

    # @classmethod
    # def get_generators(cls, gm:'GenomeManager', seq_specs):
        
    #     sg = lambda chrom, start, end: gm.get_sequence(chrom, start, end)
        
    #     # fg_streams = {}
    #     # for sp in seq_specs:
            
    #     #     needs = needs_features(sp)
            
    #     #     for need in needs:
    #     #         if need not in fg_streams:
    #     #             fg_streams[need] = gm.annotations.streams.get(need, None)
        
    #     fg = lambda chrom, start, end, feature_types: gm.annotations.stream_by_types(feature_types, chrom, start, end)
        
    #     return sg, fg