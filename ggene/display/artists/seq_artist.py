
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES
from ggene.draw.colors import Colors
from ggene.draw import highlight
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

from ggene.seqs import compare, bio

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState
    
@dataclass
class SeqArtistParams(BaseArtistParams):
    display_width = 256
    display_height = 6
    seqa_name:str = "ref" # "ref", "alt"
    seqb_name:str = "alt"
    
    display_strategy:str = "crop" # "fold"
    show_marker_line:bool = True
    highlight_details:bool = False
    highlight_match:bool = False
    highlight_rcmatch:bool = True
    
    color_bg:bool = False
    color:int = 71
    rcolor:int = 54
    
    # keep margins and related attributes as instances specific data?
    

class SeqArtist(BaseArtist):
    
    _lmargin = 4 
    _rmargin = 4
    
    def __init__(self, name, params:SeqArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        logger.debug(f"initialized artist {type(self)} with params {self.params}")
        
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        result_rows = []
        
        # seqb = window.alt_seq
        
        feats, spans = self.get_detail_features(window.features, window.start_ref)
        # if self.params.display_strategy == "crop":
        #     seqa_len = len(seqa)
        #     seqb_len = len(seqb)
        #     if seqa_len > self.params.display_width:
        #         seqa = seqa[(seqa_len - self.params.display_width)//2:(seqa_len + self.params.display_width)//2]
        #         seqb = seqb[(seqb_len - self.params.display_width)//2:(seqb_len + self.params.display_width)//2]
        
        # seqa = self.crop_seq(window.ref_seq, self.params.display_width)
        seqa = self.get_seq(self.params.seqa_name, window)
        
        if self.params.seqb_name:
            seqb = self.get_seq(self.params.seqb_name, window)
            seq_lines = self.render_dual_seq(seqa, seqb, feats, spans, state.position)
        else:    
            seq_lines = self.render_single_seq(seqa)
        
        if self.params.display_strategy == "fold":
            seq_lines = self.break_lines(seq_lines, self.params.display_width)
        
        result_rows = seq_lines
        result_rows = self.join_rows(result_rows)
        
        return result_rows
    
    def render_single_seq(self, seq):
        
        return [seq]
    
    def render_dual_seq(self, seqa, seqb, feats, spans, position):
        
        seq_lines = self.get_display_sequences(seqa, seqb, spans, self.params.highlight_details, self.params.highlight_match, self.params.highlight_rcmatch)
        
        if self.params.show_marker_line:
            marker_line = self.get_marker_line(seqa, seqb, feats, position)
            seq_lines.append(marker_line)
        
        return seq_lines
    
    def crop_seq(self, seq, display_width):
        
        seq_len = len(seq)
        if seq_len > self.params.display_width:
            seq = seq[(seq_len - display_width)//2:(seq_len + display_width)//2]
        
        return seq
    
    def get_detail_features(self, features, window_start):
        
        ftypes = ["variant","motif"]
        
        feats = []
        spans = {ft:[] for ft in ftypes}
        
        for f in features:
            
            if f.feature_type in ftypes:
                
                feats.append(f)
                spans[f.feature_type].append((f.start - window_start, f.end - window_start))
            
        return feats, spans
    
    def get_display_sequences(self, seqa, seqb, feat_spans, highlight_details, highlight_match, highlight_rcmatch):
        
        rcseqb = bio.reverse_complement(seqb)
        
        hseqa = highlight.ColoredSequence(seqa)
        hseqb = highlight.ColoredSequence(seqb)
        hrcseqb = highlight.ColoredSequence(rcseqb)
        
        mc = highlight.ColorSpec(100)
        rcmc = highlight.ColorSpec(125)
        mmc = highlight.ColorSpec(250)
        
        if highlight_details:
            hseqa.color_features(feat_spans)
            hseqb.color_features(feat_spans)
            return [hseqa.render(), hseqb.render()]
            
        elif highlight_match or highlight_rcmatch:
            
            outlist = [hseqa]
            
            if highlight_match:
                hseqa.color_matching_positions(seqb, mc, mismatch_color = mmc)
                hseqb.color_matching_positions(seqa, mc, mismatch_color = mmc)
                outlist.append(hseqb)
            
            if highlight_rcmatch:
                hseqa.color_matching_positions(rcseqb, rcmc, mismatch_color = mmc)
                hrcseqb.color_matching_positions(seqa, rcmc, mismatch_color = mmc)
                outlist.append(hrcseqb)
            
            return [cseq.render() for cseq in reversed(outlist)]
        
        else:
            return [seqa, seqb]
    
    def get_marker_line(self, seqa, seqb, features, position):
        
        motifs = [f for f in features if f.feature_type == "motif"]
        # variants = [f for f in features if f.feature_type == "variant"]
        
        variant_positions = []
        for i in range(min(len(seqa), len(seqb))):
            # A position is a variant if bases differ or if there's a gap
            is_variant = (seqa[i] != seqb[i] or 
                         seqa[i] == '-' or seqb[i] == '-')
            variant_positions.append(is_variant)
        
        marker_line = [" "]
        
        for i in range(len(seqa)):
            # Check if position is in any motif
            motif_found = None
            for motif in motifs:
                window_start = position
                motif_start = max(0, motif['start'] - window_start)
                motif_end = min(len(seqa), motif['end'] - window_start)
                
                if motif_start <= i < motif_end:
                    motif_found = motif
                    break
            
            if motif_found:
                # Choose color based on motif type
                motif_name = motif_found.get('info', {}).get('name', motif_found.get('name', ''))
                strand = motif_found.get('strand', '+')
                
                if 'splice' in motif_name.lower():
                    color = FColors.MOTIF
                elif 'tata' in motif_name.lower():
                    color = FColors.MOTIF
                elif 'cpg' in motif_name.lower():
                    color = FColors.MOTIF
                elif 'polya' in motif_name.lower():
                    color = FColors.MOTIF
                elif 'run' in motif_name.lower():
                    color = FColors.HIGHLIGHT
                    logger.debug(f"extra motif sequence: {motif_found.get("attributes",{}).get("sequence")}")
                else:
                    color = FColors.MOTIF
                    
                marker_line.append(" ")
            else:
                marker_line.append(" ")
            
        marker_sym = f"{FColors.DIM}*{FColors.RESET}"
        # Variant indicator line
        for i in range(len(marker_line)):
            if i < len(variant_positions) and variant_positions[i]:
                marker_line[i] = marker_sym
        
        return "".join(marker_line)
    
    def get_seq(self, seq_name, window, window2 = None):
        
        sn = seq_name.rstrip("2")
        do_second = seq_name.endswith("2")
        
        out_seq = ""
        
        if sn == "ref":
            if do_second:
                out_seq = window2.ref_seq
            else:
                out_seq = window.ref_seq
        elif sn == "alt":
            if do_second:
                out_seq = window2.alt_seq
            else:
                out_seq = window.alt_seq
        
        if self.params.display_strategy == "crop":
            out_seq = self.crop_seq(out_seq, self.params.display_width)
            
        return out_seq
        
    def add_seq_labels(self, seq_lines, seq_len, num_match, num_rcmatch):
        
        row_frm= self.row_frm
        # like A to B-AS .. ?
        seq_lines.insert(0, row_frm.format("",f"Duplex A-B with {num_match}/{seq_len} matching (avg {int(0.25*seq_len)})", ""))
        
        # like A to a transcript of B
        seq_lines.append(row_frm.format("", f"Duplex A-RC(B) with {num_rcmatch}/{seq_len} matching (avg {int(3/16*seq_len)})", ""))
        
        return seq_lines
