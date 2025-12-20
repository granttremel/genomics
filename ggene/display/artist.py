

from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace
import logging

from ggene import draw
from ggene.draw.format import format_genomic
from ggene.draw.ruler import Ruler
from ggene.draw.chars import HLINES
from ggene.seqs import compare, bio

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.genome_browser import BrowserState

logger = logging.getLogger(__name__)
logger.setLevel("WARNING")
# logger.setLevel("DEBUG")

@dataclass
class ArtistParams:
    display_length:int = 256
    display_height:int = 16
    pass
    
class Artist:
    
    def __init__(self, params:ArtistParams, **kwargs):
        self.params:ArtistParams = params
        
    def set_params(self, **options):
        optdict = {}
        for k,v in options.items():
            if hasattr(self.params, k):
                optdict[k] = v
        
        self.params = replace(self.params, **optdict)
    
    def __call__(self, state:'BrowserState', window:'GenomeWindow') -> List[str]:
        return self.render(state, window)
    
    def render(self, state:'BrowserState', window:'GenomeWindow') -> List[str]:
        pass
    
    @property
    def row_frm(self):
        return "{:>%s}{}{:<%s}" % (str(self._lmargin), str(self._rmargin))
    
    def get_row_labels(self, num_b4, lrb4, num_mid, lrmid, num_a4, lra4):
        rllabels = [[["",""] for ii in range(i)] for i in [num_b4, num_mid, num_a4]]
        
        if num_b4 > 0:
            rllabels[0][0] = list(lrb4)
        if num_mid > 0:
            rllabels[1][0] = list(lrmid)
        if num_a4 > 0:
            rllabels[-1][0] = list(lra4)
        return rllabels
        
    def apply_row_labels(self, b4_rows, mid_rows, a4_rows, rllabels, row_frm = ""):
        
        result_rows = []
        if not row_frm:
            # row_frm = "{:>%s}{}{:<%s}" % (str(self._lmargin), str(self._rmargin))
            row_frm= self.row_frm
        
        for eachlbls, eachrows in zip(rllabels, [b4_rows, mid_rows, a4_rows]):
            for (llbl, rlbl), row in zip(eachlbls, eachrows):
                if row:
                    row.append("\x1b[0m")
                    result_rows.append(row_frm.format(llbl,"".join(row), rlbl))
        
        return result_rows
        
    
    def make_header(self, chrom, start, end, window_size, show_reverse_strand = False, show_rna = False):
        
        lines = []
        
        pos = (start + end)//2
        strand = '-' if show_reverse_strand else '+'
        nctype = "RNA" if show_rna else "DNA"
        
        lines.append(f"Chromosome {chrom} | Position {pos} | Window {window_size} | ({strand}) strand | {nctype}")
        
        return lines
    
    def fill_lines(self, lines, display_height):
        
        rem = display_height - len(lines)
        
        b4 = a4 = rem//2
        
        if rem %2 == 1:
            b4 += 1
        
        filled_lines = []
        
        for i in range(b4):
            filled_lines.append("\n")
        filled_lines.extend(lines)
        for j in range(a4):
            filled_lines.append("\n")
        
        return filled_lines
    
    def make_footer(self, output_lines, num_rows):
        
        for i in range(len(output_lines), num_rows):
            output_lines.append("\n")
        
        return output_lines
    
    def get_ruler(self, start, end, display_length):
    
        Mpart = 1e5*(start // 1e5)
        
        def format_xlabel(x):
            
            if x==start:
                return format_genomic(x)
            else:
                dx = int(x - Mpart)
                return format_genomic(dx)
                
        rlr = Ruler(start, end, num_cols = display_length, num_labels = 11, num_ticks = 0, num_minor_ticks = 2*display_length, formatter = format_xlabel)
        rlr.render()
        
        rlr_rows = rlr.get_rows()
        return [rlr_rows]
    
    def get_feature_data(self, feat, start, display_length, display_rev = False):
        
        f_type = feat.feature_type if not "motif_len" in feat.attributes else feat.attributes.get("family_name","?")
        
        f_start = max(0, feat.attributes.get("hmm-start", feat.attributes.get("start", 0) - start))
        f_end = min(display_length, feat.attributes.get("hmm-end", feat.attributes.get("end", 0) - start))
        f_len = feat.attributes.get("motif_len", f_end - f_start)
        rev = (feat.strand == "-")^display_rev
        name = feat.name if feat.name else "(?)"
        
        return f_type, f_start, f_end, f_len, rev, name

class ProxyArtist(Artist):
    
    def __init__(self, params:ArtistParams = {}, func_handle = None, **kwargs):
        
        self.params = ArtistParams(**params)
        self.func_handle = func_handle
    
    def update_handle(self, new_func_handle:callable):
        
        if callable(new_func_handle):
            self.func_handle = new_func_handle
    
    def render(self, state:'BrowserState', window:'GenomeWindow'):
        return self.func_handle(state, window)
    
    


@dataclass
class LineArtistParams(ArtistParams):
    display_length:int = 256
    feature_types:Tuple[str] = tuple(["gene","exon","CDS","dfam_hit"])
    display_height:int = 8
    show_ruler:bool = False
    show_alt:bool = False


class LineArtist(Artist):
    
    _lmargin = 4
    _rmargin = 4
    
    def __init__(self, params:LineArtistParams, **kwargs):
        self.params:LineArtistParams = params
        
        logger.debug(f"initialized artist {type(self)} with params {self.params}")
        
    def render(self, state:'BrowserState', window:'GenomeWindow'):
        
        out_lines = []
        
        start = window.start_alt if self.params.show_alt else window.start_ref
        end = window.end_alt if self.params.show_alt else window.end_ref
        
        display_feats = self.collect_features(window)
        display_lines = self.build_display_lines(display_feats, start, end, display_length = self.params.display_length, show_ruler = self.params.show_ruler)
        out_lines.extend(display_lines)
        
        logger.debug(f"built {len(display_lines)} display lines")
        
        out_lines = self.fill_lines(out_lines, self.params.display_height)
        
        return out_lines
    
    def collect_features(self, window):
        ses = set()
        feats = []
        for feat in window.features:
            if feat and feat.feature_type in self.params.feature_types:
                
                sekey = (feat.start, feat.end, feat.feature_type)
                if sekey in ses:
                    continue
                ses.add(sekey)
                
                feats.append(feat)
                
        logger.debug(f"collected {len(feats)} features from window out of {len(window.features)}")
        return feats
    
    def build_display_lines(self, feature_insts, start, end, display_length = 256, show_ruler = False, display_rev = False, max_display_height = 8):
        """
        Build display lines for repeat instances with dynamic row allocation.

        Each repeat shows:
        - Full HMM consensus extent in gray dashes
        - Actual matched region (f.start to f.end) in color with arrowheads
        """

        logger.debug(f"enter build_display_lines with show_ruler={show_ruler}")

        scale = display_length / (end - start)

        head, body, tail = HLINES.get("head"), HLINES.get("body"), HLINES.get("tail")

        # Gray color for HMM consensus extent
        gray_color = "\x1b[38;5;240m"

        # Track which columns are occupied in each row (for conflict detection)
        # row_occupancy[row_idx] = list of (start_col, end_col) tuples
        b4_row_occupancy = []
        a4_row_occupancy = []

        # Store drawing instructions for each repeat: (row_idx, repeat_obj)
        row_assignments = []

        # Assign each repeat to the topmost available row
        for f in feature_insts:
            # Calculate full HMM extent in display coordinates
            
            f_type, f_start, f_end, f_len, rev, name = self.get_feature_data(f, start, display_length, display_rev = display_rev)

            if rev:
                tgt_row_occupancy = b4_row_occupancy
            else:
                tgt_row_occupancy = a4_row_occupancy

            # Full HMM span in genomic coordinates
            hmm_full_start = f.start - f_start
            hmm_full_end = f.end + (f_len - f_end)

            # Convert to display coordinates
            hmm_full_start_d = max(0, int(scale * (hmm_full_start - start)))
            hmm_full_end_d = min(display_length - 1, int(scale * (hmm_full_end - start)))

            # Find the first row where this repeat fits
            assigned_row = None
            for row_idx, occupancy_list in enumerate(tgt_row_occupancy):
                # Check if this row has space
                has_conflict = False
                for occupied_start, occupied_end in occupancy_list:
                    # Check for overlap (with 1 char buffer for readability)
                    if not (hmm_full_end_d < occupied_start - 1 or hmm_full_start_d > occupied_end + 1):
                        has_conflict = True
                        break

                if not has_conflict:
                    assigned_row = row_idx
                    break

            # If no existing row fits, create a new one
            if assigned_row is None:
                assigned_row = len(tgt_row_occupancy)
                tgt_row_occupancy.append([])

            # Mark this region as occupied
            tgt_row_occupancy[assigned_row].append((hmm_full_start_d, hmm_full_end_d))
            row_assignments.append((assigned_row, f))
        
        mid_rows = []
        if show_ruler:
            rlr_rows = self.get_ruler(start, end, display_length)
            mid_rows.extend(rlr_rows)
            # logger.debug(f"rlr has len {len(rlr_rows)} type {type(rlr_rows[0])}, len {len(rlr_rows[0])}")
            
        # if mid_rows:
            # logger.debug(f"mid has len {len(mid_rows)} type {type(mid_rows[0])}, len {len(mid_rows[0])}")
        
        # Create rows based on how many we need
        b4_rows = [[" "] * display_length for _ in range(max(1,len(b4_row_occupancy)))]
        a4_rows = [[" "] * display_length for _ in range(max(1,len(a4_row_occupancy)))]

        # Draw each repeat
        for row_idx, f in row_assignments:
            # Get repeat attributes
            
            f_type, f_start, f_end, f_len, rev, name = self.get_feature_data(f, start, display_length, display_rev = display_rev)
            
            if rev:
                tgt_rows = b4_rows
            else:
                tgt_rows = a4_rows
            
            # Calculate positions
            # Full HMM extent
            hmm_full_start = f.start - f_start
            hmm_full_end = f.end + (f_len - f_end)
            hmm_full_start_d = max(0, int(scale * (hmm_full_start - start)))
            hmm_full_end_d = min(display_length - 1, int(scale * (hmm_full_end - start)))

            # Actual match extent
            match_start_d = max(0, int(scale * (f.start - start)))
            match_end_d = min(display_length - 1, int(scale * (f.end - start)))

            # Get color for this repeat
            col = self.get_display_color(name = name, ftype = f_type)

            # Draw full HMM extent in gray dashes
            for i in range(hmm_full_start_d, hmm_full_end_d + 1):
                if tgt_rows[row_idx][i] == " ":
                    tgt_rows[row_idx][i] = gray_color + "·"
            
            # head = "<>"
            # tail = "┤├"
            # body = "─"
            
            # Draw actual match region in color
            # Start arrowhead
            if match_start_d < display_length:
                tgt_rows[row_idx][match_start_d] = col + (head[0] if rev else tail[1])

            # End arrowhead
            if match_end_d < display_length:
                tgt_rows[row_idx][match_end_d] = (tail[0] if rev else head[1]) + "\x1b[0m"

            # Fill in the middle with dashes
            for i in range(match_start_d + 1, match_end_d):
                if i < display_length:
                    tgt_rows[row_idx][i] = body

            # Overlay repeat name in the middle
            nname = f"{name} ({f_type})"
            name_len = len(nname)
            name_start = match_start_d + int((match_end_d - match_start_d) / 2 - name_len / 2)

            for i, char in enumerate(nname):
                name_pos = name_start + i
                if match_start_d <= name_pos <= match_end_d and name_pos < display_length:
                    tgt_rows[row_idx][name_pos] = char

        num_rows = len(b4_row_occupancy) + len(mid_rows) + len(a4_row_occupancy)
        # Add color reset at end of each row
        result_rows = []
        
        rllabels = self.get_row_labels(len(b4_rows), ("3'-","-5'"), len(mid_rows), ("",""), len(a4_rows), ("5'-","-3'"))
        
        # row_frm = "{:>%s}{}{:<%s}" % (str(self._margin), str(self._margin))
        result_rows = self.apply_row_labels(b4_rows, mid_rows, a4_rows, rllabels)
        return result_rows
    
    def get_display_color(self, f=None, name = "", ftype = ""):
        
        if not name:
            name = f.get("name", "?")
        if not ftype:
           ftype = f.feature_type if not "motif_len" in f.attributes else f.attributes.get("family_name","?")
        
        if not ftype and not name:
            return "\x1b[38;5;246m"
        
        if ftype:
            if hasattr(draw.Colors, ftype.upper()):
                return str(getattr(draw.Colors, ftype.upper()))
        
        if name:
            if "Alu" in name:
                if "AluY" in name:
                    return "\x1b[38;5;142m"
                elif "AluJ" in name:
                    return "\x1b[38;5;136m"
                elif "AluS" in name:
                    return "\x1b[38;5;130m"
                else:
                    return "\x1b[38;5;124m"
            elif "L1" in name:
                return "\x1b[38;5;121m"
            elif "SVA" in name:
                if "SVA_A" in name:
                    return "\x1b[38;5;93m"
                elif "SVA_B" in name:
                    return "\x1b[38;5;92m"
                elif "SVA_C" in name:
                    return "\x1b[38;5;91m"
                elif "SVA_D" in name:
                    return "\x1b[38;5;90m"
                else:
                    return "\x1b[38;5;89m"
            else:
                return "\x1b[38;5;130m"
        
        return "\x1b[0m"

class SeqArtistParams(ArtistParams):
    show_alt:bool = False
    color:int = 71
    rcolor:int = 54
    display_length = 256
    display_height = 42
    show_seqa_bars:bool = True
    show_seqb_bars:bool = True
    color_bg:bool = False
    
    

class SeqArtist(Artist):
    
    _lmargin = 12 
    _rmargin = 4
    
    def __init__(self, params:SeqArtistParams, **kwargs):
        self.params:SeqArtistParams = params
        logger.debug(f"initialized artist {type(self)} with params {self.params}")
        
    def render(self, state:'BrowserState', window:'GenomeWindow', state2:'BrowserState' = None, window2:'GenomeWindow' = None):
        
        result_rows = []
        
        
        if self.params.show_alt:
            tgt_seqb = window.alt_seq
            tgt_seqa = window2.alt_seq
        else:
            tgt_seqb = window.ref_seq
            tgt_seqa = window2.ref_seq
        rcseqb = bio.reverse_complement(tgt_seqb)
        
        chunk_sz = self.params.display_length
        num_chks = len(tgt_seqa) // chunk_sz + 1
        
        for n in range(num_chks):
            
            tgt_subseqa = tgt_seqa[n*chunk_sz:(n+1)*chunk_sz]
            tgt_subseqb = tgt_seqb[n*chunk_sz:(n+1)*chunk_sz]
            
            seq_lines = self.get_display_sequences(tgt_subseqa, tgt_subseqb)
            seq_lines = self.format_lines(seq_lines)
            seq_lines.append("")
            
            result_rows.extend(seq_lines)
        
        # a ruler would be nice
            
        ab_match, _ = compare.count_matching(tgt_seqa, tgt_seqb, err_tol = None)
        ab_rcmatch, _ = compare.count_matching(tgt_seqa, rcseqb, err_tol = None)
        result_rows = self.add_seq_labels(result_rows, len(tgt_seqa), ab_match, ab_rcmatch)
        
        do_seqb = self.params.show_seqb_bars^self.params.show_seqa_bars 
        do_seqa = (self.params.show_seqa_bars)^self.params.show_seqb_bars 
        do_both = self.params.show_seqa_bars and self.params.show_seqb_bars
        
        bar_lines = self.add_bars(state, window, state2, window2, do_seqa = do_seqa, do_both = do_both)
        
        result_rows.extend(bar_lines)

        
        result_rows = self.fill_lines(result_rows, self.params.display_height)
        
        return result_rows
    
    
    def add_bars(self, state, window, state2, window2, do_seqa = False, do_both = True):
        
        all_bar_lines = []
        if not do_seqa or do_both:
            bars = LineArtist(LineArtistParams(show_ruler = True))
            bar_lines = bars.render(state2, window2)
            all_bar_lines.extend(bar_lines)
            logger.debug("doing seq b")
            
        if do_seqa or do_both:
            bars = LineArtist(LineArtistParams(show_ruler = True))
            bar_lines = bars.render(state, window)
            all_bar_lines.extend(bar_lines)
            logger.debug("doing seq a")
        
        return all_bar_lines
    
    def format_lines(self, seq_lines):
        
        rllabels = self.get_row_labels(1, ("SeqB 5'-","-3'"), 1, ("SeqA 5'-","-3'"), 1, ("RC SeqB 3'-","-5'"))
        result_rows = self.apply_row_labels(*seq_lines, rllabels)
        return result_rows
    
    def add_seq_labels(self, seq_lines, seq_len, num_match, num_rcmatch):
        
        row_frm= self.row_frm
        # like A to B-AS .. ?
        seq_lines.insert(0, row_frm.format("",f"Duplex A-B with {num_match}/{seq_len} matching (avg {int(0.25*seq_len)})", ""))
        
        # like A to a transcript of B
        seq_lines.append(row_frm.format("", f"Duplex A-RC(B) with {num_rcmatch}/{seq_len} matching (avg {int(3/16*seq_len)})", ""))
        
        return seq_lines
    
    def get_display_sequences(self, seqi, seqj, color = 71, rcolor = 54):
        afwd, bfwd, brev = draw.highlight_matching(seqi, seqj, colors = (color, rcolor), do_both = True, suppress = True, color_bg = self.params.color_bg)
        return [[[bfwd]], [[afwd]], [[brev]]]
    
    



class ZoomArtist(Artist):
    
    pass

@dataclass
class MultiArtistParams:
    num_spaces:int = 1
    border_char:str = ""
    border_len:int = 80
    

class MultiArtist:
    
    def __init__(self, *artists):
        
        self.artists:List[Artist] = list(artists)
        
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
        
            