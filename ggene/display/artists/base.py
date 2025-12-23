

from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace
import logging

from ggene.draw import Colors, format_genomic, Ruler
from ggene.draw.colors import visible_slice, visible_len
from ggene.display.colors import FColors

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

logger = logging.getLogger(__name__)
logger.setLevel("WARNING")
# logger.setLevel("DEBUG")

@dataclass
class BaseArtistParams:
    display_width:int = 256
    display_height:int = 6
    
class BaseArtist:
    
    _lmargin = 4
    _rmargin = 4
    
    def __init__(self, name:str, params:BaseArtistParams, **kwargs):
        self.name = name
        self.params:BaseArtistParams = params
        self.top_label = kwargs.get("top_label","")
        self.left_label = kwargs.get("left_label","")
        self.right_label = kwargs.get("right_label","")
    
    @property
    def display_height(self):
        return self.params.display_height
    
    @property
    def display_width(self):
        return self.params.display_width
    
    def update(self, **options):
        self.set_params(**options)
    
    
    def set_params(self, **options):
        optdict = {}
        for k,v in options.items():
            if hasattr(self.params, k):
                optdict[k] = v
        
        self.params = replace(self.params, **optdict)
    
    def set_margins(self, lmargin = None, rmargin = None):
        
        if lmargin:
            self._lmargin = lmargin
        if rmargin:
            self._rmargin = rmargin
    
    def __call__(self, state:'BrowserState', window:'GenomeWindow', **kwargs) -> List[str]:
        return self.render(state, window)
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs) -> List[str]:
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
        
    def label_rows(self, b4_rows, mid_rows, a4_rows, rllabels, row_frm = ""):
        
        result_rows = []
        if not row_frm:
            row_frm= self.row_frm
        
        for eachlbls, eachrows in zip(rllabels, [b4_rows, mid_rows, a4_rows]):
            for (llbl, rlbl), row in zip(eachlbls, eachrows):
                if row:
                    if isinstance(row, str):
                        row += Colors.RESET
                    elif isinstance(row, list):
                        row.append(Colors.RESET)
                    result_rows.append(row_frm.format(llbl,"".join(row), rlbl))
        
        return result_rows
        
    def format_rows(self, rows, row_frm = ""):
        
        result_rows = []
        if not row_frm:
            row_frm= self.row_frm
        
        for row in rows:
            if row:
                if isinstance(row, str):
                    row += Colors.RESET
                elif isinstance(row, list):
                    row.append(Colors.RESET)
                result_rows.append(row_frm.format("", "".join(row), ""))
        
        return result_rows
    
    def join_rows(self, rows):
        
        result_rows = []
        for row in rows:
            if isinstance(row, list):
                result_rows.append("".join(row))
            elif isinstance(row, str):
                result_rows.append(row)
        
        return result_rows
    
    @classmethod
    def make_header(cls, chrom, start, window_size, show_reverse_strand = False, show_rna = False):
        
        lines = []

        pos = start + window_size//2
        strand = '-' if show_reverse_strand else '+'
        nctype = "RNA" if show_rna else "DNA"
        lines.append(f"{Colors.BOLD}Chromosome {chrom} | Position {pos:,} | Window {window_size}bp | ({strand}) strand | {nctype}{Colors.RESET}")
        lines.append("-"*80)
        
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
    
    def pad_rows(self, rows:List[str], to_length = None):
        
        if not to_length:
            to_length =  max([visible_len(r) for r in rows])
        
        out_rows = []
        for row in rows:
            len_diff = to_length - visible_len(row)
            prow = row.ljust(len(row) + len_diff)
            out_rows.append(prow)
        
        return out_rows
    
    @classmethod
    def make_footer(cls, output_lines, num_rows, footer_text = ""):
        
        out_lines = []
        
        for i in range(len(output_lines), num_rows):
            out_lines.append("\n")
        
        out_lines.append("-"*80)
        
        if footer_text:
            out_lines.append(footer_text)
        
        return out_lines
    
    def break_lines(self, lines, max_len):
        
        max_line = max([len(line) for line in lines])
        
        num_chunks = max_line // max_len + 1
        chunksz = max_len
        
        if num_chunks > 1:
            
            outlines = []
            
            for n in range(num_chunks):
                for eachline in lines:
                    bline = visible_slice(eachline, chunksz*n, stop = chunksz*(n+1))
                    outlines.append(bline)
            
        else:
            outlines = lines
        
        return outlines
        
    
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
    
    def collect_features(self, features, feature_types, unique_spans = False, unique_startends = False, rm_only_simple = True):
        
        spans = set()
        starts = set()
        ends = set()
        
        feats = []
        for feat in features:
            if feat and feat.feature_type in feature_types:
                
                skey = (feat.start, feat.feature_type, feat.strand)
                ekey = (feat.end, feat.feature_type, feat.strand)
                if unique_startends and skey in starts or ekey in ends:
                    continue
                
                span_key = (feat.start, feat.end, feat.feature_type, feat.strand)
                if unique_spans and span_key in spans:
                    continue
                
                if feat.source == "RepeatMasker" and feat.attributes.get("type") != "Simple_repeat":
                    continue
                
                starts.add(skey)
                ends.add(ekey)
                spans.add(span_key)
                
                feats.append(feat)

        return feats
    
    def get_feature_data(self, feat, start, display_length, display_rev = False):
        
        f_type = feat.feature_type if not "motif_len" in feat.attributes else feat.attributes.get("family_name","?")
        
        f_start = max(0, feat.attributes.get("hmm-start", feat.attributes.get("start", 0) - start))
        f_end = min(display_length, feat.attributes.get("hmm-end", feat.attributes.get("end", 0) - start))
        f_len = feat.attributes.get("motif_len", f_end - f_start)
        rev = (feat.strand == "-")^display_rev
        name = feat.name if feat.name else "(?)"
        
        return f_type, f_start, f_end, f_len, rev, name
        
    def get_display_color(self, feature):
        return FColors.get_feature_color(feature)

class ProxyArtist(BaseArtist):
    
    def __init__(self, name, params:BaseArtistParams = {}, func_handle = None, **kwargs):
        
        self.params = BaseArtistParams(name, **params)
        self.func_handle = func_handle
    
    def update_handle(self, new_func_handle:callable):
        
        if callable(new_func_handle):
            self.func_handle = new_func_handle
    
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        return self.func_handle(state, window)

