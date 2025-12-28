
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING
from dataclasses import dataclass, field, replace

from ggene.draw.chars import HLINES, FSYMS
from ggene.draw.colors import Colors
from ggene.display.colors import FColors
from ggene.display.artists.base import BaseArtistParams, BaseArtist, logger

if TYPE_CHECKING:
    from ggene.database.genome_iterator import GenomeWindow
    from ggene.browser.genome_browser import BrowserState

@dataclass
class LineArtistParams(BaseArtistParams):
    display_width:int = 256
    use_global_features:bool = True
    feature_types:Tuple[str] = ("gene","exon","CDS","dfam_hit","repeat")
    rm_only_simple_repeats:bool = True
    display_height:int = 8
    show_ruler:bool = False
    show_alt:bool = False
    # Label track settings
    orientation:str = "out" # "in"
    show_label_tracks:bool = True       # Enable label tracks for small features
    min_label_width:int = 20            # Min display chars to show inline label
    min_arrow_width:int = 3             # Below this, use single arrow char
    arrowhead_style:str = "head" # head_filled, head_triangle, head_harpoon, etc.


class LineArtist(BaseArtist):
    
    _lmargin = 4
    _rmargin = 4
    
    def __init__(self, name, params:LineArtistParams, **kwargs):
        super().__init__(name, params, **kwargs)
        
        
    def render(self, state:'BrowserState', window:'GenomeWindow', **kwargs):
        
        out_lines = []
        
        start = window.start_alt if self.params.show_alt else window.start_ref
        end = window.end_alt if self.params.show_alt else window.end_ref
        
        if self.params.use_global_features and state.feature_types:
            fts = state.feature_types
        else:
            fts = self.params.feature_types
        
        display_feats = self.collect_features(window.features, fts, unique_startends = True)
        display_lines = self.build_display_lines(display_feats, start, end, self.params.display_width)
        out_lines.extend(display_lines)
        
        out_lines = self.join_rows(out_lines)
        
        return out_lines

    
    def build_display_lines(self, feature_insts, start, end, display_width, display_rev=False):
        """
        Build display lines for features with dynamic row allocation.

        Features are displayed with:
        - Full HMM consensus extent in gray dashes (for dfam_hit)
        - Actual matched region in color with arrowheads
        - Label tracks for small features
        - Vertical connectors from labels to features
        - Larger features (genes, transcripts) closer to central ruler
        
        ooh okay what if.. "nested" features get alternating colors .. ?
        
        """
        logger.debug(f"enter build_display_lines with show_ruler={self.params.show_ruler}")

        scale = display_width / (end - start)

        # Separate features by strand
        rev_features = []
        fwd_features = []

        for f in feature_insts:
            f_type, f_start, f_end, f_len, is_rev, name = self.get_feature_data(
                f, start, display_width, display_rev=display_rev
            )
            if is_rev:
                rev_features.append(f)
            else:
                fwd_features.append(f)

        # Sort features by size (largest first) so they get rows closer to ruler
        rev_features = self._sort_by_size(rev_features, start, scale)
        fwd_features = self._sort_by_size(fwd_features, start, scale)

        # Build each strand track independently
        # Returns: (feature_rows, label_rows) ordered from near-ruler to far-from-ruler
        rev_feat_rows, rev_label_rows = self._build_strand_track(
            rev_features, start, end, display_width, scale, display_rev, is_reverse=True
        )
        fwd_feat_rows, fwd_label_rows = self._build_strand_track(
            fwd_features, start, end, display_width, scale, display_rev, is_reverse=False
        )

        # Build ruler rows
        mid_rows = []
        if self.params.show_ruler:
            rlr_rows = self.get_ruler(start, end, display_width)
            mid_rows.extend(rlr_rows)

        # Assemble rows:
        # Reverse: label_rows (far) + feat_rows (near) -> flip so near is at bottom
        # Forward: feat_rows (near) + label_rows (far) -> near is at top
        rev_combined = rev_feat_rows + rev_label_rows  # far to near

        fwd_combined = fwd_feat_rows + fwd_label_rows  # near to far (near at top)

        if self.params.orientation == "out":
            rev_combined.reverse()  # now near to far (near at bottom, adjacent to ruler)
            # Center/pad to display height
            rev_combined, fwd_combined = self._center_tracks(rev_combined, fwd_combined)
        elif self.params.orientation == "in":
            fwd_combined.reverse()
            # fwd_combined, rev_combined = self._center_tracks(fwd_combined, rev_combined)
            fwd_combined, rev_combined = self._center_tracks(rev_combined, fwd_combined)

        result_rows = rev_combined + mid_rows + fwd_combined
        return result_rows

    def _sort_by_size(self, features, start, scale):
        """Sort features by display size, largest first."""
        def get_display_width(f):
            match_start_d = max(0, int(scale * (f.start - start)))
            match_end_d = int(scale * (f.end - start))
            return match_end_d - match_start_d

        return sorted(features, key=get_display_width, reverse=True)

    def _build_strand_track(self, features, disp_start, disp_end, display_width, scale,
                            display_rev, is_reverse):
        """
        Build feature and label rows for a single strand.

        Args:
            features: List of features (already sorted by size)
            start, end: Genomic coordinates
            display_width: Width in characters
            scale: Pixels per base
            display_rev: Whether display is reversed
            is_reverse: Whether this is the reverse strand

        Returns:
            (feature_rows, label_rows) - both ordered from near-ruler to far-from-ruler
        """
        # Get drawing characters
        arrowhead_key = getattr(self.params, 'arrowhead_style', 'head')
        head = HLINES.get(arrowhead_key, HLINES.get("head"))
        body = HLINES.get("body")
        tail = HLINES.get("tail")
        small = HLINES.get("small")
        tiny = HLINES.get("tiny")
        cont = HLINES.get("cont")
        vline = HLINES.get("vline", "│")
        vline_up = HLINES.get("vline_up", "╵")
        vline_down = HLINES.get("vline_down", "╷")

        gray_color = "\x1b[38;5;240m"
        dim_color = "\x1b[2m"

        min_label_width = getattr(self.params, 'min_label_width', 0)
        min_arrow_width = getattr(self.params, 'min_arrow_width', 3)
        show_label_tracks = getattr(self.params, 'show_label_tracks', True)

        # Track row occupancy
        feat_occupancy = []   # Feature rows (near ruler)
        label_occupancy = []  # Label rows (far from ruler)

        # Store assignments for drawing
        feat_assignments = []   # (row_idx, feature, needs_label)
        label_assignments = []  # (label_row, feat_row, feature, label_start, connector_col, label_text)
        
        shorter = False
        if len(features) > 10:
            shorter = True

        # Assign features to rows (larger features processed first, get row 0)
        for f in features:
            f_type, f_start, f_end, f_len, _, name = self.get_feature_data(
                f, disp_start, display_width, display_rev=display_rev
            )

            # Calculate display coordinates
            hmm_start_d, hmm_end_d, match_start_d, match_end_d = self._get_display_coords(
                f, f_start, f_end, f_len, disp_start, scale, display_width
            )
            
            hmm_start_d = match_start_d
            hmm_end_d = match_end_d
            
            # feature_width = match_end_d - match_start_d
            feature_width = match_end_d - match_start_d
            needs_label = show_label_tracks and feature_width < min_label_width

            # Assign to feature row
            feat_row = self._find_available_row(feat_occupancy, hmm_start_d, hmm_end_d)
            feat_occupancy[feat_row].append((hmm_start_d, hmm_end_d))
            feat_assignments.append((feat_row, f, needs_label))

            # Assign label if needed
            if needs_label:
                label_text = self.format_feature_name(f, f_type, name, f.strand, display_rev=display_rev, shorter = True)
                label_len = len(label_text)
                label_start = max(0, match_start_d + (feature_width // 2) - (label_len // 2))
                label_end = min(display_width - 1, label_start + label_len)

                label_row = self._find_available_row(label_occupancy, label_start, label_end)
                label_occupancy[label_row].append((label_start, label_end))

                connector_col = match_start_d + (feature_width // 2)
                label_assignments.append((label_row, feat_row, f, label_start, connector_col, label_text))

        # Create row arrays
        num_feat_rows = max(1, len(feat_occupancy))
        num_label_rows = len(label_occupancy)
        feat_rows = [[" "] * display_width for _ in range(num_feat_rows)]
        label_rows = [[" "] * display_width for _ in range(num_label_rows)]

        # Draw features
        for feat_row, f, needs_label in feat_assignments:
            self._draw_feature(
                feat_rows[feat_row], f, disp_start, scale, display_width, display_rev,
                is_reverse, needs_label, head, body, tail, small, tiny, cont, gray_color, min_arrow_width
            )

        # Draw labels and connectors
        for label_row, feat_row, f, label_start, connector_col, label_text in label_assignments:
            self._draw_label_and_connector(
                feat_rows, label_rows, label_row, feat_row,
                f, label_start, connector_col, label_text,
                vline, vline_up, vline_down, dim_color, display_width
            )
        
        # if self.params.orientation == "in":
        #     feat_rows = list(reversed(feat_rows))
        #     label_rows = list(reversed(label_rows))

        return feat_rows, label_rows

    # def _get_display_coords(self, f, f_start, f_end, f_len, start, scale, display_width):
    #     """Calculate display coordinates for a feature."""
    #     hmm_full_start = f.start - f_start
    #     hmm_full_end = f.end + (f_len - f_end)
    #     hmm_start_d = max(0, int(scale * (hmm_full_start - start)))
    #     hmm_end_d = min(display_width - 1, int(scale * (hmm_full_end - start)))
    #     match_start_d = max(0, int(scale * (f.start - start)))
    #     match_end_d = min(display_width - 1, int(scale * (f.end - start)))
    #     return hmm_start_d, hmm_end_d, match_start_d, match_end_d
    
    def _get_display_coords(self, f, f_start, f_end, f_len, disp_start, scale, display_width):
        """Calculate display coordinates for a feature."""
        # hmm_full_start = f.start - f.get("hmm_start", 0)
        # hmm_full_end = f.start + (f_len - f.get("hmm_end", 0))
        hmm_full_start = f.start
        hmm_full_end = f.end
        hmm_start_d = max(0, int(scale * (hmm_full_start - disp_start)))
        hmm_end_d = min(display_width - 1, int(scale * (hmm_full_end - disp_start)))
        match_start_d = max(0, int(scale * (f.start - disp_start)))
        match_end_d = min(display_width - 1, int(scale * (f.end - disp_start)))
        return hmm_start_d, hmm_end_d, match_start_d, match_end_d

    def _draw_feature(self, row, f, start, scale, display_width, display_rev,
                      is_reverse, needs_label, head, body, tail, small, tiny, cont, gray_color, min_arrow_width):
        """Draw a single feature onto a row."""
        f_type, f_start, f_end, f_len, _, name = self.get_feature_data(
            f, start, display_width, display_rev=display_rev
        )
        hmm_start_d, hmm_end_d, match_start_d, match_end_d = self._get_display_coords(
            f, f_start, f_end, f_len, start, scale, display_width
        )

        # Draw HMM extent in gray
        for i in range(hmm_start_d, hmm_end_d + 1):
            if 0 <= i < display_width and row[i] == " ":
                row[i] = gray_color + "⋯" # "·"
        
        # hmm_start_d = match_start_d
        # hmm_end_d = match_end_d
        
        feature_width = match_end_d - match_start_d
        col = self.get_display_color(f) + FColors.BOLD
        
        # Very small: single arrow
        if feature_width < min_arrow_width:
            if f.end - f.start <= 1:
                arrow_char = tiny[0] if is_reverse else tiny[1]
            else:
                arrow_char = small[0] if is_reverse else small[1]
            mid_pos = match_start_d + feature_width // 2
            if 0 <= mid_pos < display_width:
                row[mid_pos] = col + arrow_char + Colors.RESET
        else:
            # Draw line with arrowheads
            if 0 < match_start_d < display_width:
                row[match_start_d] = col + (head[0] if is_reverse else tail[1])
            elif match_start_d == 0:
                row[match_start_d] = col + (cont[0] if is_reverse else cont[1])
            
            if 0 < match_end_d < display_width-1:
                row[match_end_d] = (tail[0] if is_reverse else head[1]) + Colors.RESET
            elif match_end_d == display_width-1:
                row[match_end_d] = (cont[0] if is_reverse else cont[1]) + Colors.RESET
            
            for i in range(match_start_d + 1, match_end_d):
                if 0 <= i < display_width:
                    row[i] = body[0]
            
            for j in range(hmm_start_d, hmm_end_d):
                if 0 <= j < display_width:
                    if not row[j]:
                        row[j] = "⋯"
                
                pass

            # Overlay name if not using label track
            if not needs_label:
                nname = self.format_feature_name(f, f_type, name, f.strand, display_rev=display_rev)
                name_len = len(nname)
                name_start = match_start_d + (feature_width // 2) - (name_len // 2)
                for i, char in enumerate(nname):
                    pos = name_start + i
                    if match_start_d <= pos <= match_end_d and 0 <= pos < display_width:
                        row[pos] = char

    def _draw_label_and_connector(self, feat_rows, label_rows, label_row, feat_row,
                                   f, label_start, connector_col, label_text,
                                   vline, vline_up, vline_down, dim_color, display_width):
        """Draw label text and vertical connector between label and feature."""
        # Draw label text
        for i, char in enumerate(label_text):
            pos = label_start + i
            if 0 <= pos < display_width:
                label_rows[label_row][pos] = dim_color + char

        if not (0 <= connector_col < display_width):
            return

        # Draw connector from feature row through to label row
        # Connector goes through: feat_rows[feat_row+1:] then label_rows[0:label_row]
        # then into label_row itself

        # Through feature rows (below the feature, toward labels)
        for row_idx in range(feat_row + 1, len(feat_rows)):
            if feat_rows[row_idx][connector_col] == " ":
                feat_rows[row_idx][connector_col] = dim_color + vline + Colors.RESET

        # Through label rows (above the target label row)
        for row_idx in range(0, label_row):
            if label_rows[row_idx][connector_col] == " ":
                label_rows[row_idx][connector_col] = dim_color + vline + Colors.RESET

        # Connector endpoint in label row
        if label_rows[label_row][connector_col] == " ":
            label_rows[label_row][connector_col] = dim_color + vline_down + Colors.RESET

    def _find_available_row(self, occupancy_list, start_d, end_d):
        """Find first available row without overlap, creating new row if needed."""
        for row_idx, occupied_regions in enumerate(occupancy_list):
            has_conflict = False
            for occ_start, occ_end in occupied_regions:
                if not (end_d < occ_start - 1 or start_d > occ_end + 1):
                    has_conflict = True
                    break
            if not has_conflict:
                return row_idx

        new_row = len(occupancy_list)
        occupancy_list.append([])
        return new_row

    def _center_tracks(self, rev_rows, fwd_rows):
        """Pad tracks to fill half the display height each."""
        fig_height = self.params.display_height
        half_height = fig_height // 2

        # Pad reverse rows (padding goes at top, features at bottom near ruler)
        if len(rev_rows) < half_height:
            padding = [[" "] * len(rev_rows[0]) if rev_rows else [" "]
                       for _ in range(half_height - len(rev_rows))]
            rev_rows = padding + rev_rows

        # Pad forward rows (features at top near ruler, padding at bottom)
        if len(fwd_rows) < half_height:
            width = len(fwd_rows[0]) if fwd_rows else 1
            padding = [[" "] * width for _ in range(half_height - len(fwd_rows))]
            fwd_rows = fwd_rows + padding

        return rev_rows, fwd_rows
    
    def get_hmm_extent(self, f):
        
        if f.feature_type != "dfam_hit":
            return f.start, f.start, f.end, f.end
        else:
            start_hmm = f.start - f.attributes.get("start_hmm", f.start)
            end_hmm = f.attributes.get("end_hmm", f.end) + (f.end - f.start)
            return start_hmm, f.start, f.end, end_hmm
    
    def format_feature_name(self, f, feat_type, feat_name, feat_strand, display_rev = False, shorter = False):
        
        paren = []
        
        name = "?"
        
        if feat_type in ['gene','transcript', 'exon', 'CDS']:    
            
            if feat_type == 'gene':
                name = feat_name
                if not display_rev and feat_strand == '-':
                    paren.append("-")
                elif display_rev and feat_strand == "+":
                    paren.append("+")
            
            elif feat_type == 'transcript':
                name = self._format_transcript(f)
            elif feat_type == 'exon':
                name = self._format_exon(f)
            elif feat_type == 'CDS':
                name = self._format_CDS(f)
        
        elif feat_type == "variant":
            return self._format_var(f)
        elif feat_type == "dfam_hit":
            name = feat_name
        else:  
            name = feat_name
            paren.append(feat_type)
            
        if shorter:
            paren = []
            
        if paren:
            return f" {name} ({", ".join(paren)}) "
        else:
            return f" {name} "
    
    def _format_gene(self, f):
        
        bt = f.attributes.get("gene_biotype","")
        
        if bt == "protein_coding":
            return f.name
        elif bt == "lncRNA":
            return "lncRNA"
        else:
            return f"{f.id} ({bt})"
        
        pass
    
    def _format_transcript(self, f):
        
        gn = f.attributes.get("gene_name", "?")
        tname = f.attributes.get("transcript_name","?")
        
        return tname
    
    def _format_exon(self, f):
        
        tn = f.attributes.get("transcript_name","?")
        exnum = f.attributes.get("exon_number", -1)
        
        return f"{tn}, exon {exnum}"
    
    def _format_CDS(self, f):
        
        tn = f.attributes.get("transcript_name","?")
        exnum = f.attributes.get("exon_number", -1)
        
        return f"{tn}, CDS {exnum}"
    
    def _format_var(self, f):
        
        ref = f.attributes.get("ref","")
        alt = f.attributes.get("alt",[""])[0]
        gt =f.attributes.get("genotype","")
        
        return f"{ref}->{gt}"
        