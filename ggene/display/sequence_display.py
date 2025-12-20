"""
Sequence display module for genome browser.

This module handles the rendering of DNA/RNA sequences with variant highlighting,
codon annotations, and motif underlining.
"""

from typing import List, Dict, Tuple, Optional, Any
from ggene import draw
from ggene.display.colors import Colors
from ggene.processing.coordinate_mapper import CoordinateMapper, SequenceRenderer
from ggene.processing.sequence_processor import SequenceProcessor
from ggene.database.annotations import UFeature
import logging

logger = logging.getLogger(__name__)


class SequenceDisplay:
    """Renders genomic sequences with annotations for terminal display."""

    def __init__(self, render_state, window_size: int = 240, ):
        """Initialize the sequence display renderer.

        Args:
            window_size: Size of the display window
        """
        self.rstate = render_state
        self.window_size = window_size
        self.coord_mapper = CoordinateMapper(window_size)
        self.seq_renderer = SequenceRenderer()
        self.seq_processor = SequenceProcessor()

    def render_sequences(self, ref_seq: str, personal_seq: str,
                         features: List = None, state=None) -> Tuple[List[str], int]:
        """Render reference and personal sequences with annotations.

        Args:
            ref_seq: Reference sequence
            personal_seq: Personal/alternative sequence
            features: List of genomic features
            state: Browser state object

        Returns:
            Tuple of (list of display lines, number of sequences displayed)
        """
        if not ref_seq and not personal_seq:
            return [], 0

        # Update settings based on state
        if state:
            self.coord_mapper.update_settings(
                window_size=state.window_size,
                show_reverse_strand=state.show_reverse_strand
            )
            self.seq_renderer.update_settings(
                show_reverse_strand=state.show_reverse_strand,
                show_rna=state.show_rna
            )

        # Process sequences
        aligned_ref = ref_seq
        aligned_pers = personal_seq

        # Detect variants
        variant_positions = self._detect_variant_positions(aligned_ref, aligned_pers)

        # Extract codon and motif features
        codons = self._extract_codon_features(features, state)
        motifs = self._extract_motif_features(features, state)

        # Generate display lines
        lines = []

        # Add ruler
        ruler_lines = self._generate_ruler(state, aligned_ref)
        if ruler_lines:
            lines.extend(ruler_lines)

        # Add reference sequence
        if aligned_ref:
            ref_line = self._render_sequence_line(
                aligned_ref, variant_positions, codons, motifs, [],
                is_reference=True, label="Ref"
            )
            lines.append(ref_line)

        # Add personal sequence
        if aligned_pers:
            pers_line = self._render_sequence_line(
                aligned_pers, variant_positions, codons, motifs, [],
                is_reference=False, label="Alt"
            )
            lines.append(pers_line)

        marker_row = self._generate_marker_row(variant_positions, motifs, len(aligned_ref))
        lines.append(marker_row)

        return lines, 2 if aligned_ref and aligned_pers else 1

    def _get_margin(self, lbl = "", color=""):
        fill = self.rstate.left_margin - len(lbl)
        lbl_fill = lbl + " "*fill
        margin_str = f"{Colors.BOLD}{color}{lbl_fill}{Colors.RESET}"
        return margin_str
    
    def _detect_variant_positions(self, ref_seq: str, alt_seq: str) -> List[bool]:
        
        variant_positions = []
        min_len = min(len(ref_seq), len(alt_seq))

        for i in range(min_len):
            # A position is a variant if bases differ or if there's a gap
            is_variant = (ref_seq[i] != alt_seq[i] or
                          ref_seq[i] == '-' or alt_seq[i] == '-')
            variant_positions.append(is_variant)

        return variant_positions

    def _extract_codon_features(self, features: List, state) -> Dict[str, List[Tuple[int, int]]]:
        
        start_codons = []
        stop_codons = []

        if not features or not state:
            return {'start': start_codons, 'stop': stop_codons}

        window_start = state.position

        for feature in features:
            ftype = self._get_feature_type(feature)

            if ftype == 'start_codon':
                positions = self._map_feature_to_display(feature, window_start, state)
                if positions:
                    start_codons.append(positions)

            elif ftype == 'stop_codon':
                positions = self._map_feature_to_display(feature, window_start, state)
                if positions:
                    stop_codons.append(positions)

        return {'start': start_codons, 'stop': stop_codons}

    def _extract_motif_features(self, features: List, state) -> List[Dict]:
        
        motifs = []

        if not features or not state:
            return motifs

        window_start = state.position

        for feature in features:
            ftype = self._get_feature_type(feature)

            if ftype == 'motif':
                positions = self._map_feature_to_display(feature, window_start, state)
                if positions:
                    motif_info = {
                        'start': positions[0],
                        'end': positions[1],
                        'type': getattr(feature, 'name', 'unknown')
                    }
                    motifs.append(motif_info)

        return motifs

    def _map_feature_to_display(self, feature, window_start: int, state) -> Optional[Tuple[int, int]]:
        
        # Calculate raw positions relative to window
        raw_start = max(0, feature.start - window_start)
        raw_end = min(state.window_size, feature.end - window_start + 1)

        # Apply strand-aware rendering
        if state.show_reverse_strand:
            display_start = state.window_size - raw_end
            display_end = state.window_size - raw_start
        else:
            display_start = raw_start
            display_end = raw_end

        # Ensure positions are within bounds
        display_start = max(0, min(display_start, state.window_size))
        display_end = max(0, min(display_end, state.window_size))

        if display_start < display_end:
            return (display_start, display_end)
        return None

    def _render_sequence_line(self, seq: str, variants: List[bool],
                               codons: Dict, motifs: List, highlights, 
                               is_reference: bool, label: str) -> str:
        
        # Build colored sequence
        colored_seq = []
        
        marg = self._get_margin(label)
        
        # v_spans = [v.get("")]
        # start_cd_spans = {(cd.start, cd.end):Colors.START_CODON for cd in codons.items()}
        # stop_cd_spans = {(cd.start, cd.end):Colors.START_CODON for cd in codons.items()}
        
        m_spans = {(m.get("start"), m.get("end")):Colors.MOTIF for m in motifs}
        
        col_seq = draw.highlight_sequence_by_span(seq, m_spans)
        # starts = 
        
        # for i, base in enumerate(seq):
        #     # Check if position is in a codon
        #     in_start_codon = any(start <= i < end for start, end in codons.get('start', []))
        #     in_stop_codon = any(start <= i < end for start, end in codons.get('stop', []))
            
            
            
        #     # Apply coloring based on priority
        #     if i < len(variant_positions) and variant_positions[i]:
        #         # Variant coloring
        #         if base == '-':
        #             colored_seq.append(f"{Colors.DELETION}-{Colors.RESET}")
        #         elif is_reference and seq[i] != '-':
        #             colored_seq.append(f"{Colors.SNP}{base}{Colors.RESET}")
        #         else:
        #             colored_seq.append(f"{Colors.INSERTION}{base}{Colors.RESET}")
        #     elif in_start_codon:
        #         # Start codon coloring (green background)
        #         colored_seq.append(f"{Colors.START_CODON}{base}{Colors.RESET}")
        #     elif in_stop_codon:
        #         # Stop codon coloring (red background)
        #         colored_seq.append(f"{Colors.STOP_CODON}{base}{Colors.RESET}")
        #     else:
        #         # Normal base
        #         colored_seq.append(base)

        # Format line with label
        return f"{marg}{col_seq}"


    def _generate_ruler(self, state, aligned_ref):
        rlines = []
        # Position ruler (adjusted for alignment)
        
        strand = '-' if state.show_reverse_strand else '+'
        ruler = self._get_margin("Position:", color=Colors.POSITION)
        pos_counter = 0  # Track actual genomic positions despite gaps
        for i in range(len(aligned_ref)):
            if aligned_ref[i] != '-':  # Only count non-gap positions
                if pos_counter % 10 == 0:
                    ruler += "|"
                elif pos_counter % 5 == 0:
                    ruler += strand
                else:
                    ruler += "."
                pos_counter += 1
            else:
                ruler += " "  # Space for gaps
        rlines.append(ruler)
        return rlines
    
    def _generate_marker_row(self, variant_positions:List[bool], motifs:List[Dict], seq_len:int):
        
        marker_line = [" "]*seq_len
        # Add variant markers
        if variant_positions:
            var_line = self._generate_variant_markers(variant_positions)

        # Add motif underlines
        if motifs:
            motif_line = self._generate_motif_underlines(motifs, seq_len)
        
        for i in range(seq_len):
            if var_line[i]:
                marker_line[i] = var_line[i]
            elif motif_line[i]:
                marker_line[i] = motif_line[i]
        
        return f"{self.rstate._get_margin()}{"".join(marker_line)}"
        

    def _generate_variant_markers(self, variant_positions: List[bool]) -> str:
        """Generate line showing variant positions.

        Args:
            variant_positions: List of variant flags

        Returns:
            Marker line string
        """
        markers = []
        for is_variant in variant_positions:
            markers.append('*' if is_variant else ' ')

        return markers

    def _generate_motif_underlines(self, motifs: List[Dict], seq_len: int) -> Optional[str]:
        """Generate underline indicators for motifs.

        Args:
            motifs: List of motif dictionaries
            seq_len: Length of sequence

        Returns:
            Motif underline string or None
        """
        if not motifs:
            return None

        underlines = [' '] * seq_len
        
        left, fill, right = list(self.rstate.arrow_disp)
        
        for motif in motifs:
            mstrand = motif.get("strand", "+")
            start = motif.get('start', 0)
            end = motif.get('end', seq_len)
            for i in range(start, min(end, seq_len)):
                underlines[i] = fill
            if mstrand == "+":
                underlines[end - 1] = right
            else:
                underlines[start] = left

        return underlines

    def _get_feature_type(self, feature) -> str:
        """Get the type of a feature.

        Args:
            feature: Feature object

        Returns:
            Feature type string
        """
        if isinstance(feature, UFeature):
            return feature.feature_type
        elif hasattr(feature, 'feature_type'):
            return feature.feature_type
        elif isinstance(feature, dict):
            return feature.get('feature_type', 'unknown')
        return 'unknown'