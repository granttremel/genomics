"""
Main display renderer for genome browser.

This module coordinates all display components to produce the complete
browser view.
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
from ggene import draw
import numpy as np

from ggene.display.sequence_display import SequenceDisplay
from ggene.display.feature_display import FeatureDisplay
from ggene.display.amino_acid_display import AminoAcidDisplay
from ggene.display.colors import Colors
from ggene.processing.feature_processor import FeatureProcessor
import logging

logger = logging.getLogger(__name__)

@dataclass
class RenderState:
    left_margin:int = 10
    num_rows:int = 24
    
    fg_color:int = 128
    bg_color:int = 234
    fg_color_rc: int = 122
    data_keys:List[Any] = field(default_factory = list)
    
    show_features: bool = True
    show_quality: bool = False
    show_amino_acids: bool = False
    show_reverse_strand: bool = False
    show_rna: bool = False
    show_data: bool = True
    
    arrow_disp = '↽-⇀'
    feature_disp = '|-|'
    
    
    def _get_margin(self, lbl = "", color=""):
        fill = self.left_margin - len(lbl)
        lbl_fill = lbl + " "*fill
        return f"{Colors.BOLD}{color}{lbl_fill}{Colors.RESET}"
    # plus more


class DisplayRenderer:
    """Coordinates all display operations for the genome browser."""

    def __init__(self, genome_manager, window_size: int = 240):
        """Initialize the display renderer.

        Args:
            genome_manager: GenomeManager instance
            window_size: Default window size
        """
        self.gm = genome_manager
        self.window_size = window_size

        bg, fg = draw.get_color_scheme("test")
        self.rstate = RenderState(bg_color = bg, fg_color = fg, fg_color_rc = fg - 6)
        
        # Initialize sub-renderers
        self.seq_display = SequenceDisplay(self.rstate, window_size)
        self.feature_display = FeatureDisplay(self.rstate, window_size)
        self.aa_display = AminoAcidDisplay()
        self.feature_processor = FeatureProcessor()
        
    
    def update_render_state(self, **kwargs):
        for k in kwargs:
            if hasattr(self.rstate, k):
                setattr(self.rstate, k, kwargs[k])
        return self.rstate

    def render_full_view(self, window, state, second_window=None,
                         second_state=None) -> Dict[str, List[str]]:
        all_lines = {}

        # Process features
        features = self.feature_processor.get_all_features(
            window, state, self.gm
        )

        # Render header
        header_lines = self._render_header(state, second_state)
        if header_lines:
            all_lines['header'] = header_lines

        # Render primary sequences
        if window.ref_seq or window.alt_seq:
            seq_lines, n_seqs = self.seq_display.render_sequences(
                window.ref_seq, window.alt_seq, features, state
            )
            if seq_lines:
                all_lines['sequences'] = seq_lines

        # Render secondary sequences if provided
        if second_window and second_state:
            second_features = self.feature_processor.get_all_features(
                second_window, second_state, self.gm
            )
            seq2_lines, n_seqs2 = self.seq_display.render_sequences(
                second_window.ref_seq, second_window.alt_seq,
                second_features, second_state
            )
            if seq2_lines:
                all_lines['sequences2'] = seq2_lines

        # Render features
        if self.rstate.show_features and features:
            feature_lines = self.feature_display.render_features(
                features, window.variant_features, state
            )
            if feature_lines:
                all_lines['features'] = feature_lines

        # Render amino acids
        if self.rstate.show_amino_acids and not second_state:
            aa_lines = self.aa_display.render_amino_acid_changes(
                window.ref_seq, window.alt_seq, features, state
            )
            if aa_lines:
                all_lines['amino_acids'] = aa_lines

        # Render data visualizations if available
        if hasattr(self.rstate, 'show_data') and self.rstate.show_data:
            data_lines = self._render_data_visualization(window, state)
            if data_lines:
                all_lines['data'] = data_lines

        line_count = sum([len(v)+1 for k,v in all_lines.items()])
        # Render footer
        footer_lines = self._render_footer(state, lines_used = line_count)
        if footer_lines:
            all_lines['footer'] = footer_lines

        return all_lines

    def render_interleaved(self, all_lines: Dict[str, List[str]]) -> str:
        """Interleave display sections into final output.

        Args:
            all_lines: Dictionary of section lines

        Returns:
            Complete display string
        """
        output_lines = []

        # Define section order
        section_order = [
            'header',
            'sequences',
            'sequences2',
            'data',
            'features',
            'amino_acids',
            'footer'
        ]

        # Add sections in order
        for section in section_order:
            if section in all_lines:
                output_lines.extend(all_lines[section])

        return '\n'.join(output_lines)

    def _get_margin(self, lbl = "", color=""):
        fill = self.rstate.left_margin - len(lbl)
        lbl_fill = lbl + " "*fill
        return f"{Colors.BOLD}{color}{lbl_fill}{Colors.RESET}"
        
    def _render_header(self, state, second_state=None) -> List[str]:
        """Render header information.

        Args:
            state: Primary browser state
            second_state: Optional secondary state

        Returns:
            List of header lines
        """
        lines = []

        # Title line
        title = f"Chromosome {state.chrom} | Position {state.position:,}"
        if self.rstate.show_reverse_strand:
            title += " (reverse strand)"
        if self.rstate.show_rna:
            title += " [RNA mode]"

        lines.append(f"{Colors.BOLD}{title}{Colors.RESET}")

        # Secondary position if comparing
        if second_state:
            sec_title = f"Compare: Chr{second_state.chrom} | Pos {second_state.position:,}"
            lines.append(f"{Colors.DIM}{sec_title}{Colors.RESET}")

        # Separator
        lines.append("=" * min(self.window_size, 80))

        return lines

    def _render_footer(self, state, lines_used = 24) -> List[str]:
        lines = []
        
        space = self.rstate.num_rows - lines_used
        lines.extend(["\n"]*space)
        
        # Separator
        lines.append("-" * min(self.window_size, 80))
        
        # Current gene if cached
        if hasattr(state, 'current_gene') and state.current_gene:
            gene_info = f"Gene: {state.current_gene}"
            lines.append(f"{Colors.GENE}{gene_info}{Colors.RESET}")

        # Navigation help
        nav_help = "Navigation: ← → (move) | ↑ ↓ (jump) | g (goto) | q (quit)"
        lines.append(f"{Colors.DIM}{nav_help}{Colors.RESET}")

        return lines

    def _render_data(self, data1, bg, fg, label="", minval=None, maxval=None, bit_depth = 16):
        
        data_rndr = draw.scalar_to_text_nb(data1, minval = minval, maxval =maxval, fg_color = fg, bg_color = bg, bit_depth = bit_depth)
        
        data_full = data_rndr
        num_lines = len(data_full)
        
        outlines = []
        
        for i in range(num_lines):
            
            if i == num_lines//2:
                marg_str = self._get_margin(label)
            else:
                marg_str = self._get_margin()
            
            outlines.append(f"{marg_str}{data_full[i]}")
        return outlines

    def _render_data_stacked(self, data1, data2, bg, fg1, fg2, label1="", label2="", minval=None, maxval=None, bit_depth = 16):
        
        data1_rndr = draw.scalar_to_text_nb(data1, minval = minval, maxval =maxval, fg_color = fg1, bg_color = bg, bit_depth = bit_depth)
        data2_rndr = draw.scalar_to_text_nb(data2, minval = minval, maxval = maxval, fg_color = fg2, bg_color = bg, flip = True, bit_depth = bit_depth)
        
        data_full = data1_rndr + data2_rndr
        num_lines = len(data_full)
        
        outlines = []
        
        for i in range(num_lines):
            
            if i == num_lines//2-1:
                marg_str = self._get_margin(label1)
            elif i==num_lines//2:
                marg_str = self._get_margin(label2)
            else:
                marg_str = self._get_margin()
            
            outlines.append(f"{marg_str}{data_full[i]}")
        return outlines

    def _render_stats_table(self, data, label1, data2=None, label2=""):
        outlines = []
        marg = self._get_margin(label1)
        outlines.append("{}{:<10}{:<10}{:<10}{:<10}".format(*[self._get_margin(""),"Mean", "SD", "Min", "Max"]))
        outlines.append(f"{marg}{np.mean(data):<10.2g}{np.std(data):<10.2g}{min(data):<10.2g}{max(data):<10.2g}")
        if data2 is not None:
            marg2 = self._get_margin(label2)
            outlines.append(f"{marg2:<10}{np.mean(data2):<10.2g}{np.std(data2):<10.2g}{min(data2):<10.2g}{max(data2):<10.2g}")
        return outlines

    def _render_data_visualization(self, window, state, **kwargs) -> List[str]:
        """Render data visualizations (correlation, etc.).

        Args:
            window: Window data
            state: Browser state

        Returns:
            List of data visualization lines
        """
        bit_depth = kwargs.get("bit_depth", 16)
        show_stats = kwargs.get("show_stats", True)
        
        lines = []
        stacked = kwargs.get("stacked", True) # stack alternating plots
        
        data_keys = self.rstate.data_keys
        labels = [d + ":" for d in data_keys]
        
        bg, fg, fg2 = self.rstate.bg_color, self.rstate.fg_color, self.rstate.fg_color_rc
        
        if stacked:
            for i in range(0, len(data_keys), 2):
                dk1 = data_keys[i]
                dk2 = data_keys[i+1]
                d1 = state._data_cache.get(dk1)
                d2 = state._data_cache.get(dk2)
                newlines = self._render_data_stacked(d1, d2, bg, fg, fg2, label1 = labels[i], label2 = labels[i+1], bit_depth = bit_depth)
                lines.extend(newlines)
                if show_stats:
                    statlines = self._render_stats_table(d1, labels[i], data2=d2, label2 = labels[i+1])
                    newlines.extend(statlines)
                
        else:
            for i in range(len(data_keys)):
                dk1 = data_keys[i]
                d1 = state._data_cache.get(dk1)
                newlines = self._render_data(d1, bg, fg, label = labels[i], bit_depth = bit_depth)
                lines.extend(newlines)
                if show_stats:
                    statlines = self._render_stats_table(d1, labels[i])
                    newlines.extend(statlines)

        return lines

    def render_position_info(self, state) -> str:
        """Render current position information.

        Args:
            state: Browser state

        Returns:
            Position information string
        """
        info_parts = []

        # Basic position
        info_parts.append(f"Chr{state.chrom}:{state.position:,}")

        # Window size
        info_parts.append(f"Window: {state.window_size}bp")

        # Strand
        strand = "reverse" if self.rstate.show_reverse_strand else "forward"
        info_parts.append(f"Strand: {strand}")

        # Current gene if available
        if hasattr(state, 'current_gene') and state.current_gene:
            info_parts.append(f"Gene: {state.current_gene}")

        return " | ".join(info_parts)

    def render_help(self) -> str:
        """Render help text.

        Returns:
            Help text string
        """
        help_text = """
Genome Browser Help
==================

Navigation:
  → / l        Move forward (small step)
  ← / h        Move backward (small step)
  L            Move forward (large step)
  H            Move backward (large step)
  ↑ / k        Jump to next gene
  ↓ / j        Jump to previous gene
  g            Go to position
  c            Switch chromosome

Display:
  f            Toggle features
  a            Toggle amino acids
  r            Toggle RNA mode
  R            Toggle reverse strand
  d            Toggle data visualization
  +/-          Zoom in/out

Other:
  s            Save state
  S            Load state
  q            Quit
  ?            Show this help

Feature Types:
  Gene         Blue [====]
  Transcript   Magenta [----]
  Exon         Cyan |████|
  CDS          Yellow {####}
  UTR          Gray <---->

Variant Types:
  SNP          Red
  Insertion    Green
  Deletion     Yellow
"""
        return help_text

    def update_window_size(self, new_size: int):
        """Update window size for all components.

        Args:
            new_size: New window size
        """
        self.window_size = new_size
        self.seq_display.window_size = new_size
        self.feature_display.window_size = new_size