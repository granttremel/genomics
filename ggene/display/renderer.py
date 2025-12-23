"""
Main display renderer for genome browser.

This module coordinates all display components to produce the complete
browser view with 2D layout capabilities.
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any, Union, TYPE_CHECKING
from ggene.draw.colors import Colors

from ggene.display.colors import FColors

from ggene.display.artists import BaseArtist, BaseArtistParams
from ggene.display.layout import GridLayout, RowBuilder, FullWidthRow, Panel, Label, Alignment, VAlignment

if TYPE_CHECKING:
    from ggene.browser.genome_browser import BrowserState
    from ggene.database.genome_iterator import GenomeWindow

import logging

logger = logging.getLogger(__name__)

@dataclass
class RenderParams:
    display_width:int = 256
    display_height:int = 40
    left_margin:int = 4
    right_margin:int = 4
    row_marker: Optional[str] = None
    
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
    
    def _get_margin(self, lbl = "", color=""):
        fill = self.left_margin - len(lbl)
        lbl_fill = lbl + " "*fill
        return f"{Colors.BOLD}{color}{lbl_fill}{Colors.RESET}"
    
    # plus more


class ArtistRenderer:
    """
    2D layout renderer for genome browser.

    Arranges artists in a flexible grid layout with support for:
    - Multi-row, multi-column layouts
    - Labels (top/bottom/left/right)
    - Margins and spacing
    - Headers and footers

    Example usage:
        # Create renderer
        renderer = ArtistRenderer(display_width=256, display_height=40)

        # Add header
        renderer.set_header("Chromosome {chrom} | Position {position}")

        # Add a row with multiple artists side-by-side
        row1 = renderer.add_row(height=6)
        row1.add_artist(minimap1, width=64, top_label="GC Content")
        row1.add_artist(minimap2, width=64, top_label="Coverage")
        row1.add_artist(minimap3, width=64, top_label="Repeats")

        # Add full-width artist
        renderer.add_full_width_row(sequence_artist, height=4,
                                     left_label="5'", right_label="3'")

        # Add footer
        renderer.set_footer("[Navigation help text]")

        # Render
        lines = renderer.render(state, window)
    """

    def __init__(self, display_width: int = 256, display_height: int = 40, **kwargs):
        """
        Initialize the artist renderer.

        Args:
            display_width: Total terminal width
            display_height: Total terminal height
            **kwargs: Additional render state parameters
        """
        self.display_width = display_width
        self.display_height = display_height

        # Render state for compatibility
        rparams = RenderParams(
            display_width=display_width,
            display_height=display_height,
            **kwargs
        )
        self.params = rparams

        # Grid layout system
        self.layout = GridLayout(display_width, display_height, renderer = self)

        # Header and footer
        self._header_template: Optional[str] = None
        self._footer_template: Optional[str] = None
        self._header_height = 2  # Lines reserved for header
        self._footer_height = 2  # Lines reserved for footer

        # Current row index for sequential adding
        self._next_row = 0

        # Track artists for backward compatibility
        self.artists: Dict[str, BaseArtist] = {}
        self.artist_order: List[str] = []

    def set_header(self, template: str, height: int = 2):
        """
        Set header template.

        Template can include {chrom}, {position}, {window_size}, etc.
        which will be filled from state during render.

        Args:
            template: Header text template
            height: Number of lines for header
        """
        self._header_template = template
        self._header_height = height

    def set_footer(self, template: str, height: int = 2):
        """
        Set footer template.

        Args:
            template: Footer text template
            height: Number of lines for footer
        """
        self._footer_template = template
        self._footer_height = height

    def add_row(self, height: int, fixed_height: bool = False,
                valign: VAlignment = VAlignment.TOP) -> RowBuilder:
        """
        Add a new row for horizontal artist placement.

        Args:
            height: Height of the row in lines
            fixed_height: If True, panels are padded/trimmed to exact height
            valign: Vertical alignment when fixed_height is True (TOP, CENTER, BOTTOM)

        Returns:
            RowBuilder for adding artists to this row
        """
        row_builder = RowBuilder(self.layout, self._next_row, height,
                                 fixed_height=fixed_height, valign=valign)
        self._next_row += 1
        return row_builder

    def add_full_width_row(self, artist: 'BaseArtist', height: int,
                           fixed_height: bool = False,
                           valign: VAlignment = VAlignment.TOP,
                           **label_kwargs) -> Panel:
        """
        Add a full-width row with a single artist.

        Args:
            artist: Artist to add
            height: Height of the row
            fixed_height: If True, panel is padded/trimmed to exact height
            valign: Vertical alignment when fixed_height is True
            **label_kwargs: Label parameters (top_label, bottom_label, etc.)

        Returns:
            The created Panel
        """
        FullWidthRow(self.layout, self._next_row, height, artist,
                     fixed_height=fixed_height, valign=valign, **label_kwargs)
        self._next_row += 1

        # Track artist
        self.artists[artist.name] = artist
        self.artist_order.append(artist.name)

        return self.layout.cells[-1].panel

    def add_artist_at(self, artist: 'BaseArtist', row: int, col: int,
                     width: int, height: int,
                     row_span: int = 1, col_span: int = 1,
                     **panel_kwargs) -> Panel:
        """
        Add an artist at a specific grid position.

        Args:
            artist: Artist to add
            row: Row index
            col: Column index
            width: Panel width
            height: Panel height
            row_span: Rows to span
            col_span: Columns to span
            **panel_kwargs: Additional Panel parameters

        Returns:
            The created Panel
        """
        panel = Panel(artist=artist, width=width, height=height, **panel_kwargs)
        self.layout.add_cell(panel, row, col, row_span, col_span)

        # Track artist
        # name = f"{type(artist).__name__}{len(self.artists)}"
        self.artists[artist.name] = artist
        self.artist_order.append(artist.name)

        return panel

    def add_artist(self, artist: 'BaseArtist') -> Panel:
        """
        Add an artist (legacy API - adds as full-width row).

        Args:
            artist: Artist to add
            artist_name: Optional name for artist
            height: Height for the artist's row

        Returns:
            The created Panel
        """
        # if not artist_name:
            # artist_name = f"{type(artist).__name__}{len(self.artists)}"
        artist_name = artist.name

        # panel = self.add_full_width_row(artist, height)

        self.artists[artist_name] = artist
        if artist_name not in self.artist_order:
            self.artist_order.append(artist_name)

        # Set margins for compatibility
        artist.set_margins(self.params.left_margin, self.params.right_margin)

        # return panel

    def update(self, **kwargs):
        rsd = self.params.__dict__.copy()
        for k, v in kwargs.items():
            if hasattr(self.params, k):
                rsd[k] = v    
        
        new_rparams = RenderParams(**rsd)
        self.params = new_rparams

    def update_artist(self, artist_name: Optional[str] = None,
              artist_index: Optional[int] = None, **kwargs):
        """Update an artist's parameters."""
        artist = self.get_artist(index=artist_index, name=artist_name)
        if artist:
            artist.update(**kwargs)

    def update_all(self, **options):
        
        for an, artist in self.artists.items():
            artist.update(**options)
            
    def update_by_type(self, artist_type = "LineArtist", **options):
        
        for artist in self.artists.values():
            if type(artist).__name__ == artist_type:
                artist.update(**options)
    
    def iter_artists(self, artist_type = None):
        
        for nm, artist in self.artists.items():
            
            if artist_type and not type(artist).__name__ == artist_type:
                continue
            
            yield artist
    
    def render(self, state: 'BrowserState', window: 'GenomeWindow',
              footer_text: str = "") -> List[str]:
        """
        Render the complete view.

        Args:
            state: Browser state
            window: Genome window
            footer_text: Footer text (overrides template if provided)

        Returns:
            List of lines for complete terminal view
        """
        all_lines = []
        row_border = None
        if self.params.row_marker:
            row_border = self.get_row_border(self.params.row_marker)

        # Reserve space for header/footer in layout
        content_height = (self.display_height -
                         self._header_height - self._footer_height)

        # Update layout total height
        old_height = self.layout.total_height
        self.layout.total_height = content_height

        # Render header
        if self._header_template:
            header_lines = self._render_header_frm(state)
        else:
            header_lines = self._render_header(state)
        all_lines.extend(header_lines)

        # Render grid content
        content_lines = self.layout.render(state, window, row_border = row_border)
        all_lines.extend(content_lines)

        # Render footer
        footer_lines = self._render_footer(state, footer_text)
        all_lines.extend(footer_lines)

        # Restore layout height
        self.layout.total_height = old_height

        return all_lines

    def get_row_border(self, row_border_char):
        
        if not row_border_char:
            return
        
        border_len = self.display_width - self.params.left_margin - self.params.right_margin
        
        num_space = 8
        row_border = " "*num_space + row_border_char * (border_len - 2*num_space) + " "*num_space
        row_border = FColors.DIM + row_border.center(len(row_border) + self.params.left_margin + self.params.right_margin) + FColors.RESET
        
        return row_border

    def _render_header(self, state: 'BrowserState') -> List[str]:
        """Render header with template substitution."""
        return BaseArtist.make_header(state.chrom, state.position, state.window_size, state.show_reverse_strand, state.show_rna)
    
    def _render_header_frm(self, state):
        lines = []

        # Format template
        header_text = self._header_template.format(
            chrom=state.chrom,
            position=state.position,
            window_size=getattr(state, 'window_size', 0),
            strand='-' if getattr(state, 'show_reverse_strand', False) else '+',
            nctype='RNA' if getattr(state, 'show_rna', False) else 'DNA'
        )

        lines.append(f"\x1b[1m{header_text}\x1b[0m")

        # Separator
        if self._header_height > 1:
            lines.append("-" * min(self.display_width, 80))

        # Pad to header height
        while len(lines) < self._header_height:
            lines.append("")

        return lines

    def _render_footer(self, state: 'BrowserState',
                      footer_text: str = "") -> List[str]:
        """Render footer."""
        
        lines = []
        
        # Separator
        lines.append("-" * min(self.display_width, 80))
        
        # Footer text
        if footer_text:
            formatted = FColors(effect="dim").format_string(footer_text)
            lines.append(formatted)
        elif self._footer_template:
            lines.append(FColors(effect="dim").format_string(self._footer_template))

        # Pad to footer height
        while len(lines) < self._footer_height-1:
            lines.append("+")
        return lines

    def get_artist(self, index: Optional[int] = None,
                  name: Optional[str] = None) -> Optional['BaseArtist']:
        """Get an artist by name or index."""
        if index is not None:
            if 0 <= index < len(self.artist_order):
                name = self.artist_order[index]
                return self.artists.get(name)
        elif name:
            return self.artists.get(name)
        return None

    def __getitem__(self, ind: Union[str, int]) -> Optional['BaseArtist']:
        """Get artist by index or name."""
        if isinstance(ind, str):
            return self.get_artist(name=ind)
        elif isinstance(ind, int):
            return self.get_artist(index=ind)
        return None
