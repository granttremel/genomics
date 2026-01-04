"""
2D layout system for arranging artists in terminal space.

Provides a grid-based layout manager for organizing multiple artists
into rows and columns with labels, margins, and borders.
"""

from typing import List, Dict, Optional, Union, Tuple, TYPE_CHECKING
from dataclasses import dataclass, field
from enum import Enum
import logging
import time

from ggene.display.colors import FColors

if TYPE_CHECKING:
    from ggene.display.artists.base import BaseArtist
    from ggene.browser.genome_browser import BrowserState
    from ggene.database.genome_iterator import GenomeWindow

logger = logging.getLogger(__name__)


class Alignment(Enum):
    """Label alignment options."""
    LEFT = "left"
    CENTER = "center"
    RIGHT = "right"


class VAlignment(Enum):
    """Vertical alignment options for panel content."""
    TOP = "top"
    CENTER = "center"
    BOTTOM = "bottom"


@dataclass
class Label:
    """A text label with positioning and formatting."""
    text: str
    alignment: Alignment = Alignment.LEFT
    color: str = ""
    bold: bool = False
    dim: bool = False

    def format(self, width: int) -> str:
        """Format label to fit within width."""
        text = self.text[:width]  # Truncate if too long

        # Apply color/formatting
        formatted = text
        if self.color:
            formatted = f"{self.color}{formatted}\x1b[0m"
        if self.bold:
            formatted = f"\x1b[1m{formatted}\x1b[0m"
        if self.dim:
            formatted = f"\x1b[2m{formatted}\x1b[0m"

        # Apply alignment
        if self.alignment == Alignment.LEFT:
            return formatted.ljust(width + (len(formatted) - len(text)))
        elif self.alignment == Alignment.CENTER:
            padding = (width - len(text)) // 2
            return " " * padding + formatted + " " * (width - len(text) - padding)
        else:  # RIGHT
            return formatted.rjust(width + (len(formatted) - len(text)))


@dataclass
class Panel:
    """
    A panel wraps an artist with labels and spacing.

    Handles all the formatting around an artist including:
    - Top/bottom labels
    - Left/right margins with optional labels
    - Borders (optional)
    - Padding
    - Vertical alignment within fixed height
    """
    artist: 'BaseArtist'
    width: int  # Width in characters
    height: int  # Height in lines
    enabled:bool = True

    # Labels
    top_label: Optional[Union[str, Label]] = None
    bottom_label: Optional[Union[str, Label]] = None
    left_label: Optional[Union[str, Label]] = None
    right_label: Optional[Union[str, Label]] = None

    # Margins (in characters)
    left_margin: int = 4
    right_margin: int = 0
    top_margin: int = 0
    bottom_margin: int = 0

    # Border
    show_border: bool = False
    border_char: str = "│"

    # Vertical alignment within fixed height
    valign: VAlignment = VAlignment.TOP
    fixed_height: bool = False  # If True, pad/trim to exact height

    def __post_init__(self):
        """Convert string labels to Label objects."""
        for attr in ['top_label', 'bottom_label', 'left_label', 'right_label']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Label(value))

    @property
    def content_width(self) -> int:
        """Width available for artist content."""
        return self.width - self.left_margin - self.right_margin

    @property
    def content_height(self) -> int:
        """Height available for artist content."""
        margins = self.top_margin + self.bottom_margin
        labels = 0
        if self.top_label:
            labels += 1
        if self.bottom_label:
            labels += 1
        return max(1, self.height - margins - labels)

    def render(self, state: 'BrowserState', window: 'GenomeWindow') -> List[str]:
        """
        Render the panel with artist content and all labels/margins.

        Returns list of lines (strings) representing the full panel.
        """
        lines = []
        
        if not self.enabled:
            return lines

        # Top margin
        for _ in range(self.top_margin):
            lines.append(" " * self.width)

        # Top label
        if self.top_label:
            lines.append(self._format_horizontal_label(self.top_label))

        # Get artist content
        artist_lines = self.artist.render(state, window)

        # Format content lines with margins
        formatted_content = []
        for i, line in enumerate(artist_lines):
            formatted_line = self._format_content_line(line, i, len(artist_lines))
            formatted_content.append(formatted_line)

        # Apply vertical alignment if fixed_height is enabled
        if self.fixed_height:
            formatted_content = self._apply_valign(formatted_content)

        lines.extend(formatted_content)

        # Bottom label
        if self.bottom_label:
            lines.append(self._format_horizontal_label(self.bottom_label))

        # Bottom margin
        for _ in range(self.bottom_margin):
            lines.append(" " * self.width)

        return lines

    def _apply_valign(self, content_lines: List[str]) -> List[str]:
        """Apply vertical alignment to content within fixed height."""
        target_height = self.content_height
        current_height = len(content_lines)

        if current_height >= target_height:
            # Trim if needed
            return content_lines[:target_height]

        # Calculate padding
        padding_total = target_height - current_height
        empty_line = " " * self.width

        if self.valign == VAlignment.TOP:
            # Content at top, padding at bottom
            return content_lines + [empty_line] * padding_total
        elif self.valign == VAlignment.BOTTOM:
            # Padding at top, content at bottom
            return [empty_line] * padding_total + content_lines
        else:  # CENTER
            # Split padding between top and bottom
            padding_top = padding_total // 2
            padding_bottom = padding_total - padding_top
            return [empty_line] * padding_top + content_lines + [empty_line] * padding_bottom

    def _fit_content_height(self, lines: List[str]) -> List[str]:
        """Fit content to content_height by padding or trimming."""
        target = self.content_height
        current = len(lines)

        if current < target:
            # Pad with empty lines
            padding = [" " * self.content_width for _ in range(target - current)]
            return lines + padding
        elif current > target:
            # Trim (keep first content_height lines)
            return lines[:target]
        return lines

    def _format_content_line(self, line: str, line_idx: int, total_lines: int) -> str:
        """Format a single content line with margins and labels."""
        # Strip ANSI and measure actual visible length
        visible_len = self._visible_length(line)

        # Pad or trim to content_width
        # if visible_len < self.content_width:
        #     line = line + " " * (self.content_width - visible_len)
        # elif visible_len > self.content_width:
        #     # Trim (crude, doesn't account for ANSI mid-sequence)
        #     line = line[:self.content_width]

        # Left margin/label
        left = " " * self.left_margin
        if self.left_label:
            # Show label at middle line
            if line_idx == total_lines // 2:
                left = self.left_label.format(self.left_margin)
            else:
                left = " " * self.left_margin

        # Right margin/label
        right = " " * self.right_margin
        if self.right_label:
            if line_idx == total_lines // 2:
                right = self.right_label.format(self.right_margin)
            else:
                right = " " * self.right_margin

        return f"{left}{line}{right}"

    def _format_horizontal_label(self, label: Label) -> str:
        """Format a top or bottom label with margins."""
        # Label should be formatted to content width, then padded with margins
        left = " " * self.left_margin
        right = " " * self.right_margin
        label_text = label.format(self.content_width)
        return f"{FColors.BOLD}{left}{label_text}{right}{FColors.RESET}"

    def _visible_length(self, text: str) -> int:
        """Calculate visible length of text (excluding ANSI codes)."""
        import re
        ansi_escape = re.compile(r'\x1b\[[0-9;]*m')
        return len(ansi_escape.sub('', text))


@dataclass
class GridCell:
    """A cell in the grid layout."""
    row: int
    col: int
    row_span: int = 1
    col_span: int = 1
    panel: Optional[Panel] = None


class GridLayout:
    """
    2D grid layout manager.

    Divides terminal space into rows and columns, placing panels in cells.
    Automatically calculates dimensions based on available space.

    Column widths are row-specific, allowing different column configurations
    per row (e.g., 3 equal columns in row 0, full width in row 1).
    """

    def __init__(self, total_width: int, total_height: int, renderer = None):
        """
        Initialize grid layout.

        Args:
            total_width: Total terminal width
            total_height: Total terminal height
        """
        self.total_width = total_width
        self.total_height = total_height
        self._renderer = renderer

        # Grid structure
        self.cells: List[GridCell] = []
        self.row_heights: Dict[int, int] = {}  # row_idx -> height

        # Row-specific column widths: row_idx -> {col_idx -> width}
        self.row_col_widths: Dict[int, Dict[int, int]] = {}

        # Row-specific auto columns: row_idx -> [col_indices]
        self.row_auto_cols: Dict[int, List[int]] = {}

        # Auto-sizing for rows
        self.auto_rows: List[int] = []  # Rows with auto height

    def add_cell(self, panel: Panel, row: int, col: int,
                row_span: int = 1, col_span: int = 1) -> GridCell:
        """
        Add a panel to the grid.

        Args:
            panel: Panel to add
            row: Row index (0-based)
            col: Column index (0-based)
            row_span: Number of rows to span
            col_span: Number of columns to span

        Returns:
            The created GridCell
        """
        cell = GridCell(row=row, col=col, row_span=row_span,
                       col_span=col_span, panel=panel)
        self.cells.append(cell)

        # Initialize row structures if needed
        if row not in self.row_col_widths:
            self.row_col_widths[row] = {}
        if row not in self.row_auto_cols:
            self.row_auto_cols[row] = []

        if self._renderer:
            self._renderer.add_artist(panel.artist)

        return cell

    def set_row_height(self, row: int, height: int):
        """Set fixed height for a row."""
        self.row_heights[row] = height
        if row in self.auto_rows:
            self.auto_rows.remove(row)

    def set_col_width(self, row: int, col: int, width: int):
        """Set fixed width for a column in a specific row."""
        if row not in self.row_col_widths:
            self.row_col_widths[row] = {}
        self.row_col_widths[row][col] = width

        # Remove from auto cols if present
        if row in self.row_auto_cols and col in self.row_auto_cols[row]:
            self.row_auto_cols[row].remove(col)

    def set_auto_row(self, row: int):
        """Mark row for automatic height distribution."""
        if row not in self.auto_rows:
            self.auto_rows.append(row)

    def set_auto_col(self, row: int, col: int):
        """Mark column for automatic width distribution within a row."""
        if row not in self.row_auto_cols:
            self.row_auto_cols[row] = []
        if col not in self.row_auto_cols[row]:
            self.row_auto_cols[row].append(col)

    def calculate_layout(self):
        """Calculate final dimensions for all rows and columns."""
        # Get max row index
        max_row = max((cell.row + cell.row_span - 1 for cell in self.cells), default=0)

        # Calculate row heights
        fixed_height = sum(self.row_heights.values())
        remaining_height = self.total_height - fixed_height
        auto_row_height = remaining_height // max(len(self.auto_rows), 1)

        for row in range(max_row + 1):
            if row not in self.row_heights:
                self.row_heights[row] = auto_row_height

        # Calculate column widths for each row independently
        for row in range(max_row + 1):
            if row not in self.row_col_widths:
                self.row_col_widths[row] = {}

            # Get fixed widths for this row
            fixed_width = sum(self.row_col_widths[row].values())
            remaining_width = self.total_width - fixed_width

            # Get auto columns for this row
            auto_cols = self.row_auto_cols.get(row, [])
            if auto_cols:
                auto_col_width = remaining_width // len(auto_cols)
                for col in auto_cols:
                    self.row_col_widths[row][col] = auto_col_width

        # Update panel dimensions
        for cell in self.cells:
            if cell.panel:
                # Get width from row-specific column widths
                width = sum(self.row_col_widths.get(cell.row, {}).get(cell.col + i, 0)
                           for i in range(cell.col_span))
                height = sum(self.row_heights.get(cell.row + i, 0)
                            for i in range(cell.row_span))

                cell.panel.width = width
                cell.panel.height = height

    def render(self, state: 'BrowserState', window: 'GenomeWindow', row_border = None) -> List[str]:
        """
        Render the entire grid.

        Returns list of lines representing the full terminal view.
        """
        # Ensure layout is calculated and artist params updated
        self.calculate_layout()
        self._update_artist_params()

        # Organize cells by row
        rows_data: Dict[int, List[GridCell]] = {}
        for cell in self.cells:
            if cell.row not in rows_data:
                rows_data[cell.row] = []
            rows_data[cell.row].append(cell)

        # Sort cells in each row by column
        for row_cells in rows_data.values():
            row_cells.sort(key=lambda c: c.col)

        # Render row by row
        all_lines = []
        max_row = max(self.row_heights.keys(), default=0)

        for row_idx in range(max_row + 1):
            if row_idx not in rows_data:
                # Empty row
                row_height = self.row_heights.get(row_idx, 0)
                for _ in range(row_height):
                    all_lines.append(" " * self.total_width)
                continue

            # Render all cells in this row
            row_cells = rows_data[row_idx]
            cells_rendered = []
            cell_widths = []
            for cell in row_cells:
                if cell.panel:
                    t0 = time.perf_counter()
                    rendered = cell.panel.render(state, window)
                    dt = time.perf_counter() - t0
                    logger.debug(f"spent {1000*dt:0.1f}ms on panel with artist {cell.panel.artist.name}")
                    cells_rendered.append(rendered)
                    cell_widths.append(cell.panel.width)

            # Get max height of row (from panels or rendered content)
            row_height = self.row_heights.get(row_idx, 0)
            max_rendered_height = max((len(lines) for lines in cells_rendered), default=0)
            row_height = max(row_height, max_rendered_height)

            # Combine cells horizontally for each line
            for line_idx in range(row_height):
                line_parts = []
                for i, cell_lines in enumerate(cells_rendered):
                    cell_width = cell_widths[i] if i < len(cell_widths) else 0
                    if line_idx < len(cell_lines):
                        line = cell_lines[line_idx]
                        # Pad line to cell width for alignment
                        line_parts.append(self._pad_to_width(line, cell_width))
                    else:
                        # Pad with spaces
                        line_parts.append(" " * cell_width)

                all_lines.append("".join(line_parts))
            if row_border and row_idx < max_row:
                all_lines.append(row_border)

        return all_lines



    def _update_artist_params(self):
        """Update artist parameters with calculated dimensions."""
        for cell in self.cells:
            if cell.panel and cell.panel.artist:
                artist = cell.panel.artist
                # Update artist's display width/height from panel content dimensions
                content_width = cell.panel.content_width
                content_height = cell.panel.content_height
                if hasattr(artist, 'params'):
                    if hasattr(artist.params, 'display_width'):
                        artist.params.display_width = content_width
                    if hasattr(artist.params, 'display_height'):
                        artist.params.display_height = content_height
                if hasattr(artist, 'set_params'):
                    artist.set_params(display_width=content_width, display_height=content_height)

    def _pad_to_width(self, text: str, target_width: int) -> str:
        """Pad text to target width, accounting for ANSI codes."""
        visible_len = self._visible_length(text)
        if visible_len < target_width:
            return text + " " * (target_width - visible_len)
        return text

    def _visible_length(self, text: str) -> int:
        """Calculate visible length of text (excluding ANSI codes)."""
        import re
        ansi_escape = re.compile(r'\x1b\[[0-9;]*m')
        return len(ansi_escape.sub('', text))


class RowBuilder:
    """
    Helper for building rows with artists placed horizontally.

    Simpler API than GridLayout for common case of horizontal stacking.
    """

    def __init__(self, layout: GridLayout, row_idx: int, height: int,
                 fixed_height: bool = False, valign: VAlignment = VAlignment.TOP):
        """
        Initialize row builder.

        Args:
            layout: Parent GridLayout
            row_idx: Row index in the grid
            height: Height for this row
            fixed_height: If True, panels use fixed height with padding
            valign: Vertical alignment when fixed_height is True
        """
        self.layout = layout
        self.row_idx = row_idx
        self.height = height
        self.next_col = 0
        self.fixed_height = fixed_height
        self.valign = valign

        layout.set_row_height(row_idx, height)

    def add_artist(self, artist: 'BaseArtist',
                   width: Union[int, str] = "auto",
                   **label_kwargs) -> Panel:
        """
        Add an artist to this row.

        Args:
            artist: Artist to add
            width: Width in characters, or "auto" for equal distribution
            **label_kwargs: Label arguments (top_label, bottom_label, etc.)

        Returns:
            The created Panel
        """
        # Handle width - use row-specific column width
        if width == "auto":
            self.layout.set_auto_col(self.row_idx, self.next_col)
            actual_width = 0  # Will be calculated later
        else:
            actual_width = width
            self.layout.set_col_width(self.row_idx, self.next_col, width)

        # Create panel with vertical alignment settings
        panel = Panel(
            artist=artist,
            width=actual_width,
            height=self.height,
            fixed_height=self.fixed_height,
            valign=self.valign,
            **label_kwargs
        )

        # Add to grid
        self.layout.add_cell(panel, self.row_idx, self.next_col)

        self.next_col += 1
        return panel


class FullWidthRow:
    """Helper for creating a full-width row with a single artist."""

    def __init__(self, layout: GridLayout, row_idx: int, height: int,
                 artist: 'BaseArtist', fixed_height: bool = False,
                 valign: VAlignment = VAlignment.TOP, **label_kwargs):
        """Create a full-width row."""
        self.layout = layout
        self.row_idx = row_idx

        # Create full-width panel
        panel = Panel(
            artist=artist,
            width=layout.total_width,
            height=height,
            fixed_height=fixed_height,
            valign=valign,
            **label_kwargs
        )

        layout.set_row_height(row_idx, height)
        layout.set_col_width(row_idx, 0, layout.total_width)  # Row-specific
        layout.add_cell(panel, row_idx, 0)
