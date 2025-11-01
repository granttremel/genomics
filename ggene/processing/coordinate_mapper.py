"""
Coordinate transformation utilities for genome browser.

This module handles coordinate system transformations between reference,
display, and strand-aware coordinates.
"""

from ggene.seqs.bio import reverse_complement, to_rna


class CoordinateMapper:
    """Handles coordinate transformations for genomic display."""

    def __init__(self, window_size: int, show_reverse_strand: bool = False):
        """Initialize the coordinate mapper.

        Args:
            window_size: Size of the display window
            show_reverse_strand: Whether to show reverse strand
        """
        self.window_size = window_size
        self.show_reverse_strand = show_reverse_strand

    def render_pos(self, start: int, end: int) -> tuple[int, int]:
        """Transform positions for display based on strand.

        When showing reverse strand, coordinates are flipped.

        Args:
            start: Start position in window coordinates
            end: End position in window coordinates

        Returns:
            Tuple of (display_start, display_end)
        """
        if self.show_reverse_strand:
            return self.window_size - end, self.window_size - start
        else:
            return start, end

    def map_feature_to_window(self, feature_start: int, feature_end: int,
                              window_start: int) -> tuple[int, int]:
        """Map feature coordinates to window coordinates.

        Args:
            feature_start: Feature start in genomic coordinates
            feature_end: Feature end in genomic coordinates
            window_start: Window start in genomic coordinates

        Returns:
            Tuple of (window_start, window_end) coordinates
        """
        # Convert to window-relative coordinates
        rel_start = max(0, feature_start - window_start)
        rel_end = min(self.window_size - 1, feature_end - window_start)

        # Apply strand transformation if needed
        return self.render_pos(rel_start, rel_end)

    def window_to_genomic(self, window_pos: int, window_start: int) -> int:
        """Convert window position to genomic coordinate.

        Args:
            window_pos: Position within the window (0-based)
            window_start: Genomic coordinate of window start

        Returns:
            Genomic coordinate
        """
        if self.show_reverse_strand:
            # When showing reverse strand, positions are flipped
            actual_pos = self.window_size - window_pos - 1
        else:
            actual_pos = window_pos

        return window_start + actual_pos

    def genomic_to_window(self, genomic_pos: int, window_start: int) -> int:
        """Convert genomic coordinate to window position.

        Args:
            genomic_pos: Genomic coordinate
            window_start: Genomic coordinate of window start

        Returns:
            Window position (0-based), or -1 if outside window
        """
        rel_pos = genomic_pos - window_start

        if rel_pos < 0 or rel_pos >= self.window_size:
            return -1  # Position outside window

        if self.show_reverse_strand:
            return self.window_size - rel_pos - 1
        else:
            return rel_pos

    def update_settings(self, window_size: int = None,
                        show_reverse_strand: bool = None):
        """Update mapper settings.

        Args:
            window_size: New window size (if provided)
            show_reverse_strand: New reverse strand setting (if provided)
        """
        if window_size is not None:
            self.window_size = window_size
        if show_reverse_strand is not None:
            self.show_reverse_strand = show_reverse_strand


class SequenceRenderer:
    """Handles sequence rendering transformations."""

    def __init__(self, show_reverse_strand: bool = False, show_rna: bool = False):
        """Initialize the sequence renderer.

        Args:
            show_reverse_strand: Whether to show reverse complement
            show_rna: Whether to show RNA (U instead of T)
        """
        self.show_reverse_strand = show_reverse_strand
        self.show_rna = show_rna

    def render_seq(self, seq: str) -> str:
        """Transform sequence for display.

        Args:
            seq: Input DNA sequence

        Returns:
            Transformed sequence based on settings
        """
        if self.show_reverse_strand:
            seq = reverse_complement(seq)
        if self.show_rna:
            seq = to_rna(seq)
        return seq

    def update_settings(self, show_reverse_strand: bool = None,
                        show_rna: bool = None):
        """Update renderer settings.

        Args:
            show_reverse_strand: New reverse strand setting (if provided)
            show_rna: New RNA display setting (if provided)
        """
        if show_reverse_strand is not None:
            self.show_reverse_strand = show_reverse_strand
        if show_rna is not None:
            self.show_rna = show_rna