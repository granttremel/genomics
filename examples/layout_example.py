"""
Example demonstrating the 2D layout system for genome browser.

Shows how to arrange multiple artists in rows and columns with labels.
"""

from ggene.display.renderer import ArtistRenderer
from ggene.display.layout import Label, Alignment
from ggene.display.artists.base import BaseArtist, BaseArtistParams


# Example 1: Simple row-based layout with side-by-side artists
def example_row_layout():
    """Create a layout with artists arranged horizontally."""

    # Create renderer
    renderer = ArtistRenderer(display_width=256, display_height=40)

    # Set header and footer
    renderer.set_header("Chromosome {chrom} | Position {position} | ({strand}) {nctype}")
    renderer.set_footer("[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]")

    # Add a row with 3 minimaps side-by-side
    row1 = renderer.add_row(height=6)
    row1.add_artist(
        gc_content_artist,
        width=85,
        top_label="GC Content",
        left_margin=4,
        right_margin=1
    )
    row1.add_artist(
        coverage_artist,
        width=85,
        top_label="Coverage",
        left_margin=1,
        right_margin=1
    )
    row1.add_artist(
        repeats_artist,
        width=85,
        top_label="Repeat Density",
        left_margin=1,
        right_margin=4
    )

    # Add another row with 2 wider artists
    row2 = renderer.add_row(height=6)
    row2.add_artist(correlation_artist, width=128, top_label="Autocorrelation")
    row2.add_artist(conservation_artist, width=128, top_label="Conservation")

    # Add full-width sequence display with strand labels
    renderer.add_full_width_row(
        sequence_artist,
        height=4,
        left_label=Label("5'", alignment=Alignment.RIGHT, bold=True),
        right_label=Label("3'", alignment=Alignment.LEFT, bold=True)
    )

    # Add full-width feature display
    renderer.add_full_width_row(
        feature_artist,
        height=8,
        top_label="Genomic Features"
    )

    return renderer


# Example 2: Grid-based layout with precise positioning
def example_grid_layout():
    """Create a layout using explicit grid positions."""

    renderer = ArtistRenderer(display_width=256, display_height=40)

    # Header/footer
    renderer.set_header("Chromosome {chrom} | Position {position}")
    renderer.set_footer("[Navigation: ←/→ arrows]")

    # Position artists at specific grid cells
    # Row 0, Col 0: Large overview map (spans 2 columns)
    renderer.add_artist_at(
        overview_artist,
        row=0, col=0,
        width=170, height=8,
        col_span=2,
        top_label="Chromosome Overview",
        left_margin=4
    )

    # Row 0, Col 2: Statistics panel
    renderer.add_artist_at(
        stats_artist,
        row=0, col=2,
        width=82, height=8,
        top_label="Statistics",
        left_margin=2,
        right_margin=4
    )

    # Row 1: Full-width sequence (spans all 3 columns)
    renderer.add_artist_at(
        sequence_artist,
        row=1, col=0,
        width=256, height=4,
        col_span=3,
        left_label="5'",
        right_label="3'"
    )

    # Row 2: Features
    renderer.add_artist_at(
        feature_artist,
        row=2, col=0,
        width=256, height=8,
        col_span=3,
        top_label="Annotations"
    )

    return renderer


# Example 3: Auto-width distribution
def example_auto_width():
    """Demonstrate automatic width distribution."""

    renderer = ArtistRenderer(display_width=256, display_height=40)

    # Header
    renderer.set_header("Multi-track view")

    # Add row with auto-width artists (will divide space equally)
    row = renderer.add_row(height=6)

    # Each artist gets equal width automatically
    for i, (artist, label) in enumerate([
        (track1_artist, "Track 1"),
        (track2_artist, "Track 2"),
        (track3_artist, "Track 3"),
        (track4_artist, "Track 4")
    ]):
        row.add_artist(
            artist,
            width="auto",  # Auto-distribute available width
            top_label=label,
            left_margin=1,
            right_margin=1
        )

    # Sequence row
    renderer.add_full_width_row(sequence_artist, height=3)

    return renderer


# Example 4: Complex multi-row layout (like your current view)
def example_complex_layout():
    """Recreate something like the view.txt layout."""

    renderer = ArtistRenderer(display_width=256, display_height=50)

    # Header
    renderer.set_header("Chromosome {chrom} | Position {position} | Window {window_size} | ({strand}) {nctype}")

    # Row 0: Three small minimaps side-by-side
    row0 = renderer.add_row(height=5)
    for artist, label in [
        (minimap1, "Signal 1"),
        (minimap2, "Signal 2"),
        (minimap3, "Signal 3")
    ]:
        row0.add_artist(artist, width="auto", top_label=label)

    # Row 1: Ruler
    renderer.add_full_width_row(ruler_artist, height=1)

    # Row 2: Another set of minimaps
    row2 = renderer.add_row(height=5)
    for artist, label in [
        (minimap4, "Signal 4"),
        (minimap5, "Signal 5"),
        (minimap6, "Signal 6")
    ]:
        row2.add_artist(artist, width="auto", top_label=label)

    # Row 3: Another ruler
    renderer.add_full_width_row(ruler_artist, height=1)

    # Row 4: Larger visualization
    renderer.add_full_width_row(
        large_viz_artist,
        height=6,
        top_label="Detailed Metrics"
    )

    # Row 5: Sequences (ref and alt)
    renderer.add_full_width_row(
        sequence_artist,
        height=3,
        left_label=Label("5'", bold=True),
        right_label=Label("3'", bold=True)
    )

    # Row 6: Another metrics row
    renderer.add_full_width_row(metrics_artist, height=6)

    # Row 7: Features with strand indicators
    renderer.add_full_width_row(
        feature_artist,
        height=5,
        top_label=Label("Genomic Features", alignment=Alignment.CENTER),
        bottom_label="Ruler",
        left_label="3'-",
        right_label="-5'"
    )

    # Row 8: Gene annotations
    renderer.add_full_width_row(
        gene_artist,
        height=4,
        left_label="5'-",
        right_label="-3'"
    )

    # Footer
    renderer.set_footer("[←/→: move | ↑/↓: jump | g: goto | ?: help | q: quit]")

    return renderer


# Example 5: Using labels with colors and alignment
def example_styled_labels():
    """Demonstrate label styling options."""

    renderer = ArtistRenderer(display_width=256, display_height=30)

    # Header
    renderer.set_header("Styled Layout Example")

    # Row with styled labels
    renderer.add_full_width_row(
        artist,
        height=10,
        top_label=Label(
            "Forward Strand",
            alignment=Alignment.CENTER,
            color="\x1b[38;5;82m",  # Green
            bold=True
        ),
        left_label=Label(
            "5'→",
            alignment=Alignment.RIGHT,
            bold=True
        ),
        right_label=Label(
            "→3'",
            alignment=Alignment.LEFT,
            bold=True
        ),
        bottom_label=Label(
            "Position ruler",
            alignment=Alignment.CENTER,
            dim=True
        ),
        left_margin=6,
        right_margin=6,
        top_margin=1,
        bottom_margin=1
    )

    return renderer


# Usage in browser
def setup_browser_renderer(browser_state):
    """Example of setting up renderer in actual browser."""

    renderer = ArtistRenderer(
        display_width=browser_state.window_size,
        display_height=browser_state.display_height
    )

    # Configure header
    renderer.set_header(
        "Chromosome {chrom} | Position {position} | Window {window_size} | ({strand}) {nctype}"
    )

    # Add your artists based on what you want to show
    if browser_state.show_minimaps:
        row = renderer.add_row(height=6)
        row.add_artist(gc_artist, width="auto", top_label="GC%")
        row.add_artist(cov_artist, width="auto", top_label="Coverage")
        row.add_artist(rep_artist, width="auto", top_label="Repeats")

    # Always show sequence
    renderer.add_full_width_row(
        sequence_artist,
        height=3,
        left_label="5'",
        right_label="3'"
    )

    # Features if enabled
    if browser_state.show_features:
        renderer.add_full_width_row(
            feature_artist,
            height=6,
            top_label="Features"
        )

    # Footer
    renderer.set_footer("[Navigation help]")

    return renderer


if __name__ == "__main__":
    # Print examples
    print("See examples above for different layout configurations!")
    print("\nKey features:")
    print("  - Row-based layout: renderer.add_row(height).add_artist(...)")
    print("  - Full-width rows: renderer.add_full_width_row(artist, height, ...)")
    print("  - Grid positioning: renderer.add_artist_at(artist, row, col, width, height)")
    print("  - Labels: top_label, bottom_label, left_label, right_label")
    print("  - Auto-width: width='auto' for equal distribution")
    print("  - Styled labels: Label(text, alignment, color, bold, dim)")
