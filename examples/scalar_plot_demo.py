"""
Demo script showing how to use the ScalarPlot class.
"""

import sys
import math
sys.path.insert(0, '/home/grant/genomics-prj')

from ggene.scalar_plot import ScalarPlot

def demo_basic():
    """Basic usage example."""
    print("=" * 60)
    print("BASIC PLOT")
    print("=" * 60)

    # Create some sample data
    data = [math.sin(x/10) * 5 + 10 for x in range(80)]

    # Create and show a basic plot
    plot = ScalarPlot(data)
    plot.show()
    print()


def demo_color_schemes():
    """Demonstrate different color schemes."""
    print("=" * 60)
    print("COLOR SCHEMES")
    print("=" * 60)

    data = [math.sin(x/10) * 5 + 10 for x in range(80)]
    plot = ScalarPlot(data)

    for scheme in ['vscode', 'blue', 'gray', 'icy', 'foggy']:
        print(f"\n{scheme.upper()}:")
        plot.set_color_scheme(scheme)
        plot.show()
    print()


def demo_with_ruler():
    """Demonstrate adding a ruler."""
    print("=" * 60)
    print("WITH RULER")
    print("=" * 60)

    data = [math.sin(x/10) * 5 + 10 for x in range(80)]
    plot = ScalarPlot(data, bit_depth=16)
    plot.set_color_scheme('vscode')
    plot.enable_ruler(xmin=0, xmax=1000, genomic=True)
    plot.show()
    print()


def demo_distribution():
    """Demonstrate distribution mode."""
    print("=" * 60)
    print("DISTRIBUTION MODE")
    print("=" * 60)

    dist_data = {
        'A': 10,
        'C': 15,
        'G': 12,
        'T': 13,
    }

    plot = ScalarPlot(dist_data, labels=True, edges=True)
    plot.set_color_scheme('blue')
    plot.show()
    print()


def demo_mid_mode():
    """Demonstrate mid mode (centered around zero)."""
    print("=" * 60)
    print("MID MODE (centered)")
    print("=" * 60)

    # Data that oscillates around zero
    data = [math.sin(x/5) * 3 for x in range(60)]

    plot = ScalarPlot(data, mode='mid')
    plot.set_color_scheme('dusty')
    plot.show()
    print()


def demo_method_chaining():
    """Demonstrate method chaining."""
    print("=" * 60)
    print("METHOD CHAINING")
    print("=" * 60)

    data = [math.exp(-((x-40)/10)**2) * 10 for x in range(80)]

    # Chain multiple operations
    plot = (ScalarPlot(data)
            .set_color_scheme('vscode')
            .set_bit_depth(24)
            .enable_range_labels(True, '0.1f')
            .enable_ruler(xmin=0, xmax=800, num_labels=5))

    plot.show()
    print()


def demo_rerendering():
    """Demonstrate re-rendering with different options."""
    print("=" * 60)
    print("RE-RENDERING")
    print("=" * 60)

    data = [x**0.5 * 2 for x in range(1, 81)]
    plot = ScalarPlot(data, bit_depth=8)

    print("Original (8-bit):")
    plot.show()

    print("\nIncreased to 16-bit:")
    plot.set_bit_depth(16)
    plot.show()

    print("\nWith different colors:")
    plot.set_color(fg=220, bg=232)
    plot.show()

    print("\nFlipped:")
    plot.set_flip(True)
    plot.show()
    print()


if __name__ == '__main__':
    demo_basic()
    demo_color_schemes()
    demo_with_ruler()
    demo_distribution()
    demo_mid_mode()
    demo_method_chaining()
    demo_rerendering()
