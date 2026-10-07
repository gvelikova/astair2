"""Lazy access to matplotlib, which is only needed when plots are requested."""

import logging

logs = logging.getLogger(__name__)


def pyplot():
    """matplotlib.pyplot configured for file output, or None if matplotlib is not installed."""
    try:
        import matplotlib
    except ImportError:
        logs.error(
            "Matplotlib was not found, visualisation output is not supported. Install astair2[plot]."
        )
        return None
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.ioff()
    return plt


def plot_colors(colors, default):
    """Colours given on the command line as 'color1,color2,...' arrive as a list of characters."""
    return colors if colors == default else "".join(colors).split(",")
