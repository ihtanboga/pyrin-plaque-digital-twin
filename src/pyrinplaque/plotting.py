"""Reproducible plotting helpers for Day-3 cell-state figures. Publication-grade defaults."""
import matplotlib as mpl
import matplotlib.pyplot as plt

def apply_style(sizes=(8, 7, 6)):
    base, mid, small = sizes
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": mid, "xtick.labelsize": small, "ytick.labelsize": small,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.titlelocation": "left", "legend.frameon": False,
    })

def set_frame(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
