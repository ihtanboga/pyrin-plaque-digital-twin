"""Shared publication style for revision figures (Elsevier column widths, 600-dpi TIFF + vector PDF).

Palette: categorical slots validated with the dataviz validator (adjacent CVD dE >= 9.1,
normal-vision dE >= 19.6 on a white surface). Identity is never colour-alone: panels
carry direct labels or legends.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
WIDTH = {"single": 90 * MM, "onehalf": 140 * MM, "double": 190 * MM}
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "magenta": "#e87ba4", "green": "#008300", "violet": "#4a3aa7", "red": "#e34948",
     "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#8a8984", "grid": "#e4e3df", "gray": "#9a9a96",
     "midgray": "#f0efec"}
# semantic assignments used across figures
MODULE = {"PYRIN_BACKBONE": C["blue"], "NLRP3_BACKBONE": C["orange"],
          "EFFECTOR_CYTOKINE_ARM": C["aqua"], "EFFECTOR_LYSIS_ARM": C["magenta"]}
LINEAGE = {"Monocyte": C["orange"], "Macrophage": C["blue"], "Dendritic cell": C["aqua"],
           "Neutrophil": C["ink2"], "Mixed": "#c4c3be"}
DIVERGING = ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f4b3b2", "#e66767", "#a83232"]


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": C["ink2"],
        "axes.labelcolor": C["ink"], "xtick.color": C["ink2"], "ytick.color": C["ink2"],
        "axes.titlelocation": "left", "axes.titleweight": "bold", "legend.frameon": False,
        "savefig.dpi": 600, "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.grid": False, "lines.linewidth": 1.0, "patch.linewidth": 0,
    })


def panel_label(ax, letter, dx=-0.12, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right",
            color=C["ink"])


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"{name}.tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
