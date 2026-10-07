#!/usr/bin/env python
"""Revised Figure 1 (schematic, drawn programmatically): sensor-proximal backbones, the shared
caspase-1 node and the separable cytokine and lysis arms; analytic workflow with inferential roles.

Usage (repo root): python -I scripts/22_fig_schematic.py
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, "src")
from pyrinplaque import figstyle as F

OUT = "results/revision/figures"


def box(ax, x, y, w, h, title, body=None, ec=F.C["ink2"], fc="white", tsize=6.6, bsize=5.8, tcolor=None, lw=0.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012", ec=ec, fc=fc, lw=lw))
    if body:
        ax.text(x + w / 2, y + h * 0.76, title, ha="center", va="center", fontsize=tsize, fontweight="bold",
                color=tcolor or F.C["ink"])
        ax.text(x + w / 2, y + h * 0.36, body, ha="center", va="center", fontsize=bsize, color=F.C["ink2"], linespacing=1.25)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=tsize, fontweight="bold",
                color=tcolor or F.C["ink"], linespacing=1.2)


def arrow(ax, p0, p1, color=F.C["ink2"], lw=0.9, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=7, color=color, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=1, shrinkB=1))


def main():
    F.apply()
    fig = plt.figure(figsize=(F.WIDTH["double"], 150 * F.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.015, 0.975, "A", fontsize=10, fontweight="bold", va="top")
    ax.text(0.045, 0.972, "Sensor-proximal programs, the shared caspase-1 node and separable effector arms",
            fontsize=7.5, va="top", color=F.C["ink"])
    # sensors
    box(ax, 0.04, 0.795, 0.40, 0.125, "PYRIN regulatory backbone (9 genes)",
        "MEFV · RHOA · RAC1 · CDC42 · PKN1\nPKN2 · YWHAB · YWHAZ · PSTPIP1\n"
        "RhoA–PKN phosphorylation and 14-3-3 binding restrain pyrin", ec=F.C["blue"], tcolor=F.C["blue"], lw=1.2)
    box(ax, 0.56, 0.795, 0.40, 0.125, "NLRP3 comparator backbone (5 genes)",
        "NLRP3 · NEK7 · TXNIP · NFKB1 · RELA\n\nNEK7 licensing; TXNIP/NF-κB priming context",
        ec=F.C["orange"], tcolor=F.C["orange"], lw=1.2)
    box(ax, 0.41, 0.665, 0.18, 0.055, "ASC (PYCARD)")
    box(ax, 0.41, 0.565, 0.18, 0.055, "Caspase-1 (CASP1)")
    arrow(ax, (0.24, 0.795), (0.45, 0.722))
    arrow(ax, (0.76, 0.795), (0.55, 0.722))
    arrow(ax, (0.50, 0.665), (0.50, 0.622))
    ax.text(0.50, 0.765, "shared downstream genes (PYCARD, CASP1, GSDMD, IL1B, IL18)\nexcluded from discriminative backbone scores",
            ha="center", va="center", fontsize=5.6, style="italic", color=F.C["ink2"])
    # arms
    box(ax, 0.05, 0.395, 0.38, 0.12, "Cytokine arm (CASP1, IL1B, IL18)",
        "caspase-1 matures pro-IL-1β (priming-inducible) and\npro-IL-18 (constitutive); release can occur from\nliving cells through GSDMD pores",
        ec=F.C["aqua"], tcolor="#137a55", lw=1.2)
    box(ax, 0.57, 0.395, 0.38, 0.12, "Lysis arm (GSDMD, GSDME, NINJ1)",
        "GSDMD-N pores; NINJ1-dependent plasma-membrane\nrupture; GSDME (caspase-3) as an inflammasome-\nindependent route",
        ec=F.C["magenta"], tcolor="#b8456f", lw=1.2)
    arrow(ax, (0.44, 0.565), (0.30, 0.517))
    arrow(ax, (0.56, 0.565), (0.70, 0.517))
    box(ax, 0.62, 0.575, 0.33, 0.06, "Non-canonical caspase-4/5 → GSDMD (reported separately)",
        ec=F.C["gray"], tsize=5.6, lw=0.7)
    arrow(ax, (0.75, 0.575), (0.73, 0.517), color=F.C["gray"])
    arrow(ax, (0.43, 0.43), (0.57, 0.43), style="<|-|>", color=F.C["gray"])
    ax.text(0.50, 0.448, "separable", ha="center", fontsize=5.4, color=F.C["ink2"])
    # B: workflow
    ax.text(0.015, 0.355, "B", fontsize=10, fontweight="bold", va="top")
    ax.text(0.045, 0.352, "Analytic workflow; each layer has a distinct inferential role", fontsize=7.5, va="top")
    steps = [("Mechanism-guided\nmodules", "prespecified;\nshared genes\nexcluded"),
             ("GSE159677\nscRNA-seq", "QC, doublets,\nambient RNA,\ndepth controls"),
             ("Myeloid\nsub-clusters", "monocyte, macro-\nphage, DC,\nneutrophil"),
             ("Core vs adjacent", "3 patients;\npatient-level,\ndescriptive"),
             ("GSE120521 bulk", "4 paired plaques;\ncomposition-\nadjusted"),
             ("Two-signal\nmodel", "dimensionless;\nnot calibrated"),
             ("Caveat-aware\nprioritisation", "expert-weighted\nevidence minus\npenalties")]
    w, gap, x0 = 0.116, 0.0247, 0.02
    for i, (t, b) in enumerate(steps):
        x = x0 + i * (w + gap)
        box(ax, x, 0.13, w, 0.16, t, b, tsize=6.0, bsize=5.3, ec=F.C["blue"] if i in (1, 2) else F.C["ink2"])
        if i < len(steps) - 1:
            arrow(ax, (x + w, 0.21), (x + w + gap, 0.21))
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.075, boxstyle="round,pad=0.004,rounding_size=0.01",
                                ec="none", fc=F.C["midgray"]))
    ax.text(0.50, 0.0575, "Interpretation is constrained by sparse MEFV detection, three single-cell donors, four bulk pairs, "
            "cell-composition confounding,\nabsent genotype data and an uncalibrated model; transcript scores do not measure "
            "inflammasome activation.", ha="center", va="center", fontsize=5.8, color=F.C["ink"])
    F.save(fig, OUT, "Figure1_schematic")


if __name__ == "__main__":
    main()
