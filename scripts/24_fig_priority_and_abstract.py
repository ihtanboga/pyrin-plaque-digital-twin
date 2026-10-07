#!/usr/bin/env python
"""Revised Figure 6 (caveat-aware priority score) and the graphical abstract, both drawn from
the revision tables.

Usage (repo root): python -I scripts/24_fig_priority_and_abstract.py
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, "src")
from pyrinplaque import figstyle as F

TAB = "results/revision/tables"
OUT = "results/revision/figures"
SHORT = {"L2-01": "MEFV/pyrin threshold axis", "L1-01": "Myeloid-compartment pyrin-backbone state",
         "L2-06": "NLRP3 (comparator)", "L1-02": "Classical-monocyte MEFV-high subset",
         "L2-02": "RHOA-RAC1-CDC42 context", "L2-03": "PKN1/PKN2 / 14-3-3 regulation",
         "L1-03": "T/NK NLRP3-backbone compartment (comparator)", "L2-04": "PSTPIP1-pyrin interaction",
         "L1-05": "Unstable-plaque bulk pyrin signal", "L2-05": "Shared effector arms (cytokine, lysis)",
         "L1-04": "Plaque-core MEFV signal (reversed)"}


def fig6():
    t = pd.read_csv(f"{TAB}/R41_priority_table_revised.csv")
    t = t.sort_values("final_priority_score")
    y = np.arange(len(t))
    fig, ax = plt.subplots(figsize=(F.WIDTH["double"], 95 * F.MM))
    ax.barh(y, t.raw_positive_score, color=F.C["blue"], height=0.6, label="Weighted positive evidence")
    ax.barh(y, -t.total_penalty, color=F.C["red"], height=0.6, label="Weighted penalties")
    ax.scatter(t.final_priority_score, y, marker="D", s=22, color=F.C["ink"], zorder=3, label="Final score")
    for yi, v in zip(y, t.final_priority_score):
        ax.text(max(v, 0) + 1.5, yi + 0.33, f"{v:.1f}", fontsize=5.6, color=F.C["ink"])
    for x, lab in ((35, "moderate ≥ 35"), (55, "high ≥ 55")):
        ax.axvline(x, color=F.C["ink2"], lw=0.7, ls=(0, (3, 2)))
        ax.text(x + 0.8, len(t) - 0.35, lab, fontsize=5.6, color=F.C["ink2"])
    ax.axvline(0, color=F.C["ink2"], lw=0.7)
    ax.set_yticks(y, [SHORT[c] for c in t.candidate_id], fontsize=6)
    ax.set_xlabel("Points (expert-defined research prioritisation; not a clinical or target ranking)")
    ax.set_xlim(-42, 82)
    ax.legend(fontsize=5.8, loc="lower right")
    F.save(fig, OUT, "Figure6_priority")


def box(ax, x, y, w, h, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.015", fc=fc, ec=ec, lw=1.0))


def graphical_abstract():
    lin = pd.read_csv(f"{TAB}/R12b_myeloid_lineages.csv")
    lin = lin[lin.definition == "lineage"].set_index("lineage")
    bulk = pd.read_csv(f"{TAB}/R22_bulk_paired_results.csv").set_index("module")
    comp = pd.read_csv(f"{TAB}/R05_compartment_summary_revised.csv", index_col=0)
    fig = plt.figure(figsize=(13.28 / 2.54 * 2, 5.31 / 2.54 * 2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.95, "Pyrin (MEFV) versus NLRP3 programs in human carotid plaque: a hypothesis-generating map",
            ha="center", fontsize=10.5, fontweight="bold", color=F.C["ink"])
    cols = [(0.015, "Data and controls", F.C["ink2"]), (0.265, "Where is MEFV?", F.C["orange"]),
            (0.515, "What drives bulk signals?", F.C["blue"]), (0.765, "Two-signal model", F.C["aqua"])]
    for x, title, colr in cols:
        box(ax, x, 0.10, 0.22, 0.74, "white", colr)
        ax.text(x + 0.11, 0.79, title, ha="center", fontsize=9, fontweight="bold", color=colr)
    ax.text(0.125, 0.47, "GSE159677 scRNA-seq\n3 patients, core + adjacent\n"
            f"{int(comp.n_cells.sum()):,} singlets\n\nGSE120521 bulk RNA-seq\n4 stable/unstable pairs\n\n"
            "Doublets, ambient RNA,\nsequencing depth controlled;\nsample labels verified", ha="center", va="center",
            fontsize=7.2, color=F.C["ink"], linespacing=1.25)
    # MEFV by lineage mini bar chart
    a1 = fig.add_axes([0.30, 0.30, 0.17, 0.40])
    order = ["Monocyte", "Dendritic cell", "Macrophage"]
    vals = [lin.loc[k, "MEFV_pct"] for k in order]
    a1.barh(range(3)[::-1], vals, color=[F.LINEAGE[k] for k in order], height=0.6)
    for i, v in zip(range(3)[::-1], vals):
        a1.text(v + 0.2, i, f"{v:.1f}%", va="center", fontsize=7)
    a1.set_yticks(range(3)[::-1], order, fontsize=7)
    a1.set_xlim(0, 10)
    a1.set_xlabel("MEFV+ cells (%)", fontsize=7)
    a1.tick_params(labelsize=6.5)
    ax.text(0.375, 0.17, "MEFV is sparse and monocyte-enriched;\nrare in TREM2+ lipid-associated macrophages",
            ha="center", fontsize=6.6, color=F.C["ink"])
    # bulk dumbbell
    a2 = fig.add_axes([0.585, 0.33, 0.12, 0.36])
    mods = [("EFFECTOR_CYTOKINE_ARM", "Cytokine arm"), ("PYRIN_BACKBONE", "Pyrin backbone"), ("NLRP3_BACKBONE", "NLRP3 backbone")]
    for i, (k, lab) in enumerate(mods[::-1]):
        a2.plot([bulk.loc[k, "resid_mean_paired_diff"], bulk.loc[k, "mean_paired_diff"]], [i, i], color=F.C["grid"], lw=3)
        a2.scatter(bulk.loc[k, "mean_paired_diff"], i, color=F.C["blue"], s=18, zorder=3)
        a2.scatter(bulk.loc[k, "resid_mean_paired_diff"], i, color=F.C["orange"], s=18, zorder=3)
    a2.set_yticks(range(3), [m[1] for m in mods[::-1]], fontsize=6.5)
    a2.axvline(0, color=F.C["muted"], lw=0.6)
    a2.set_xlabel("Unstable − stable", fontsize=6.5)
    a2.tick_params(labelsize=6)
    ax.text(0.625, 0.17, "Blue: unadjusted; orange: adjusted for myeloid\ncontent. Bulk signals track myeloid abundance",
            ha="center", fontsize=6.6, color=F.C["ink"])
    ax.text(0.875, 0.47, "Signal 1 (priming)\n→ IL-1β magnitude and\nIL-1β:IL-18 ratio\n\n"
            "Signal 2 / pyrin gate\n→ timing of IL-1β, IL-18\nand lysis\n\nDimensionless, uncalibrated;\n"
            "generates testable predictions", ha="center", va="center", fontsize=7.2, color=F.C["ink"], linespacing=1.25)
    ax.text(0.5, 0.035, "Transcript scores do not measure inflammasome activation; findings prioritise cell-type-resolved, "
            "protein-level and genotype-aware validation.", ha="center", fontsize=6.8, style="italic", color=F.C["ink2"])
    for x0 in (0.237, 0.487, 0.737):
        ax.add_patch(FancyArrowPatch((x0, 0.47), (x0 + 0.026, 0.47), arrowstyle="-|>", mutation_scale=10, color=F.C["ink2"]))
    F.save(fig, OUT, "Graphical_abstract")


def main():
    F.apply()
    fig6()
    graphical_abstract()


if __name__ == "__main__":
    main()
