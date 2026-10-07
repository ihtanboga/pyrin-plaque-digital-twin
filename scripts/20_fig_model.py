#!/usr/bin/env python
"""Revised Figure 5: two-signal pyrin model (v2) — scenario trajectories, matched 20% perturbations
and LHS partial rank correlations. Reads R30-R33 and data/processed/ode2_scenario_trajectories.csv.

Usage (repo root): python -I scripts/20_fig_model.py
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.insert(0, "src")
from pyrinplaque import figstyle as F

TAB = "results/revision/tables"
OUT = "results/revision/figures"
SCEN = [("baseline", "Baseline", F.C["ink2"]), ("pyrin_sensitized_low_threshold", "Lower pyrin threshold", F.C["blue"]),
        ("high_inflammatory_priming", "Higher priming (signal 1)", F.C["orange"]),
        ("inhibited_or_delayed_activation", "Inhibited/delayed gate", F.C["aqua"])]
STATES = [("P_active", "Active pyrin"), ("IL1B_out", "Extracellular IL-1β"),
          ("IL18_out", "Extracellular IL-18"), ("lysed", "Lysed fraction")]
OUTS = [("IL1B_AUC", "IL-1β\nAUC"), ("IL18_AUC", "IL-18\nAUC"), ("IL1B_to_IL18_AUC_ratio", "IL-1β:\nIL-18"),
        ("time_to_50pct_peak_IL1B", "t½\nIL-1β"), ("time_to_50pct_final_lysis", "t½\nlysis")]
PERT = [("threshold_minus20pct", "Threshold −20%"), ("signal2_plus20pct", "Signal 2 +20%"),
        ("priming_plus20pct", "Priming +20%"), ("inhibition_0p2", "Inhibition 0.2")]
PAR = {"priming_input": "Priming (signal 1)", "signal2_input": "Signal 2 (RhoA-inactivating)",
       "pyrin_activation_threshold": "Pyrin activation threshold", "inhibition_factor": "Inhibition factor",
       "pyrin_activation_rate": "Pyrin activation rate", "proIL1B_priming_synthesis": "pro-IL-1β synthesis",
       "proIL18_constitutive_synthesis": "pro-IL-18 synthesis", "CASP1_activation_rate": "Caspase-1 activation",
       "GSDMD_cleavage_rate": "GSDMD cleavage", "pore_release_rate": "Pore release", "lysis_rate": "Lysis rate",
       "lysis_pore_threshold": "Lysis pore threshold", "decay": "Active-species decay", "clearance": "Cytokine clearance"}
PRCC_OUT = [("IL1B_AUC", "IL-1β\nAUC"), ("IL18_AUC", "IL-18\nAUC"), ("IL1B_to_IL18_AUC_ratio", "IL-1β:\nIL-18"),
            ("time_to_50pct_peak_IL1B", "t½\nIL-1β"), ("lysed_fraction_end", "Final\nlysis"),
            ("time_to_50pct_final_lysis", "t½\nlysis")]
DIV = LinearSegmentedColormap.from_list("div", F.DIVERGING)


def main():
    F.apply()
    tr = pd.read_csv("data/processed/ode2_scenario_trajectories.csv")
    mp = pd.read_csv(f"{TAB}/R31_ode2_matched_perturbations.csv").set_index("perturbation")
    pr = pd.read_csv(f"{TAB}/R33_ode2_global_prcc.csv")
    fig = plt.figure(figsize=(F.WIDTH["double"], 130 * F.MM))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.35], hspace=0.62, wspace=0.42, left=0.07, right=0.98)
    for j, (col, lab) in enumerate(STATES):
        ax = fig.add_subplot(gs[0, j])
        for key, name, colr in SCEN:
            d = tr[tr.scenario == key]
            ax.plot(d.t, d[col], color=colr, lw=1.1, label=name)
        ax.set_xlim(0, 30)
        ax.set_title(lab, fontsize=7, fontweight="normal")
        ax.set_xlabel("Dimensionless time")
        if j == 0:
            F.panel_label(ax, "A", dx=-0.3)
            ax.set_ylabel("Dimensionless level")
    fig.axes[0].legend(loc="upper right", fontsize=5.8, handlelength=1.4)

    ax = fig.add_axes([0.13, 0.08, 0.27, 0.36])
    M = np.array([[mp.loc[p, f"{o}_pct_change"] for o, _ in OUTS] for p, _ in PERT])
    lim = max(5, np.nanmax(np.abs(M)))
    im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:+.1f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(M[i, j]) > 0.6 * lim else F.C["ink"])
    ax.set_xticks(range(len(OUTS)), [o[1] for o in OUTS])
    ax.set_yticks(range(len(PERT)), [p[1] for p in PERT])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Equal-size (20%) perturbations:\n% change from baseline", fontsize=7, fontweight="normal")
    F.panel_label(ax, "B", dx=-0.32)
    cb = fig.colorbar(im, cax=fig.add_axes([0.41, 0.08, 0.008, 0.36]))
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("% change", fontsize=6)

    ax = fig.add_axes([0.63, 0.08, 0.30, 0.36])
    P = pr.pivot_table(index="parameter", columns="output", values="PRCC").reindex(list(PAR))[[o for o, _ in PRCC_OUT]]
    im = ax.imshow(P.values, cmap=DIV, norm=TwoSlopeNorm(0, -1, 1), aspect="auto")
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            v = P.values[i, j]
            if abs(v) >= 0.3:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.3,
                        color="white" if abs(v) > 0.65 else F.C["ink"])
    ax.set_xticks(range(len(PRCC_OUT)), [o[1] for o in PRCC_OUT], fontsize=5.8)
    ax.set_yticks(range(len(PAR)), [PAR[k] for k in P.index], fontsize=5.8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Global sensitivity: PRCC\n(LHS, n = 1,000; |PRCC| ≥ 0.3 printed)", fontsize=7, fontweight="normal")
    F.panel_label(ax, "C", dx=-0.42)
    cb = fig.colorbar(im, cax=fig.add_axes([0.94, 0.08, 0.008, 0.36]))
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("PRCC", fontsize=6)
    F.save(fig, OUT, "Figure5_two_signal_model")


if __name__ == "__main__":
    main()
