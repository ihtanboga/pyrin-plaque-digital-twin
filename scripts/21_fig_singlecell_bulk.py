#!/usr/bin/env python
"""Revised Figures 2-4 (single-cell compartments + depth; myeloid subsets; corrected regional
and bulk triangulation). Reads results/revision/tables/R0x-R2x and the processed h5ad files.

Usage (repo root): python -I scripts/21_fig_singlecell_bulk.py
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, "src")
import anndata as ad
from pyrinplaque import figstyle as F

TAB = "results/revision/tables"
OUT = "results/revision/figures"
CT = ["Macrophage/Myeloid", "Endothelial", "T/NK", "SMC/Fibroblast", "Mast", "B/Plasma"]
CT_LAB = {"Macrophage/Myeloid": "Myeloid", "Endothelial": "Endothelial", "T/NK": "T/NK",
          "SMC/Fibroblast": "SMC/fibroblast", "Mast": "Mast", "B/Plasma": "B/plasma"}
PAT_MK = {"Patient 1": "o", "Patient 2": "s", "Patient 3": "^"}
DIV = LinearSegmentedColormap.from_list("div", F.DIVERGING)
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])


def tidy_heat(ax):
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def fig2(s):
    o = s.obs
    comp = pd.read_csv(f"{TAB}/R05_compartment_summary_revised.csv", index_col=0)
    pat = pd.read_csv(f"{TAB}/R06b_patient_compartment_values.csv")
    strat = pd.read_csv(f"{TAB}/R07c_depth_stratified_detection.csv")
    glm = pd.read_csv(f"{TAB}/R07b_cloglog_glm_depth_offset.csv")
    fig = plt.figure(figsize=(F.WIDTH["double"], 165 * F.MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.05, 1], hspace=0.75, wspace=0.55, left=0.1, right=0.97, top=0.96, bottom=0.07)
    # A: MEFV detection
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(CT))[::-1]
    ax.barh(y, comp.loc[CT, "MEFV_pct"], color=F.C["blue"], height=0.62, label="Observed (singlets)")
    for p, mk in PAT_MK.items():
        v = pat[pat.patient == p].set_index("cell_type").reindex(CT)["MEFV_pct"]
        ax.scatter(v, y, marker=mk, s=9, facecolor="white", edgecolor=F.C["ink"], lw=0.6, zorder=3, label=p)
    ax.scatter(comp.loc[CT, "MEFV_pct_SoupX"], y + 0.0, marker="|", s=40, color=F.C["orange"], lw=1.2, zorder=4, label="SoupX-corrected")
    ax.scatter(comp.loc[CT, "MEFV_pct_DecontX"], y + 0.0, marker="|", s=40, color=F.C["aqua"], lw=1.2, zorder=4, label="DecontX-corrected")
    ax.set_yticks(y, [CT_LAB[c] for c in CT])
    ax.set_xlabel("Cells with detectable MEFV (%)")
    for yi, c in zip(y, CT):
        ax.text(comp.loc[c, "MEFV_pct"] + 0.25, yi, f"{comp.loc[c, 'MEFV_pct']:.2f}%", va="center", fontsize=5.8, color=F.C["ink2"])
    ax.legend(fontsize=5.5, loc="lower right", handletextpad=0.3, borderaxespad=0.2)
    ax.set_title("MEFV detection by compartment", fontweight="normal")
    F.panel_label(ax, "A", dx=-0.32)
    # B: standardized compartment means (cell-level z, the Table 2 definition)
    ax = fig.add_subplot(gs[0, 1])
    cols = {"PYRIN_BACKBONE_score": "PYRIN\nbackbone", "NLRP3_BACKBONE_score": "NLRP3\nbackbone",
            "EFFECTOR_CYTOKINE_ARM_score": "Cytokine\narm", "EFFECTOR_LYSIS_ARM_score": "Lysis\narm"}
    Z = pd.DataFrame({lab: ((o[c] - o[c].mean()) / o[c].std()).groupby(o.cell_type, observed=True).mean()
                      for c, lab in cols.items()}).reindex(CT)
    Z.insert(2, "Specificity\ndelta", comp.loc[CT, "delta_cell_level_z"].values)
    im = ax.imshow(Z.values, cmap=DIV, norm=TwoSlopeNorm(0, -1.5, 1.5), aspect="auto")
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            v = Z.values[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.8, color="white" if abs(v) > 1.0 else F.C["ink"])
    ax.set_xticks(range(Z.shape[1]), Z.columns, fontsize=5.8)
    ax.set_yticks(range(len(CT)), [CT_LAB[c] for c in CT])
    tidy_heat(ax)
    ax.axvline(1.5, color="white", lw=2)
    ax.axvline(2.5, color="white", lw=2)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("Mean of cell-level z", fontsize=6)
    ax.set_title("Compartment means (cell-level standardisation)", fontweight="normal")
    F.panel_label(ax, "B", dx=-0.3)
    # C: violins
    ax = fig.add_subplot(gs[1, :])
    pos = np.arange(len(CT))
    for k, (col, colr, off) in enumerate((("PYRIN_BACKBONE_score", F.C["blue"], -0.18), ("NLRP3_BACKBONE_score", F.C["orange"], 0.18))):
        data = [o.loc[o.cell_type == c, col].values for c in CT]
        vp = ax.violinplot(data, positions=pos + off, widths=0.34, showextrema=False, showmedians=False)
        for b in vp["bodies"]:
            b.set_facecolor(colr)
            b.set_alpha(0.75)
            b.set_edgecolor("none")
        med = [np.median(d) for d in data]
        ax.hlines(med, pos + off - 0.1, pos + off + 0.1, color=F.C["ink"], lw=1.0)
    ax.axhline(0, color=F.C["grid"], lw=0.6, zorder=0)
    ax.set_xticks(pos, [CT_LAB[c] for c in CT])
    ax.set_ylabel("Module score")
    ax.legend(handles=[Line2D([], [], color=F.C["blue"], lw=5, alpha=0.75, label="PYRIN_BACKBONE"),
                       Line2D([], [], color=F.C["orange"], lw=5, alpha=0.75, label="NLRP3_BACKBONE"),
                       Line2D([], [], color=F.C["ink"], lw=1, label="Median")], ncol=3, loc="upper right", fontsize=6)
    ax.set_title("Cell-level module-score distributions (singlets)", fontweight="normal")
    F.panel_label(ax, "C", dx=-0.04)
    # D: depth-stratified detection
    ax = fig.add_subplot(gs[2, 0])
    show = {"Macrophage/Myeloid": F.C["blue"], "Endothelial": F.C["orange"], "SMC/Fibroblast": F.C["aqua"], "T/NK": F.C["ink2"]}
    for c, colr in show.items():
        d = strat[strat.cell_type == c].sort_values("nUMI_quintile")
        ax.plot(d.nUMI_quintile, d.MEFV_pct, marker="o", ms=3, color=colr, lw=1.1, label=CT_LAB[c])
    for c in show:
        d = strat[strat.cell_type == c].sort_values("nUMI_quintile")
        for _, r in d.iterrows():
            if r.n < 300 and r.MEFV_pct > 0:
                ax.annotate(f"n={int(r.n)}", (r.nUMI_quintile, r.MEFV_pct), xytext=(4, -2), textcoords="offset points",
                            fontsize=5.2, color=F.C["ink2"])
    ax.set_ylabel("MEFV detection (%)")
    ax.set_xlabel("Global nUMI quintile (Q1 lowest depth)")
    ax.legend(fontsize=5.6, loc="upper left")
    ax.set_title("Detection at matched sequencing depth", fontweight="normal")
    F.panel_label(ax, "D", dx=-0.3)
    # E: cloglog rate ratios
    ax = fig.add_subplot(gs[2, 1])
    g = glm[glm.model.str.contains("sample FE")].reset_index(drop=True)
    labels = [r.contrast.replace("Macrophage/Myeloid", "Myeloid").replace("SMC/Fibroblast", "SMC/fibroblast")
              .replace("myeloid vs non-myeloid", "Myeloid vs all non-myeloid") for r in g.itertuples()]
    yy = np.arange(len(g))[::-1]
    ax.errorbar(g.rate_ratio_per_UMI, yy, xerr=[g.rate_ratio_per_UMI - g.ci95_lo, g.ci95_hi - g.rate_ratio_per_UMI],
                fmt="o", ms=3.5, color=F.C["ink"], ecolor=F.C["ink2"], elinewidth=0.9, capsize=1.5)
    ax.axvline(1, color=F.C["muted"], lw=0.7)
    ax.set_xscale("log")
    ax.set_yticks(yy, labels, fontsize=6)
    ax.set_xlabel("MEFV detection rate ratio per UMI (log scale)")
    for xi, hi, yi in zip(g.rate_ratio_per_UMI, g.ci95_hi, yy):
        ax.text(hi * 1.25, yi, f"{xi:.1f}", va="center", fontsize=5.6, color=F.C["ink2"])
    ax.set_xlim(0.08, 400)
    ax.set_title("cloglog GLM, offset log(nUMI), sample fixed effects", fontweight="normal")
    F.panel_label(ax, "E", dx=-0.55)
    F.save(fig, OUT, "Figure2_compartments_depth")


def fig3(m):
    o = m.obs
    sub = pd.read_csv(f"{TAB}/R12a_myeloid_subsets.csv")
    co = pd.read_csv(f"{TAB}/R12d_MEFV_coexpression_myeloid.csv")
    fig = plt.figure(figsize=(F.WIDTH["double"], 150 * F.MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1, 1.15], hspace=0.55, wspace=0.62,
                          left=0.06, right=0.90, top=0.95, bottom=0.08)
    # A: UMAP by lineage with subset labels
    ax = fig.add_subplot(gs[0, 0])
    U = m.obsm["X_umap"]
    for lin, colr in F.LINEAGE.items():
        k = (o.lineage == lin).values
        ax.scatter(U[k, 0], U[k, 1], s=0.6, color=colr, lw=0, rasterized=True, label=f"{lin} ({k.sum():,})")
    for subset, g in o.groupby("subset", observed=True):
        if len(g) < 40:
            continue
        cx, cy = np.median(U[(o.subset == subset).values], axis=0)
        lab = {"TREM2+ lipid-associated macrophage": "TREM2+ LAM", "LYVE1+ resident-like macrophage": "LYVE1+ mac.",
               "Mixed (T-cell transcripts)": "Mixed (T)", "Mixed (SMC transcripts)": "Mixed (SMC)",
               "Non-classical monocyte": "Non-classical\nmonocyte", "Classical monocyte": "Classical\nmonocyte",
               "C1Q+ macrophage": "C1Q+ mac."}.get(subset, subset)
        ax.text(cx, cy, lab, fontsize=5, ha="center", va="center", color=F.C["ink"],
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=8, fontsize=5.5, loc="upper left", bbox_to_anchor=(0.0, -0.06), ncol=3, handletextpad=0.1,
              borderaxespad=0, columnspacing=0.8)
    ax.set_title("Myeloid sub-clusters (Harmony)", fontweight="normal")
    F.panel_label(ax, "A", dx=-0.06)
    # B: MEFV by subset
    ax = fig.add_subplot(gs[0, 1])
    agg = sub.groupby("subset").apply(lambda g: pd.Series({"n": g.n_cells.sum(),
                                                            "MEFV_pct": np.average(g.MEFV_pct, weights=g.n_cells),
                                                            "lineage": g.lineage.iloc[0]}), include_groups=False)
    agg = agg[agg.n >= 40].sort_values("MEFV_pct")
    yy = np.arange(len(agg))
    ax.barh(yy, agg.MEFV_pct, color=[F.LINEAGE[l] for l in agg.lineage], height=0.65)
    per_pat = o.groupby(["subset", "patient"], observed=True)["MEFV_detected"].agg(["mean", "size"]).reset_index()
    for p, mk in PAT_MK.items():
        d = per_pat[(per_pat.patient == p) & (per_pat["size"] >= 20)].set_index("subset")
        v = [100 * d.loc[s_, "mean"] if s_ in d.index else np.nan for s_ in agg.index]
        ax.scatter(v, yy, marker=mk, s=9, facecolor="white", edgecolor=F.C["ink"], lw=0.6, zorder=3, label=f"{p} (≥20 cells)")
    for yi, (s_, r) in enumerate(agg.iterrows()):
        ax.text(1.02, yi, f"{r.MEFV_pct:.1f}%  n={int(r.n):,}", va="center", ha="left", fontsize=5.4,
                color=F.C["ink2"], transform=ax.get_yaxis_transform())
    ax.set_xlim(0, 24)
    ax.set_yticks(yy, agg.index, fontsize=5.8)
    ax.set_xlabel("Cells with detectable MEFV (%)")
    ax.legend(fontsize=5.3, loc="lower right")
    ax.set_title("MEFV detection by myeloid subset", fontweight="normal")
    F.panel_label(ax, "B", dx=-0.62)
    # C: marker dot plot
    ax = fig.add_subplot(gs[1, 0])
    markers = ["MEFV", "CD14", "FCN1", "S100A8", "VCAN", "FCGR3A", "IL1B", "C1QA", "APOE", "TREM2", "SPP1", "LYVE1", "CD1C", "CLEC9A", "LAMP3", "CSF3R"]
    order = list(agg.index[::-1])
    X = m[:, markers].layers["counts"]
    Xd = pd.DataFrame(np.asarray(X.todense()), columns=markers, index=o.index)
    lx = pd.DataFrame(np.asarray(m[:, markers].X.todense()), columns=markers, index=o.index)
    frac = (Xd > 0).groupby(o.subset, observed=True).mean().reindex(order)
    mean = lx.groupby(o.subset, observed=True).mean().reindex(order)
    zs = (mean - mean.mean()) / mean.std()
    for i, s_ in enumerate(order):
        ax.scatter(np.arange(len(markers)), [i] * len(markers), s=frac.loc[s_].values * 40 + 0.5,
                   c=zs.loc[s_].values, cmap=SEQ, vmin=-1, vmax=2, lw=0)
    ax.set_xticks(range(len(markers)), markers, rotation=90, fontsize=5.6)
    ax.set_yticks(range(len(order)), order, fontsize=5.6)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title("Marker detection (size) and scaled mean (colour)", fontweight="normal")
    sm_ = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(-1, 2))
    cb = fig.colorbar(sm_, ax=ax, fraction=0.03, pad=0.17, location="right")
    cb.ax.tick_params(labelsize=5)
    cb.set_label("Scaled mean (z)", fontsize=5.5)
    for fr in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=fr * 40 + 0.5, color=F.C["ink2"], label=f"{int(fr * 100)}%")
    ax.legend(title="% cells", fontsize=5, title_fontsize=5.5, loc="upper left", bbox_to_anchor=(1.0, 1.0), labelspacing=0.8)
    F.panel_label(ax, "C", dx=-0.42)
    # D: co-expression in MEFV+ vs MEFV- myeloid cells
    ax = fig.add_subplot(gs[1, 1])
    co = co.set_index("marker").loc[["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "SELL", "LYZ", "CSF3R", "C1QA", "APOE", "TREM2", "CD68"]]
    yy = np.arange(len(co))[::-1]
    ax.hlines(yy, co.pct_in_MEFV_neg, co.pct_in_MEFV_pos, color=F.C["grid"], lw=2.2)
    ax.scatter(co.pct_in_MEFV_neg, yy, s=14, color=F.C["blue"], zorder=3, label="MEFV− myeloid cells")
    ax.scatter(co.pct_in_MEFV_pos, yy, s=14, color=F.C["orange"], zorder=3, label="MEFV+ myeloid cells")
    for yi, (g, r) in zip(yy, co.iterrows()):
        ax.text(101, yi, f"OR {r.odds_ratio:.1f}", va="center", fontsize=5.3, color=F.C["ink2"])
    ax.set_yticks(yy, co.index, fontsize=5.8)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Cells with marker detected (%)")
    ax.legend(fontsize=5.5, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, borderaxespad=0.2)
    ax.set_title(" ", fontweight="normal")
    F.panel_label(ax, "D", dx=-0.2)
    F.save(fig, OUT, "Figure3_myeloid_subsets")


def paired(ax, df, ylab, labels=("PA", "AC"), fmt="{:+.2f}"):
    for p, r in df.iterrows():
        ax.plot([0, 1], [r[labels[0]], r[labels[1]]], color=F.C["ink2"], lw=0.9, marker=PAT_MK.get(p, "o"), ms=3.2,
                mfc="white", mec=F.C["ink"], mew=0.6)
        yv = r[labels[1]]
        rng = np.nanmax(df[labels[1]]) - np.nanmin(df[labels[1]]) or 1
        close = [q for q in df.index if q != p and abs(df.loc[q, labels[1]] - yv) < 0.04 * rng]
        nud = 0.05 * rng * (1 if any(q > p for q in close) else -1 if close else 0)
        ax.text(1.06, yv - nud, p.replace("Patient ", "P"), fontsize=5.3, va="center", color=F.C["ink2"])
    d = (df[labels[1]] - df[labels[0]])
    ax.set_xticks([0, 1], labels)
    ax.set_xlim(-0.3, 1.35)
    ax.set_ylabel(ylab, fontsize=6)
    ax.set_title(f"mean Δ {fmt.format(d.mean())}; {labels[1]} higher {int((d > 0).sum())}/{len(d)}",
                 fontsize=5.6, fontweight="normal", color=F.C["ink2"], loc="left")


def fig4(m):
    o = m.obs
    reg = {}
    for lab, col, fn in (("MEFV detection (%)", "MEFV_detected", lambda v: 100 * v.mean()),
                         ("Monocyte-lineage fraction (%)", "lineage", lambda v: 100 * (v == "Monocyte").mean()),
                         ("PYRIN_BACKBONE", "PYRIN_BACKBONE_score", np.mean),
                         ("Specificity delta", "delta_cell", np.mean)):
        reg[lab] = o.groupby(["patient", "region"], observed=True)[col].agg(fn).unstack()
    bulk = pd.read_csv(f"{TAB}/R21_bulk_module_scores.csv", index_col=0)
    res = pd.read_csv(f"{TAB}/R22_bulk_paired_results.csv").set_index("module")
    rec = pd.read_csv(f"{TAB}/R23_bulk_effect_vs_myeloid_enrichment.csv").set_index("module")
    fig = plt.figure(figsize=(F.WIDTH["double"], 170 * F.MM))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.2], hspace=1.0, wspace=0.8, left=0.07, right=0.97, top=0.93, bottom=0.07)
    for j, (lab, df) in enumerate(reg.items()):
        ax = fig.add_subplot(gs[0, j])
        paired(ax, df, lab, fmt="{:+.1f}" if "%" in lab else "{:+.3f}")
        if j == 0:
            F.panel_label(ax, "A", dx=-0.45, dy=1.22)
            ax.text(0, 1.24, "Myeloid compartment, adjacent (PA) vs core (AC), n = 3 patients", transform=ax.transAxes,
                    fontsize=7, va="bottom")
    mods = [("PYRIN_BACKBONE", "PYRIN_BACKBONE"), ("NLRP3_BACKBONE", "NLRP3_BACKBONE"),
            ("EFFECTOR_CYTOKINE_ARM", "Cytokine arm"), ("EFFECTOR_LYSIS_ARM", "Lysis arm")]
    for j, (mod, lab) in enumerate(mods):
        ax = fig.add_subplot(gs[1, j])
        df = bulk.pivot_table(index="patient", columns="status", values=mod)
        paired(ax, df, f"{lab} score", labels=("stable", "unstable"))
        ax.set_title(ax.get_title(loc="left") + f"; d$_z$ {res.loc[mod, 'cohens_dz']:.2f}", fontsize=5.6,
                     fontweight="normal", color=F.C["ink2"], loc="left")
        if j == 0:
            F.panel_label(ax, "B", dx=-0.45, dy=1.22)
            ax.text(0, 1.24, "Bulk GSE120521, stable vs unstable region, n = 4 plaques", transform=ax.transAxes,
                    fontsize=7, va="bottom")
    # C: raw vs residualized
    ax = fig.add_subplot(gs[2, 0:2])
    order = ["MYELOID_MARKER", "EFFECTOR_CYTOKINE_ARM", "NONCANONICAL_CASP4_5", "PYRIN_FULL", "PYRIN_BACKBONE",
             "EFFECTOR_LYSIS_ARM", "NLRP3_FULL", "NLRP3_BACKBONE"]
    nice = {"MYELOID_MARKER": "Myeloid marker", "EFFECTOR_CYTOKINE_ARM": "Cytokine arm", "NONCANONICAL_CASP4_5": "Caspase-4/5",
            "PYRIN_FULL": "PYRIN_FULL", "PYRIN_BACKBONE": "PYRIN_BACKBONE", "EFFECTOR_LYSIS_ARM": "Lysis arm",
            "NLRP3_FULL": "NLRP3_FULL", "NLRP3_BACKBONE": "NLRP3_BACKBONE"}
    yy = np.arange(len(order))[::-1]
    raw = res.loc[order, "mean_paired_diff"]
    rr = res.loc[order, "resid_mean_paired_diff"].fillna(raw)
    ax.hlines(yy, rr, raw, color=F.C["grid"], lw=2.2)
    ax.scatter(raw, yy, s=16, color=F.C["blue"], zorder=3, label="Unadjusted")
    ax.scatter(rr[order[1:]], yy[1:], s=16, color=F.C["orange"], zorder=3, label="Myeloid-residualised")
    ax.axvline(0, color=F.C["muted"], lw=0.7)
    ax.set_yticks(yy, [nice[k] for k in order], fontsize=6)
    ax.set_xlabel("Mean paired difference, unstable − stable (z units)")
    ax.legend(fontsize=5.6, loc="lower right")
    ax.set_title("Composition adjustment applied to every module", fontweight="normal")
    F.panel_label(ax, "C", dx=-0.3)
    # D: effect vs sc myeloid enrichment
    ax = fig.add_subplot(gs[2, 2:4])
    r2 = rec.drop(index=["LEGACY_GENERIC_8GENE"])
    ax.scatter(r2.myeloid_log2_enrichment_mean, r2.mean_paired_diff, s=18, color=F.C["blue"], zorder=3)
    off = {"PYRIN_BACKBONE": (4, -9), "NLRP3_FULL": (5, -9), "EFFECTOR_LYSIS_ARM": (-8, 6), "PYRIN_FULL": (4, 3),
           "NLRP3_BACKBONE": (4, 2), "NONCANONICAL_CASP4_5": (4, 2), "EFFECTOR_CYTOKINE_ARM": (-18, 6), "MYELOID_MARKER": (-30, 6)}
    for k, r in r2.iterrows():
        ax.annotate(nice[k], (r.myeloid_log2_enrichment_mean, r.mean_paired_diff), xytext=off.get(k, (3, 2)),
                    textcoords="offset points", fontsize=5.3, color=F.C["ink2"])
    from scipy import stats
    rho = stats.spearmanr(r2.myeloid_log2_enrichment_mean, r2.mean_paired_diff)
    ax.text(0.02, 0.97, f"Spearman ρ = {rho.statistic:.2f} (8 modules)", transform=ax.transAxes, fontsize=5.8, va="top", color=F.C["ink2"])
    ax.set_xlabel("Single-cell myeloid enrichment of module genes\n(mean log$_2$ ratio, myeloid vs other cells)")
    ax.set_ylabel("Bulk unstable − stable (z units)")
    ax.set_title("Bulk effect sizes track module myeloid specificity", fontweight="normal")
    F.panel_label(ax, "D", dx=-0.22)
    F.save(fig, OUT, "Figure4_regional_bulk")


def main():
    F.apply()
    s = ad.read_h5ad("data/processed/gse159677_singlets.h5ad", backed=None)
    fig2(s)
    del s
    m = ad.read_h5ad("data/processed/gse159677_myeloid.h5ad")
    fig3(m)
    fig4(m)


if __name__ == "__main__":
    main()
