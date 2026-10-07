#!/usr/bin/env python
"""Revised supplementary figures.
  S1  UMAP of all singlets (compartments) with module-score overlays
  S2  inflammasome-gene dot plot by compartment (singlets)
  S3  doublet and ambient-RNA controls
  S4  sequencing-depth controls (nUMI of MEFV+ vs MEFV- myeloid cells; fixed-depth detection)
  S5  bulk gene-level heatmap (split effector arms; GSDME recovered from DFNA5)
  S6  myeloid sub-cluster signature heatmap
  S7  expression-matched random gene-set null
  S8  superseded v1 model: structural symmetry of priming and threshold

Usage (repo root): python -I scripts/23_fig_supplement.py [S1 S2 ...]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.insert(0, "src")
import anndata as ad
from pyrinplaque import figstyle as F

TAB = "results/revision/tables"
OUT = "results/revision/figures/supplementary"
CT = ["Macrophage/Myeloid", "Endothelial", "T/NK", "SMC/Fibroblast", "Mast", "B/Plasma"]
CT_LAB = {"Macrophage/Myeloid": "Myeloid", "Endothelial": "Endothelial", "T/NK": "T/NK",
          "SMC/Fibroblast": "SMC/fibroblast", "Mast": "Mast", "B/Plasma": "B/plasma"}
DIV = LinearSegmentedColormap.from_list("div", F.DIVERGING)
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
CT_COL = {"Macrophage/Myeloid": F.C["blue"], "Endothelial": F.C["orange"], "T/NK": F.C["aqua"],
          "SMC/Fibroblast": "#b9b8b3", "Mast": F.C["ink2"], "B/Plasma": "#7d7c77"}


def s1(s):
    o = s.obs
    U = s.obsm["X_umap"]
    fig, axs = plt.subplots(1, 4, figsize=(F.WIDTH["double"], 52 * F.MM))
    ax = axs[0]
    for c in CT:
        k = (o.cell_type == c).values
        ax.scatter(U[k, 0], U[k, 1], s=0.15, color=CT_COL[c], lw=0, rasterized=True)
        cx, cy = np.median(U[k], axis=0)
        ax.text(cx, cy, CT_LAB[c], fontsize=5.5, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
    ax.set_title("Compartments (singlets)", fontweight="normal")
    for ax, col, lab in zip(axs[1:], ["PYRIN_BACKBONE_score", "NLRP3_BACKBONE_score", "EFFECTOR_CYTOKINE_ARM_score"],
                            ["PYRIN_BACKBONE", "NLRP3_BACKBONE", "Cytokine arm"]):
        v = o[col].values
        lo, hi = np.percentile(v, [1, 99])
        order = np.argsort(v)
        sc_ = ax.scatter(U[order, 0], U[order, 1], c=v[order], s=0.15, cmap=SEQ, vmin=lo, vmax=hi, lw=0, rasterized=True)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.01)
        cb.ax.tick_params(labelsize=5)
        ax.set_title(lab, fontweight="normal")
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[:].set_visible(False)
    F.save(fig, OUT, "FigureS1_umap_modules")


def s2(s):
    genes = [("Pyrin backbone", ["MEFV", "RHOA", "RAC1", "CDC42", "PKN1", "PKN2", "YWHAB", "YWHAZ", "PSTPIP1"]),
             ("NLRP3 backbone", ["NLRP3", "NEK7", "TXNIP", "NFKB1", "RELA"]),
             ("Shared / effector", ["PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "GSDME", "NINJ1", "CASP4", "CASP5"])]
    allg = [g for _, gs in genes for g in gs]
    o = s.obs
    cnt = pd.DataFrame(np.asarray(s[:, allg].layers["counts"].todense()) > 0, columns=allg, index=o.index)
    ex = pd.DataFrame(np.asarray(s[:, allg].X.todense()), columns=allg, index=o.index)
    frac = cnt.groupby(o.cell_type, observed=True).mean().reindex(CT)
    mean = ex.groupby(o.cell_type, observed=True).mean().reindex(CT)
    z = (mean - mean.mean()) / mean.std()
    fig, ax = plt.subplots(figsize=(F.WIDTH["double"], 55 * F.MM))
    for i, c in enumerate(CT):
        ax.scatter(np.arange(len(allg)), [i] * len(allg), s=frac.loc[c].values * 60 + 0.3, c=z.loc[c].values,
                   cmap=SEQ, vmin=-1, vmax=2, lw=0)
    ax.set_xticks(range(len(allg)), allg, rotation=90, fontsize=6)
    ax.set_yticks(range(len(CT)), [CT_LAB[c] for c in CT])
    ax.set_ylim(len(CT) - 0.4, -0.6)
    x = 0
    for lab, gs in genes:
        ax.text(x + (len(gs) - 1) / 2, -1.0, lab, ha="center", fontsize=6.2, color=F.C["ink"])
        if x:
            ax.axvline(x - 0.5, color=F.C["grid"], lw=0.8)
        x += len(gs)
    for fr in (0.05, 0.25, 0.5, 1.0):
        ax.scatter([], [], s=fr * 60 + 0.3, color=F.C["ink2"], label=f"{int(fr * 100)}%")
    ax.legend(title="% cells", fontsize=5.3, title_fontsize=5.5, loc="upper left", bbox_to_anchor=(1.0, 1.0), labelspacing=1)
    F.save(fig, OUT, "FigureS2_gene_dotplot")


def s3():
    d = pd.read_csv(f"{TAB}/R03a_doublets_by_sample.csv", index_col=0)
    a = pd.read_csv(f"{TAB}/R09a_ambient_MEFV.csv", index_col=0).reindex(CT)
    fig, axs = plt.subplots(1, 2, figsize=(F.WIDTH["double"], 60 * F.MM), gridspec_kw={"wspace": 0.45})
    ax = axs[0]
    x = np.arange(len(d))
    ax.bar(x - 0.18, d.scDblFinder_rate_pct, width=0.34, color=F.C["blue"], label="scDblFinder (primary)")
    ax.bar(x + 0.18, d.scrublet_rate_pct, width=0.34, color=F.C["orange"], label="Scrublet")
    ax.set_xticks(x, [f"{i}\n(n={n:,})" for i, n in zip(d.index, d.n_QC)], fontsize=5.6)
    ax.set_ylabel("Predicted doublets (%)")
    ax.legend(fontsize=5.6)
    ax.set_title("Doublet rates by sample", fontweight="normal")
    F.panel_label(ax, "A", dx=-0.15)
    ax = axs[1]
    y = np.arange(len(CT))[::-1]
    for col, colr, lab, off in (("observed_MEFV_pct", F.C["blue"], "Observed", 0.27), ("MEFV_pct_SoupX", F.C["orange"], "SoupX", 0.09),
                                ("MEFV_pct_DecontX", F.C["aqua"], "DecontX", -0.09),
                                ("expected_ambient_MEFV_pct", F.C["ink2"], "Expected from ambient RNA", -0.27)):
        ax.barh(y + off, a[col], height=0.17, color=colr, label=lab)
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_yticks(y, [CT_LAB[c] for c in CT])
    ax.set_xlabel("Cells with detectable MEFV (%) (symlog)")
    ax.legend(fontsize=5.6, loc="lower right")
    ax.set_title("Ambient-RNA controls", fontweight="normal")
    F.panel_label(ax, "B", dx=-0.25)
    F.save(fig, OUT, "FigureS3_doublet_ambient")


def s4(s):
    o = s.obs
    r08 = pd.read_csv(f"{TAB}/R08_fixed_depth_2000UMI.csv", index_col=0).reindex(CT)
    my = o[o.cell_type == "Macrophage/Myeloid"]
    fig, axs = plt.subplots(1, 2, figsize=(F.WIDTH["double"], 60 * F.MM), gridspec_kw={"wspace": 0.4})
    ax = axs[0]
    smp = sorted(my["sample"].unique())
    for i, sm in enumerate(smp):
        g = my[my["sample"] == sm]
        for k, (lab, colr, off) in enumerate((("MEFV−", F.C["blue"], -0.18), ("MEFV+", F.C["orange"], 0.18))):
            v = np.log10(g.loc[g.MEFV_detected == (lab == "MEFV+"), "total_counts"])
            if len(v) < 3:
                continue
            vp = ax.violinplot([v], positions=[i + off], widths=0.32, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor(colr)
                b.set_alpha(0.75)
            ax.hlines(np.median(v), i + off - 0.1, i + off + 0.1, color=F.C["ink"], lw=0.9)
    ax.set_xticks(range(len(smp)), smp, fontsize=6)
    ax.set_ylabel("log$_{10}$ UMIs per cell")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=F.C["blue"], lw=5, alpha=.75, label="MEFV−"),
                       Line2D([], [], color=F.C["orange"], lw=5, alpha=.75, label="MEFV+")], fontsize=5.6)
    ax.set_title("Myeloid cells: library size by MEFV status", fontweight="normal")
    F.panel_label(ax, "A", dx=-0.15)
    ax = axs[1]
    y = np.arange(len(CT))[::-1]
    ax.barh(y, r08.MEFV_pct_at_2000UMI, color=F.C["blue"], height=0.6)
    for yi, c in zip(y, CT):
        ax.text(r08.loc[c, "MEFV_pct_at_2000UMI"] + 0.05, yi, f"{r08.loc[c, 'MEFV_pct_at_2000UMI']:.2f}%  (n={int(r08.loc[c, 'n_cells_ge2000UMI']):,})",
                va="center", fontsize=5.5, color=F.C["ink2"])
    ax.set_yticks(y, [CT_LAB[c] for c in CT])
    ax.set_xlabel("MEFV detection after downsampling every cell to 2,000 UMIs (%)")
    ax.set_title("Fixed-depth detection", fontweight="normal")
    F.panel_label(ax, "B", dx=-0.25)
    F.save(fig, OUT, "FigureS4_depth_controls")


def s5():
    df = pd.read_excel("data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx")
    df["name"] = df["name"].replace({"DFNA5": "GSDME"})
    cols = [c for c in df.columns if "FPKM" in c]
    e = np.log2(df.dropna(subset=["name"]).groupby("name")[cols].max() + 1)
    groups = [("Pyrin backbone", ["MEFV", "RHOA", "RAC1", "CDC42", "PKN1", "PKN2", "YWHAB", "YWHAZ", "PSTPIP1"]),
              ("NLRP3 backbone", ["NLRP3", "NEK7", "TXNIP", "NFKB1", "RELA"]),
              ("Cytokine arm", ["CASP1", "IL1B", "IL18"]), ("Lysis arm", ["GSDMD", "GSDME", "NINJ1"]),
              ("Caspase-4/5", ["CASP4", "CASP5"]), ("Adaptor", ["PYCARD"]),
              ("Myeloid marker", ["LYZ", "CD68", "CD14", "FCGR3A", "LST1", "C1QA", "C1QB", "C1QC", "MS4A7", "TYROBP"])]
    genes = [g for _, gs in groups for g in gs]
    order = ["MB1_FPKM_stable", "MB3_FPKM_stable", "MB5_FPKM_stable", "MB9_FPKM_stable",
             "MB2_FPKM_unstable", "MB4_FPKM_unstable", "MB6_FPKM_unstable", "MB10_FPKM_unstable"]
    Z = e.loc[genes, order]
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1), axis=0)
    fig, ax = plt.subplots(figsize=(F.WIDTH["onehalf"], 150 * F.MM))
    im = ax.imshow(Z.values, cmap=DIV, norm=TwoSlopeNorm(0, -2, 2), aspect="auto")
    ax.set_yticks(range(len(genes)), genes, fontsize=5.8)
    ax.set_xticks(range(8), [f"P{i} {s}" for s in ("stable", "unstable") for i in (1, 2, 3, 4)], rotation=45, ha="right", fontsize=6)
    ax.axvline(3.5, color="white", lw=2)
    y = 0
    for lab, gs in groups:
        if y:
            ax.axhline(y - 0.5, color="white", lw=2)
        ax.text(8.0, y + (len(gs) - 1) / 2, lab, fontsize=6, va="center")
        y += len(gs)
    ax.tick_params(length=0)
    ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.25)
    cb.set_label("Row z score, log$_2$(FPKM+1)", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    F.save(fig, OUT, "FigureS5_bulk_gene_heatmap")


def s6():
    z = pd.read_csv(f"{TAB}/R12a2_subcluster_signature_z.csv", index_col=0)
    sub = pd.read_csv(f"{TAB}/R12a_myeloid_subsets.csv").set_index("subcluster")
    z.columns = [c.replace("sig_", "") for c in z.columns]
    z.index = [f"{i}: {sub.loc[int(i), 'subset']} (n={int(sub.loc[int(i), 'n_cells'])})" for i in z.index]
    fig, ax = plt.subplots(figsize=(F.WIDTH["double"], 95 * F.MM))
    im = ax.imshow(z.values, cmap=DIV, norm=TwoSlopeNorm(0, -3, 3), aspect="auto")
    ax.set_xticks(range(z.shape[1]), z.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(z.shape[0]), z.index, fontsize=5.8)
    ax.tick_params(length=0)
    ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Signature score, z across sub-clusters", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    F.save(fig, OUT, "FigureS6_myeloid_signatures")


def s7():
    r = pd.read_csv(f"{TAB}/R10a_random_geneset_null.csv", index_col=0).reindex(CT)
    fig, ax = plt.subplots(figsize=(F.WIDTH["onehalf"], 60 * F.MM))
    y = np.arange(len(CT))[::-1]
    ax.errorbar(r.null_mean, y, xerr=1.96 * r.null_sd, fmt="o", ms=3, color=F.C["ink2"], ecolor=F.C["muted"], label="Null mean ± 1.96 SD (200 sets)")
    ax.scatter(r.observed, y, s=20, color=F.C["blue"], zorder=3, label="Observed PYRIN_BACKBONE")
    for yi, c in zip(y, CT):
        ax.text(r.loc[c, "observed"] + 0.03, yi + 0.18, f"{r.loc[c, 'empirical_percentile']:.0f}th pct", fontsize=5.4, color=F.C["ink2"])
    ax.set_yticks(y, [CT_LAB[c] for c in CT])
    ax.set_xlabel("Mean module score")
    ax.legend(fontsize=5.6, loc="lower right")
    F.save(fig, OUT, "FigureS7_random_null")


def s8():
    t = pd.read_csv(f"{TAB}/R34_ode1_symmetry_check.csv")
    fig, ax = plt.subplots(figsize=(F.WIDTH["onehalf"], 55 * F.MM))
    y = np.arange(len(t))[::-1]
    ax.barh(y, t.time_to_50pct_IL1B, color=[F.C["ink2"], F.C["blue"], F.C["orange"], F.C["orange"], F.C["blue"]], height=0.6)
    for yi, (_, r) in zip(y, t.iterrows()):
        ax.text(r.time_to_50pct_IL1B + 0.1, yi, f"t½ {r.time_to_50pct_IL1B:.1f}; peak IL-1β {r.peak_IL1B:.3f}",
                va="center", fontsize=5.5, color=F.C["ink2"])
    ax.set_yticks(y, t.case, fontsize=6)
    ax.set_xlabel("Time to 50% of peak IL-1β (v1 model, dimensionless)")
    ax.set_xlim(0, 11)
    ax.set_title("v1 model: equal shifts of priming and threshold give identical outputs", fontweight="normal")
    F.save(fig, OUT, "FigureS8_v1_symmetry")


def main():
    F.apply()
    want = set(sys.argv[1:]) or {"S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"}
    if want & {"S1", "S2", "S4"}:
        s = ad.read_h5ad("data/processed/gse159677_singlets.h5ad")
        if "X_umap" not in s.obsm:
            raise SystemExit("singlets h5ad lacks X_umap")
        "S1" in want and s1(s)
        "S2" in want and s2(s)
        "S4" in want and s4(s)
    "S3" in want and s3()
    "S5" in want and s5()
    "S6" in want and s6()
    "S7" in want and s7()
    "S8" in want and s8()


if __name__ == "__main__":
    main()
