#!/usr/bin/env python
"""Revision step 5: myeloid sub-clustering (macrophage / monocyte / DC / neutrophil subsets),
MEFV detection by subset, monocyte-marker co-expression, regional composition, and the
corrected paired core (AC) vs adjacent (PA) comparison.

Input : data/processed/gse159677_singlets.h5ad (scripts/14_sc_revised_analysis.py)
Output: results/revision/tables/R12*, R13*; data/processed/gse159677_myeloid.h5ad

Usage (repo root): python -I scripts/15_sc_myeloid_subsets_and_regions.py
"""
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, "src")
import scanpy as sc
import anndata as ad
import harmonypy
import statsmodels.api as sm

SEED = 0
IN = "data/processed/gse159677_singlets.h5ad"
OUT = "data/processed/gse159677_myeloid.h5ad"
TAB = "results/revision/tables"
SIG = {
    "Mono_classical": ["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "SELL"],
    "Mono_nonclassical": ["FCGR3A", "CX3CR1", "CDKN1C", "LST1", "MS4A7"],
    "Mono_inflammatory_IL1B": ["IL1B", "CXCL8", "TNF", "CCL3", "CCL4", "NLRP3", "EREG", "PTGS2"],
    "Mac_C1Q": ["C1QA", "C1QB", "C1QC", "APOE", "CD163"],
    "Mac_TREM2_lipid": ["TREM2", "SPP1", "GPNMB", "LPL", "CD9", "FABP5", "APOC1", "ACP5"],
    "Mac_LYVE1_resident": ["LYVE1", "F13A1", "FOLR2", "SELENOP", "MRC1", "STAB1"],
    "cDC1": ["CLEC9A", "XCR1", "CADM1", "IDO1"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "CD1E"],
    "mregDC_LAMP3": ["LAMP3", "CCR7", "FSCN1", "CCL19", "CCL22"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "JCHAIN", "TCF4"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CXCR2", "CMTM2", "G0S2", "S100P", "PROK2"],
    "Proliferating": ["MKI67", "TOP2A", "STMN1"],
    "T_contaminant": ["CD3E", "CD3D", "TRAC"],
    "SMC_contaminant": ["ACTA2", "TAGLN", "MYH11"],
    "EC_contaminant": ["PECAM1", "VWF", "CDH5"],
}
NAME = {"Mono_classical": "Classical monocyte", "Mono_nonclassical": "Non-classical monocyte",
        "Mono_inflammatory_IL1B": "IL1B+ inflammatory monocyte-derived", "Mac_C1Q": "C1Q+ macrophage",
        "Mac_TREM2_lipid": "TREM2+ lipid-associated macrophage", "Mac_LYVE1_resident": "LYVE1+ resident-like macrophage",
        "cDC1": "cDC1", "cDC2": "cDC2", "mregDC_LAMP3": "LAMP3+ mregDC", "pDC": "pDC", "Neutrophil": "Neutrophil",
        "Proliferating": "Proliferating myeloid", "T_contaminant": "Mixed (T-cell transcripts)",
        "SMC_contaminant": "Mixed (SMC transcripts)", "EC_contaminant": "Mixed (endothelial transcripts)"}
LINEAGE = {"Classical monocyte": "Monocyte", "Non-classical monocyte": "Monocyte",
           "IL1B+ inflammatory monocyte-derived": "Monocyte", "C1Q+ macrophage": "Macrophage",
           "TREM2+ lipid-associated macrophage": "Macrophage", "LYVE1+ resident-like macrophage": "Macrophage",
           "Proliferating myeloid": "Macrophage", "cDC1": "Dendritic cell", "cDC2": "Dendritic cell",
           "LAMP3+ mregDC": "Dendritic cell", "pDC": "Dendritic cell", "Neutrophil": "Neutrophil",
           "Mixed (T-cell transcripts)": "Mixed", "Mixed (SMC transcripts)": "Mixed",
           "Mixed (endothelial transcripts)": "Mixed"}
COEXP = ["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "LYZ", "SELL", "CSF3R", "FCGR3B", "C1QA", "APOE", "TREM2", "CD68"]


def dummies(x, categories=None):
    cat = pd.Categorical(np.asarray(x).astype(str), categories=categories)
    d = pd.get_dummies(cat, drop_first=True).astype(float)
    d.columns = [str(c) for c in d.columns]
    d.index = x.index
    return d


def pct(v):
    return 100 * np.mean(v)


def haldane_or(a, b, c, d):
    a, b, c, d = a + .5, b + .5, c + .5, d + .5
    lo = np.log(a * d / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return np.exp(lo), np.exp(lo - 1.96 * se), np.exp(lo + 1.96 * se)


def zdelta_global(obs, by):
    p, n = obs["PYRIN_BACKBONE_score"], obs["NLRP3_BACKBONE_score"]
    return ((p - p.mean()) / p.std() - (n - n.mean()) / n.std()).groupby([obs[b] for b in by], observed=True).mean()


def main():
    s = ad.read_h5ad(IN)
    # global cell-level z (all singlets) used for the specificity delta everywhere
    s.obs["delta_cell"] = ((s.obs.PYRIN_BACKBONE_score - s.obs.PYRIN_BACKBONE_score.mean()) / s.obs.PYRIN_BACKBONE_score.std()
                           - (s.obs.NLRP3_BACKBONE_score - s.obs.NLRP3_BACKBONE_score.mean()) / s.obs.NLRP3_BACKBONE_score.std())
    m = s[s.obs.cell_type == "Macrophage/Myeloid"].copy()
    for k in ("soupx", "decontx"):
        del m.layers[k]
    sc.pp.highly_variable_genes(m, n_top_genes=2000, flavor="seurat", batch_key="sample")
    mh = m[:, m.var.highly_variable].copy()
    sc.pp.scale(mh, max_value=10)
    sc.tl.pca(mh, n_comps=30, svd_solver="arpack", random_state=SEED)
    ho = harmonypy.run_harmony(mh.obsm["X_pca"], m.obs[["sample"]], "sample", random_state=SEED, verbose=False)
    m.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr).T
    m.obsm["X_pca"] = mh.obsm["X_pca"]
    sc.pp.neighbors(m, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(m, random_state=SEED)
    sc.tl.leiden(m, resolution=0.6, random_state=SEED, flavor="igraph", n_iterations=2, directed=False, key_added="sub")
    for k, g in SIG.items():
        sc.tl.score_genes(m, [x for x in g if x in m.var_names], score_name=f"sig_{k}", random_state=SEED, use_raw=False)
    sigc = [f"sig_{k}" for k in SIG]
    cm = m.obs.groupby("sub", observed=True)[sigc].mean()
    z = (cm - cm.mean()) / cm.std()
    best = z.idxmax(1).str.replace("sig_", "")
    m.obs["subset"] = m.obs["sub"].map(best.map(NAME)).astype(str)
    m.obs["lineage"] = m.obs["subset"].map(LINEAGE)
    # cell-level lineage call (cluster-independent sensitivity)
    lin_sig = {"Monocyte": ["sig_Mono_classical", "sig_Mono_nonclassical"],
               "Macrophage": ["sig_Mac_C1Q", "sig_Mac_TREM2_lipid", "sig_Mac_LYVE1_resident"],
               "Dendritic cell": ["sig_cDC1", "sig_cDC2", "sig_mregDC_LAMP3", "sig_pDC"], "Neutrophil": ["sig_Neutrophil"]}
    L = pd.DataFrame({k: m.obs[v].max(1) for k, v in lin_sig.items()})
    m.obs["lineage_cell_level"] = np.where(L.max(1) > 0, L.idxmax(1), "Unassigned")

    o = m.obs
    counts = m.layers["counts"]
    gi = {g: i for i, g in enumerate(m.var_names)}
    for g in COEXP:
        o[f"det_{g}"] = np.asarray(counts[:, gi[g]].todense()).ravel() > 0
    # R12a subset summary
    rows = []
    for (sub, subset), g in o.groupby(["sub", "subset"], observed=True):
        per_pat = g.groupby("patient", observed=True)["MEFV_detected"].mean() * 100
        rows.append({"subcluster": sub, "subset": subset, "lineage": LINEAGE[subset], "n_cells": len(g),
                     "n_samples": g["sample"].nunique(), "max_single_sample_fraction": g["sample"].value_counts(normalize=True).max(),
                     "AC_fraction": (g.region == "AC").mean(), "MEFV_pct": pct(g.MEFV_detected),
                     "MEFV_pct_SoupX": pct(g.MEFV_detected_soupx), "MEFV_pct_DecontX": pct(g.MEFV_detected_decontx),
                     "MEFV_pct_by_patient": ";".join(f"{k}:{v:.1f}" for k, v in per_pat.items()),
                     "PYRIN_BACKBONE_mean": g.PYRIN_BACKBONE_score.mean(), "NLRP3_BACKBONE_mean": g.NLRP3_BACKBONE_score.mean(),
                     "delta_cell_level_z": g.delta_cell.mean(), "median_nUMI": g.total_counts.median(),
                     "top_signature_z": float(z.loc[sub].max()), "second_signature": z.loc[sub].drop(z.loc[sub].idxmax()).idxmax().replace("sig_", "")})
    r12a = pd.DataFrame(rows).sort_values("MEFV_pct", ascending=False)
    r12a.to_csv(f"{TAB}/R12a_myeloid_subsets.csv", index=False)
    z.to_csv(f"{TAB}/R12a2_subcluster_signature_z.csv")
    # lineage-level summary (cluster-based and cell-level)
    lin_rows = []
    for col in ("lineage", "lineage_cell_level"):
        for lin, g in o.groupby(col, observed=True):
            lin_rows.append({"definition": col, "lineage": lin, "n_cells": len(g), "AC_fraction": (g.region == "AC").mean(),
                             "MEFV_pos": int(g.MEFV_detected.sum()), "MEFV_pct": pct(g.MEFV_detected),
                             "MEFV_pct_SoupX": pct(g.MEFV_detected_soupx), "MEFV_pct_DecontX": pct(g.MEFV_detected_decontx),
                             "PYRIN_BACKBONE_mean": g.PYRIN_BACKBONE_score.mean(), "NLRP3_BACKBONE_mean": g.NLRP3_BACKBONE_score.mean(),
                             "median_nUMI": g.total_counts.median()})
    r12l = pd.DataFrame(lin_rows)
    r12l.to_csv(f"{TAB}/R12b_myeloid_lineages.csv", index=False)
    # per sample x lineage
    ps = o.groupby(["patient", "region", "lineage"], observed=True).agg(n=("MEFV_detected", "size"), MEFV_pos=("MEFV_detected", "sum"))
    ps["MEFV_pct"] = 100 * ps.MEFV_pos / ps.n
    ps.to_csv(f"{TAB}/R12c_sample_by_lineage_MEFV.csv")
    # R12d co-expression: MEFV+ vs MEFV- myeloid cells
    pos, neg = o[o.MEFV_detected], o[~o.MEFV_detected]
    co = []
    for g in COEXP:
        a_, b_ = int(pos[f"det_{g}"].sum()), int((~pos[f"det_{g}"]).sum())
        c_, d_ = int(neg[f"det_{g}"].sum()), int((~neg[f"det_{g}"]).sum())
        orr, lo, hi = haldane_or(a_, b_, c_, d_)
        co.append({"marker": g, "pct_in_MEFV_pos": 100 * a_ / (a_ + b_), "pct_in_MEFV_neg": 100 * c_ / (c_ + d_),
                   "odds_ratio": orr, "or_ci95_lo": lo, "or_ci95_hi": hi, "fisher_p": stats.fisher_exact([[a_, b_], [c_, d_]])[1]})
    pd.DataFrame(co).to_csv(f"{TAB}/R12d_MEFV_coexpression_myeloid.csv", index=False)
    # R12e GLM within myeloid: lineage rate ratios per UMI (sample FE)
    d = o[o.lineage.isin(["Monocyte", "Macrophage", "Dendritic cell"])].copy()
    X = dummies(d["lineage"], ["Macrophage", "Monocyte", "Dendritic cell"])
    X = sm.add_constant(X.join(dummies(d["sample"])))
    f = sm.GLM(d.MEFV_detected.astype(float), X, family=sm.families.Binomial(link=sm.families.links.CLogLog()),
               offset=np.log(d.total_counts)).fit()
    glm = [{"contrast": f"{t} vs Macrophage", "rate_ratio_per_UMI": np.exp(f.params[t]),
            "ci95_lo": np.exp(f.params[t] - 1.96 * f.bse[t]), "ci95_hi": np.exp(f.params[t] + 1.96 * f.bse[t]),
            "p": f.pvalues[t]} for t in ("Monocyte", "Dendritic cell")]
    pd.DataFrame(glm).to_csv(f"{TAB}/R12e_lineage_cloglog_glm.csv", index=False)
    # R12f mapping of the v1 exploratory MEFV-high Leiden cluster
    gl = s.obs.groupby("leiden", observed=True).agg(n=("leiden", "size"), MEFV_pct=("MEFV_detected", pct),
                                                     ct=("cell_type", "first"))
    top = gl.sort_values("MEFV_pct", ascending=False).index[0]
    mp = o[o.leiden == top]["subset"].value_counts()
    pd.DataFrame({"v1_style_leiden_cluster": top, "cluster_MEFV_pct": gl.loc[top, "MEFV_pct"], "subset": mp.index,
                  "n_cells": mp.values, "fraction": mp.values / mp.sum()}).to_csv(f"{TAB}/R12f_v1_MEFVhigh_cluster_mapping.csv", index=False)

    # ---------------- R13 corrected paired AC vs PA (myeloid compartment) ----------------
    metrics = {"MEFV detection (%)": ("MEFV_detected", pct), "PYRIN_BACKBONE": ("PYRIN_BACKBONE_score", np.mean),
               "NLRP3_BACKBONE": ("NLRP3_BACKBONE_score", np.mean), "Specificity delta (cell-level z)": ("delta_cell", np.mean),
               "Cytokine arm (CASP1/IL1B/IL18)": ("EFFECTOR_CYTOKINE_ARM_score", np.mean),
               "Lysis arm (GSDMD/GSDME/NINJ1)": ("EFFECTOR_LYSIS_ARM_score", np.mean),
               "Monocyte-lineage fraction (%)": ("lineage", lambda v: 100 * np.mean(np.asarray(v) == "Monocyte")),
               "TREM2+ lipid-associated macrophage fraction (%)": ("subset", lambda v: 100 * np.mean(np.asarray(v) == "TREM2+ lipid-associated macrophage"))}
    rows = []
    for lab, (col, fn) in metrics.items():
        vals = {(p, r): fn(g[col]) for (p, r), g in o.groupby(["patient", "region"], observed=True)}
        pats = sorted({p for p, _ in vals})
        diffs = np.array([vals[(p, "AC")] - vals[(p, "PA")] for p in pats])
        rows.append({"metric": lab, "PA_mean_of_patients": np.mean([vals[(p, "PA")] for p in pats]),
                     "AC_mean_of_patients": np.mean([vals[(p, "AC")] for p in pats]),
                     "mean_AC_minus_PA": diffs.mean(), "n_patients_AC_higher": int((diffs > 0).sum()),
                     "per_patient_AC_minus_PA": ";".join(f"{x:+.3f}" for x in diffs)})
    pd.DataFrame(rows).to_csv(f"{TAB}/R13a_myeloid_AC_vs_PA_corrected.csv", index=False)
    # composition-standardized MEFV detection: region-specific lineage rates applied to the pooled lineage mix
    std_rows = []
    mix = o[o.lineage.isin(["Monocyte", "Macrophage", "Dendritic cell"])].lineage.value_counts(normalize=True)
    for (p, r), g in o.groupby(["patient", "region"], observed=True):
        rates = g.groupby("lineage", observed=True)["MEFV_detected"].mean()
        comp_ = g.lineage.value_counts(normalize=True)
        std_rows.append({"patient": p, "region": r, "crude_MEFV_pct": pct(g.MEFV_detected),
                         "lineage_standardized_MEFV_pct": 100 * sum(mix[k] * rates.get(k, np.nan) for k in mix.index),
                         **{f"rate_{k}": 100 * rates.get(k, np.nan) for k in mix.index},
                         **{f"frac_{k}": comp_.get(k, 0.0) for k in ("Monocyte", "Macrophage", "Dendritic cell", "Neutrophil", "Mixed")}})
    pd.DataFrame(std_rows).to_csv(f"{TAB}/R13b_region_MEFV_lineage_standardized.csv", index=False)
    # regional GLM within myeloid: AC vs PA, with and without lineage adjustment
    g_rows = []
    dm = o[o.lineage.isin(["Monocyte", "Macrophage", "Dendritic cell"])].copy()
    for spec, adj in (("AC vs PA | patient FE + offset", False), ("AC vs PA | patient FE + lineage + offset", True)):
        X = pd.DataFrame({"AC": (dm.region == "AC").astype(float)}, index=dm.index).join(dummies(dm["patient"]))
        if adj:
            X = X.join(dummies(dm["lineage"], ["Macrophage", "Monocyte", "Dendritic cell"]))
        f = sm.GLM(dm.MEFV_detected.astype(float), sm.add_constant(X),
                   family=sm.families.Binomial(link=sm.families.links.CLogLog()), offset=np.log(dm.total_counts)).fit()
        g_rows.append({"model": spec, "rate_ratio_AC_vs_PA": np.exp(f.params["AC"]),
                       "ci95_lo": np.exp(f.params["AC"] - 1.96 * f.bse["AC"]),
                       "ci95_hi": np.exp(f.params["AC"] + 1.96 * f.bse["AC"]), "p": f.pvalues["AC"]})
    pd.DataFrame(g_rows).to_csv(f"{TAB}/R13c_region_glm_lineage_adjusted.csv", index=False)
    # all compartments paired (Table S12 analogue)
    rows = []
    so = s.obs
    for ct, gct in so.groupby("cell_type", observed=True):
        for lab, (col, fn) in list(metrics.items())[:6]:
            vals = {(p, r): fn(g[col]) for (p, r), g in gct.groupby(["patient", "region"], observed=True) if len(g)}
            pats = sorted({p for p, _ in vals if (p, "AC") in vals and (p, "PA") in vals})
            diffs = np.array([vals[(p, "AC")] - vals[(p, "PA")] for p in pats])
            rows.append({"cell_type": ct, "metric": lab, "n_patients": len(pats), "mean_AC_minus_PA": diffs.mean(),
                         "n_AC_higher": int((diffs > 0).sum()), "per_patient": ";".join(f"{x:+.3f}" for x in diffs)})
    pd.DataFrame(rows).to_csv(f"{TAB}/R13d_all_compartments_AC_vs_PA_corrected.csv", index=False)
    m.write(OUT)
    pd.set_option("display.width", 250)
    print(r12a.drop(columns=["MEFV_pct_by_patient"]).round(3).to_string(index=False))
    print(r12l.round(2).to_string(index=False))
    print(pd.DataFrame(co).round(3).to_string(index=False))
    print(pd.DataFrame(glm).round(3).to_string(index=False))
    print(pd.read_csv(f"{TAB}/R12f_v1_MEFVhigh_cluster_mapping.csv").round(3).to_string(index=False))
    print(pd.read_csv(f"{TAB}/R13a_myeloid_AC_vs_PA_corrected.csv").round(3).to_string(index=False))
    print(pd.DataFrame(std_rows).round(2).to_string(index=False))
    print(pd.DataFrame(g_rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
