#!/usr/bin/env python
"""Revision step 4: revised single-cell analysis set and robustness analyses (GSE159677).

Primary analysis set = QC cells (v1 QC, v1 clustering/annotation reproduced) minus
scDblFinder doublets, with sample labels corrected (see 10_sc_reproduce_v1.py).

Outputs (results/revision/tables):
  R03 doublets per sample / compartment (scDblFinder, Scrublet)
  R04 per-sample myeloid MEFV detection (all QC, singlets, SoupX, DecontX)
  R05 compartment summary (revised Table 2) incl. specificity delta under three definitions
  R06 patient- and sample-level specificity delta
  R07 depth: nUMI MEFV+ vs MEFV-, cloglog GLM with log(nUMI) offset, depth-stratified detection
  R08 fixed-depth (2,000 UMI) downsampling: detection and module scores by compartment
  R09 ambient: observed vs expected-ambient MEFV detection; SoupX/DecontX-corrected detection and scores
  R10 expression-matched random gene-set null (200 sets) and leave-one-sample-out
  R11 module scores by compartment incl. split effector arms
Writes data/processed/gse159677_singlets.h5ad.

Usage (repo root): python -I scripts/14_sc_revised_analysis.py <ambient_work_dir>
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, "src")
import scanpy as sc
import anndata as ad
import statsmodels.api as sm
from pyrinplaque import modules as M

SEED = 0
AMB = sys.argv[1]
H5AD = "data/processed/gse159677_v1repro.h5ad"
OUT = "data/processed/gse159677_singlets.h5ad"
TAB = "results/revision/tables"
META = pd.read_csv("data/metadata/gse159677_sample_metadata.csv", dtype={"aggr_suffix": str})
CT_ORDER = ["Macrophage/Myeloid", "Endothelial", "T/NK", "SMC/Fibroblast", "Mast", "B/Plasma"]
SPLIT = {"EFFECTOR_CYTOKINE_ARM": ["CASP1", "IL1B", "IL18"],
         "EFFECTOR_LYSIS_ARM": ["GSDMD", "GSDME", "NINJ1"],
         "NONCANONICAL_CASP4_5": ["CASP4", "CASP5"]}


def score_sets():
    s = M.scoring_sets("config/gene_modules.yaml")
    s = {"PYRIN_FULL": s["PYRIN_FULL"], "PYRIN_BACKBONE": s["PYRIN_BACKBONE"], "NLRP3_FULL": s["NLRP3_FULL"],
         "NLRP3_BACKBONE": s["NLRP3_BACKBONE"], "LEGACY_GENERIC_8GENE": s["GENERIC_PYROPTOSIS"]}
    s.update(SPLIT)
    return s


def lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1
    Xn = sp.diags(1e4 / tot) @ X
    Xn = sp.csr_matrix(Xn, dtype=np.float32)
    Xn.data = np.log1p(Xn.data)
    return Xn


def score(adata, sets, prefix=""):
    for name, genes in sets.items():
        g = [x for x in genes if x in adata.var_names]
        sc.tl.score_genes(adata, g, score_name=f"{prefix}{name}_score", ctrl_size=50, n_bins=25,
                          random_state=SEED, use_raw=False)


def dummies(x, categories=None):
    cat = pd.Categorical(np.asarray(x).astype(str), categories=categories)
    d = pd.get_dummies(cat, drop_first=True).astype(float)
    d.columns = [str(c) for c in d.columns]
    d.index = x.index
    return d


def zdelta(obs, by, p="PYRIN_BACKBONE_score", n="NLRP3_BACKBONE_score"):
    zp = (obs[p] - obs[p].mean()) / obs[p].std()
    zn = (obs[n] - obs[n].mean()) / obs[n].std()
    return (zp - zn).groupby([obs[b] for b in by], observed=True).mean()


def main():
    os.makedirs(TAB, exist_ok=True)
    CACHE = "data/processed/gse159677_qc_controls.h5ad"
    if os.path.exists(CACHE):
        a = ad.read_h5ad(CACHE)
    else:
        a = stage1()
        a.write(CACHE)
    stage2(a)


def stage1():
    a = ad.read_h5ad(H5AD)
    a.raw = None
    # ---------------- doublets ----------------
    dbl, cont = [], []
    for _, m in META.iterrows():
        d = os.path.join(AMB, m.gsm)
        x = pd.read_csv(os.path.join(d, "doublets.csv")).set_index("barcode")
        dbl.append(x)
        cont.append(pd.read_csv(os.path.join(d, "decontx_contamination.csv")).set_index("barcode"))
    dbl, cont = pd.concat(dbl), pd.concat(cont)
    a.obs = a.obs.join(dbl).join(cont)
    # Scrublet per sample (expected rate 0.8% per 1,000 recovered cells)
    a.obs["scrublet_score"], a.obs["scrublet_doublet"] = np.nan, False
    for s_ in a.obs["sample"].unique():
        sub = ad.AnnData(a[a.obs["sample"] == s_].layers["counts"].copy())
        sub.obs_names = a.obs_names[a.obs["sample"] == s_]
        rate = min(0.25, 0.008 * sub.n_obs / 1000)
        sc.pp.scrublet(sub, expected_doublet_rate=rate, random_state=SEED, verbose=False)
        a.obs.loc[sub.obs_names, "scrublet_score"] = sub.obs["doublet_score"].values
        a.obs.loc[sub.obs_names, "scrublet_doublet"] = sub.obs["predicted_doublet"].astype(bool).values
    a.obs["is_doublet"] = a.obs["scDblFinder_class"].eq("doublet")
    r03a = a.obs.groupby("sample").agg(n_QC=("is_doublet", "size"), scDblFinder_doublets=("is_doublet", "sum"),
                                       scrublet_doublets=("scrublet_doublet", "sum"))
    r03a["scDblFinder_rate_pct"] = 100 * r03a.scDblFinder_doublets / r03a.n_QC
    r03a["scrublet_rate_pct"] = 100 * r03a.scrublet_doublets / r03a.n_QC
    r03b = a.obs.groupby("cell_type").agg(n_QC=("is_doublet", "size"), scDblFinder_doublets=("is_doublet", "sum"),
                                          scrublet_doublets=("scrublet_doublet", "sum"),
                                          MEFV_pct_all=("MEFV_detected", lambda v: 100 * v.mean()))
    r03b["MEFV_pct_doublets"] = a.obs[a.obs.is_doublet].groupby("cell_type")["MEFV_detected"].mean() * 100
    r03b["MEFV_pct_singlets"] = a.obs[~a.obs.is_doublet].groupby("cell_type")["MEFV_detected"].mean() * 100
    agree = pd.crosstab(a.obs.is_doublet, a.obs.scrublet_doublet)
    r03a.to_csv(f"{TAB}/R03a_doublets_by_sample.csv")
    r03b.to_csv(f"{TAB}/R03b_doublets_by_compartment.csv")
    agree.to_csv(f"{TAB}/R03c_doublet_method_agreement.csv")

    # ---------------- ambient-corrected MEFV ----------------
    gene_idx = {g: i for i, g in enumerate(a.var_names)}
    corr = {}
    for meth in ("soupx", "decontx"):
        mats = []
        for _, m in META.iterrows():
            d = os.path.join(AMB, m.gsm)
            cells = pd.read_csv(os.path.join(d, "cells.csv"))["barcode"].astype(str).values
            X = sp.csr_matrix(sio.mmread(os.path.join(d, f"{meth}_counts.mtx")).T)  # cells x genes
            mats.append(pd.Series(range(len(cells)), index=cells).to_frame("i").assign(blk=len(mats)))
            corr.setdefault(meth, []).append((cells, X))
    for meth in ("soupx", "decontx"):
        cells = np.concatenate([c for c, _ in corr[meth]])
        X = sp.vstack([x for _, x in corr[meth]]).tocsr()
        order = pd.Index(cells).get_indexer(a.obs_names)
        assert (order >= 0).all()
        a.layers[meth] = X[order]
        a.obs[f"MEFV_detected_{meth}"] = np.asarray(a.layers[meth][:, gene_idx["MEFV"]].todense()).ravel() >= 1
    return a


def stage2(a):
    # ---------------- singlet analysis set ----------------
    s = a[~a.obs.is_doublet].copy()
    s.X = lognorm(s.layers["counts"])
    sets = score_sets()
    score(s, sets)
    s.obs["log_nUMI"] = np.log(s.obs["total_counts"])
    o = s.obs

    # R04 per-sample myeloid detection
    rows = []
    for (gsm, pat, reg), g in a.obs.groupby(["gsm", "patient", "region"], observed=True):
        my = g[g.cell_type == "Macrophage/Myeloid"]
        mys = my[~my.is_doublet]
        rows.append({"gsm": gsm, "patient": pat, "region": reg, "n_QC_cells": len(g),
                     "n_myeloid_QC": len(my), "MEFV_pos_myeloid_QC": int(my.MEFV_detected.sum()),
                     "MEFV_pct_myeloid_QC": 100 * my.MEFV_detected.mean(),
                     "n_myeloid_singlets": len(mys), "MEFV_pos_myeloid_singlets": int(mys.MEFV_detected.sum()),
                     "MEFV_pct_myeloid_singlets": 100 * mys.MEFV_detected.mean(),
                     "MEFV_pct_myeloid_singlets_SoupX": 100 * mys.MEFV_detected_soupx.mean(),
                     "MEFV_pct_myeloid_singlets_DecontX": 100 * mys.MEFV_detected_decontx.mean(),
                     "median_nUMI_myeloid_singlets": float(mys.total_counts.median())})
    r04 = pd.DataFrame(rows).sort_values(["patient", "region"])
    pooled = {}
    for reg in ("AC", "PA"):
        x = r04[r04.region == reg]
        pooled[reg] = {"pooled_pct": 100 * x.MEFV_pos_myeloid_singlets.sum() / x.n_myeloid_singlets.sum(),
                       "mean_of_patient_pct": x.MEFV_pct_myeloid_singlets.mean()}
    pooled["ALL"] = {"pooled_pct": 100 * r04.MEFV_pos_myeloid_singlets.sum() / r04.n_myeloid_singlets.sum(),
                     "mean_of_patient_pct": r04.MEFV_pct_myeloid_singlets.mean()}
    r04.to_csv(f"{TAB}/R04_per_sample_myeloid_MEFV.csv", index=False)
    pd.DataFrame(pooled).T.to_csv(f"{TAB}/R04b_pooled_vs_mean_of_samples.csv")

    # R05 compartment summary + delta definitions
    comp = o.groupby("cell_type", observed=True).agg(
        n_cells=("cell_type", "size"), MEFV_pct=("MEFV_detected", lambda v: 100 * v.mean()),
        MEFV_pct_SoupX=("MEFV_detected_soupx", lambda v: 100 * v.mean()),
        MEFV_pct_DecontX=("MEFV_detected_decontx", lambda v: 100 * v.mean()),
        PYRIN_BACKBONE_mean=("PYRIN_BACKBONE_score", "mean"), NLRP3_BACKBONE_mean=("NLRP3_BACKBONE_score", "mean"),
        median_nUMI=("total_counts", "median")).reindex(CT_ORDER)
    comp["delta_cell_level_z"] = zdelta(o, ["cell_type"]).reindex(CT_ORDER)
    zc = lambda v: (v - v.mean()) / v.std(ddof=1)
    comp["delta_compartment_mean_z"] = zc(comp.PYRIN_BACKBONE_mean) - zc(comp.NLRP3_BACKBONE_mean)
    comp["check_weighted_sum_cell_level"] = float((comp.n_cells * comp.delta_cell_level_z).sum())
    pat = zdelta(o, ["patient", "cell_type"]).unstack()
    comp["delta_patient_mean"] = pat.mean().reindex(CT_ORDER)
    comp["delta_patient_min"] = pat.min().reindex(CT_ORDER)
    comp["delta_patient_max"] = pat.max().reindex(CT_ORDER)
    comp["rank_delta_by_patient"] = ";".join(f"{p}:{pat.loc[p].idxmax()}" for p in pat.index)
    comp.to_csv(f"{TAB}/R05_compartment_summary_revised.csv")
    # R06 patient/sample-level values
    samp = o.groupby(["patient", "region", "cell_type"], observed=True).agg(
        n=("cell_type", "size"), PYRIN_BACKBONE=("PYRIN_BACKBONE_score", "mean"),
        NLRP3_BACKBONE=("NLRP3_BACKBONE_score", "mean"),
        MEFV_pct=("MEFV_detected", lambda v: 100 * v.mean()))
    samp["delta_cell_level_z"] = zdelta(o, ["patient", "region", "cell_type"])
    samp.to_csv(f"{TAB}/R06a_sample_compartment_values.csv")
    patv = o.groupby(["patient", "cell_type"], observed=True).agg(
        n=("cell_type", "size"), PYRIN_BACKBONE=("PYRIN_BACKBONE_score", "mean"),
        NLRP3_BACKBONE=("NLRP3_BACKBONE_score", "mean"), MEFV_pct=("MEFV_detected", lambda v: 100 * v.mean()))
    patv["delta_cell_level_z"] = zdelta(o, ["patient", "cell_type"])
    patv["PYRIN_BACKBONE_rank_within_patient"] = patv.groupby(level=0)["PYRIN_BACKBONE"].rank(ascending=False)
    patv["delta_rank_within_patient"] = patv.groupby(level=0)["delta_cell_level_z"].rank(ascending=False)
    patv.to_csv(f"{TAB}/R06b_patient_compartment_values.csv")

    # R07 depth
    rows = []
    for lab, g in (("all compartments", o), ("myeloid", o[o.cell_type == "Macrophage/Myeloid"])):
        pos, neg = g[g.MEFV_detected], g[~g.MEFV_detected]
        rows.append({"subset": lab, "n_MEFV_pos": len(pos), "n_MEFV_neg": len(neg),
                     "median_nUMI_pos": pos.total_counts.median(), "median_nUMI_neg": neg.total_counts.median(),
                     "median_nGenes_pos": pos.n_genes_by_counts.median(), "median_nGenes_neg": neg.n_genes_by_counts.median(),
                     "mannwhitney_p_nUMI": stats.mannwhitneyu(pos.total_counts, neg.total_counts).pvalue})
    for smp, g in o[o.cell_type == "Macrophage/Myeloid"].groupby("sample", observed=True):
        pos, neg = g[g.MEFV_detected], g[~g.MEFV_detected]
        rows.append({"subset": f"myeloid {smp}", "n_MEFV_pos": len(pos), "n_MEFV_neg": len(neg),
                     "median_nUMI_pos": pos.total_counts.median(), "median_nUMI_neg": neg.total_counts.median(),
                     "median_nGenes_pos": pos.n_genes_by_counts.median(), "median_nGenes_neg": neg.n_genes_by_counts.median(),
                     "mannwhitney_p_nUMI": stats.mannwhitneyu(pos.total_counts, neg.total_counts).pvalue if len(pos) else np.nan})
    pd.DataFrame(rows).to_csv(f"{TAB}/R07a_nUMI_MEFVpos_vs_neg.csv", index=False)
    # cloglog binomial GLM with log(nUMI) offset = Poisson detection model; exp(beta) = rate ratio per UMI
    glm_rows = []
    d = o[o.cell_type.isin(["Macrophage/Myeloid", "Endothelial", "T/NK", "SMC/Fibroblast"])].copy()
    d["ct"] = pd.Categorical(d.cell_type, categories=["T/NK", "Macrophage/Myeloid", "Endothelial", "SMC/Fibroblast"])
    for spec, extra in (("compartment + sample FE + offset log(nUMI)", True), ("compartment + offset log(nUMI)", False)):
        X = dummies(d["cell_type"], ["T/NK", "Macrophage/Myeloid", "Endothelial", "SMC/Fibroblast"])
        if extra:
            X = X.join(dummies(d["sample"]))
        X = sm.add_constant(X)
        fam = sm.families.Binomial(link=sm.families.links.CLogLog())
        f = sm.GLM(d.MEFV_detected.astype(float), X, family=fam, offset=d.log_nUMI).fit()
        for term in ["Macrophage/Myeloid", "Endothelial", "SMC/Fibroblast"]:
            b, se = f.params[term], f.bse[term]
            glm_rows.append({"model": spec, "contrast": f"{term} vs T/NK", "rate_ratio_per_UMI": np.exp(b),
                             "ci95_lo": np.exp(b - 1.96 * se), "ci95_hi": np.exp(b + 1.96 * se), "p": f.pvalues[term]})
    # myeloid vs all non-myeloid (all compartments, incl. zero-event ones)
    d2 = o.copy()
    X = sm.add_constant(pd.DataFrame({"myeloid": (d2.cell_type == "Macrophage/Myeloid").astype(float)}, index=d2.index)
                        .join(dummies(d2["sample"])))
    f = sm.GLM(d2.MEFV_detected.astype(float), X, family=sm.families.Binomial(link=sm.families.links.CLogLog()),
               offset=d2.log_nUMI).fit()
    b, se = f.params["myeloid"], f.bse["myeloid"]
    glm_rows.append({"model": "myeloid vs all other + sample FE + offset log(nUMI)", "contrast": "myeloid vs non-myeloid",
                     "rate_ratio_per_UMI": np.exp(b), "ci95_lo": np.exp(b - 1.96 * se), "ci95_hi": np.exp(b + 1.96 * se),
                     "p": f.pvalues["myeloid"]})
    # region within myeloid (corrected labels), patient FE
    dm = o[o.cell_type == "Macrophage/Myeloid"].copy()
    X = sm.add_constant(pd.DataFrame({"AC": (dm.region == "AC").astype(float)}, index=dm.index)
                        .join(dummies(dm["patient"])))
    f = sm.GLM(dm.MEFV_detected.astype(float), X, family=sm.families.Binomial(link=sm.families.links.CLogLog()),
               offset=dm.log_nUMI).fit()
    b, se = f.params["AC"], f.bse["AC"]
    glm_rows.append({"model": "myeloid only: region + patient FE + offset log(nUMI)", "contrast": "AC vs PA",
                     "rate_ratio_per_UMI": np.exp(b), "ci95_lo": np.exp(b - 1.96 * se), "ci95_hi": np.exp(b + 1.96 * se),
                     "p": f.pvalues["AC"]})
    pd.DataFrame(glm_rows).to_csv(f"{TAB}/R07b_cloglog_glm_depth_offset.csv", index=False)
    # depth-stratified detection (global nUMI quintiles)
    o["nUMI_quintile"] = pd.qcut(o.total_counts, 5, labels=[f"Q{i}" for i in range(1, 6)])
    strat = o.groupby(["nUMI_quintile", "cell_type"], observed=True).agg(
        n=("MEFV_detected", "size"), MEFV_pct=("MEFV_detected", lambda v: 100 * v.mean()),
        nUMI_median=("total_counts", "median")).reset_index()
    strat.to_csv(f"{TAB}/R07c_depth_stratified_detection.csv", index=False)

    # R08 fixed-depth downsampling to 2,000 UMI
    keep = o.total_counts >= 2000
    dsa = ad.AnnData(s[keep].layers["counts"].copy(), obs=o.loc[keep, ["cell_type", "sample", "patient", "region"]].copy(),
                     var=pd.DataFrame(index=s.var_names))
    sc.pp.downsample_counts(dsa, counts_per_cell=2000, random_state=SEED)
    dsa.obs["MEFV_detected_ds"] = np.asarray(dsa[:, "MEFV"].X.todense()).ravel() > 0
    dsa.X = lognorm(dsa.X)
    score(dsa, {k: sets[k] for k in ("PYRIN_BACKBONE", "NLRP3_BACKBONE")}, prefix="ds_")
    r08 = dsa.obs.groupby("cell_type", observed=True).agg(
        n_cells_ge2000UMI=("cell_type", "size"), MEFV_pct_at_2000UMI=("MEFV_detected_ds", lambda v: 100 * v.mean()),
        PYRIN_BACKBONE_mean_ds=("ds_PYRIN_BACKBONE_score", "mean"),
        NLRP3_BACKBONE_mean_ds=("ds_NLRP3_BACKBONE_score", "mean")).reindex(CT_ORDER)
    r08["delta_cell_level_z_ds"] = zdelta(dsa.obs, ["cell_type"], "ds_PYRIN_BACKBONE_score",
                                          "ds_NLRP3_BACKBONE_score").reindex(CT_ORDER)
    r08.to_csv(f"{TAB}/R08_fixed_depth_2000UMI.csv")
    # depth-adjusted module scores (OLS with log nUMI + sample)
    ols_rows = []
    for mod in ("PYRIN_BACKBONE_score", "NLRP3_BACKBONE_score"):
        dd = o.copy()
        X = dummies(dd["cell_type"], CT_ORDER)
        X = sm.add_constant(X.join(dummies(dd["sample"])).assign(log_nUMI=dd.log_nUMI))
        f = sm.OLS(dd[mod], X).fit()
        for t in CT_ORDER[1:]:
            ols_rows.append({"module": mod, "contrast": f"{t} minus Macrophage/Myeloid (adjusted for log nUMI + sample)",
                             "estimate": f.params[t], "ci95_lo": f.conf_int().loc[t, 0], "ci95_hi": f.conf_int().loc[t, 1]})
        ols_rows.append({"module": mod, "contrast": "log_nUMI slope", "estimate": f.params["log_nUMI"],
                         "ci95_lo": f.conf_int().loc["log_nUMI", 0], "ci95_hi": f.conf_int().loc["log_nUMI", 1]})
    pd.DataFrame(ols_rows).to_csv(f"{TAB}/R08b_depth_adjusted_module_scores.csv", index=False)

    # R09 ambient
    amb = pd.read_csv(os.path.join(AMB, "ambient_summary.csv")).set_index("gsm")
    soup = pd.read_csv(os.path.join(AMB, "soup_profiles_counts.csv"), index_col=0)
    pm = (soup.loc["MEFV"] / soup.sum()).rename("p_soup_MEFV")
    gsm_s = o.gsm.astype(str)
    o["rho_soupx"] = gsm_s.map(amb.soupx_rho_used).astype(float)
    o["lambda_ambient_MEFV"] = o.rho_soupx * o.total_counts.astype(float) * gsm_s.map(pm).astype(float)
    o["p_ambient_MEFV_ge1"] = 1 - np.exp(-o.lambda_ambient_MEFV)
    r09 = o.groupby("cell_type", observed=True).agg(
        observed_MEFV_pct=("MEFV_detected", lambda v: 100 * v.mean()),
        expected_ambient_MEFV_pct=("p_ambient_MEFV_ge1", lambda v: 100 * v.mean()),
        MEFV_pct_SoupX=("MEFV_detected_soupx", lambda v: 100 * v.mean()),
        MEFV_pct_DecontX=("MEFV_detected_decontx", lambda v: 100 * v.mean()),
        decontX_contamination_median=("decontX_contamination", "median")).reindex(CT_ORDER)
    r09.to_csv(f"{TAB}/R09a_ambient_MEFV.csv")
    amb.join(pm.rename_axis("gsm")).to_csv(f"{TAB}/R09b_ambient_by_sample.csv")
    for meth in ("soupx", "decontx"):
        t = ad.AnnData(lognorm(s.layers[meth]), obs=o[["cell_type"]].copy(), var=pd.DataFrame(index=s.var_names))
        score(t, {k: sets[k] for k in ("PYRIN_BACKBONE", "NLRP3_BACKBONE")}, prefix=f"{meth}_")
        o[f"{meth}_PYRIN_BACKBONE_score"] = t.obs[f"{meth}_PYRIN_BACKBONE_score"].values
        o[f"{meth}_NLRP3_BACKBONE_score"] = t.obs[f"{meth}_NLRP3_BACKBONE_score"].values
    r09c = o.groupby("cell_type", observed=True)[[c for c in o.columns if c.endswith("BACKBONE_score")]].mean().reindex(CT_ORDER)
    r09c["delta_cell_level_z_soupx"] = zdelta(o, ["cell_type"], "soupx_PYRIN_BACKBONE_score", "soupx_NLRP3_BACKBONE_score")
    r09c["delta_cell_level_z_decontx"] = zdelta(o, ["cell_type"], "decontx_PYRIN_BACKBONE_score", "decontx_NLRP3_BACKBONE_score")
    r09c.to_csv(f"{TAB}/R09c_ambient_corrected_scores.csv")

    # R10 random gene-set null (expression-bin matched) + leave-one-sample-out
    rng = np.random.default_rng(SEED)
    bb = [g for g in sets["PYRIN_BACKBONE"] if g in s.var_names]
    gm = np.asarray(s.X.mean(0)).ravel()
    ranks = pd.Series(gm, index=s.var_names).rank(method="min")
    bins = pd.cut(ranks, 25, labels=False)
    excl = set(g for v in M.scoring_sets("config/gene_modules.yaml").values() for g in v)
    pool = {b: [g for g in bins.index[bins == b] if g not in excl] for b in range(25)}
    obs_mean = o.groupby("cell_type", observed=True)["PYRIN_BACKBONE_score"].mean()
    null = []
    for i in range(200):
        rs = [rng.choice(pool[bins[g]]) for g in bb]
        sc.tl.score_genes(s, rs, score_name="_rnd", ctrl_size=50, n_bins=25, random_state=SEED, use_raw=False)
        null.append(s.obs.groupby("cell_type", observed=True)["_rnd"].mean())
    N = pd.DataFrame(null)
    r10 = pd.DataFrame({"observed": obs_mean, "null_mean": N.mean(), "null_sd": N.std(),
                        "empirical_percentile": [100 * (N[c] < obs_mean[c]).mean() for c in obs_mean.index],
                        "excess_over_null": obs_mean - N.mean()}).reindex(CT_ORDER)
    r10["excess_over_null_sd_units"] = r10.excess_over_null / r10.null_sd
    r10.to_csv(f"{TAB}/R10a_random_geneset_null.csv")
    loo = []
    for smp in o["sample"].unique():
        mm = o[o["sample"] != smp].groupby("cell_type", observed=True)["PYRIN_BACKBONE_score"].mean()
        loo.append({"excluded_sample": smp, "top_compartment": mm.idxmax(),
                    **{f"mean_{k}": v for k, v in mm.items()}})
    pd.DataFrame(loo).to_csv(f"{TAB}/R10b_leave_one_sample_out.csv", index=False)

    # R11 all module scores by compartment
    cols = [f"{k}_score" for k in sets]
    o.groupby("cell_type", observed=True)[cols].mean().reindex(CT_ORDER).to_csv(f"{TAB}/R11_module_scores_by_compartment.csv")
    s.obs = o.drop(columns=["_rnd"], errors="ignore")
    s.write(OUT)
    pd.set_option("display.width", 250)
    print(r04.round(2).to_string(index=False))
    print(pd.DataFrame(pooled).T.round(2)); print(comp.round(3).to_string())
    print(pd.DataFrame(glm_rows).round(3).to_string(index=False)); print(r08.round(3).to_string())
    print(r09.round(3).to_string()); print(r10.round(3).to_string())


if __name__ == "__main__":
    main()
