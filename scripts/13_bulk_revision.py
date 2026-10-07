#!/usr/bin/env python
"""Revision: paired bulk analysis of GSE120521 with
  * the former 'generic pyroptosis' module split into a cytokine arm (CASP1, IL1B, IL18)
    and a lysis arm (GSDMD, GSDME, NINJ1); caspase-4/5 reported separately,
  * GSDME recovered under its legacy symbol DFNA5,
  * t-distribution 95% CIs for the four paired differences (bootstrap CIs dropped),
  * myeloid-marker residualization applied to every module,
  * an explicit reconciliation of bulk effect sizes with each module's single-cell
    myeloid enrichment (GSE159677).

Usage (repo root): python -I scripts/13_bulk_revision.py
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "src")
from pyrinplaque import modules as M

FPKM = "data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx"
META = "data/metadata/gse120521_sample_metadata.csv"
H5AD = "data/processed/gse159677_v1repro.h5ad"
TAB = "results/revision/tables"
ALIAS = {"DFNA5": "GSDME"}
MYELOID = ["LYZ", "CD68", "CD14", "FCGR3A", "LST1", "C1QA", "C1QB", "C1QC", "MS4A7", "TYROBP"]


def module_sets():
    s = M.scoring_sets("config/gene_modules.yaml")
    sets = {
        "PYRIN_BACKBONE": s["PYRIN_BACKBONE"],
        "PYRIN_FULL": s["PYRIN_FULL"],
        "NLRP3_BACKBONE": s["NLRP3_BACKBONE"],
        "NLRP3_FULL": s["NLRP3_FULL"],
        "EFFECTOR_CYTOKINE_ARM": ["CASP1", "IL1B", "IL18"],
        "EFFECTOR_LYSIS_ARM": ["GSDMD", "GSDME", "NINJ1"],
        "NONCANONICAL_CASP4_5": ["CASP4", "CASP5"],
        "LEGACY_GENERIC_8GENE": s["GENERIC_PYROPTOSIS"],
        "MYELOID_MARKER": MYELOID,
    }
    return sets


def load_logexpr():
    df = pd.read_excel(FPKM)
    df["name"] = df["name"].replace(ALIAS)
    cols = [c for c in df.columns if "FPKM" in c]
    expr = df.dropna(subset=["name"]).groupby("name")[cols].max()
    return np.log2(expr + 1.0), cols


def module_score(logexpr, genes):
    present = [g for g in genes if g in logexpr.index]
    sub = logexpr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


def paired_summary(diff):
    diff = np.asarray(diff, float)
    n = len(diff)
    m, sd = diff.mean(), diff.std(ddof=1)
    se = sd / np.sqrt(n)
    tcrit = stats.t.ppf(0.975, n - 1)
    p = stats.ttest_1samp(diff, 0).pvalue
    return {"mean_paired_diff": m, "n_pairs_positive": int((diff > 0).sum()), "cohens_dz": m / sd if sd > 0 else np.nan,
            "t_ci95_lo": m - tcrit * se, "t_ci95_hi": m + tcrit * se, "paired_t_p": p,
            "per_pair_diffs": ";".join(f"{d:+.3f}" for d in diff)}


def myeloid_enrichment_index(sets):
    """Mean over module genes of log2 ratio of mean CP10k expression in myeloid vs all other
    QC cells (GSE159677, original QC/annotation)."""
    import anndata as ad
    a = ad.read_h5ad(H5AD, backed="r")
    myl = (a.obs["cell_type"] == "Macrophage/Myeloid").to_numpy()
    rows = {}
    for name, genes in sets.items():
        g = [x for x in genes if x in a.var_names]
        X = a[:, g].to_memory().layers["counts"]
        tot = a.obs["total_counts"].to_numpy()[:, None]
        cp = np.asarray(X.todense()) / tot * 1e4
        lr = np.log2((cp[myl].mean(0) + 0.01) / (cp[~myl].mean(0) + 0.01))
        rows[name] = {"myeloid_log2_enrichment_mean": float(np.mean(lr)), "genes_sc": ";".join(g),
                      "per_gene": ";".join(f"{x}:{v:+.2f}" for x, v in zip(g, lr))}
    return pd.DataFrame(rows).T


def main():
    sets = module_sets()
    logexpr, cols = load_logexpr()
    meta = pd.read_csv(META).set_index("column")
    scores, avail = {}, []
    for name, genes in sets.items():
        s, present = module_score(logexpr, genes)
        scores[name] = s
        for g in genes:
            avail.append({"module": name, "gene": g, "present": g in present,
                          "mean_log2FPKM": float(logexpr.loc[g].mean()) if g in present else np.nan})
    S = pd.DataFrame(scores).join(meta[["sample_id", "patient", "status", "pair_id"]])
    pd.DataFrame(avail).to_csv(f"{TAB}/R20_bulk_gene_availability.csv", index=False)

    # myeloid residualization for every module (OLS across the 8 samples, as in v1)
    x = S["MYELOID_MARKER"].to_numpy()
    rows = []
    for name in sets:
        y = S[name].to_numpy()
        piv = S.pivot_table(index="patient", columns="status", values=name)
        diff = (piv["unstable"] - piv["stable"]).to_numpy()
        row = {"module": name, "n_genes_present": int(sum(g in logexpr.index for g in sets[name])),
               **paired_summary(diff)}
        if name != "MYELOID_MARKER":
            b1, b0 = np.polyfit(x, y, 1)
            S[f"{name}_resid"] = y - (b0 + b1 * x)
            r_xy = stats.pearsonr(x, y)
            pr = S.pivot_table(index="patient", columns="status", values=f"{name}_resid")
            rd = (pr["unstable"] - pr["stable"]).to_numpy()
            rs = paired_summary(rd)
            row.update({"r_with_myeloid": r_xy.statistic, "r_p": r_xy.pvalue,
                        "resid_mean_paired_diff": rs["mean_paired_diff"],
                        "resid_n_pairs_positive": rs["n_pairs_positive"], "resid_paired_t_p": rs["paired_t_p"],
                        "resid_per_pair_diffs": rs["per_pair_diffs"],
                        "attenuation_pct": 100 * (1 - rs["mean_paired_diff"] / row["mean_paired_diff"])})
        rows.append(row)
    R = pd.DataFrame(rows)
    S.to_csv(f"{TAB}/R21_bulk_module_scores.csv")
    R.to_csv(f"{TAB}/R22_bulk_paired_results.csv", index=False)

    # reconciliation: bulk effect size vs single-cell myeloid enrichment of each module
    E = myeloid_enrichment_index({k: v for k, v in sets.items()})
    E.index.name = "module"
    rec = R.set_index("module")[["mean_paired_diff", "cohens_dz", "attenuation_pct"]].join(E)
    rec["myeloid_log2_enrichment_mean"] = rec["myeloid_log2_enrichment_mean"].astype(float)
    sub = rec.drop(index=["LEGACY_GENERIC_8GENE"])
    rho = stats.spearmanr(sub["myeloid_log2_enrichment_mean"], sub["mean_paired_diff"])
    rec.to_csv(f"{TAB}/R23_bulk_effect_vs_myeloid_enrichment.csv")
    pd.set_option("display.width", 250)
    print(R[["module", "n_genes_present", "mean_paired_diff", "n_pairs_positive", "cohens_dz", "t_ci95_lo",
             "t_ci95_hi", "paired_t_p", "r_with_myeloid", "resid_mean_paired_diff", "resid_n_pairs_positive",
             "attenuation_pct"]].round(3).to_string(index=False))
    print(rec[["mean_paired_diff", "myeloid_log2_enrichment_mean"]].round(3).to_string())
    print(f"Spearman(bulk paired diff, sc myeloid enrichment) across modules: rho={rho.statistic:.3f}, p={rho.pvalue:.3f}, n={len(sub)}")


if __name__ == "__main__":
    main()
