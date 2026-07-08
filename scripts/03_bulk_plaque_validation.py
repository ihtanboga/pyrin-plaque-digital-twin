#!/usr/bin/env python
"""Day-4 Part A: bulk plaque validation on GSE120521 (paired stable vs unstable).
Triangulation only — cannot prove cell-type localization, mechanism, or causality.
Config-driven module scoring; shared downstream genes excluded from backbones.

Usage: python scripts/03_bulk_plaque_validation.py
"""
import sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "src")
from pyrinplaque import bulk as B, modules as M

FPKM = "data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx"
META = "data/metadata/gse120521_sample_metadata.csv"
MYELOID = ["LYZ","CD68","CD14","FCGR3A","LST1","C1QA","C1QB","C1QC","MS4A7","TYROBP"]

def main():
    sets = M.scoring_sets("config/gene_modules.yaml")
    sets["MYELOID_MARKER"] = MYELOID
    expr, scols = B.load_fpkm(FPKM)
    logexpr = B.log_transform(expr)
    scores, detected = B.score_all(logexpr, sets)
    meta = pd.read_csv(META).set_index("column")
    scores = scores.join(meta[["sample_id","patient","status","pair_id"]])
    scores.to_csv("results/tables/bulk_module_scores.csv")

    # paired stable-vs-unstable
    rows = []
    for m in [c for c in sets]:
        piv = scores.pivot_table(index="patient", columns="status", values=m)
        diff = (piv["unstable"] - piv["stable"]).values
        t, p = stats.ttest_rel(piv["unstable"], piv["stable"])
        rows.append({"module": m, "mean_stable": piv["stable"].mean(), "mean_unstable": piv["unstable"].mean(),
                     "mean_paired_diff": diff.mean(), "n_unstable_higher": int((diff > 0).sum()),
                     "cohens_d_paired": B.cohens_d_paired(diff), "paired_t_p": float(p)})
    pd.DataFrame(rows).to_csv("results/tables/bulk_status_tests.csv", index=False)

    # myeloid confounding: residualize PYRIN_BACKBONE on MYELOID_MARKER
    x, y = scores["MYELOID_MARKER"].values, scores["PYRIN_BACKBONE"].values
    b1, b0 = np.polyfit(x, y, 1)
    scores["PYRIN_BACKBONE_resid"] = y - (b0 + b1 * x)
    print("bulk validation done:", scores.shape)

if __name__ == "__main__":
    main()
