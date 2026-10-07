#!/usr/bin/env python
"""Revision step 1: reproduce the v1 single-cell pipeline exactly (same QC, clustering,
annotation and module scoring as scripts/02_score_plaque_cell_states.py), with ONE
correction: aggregate barcode suffixes are mapped to samples using the depositors'
GSM4837528_Aggregated.Sample.Meta.txt (verified by exact barcode matching against each
sample's molecule_info.h5), not by GSM accession order. v1 swapped AC and PA within
every patient.

Writes data/processed/gse159677_v1repro.h5ad and results/revision/tables/
R01_v1_reproduction_check.csv comparing the reproduced compartment summaries with the
v1 derived tables.

Usage (from repo root):
  python -I scripts/10_sc_reproduce_v1.py <filtered_feature_bc_matrix_dir>
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
import scanpy as sc
from pyrinplaque import modules as M, scoring as S

SEED = 0
np.random.seed(SEED)
MTX = sys.argv[1] if len(sys.argv) > 1 else "data/raw/GSE159677/aggregate/outs/filtered_feature_bc_matrix"
META = "data/metadata/gse159677_sample_metadata.csv"
OUT_H5AD = "data/processed/gse159677_v1repro.h5ad"
TAB = "results/revision/tables"

MARKERS = {
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "SMC/Fibroblast": ["ACTA2", "MYH11", "TAGLN", "DCN", "LUM", "COL1A1"],
    "Macrophage/Myeloid": ["CD68", "LYZ", "CD14", "AIF1", "ITGAM", "C1QA", "C1QB"],
    "T/NK": ["CD3D", "CD3E", "TRAC", "CD8A", "NKG7", "GNLY"],
    "B/Plasma": ["CD79A", "MS4A1", "CD19", "IGHG1", "MZB1"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2"]}


def annotate(adata):
    for ct, gs in MARKERS.items():
        sc.tl.score_genes(adata, [g for g in gs if g in adata.var_names],
                          score_name=f"sig_{ct}", random_state=SEED)
    cols = [f"sig_{ct}" for ct in MARKERS]
    cl_mean = adata.obs.groupby("leiden", observed=True)[cols].mean()
    assign = {cl: cols[int(np.argmax(cl_mean.loc[cl].values))].replace("sig_", "") for cl in cl_mean.index}
    adata.obs["cell_type"] = adata.obs["leiden"].map(assign).astype(str).values


def main():
    os.makedirs(TAB, exist_ok=True)
    adata = sc.read_10x_mtx(MTX, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    meta = pd.read_csv(META, dtype={"aggr_suffix": str}).set_index("aggr_suffix")
    suf = adata.obs_names.str.split("-").str[-1]
    for col in ("gsm", "patient", "region"):
        adata.obs[col] = suf.map(meta[col]).values
    adata.obs["sample"] = (adata.obs["patient"].str.replace("Patient ", "P") + "_" + adata.obs["region"]).values
    assert adata.obs["gsm"].notna().all()

    # --- v1 QC (unchanged) ---
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    n_loaded = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[(adata.obs.pct_counts_mt < 20) & (adata.obs.n_genes_by_counts < 6000)].copy()
    adata.layers["counts"] = adata.X.copy()
    cpc = np.asarray(adata.X.sum(1)).ravel()
    cpc[cpc == 0] = 1
    Xn = adata.X.multiply(1e4 / cpc[:, None]).tocsr().astype("float32")
    Xn.data = np.log1p(Xn.data)
    adata.X = Xn
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    ah = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(ah, max_value=10)
    sc.tl.pca(ah, n_comps=30, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca"] = ah.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.6, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    annotate(adata)
    sets = M.scoring_sets("config/gene_modules.yaml")
    S.score_modules(adata, sets, seed=SEED)
    mefv = np.asarray((adata.layers["counts"][:, adata.var_names.get_loc("MEFV")] > 0).todense()).ravel()
    adata.obs["MEFV_detected"] = mefv
    adata.uns["v1repro"] = {"n_loaded": int(n_loaded)}
    adata.write(OUT_H5AD)

    # --- reproduction check against v1 derived tables ---
    v1 = pd.read_csv("results/tables/cell_state_summary.csv").set_index("cell_type")
    v1d = pd.read_csv("results/tables/pyrin_specificity_delta.csv").set_index("cell_type")
    o = adata.obs
    rep = o.groupby("cell_type").agg(n_cells=("cell_type", "size"),
                                     PYRIN_BACKBONE_mean=("PYRIN_BACKBONE_score", "mean"),
                                     NLRP3_BACKBONE_mean=("NLRP3_BACKBONE_score", "mean"),
                                     GENERIC_PYROPTOSIS_mean=("GENERIC_PYROPTOSIS_score", "mean"),
                                     MEFV_pct=("MEFV_detected", lambda x: 100 * x.mean()))
    rep["pyrin_specificity_delta_cellz"] = S.specificity_delta(o)
    rows = []
    for ct in rep.index:
        for col, v1col, src in (("n_cells", "n_cells", v1), ("PYRIN_BACKBONE_mean", "PYRIN_BACKBONE_mean", v1),
                                ("NLRP3_BACKBONE_mean", "NLRP3_BACKBONE_mean", v1),
                                ("GENERIC_PYROPTOSIS_mean", "GENERIC_PYROPTOSIS_mean", v1),
                                ("MEFV_pct", "MEFV_pct", v1d),
                                ("pyrin_specificity_delta_cellz", "pyrin_specificity_delta", v1d)):
            a, b = float(rep.loc[ct, col]), float(src.loc[ct, v1col]) if ct in src.index else np.nan
            rows.append({"cell_type": ct, "metric": col, "reproduced": a, "v1_reported": b,
                         "abs_diff": abs(a - b)})
    chk = pd.DataFrame(rows)
    chk.to_csv(f"{TAB}/R01_v1_reproduction_check.csv", index=False)
    print("cells loaded", n_loaded, "| after QC", adata.n_obs, "| genes", adata.n_vars,
          "| leiden clusters", adata.obs.leiden.nunique())
    print(chk.groupby("metric").abs_diff.max())
    print(o.groupby(["sample", "cell_type"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
