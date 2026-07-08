#!/usr/bin/env python
"""Day-3: score PYRIN/NLRP3/generic pyroptosis programs in GSE159677 plaque scRNA-seq.
Runs from clean config files; no hard-coded gene lists. See reports/day3_dataset_qc_summary.md.

Usage: NUMBA_CACHE_DIR=/tmp/numba_cache python scripts/02_score_plaque_cell_states.py
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import scanpy as sc
from pyrinplaque import modules as M, scoring as S, provenance as P

SEED = 0
np.random.seed(SEED)
RAW = "data/raw/GSE159677/aggregate/outs/filtered_feature_bc_matrix"

def main():
    adata = sc.read_10x_mtx(RAW, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    # sample id from 10x-aggr barcode suffix -> GSM (aggr order)
    aggr = ["GSM4837523","GSM4837524","GSM4837525","GSM4837526","GSM4837527","GSM4837528"]
    suf = adata.obs_names.str.split("-").str[-1]
    adata.obs["gsm"] = suf.map({str(i+1): g for i, g in enumerate(aggr)}).values
    meta = pd.read_csv("data/metadata/gse159677_sample_metadata.csv").set_index("gsm")
    adata.obs["region"] = adata.obs["gsm"].map(meta["location"]).map(
        {"atherosclerotic core":"AC","proximal adjacent":"PA"}).values
    # QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[(adata.obs.pct_counts_mt < 20) & (adata.obs.n_genes_by_counts < 6000)].copy()
    adata.layers["counts"] = adata.X.copy()
    # CP10k + log1p (manual to avoid numba-parallel edge case under some builds)
    cpc = np.asarray(adata.X.sum(1)).ravel(); cpc[cpc == 0] = 1
    Xn = adata.X.multiply(1e4 / cpc[:, None]).tocsr().astype("float32")
    Xn.data = np.log1p(Xn.data); adata.X = Xn; adata.raw = adata
    # cluster + provisional annotation (published per-cell labels not in GEO supp)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    ah = adata[:, adata.var.highly_variable].copy(); sc.pp.scale(ah, max_value=10)
    sc.tl.pca(ah, n_comps=30, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca"] = ah.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.6, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    _annotate(adata)
    # module scoring from frozen config
    sets = M.scoring_sets("config/gene_modules.yaml")
    S.score_modules(adata, sets, seed=SEED)
    adata.write("data/processed/gse159677_scored.h5ad")
    prov = P.run_provenance(
        inputs=["data/raw/GSE159677/GSE159677_aggregate_filtered.tar.gz", "config/gene_modules.yaml"],
        params={"seed": SEED, "min_genes": 200, "min_cells": 3, "pct_mt_max": 20,
                "norm": "CP10k+log1p", "n_hvg": 2000, "n_pcs": 30, "leiden_res": 0.6,
                "score_method": "scanpy.score_genes ctrl_size=50 n_bins=25", "rank_score": "mean gene percentile"},
        outputs=["data/processed/gse159677_scored.h5ad"])
    P.write_provenance(prov, "audit/day3_provenance.json")
    print("done:", adata.shape)

MARKERS = {
 "Endothelial": ["PECAM1","VWF","CDH5","CLDN5"],
 "SMC/Fibroblast": ["ACTA2","MYH11","TAGLN","DCN","LUM","COL1A1"],
 "Macrophage/Myeloid": ["CD68","LYZ","CD14","AIF1","ITGAM","C1QA","C1QB"],
 "T/NK": ["CD3D","CD3E","TRAC","CD8A","NKG7","GNLY"],
 "B/Plasma": ["CD79A","MS4A1","CD19","IGHG1","MZB1"],
 "Mast": ["TPSAB1","CPA3","MS4A2"]}

def _annotate(adata):
    for ct, gs in MARKERS.items():
        sc.tl.score_genes(adata, [g for g in gs if g in adata.var_names],
                          score_name=f"sig_{ct}", random_state=SEED)
    cols = [f"sig_{ct}" for ct in MARKERS]
    cl_mean = adata.obs.groupby("leiden")[cols].mean()
    assign = {cl: cols[int(np.argmax(cl_mean.loc[cl].values))].replace("sig_","") for cl in cl_mean.index}
    adata.obs["cell_type"] = adata.obs["leiden"].map(assign).values

if __name__ == "__main__":
    main()
