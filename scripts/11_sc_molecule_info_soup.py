#!/usr/bin/env python
"""Revision step 2: per-sample empty-droplet (ambient RNA) profiles from molecule_info.h5.

For each GSE159677 sample this script
  * verifies which aggregate barcode suffix the sample's cell barcodes map to,
  * computes full-depth UMI totals for cell barcodes (to document the aggr read subsampling),
  * builds the ambient "soup" profile from non-cell droplets with 1-100 UMIs (SoupX default),
  * writes a genes x droplets background matrix (droplets with 10-100 UMIs) for DecontX,
  * exports the QC-passing cells' aggregate counts + Leiden clusters for the R step.

Usage (repo root): python -I scripts/11_sc_molecule_info_soup.py <raw_dir> <work_dir>
"""
import os
import sys
import glob
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import anndata as ad

RAW = sys.argv[1]
WORK = sys.argv[2]
TAB = "results/revision/tables"
META = pd.read_csv("data/metadata/gse159677_sample_metadata.csv", dtype={"aggr_suffix": str})
H5AD = "data/processed/gse159677_v1repro.h5ad"


def main():
    os.makedirs(WORK, exist_ok=True)
    adata = ad.read_h5ad(H5AD)
    agg_seq = pd.Series(adata.obs_names.str.rsplit("-", n=1).str[0], index=adata.obs_names)
    agg_suf = pd.Series(adata.obs_names.str.rsplit("-", n=1).str[-1], index=adata.obs_names)
    genes_agg = pd.Index(adata.var_names)
    rows, soup_cols = [], {}
    for _, m in META.iterrows():
        path = glob.glob(os.path.join(RAW, f"{m.gsm}_*_moleculeinfo.h5"))[0]
        with h5py.File(path, "r") as f:
            bc_idx = f["barcode_idx"][:]
            ft_idx = f["feature_idx"][:]
            barcodes = f["barcodes"][:].astype(str)
            fnames = f["features/name"][:].astype(str)
            pf = f["barcode_info/pass_filter"][:, 0]
        nb, nf = len(barcodes), len(fnames)
        tot = np.bincount(bc_idx, minlength=nb)
        is_cell = np.zeros(nb, bool)
        is_cell[pf] = True
        cells = set(barcodes[pf])
        # suffix verification against the aggregate (all QC cells of that suffix)
        overlap = {s: int(agg_seq[agg_suf == s].isin(cells).sum()) for s in sorted(agg_suf.unique())}
        n_qc = int((agg_suf == m.aggr_suffix).sum())
        # full-depth UMI totals for this sample's QC cells vs aggregate (subsampled) totals
        sel = agg_suf == m.aggr_suffix
        bc2tot = pd.Series(tot[pf], index=barcodes[pf])
        full = bc2tot.reindex(agg_seq[sel].values).to_numpy(float)
        aggtot = adata.obs.loc[sel, "total_counts"].to_numpy(float)
        # soup profile (SoupX default range: droplets with <=100 UMIs that are not cells)
        soup_bc = (tot > 0) & (tot <= 100) & ~is_cell
        soup = np.bincount(ft_idx[soup_bc[bc_idx]], minlength=nf).astype(float)
        soup_s = pd.Series(soup, index=fnames).groupby(level=0).sum()
        soup_cols[m.gsm] = soup_s
        # background matrix for DecontX (droplets 10-100 UMIs)
        bg_bc = (tot >= 10) & (tot <= 100) & ~is_cell
        mol = bg_bc[bc_idx]
        cols = np.flatnonzero(bg_bc)
        colmap = np.full(nb, -1, np.int64)
        colmap[cols] = np.arange(len(cols))
        B = sp.coo_matrix((np.ones(int(mol.sum()), np.float32), (ft_idx[mol], colmap[bc_idx[mol]])),
                          shape=(nf, len(cols))).tocsr()
        B.sum_duplicates()
        # align background genes to aggregate genes (sum duplicated symbols)
        Bdf_genes = pd.Index(fnames)
        keep = Bdf_genes.isin(genes_agg)
        Bk = B[keep]
        g_k = Bdf_genes[keep]
        # collapse duplicate symbols by summing rows
        codes, uniq = pd.factorize(g_k)
        agg_map = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        Bc = (agg_map @ Bk).tocsr()
        order = pd.Index(uniq).get_indexer(genes_agg)
        Bfull = sp.vstack([Bc, sp.csr_matrix((1, Bc.shape[1]))]).tocsr()[np.where(order < 0, Bc.shape[0], order)]
        sdir = os.path.join(WORK, m.gsm)
        os.makedirs(sdir, exist_ok=True)
        sio.mmwrite(os.path.join(sdir, "background.mtx"), Bfull.astype(np.float32))
        # QC-cell aggregate counts for this sample (genes x cells) + clusters
        sub = adata[sel]
        sio.mmwrite(os.path.join(sdir, "cells_counts.mtx"), sub.layers["counts"].T.tocsr().astype(np.float32))
        pd.DataFrame({"barcode": sub.obs_names, "leiden": sub.obs["leiden"].astype(str).values,
                      "cell_type": sub.obs["cell_type"].values}).to_csv(os.path.join(sdir, "cells.csv"), index=False)
        pd.Series(genes_agg).to_csv(os.path.join(sdir, "genes.csv"), index=False, header=["gene"])
        rows.append({"gsm": m.gsm, "patient": m.patient, "region": m.region, "aggr_suffix": m.aggr_suffix,
                     "n_cell_barcodes_cellranger": int(len(pf)), "n_QC_cells_in_aggregate_suffix": n_qc,
                     "barcode_overlap_by_suffix": ";".join(f"-{k}:{v}" for k, v in overlap.items()),
                     "median_UMI_full_depth": float(np.nanmedian(full)), "median_UMI_aggregate": float(np.median(aggtot)),
                     "aggregate_to_full_UMI_ratio_median": float(np.nanmedian(aggtot / full)),
                     "n_soup_droplets_1_100": int(soup_bc.sum()), "soup_UMIs": float(soup.sum()),
                     "n_background_droplets_10_100": int(len(cols)),
                     "MEFV_soup_fraction": float(soup_s.get("MEFV", 0) / soup_s.sum())})
        print(rows[-1], flush=True)
        del bc_idx, ft_idx
    pd.DataFrame(rows).to_csv(f"{TAB}/R02_sample_barcode_verification_and_soup.csv", index=False)
    soup_df = pd.DataFrame(soup_cols)
    soup_df.to_csv(os.path.join(WORK, "soup_profiles_counts.csv"))


if __name__ == "__main__":
    main()
