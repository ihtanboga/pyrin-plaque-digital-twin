# Fig — Module score UMAP

- **Output file:** results/figures/module_score_umap.png
- **Code cell / function:** `UMAP scatter module score`
- **Data source:** GSE159677 (verified primary single-cell; Nat Commun/PMID 36224302). Integrated atlas NOT used.
- **Input:** data/processed/gse159677_scored.h5ad (49,380 cells × 23,578 genes, 6 samples / 3 patients, AC+PA)
- **Generating script:** scripts/02_score_plaque_cell_states.py + figure code in notebooks/02_score_plaque_cell_states.ipynb
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13); NUMBA_CACHE_DIR=/tmp/numba_cache
- **Random seed:** 0
- **Cell labels:** provisional, marker-based (no published per-cell labels in GEO supp) — labelled provisional in every panel.

**What it shows:** UMAP (computed here; not published) colored by cell type and by PYRIN_BACKBONE, NLRP3_BACKBONE, GENERIC_PYROPTOSIS per-cell scores.

**Interpretation:** Generic pyroptosis and pyrin-backbone scores are visibly enriched in the myeloid compartment; NLRP3-backbone is more diffuse with T/NK weight.

**Limitations:** UMAP embedding is provisional (computed from this run's PCA, seed 0). Per-cell scores are noisy; interpret compartment-level, not per-cell. Relates to C-B-03, C-B-04.
