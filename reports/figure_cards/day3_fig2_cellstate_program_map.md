# Fig 2 — Cell-state program map (heatmap)

- **Output file:** results/figures/fig2_cellstate_program_map.png
- **Code cell / function:** `build_fig2 heatmap cell`
- **Data source:** GSE159677 (verified primary single-cell; Nat Commun/PMID 36224302). Integrated atlas NOT used.
- **Input:** data/processed/gse159677_scored.h5ad (49,380 cells × 23,578 genes, 6 samples / 3 patients, AC+PA)
- **Generating script:** scripts/02_score_plaque_cell_states.py + figure code in notebooks/02_score_plaque_cell_states.ipynb
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13); NUMBA_CACHE_DIR=/tmp/numba_cache
- **Random seed:** 0
- **Cell labels:** provisional, marker-based (no published per-cell labels in GEO supp) — labelled provisional in every panel.

**What it shows:** Heatmap of PYRIN_BACKBONE, NLRP3_BACKBONE, GENERIC_PYROPTOSIS mean scores, z-scored across cell types.

**Interpretation:** PYRIN-backbone program score is highest in Macrophage/Myeloid cells (z=+1.73) and lowest in B/Plasma; NLRP3-backbone is highest in T/NK. Pyrin and NLRP3 programs peak in DIFFERENT cell types — a clean side-by-side separation.

**Limitations:** Program score is not causal activity. z-scoring is across only 6 cell types (small). Shared downstream genes excluded from both backbones. Relates to claims C-A-01, C-A-07, C-B-01.
