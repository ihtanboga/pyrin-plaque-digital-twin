# Fig 3 — PYRIN vs NLRP3 backbone violins

- **Output file:** results/figures/day3_pyrin_vs_nlrp3_violin.png
- **Code cell / function:** `violin pyrin nlrp3`
- **Data source:** GSE159677 (verified primary single-cell; Nat Commun/PMID 36224302). Integrated atlas NOT used.
- **Input:** data/processed/gse159677_scored.h5ad (49,380 cells × 23,578 genes, 6 samples / 3 patients, AC+PA)
- **Generating script:** scripts/02_score_plaque_cell_states.py + figure code in notebooks/02_score_plaque_cell_states.ipynb
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13); NUMBA_CACHE_DIR=/tmp/numba_cache
- **Random seed:** 0
- **Cell labels:** provisional, marker-based (no published per-cell labels in GEO supp) — labelled provisional in every panel.

**What it shows:** Side-by-side violins of PYRIN_BACKBONE and NLRP3_BACKBONE scores by cell type (medians marked).

**Interpretation:** PYRIN-backbone median highest in Macrophage/Myeloid; NLRP3-backbone highest in T/NK. Always reported side-by-side per carry-forward rule 3.

**Limitations:** Distributions overlap substantially; medians differ modestly. Not causal. Relates to C-A-01, C-B-01.
