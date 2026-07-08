# Fig 4 — Inflammasome gene dotplot

- **Output file:** results/figures/day3_gene_dotplot.png
- **Code cell / function:** `dotplot genes cell type`
- **Data source:** GSE159677 (verified primary single-cell; Nat Commun/PMID 36224302). Integrated atlas NOT used.
- **Input:** data/processed/gse159677_scored.h5ad (49,380 cells × 23,578 genes, 6 samples / 3 patients, AC+PA)
- **Generating script:** scripts/02_score_plaque_cell_states.py + figure code in notebooks/02_score_plaque_cell_states.ipynb
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13); NUMBA_CACHE_DIR=/tmp/numba_cache
- **Random seed:** 0
- **Cell labels:** provisional, marker-based (no published per-cell labels in GEO supp) — labelled provisional in every panel.

**What it shows:** Dotplot of pyrin-backbone, NLRP3-backbone and shared-downstream genes by cell type (dot size = % expressing, color = scaled mean).

**Interpretation:** MEFV is sparse across all cell types (largest/darkest dot in Macrophage/Myeloid). Shared downstream genes (PYCARD/CASP1/IL1B) are broadly expressed in myeloid cells; the dashed dividers mark the shared block, which is excluded from discriminative scoring.

**Limitations:** MEFV dropout limits interpretation; shown for transparency. Relates to C-A-05, C-A-08, C-C shared-node claims.
