# Fig 5 — MEFV detection

- **Output file:** results/figures/day3_mefv_detection.png
- **Code cell / function:** `mefv detection bar`
- **Data source:** GSE159677 (verified primary single-cell; Nat Commun/PMID 36224302). Integrated atlas NOT used.
- **Input:** data/processed/gse159677_scored.h5ad (49,380 cells × 23,578 genes, 6 samples / 3 patients, AC+PA)
- **Generating script:** scripts/02_score_plaque_cell_states.py + figure code in notebooks/02_score_plaque_cell_states.ipynb
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13); NUMBA_CACHE_DIR=/tmp/numba_cache
- **Random seed:** 0
- **Cell labels:** provisional, marker-based (no published per-cell labels in GEO supp) — labelled provisional in every panel.

**What it shows:** MEFV detection rate by cell type (left) and by sample (right; red=atherosclerotic core, blue=proximal adjacent).

**Interpretation:** MEFV detected in only 0.97% of cells overall, concentrated in Macrophage/Myeloid (4.3%). Reported up front per carry-forward rule 4. No single-donor drives the signal (donor-influence control: all flags False).

**Limitations:** MEFV is dropout-prone in 3' scRNA-seq; low detection is expected and does NOT imply absence of pyrin biology. Patient 3 shows modestly higher detection but does not drive cell-type results. Relates to C-A-05.
