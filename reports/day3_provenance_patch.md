# Day 3.5 — Provenance Patch

**Why:** Day-3 `audit/day3_provenance.json` was written before scanpy was imported in that kernel, so its `packages` field was empty. Claude Science artifacts must be defensible by exact code + environment + provenance, so this patch backfills package versions and records code/output checksums.

## 1. Environment (conda env `sc`)
| Package | Version |
|---------|---------|
| python | 3.13.14 |
| scanpy | 1.12.1 |
| anndata | 0.12.19 |
| pandas | 2.3.3 |
| numpy | 2.4.6 |
| scipy | 1.18.0 |
| matplotlib | 3.11.0 |
| seaborn | 0.13.2 |
| scikit-learn | 1.9.0 |
| numba | 0.65.1 |
| leidenalg | 0.12.0 |
| python-igraph | 1.0.0 |
| umap-learn | 0.5.12 |
| pyarrow | 24.0.0 |
| pyyaml | 6.0.3 |

Full freeze: `audit/day3_environment_lock.txt` (62 pkgs). Pinned key deps: `requirements-freeze.txt`.

## 2. NUMBA_CACHE_DIR fix (kept documented)
The sandbox cannot write `__pycache__` into the read-only conda `site-packages`, so scanpy's numba-cached kernels raise `cannot cache function … no locator available`. **Fix:** set `NUMBA_CACHE_DIR=/tmp/numba_cache` **before importing scanpy** (JIT stays enabled). This is applied at the top of `scripts/02_score_plaque_cell_states.py` and the notebook, and must be re-applied on every fresh kernel in env `sc`.

## 3. Script checksums (sha256)
- `scripts/02_score_plaque_cell_states.py` — `8325332bbd4a08dc…`
- `src/pyrinplaque/modules.py` — `deb97076e5fd3fcb…`
- `src/pyrinplaque/scoring.py` — `fcd1d91551a3e38e…`
- `src/pyrinplaque/plotting.py` — `ff7479e46d60dd0d…`
- `src/pyrinplaque/provenance.py` — `1fad5f366cdfe36a…`

## 4. Figure checksums (sha256)
- `fig2_cellstate_program_map.png` — `0d6adf028ddc8ec6…`
- `module_score_umap.png` — `d2795b1245a377dd…`
- `day3_pyrin_vs_nlrp3_violin.png` — `1deaf5728c4f3a96…`
- `day3_gene_dotplot.png` — `181a2263b4f792f6…`
- `day3_mefv_detection.png` — `8193b9f1d3cbf130…`

## 5. Table checksums (sha256)
- `cell_state_scores.parquet` — `3f7ecc815cd91e90…`
- `cell_state_summary.csv` — `a7472737cbc805fa…`
- `pyrin_specificity_delta.csv` — `72c2d5b20855460c…`
- `day3_gene_availability.csv` — `f80cd6871ec20486…`
- `day3_gene_detection_by_celltype.csv` — `922367c2cdd3b50b…`
- `day3_random_geneset_control.csv` — `624f55e6fe90dfd7…`
- `day3_donor_influence.csv` — `2cd2aad39a98bf86…`
- `day3_shared_downstream_exclusion.csv` — `0ddc8961f2dd66ef…`
- `day3_cellstate_ranking.csv` — `1a0aff261d427376…`

## 6. Data provenance (unchanged from Day 3)
- Input: `GSE159677_aggregate_filtered.tar.gz` (sha256 `8cd0b7f84725e0c8…`), `config/gene_modules.yaml` (`f1078086783422bd…`).
- Atlas NOT used (`atlas_used=False`).
- Seed 0; CP10k+log1p; HVG 2000; 30 PCs; Leiden res 0.6; scanpy `score_genes` ctrl_size=50 n_bins=25.
