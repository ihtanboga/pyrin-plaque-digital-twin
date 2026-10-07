# Revision (Cytokine, major revision) — reproducible re-analysis

This revision rebuilds the single-cell analysis from the raw GEO deposits, adds the controls requested
by the reviewers, and regenerates every table and figure from code. Manuscript numbers are filled
programmatically from `results/revision/tables` (see `tools/build_manuscript.py`).

## Corrections disclosed in the revised manuscript
* **Sample labels (GSE159677).** Aggregate barcode suffixes map to samples as given in the depositors'
  `GSM4837528_Aggregated.Sample.Meta.txt` (PA = -1/-3/-5, AC = -2/-4/-6), verified by exact barcode
  matching against each sample's `molecule_info.h5`. The original script `02_score_plaque_cell_states.py`
  assumed GSM accession order and therefore interchanged AC and PA within each patient. See
  `data/metadata/gse159677_sample_metadata.csv` (verified) and `..._ORIGINAL_ERRONEOUS.csv`.
* **Specificity delta.** Defined at the cell level (z across all cells; compartment mean of z_P − z_N);
  compartment-mean standardization is reported for comparison.
* **GSDME** is annotated as `DFNA5` in the GSE120521 FPKM matrix and is now mapped.
* **Dynamic model.** In the original model priming and threshold entered only through their difference
  (`R34_ode1_symmetry_check.csv`); replaced by the two-signal model in `src/pyrinplaque/ode_twin_v2.py`.

## Pipeline (run from this directory; Python with `-I`)
| Step | Script | Main outputs |
|---|---|---|
| 1 | `scripts/10_sc_reproduce_v1.py <filtered_matrix_dir>` | `gse159677_v1repro.h5ad`, R01 |
| 2 | `scripts/11_sc_molecule_info_soup.py <raw_dir> <work_dir>` | per-sample counts, empty-droplet profiles, R02 |
| 3 | `Rscript scripts/12_sc_doublets_ambient.R <work_dir>` | scDblFinder, SoupX, DecontX outputs |
| 4 | `scripts/14_sc_revised_analysis.py <work_dir>` | singlet set, R03–R11 |
| 5 | `scripts/15_sc_myeloid_subsets_and_regions.py` | myeloid subsets, corrected AC–PA, R12–R13 |
| 6 | `scripts/13_bulk_revision.py` | bulk split arms, residualization, R20–R23 |
| 7 | `scripts/05b_run_ode_v2.py` | two-signal model, R30–R34 |
| 8 | `scripts/16_priority_revision.py` | revised priority score and change log, R40–R42 |
| 9 | `scripts/20–24_*.py` | Figures 1–6, Supplementary Figures S1–S8, graphical abstract |
| 10 | `tools/build_manuscript.py`, `tools/build_response.py`, `tools/build_supplement.py` | manuscript (clean + tracked), response letter, supplement |

Environments: Python 3.13 (scanpy 1.12.1, anndata 0.12.19, numpy 2.4.6, scipy 1.18.0, pandas 2.3.3,
statsmodels 0.15.0, harmonypy 0.0.10, scikit-image) and R 4.5 (scDblFinder 1.24.10, SoupX 1.6.2,
celda 1.26.0). Raw data (`data/raw`) and processed objects (`data/processed`) are not committed.
