# Fig — Myeloid confounding check

- **Output file:** results/figures/day4_bulk_myeloid_confounding.png
- **Data source:** GSE120521 (verified public bulk carotid plaque RNA-seq, PMID pending; 4 patients × stable/unstable, FPKM). Triangulation only.
- **Input:** data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx (log2(FPKM+1))
- **Generating script:** scripts/03_bulk_plaque_validation.py
- **Environment:** conda env `sc` (scanpy 1.12.1, pandas 2.3.3, python 3.13.14); seed 0
- **Scoring:** mean per-gene z-score across samples, from frozen config/gene_modules.yaml (no hard-coded lists); shared downstream genes excluded from backbones.

**What it shows:** PYRIN_BACKBONE vs MYELOID_MARKER score (canonical myeloid markers LYZ,CD68,CD14,FCGR3A,LST1,C1QA/B/C,MS4A7,TYROBP), colored by status.

**Interpretation:** PYRIN_BACKBONE correlates with myeloid abundance (Pearson r=0.675, p=0.066). After exploratory residualization on MYELOID_MARKER, the stable-vs-unstable pyrin difference collapses (mean paired residual diff +0.140, 3/4, paired t p=0.69). Much of the bulk pyrin-backbone rise is explained by increased myeloid content in unstable plaque, NOT necessarily per-cell pyrin biology.

**Limitations:** n=8 samples; residualization is exploratory and underpowered; confounding is real and stated honestly.

**Interpretation rule:** Bulk data is triangulation only — it cannot prove cell-type localization, mechanism, or causality.
