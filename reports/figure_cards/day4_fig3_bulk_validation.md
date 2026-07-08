# Fig 3 — Bulk validation (paired program scores by status)

- **Output file:** results/figures/fig3_bulk_validation.png
- **Data source:** GSE120521 (verified public bulk carotid plaque RNA-seq, PMID pending; 4 patients × stable/unstable, FPKM). Triangulation only.
- **Input:** data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx (log2(FPKM+1))
- **Generating script:** scripts/03_bulk_plaque_validation.py
- **Environment:** conda env `sc` (scanpy 1.12.1, pandas 2.3.3, python 3.13.14); seed 0
- **Scoring:** mean per-gene z-score across samples, from frozen config/gene_modules.yaml (no hard-coded lists); shared downstream genes excluded from backbones.

**What it shows:** Paired within-patient (n=4) PYRIN_BACKBONE, NLRP3_BACKBONE, GENERIC_PYROPTOSIS in stable vs unstable plaque regions.

**Interpretation:** All three programs are higher in unstable regions (PYRIN_BACKBONE 4/4 patients, paired Cohen's d=1.44; NLRP3_BACKBONE 3/4, d=0.86; GENERIC 4/4, d=1.84). Directionally consistent with the Day-3 hypothesis that pyrin-associated programs rise with plaque instability.

**Limitations:** n=4 patients; FPKM only; bulk (no cell-type resolution); the pyrin rise is NOT independent of myeloid abundance (see confounding figure). p-values secondary and underpowered.

**Interpretation rule:** Bulk data is triangulation only — it cannot prove cell-type localization, mechanism, or causality.
