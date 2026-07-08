# Fig — Bulk module-gene heatmap

- **Output file:** results/figures/day4_bulk_gene_heatmap.png
- **Data source:** GSE120521 (verified public bulk carotid plaque RNA-seq, PMID pending; 4 patients × stable/unstable, FPKM). Triangulation only.
- **Input:** data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx (log2(FPKM+1))
- **Generating script:** scripts/03_bulk_plaque_validation.py
- **Environment:** conda env `sc` (scanpy 1.12.1, pandas 2.3.3, python 3.13.14); seed 0
- **Scoring:** mean per-gene z-score across samples, from frozen config/gene_modules.yaml (no hard-coded lists); shared downstream genes excluded from backbones.

**What it shows:** Row z-scored log2(FPKM+1) of PYRIN-backbone, NLRP3-backbone and shared-downstream genes across all 8 samples (ordered stable then unstable).

**Interpretation:** Unstable samples (right) show broadly higher expression across most inflammasome genes; MEFV is modestly higher in most unstable samples.

**Limitations:** n=4; per-gene noise high at this sample size; visual, not inferential.

**Interpretation rule:** Bulk data is triangulation only — it cannot prove cell-type localization, mechanism, or causality.
