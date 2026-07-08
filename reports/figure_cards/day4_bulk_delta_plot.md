# Fig — Pyrin-minus-NLRP3 backbone delta

- **Output file:** results/figures/day4_bulk_delta_plot.png
- **Data source:** GSE120521 (verified public bulk carotid plaque RNA-seq, PMID pending; 4 patients × stable/unstable, FPKM). Triangulation only.
- **Input:** data/raw/GSE120521/GSE120521_Athero_RNAseq_FPKM.xlsx (log2(FPKM+1))
- **Generating script:** scripts/03_bulk_plaque_validation.py
- **Environment:** conda env `sc` (scanpy 1.12.1, pandas 2.3.3, python 3.13.14); seed 0
- **Scoring:** mean per-gene z-score across samples, from frozen config/gene_modules.yaml (no hard-coded lists); shared downstream genes excluded from backbones.

**What it shows:** PYRIN_BACKBONE minus NLRP3_BACKBONE score, paired stable→unstable per patient.

**Interpretation:** The pyrin-minus-NLRP3 delta increases in unstable regions in 4/4 patients (mean Δ=+0.364) — pyrin-backbone rises disproportionately versus NLRP3-backbone with instability.

**Limitations:** n=4, descriptive; does not establish pyrin-specific mechanism; shared downstream excluded from both backbones.

**Interpretation rule:** Bulk data is triangulation only — it cannot prove cell-type localization, mechanism, or causality.
