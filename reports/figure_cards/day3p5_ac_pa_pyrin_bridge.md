# Figure card — Day 3.5 AC-vs-PA pyrin bridge

- **Output:** results/figures/day3p5_ac_pa_pyrin_bridge.png
- **Input:** data/processed/gse159677_scored.h5ad (existing scored object; no re-analysis)
- **Code:** Day-3.5 aggregation cell (paired within-patient); tables day3p5_ac_pa_compartment_scores.csv, day3p5_ac_pa_pairwise_tests.csv
- **Environment:** conda env `sc` (scanpy 1.12.1, python 3.13.14); seed 0
- **Data source:** GSE159677 (verified primary scRNA-seq); atlas NOT used
- **Shows:** Paired within-patient (n=3) Macrophage/Myeloid PYRIN-backbone score, pyrin specificity delta, and MEFV detection, PA adjacent → AC core.
- **Interpretation:** MEFV detection rises in AC core in all 3 patients; composite PYRIN_BACKBONE score and pyrin_specificity_delta do not (delta falls because NLRP3_BACKBONE rises in core).
- **Limitations:** n=3 patients, descriptive only (no p-values); provisional cell labels; MEFV sparse/dropout-prone; internal bridge, NOT independent validation.
- **Relates to:** Day-3 myeloid pyrin-permissiveness finding; claim cards C-A-05, C-A-08 (MEFV/myeloid).
