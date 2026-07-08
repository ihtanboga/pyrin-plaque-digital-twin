# Day 5 — Caveat Ledger

Every caveat from Days 1–4 and exactly how it constrains Day-5 scoring. Nothing here is hidden; each maps to a named penalty in `config/scoring_weights.yaml` and a line in every affected target card.

| ID | Day | Layer | Severity | Effect on priority | Representation |
|----|-----|-------|----------|--------------------|----------------|
| CV01 | Day3 | single-cell | **high** | confidence downgrade | MEFV detection sparse (~0.97% overall, 4.3% myeloid); MEFV_detection_score capped, sparse_expression_penalty applied to any MEFV-based claim |
| CV02 | Day3 | single-cell | **medium** | penalty | Provisional marker-based cell labels (no published per-cell labels); provisional_annotation_penalty on all compartment claims |
| CV03 | Day3 | single-cell | **high** | penalty | Small n=3 scRNA-seq patients; small_n_penalty; no patient-level inference |
| CV04 | Day3 | single-cell | **medium** | penalty | Random matched gene-set control modest (~82nd pctile, p≈0.18); random_control_penalty; signal 'suggestive not definitive' |
| CV05 | Day3 | single-cell | **medium** | confidence downgrade | Small Leiden clusters (leiden16 n=138, leiden4 n=160) low donor support; low-confidence grade for substate claims |
| CV06 | Day3.5 | AC-PA bridge | **high** | penalty | AC-vs-PA mixed: MEFV up in core 3/3 but PYRIN_BACKBONE not consistently up (1/3), specificity delta DOWN in AC (NLRP3 rises in core); ac_pa_mixed_signal_penalty strong on any 'plaque core' claim |
| CV07 | Day4 | bulk | **high** | penalty | Bulk myeloid confounding: PYRIN_BACKBONE unstable>stable 4/4 (d=1.44) collapses to +0.14 (p≈0.69) after myeloid residualization; myeloid_confounding_penalty strong on all bulk-derived support |
| CV08 | Day4 | bulk | **medium** | penalty | GSE120521 n=4 patients, FPKM-only (no raw counts, no count-based modeling); small_n_penalty on bulk |
| CV09 | Day1-5 | all | **high** | exclusion | No MEFV genotype / FMF stratification; cannot link expression to pyrin genetic status; no_genotype_penalty; exclude any genotype-conditioned claim |
| CV10 | Day4 | ODE | **high** | confidence downgrade | ODE dimensionless, uncalibrated, hypothesis-generating; ODE_uncalibrated_penalty; ODE_leverage_score contributes only qualitative ordering, never magnitude |
| CV11 | Day1-2 | literature/all | **medium** | exclusion | Shared downstream PYCARD/CASP1/GSDMD/IL1B/IL18 non-discriminative; shared_downstream_non_specificity_penalty; these nodes get ZERO pyrin-specificity credit |
| CV12 | Day2 | literature | **low** | none | Integrated atlas (PMID 40931012) is literature precedent only, NOT a data source; no atlas-dependent data claim |

## Scoring rules that follow directly from this ledger
1. **Bulk support is discounted hard.** Any candidate whose support is bulk-derived (unstable-plaque pyrin signal) carries `myeloid_confounding_penalty` (CV07) — the raw 4/4 effect is largely myeloid-abundance-driven.
2. **'Plaque core' claims are discounted.** AC-vs-PA is mixed (CV06): MEFV detection rises in core but the composite pyrin program does not, and specificity delta falls. Any core-localization claim carries `ac_pa_mixed_signal_penalty`.
3. **Shared downstream genes get zero pyrin-specificity credit** (CV11) — they score only as shared pathway-output comparators.
4. **NLRP3 is a saturated comparator** — it may score high biologically but is labeled 'not novel / saturated comparator' and never as a pyrin-specific novelty.
5. **MEFV-based claims are confidence-downgraded** for sparsity (CV01) and lack of genotype (CV09).
6. **ODE contributes qualitative leverage only** (CV10) — never magnitude, never clinical prediction.

## Non-negotiable language bans (carry-forward rule 10)
No 'first ever'; no causality; no clinical prediction; no 'MEFV causes plaque instability'; no 'colchicine acts via pyrin here'; no 'validated target'.