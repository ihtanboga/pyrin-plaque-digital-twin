# Target Card 2: Macrophage/Myeloid pyrin-permissive state

- **Candidate ID:** L1-01  ·  **Type:** cell_state  ·  **Axis:** pyrin cell-state localization
- **Final priority score:** **33.1**  (raw positive 66.0 − penalties 32.9)  ·  **Confidence grade: low**
- **Rank:** 2 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
Day-3: PYRIN_BACKBONE highest in Macrophage/Myeloid (specificity delta +0.74); survives shared-downstream exclusion (myeloid #1 in FULL and BACKBONE); consistent across all 6 donors.

## Strongest counter-evidence / penalties
Random matched gene-set control modest (82nd pctile, p≈0.18); provisional labels; n=3 patients; bulk support confounded by myeloid abundance.

## Positive score components
| Component | Points |
|-----------|--------|
| pyrin_specificity_score | +15.0 |
| singlecell_compartment_score | +13.5 |
| literature_support_score | +10.5 |
| bulk_directional_support_score | +7.0 |
| MEFV_detection_score | +6.0 |
| shared_downstream_exclusion_pass_score | +5.0 |
| clinical_alignment_score | +4.0 |
| druggability_or_intervention_plausibility_score | +3.0 |
| ODE_leverage_score | +2.0 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| sparse_expression_penalty | -5.6 |
| myeloid_confounding_penalty | -5.4 |
| provisional_annotation_penalty | -4.0 |
| random_control_penalty | -4.0 |
| no_genotype_penalty | -4.0 |
| ac_pa_mixed_signal_penalty | -2.4 |
| ODE_uncalibrated_penalty | -1.5 |

## Druggability / clinical alignment
Representative target **MEFV**: ChEMBL target found=1, Open Targets tractable buckets=4. Clinical: colchicine (FMF standard of care — pyrin pathway; NOT proven mechanism in plaque). (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
Targeted MEFV/pyrin RNA-FISH or protein staining in myeloid cells across more donors; genotype-stratified plaque cohort.

## Safe claim (approved language)
> Human plaque myeloid compartments show a pyrin-regulatory permissiveness pattern that may justify targeted experimental validation of the MEFV/pyrin threshold axis.

## Claims to avoid (banned)
causal pyrin activation; MEFV causes plaque instability; validated target

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
