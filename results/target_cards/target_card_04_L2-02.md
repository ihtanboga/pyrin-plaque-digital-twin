# Target Card 4: RHOA-RAC1-CDC42 regulatory context

- **Candidate ID:** L2-02  ·  **Type:** mechanism_node  ·  **Axis:** pyrin upstream regulation
- **Final priority score:** **30.6**  (raw positive 55.5 − penalties 24.9)  ·  **Confidence grade: low**
- **Rank:** 4 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
RhoA-pyrin regulation is well-supported pyrin-specific biology (RhoA inactivation activates pyrin); part of the backbone; broadly expressed so not dropout-limited.

## Strongest counter-evidence / penalties
RHOA/RAC1/CDC42 are pleiotropic housekeeping GTPases — expression not pyrin-specific; RAC1/CDC42 not directly cited in top pyrin claims.

## Positive score components
| Component | Points |
|-----------|--------|
| pyrin_specificity_score | +14.0 |
| literature_support_score | +10.5 |
| singlecell_compartment_score | +9.0 |
| shared_downstream_exclusion_pass_score | +5.0 |
| bulk_directional_support_score | +5.0 |
| druggability_or_intervention_plausibility_score | +4.0 |
| MEFV_detection_score | +3.0 |
| clinical_alignment_score | +3.0 |
| ODE_leverage_score | +2.0 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| myeloid_confounding_penalty | -4.5 |
| random_control_penalty | -3.5 |
| no_genotype_penalty | -2.8 |
| ac_pa_mixed_signal_penalty | -2.4 |
| provisional_annotation_penalty | -2.0 |
| sparse_expression_penalty | -1.6 |
| ODE_uncalibrated_penalty | -1.5 |
| shared_downstream_non_specificity_penalty | -0.6 |

## Druggability / clinical alignment
Representative target **RHOA**: ChEMBL target found=1, Open Targets tractable buckets=9. Clinical: no direct clinical agent identified. (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
RhoA-modulation in myeloid pyrin reporter assays.

## Safe claim (approved language)
> The RhoA/Rac/Cdc42 axis is a pyrin-regulatory context worth probing, acknowledging these GTPases are pleiotropic.

## Claims to avoid (banned)
pyrin-specific expression marker; causal in plaque

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
