# Target Card 1: MEFV / pyrin threshold axis

- **Candidate ID:** L2-01  ·  **Type:** mechanism_node  ·  **Axis:** pyrin sensor
- **Final priority score:** **37.9**  (raw positive 68.5 − penalties 30.6)  ·  **Confidence grade: moderate**
- **Rank:** 1 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
Strong pyrin-specific literature (PYD/B30.2 architecture, PSTPIP1, threshold biology); MEFV is the defining pyrin sensor; survives shared-downstream exclusion.

## Strongest counter-evidence / penalties
MEFV sparse in scRNA-seq; ODE sensitivity shows generic priming outweighs the pyrin threshold lever; no genotype data.

## Positive score components
| Component | Points |
|-----------|--------|
| pyrin_specificity_score | +19.0 |
| literature_support_score | +13.5 |
| singlecell_compartment_score | +10.5 |
| MEFV_detection_score | +7.0 |
| shared_downstream_exclusion_pass_score | +5.0 |
| bulk_directional_support_score | +5.0 |
| clinical_alignment_score | +4.0 |
| ODE_leverage_score | +2.5 |
| druggability_or_intervention_plausibility_score | +2.0 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| sparse_expression_penalty | -5.6 |
| myeloid_confounding_penalty | -4.5 |
| no_genotype_penalty | -4.0 |
| random_control_penalty | -3.5 |
| ac_pa_mixed_signal_penalty | -3.2 |
| provisional_annotation_penalty | -2.0 |
| ODE_uncalibrated_penalty | -1.8 |

## Druggability / clinical alignment
Representative target **MEFV**: ChEMBL target found=1, Open Targets tractable buckets=4. Clinical: colchicine (FMF standard of care — pyrin pathway; NOT proven mechanism in plaque). (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
MEFV genotype-stratified plaque cohort; pyrin-threshold perturbation in myeloid models.

## Safe claim (approved language)
> The MEFV/pyrin threshold axis is the most pyrin-specific mechanistic node and a priority for hypothesis-driven validation.

## Claims to avoid (banned)
MEFV causes atherosclerosis; threshold proven; druggable target

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
