# Target Card 5: PKN1/PKN2 / 14-3-3 inhibitory regulation

- **Candidate ID:** L2-03  ·  **Type:** mechanism_node  ·  **Axis:** pyrin inhibitory regulation
- **Final priority score:** **29.75**  (raw positive 55.75 − penalties 26.0)  ·  **Confidence grade: low**
- **Rank:** 5 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
PKN1/PKN2 phosphorylate pyrin; 14-3-3 (YWHAB/YWHAZ) binding keeps pyrin inhibited — a well-defined pyrin-specific ON/OFF switch mapping to the ODE inhibition/threshold levers.

## Strongest counter-evidence / penalties
Low MEFV context; ODE shows the inhibition lever is secondary to priming; no direct plaque perturbation data.

## Positive score components
| Component | Points |
|-----------|--------|
| pyrin_specificity_score | +16.0 |
| literature_support_score | +11.2 |
| singlecell_compartment_score | +7.5 |
| shared_downstream_exclusion_pass_score | +5.0 |
| druggability_or_intervention_plausibility_score | +5.0 |
| bulk_directional_support_score | +4.0 |
| ODE_leverage_score | +3.0 |
| MEFV_detection_score | +2.0 |
| clinical_alignment_score | +2.0 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| myeloid_confounding_penalty | -4.5 |
| random_control_penalty | -3.5 |
| no_genotype_penalty | -2.8 |
| sparse_expression_penalty | -2.4 |
| ac_pa_mixed_signal_penalty | -2.4 |
| provisional_annotation_penalty | -2.0 |
| ODE_uncalibrated_penalty | -1.8 |
| shared_downstream_non_specificity_penalty | -0.6 |

## Druggability / clinical alignment
Representative target **PKN1**: ChEMBL target found=1, Open Targets tractable buckets=10. Clinical: no direct clinical agent identified. (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
PKN inhibitor / 14-3-3 disruption in myeloid pyrin assays; map to ODE inhibition_factor.

## Safe claim (approved language)
> The PKN/14-3-3 inhibitory switch is a mechanistically specific pyrin regulatory node suited to perturbation experiments.

## Claims to avoid (banned)
druggable in plaque; controls instability

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
