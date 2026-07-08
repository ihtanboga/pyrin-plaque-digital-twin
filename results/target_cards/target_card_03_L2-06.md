# Target Card 3: NLRP3 (saturated comparator)

- **Candidate ID:** L2-06  ·  **Type:** mechanism_node  ·  **Axis:** NLRP3 comparator
- **Final priority score:** **32.0**  (raw positive 51.0 − penalties 19.0)  ·  **Confidence grade: low**
- **Rank:** 3 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
NLRP3-in-atherosclerosis is strongly supported (CANTOS, colchicine trials, extensive mechanism).

## Strongest counter-evidence / penalties
Saturated, well-studied comparator — NOT novel, NOT pyrin-specific; included only to contextualize pyrin.

## Positive score components
| Component | Points |
|-----------|--------|
| literature_support_score | +13.5 |
| singlecell_compartment_score | +9.0 |
| clinical_alignment_score | +9.0 |
| druggability_or_intervention_plausibility_score | +7.0 |
| bulk_directional_support_score | +6.0 |
| shared_downstream_exclusion_pass_score | +5.0 |
| ODE_leverage_score | +1.5 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| myeloid_confounding_penalty | -2.7 |
| provisional_annotation_penalty | -2.0 |
| no_genotype_penalty | -2.0 |
| ac_pa_mixed_signal_penalty | -1.6 |
| random_control_penalty | -1.5 |
| shared_downstream_non_specificity_penalty | -1.5 |
| ODE_uncalibrated_penalty | -0.9 |
| sparse_expression_penalty | -0.8 |

## Druggability / clinical alignment
Representative target **NLRP3**: ChEMBL target found=1, Open Targets tractable buckets=0. Clinical: colchicine (approved; COLCOT NCT02551094 n=4745, LoDoCo2 NCT03048825 n=7264 — NLRP3-associated, non-specific); MCC950 preclinical SM inhibitor. (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
None as novelty; comparator baseline only.

## Safe claim (approved language)
> NLRP3 is a well-validated inflammasome comparator in atherosclerosis and frames what pyrin-specific work must add beyond.

## Claims to avoid (banned)
novel; pyrin-specific; first report

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
