# Target Card 6: T/NK NLRP3-backbone-high compartment (comparator)

- **Candidate ID:** L1-03  ·  **Type:** cell_state  ·  **Axis:** NLRP3 comparator
- **Final priority score:** **29.2**  (raw positive 47.0 − penalties 17.8)  ·  **Confidence grade: low**
- **Rank:** 6 of 11 candidates

> **This is a hypothesis-prioritization card for future validation. It is NOT a validated target, NOT a drug recommendation, and NOT a causal or clinical claim.**

## Strongest supporting evidence
Day-3: NLRP3_BACKBONE highest in T/NK; NLRP3-in-atherosclerosis is well-established (CANTOS etc.).

## Strongest counter-evidence / penalties
NLRP3-in-atherosclerosis is a saturated, well-studied comparator — NOT novel and NOT pyrin-specific.

## Positive score components
| Component | Points |
|-----------|--------|
| literature_support_score | +12.0 |
| singlecell_compartment_score | +12.0 |
| clinical_alignment_score | +6.0 |
| druggability_or_intervention_plausibility_score | +6.0 |
| shared_downstream_exclusion_pass_score | +5.0 |
| bulk_directional_support_score | +5.0 |
| ODE_leverage_score | +1.0 |

## Penalty components (caveats applied)
| Penalty | Points |
|---------|--------|
| small_n_penalty | -6.0 |
| provisional_annotation_penalty | -4.0 |
| no_genotype_penalty | -2.0 |
| myeloid_confounding_penalty | -1.8 |
| random_control_penalty | -1.5 |
| ODE_uncalibrated_penalty | -0.9 |
| sparse_expression_penalty | -0.8 |
| ac_pa_mixed_signal_penalty | -0.8 |

## Druggability / clinical alignment
Representative target **NLRP3**: ChEMBL target found=1, Open Targets tractable buckets=0. Clinical: colchicine (approved; COLCOT NCT02551094 n=4745, LoDoCo2 NCT03048825 n=7264 — NLRP3-associated, non-specific); MCC950 preclinical SM inhibitor. (Tractability ≠ pyrin-selective druggability.)

## Recommended next validation
None as novelty; use only as a comparator baseline.

## Safe claim (approved language)
> NLRP3-backbone localizes to a distinct (T/NK-weighted) compartment and serves as a saturated comparator to pyrin.

## Claims to avoid (banned)
novel finding; pyrin-specific; first report

## Caveats carried (from caveat ledger)
See reports/day5_caveat_ledger.md. Relevant: sparse MEFV (CV01), provisional labels (CV02), small n (CV03/CV08), random-control modest (CV04), AC/PA mixed (CV06), myeloid confounding (CV07), no genotype (CV09), ODE uncalibrated (CV10), shared-downstream non-specific (CV11).
