# Day 5 — Druggability & Clinical-Alignment Note

**Purpose:** tractability/clinical context for pyrin-axis and comparator targets, as a *plausibility* input to the priority score. **This is NOT a drug recommendation, NOT evidence of efficacy in atherosclerosis, and NOT proof that any agent acts via pyrin in plaque.**

## Sources
- **ChEMBL** (`target_search`, `compound_search`) — target existence + known compounds.
- **Open Targets Platform** (`open_targets_graphql` tractability) — small-molecule / antibody tractability buckets.
- **Clinical trials** (verified Day-2): CANTOS (canakinumab, NCT01327846, n=10066), COLCOT (colchicine, NCT02551094, n=4745), LoDoCo2 (colchicine, NCT03048825, n=7264) — all COMPLETED.

## Key findings (results/tables/druggability_clinical_alignment.csv)
- **MEFV (pyrin sensor):** ChEMBL target exists; Open Targets returns 4 tractability buckets. Colchicine is FMF standard-of-care acting on the pyrin pathway — but **colchicine is non-specific and its atherosclerosis benefit (COLCOT/LoDoCo2) is NOT proven to be pyrin-mediated.**
- **Pyrin backbone (RHOA/RAC1/CDC42/PKN1/PKN2):** broadly tractable (7–10 buckets) as proteins, but these are **pleiotropic GTPases/kinases with NO pyrin-specific clinical agent** — tractability ≠ pyrin-selective druggability.
- **PSTPIP1:** no ChEMBL target, limited tractability — a protein-interaction node, not a classical drug target.
- **NLRP3 (comparator):** Open Targets returned 0 tractability buckets for this Ensembl ID, but well-known preclinical SM inhibitors exist (e.g., MCC950, ChEMBL3183703). NLRP3 is the **saturated comparator** — strong clinical momentum (colchicine trials) but not pyrin-specific.
- **Shared downstream (CASP1/GSDMD/IL1B/IL18):** most tractable and clinically advanced (canakinumab approved; VX-765, IL-18BP investigational) — but by construction **NOT pyrin-specific**; IL-1β inhibition benefit (CANTOS) does not attribute to pyrin.

## How this feeds the priority score
- `druggability_or_intervention_plausibility_score` and `clinical_alignment_score` are set from this table.
- Shared-downstream nodes score high on druggability/clinical alignment but receive **zero pyrin-specificity credit** (CV11) — high tractability does not rescue a non-specific node.
- No target is called "druggable in plaque" or "validated." Language is plausibility-only.
