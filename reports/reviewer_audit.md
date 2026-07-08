# PYRIN-PLAQUE — Reviewer-Agent Audit History

A cumulative record of every adversarial reviewer-agent gate and the corrections it produced. All gates run on ground-truth computed facts, not prose.

## Gate results
| Day | Gate | Result | Key flags → action |
|-----|------|--------|--------------------|
| 1 | Data gate | PROCEED_WITH_CAVEATS (5 PASS / 2 FLAG) | Shared-downstream circularity → excluded from discriminative scoring; integrated atlas unverified → excluded as data source |
| 2 | Claim gate | 56 PASS / 1 FLAG of 57 claims | FMF-premature-atherosclerosis claim reworded to association-only; one observational claim confidence-downgraded |
| 3 | Single-cell figure/data gate | PROCEED_WITH_CAVEATS (9 PASS / 1 FLAG) | Small Leiden clusters low donor support → flagged low-confidence |
| 3.5 | (patch) | provenance + AC/PA bridge | Provenance backfilled; AC/PA mixed signal documented honestly |
| 4 | Bulk/ODE gate | PROCEED_WITH_CAVEATS (11 PASS / 2 FLAG) | Myeloid confounding collapses bulk pyrin effect; AC/PA mixed signal — both carried forward |
| 5 | Priority/ODE/target-card gate | PROCEED_WITH_CAVEATS (13 PASS / 0 FLAG) | Priming>threshold verdict stated; nothing graded high confidence |
| 6 | Final cross-artifact audit | PROCEED (13 PASS / 0 FLAG; 0 blocker/major/minor) | GSDME + GSE120521 pairing consistency confirmed |

## Corrections made because of reviewer agents
1. **Shared downstream exclusion** — PYCARD/CASP1/GSDMD/IL1B/IL18 never used as discriminative pyrin-vs-NLRP3 features.
2. **Integrated atlas excluded** as a data source (literature precedent only).
3. **Provisional labels** — single-cell annotations labeled provisional marker-based throughout.
4. **Myeloid confounding carried forward** — bulk pyrin signal collapses under myeloid residualization (+0.72→+0.14).
5. **AC/PA mixed signal carried forward** — core MEFV rise vs non-rising composite pyrin program.
6. **Generic priming > pyrin threshold** — stated from ODE sensitivity.
7. **No high-confidence target** — priority score caps at moderate.

## External auditor corrections (independent of reviewer gates)
- **GSDME** reassigned to GENERIC_PYROPTOSIS module (all 14/14 PYRIN_MODULE genes present in bulk).
- **GSE120521 pairing** reworded: stable/unstable status + paired structure verified; exact GSM-to-column crosswalk documented as an assumption.
- Day-3 donor-robustness and file-count prose figures corrected.

## Standing conclusion
No unsupported first-person claim, no first-ever language, no causal or clinical claim survives in the deliverables. Null and confounded results are preserved, not hidden. Day 7 packaging may proceed.
