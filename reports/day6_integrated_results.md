# PYRIN-PLAQUE Digital Twin — Day 6 Integrated Results

*Auditable public-data workbench mapping MEFV/pyrin inflammasome biology to human atherosclerotic plaque cell states, bulk plaque instability, and a simple pyroptosis digital twin.*

> Language convention: this report uses **suggestive / hypothesis-generating / triangulation / caveat-aware / future validation**. It deliberately avoids *validated, causal, first-ever, clinical prediction, drug target, treatment recommendation*.

## 1. Project question
Do MEFV/pyrin inflammasome programs converge on specific human atherosclerotic plaque cell states — distinct from the well-studied generic NLRP3 program — and can a simple pyrin-adapted pyroptosis model generate testable, honestly-caveated IL-1β/IL-18 hypotheses, all through an auditable public-data + Claude Science workflow?

## 2. What was already known / novelty risk
NLRP3-inflammasome involvement in atherosclerosis is **well established** (CANTOS/canakinumab, colchicine trials). Combining single-cell + bulk + Mendelian-randomization/omics + AI agents is now **common** in cardiovascular omics (Day-2 method-precedent matrix: 6 atherosclerosis papers combine ≥3 such methods). **Therefore we make no novelty claim for the pipeline itself.** The novelty risk is overclaiming; the mitigation is conservative framing plus reviewer-agent gates.

## 3. What PYRIN-PLAQUE adds
A **conservative, reproducible contribution**: (a) an explicit **pyrin-vs-NLRP3 backbone separation** that excludes shared downstream effectors from discriminative scoring; (b) **cell-state-resolved** pyrin-permissiveness scoring in verified public plaque single-cell data; (c) a **caveat-aware priority score** that penalizes confounded and mixed signals rather than hiding them; (d) a **dimensionless digital-twin** that generates IL-1β/IL-18 hypotheses and, tellingly, shows generic priming outweighing the pyrin-specific threshold; (e) end-to-end **auditable Claude Science provenance** with six reviewer-agent gates.

## 4. Data sources actually used
- **GSE159677** (verified) — primary single-cell: 6 samples, 3 patients × atherosclerotic core (AC) + proximal-adjacent (PA); 51,981 cells (49,380 post-QC), 10x, calcified carotid plaque.
- **GSE120521** (verified) — primary bulk triangulation: 8 samples, 4 patients × stable/unstable regions, FPKM only. Stable/unstable status and paired structure verified; exact GSM-to-column crosswalk documented as an assumption.
- **DS_SC_ATLAS_INTEGRATED** — **NOT used as a data source** (NEEDS_REVIEW). The integrated atlas (Nat Commun 2025, PMID 40931012) is cited as **literature method-precedent only**.
- Gene modules frozen v0 (22 symbols, MyGene-validated): PYRIN_MODULE (14), NLRP3_COMPARISON_MODULE (10), GENERIC_PYROPTOSIS_MODULE (8); backbones exclude the 5 shared downstream genes (PYCARD, CASP1, GSDMD, IL1B, IL18).

## 5. Single-cell result (GSE159677)
**Suggestive myeloid pyrin-permissiveness.** MEFV detection is **sparse overall (~0.97%)** but **myeloid-enriched (Macrophage/Myeloid 4.3%)**. PYRIN_BACKBONE score is **highest in Macrophage/Myeloid** (pyrin_specificity_delta = z(PYRIN_BACKBONE) − z(NLRP3_BACKBONE) = **+0.74**, highest of all compartments), and the signal **survives shared-downstream exclusion** (myeloid ranks #1 in both FULL and BACKBONE). **NLRP3_BACKBONE peaks in a different (T/NK-weighted) compartment.** Caveats: the random matched gene-set control is **modest (~82nd percentile, p≈0.18) — suggestive, not definitive**; cell labels are **provisional marker-based** (no published per-cell labels); n=3 patients; no single-donor artifact (all `single_donor_flag=False`).

## 6. AC-vs-PA bridge result
**Mixed, honestly reported.** MEFV **detection rises in plaque core (AC) in 3/3 patients** (+2.97 pp, paired). But the composite **PYRIN_BACKBONE does not consistently rise** in AC (1/3), and **pyrin_specificity_delta falls** in AC (0/3) because **NLRP3_BACKBONE rises fastest in the core**. Therefore the AC-core signal is **MEFV-detection-driven and NLRP3-competed, not a clean pyrin-program gain**. This is an internal bridge, **not** independent validation.

## 7. Bulk validation result (GSE120521)
**Directionally supportive but myeloid-abundance-confounded.** PYRIN_BACKBONE is **higher in unstable plaque in 4/4 paired patients** (paired diff +0.72, Cohen's d=1.44). However, **residualizing on a myeloid-marker score collapses the difference to +0.14 (paired t p≈0.69)**, and PYRIN_BACKBONE correlates with myeloid abundance (r=0.68). Therefore **much of the bulk pyrin signal reflects increased myeloid content in unstable plaque, not per-cell pyrin biology**. Bulk is **triangulation only** (n=4, FPKM-only, no cell-type resolution).

## 8. ODE / digital-twin result
The pyrin-adapted pyroptosis ODE is **dimensionless, uncalibrated, hypothesis-generating — not a clinical prediction model**. Across 4 scenarios and global sensitivity (LHS N=400), **generic `priming_input` dominates the pyrin-specific `pyrin_activation_threshold`** on cytokine timing (mean|ρ| 0.34 vs 0.19; IL-1β-timing ρ −0.56 vs +0.31). **Therefore the ODE supports a generic inflammatory-priming explanation more than a pyrin-specific threshold explanation** — consistent with the bulk myeloid-confounding finding. It nonetheless yields concrete, falsifiable IL-1β/IL-18 kinetic hypotheses for future calibrated work.

## 9. Caveat-aware priority score
Fully decomposable (final = Σ positive·weight − Σ penalty·weight), config-driven, 11 candidates. **Top: MEFV/pyrin threshold axis — 37.9, moderate confidence only.** The Macrophage/Myeloid pyrin-permissive state is **low confidence (33.1)**. **No candidate reaches the high band (≥55).** The myeloid-confounded bulk signal (L1-05, 21.6) and the mixed AC-core signal (L1-04, 11.7) are **penalty-suppressed to the bottom** by `myeloid_confounding_penalty` and `ac_pa_mixed_signal_penalty`.

## 10. Druggability / clinical alignment
**Tractability ≠ pyrin-selective druggability.** Pyrin-backbone proteins (RhoA/Rac/Cdc42/PKN) are broadly tractable but are pleiotropic and have **no pyrin-selective clinical agent**. **Colchicine's atherosclerosis benefit (COLCOT/LoDoCo2) does not prove a pyrin mechanism in plaque.** IL-1β/CASP1/GSDMD/IL18 are the most drug-advanced (canakinumab approved) but are **shared downstream, not pyrin-specific** — they receive zero pyrin-specificity credit.

## 11. Reviewer-agent audit history
Six adversarial gates, each on ground-truth computed facts: **Day-1 data gate** (PROCEED_WITH_CAVEATS; flagged shared-downstream circularity + atlas), **Day-2 claim gate** (56 PASS / 1 FLAG of 57 claims), **Day-3 single-cell figure gate** (9 PASS / 1 FLAG), **Day-4 bulk/ODE gate** (11 PASS / 2 FLAG: myeloid confounding + AC/PA mixed), **Day-5 priority/ODE/target-card gate** (13 PASS / 0 FLAG), **Day-6 final audit** (this report). Two external auditor corrections (GSDME attribution; GSE120521 pairing wording) were applied.

## 12. Final safe interpretation
**PYRIN-PLAQUE identifies a suggestive myeloid pyrin-permissiveness pattern and prioritizes the MEFV/pyrin threshold axis for future validation, while showing that much of the bulk signal is myeloid-abundance-confounded and that generic inflammatory priming may dominate over pyrin-specific threshold effects.** It is an **auditable hypothesis engine**, not a validated-target or clinical claim.

## 13. Claims to avoid
MEFV causes atherosclerosis; validated pyrin target; colchicine proves pyrin mechanism; single-cell expression proves causality; ODE predicts clinical outcomes; "first ever"; drug recommendation; patient-stratification tool.
