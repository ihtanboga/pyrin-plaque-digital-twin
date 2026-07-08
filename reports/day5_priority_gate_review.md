# Day 5 — Reviewer-Agent Priority + ODE + Target-Card Gate

**Project:** PYRIN-PLAQUE Digital Twin · **Gate:** Day-5 (caveat ledger + 4-scenario ODE + sensitivity + caveat-aware priority score + 6 target cards + druggability).
**Reviewer:** adversarial reviewer agent (`host.llm`, reasoning model) — audits computed results only.

| # | Check | Status | Finding |
|---|-------|--------|---------|
| 1 | Are all Day1-4 caveats explicitly carried into a caveat ledger (not hidden)? | **PASS** | 12 caveats (CV01-CV12) logged, all 10 mandatory carry-forwards present, each mapped to a named penalty and echoed in every target card |
| 2 | Were all 4 ODE scenarios run with readouts? | **PASS** | baseline, pyrin_sensitized_low_threshold, high_inflammatory_priming, inhibited_or_delayed_activation all run with readouts.csv and fig4 twin plot |
| 3 | Was ODE sensitivity (OAT + global) run and ranked? | **PASS** | OAT (+/-20%) plus LHS global (N=400, Spearman) run and ranked: decay > priming_input (0.34) > threshold (0.19), tornado figure produced |
| 4 | Is the priming-vs-threshold verdict stated per the reviewer rule? | **PASS** | Verdict explicitly states priming_input dominates threshold (rho -0.56 vs 0.31) and concludes ODE favors generic priming over a pyrin-specific explanation, recorded on fig4/tornado cards |
| 5 | Is the priority score fully decomposable (positive - penalties) and config-driven? | **PASS** | final = sum(positive*w) - sum(penalty*w) in priority.py, weights externalized in scoring_weights.yaml, 198-row component table backs every score |
| 6 | Do penalties actually discount the confounded/mixed signals (bulk myeloid, AC/PA)? | **PASS** | myeloid_confounding_penalty (w=9) and ac_pa_mixed_signal_penalty (w=8) applied at full strength, dropping L1-05 bulk to 21.6 and L1-04 AC-core to 11.7 (lowest of all candidates) |
| 7 | Does nothing get over-graded (no false 'high' confidence)? | **PASS** | Top score is L2-01 at 37.9 labeled 'moderate'; nothing reaches the >=55 high band, confirming conservative-by-design scoring |
| 8 | Do target cards include BOTH supporting evidence AND penalties/counter-evidence? | **PASS** | All 6 target cards carry strongest_supporting_evidence, strongest_counterevidence, and a penalty component table alongside safe_claim/claims_to_avoid |
| 9 | Is NLRP3 labeled a saturated non-novel comparator? | **PASS** | Both L2-06 (NLRP3) and L1-03 (T/NK NLRP3 compartment) explicitly labeled 'saturated comparator / not novel / not pyrin-specific' |
| 10 | Do shared-downstream genes get zero pyrin-specificity credit? | **PASS** | L2-05 (PYCARD/CASP1/GSDMD/IL1B/IL18) has pyrin_specificity_score=0 with shared_downstream_non_specificity_penalty applied at full strength |
| 11 | Is druggability framed as plausibility only (no efficacy/drug claims)? | **PASS** | druggability_clinical_alignment.csv note explicitly separates tractability from pyrin-selective efficacy; states MEFV/colchicine link is not proven pyrin-mediated in plaque |
| 12 | Were frozen gene modules left unchanged? | **PASS** | gene_modules.yaml confirmed unmodified during Day 5 |
| 13 | Should Day 6 proceed? | **PASS** | All process gates satisfied with caveats fully retained and visible; safe to proceed to Day 6 provided the ledger, penalty structure, and banned-claims list carry forward unchanged |

## Verdict: **PROCEED_WITH_CAVEATS**  (PASS 13 · FLAG 0 · FAIL 0)

### Key flags carried to Day 6
- No candidate clears the 'high' confidence band (top = 37.9, moderate) — Day 6 must not upgrade this language
- ODE sensitivity supports a generic inflammatory-priming explanation over a pyrin-specific threshold mechanism — this framing must persist into Day 6 narrative
- Lowest-ranked candidates (bulk myeloid, AC/PA-core) are penalty-suppressed for confounding, not disproven — Day 6 should not silently drop them without noting why
- MEFV/colchicine clinical trial data (CANTOS/COLCOT/LoDoCo2) remain non-pyrin-specific; any Day 6 clinical framing must retain the explicit non-causal caveat

### Headline honest findings
- **ODE supports generic inflammatory-priming over pyrin-specific threshold.** Global sensitivity: priming_input (mean|ρ|=0.35) outranks pyrin_activation_threshold (mean|ρ|=0.19) on cytokine timing — stated explicitly per the reviewer rule. This aligns with the Day-4 myeloid confounding finding.
- **Priority score is conservative by design:** the top candidate (MEFV/pyrin threshold axis, 37.9) is only 'moderate'; nothing reaches 'high'. Penalties push the myeloid-confounded bulk signal and the mixed AC/PA core signal to the bottom.
- **Every target card shows both support and penalties**; NLRP3 is a labeled saturated comparator; shared-downstream nodes get zero pyrin-specificity credit; druggability is plausibility-only.