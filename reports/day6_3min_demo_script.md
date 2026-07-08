# Day 6 — 3-Minute Demo Script

**[0:00–0:25 · Problem]**
Inflammation drives atherosclerosis, and the NLRP3 inflammasome is already a saturated, well-studied story — CANTOS, colchicine trials. But pyrin, encoded by MEFV, is a *mechanistically distinct* inflammasome sensor with its own regulatory logic, and it has not been well mapped in human plaque. PYRIN-PLAQUE asks: does a pyrin program localize to specific plaque cell states, distinct from generic NLRP3?

**[0:25–0:55 · Claude Science workflow]**
Everything here is public data, run through Claude Science. MCP connectors pulled GEO datasets, PubMed literature, ChEMBL and Open Targets druggability, and clinical trials. We built a literature evidence graph of 57 reviewer-gated claims, then analyzed two verified public datasets. **Claude Science reviewer agents flagged and corrected overclaims at every stage** — six adversarial gates in total.

**[0:55–1:35 · Single-cell finding]**
In the verified single-cell plaque dataset GSE159677 — nearly 50,000 cells — we scored a pyrin backbone against an NLRP3 backbone, deliberately excluding shared downstream effectors so the comparison is real. The result is a *suggestive myeloid pyrin-permissiveness*: MEFV itself is sparse, under one percent overall, but it's enriched in macrophage/myeloid cells at over four percent, and the pyrin backbone score is highest there. The NLRP3 comparator peaks in a different, T/NK-weighted compartment. Honestly: a matched random gene-set control was only modest, and the cell labels are provisional. This is suggestive, not definitive.

**[1:35–2:00 · Bulk + caveat]**
We triangulated with bulk plaque RNA-seq, GSE120521 — four patients, stable versus unstable regions, paired. The pyrin backbone is higher in unstable plaque in all four patients. But when we residualize on a myeloid-marker score, that difference collapses — from plus-point-seven-two to plus-point-one-four, no longer significant. So the bulk signal is directionally supportive but largely *myeloid-abundance-confounded*. **Null and confounded results are preserved, not hidden.**

**[2:00–2:25 · ODE]**
The digital twin is a dimensionless, hypothesis-generating pyroptosis ODE — not a clinical prediction. Its sensitivity analysis is revealing: generic inflammatory priming dominates the pyrin-specific activation threshold. So the model, honestly read, supports a generic-priming explanation more than a pyrin-specific one — the same story the bulk confounding told.

**[2:25–2:45 · Priority score]**
All of this feeds a caveat-aware priority score — fully decomposed, every penalty visible. The MEFV/pyrin threshold axis comes out on top, but only at *moderate* confidence. Nothing reaches high confidence. The myeloid-confounded bulk signal and the mixed plaque-core signal are pushed to the bottom by explicit penalties.

**[2:45–3:00 · Close]**
So what is PYRIN-PLAQUE? A suggestive myeloid pyrin-permissiveness pattern, a prioritized hypothesis for future validation, and a fully auditable trail. **Every figure has code, environment, and provenance.** This is an auditable hypothesis engine — not a validated target, and not a clinical claim.
