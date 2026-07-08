# Day 4 — Reviewer-Agent Bulk + ODE Gate

**Project:** PYRIN-PLAQUE Digital Twin · **Gate:** Day-4 (Day-3.5 patch + AC/PA bridge + GSE120521 bulk validation + minimal ODE scaffold).
**Reviewer:** adversarial reviewer agent (`host.llm`, reasoning model) — audits computed results only.

| # | Check | Status | Finding |
|---|-------|--------|---------|
| 1 | Did Day3.5 provenance patch fix missing/incomplete environment metadata? | **PASS** | Empty packages field patched with 15 exact versions, script/figure/table checksums, 62-pkg env lock, and requirements-freeze.txt. |
| 2 | Was AC-vs-PA bridge computed only from existing scored data? | **PASS** | Computed solely from gse159677_scored.h5ad with no re-clustering, re-annotation, or new download. |
| 3 | Are Day-3 provisional labels still described as provisional? | **PASS** | Labels consistently termed provisional marker-based across bridge report, figure card, and QC summary. |
| 4 | Is GSE120521 truly the bulk dataset used? | **PASS** | Verified public GSE120521 carotid plaque bulk RNA-seq FPKM file (checksum 27e608a5), 8 samples/4 patients. |
| 5 | Are stable/unstable labels and pairing documented? | **PASS** | Stable/unstable status and within-patient paired structure verified from GEO characteristics and FPKM column labels; exact GSM-to-column crosswalk not directly encoded in the FPKM file is documented as an assumption. FPKM-only caveat documented. |
| 6 | Are bulk scores computed from frozen gene modules? | **PASS** | Scores use config/gene_modules.yaml (no hard-coded lists); 31/32 module genes present (all 14/14 PYRIN_MODULE genes present; only GSDME from GENERIC_PYROPTOSIS absent). |
| 7 | Are shared downstream genes excluded from discriminative PYRIN-vs-NLRP3 claims? | **PASS** | PYCARD/CASP1/GSDMD/IL1B/IL18 excluded from both PYRIN and NLRP3 backbones. |
| 8 | Is myeloid composition confounding addressed? | **FLAG** | Residualizing on myeloid score collapses stable-vs-unstable pyrin diff from +0.72 to +0.14 (p=0.69), indicating the main bulk PYRIN effect is largely myeloid-abundance-driven; honestly reported but weakens headline finding. |
| 9 | Are small sample limitations stated clearly? | **PASS** | n=4 patients/8 samples stated throughout; effect sizes and bootstrap CIs prioritized over p-values, FPKM-only limitation noted. |
| 10 | Is bulk framed as triangulation, not proof? | **PASS** | All bulk figure cards/data notes explicitly state triangulation only, no causal/localization claims. |
| 11 | Is the ODE dimensionless and hypothesis-generating? | **PASS** | 8-state dimensionless ODE explicitly labeled hypothesis-generating, uncalibrated, non-clinical; 13/13 tests pass, no NaN/inf. |
| 12 | Are all figures linked to code and figure cards? | **PASS** | All six Day-4 figures have figure cards naming code, environment, and data source. |
| 13 | Should Day 5 proceed? | **FLAG** | Core pipeline and provenance are sound, but myeloid confounding of the primary bulk PYRIN signal and mixed AC-vs-PA bridge results must be carried forward as caveats into Day 5, not treated as resolved. |

## Verdict: **PROCEED_WITH_CAVEATS**  (PASS 11 · FLAG 2 · FAIL 0)

### Key flags carried to Day 5
- Myeloid-composition residualization collapses the primary bulk stable-vs-unstable PYRIN_BACKBONE difference (0.72 -> 0.14, paired p=0.69), showing the raw d=1.44 effect is substantially confounded by cell-type abundance rather than pure PYRIN pathway activity.
- AC-vs-PA single-cell bridge shows MEFV detection up in AC (3/3) but PYRIN_BACKBONE score not consistently higher (1/3) and pyrin-specificity delta actually lower in AC (0/3) because NLRP3_BACKBONE also rises in core — correctly framed as non-validating, but is a genuine mixed/contradictory signal that should not be glossed over in Day 5 synthesis.
- One GENERIC_PYROPTOSIS module gene (GSDME) absent from bulk dataset (31/32 module genes present; all 14/14 PYRIN_MODULE genes present) — minor but should be tracked for reproducibility.

### Headline honest findings
- **Bulk triangulation is directionally supportive but confounded:** PYRIN_BACKBONE is higher in unstable plaque in 4/4 patients (raw d=1.44), but after residualizing on a myeloid-marker score the stable-vs-unstable difference collapses (+0.14, p≈0.69) — much of the bulk pyrin rise reflects myeloid abundance, not per-cell pyrin biology. Stated explicitly.
- **AC-vs-PA single-cell bridge is mixed:** MEFV detection rises in plaque core (3/3 patients) but the composite PYRIN_BACKBONE and pyrin_specificity_delta do not (NLRP3 rises in core). Reported honestly, not oversold.
- **ODE is an explicitly dimensionless, uncalibrated, hypothesis-generating scaffold**; baseline only, 13/13 tests pass.