# Day 3.5 — Single-cell → bulk bridge (AC core vs PA adjacent)

**Source:** existing scored object `data/processed/gse159677_scored.h5ad` (no re-clustering, no re-annotation, no new download). Aggregated by patient × region × provisional cell type.

**Design:** 3 patients, each with a matched atherosclerotic core (AC) and proximal adjacent (PA) sample → paired, within-patient. n=3 ⇒ **descriptive paired statistics only; no p-values, no significance claims.**

## Main question
Within the same patients, is Macrophage/Myeloid **PYRIN_BACKBONE** or **pyrin_specificity_delta** directionally higher in plaque core (AC) than adjacent (PA)?

## Result (Macrophage/Myeloid, AC minus PA, per patient)
| Metric | mean AC | mean PA | paired diff (per-patient) | AC-higher | direction-consistent |
|--------|---------|---------|---------------------------|-----------|----------------------|
| PYRIN_BACKBONE_mean | 0.532 | 0.547 | -0.015 (-0.047;-0.014;+0.016) | 1/3 | False |
| NLRP3_BACKBONE_mean | 0.132 | 0.024 | +0.108 (+0.082;+0.064;+0.180) | 3/3 | True |
| pyrin_specificity_delta | 0.306 | 0.705 | -0.398 (-0.435;-0.254;-0.505) | 0/3 | True |
| MEFV detection (%) | 7.65 | 4.68 | +2.97 (+2.126;+3.011;+3.771) | 3/3 | True |

## Honest interpretation
- **MEFV detection is directionally higher in plaque core (AC) in all 3 patients** (+2.97 pp mean) — a consistent, encouraging within-patient signal that the myeloid pyrin *marker* concentrates where disease is most advanced.
- **The composite PYRIN_BACKBONE score is NOT higher in AC** (mixed, 1/3 patients) — the multi-gene backbone score is dominated by broadly-expressed regulators (RHOA/RAC1/CDC42/YWHA*) whose expression does not track the core/adjacent axis.
- **pyrin_specificity_delta is consistently LOWER in AC** (0/3 AC-higher) — because **NLRP3_BACKBONE rises in the core in all 3 patients** (+0.108). In advanced plaque core, the NLRP3 program increases more than the pyrin backbone, so the *relative* pyrin specificity drops even as MEFV itself goes up.

## Bottom line
This is an **internal bridge, not validation**. The Day-3 myeloid pyrin *marker* (MEFV) is directionally stronger in plaque core, but the composite pyrin-backbone program and the pyrin-vs-NLRP3 specificity delta are **not** — the core is where NLRP3 rises fastest. Any downstream claim must state that the core signal is MEFV-detection-driven and NLRP3-competed, not a clean pyrin-program gain. n=3, provisional labels, sparse MEFV — descriptive only.
