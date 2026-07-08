# Novelty audit — PYRIN-PLAQUE Digital Twin (Day 2)

## Purpose
Establish, from real literature, what is ALREADY common so the project's novelty claim stays conservative and defensible. No "first ever" language is used anywhere.

## What the literature shows is already saturated
1. **NLRP3 inflammasome in atherosclerosis is extensively studied.** Axis-B claims (n=11) document NLRP3/macrophage/endothelial pyroptosis and plaque instability as an established field. Anchor A2 (macrophage galectin-3 -> TLR4/NLRP3 pyroptosis) and A3 (HDAC3 -> endothelial NLRP3) are recent 2025 examples.
2. **The "single-cell + bulk + MR/TWAS + druggability" pipeline is now common in cardiovascular omics.** 6 atherosclerosis papers in our axis-E sample already combine 3+ of: scRNA-seq, bulk RNA-seq, Mendelian randomization/SMR, TWAS, druggability, molecular docking. This means the multi-omics + target-ranking combination is NOT itself novel.
3. **An integrated single-cell atlas of human plaques already exists** (anchor A4, Nat Commun 2025, PMID 40931012). We therefore do NOT claim to build an atlas, and we keep this as a literature precedent while its use as a DATA source remains NEEDS_REVIEW (Day-1 gate).

## What remains genuinely under-addressed (defensible contribution)
- **Pyrin/MEFV-SPECIFIC mechanism separation from NLRP3** at the level of a curated, version-pinned gene module with a pyrin-specific regulatory backbone (MEFV, RhoA/RAC1/CDC42, PKN1/2, 14-3-3, PSTPIP1) distinct from NLRP3 (axis-A + InterPro domain evidence C-A-* ). Most atherosclerosis inflammasome work treats IL-1β/GSDMD signal as generically "NLRP3".
- **Cell-state-resolved pyrin-permissiveness** in human plaque, with shared downstream effectors explicitly excluded from discriminative scoring (Day-1 caveat 1).
- **An auditable, provenance-linked workflow** where every claim/number/figure is source-linked and adversarially reviewed — a process contribution, not a biological firstness claim.

## Safest novelty statement (defensible)
"We did not identify a public, auditable workflow that (a) separates MEFV/pyrin-specific inflammasome biology from generic NLRP3 at the gene-module level, (b) maps it to human plaque cell states with shared-effector controls, and (c) couples it to a pyrin-adapted pyroptosis digital twin under claim-level provenance. Individual components — single-cell plaque atlases, NLRP3 pyroptosis studies, multi-omics target ranking — are already common; this project's contribution is the pyrin-specific, cell-state-resolved, auditable integration."

## Claims to AVOID (hard bans, enforced by reviewer gate)
- "First ever" anything.
- "MEFV causes atherosclerosis."
- "Single-cell expression proves causality."
- "Colchicine benefit proves a pyrin mechanism."
- "The ODE model predicts clinical outcomes."
- Any conflation of pyrin with NLRP3.

## Top method-precedent papers that threaten novelty (must be cited to stay honest)
- PMID 41806154 (J Cardiovasc Transl Res 2026, 4 methods): Bulk RNA-seq + WGCNA + machine learning + SMR + scRNA-seq + molecular docking drug prediction
- PMID 40555237 (Am J Hum Genet 2025, 3 methods): scRNA-seq + TWAS + cis-MR + colocalization for atherosclerotic disease gene discovery
- PMID 42032825 (Ann Hum Genet 2026, 3 methods): SMR + two-sample MR + multiomics (mQTL, sQTL) + scRNA-seq + drug prediction (DGIdb)
- PMID 42176975 (Ann Vasc Surg 2026, 3 methods): Bulk RNA-seq + SMR + scRNA-seq + spatial transcriptomics + MR-scRNA immune phenotyping
- PMID 41601048 (Ecotoxicol Environ Saf 2026, 3 methods): MR + network toxicology (PPI network, enrichment) + scRNA-seq + molecular docking
- PMID 37351287 (Front Cardiovasc Med 2023, 3 methods): TWAS (4 approaches) + SMR + scRNA-seq immune characterization in coronary atherosclerosis
- PMID 39364414 (Front Immunol 2024, 2 methods): scRNA-seq + bulk RNA-seq + WGCNA + machine learning risk model + ROC validation
- PMID 41262249 (Front Immunol 2025, 2 methods): Bulk RNA-seq + WGCNA + machine learning (SVM-RFE, LASSO, XGBoost, RF) + scRNA-seq analysis
