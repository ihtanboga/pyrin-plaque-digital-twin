# Day 3 — Dataset & QC Summary (GSE159677)

## Dataset
- **Accession:** GSE159677 — "Decoding the transcriptome of calcified atherosclerotic plaque at single-cell resolution" (PMID 36224302; platform GPL18573, Illumina NextSeq 500).
- **Why this dataset:** verified primary single-cell source for the MVP (Day-1 gate). The integrated atlas (DS_SC_ATLAS_INTEGRATED) was NOT used — it remains NEEDS_REVIEW.
- **Downloaded object:** `GSE159677_AGGREGATEMAPPED-tisCAR6samples_featurebcmatrixfiltered.tar.gz` (307.8 MB, sha256 8cd0b7f8…), a **filtered** 10x feature-barcode matrix aggregating all 6 samples. Local path: `data/raw/GSE159677/aggregate/outs/filtered_feature_bc_matrix/`.

## Samples (n=6; 3 patients × 2 regions)
| GSM | Patient | Region | Tissue |
|-----|---------|--------|--------|
| GSM4837523 | Patient 1 AC | atherosclerotic core | carotid |
| GSM4837524 | Patient 1 PA | proximal adjacent | carotid |
| GSM4837525 | Patient 2 AC | atherosclerotic core | carotid |
| GSM4837526 | Patient 2 PA | proximal adjacent | carotid |
| GSM4837527 | Patient 3 AC | atherosclerotic core | carotid |
| GSM4837528 | Patient 3 PA | proximal adjacent | carotid |

Regions: **AC** = atherosclerotic core, **PA** = patient-matched proximal adjacent carotid.

## Cell / gene counts & QC (after filtering)
- **Cells loaded:** 51,981 → **49,380** after QC (min_genes=200, min_cells=3, pct_mt<20, n_genes<6000).
- **Genes:** 33,538 → 23,578 after gene filter.
- **Median genes/cell:** 1,339 · **Median UMIs/cell:** 3,776 · **Median % mito:** 3.6%.
- Per-sample cell counts range 3,379–15,960 (see data/metadata/gse159677_sample_metadata.csv).

## Cell annotation
- **Published per-cell labels are NOT provided in the GEO supplementary files.** Annotation here is **provisional, marker-based**: Leiden clustering (res=0.6, 22 clusters) → majority marker-score label. Six compartments: T/NK (19,114), Macrophage/Myeloid (10,371), SMC/Fibroblast (10,290), Endothelial (6,324), B/Plasma (2,707), Mast (574). Every figure labels these as provisional.

## Gene availability
- **All 22 module genes present** in the matrix (0 absent).

## MEFV detection (reported up front, carry-forward rule 4)
- **Overall MEFV detection: 0.97% of cells** — sparse, as expected for a low-abundance transcript in 3′ droplet scRNA-seq (dropout-prone).
- **By cell type:** Macrophage/Myeloid 4.30% ≫ Endothelial 0.24% > SMC/Fibroblast 0.06% ≈ T/NK 0.05%; Mast & B/Plasma 0%.
- **By sample:** 0.49–2.48% (highest in Patient 3 AC). No single donor drives cell-type results (donor-influence control: all flags False).
- **Caveat:** low MEFV detection does NOT imply absence of pyrin biology; sparse detection limits per-cell inference and is interpreted at compartment level only.

## Known limitations
- **Small n** (3 patients, 6 samples) — no formal patient-level statistics; bootstrap CIs are within-cell and descriptive.
- **Calcified plaque / carotid endarterectomy context** — advanced, calcified disease; not a time course.
- **No MEFV genotype / FMF stratification** — cannot relate expression to pyrin genetic status.
- **Provisional labels** — not author-curated; compartment-level conclusions only.

## Headline result (descriptive, not causal)
- PYRIN_BACKBONE program score is highest in **Macrophage/Myeloid** cells (pyrin_specificity_delta = z(PYRIN_BB) − z(NLRP3_BB) = **+0.74**), while NLRP3_BACKBONE peaks in **T/NK** — the two programs localize to different compartments.
- A myeloid sub-state (leiden9) shows notably higher MEFV detection (**12.9%**) and the top pyrin-backbone score.
- **Sensitivity:** random expression-matched gene sets place the myeloid PYRIN_BACKBONE score at the **82nd percentile** (elevated but not extreme, p≈0.18) — reported honestly; the pyrin-backbone signal is suggestive, not definitive.
