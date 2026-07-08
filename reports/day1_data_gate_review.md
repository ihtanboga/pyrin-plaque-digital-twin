# Day 1 — Reviewer-Agent Data Gate Review

**Project:** PYRIN-PLAQUE Digital Twin
**Gate:** Day 1 Data Gate (scope lock, source registry, dataset verification, gene-module freeze)
**Date:** 2026-07-07
**Reviewer:** adversarial reviewer agent (Claude Science `host.llm`, reasoning model) — audits only; creates no new biological claims.

## Scope of this gate
This review audits ONLY Day-1 deliverables. No literature-claim extraction, no single-cell analysis, no ODE work was performed today (per instruction to run Day 1 only).

## Inputs audited (ground truth, verified this session)
- Dataset verification via `omics-archives` `geo_get_series` (real GEO metadata).
- Gene-symbol validation via `genes-ontologies` `query_genes` (mygene.info).
- Files: `manifests/data_manifest.tsv`, `config/datasets.yaml`, `config/gene_modules.yaml`, `results/tables/gene_module_overlap.csv`, `data/metadata/source_registry.csv`.

## Seven-point data gate

DAY-1 DATA GATE REVIEW — PYRIN-PLAQUE DIGITAL TWIN

1. Are all datasets real and public?
PASS — All three core datasets (GSE120521, GSE159677, GSE253902) carry status VERIFIED with explicit public-release dates and associated PubMed IDs (31339449, 36224302, 38385291); no dataset is presented as public without a corroborating verification stamp.

2. Are accessions correct/consistent with stated assay & tissue?
PASS — Each accession's assay type and tissue align internally (bulk RNA-seq / FPKM matches GSE120521's carotid stable-vs-unstable design; scRNA-seq / 10x matrices match GSE159677's carotid core+adjacent; CITE-seq / mtx+ADT matches GSE253902's carotid artery), with no cross-assay or cross-tissue mismatch in the stated facts.

3. Are any datasets restricted or license-problematic?
PASS — No restricted-access, embargo, or controlled-access flags appear on any of the three verified GEO entries; all show standard "Public on [date]" status typical of open GEO records.

4. Are all gene symbols valid?
PASS — 22/22 symbols validated, 0 not_found, all confirmed taxid 9606 via mygene.info; no orphan or species-ambiguous symbols in the panel.

5. Is pyrin kept separate from NLRP3?
PASS — A distinct PYRIN-specific backbone (MEFV, RHOA, RAC1, CDC42, PKN1, PKN2, YWHAB, YWHAZ, PSTPIP1) is enumerated separately from the NLRP3-specific set (NLRP3, NEK7, TXNIP, NFKB1, RELA) and from the generic-inflammasome set (CASP4, CASP5, GSDME); no gene appears in both the pyrin and NLRP3 lists.

6. Are any genes circular/cherry-picked without rationale?
FLAG — The shared downstream set (PYCARD, CASP1, GSDMD, IL1B, IL18) is common to both pyrin and NLRP3 arms, which is biologically expected (canonical shared inflammasome effector machinery), not evidence of circularity — but because these genes sit downstream of both backbones, any classifier or signature built using them as discriminating features between pyrin- and NLRP3-driven states would be circular. Flagging so Day-2 modeling explicitly excludes shared-downstream genes from any "pyrin vs NLRP3" discriminative feature set and uses them only as pathway-activity/output nodes.

7. Are there unresolved blockers before Day 2?
FLAG — DS_SC_ATLAS_INTEGRATED is NEEDS_REVIEW with no verified single-accession public download object; it is not a hard blocker for Day-2 work since three independently verified GEO datasets already cover bulk RNA-seq, scRNA-seq, and CITE-seq, but any Day-2 claim relying on an "integrated atlas" must not proceed until a concrete, verifiable accession/object is located — currently that claim class is unsupported.

Overall recommendation: PROCEED_WITH_CAVEATS — The three core datasets and 22-gene module pass identity, accession-consistency, licensing, and pyrin/NLRP3-separation checks cleanly, so foundational Day-1 work (data loading, QC, module scoring) can proceed; however, Day-2 plans must explicitly (a) exclude the shared downstream effector genes (PYCARD, CASP1, GSDMD, IL1B, IL18) from any pyrin-vs-NLRP3 discriminative/classification feature set to avoid circularity, and (b) refrain from any analysis or claim that depends on the unverified "integrated atlas" (DS_SC_ATLAS_INTEGRATED) until a concrete public accession is confirmed.

## Actions carried into Day 2 (from FLAGs)
1. **Specificity control is mandatory (FLAG #6):** In all pyrin-vs-NLRP3 discriminative scoring, compute a "pyrin backbone" score EXCLUDING shared downstream genes (PYCARD, CASP1, GSDMD, IL1B, IL18). Shared genes are used only as pathway-output nodes, never as discriminating features. This is already specified as a control in `05_analysis_methods.md`; the gate makes it a hard requirement.
2. **Integrated atlas gated (FLAG #7):** `DS_SC_ATLAS_INTEGRATED` remains `NEEDS_REVIEW`. No Day-2+ claim may depend on an "integrated atlas" until a concrete public accession/download object is verified. Verified primary single-cell source for the MVP is **GSE159677**; CITE-seq **GSE253902** is the stretch independent validation.

## Verdict
**PROCEED_WITH_CAVEATS** — Day-1 foundations (three verified public datasets spanning bulk / scRNA-seq / CITE-seq; 22/22 validated gene symbols; pyrin backbone cleanly separated from NLRP3) pass the gate. The two FLAGs are process constraints for Day 2, not blockers.
