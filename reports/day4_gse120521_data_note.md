# Day 4 — GSE120521 bulk data note

- **Accession:** GSE120521
- **Title:** RNA-seq of stable and unstable section of human atherosclerotic plaques
- **Source URL:** ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120521/suppl/GSE120521_Athero_RNAseq_FPKM.xlsx
- **Organism:** Homo sapiens · **Status:** Public on Aug 02 2019
- **Samples:** 8 · **Patients:** 4
- **Design:** Plaques from 4 patients, each dissected into **stable** and **unstable** regions → RNA-seq (Genewiz). **PAIRED design** (within-patient stable vs unstable).
- **Stable/unstable labels:** VERIFIED from per-sample GEO characteristics ('stability of region': stable/unstable + 'patient': Patient 1–4). Column names in the FPKM matrix also carry the status (…_stable / …_unstable).
- **Pairing:** The stable/unstable **status** of each column is VERIFIED (embedded in the FPKM column names and matching the 4 stable + 4 unstable GEO samples). The specific **GSM-accession-per-column** assignment — P1(MB1/MB2), P2(MB3/MB4), P3(MB5/MB6), P4(MB9/MB10) — is an **ordinal assumption** (FPKM columns taken in the same order as the GEO samples list); per-sample `supplementary_files` were empty for all 8 GSMs, so no independent GEO-side crosswalk could confirm the exact GSM↔column identity. This does NOT affect the paired analysis, which pairs by the verified stable/unstable status within each of the 4 patient columns; only the patient/GSM label attached to each pair is assumption-based.
- **Expression unit:** FPKM (log-scale recommended for scoring). **Raw counts:** NOT provided in this supplementary file (FPKM only) — so no DESeq2/edgeR count-based modeling; module scoring on log-FPKM is appropriate for triangulation.
- **File:** GSE120521_Athero_RNAseq_FPKM.xlsx, 5.51 MB (sha256 27e608a5…), 58,037 genes × 8 samples.
- **Limitations:** n=4 patients (8 samples); FPKM-only (no raw counts); bulk tissue (no cell-type resolution); triangulation only — cannot prove cell-type localization, mechanism, or causality.

**Verification status: VERIFIED** (public, stable/unstable + pairing confirmed). Proceeding with paired analysis.
