# Day 6 — Consistency Patch Log

Cross-artifact scan of all reports, figure cards, target cards, evidence, and audit files for the 6 mandated consistency checks. Patches applied in place; artifacts re-versioned.

## Patches applied
1. **GSDME module label** — `reports/day4_bulk_ode_gate_review.md` and `audit/day4_bulk_ode_flags.csv` reworded (v2): GSDME is a GENERIC_PYROPTOSIS module gene; **all 14/14 PYRIN_MODULE genes are present** in the bulk dataset (31/32 total module genes; only GSDME absent). Verified no remaining 'PYRIN module gene (GSDME)' phrasing anywhere.
2. **GSE120521 pairing** — canonical conservative wording applied to `reports/day4_gse120521_data_note.md`, `reports/day4_bulk_ode_gate_review.md`, and `data/metadata/gse120521_sample_metadata.csv` (new columns `status_verified`, `gsm_column_assignment`). Statement: stable/unstable status and within-patient paired structure verified from GEO characteristics + FPKM column labels; exact GSM-to-column crosswalk documented as an assumption (per-sample supplementary_files were empty). Only one version of the statement now exists in the repo.
3. **Confidence language** — scanned for 'validated target','high confidence','first ever','proven','causal','predicts clinical','drug recommendation'. All occurrences are in **banned-claims lists** (target-card `claims_to_avoid`, novelty_audit avoid-list) or **attributed external literature claims** (e.g. a cited paper's PSMA4 'causal target' finding in the evidence graph — clearly the paper's claim, not ours). No unsafe first-person claim required rewriting.
4. **NLRP3 comparator** — confirmed always 'saturated / non-novel comparator' (L2-06, L1-03); never project novelty.
5. **Shared downstream** — confirmed PYCARD/CASP1/GSDMD/IL1B/IL18 receive zero pyrin-specificity credit; GSDME is generic pyroptosis.
6. **Claude Science usage** — dedicated statement created (Part F) covering connectors, literature extraction, dataset verification, omics analysis, figures, provenance, and the six reviewer-agent gates.

## Unresolved contradictions
**None.** All six checks resolved or confirmed. The attributed-literature 'causal' language in the evidence graph (external papers' own findings) is intentionally preserved as a faithful literature record and is never used as a PYRIN-PLAQUE claim.