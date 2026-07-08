# PYRIN-PLAQUE — Artifact Index

Complete machine-readable inventory with checksums: [audit/day7_artifact_inventory.csv](audit/day7_artifact_inventory.csv) (143 files, all packaged, ~10.8 MB; raw data and the ~1.9 GB scored `.h5ad` excluded).

## Core deliverables
| Path | Type | Description |
|------|------|-------------|
| README.md | doc | Demo-ready project README |
| reports/submission_summary.md | doc | 100–200 word submission summary (+100/150-word variants) |
| reports/day6_integrated_results.md | report | 13-section integrated scientific results |
| reports/reviewer_audit.md | report | Cumulative 6-gate reviewer-agent audit history |
| reports/claude_science_usage_statement.md | report | Concrete Claude Science usage |
| LICENSE / DATA_LICENSE.md / CITATION.cff | license | MIT code license; data attribution; citation metadata |

## Config (frozen)
| Path | Description |
|------|-------------|
| config/gene_modules.yaml | Frozen v0 modules (PYRIN 14 / NLRP3 10 / GENERIC 8) |
| config/scoring_weights.yaml | Priority score weights (9 positive / 9 penalty) |
| config/ode_parameters.yaml · config/simulation_scenarios.yaml | ODE baseline + 4 scenarios |

## Evidence graph
| Path | Description |
|------|-------------|
| evidence/evidence_graph.jsonl | 57 claims + source nodes/edges |
| evidence/claim_cards.md · evidence/novelty_audit.md | Claim cards; conservative novelty audit |
| manifests/evidence_table.tsv | 57-row evidence table |

## Single-cell (GSE159677)
| Path | Description |
|------|-------------|
| results/tables/pyrin_specificity_delta.csv | z(PYRIN_BACKBONE) − z(NLRP3_BACKBONE) by cell type |
| results/tables/cell_state_summary.csv · day3_cellstate_ranking.csv | Aggregated scores; top-state ranking |
| results/tables/day3_random_geneset_control.csv · day3_donor_influence.csv | Controls |
| results/figures/fig2_cellstate_program_map.png | Cell-state program heatmap (+ card) |

## Bulk (GSE120521)
| Path | Description |
|------|-------------|
| results/tables/bulk_status_tests.csv | Paired stable/unstable module tests |
| results/figures/fig3_bulk_validation.png | Bulk + myeloid confounding (+ card) |

## ODE digital twin
| Path | Description |
|------|-------------|
| results/tables/ode_scenario_readouts.csv · ode_sensitivity_ranked_parameters.csv | 4-scenario readouts; ranked sensitivity |
| results/figures/fig4_pyrin_ode_twin.png | Scenarios + tornado (+ cards) |

## Priority & targets
| Path | Description |
|------|-------------|
| results/tables/pyrin_plaque_priority_table.csv · priority_score_components.csv | Decomposed 11-candidate priority |
| results/target_cards/ (6) | Top-6 target/cell-state cards |
| results/figures/fig5_priority_ranking.png | Decomposed ranking (+ card) |
| results/tables/druggability_clinical_alignment.csv | ChEMBL/OpenTargets/trials druggability |

## Figures & provenance
| Path | Description |
|------|-------------|
| results/figures/fig1_mechanism_architecture.png | Mechanism + workflow schematic (+ card) |
| audit/day6_final_figure_inventory.csv | 5-figure inventory with checksums |
| audit/day3_provenance.json · requirements-freeze.txt · environment.yml | Provenance + pinned env |

## Code
| Path | Description |
|------|-------------|
| src/pyrinplaque/ | modules · scoring · bulk · ode_twin · priority · plotting · provenance |
| scripts/ | 02–06 analysis + smoke_test.py |
| tests/ | 13 unit tests (pass) |
| Makefile | test / smoke / figures / tables / all |
| app/streamlit_app.py | Static demo app (optional) |

## Audit
reports/day1–day6 reviews + gate reviews; audit/*_flags.csv; day7_stale_artifact_check.csv; day7_reproducibility_status.csv.

## Excluded from package (regenerate)
| Item | Reason | Regenerate |
|------|--------|-----------|
| data/raw/ | Raw public GEO data — not redistributed | Download per DATA_LICENSE.md |
| data/processed/*.h5ad | ~1.9 GB scored object | `make all` after downloads |
