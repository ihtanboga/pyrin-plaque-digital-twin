# Day 7 — Smoke Test Report

Fast integrity check for the hackathon package (no raw data, no heavy compute).

## `make smoke` — scripts/smoke_test.py
| Check | Result |
|-------|--------|
| Imports (numpy/scipy/pandas/yaml + pyrinplaque modules) | PASS |
| config/gene_modules.yaml loads; PYRIN_MODULE = 14 genes | PASS |
| PYRIN_BACKBONE excludes shared downstream (9 genes) | PASS |
| ODE baseline integrates (finite, non-negative) | PASS |
| Expected output files present (priority table, readouts, figs, configs, README, licenses) | PASS |
| No raw data shipped (data/raw excluded at package time) | PASS |
| **Overall** | **PASS (exit 0)** |

## `make test` — pytest tests/
**13 passed** in ~1.3 s (module tests incl. backbone-excludes-shared-downstream; scoring tests; 6 ODE tests incl. inhibition smoke check).

## Scripts runnable without raw data
- `scripts/04_pyrin_ode_baseline.py` — ODE baseline (fast, <5 s)
- `scripts/05_run_ode_scenarios.py` — 4 scenarios + readouts (fast, <10 s)
- `scripts/06_priority_score.py` — priority score from handoff candidates (fast)
- `scripts/smoke_test.py` — integrity check (fast)

## Scripts requiring raw GEO downloads (documented in `make all`)
- `scripts/02_score_plaque_cell_states.py` — needs GSE159677 (~308 MB tar)
- `scripts/03_bulk_plaque_validation.py` — needs GSE120521 FPKM (~5.6 MB)

## Runtime
Smoke + tests complete in **under 5 seconds** combined on the reference environment.

## Missing artifacts
None. All expected derived tables, figures, configs, and docs present.

## Raw data status
Raw data present locally under `data/raw/` for development but **excluded from the package** (see final_package_manifest.md). The scored `.h5ad` (~1.9 GB) is not shipped; regenerate via `make all`.

## Verdict
**The final package is reproducible enough for hackathon review:** the fast path (`make smoke` + `make test`) verifies integrity in seconds; the full path (`make all`) is documented and requires only public downloads.
