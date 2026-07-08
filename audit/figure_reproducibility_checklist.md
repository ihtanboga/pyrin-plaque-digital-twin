# Day 6 — Figure Reproducibility Checklist

Every main figure traces to a named script, a pinned environment, an input, a checksum, and a figure card.

| Figure | Script | Env | Input | Checksum | Figure card |
|--------|--------|-----|-------|----------|-------------|
| Figure 1 | Day-6 Fig-1 cell | sc | schematic (frozen modules + workflow) | `e264b4833b03cbf6` | reports/figure_cards/fig1_mechanism_architecture.md |
| Figure 2 | scripts/02_score_plaque_cell_states.py | sc | gse159677_scored.h5ad | `0d6adf028ddc8ec6` | reports/figure_cards/day3_fig2_cellstate_program_map.md |
| Figure 3 | scripts/03_bulk_plaque_validation.py | sc | GSE120521 FPKM | `c765784ccb5817c3` | reports/figure_cards/day4_fig3_bulk_validation.md |
| Figure 4 | scripts/05_run_ode_scenarios.py | sc | ODE simulation | `e70573ee27568702` | reports/figure_cards/day5_fig4_pyrin_ode_twin.md |
| Figure 5 | scripts/06_priority_score.py | sc | priority tables | `42f2ae14c0b3c43c` | reports/figure_cards/day5_fig5_priority_ranking.md |

**Environment lock:** requirements-freeze.txt (14 pinned deps); audit/day3_environment_lock.txt (62 pkgs).
**Numba note:** set NUMBA_CACHE_DIR=/tmp/numba_cache before importing scanpy (documented in audit/day3_provenance.json).
**Frozen inputs:** config/gene_modules.yaml (sha256 f1078086…); GSE159677 tar (8cd0b7f8…); GSE120521 FPKM (27e608a5…).

All five main figures reproducible from clean config + named scripts. 13/13 unit tests pass.