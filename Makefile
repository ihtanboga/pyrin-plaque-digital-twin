.PHONY: help test smoke figures tables all clean

help:
	@echo "PYRIN-PLAQUE Digital Twin — Make targets"
	@echo "  make test     - run unit tests (pytest)"
	@echo "  make smoke    - fast smoke test (imports, config, gene modules, ODE baseline, output presence)"
	@echo "  make figures  - regenerate ODE figures (single-cell/bulk figs need raw data)"
	@echo "  make tables   - regenerate ODE scenario + sensitivity tables"
	@echo "  make all      - full reproduction path (documented; needs raw GEO downloads)"

test:
	NUMBA_CACHE_DIR=/tmp/numba_cache pytest tests/ -q

smoke:
	NUMBA_CACHE_DIR=/tmp/numba_cache python scripts/smoke_test.py

tables:
	NUMBA_CACHE_DIR=/tmp/numba_cache python scripts/05_run_ode_scenarios.py

figures: tables
	NUMBA_CACHE_DIR=/tmp/numba_cache python scripts/04_pyrin_ode_baseline.py

all:
	@echo "Full reproduction:"
	@echo " 1) download GSE159677 + GSE120521 (see DATA_LICENSE.md)"
	@echo " 2) python scripts/02_score_plaque_cell_states.py   # single-cell"
	@echo " 3) python scripts/03_bulk_plaque_validation.py     # bulk"
	@echo " 4) python scripts/04_pyrin_ode_baseline.py"
	@echo " 5) python scripts/05_run_ode_scenarios.py"
	@echo " 6) python scripts/06_priority_score.py"
	$(MAKE) test

clean:
	rm -rf __pycache__ */__pycache__ .pytest_cache
