# Figure card — Day-4 baseline ODE smoke test

- **Output file:** results/figures/day4_ode_baseline_smoke.png
- **Input:** config/ode_parameters.yaml (baseline), config/simulation_scenarios.yaml
- **Generating script:** scripts/04_pyrin_ode_baseline.py (function `main`); model in src/pyrinplaque/ode_twin.py
- **Environment:** conda env `sc` (scipy 1.18.0, python 3.13.14); solver LSODA, seed 0
- **Data source:** none — this is a mechanistic simulation, not data.

**What it shows:** Baseline trajectories of the 8 dimensionless state variables (P_inactive, P_active, ASC_complex, CASP1_active, GSDMD_N, IL1B_external, IL18_external, rupture_proxy) over t=0–50.

**Interpretation:** The inactive pyrin pool depletes as P_active rises early, propagating through ASC → caspase-1 → GSDMD_N and then to IL-1β / IL-18 externalization and a rupture proxy — a qualitatively sensible pyroptosis cascade with sequential peaks. This is a **smoke test** confirming the scaffold integrates and behaves monotonically through the cascade.

**Limitations:** DIMENSIONLESS, hypothesis-generating, NOT calibrated to patients, NOT a clinical prediction model. The NLRP3 pyroptosis ODE literature is used only as a structural template; this does not assert pyrin == NLRP3. Parameter values are illustrative. Full scenarios (pyrin_sensitized, intervention) and sensitivity analysis are deferred to Day 5.
