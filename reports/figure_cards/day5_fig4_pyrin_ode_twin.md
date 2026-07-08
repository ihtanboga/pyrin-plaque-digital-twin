# Figure card — Fig 4 Pyrin ODE digital twin (4 scenarios)

- **Output:** results/figures/fig4_pyrin_ode_twin.png
- **Input:** config/ode_parameters.yaml, config/simulation_scenarios.yaml
- **Code:** scripts/05_run_ode_scenarios.py; model src/pyrinplaque/ode_twin.py; readouts results/tables/ode_scenario_readouts.csv
- **Environment:** conda env `sc` (scipy 1.18.0, python 3.13.14); LSODA solver; seed 0
- **Data source:** none — mechanistic simulation.

**Shows:** P_active, IL1B_external, IL18_external, rupture_proxy across 4 scenarios: baseline, pyrin_sensitized (low threshold), high_inflammatory_priming, inhibited/delayed.

**Interpretation:** Lowering the pyrin activation threshold (pyrin_sensitized) speeds and slightly raises the cascade, but **raising priming_input shifts the cascade MORE** (peak P_active 0.29→0.48 vs 0.35; IL1B t50 7.5→6.5 vs 7.1). The distinction matters: **the generic inflammatory-priming lever outweighs the pyrin-specific threshold lever** — consistent with Day-4 bulk myeloid confounding. Inhibited/delayed scenario lowers and delays all readouts.

**Limitations:** Dimensionless, qualitative, hypothesis-generating; NOT calibrated to any MEFV variant, drug, or patient; NO clinical prediction or treatment-efficacy claim.
