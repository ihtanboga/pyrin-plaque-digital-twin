# Figure card — ODE sensitivity tornado

- **Output:** results/figures/day5_ode_sensitivity_tornado.png
- **Input:** config/ode_parameters.yaml (baseline ± ranges)
- **Code:** Day-5 OAT + LHS sensitivity cell; tables ode_oat_sensitivity.csv, ode_global_sensitivity.csv, ode_sensitivity_ranked_parameters.csv
- **Environment:** conda env `sc` (scipy 1.18.0); LHS N=400, seed 0
- **Data source:** none — mechanistic simulation.

**Shows:** Spearman ρ of each parameter vs time_to_50%_IL1B (global sensitivity, LHS N=400). Red speeds cytokine onset, blue delays.

**Interpretation:** `decay` dominates overall (scale parameter). Among mechanism levers, **priming_input (mean|ρ|=0.34, IL1B-timing ρ=-0.56) clearly outranks pyrin_activation_threshold (mean|ρ|=0.19, ρ=0.31)**. Per the Day-5 reviewer requirement: **the ODE supports a generic inflammatory-priming explanation more than a pyrin-specific threshold explanation.**

**Limitations:** Dimensionless, uncalibrated; sensitivity indices are qualitative parameter rankings, not effect magnitudes for any biological quantity.
