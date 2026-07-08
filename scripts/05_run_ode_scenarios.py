#!/usr/bin/env python
"""Day-5 Part B/C: run 4 pyrin ODE scenarios + OAT/LHS sensitivity.
Dimensionless, hypothesis-generating. NOT calibrated; NOT clinical prediction.

Usage: python scripts/05_run_ode_scenarios.py
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from pyrinplaque import ode_twin as O

def main():
    base = O.load_params("config/ode_parameters.yaml", "baseline")
    cfg = O.load_scenarios("config/simulation_scenarios.yaml")
    scen, tspan, npts = cfg["scenarios"], tuple(cfg["t_span"]), cfg["n_points"]
    order = ["baseline", "pyrin_sensitized_low_threshold",
             "high_inflammatory_priming", "inhibited_or_delayed_activation"]
    trajs, reads = [], []
    for name in order:
        p = dict(base); p.update(scen[name].get("overrides") or {})
        t, traj = O.simulate(p, t_span=tspan, n_points=npts)
        df = pd.DataFrame(traj); df["scenario"] = name; trajs.append(df)
        r = {"scenario": name}; r.update(O.scenario_readouts(t, traj)); reads.append(r)
    pd.concat(trajs, ignore_index=True).to_csv("data/processed/ode_scenario_trajectories.csv", index=False)
    pd.DataFrame(reads).set_index("scenario").to_csv("results/tables/ode_scenario_readouts.csv")
    print("scenarios done")

if __name__ == "__main__":
    main()
