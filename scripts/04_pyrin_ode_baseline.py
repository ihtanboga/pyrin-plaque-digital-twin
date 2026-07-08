#!/usr/bin/env python
"""Day-4 Part B: run baseline pyrin-adapted pyroptosis ODE and save trajectory + smoke figure.
DIMENSIONLESS, hypothesis-generating. Full scenarios + sensitivity are Day-5 work.

Usage: python scripts/04_pyrin_ode_baseline.py
"""
import sys
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "src")
from pyrinplaque import ode_twin as O

def main():
    p = O.load_params("config/ode_parameters.yaml", "baseline")
    t, traj = O.simulate(p, t_span=(0, 50), n_points=501)
    df = pd.DataFrame(traj)
    df.to_csv("data/processed/day4_ode_baseline_trajectory.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for name in O.STATE_NAMES:
        ax.plot(df["t"], df[name], label=name, lw=1.4)
    ax.set_xlabel("time (dimensionless)"); ax.set_ylabel("state value (dimensionless)")
    ax.set_title("Baseline pyrin pyroptosis ODE (dimensionless, hypothesis-generating)")
    ax.legend(ncol=2, fontsize=6)
    fig.tight_layout(); fig.savefig("results/figures/day4_ode_baseline_smoke.png", dpi=300)
    print("ODE baseline done; final rupture_proxy=%.4f" % df["rupture_proxy"].iloc[-1])

if __name__ == "__main__":
    main()
