#!/usr/bin/env python
"""Revision: two-signal pyrin model (v2) — scenarios, matched perturbations, one-at-a-time
and Latin-hypercube global sensitivity (PRCC). Also documents the v1 structural symmetry.

Dimensionless and hypothesis-generating; NOT calibrated to patients, drugs or MEFV variants.

Usage (repo root): python -I scripts/05b_run_ode_v2.py
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import qmc

sys.path.insert(0, "src")
from pyrinplaque import ode_twin_v2 as O2
from pyrinplaque import ode_twin as O1

TAB = "results/revision/tables"
PROC = "data/processed"
SEED = 0
BASE = O2.load_params("config/ode_parameters_v2.yaml")

SCENARIOS = {
    "baseline": {},
    "pyrin_sensitized_low_threshold": {"pyrin_activation_threshold": 0.30},
    "high_inflammatory_priming": {"priming_input": 1.8},
    "inhibited_or_delayed_activation": {"inhibition_factor": 0.5, "pyrin_activation_threshold": 0.70},
}
# equal relative (20%) change of each lever
MATCHED = {
    "baseline": {},
    "threshold_minus20pct": {"pyrin_activation_threshold": 0.40},
    "signal2_plus20pct": {"signal2_input": 1.20},
    "priming_plus20pct": {"priming_input": 1.20},
    "inhibition_0p2": {"inhibition_factor": 0.20},
}
KEY_OUT = ["peak_P_active", "peak_IL1B_out", "IL1B_AUC", "time_to_50pct_peak_IL1B",
           "peak_IL18_out", "IL18_AUC", "time_to_50pct_peak_IL18", "IL1B_to_IL18_AUC_ratio",
           "lysed_fraction_end", "time_to_50pct_final_lysis", "IL1B_released_via_pores_fraction"]
VARY = ["priming_input", "signal2_input", "pyrin_activation_threshold", "inhibition_factor",
        "pyrin_activation_rate", "pyrin_basal_synthesis", "pyrin_priming_synthesis", "pyrin_turnover",
        "ASC_recruitment_rate", "CASP1_activation_rate", "GSDMD_cleavage_rate",
        "cytokine_maturation_rate", "proIL1B_priming_synthesis", "proIL18_constitutive_synthesis",
        "pro_turnover", "mature_turnover", "pore_release_rate", "pore_repair", "lysis_rate",
        "lysis_pore_threshold", "decay", "clearance", "saturation"]
PRCC_OUT = ["IL1B_AUC", "IL18_AUC", "time_to_50pct_peak_IL1B", "IL1B_to_IL18_AUC_ratio",
            "lysed_fraction_end", "time_to_50pct_final_lysis", "IL1B_released_via_pores_fraction"]


def run(over):
    p = dict(BASE)
    p.update(over)
    t, tr = O2.simulate(p)
    return t, tr, O2.readouts(t, tr)


def prcc(X, y):
    """Partial rank correlation coefficients (Marino et al. 2008) with t-test p values."""
    ok = np.isfinite(y)
    X, y = X[ok], y[ok]
    R = np.column_stack([stats.rankdata(c) for c in X.T])
    ry = stats.rankdata(y)
    n, k = R.shape
    out = []
    for i in range(k):
        Z = np.column_stack([np.ones(n), np.delete(R, i, axis=1)])
        bx, *_ = np.linalg.lstsq(Z, R[:, i], rcond=None)
        by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
        rx, rr = R[:, i] - Z @ bx, ry - Z @ by
        r = float(np.corrcoef(rx, rr)[0, 1])
        df = n - 2 - (k - 1)
        tval = r * np.sqrt(df / max(1e-12, 1 - r * r))
        out.append((r, float(2 * stats.t.sf(abs(tval), df)), int(n)))
    return out


def main():
    # 1) scenarios + trajectories
    rows, trajs = [], []
    for name, over in SCENARIOS.items():
        t, tr, r = run(over)
        rows.append({"scenario": name, **r})
        df = pd.DataFrame(tr)
        df["scenario"] = name
        trajs.append(df)
    pd.DataFrame(rows).set_index("scenario").to_csv(f"{TAB}/R30_ode2_scenario_readouts.csv")
    pd.concat(trajs).to_csv(f"{PROC}/ode2_scenario_trajectories.csv", index=False)

    # 2) matched (equal-relative) perturbations, expressed as % change from baseline
    base_r = run({})[2]
    mrows = []
    for name, over in MATCHED.items():
        r = run(over)[2]
        row = {"perturbation": name, **{k: r[k] for k in KEY_OUT}}
        for k in KEY_OUT:
            b = base_r[k]
            row[f"{k}_pct_change"] = 100 * (r[k] - b) / b if b else np.nan
        mrows.append(row)
    pd.DataFrame(mrows).to_csv(f"{TAB}/R31_ode2_matched_perturbations.csv", index=False)

    # 3) one-at-a-time +-20%
    orows = []
    for par in VARY:
        for lab, f in (("minus20", 0.8), ("plus20", 1.2)):
            val = BASE[par] * f if par != "inhibition_factor" else (0.0 if f < 1 else 0.2)
            r = run({par: val})[2]
            orows.append({"parameter": par, "perturb": lab, "value": val,
                          **{f"{k}_delta": r[k] - base_r[k] for k in KEY_OUT}})
    pd.DataFrame(orows).to_csv(f"{TAB}/R32_ode2_oat_sensitivity.csv", index=False)

    # 4) global: LHS N=1000, log-uniform 0.5x-2x for all positive parameters (equal relative
    #    uncertainty); inhibition_factor uniform 0-0.5
    N = 1000
    U = qmc.LatinHypercube(d=len(VARY), seed=SEED).random(N)
    X = np.empty_like(U)
    for j, par in enumerate(VARY):
        X[:, j] = 0.5 * U[:, j] if par == "inhibition_factor" else BASE[par] * 2.0 ** (2 * U[:, j] - 1)
    res = []
    for i in range(N):
        _, _, r = run(dict(zip(VARY, X[i])))
        res.append(r)
    R = pd.DataFrame(res)
    samples = pd.concat([pd.DataFrame(X, columns=VARY), R], axis=1)
    samples["gate_open"] = samples["signal2_input"] > samples["pyrin_activation_threshold"]
    samples.to_csv(f"{PROC}/ode2_lhs_samples.csv", index=False)
    grows = []
    for out in PRCC_OUT:
        y = R[out].to_numpy(float)
        if out.startswith("time_to") or out.endswith("ratio") or out.endswith("fraction"):
            y = np.where(samples["gate_open"], y, np.nan)  # timing only defined when the gate opens
        for par, (r, p, n) in zip(VARY, prcc(X, y)):
            rho = stats.spearmanr(X[np.isfinite(y), VARY.index(par)], y[np.isfinite(y)])[0]
            grows.append({"output": out, "parameter": par, "PRCC": r, "PRCC_p": p,
                          "spearman_rho": rho, "n_sets": n})
    G = pd.DataFrame(grows)
    G.to_csv(f"{TAB}/R33_ode2_global_prcc.csv", index=False)

    # 5) v1 symmetry demonstration (superseded model)
    b1 = O1.load_params("config/ode_parameters.yaml", "baseline")
    srows = []
    for lab, over in (("v1 baseline (drive 0.5)", {}),
                      ("v1 threshold 0.5->0.3 (drive 0.7)", {"pyrin_activation_threshold": 0.3}),
                      ("v1 priming 1.0->1.2 (drive 0.7)", {"priming_input": 1.2}),
                      ("v1 priming 1.0->1.8 (drive 1.3)", {"priming_input": 1.8}),
                      ("v1 threshold 0.5->-0.3 (drive 1.3)", {"pyrin_activation_threshold": -0.3})):
        p = dict(b1)
        p.update(over)
        t, tr = O1.simulate(p)
        r = O1.scenario_readouts(t, tr)
        srows.append({"case": lab, "peak_P_active": r["peak_P_active"],
                      "peak_IL1B": r["peak_IL1B_external"], "time_to_50pct_IL1B": r["time_to_50pct_IL1B"],
                      "IL18_over_IL1B_peak": r["peak_IL18_external"] / r["peak_IL1B_external"]})
    pd.DataFrame(srows).to_csv(f"{TAB}/R34_ode1_symmetry_check.csv", index=False)

    print(pd.DataFrame(rows).set_index("scenario")[KEY_OUT].round(3).to_string())
    print(pd.DataFrame(mrows).set_index("perturbation")[[c for c in pd.DataFrame(mrows).columns if c.endswith("pct_change")]].round(1).to_string())
    print("gate open in", int(samples.gate_open.sum()), "of", N, "LHS sets")
    top = G.assign(a=G.PRCC.abs()).sort_values("a", ascending=False).groupby("output").head(4)
    print(top[["output", "parameter", "PRCC", "spearman_rho"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
