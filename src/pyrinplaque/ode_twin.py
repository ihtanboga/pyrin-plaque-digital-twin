"""Minimal pyrin-adapted pyroptosis ODE scaffold (Day-4 baseline).

DIMENSIONLESS, hypothesis-generating. NOT calibrated to patients; NOT a clinical prediction model.
The NLRP3 pyroptosis ODE literature is used as a structural TEMPLATE only — this does not assert
that pyrin and NLRP3 are mechanistically equivalent.

State vector y (8 variables):
    0 P_inactive    inactive pyrin sensor pool
    1 P_active      active pyrin sensor
    2 ASC_complex   ASC speck / adaptor complex
    3 CASP1_active  active caspase-1
    4 GSDMD_N       cleaved gasdermin-D N-terminal (pore-forming)
    5 IL1B_external externalized IL-1beta (proxy)
    6 IL18_external externalized IL-18 (proxy)
    7 rupture_proxy membrane-rupture / pyroptosis proxy
"""
import numpy as np
import yaml

STATE_NAMES = ["P_inactive", "P_active", "ASC_complex", "CASP1_active",
               "GSDMD_N", "IL1B_external", "IL18_external", "rupture_proxy"]

def load_params(path="config/ode_parameters.yaml", scenario="baseline"):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return dict(cfg[scenario])

def _sat(x, ceiling):
    """Hill-like saturation: x/(1+x/ceiling), keeps production bounded and nonnegative."""
    return x / (1.0 + x / ceiling)

def rhs(t, y, p):
    """Right-hand side of the ODE system. y nonnegative by construction of the terms."""
    P_i, P_a, ASC, C1, GN, IL1, IL18, RUP = y
    sat = p["saturation"]
    # pyrin activation gated by priming above threshold, modulated by inhibition
    drive = max(0.0, p["priming_input"] - p["pyrin_activation_threshold"])
    act = p["pyrin_activation_rate"] * drive * P_i * (1.0 - p["inhibition_factor"])
    dP_i = -act + p["decay"] * P_a
    dP_a = act - p["decay"] * P_a - p["ASC_recruitment_rate"] * _sat(P_a, sat)
    dASC = p["ASC_recruitment_rate"] * _sat(P_a, sat) - p["CASP1_activation_rate"] * _sat(ASC, sat) - p["decay"] * ASC
    dC1  = p["CASP1_activation_rate"] * _sat(ASC, sat) - p["GSDMD_cleavage_rate"] * _sat(C1, sat) - p["decay"] * C1
    dGN  = p["GSDMD_cleavage_rate"] * _sat(C1, sat) - p["decay"] * GN
    dIL1 = p["cytokine_release_rate_IL1B"] * _sat(GN, sat) - p["decay"] * IL1
    dIL18= p["cytokine_release_rate_IL18"] * _sat(GN, sat) - p["decay"] * IL18
    dRUP = p["rupture_rate"] * _sat(GN, sat) - p["decay"] * RUP
    return [dP_i, dP_a, dASC, dC1, dGN, dIL1, dIL18, dRUP]

def initial_state():
    y0 = np.zeros(len(STATE_NAMES))
    y0[0] = 1.0  # start with a full inactive pyrin pool
    return y0

def simulate(params, t_span=(0, 50), n_points=501, clip_negative=True):
    """Integrate the baseline system. Returns (t, DataFrame-ready dict of trajectories)."""
    from scipy.integrate import solve_ivp
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(rhs, t_span, initial_state(), t_eval=t_eval, args=(params,),
                    method="LSODA", rtol=1e-8, atol=1e-10)
    Y = sol.y.T
    if clip_negative:
        Y = np.clip(Y, 0.0, None)  # documented clipping: tiny negative numerical excursions -> 0
    traj = {"t": sol.t}
    for i, name in enumerate(STATE_NAMES):
        traj[name] = Y[:, i]
    return sol.t, traj


# ---- Day-5 scenario readouts ----
def _auc(t, y):
    import numpy as np
    trap = getattr(np, "trapezoid", getattr(np, "trapz", None))
    return float(trap(y, t))

def _time_to_frac(t, y, frac):
    """First time y reaches frac*max(y). Returns nan if never."""
    import numpy as np
    y = np.asarray(y); ymax = y.max()
    if ymax <= 0: return float("nan")
    thr = frac * ymax
    idx = np.argmax(y >= thr)
    return float(t[idx]) if y[idx] >= thr else float("nan")

def _time_to_abs(t, y, level):
    """First time y reaches an absolute level. nan if never."""
    import numpy as np
    y = np.asarray(y)
    idx = np.argmax(y >= level)
    return float(t[idx]) if (y >= level).any() else float("nan")

def scenario_readouts(t, traj):
    """13 readouts from a trajectory dict (keys = STATE_NAMES + 't')."""
    r = {}
    r["peak_P_active"] = float(max(traj["P_active"]))
    r["time_to_50pct_P_active"] = _time_to_frac(t, traj["P_active"], 0.5)
    r["peak_ASC_complex"] = float(max(traj["ASC_complex"]))
    r["peak_CASP1_active"] = float(max(traj["CASP1_active"]))
    r["peak_GSDMD_N"] = float(max(traj["GSDMD_N"]))
    r["peak_IL1B_external"] = float(max(traj["IL1B_external"]))
    r["peak_IL18_external"] = float(max(traj["IL18_external"]))
    r["IL1B_AUC"] = _auc(t, traj["IL1B_external"])
    r["IL18_AUC"] = _auc(t, traj["IL18_external"])
    r["time_to_50pct_IL1B"] = _time_to_frac(t, traj["IL1B_external"], 0.5)
    r["time_to_50pct_IL18"] = _time_to_frac(t, traj["IL18_external"], 0.5)
    r["time_to_rupture_proxy_0p5"] = _time_to_abs(t, traj["rupture_proxy"], 0.5 * max(traj["rupture_proxy"]))
    r["rupture_proxy_AUC"] = _auc(t, traj["rupture_proxy"])
    return r

def load_scenarios(path="config/simulation_scenarios.yaml"):
    import yaml
    with open(path) as f:
        return yaml.safe_load(f)
