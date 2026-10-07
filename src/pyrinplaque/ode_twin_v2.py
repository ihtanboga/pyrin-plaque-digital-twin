"""Two-signal pyrin inflammasome output model (revision v2).

DIMENSIONLESS, hypothesis-generating; NOT calibrated to patients, drugs or MEFV variants.

Why a v2 model
--------------
In v1 the priming input and the pyrin activation threshold entered the system only through
their difference, max(0, priming - threshold). The two levers were therefore interchangeable
(an equal shift of either produced identical trajectories), and IL-1beta, IL-18 and the
rupture proxy were driven by the same GSDMD-N term with the same decay, so they were exact
scalar multiples of one another. v2 separates the two signals and the two output arms:

* Signal 1 (priming, `priming_input`) acts on transcription only: it sets the size of the
  pro-IL-1beta pool and, more weakly, the MEFV/pyrin pool. pro-IL-18 is constitutive.
* Signal 2 (`signal2_input`, a RhoA-inactivating stimulus that relieves PKN/14-3-3
  restraint) gates pyrin activation through max(0, signal2 - threshold) * (1 - inhibition).
  The threshold is defined relative to signal 2, so the two are interchangeable by
  construction; this is stated, not hidden.
* Caspase-1 is catalytic. It matures pro-IL-1beta / pro-IL-18 (cytokine arm) and cleaves
  GSDMD (lysis arm).
* Mature cytokines leave surviving cells through GSDMD pores (sub-lytic release) and are
  released in bulk when cells lyse. Lysis is a steep (Hill) function of pore load, a
  stand-in for NINJ1-dependent plasma-membrane rupture, so cytokine release and lysis are
  separable.

State vector (population-level amounts):
    0 P_inactive     inactive (phosphorylated, 14-3-3-bound) pyrin
    1 P_active       active pyrin
    2 ASC            ASC specks
    3 CASP1          active caspase-1
    4 GSDMD_N        GSDMD-N pores
    5 proIL1B        intracellular pro-IL-1beta
    6 matIL1B_in     intracellular mature IL-1beta
    7 IL1B_out       extracellular mature IL-1beta
    8 proIL18        intracellular pro-IL-18
    9 matIL18_in     intracellular mature IL-18
   10 IL18_out       extracellular mature IL-18
   11 lysed          lysed fraction of the cell population (0-1)
   12 IL1B_pore      cumulative IL-1beta released through pores (bookkeeping)
   13 IL1B_lysis     cumulative IL-1beta released by lysis (bookkeeping)
"""
import numpy as np
import yaml

STATE_NAMES = ["P_inactive", "P_active", "ASC", "CASP1", "GSDMD_N",
               "proIL1B", "matIL1B_in", "IL1B_out",
               "proIL18", "matIL18_in", "IL18_out",
               "lysed", "IL1B_pore", "IL1B_lysis"]
IDX = {n: i for i, n in enumerate(STATE_NAMES)}


def load_params(path="config/ode_parameters_v2.yaml", scenario="baseline"):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return dict(cfg[scenario])


def _sat(x, ceiling):
    """Bounded production term x/(1+x/ceiling)."""
    return x / (1.0 + x / ceiling)


def _hill(x, k, n):
    x = max(x, 0.0)
    return x ** n / (k ** n + x ** n)


def rhs(t, y, p):
    (P_i, P_a, ASC, C1, G, proB, mB, B_out,
     pro18, m18, I18_out, L, B_pore, B_lys) = y
    sat = p["saturation"]
    alive = max(1.0 - L, 0.0)
    # lysis hazard per surviving cell (NINJ1-like step, steep in pore load)
    lam = p["lysis_rate"] * _hill(G, p["lysis_pore_threshold"], p["lysis_hill"])
    dL = lam * alive

    # signal 2 gate (threshold is relative to signal 2 by construction)
    drive = max(0.0, p["signal2_input"] - p["pyrin_activation_threshold"])
    act = p["pyrin_activation_rate"] * drive * P_i * (1.0 - p["inhibition_factor"])
    # signal 1: transcriptional supply in surviving cells
    pyrin_syn = alive * (p["pyrin_basal_synthesis"] + p["pyrin_priming_synthesis"] * p["priming_input"])
    proB_syn = alive * p["proIL1B_priming_synthesis"] * p["priming_input"]
    pro18_syn = alive * p["proIL18_constitutive_synthesis"]

    dP_i = pyrin_syn - act + p["decay"] * P_a - p["pyrin_turnover"] * P_i - lam * P_i
    dP_a = act - p["decay"] * P_a - p["ASC_recruitment_rate"] * _sat(P_a, sat) - lam * P_a
    dASC = p["ASC_recruitment_rate"] * _sat(P_a, sat) - p["decay"] * ASC - lam * ASC
    dC1 = p["CASP1_activation_rate"] * _sat(ASC, sat) - p["decay"] * C1 - lam * C1
    # caspase-1 is catalytic (not consumed)
    dG = alive * p["GSDMD_cleavage_rate"] * _sat(C1, sat) - p["pore_repair"] * G - lam * G

    cleave = p["cytokine_maturation_rate"] * _sat(C1, sat)
    pore_rel = p["pore_release_rate"] * G
    dproB = proB_syn - cleave * proB - p["pro_turnover"] * proB - lam * proB
    dmB = cleave * proB - pore_rel * mB - p["mature_turnover"] * mB - lam * mB
    dB_out = pore_rel * mB + lam * mB - p["clearance"] * B_out
    dpro18 = pro18_syn - cleave * pro18 - p["pro_turnover"] * pro18 - lam * pro18
    dm18 = cleave * pro18 - pore_rel * m18 - p["mature_turnover"] * m18 - lam * m18
    dI18_out = pore_rel * m18 + lam * m18 - p["clearance"] * I18_out
    return [dP_i, dP_a, dASC, dC1, dG, dproB, dmB, dB_out,
            dpro18, dm18, dI18_out, dL, pore_rel * mB, lam * mB]


def initial_state(p):
    """Primed steady state before signal 2 (t=0 is signal-2 onset)."""
    y0 = np.zeros(len(STATE_NAMES))
    y0[IDX["P_inactive"]] = (p["pyrin_basal_synthesis"] + p["pyrin_priming_synthesis"] * p["priming_input"]) / p["pyrin_turnover"]
    y0[IDX["proIL1B"]] = p["proIL1B_priming_synthesis"] * p["priming_input"] / p["pro_turnover"]
    y0[IDX["proIL18"]] = p["proIL18_constitutive_synthesis"] / p["pro_turnover"]
    return y0


def simulate(params, t_span=(0, 50), n_points=1001):
    from scipy.integrate import solve_ivp
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(rhs, t_span, initial_state(params), t_eval=t_eval, args=(params,),
                    method="LSODA", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = np.clip(sol.y.T, 0.0, None)  # documented clipping of tiny negative excursions
    traj = {"t": sol.t}
    for i, name in enumerate(STATE_NAMES):
        traj[name] = Y[:, i]
    return sol.t, traj


def _trap(y, t):
    return float(np.trapezoid(y, t))


def _time_to_frac(t, y, frac):
    y = np.asarray(y)
    ymax = y.max()
    if ymax <= 1e-12:
        return float("nan")
    idx = int(np.argmax(y >= frac * ymax))
    return float(t[idx])


def readouts(t, traj):
    r = {}
    r["peak_P_active"] = float(traj["P_active"].max())
    r["time_to_50pct_P_active"] = _time_to_frac(t, traj["P_active"], 0.5)
    for name, key in (("IL1B", "IL1B_out"), ("IL18", "IL18_out")):
        r[f"peak_{name}_out"] = float(traj[key].max())
        r[f"{name}_AUC"] = _trap(traj[key], t)
        r[f"time_to_50pct_peak_{name}"] = _time_to_frac(t, traj[key], 0.5)
    r["IL1B_to_IL18_AUC_ratio"] = r["IL1B_AUC"] / r["IL18_AUC"] if r["IL18_AUC"] > 0 else float("nan")
    r["lysed_fraction_end"] = float(traj["lysed"][-1])
    r["time_to_50pct_final_lysis"] = _time_to_frac(t, traj["lysed"], 0.5)
    tot = traj["IL1B_pore"][-1] + traj["IL1B_lysis"][-1]
    r["IL1B_released_via_pores_fraction"] = float(traj["IL1B_pore"][-1] / tot) if tot > 0 else float("nan")
    # IL-1beta already released before 10% of cells have lysed (sub-lytic release)
    rel = traj["IL1B_pore"] + traj["IL1B_lysis"]
    i10 = int(np.argmax(traj["lysed"] >= 0.10)) if (traj["lysed"] >= 0.10).any() else len(t) - 1
    r["IL1B_release_before_10pct_lysis_fraction"] = float(rel[i10] / tot) if tot > 0 else float("nan")
    return r
