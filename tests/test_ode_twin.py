import sys; sys.path.insert(0, "src")
import numpy as np
from pyrinplaque import ode_twin as O

def _run(inhibition=0.0):
    p = O.load_params("config/ode_parameters.yaml", "baseline")
    p["inhibition_factor"] = inhibition
    t, traj = O.simulate(p, t_span=(0, 50), n_points=501)
    return t, traj

def test_module_imports():
    assert hasattr(O, "simulate") and hasattr(O, "rhs")
    assert len(O.STATE_NAMES) == 8

def test_expected_columns():
    t, traj = _run()
    assert "t" in traj
    for s in O.STATE_NAMES:
        assert s in traj

def test_no_nan_inf():
    t, traj = _run()
    for s in O.STATE_NAMES:
        arr = np.asarray(traj[s])
        assert not np.isnan(arr).any()
        assert not np.isinf(arr).any()

def test_nonnegative():
    t, traj = _run()
    for s in O.STATE_NAMES:
        assert (np.asarray(traj[s]) >= -1e-9).all(), f"{s} negative"

def test_baseline_completes():
    t, traj = _run()
    assert len(t) == 501
    assert t[-1] == 50

def test_inhibition_reduces_cytokines():
    # tiny inhibition smoke check: higher inhibition -> less IL1B externalization at peak
    _, base = _run(inhibition=0.0)
    _, inhib = _run(inhibition=0.9)
    assert max(inhib["IL1B_external"]) < max(base["IL1B_external"]), "inhibition should lower peak IL1B"
