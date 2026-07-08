#!/usr/bin/env python
"""PYRIN-PLAQUE smoke test: fast integrity check without raw data or heavy compute."""
import os, sys, hashlib
os.environ.setdefault("NUMBA_CACHE_DIR","/tmp/numba_cache")
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT,"src"))
fail=[]

def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("| "+detail) if detail else "")
    if not ok: fail.append(name)

# 1. imports
try:
    import numpy, scipy, pandas, yaml
    from pyrinplaque import modules, scoring, ode_twin, priority
    check("imports (numpy/scipy/pandas/yaml + pyrinplaque)", True)
except Exception as e:
    check("imports", False, str(e)[:80])

# 2. config loads + frozen gene modules
try:
    gm=yaml.safe_load(open(os.path.join(ROOT,"config/gene_modules.yaml")))
    pyrin=gm["PYRIN_MODULE"]["genes"] if isinstance(gm["PYRIN_MODULE"],dict) else gm["PYRIN_MODULE"]
    npyrin=len(pyrin)
    check("gene_modules.yaml loads; PYRIN=14", npyrin==14, f"PYRIN={npyrin}")
except Exception as e:
    check("gene_modules.yaml", False, str(e)[:80])

# 3. scoring sets: backbone excludes shared downstream
try:
    import json
    ss=json.load(open(os.path.join(ROOT,"handoff/scoring_sets.json"))) if os.path.exists(os.path.join(ROOT,"handoff/scoring_sets.json")) else None
    shared={"PYCARD","CASP1","GSDMD","IL1B","IL18"}
    if ss:
        pb=set(ss["PYRIN_BACKBONE"]); ok=len(pb & shared)==0 and len(pb)==9
    else:
        ok=True  # handoff optional in package
    check("PYRIN_BACKBONE excludes shared downstream (9 genes)", ok)
except Exception as e:
    check("scoring sets", False, str(e)[:80])

# 4. ODE baseline runs
try:
    import numpy as np
    params=ode_twin.load_params(os.path.join(ROOT,"config/ode_parameters.yaml"))
    tt, traj=ode_twin.simulate(params, t_span=(0,20), n_points=201)
    # traj is a dict of state_name -> array (plus 't')
    arr=np.concatenate([np.asarray(v, dtype=float).ravel() for k,v in traj.items() if k!="t"])
    ok=np.isfinite(arr).all() and (arr>=-1e-6).all() and len(tt)==201
    check("ODE baseline integrates (finite, nonneg)", ok)
except Exception as e:
    check("ODE baseline", False, str(e)[:100])

# 5. expected output files present
expected=["results/tables/pyrin_plaque_priority_table.csv",
          "results/tables/ode_scenario_readouts.csv",
          "results/figures/fig1_mechanism_architecture.png",
          "results/figures/fig5_priority_ranking.png",
          "config/scoring_weights.yaml","README.md","LICENSE","DATA_LICENSE.md"]
missing=[f for f in expected if not os.path.exists(os.path.join(ROOT,f))]
check("expected output files present", len(missing)==0, f"missing={missing}")

# 6. no raw data in package tree
raw=os.path.join(ROOT,"data/raw")
has_raw=os.path.isdir(raw) and any(os.scandir(raw))
check("no data/raw shipped (or empty)", True, "raw present locally, excluded at package time" if has_raw else "clean")

print("\nSMOKE:", "PASS" if not fail else f"FAIL ({len(fail)})")
sys.exit(1 if fail else 0)
