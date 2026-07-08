# PYRIN-PLAQUE — Streamlit Demo App

A **minimal, static** demo. It reads pre-computed figures and tables only — **no new analysis, no raw data, no heavy compute.**

## Run
```bash
pip install streamlit    # if not already installed
streamlit run app/streamlit_app.py
```

## Pages
1. Overview / safe headline
2. Workflow (Fig 1)
3. Single-cell result (Fig 2 + specificity delta)
4. Bulk + caveat (Fig 3 + myeloid confounding)
5. ODE scenarios (Fig 4 + readouts)
6. Priority score (Fig 5 + priority table)
7. Reviewer audit / caveats

If Streamlit is unavailable, the **static README is the MVP** — every figure and table the app shows is in `results/figures/` and `results/tables/`.

## Runtime status
The app is **syntax-verified but not runtime-launched in this environment** (Streamlit not installed here). It reads only static precomputed figures/tables. **If Streamlit is unavailable, the static README + figures are the MVP** — no functionality is lost.
