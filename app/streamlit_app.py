"""PYRIN-PLAQUE Digital Twin — minimal static demo app.
No new analysis: reads pre-computed tables and figures only.
Run: streamlit run app/streamlit_app.py
"""
import os
import streamlit as st
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "results", "figures")
TAB = os.path.join(ROOT, "results", "tables")

st.set_page_config(page_title="PYRIN-PLAQUE Digital Twin", layout="wide")

SAFE_HEADLINE = (
    "PYRIN-PLAQUE identifies a **suggestive myeloid pyrin-permissiveness** pattern and "
    "prioritizes the **MEFV/pyrin threshold axis** for future validation, while showing that "
    "much of the bulk signal is **myeloid-abundance-confounded** and that **generic inflammatory "
    "priming may dominate** over pyrin-specific threshold effects. It is an **auditable hypothesis "
    "engine, not a validated-target or clinical claim.**"
)

def show_fig(name, caption):
    p = os.path.join(FIG, name)
    if os.path.exists(p):
        st.image(p, caption=caption, use_container_width=True)
    else:
        st.info(f"Figure not found in package: {name}")

def show_table(name):
    p = os.path.join(TAB, name)
    if os.path.exists(p):
        st.dataframe(pd.read_csv(p), use_container_width=True)
    else:
        st.info(f"Table not found in package: {name}")

page = st.sidebar.radio("Pages", [
    "1. Overview", "2. Workflow (Fig 1)", "3. Single-cell", "4. Bulk + caveat",
    "5. ODE scenarios", "6. Priority score", "7. Reviewer audit / caveats"])

st.sidebar.markdown("---")
st.sidebar.caption("Public data only · hypothesis-generating · auditable Claude Science workflow")

if page == "1. Overview":
    st.title("PYRIN-PLAQUE Digital Twin")
    st.markdown(SAFE_HEADLINE)
    st.markdown("**This is:** public-data-only · hypothesis-generating · reproducible · auditable.")
    st.markdown("**This is NOT:** validated target discovery · clinical prediction · drug recommendation.")
elif page == "2. Workflow (Fig 1)":
    st.header("Mechanism + Claude Science workflow")
    show_fig("fig1_mechanism_architecture.png", "Fig 1 — Pyrin vs NLRP3 architecture + auditable workflow")
elif page == "3. Single-cell":
    st.header("Single-cell cell-state program map (GSE159677)")
    show_fig("fig2_cellstate_program_map.png", "Fig 2 — PYRIN_BACKBONE highest in Macrophage/Myeloid")
    st.caption("MEFV sparse overall (~0.97%) but myeloid-enriched (4.3%). Random control modest (~82nd pctile).")
    show_table("pyrin_specificity_delta.csv")
elif page == "4. Bulk + caveat":
    st.header("Bulk validation + myeloid confounding (GSE120521)")
    show_fig("fig3_bulk_validation.png", "Fig 3 — PYRIN_BACKBONE higher in unstable 4/4, but myeloid-confounded")
    st.warning("Myeloid residualization collapses the effect: +0.72 -> +0.14 (p about 0.69). Directionally supportive but confounded.")
    show_table("bulk_status_tests.csv")
elif page == "5. ODE scenarios":
    st.header("Dimensionless pyroptosis digital twin")
    show_fig("fig4_pyrin_ode_twin.png", "Fig 4 — 4 scenarios + sensitivity")
    st.info("Generic priming dominates the pyrin threshold. Hypothesis-generating only; not a clinical prediction.")
    show_table("ode_scenario_readouts.csv")
elif page == "6. Priority score":
    st.header("Caveat-aware priority score")
    show_fig("fig5_priority_ranking.png", "Fig 5 — MEFV/pyrin threshold axis top (moderate); no high-confidence target")
    show_table("pyrin_plaque_priority_table.csv")
elif page == "7. Reviewer audit / caveats":
    st.header("Reviewer-agent audit & caveats")
    st.markdown("Six adversarial gates (Day 1-6). Null and confounded results preserved, not hidden.")
    st.markdown("- Sparse MEFV · provisional labels · n=3 single-cell / n=4 bulk (FPKM)")
    st.markdown("- Bulk myeloid-confounded · AC/PA mixed · no genotype · ODE uncalibrated")
    st.markdown("- Integrated atlas not used as data source · shared downstream non-specific")
    st.caption("Full history: reports/reviewer_audit.md")
