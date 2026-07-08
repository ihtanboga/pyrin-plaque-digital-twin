# Figure card — Fig 1 Mechanism + Claude Science workflow

- **Output:** results/figures/fig1_mechanism_architecture.png
- **Input:** none — schematic built from frozen config/gene_modules.yaml module membership + project workflow.
- **Generating script/cell:** Day-6 Fig-1 matplotlib cell (FancyBboxPatch schematic)
- **Environment:** conda env `sc` (matplotlib 3.11.0, python 3.13.14)
- **Data source:** conceptual/schematic (no dataset).

**Shows:** (A) Pyrin backbone vs NLRP3 comparator backbone converging on shared downstream (ASC/PYCARD → CASP1 → GSDMD → IL-1β/IL-18), with shared nodes explicitly labeled 'zero pyrin-specificity credit'. (B) The seven-step Claude Science workflow (connectors → evidence graph → scRNA-seq → bulk → ODE → priority score → reviewer audit ×6) with caveat-ledger and provenance banners.

**Interpretation:** Orients the reader to the pyrin-vs-NLRP3 separation logic and the auditable public-data workflow.

**Limitations:** Schematic only; not data. Backbone membership reflects the frozen v0 modules; shared genes are comparators, not pyrin-specific.

**Related claims/caveats:** shared-downstream non-specificity (CV11); NLRP3 saturated comparator; whole caveat ledger (reports/day5_caveat_ledger.md).
