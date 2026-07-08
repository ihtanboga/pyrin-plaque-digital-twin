# Figure card — Fig 5 Priority ranking (decomposed)

- **Output:** results/figures/fig5_priority_ranking.png
- **Input:** results/tables/pyrin_plaque_priority_table.csv, priority_score_components.csv
- **Code:** scripts/06_priority_score.py; src/pyrinplaque/priority.py; weights config/scoring_weights.yaml; candidates handoff/candidates.json
- **Environment:** conda env `sc` (python 3.13.14)
- **Data source:** derived from Day1-4 result tables (evidence graph, single-cell scores, bulk tests, ODE).

**Shows:** Fully decomposed caveat-aware priority score for 11 candidates (5 cell-state L1 + 6 mechanism-node L2). Green = positive contributions (stacked), red = penalties (stacked left), ◆ = final score. Dotted lines: moderate≥35, high≥55.

**Interpretation:** The **MEFV/pyrin threshold axis (L2-01, 37.9, moderate)** and the **Macrophage/Myeloid pyrin-permissive state (L1-01, 33.1)** rank highest. **Nothing reaches the 'high' band — conservative by design.** Penalties do exactly what the caveat ledger requires: the AC-core signal (L1-04) is pushed to the bottom by `ac_pa_mixed_signal_penalty`, and the unstable-bulk signal (L1-05) is heavily discounted by `myeloid_confounding_penalty`. Shared-downstream (L2-05) gets zero pyrin-specificity credit. NLRP3 (L2-06) scores respectably but is labeled a saturated comparator.

**Limitations:** This is a hypothesis-prioritization score for future validation, NOT a validated target ranking, NOT a drug-discovery score, NOT a clinical decision tool. Component values are transparent expert-assigned mappings from Day1-4 evidence, not learned weights.
