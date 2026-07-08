#!/usr/bin/env python
"""Day-5 Part D: caveat-aware Pyrin-Plaque Priority Score.
"Hypothesis prioritization for future validation" — NOT a drug-discovery ranking.
Fully decomposable; weights in config/scoring_weights.yaml; candidates in handoff/candidates.json.

Usage: python scripts/06_priority_score.py
"""
import sys, json
sys.path.insert(0, "src")
from pyrinplaque import priority as P

def main():
    weights = P.load_weights("config/scoring_weights.yaml")
    candidates = json.load(open("handoff/candidates.json"))
    pdf, cdf = P.build_table(candidates, weights)
    pdf.to_csv("results/tables/pyrin_plaque_priority_table.csv", index=False)
    cdf.to_csv("results/tables/priority_score_components.csv", index=False)
    print("priority table:", pdf.shape, "| top:", pdf.iloc[0].candidate_id, pdf.iloc[0].final_priority_score)

if __name__ == "__main__":
    main()
