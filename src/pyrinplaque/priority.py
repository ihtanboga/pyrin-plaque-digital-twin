"""Caveat-aware Pyrin-Plaque Priority Score.

"Hypothesis prioritization for future validation" — NOT a drug-discovery ranking.
Fully decomposable: final = sum(positive_i * w_i) - sum(penalty_j * w_j).
Component values (0..1) are assigned transparently per candidate from Day1-4 evidence;
weights live in config/scoring_weights.yaml.
"""
import yaml
import pandas as pd

def load_weights(path="config/scoring_weights.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)

def score_candidate(cand, weights):
    """cand: dict with 'positives' and 'penalties' sub-dicts (each component 0..1).
    Returns (raw_positive, total_penalty, final, components_df_rows)."""
    pos_w = weights["positive_components"]
    pen_w = weights["penalty_components"]
    rows = []
    raw_pos = 0.0
    for k, w in pos_w.items():
        v = float(cand["positives"].get(k, 0.0))
        pts = v * w
        raw_pos += pts
        rows.append({"candidate_id": cand["candidate_id"], "component": k, "kind": "positive",
                     "value_0_1": round(v, 3), "weight": w, "points": round(pts, 2)})
    tot_pen = 0.0
    for k, w in pen_w.items():
        v = float(cand["penalties"].get(k, 0.0))
        pts = v * w
        tot_pen += pts
        rows.append({"candidate_id": cand["candidate_id"], "component": k, "kind": "penalty",
                     "value_0_1": round(v, 3), "weight": w, "points": round(-pts, 2)})
    final = raw_pos - tot_pen
    return round(raw_pos, 2), round(tot_pen, 2), round(final, 2), rows

def grade(final, weights):
    bins = weights["confidence_grade_bins"]
    if final >= bins["high"]:
        return "high"
    if final >= bins["moderate"]:
        return "moderate"
    return "low"

def build_table(candidates, weights):
    """Returns (priority_df, components_df)."""
    prows, crows = [], []
    for c in candidates:
        raw, pen, final, comps = score_candidate(c, weights)
        crows.extend(comps)
        prows.append({
            "candidate_id": c["candidate_id"], "candidate_type": c["candidate_type"],
            "candidate_name": c["candidate_name"], "biological_axis": c["biological_axis"],
            "raw_positive_score": raw, "total_penalty": pen, "final_priority_score": final,
            "confidence_grade": grade(final, weights),
            "strongest_supporting_evidence": c["strongest_supporting_evidence"],
            "strongest_counterevidence": c["strongest_counterevidence"],
            "recommended_next_validation": c["recommended_next_validation"],
            "safe_claim": c["safe_claim"], "claims_to_avoid": c["claims_to_avoid"],
        })
    pdf = pd.DataFrame(prows).sort_values("final_priority_score", ascending=False).reset_index(drop=True)
    return pdf, pd.DataFrame(crows)
