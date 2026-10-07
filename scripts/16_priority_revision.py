#!/usr/bin/env python
"""Revision: caveat-aware priority score re-run.

1. Rebuilds the v1 candidate inputs (handoff/candidates.json was not in the repository) from
   results/tables/priority_score_components.csv and pyrin_plaque_priority_table.csv and checks
   that the v1 scores are reproduced exactly.
2. Applies a documented set of component changes driven by the revision analyses, writes the
   change log (R40), the revised priority table (R41) and components (R42), and saves the
   revised inputs to config/priority_candidates_revised.json.

Weights and grade thresholds are unchanged (config/scoring_weights.yaml).
Usage (repo root): python -I scripts/16_priority_revision.py
"""
import json
import sys
import pandas as pd

sys.path.insert(0, "src")
from pyrinplaque import priority as P

TAB = "results/revision/tables"
W = P.load_weights("config/scoring_weights.yaml")

RENAME = {
    "L1-01": ("Myeloid-compartment pyrin-backbone-high state", "pyrin cell-state localization"),
    "L1-02": ("Classical-monocyte MEFV-detection-high subset", "pyrin cell-state localization"),
    "L1-04": ("Plaque-core MEFV-detection signal (direction reversed after sample-label correction)", "plaque-core localization"),
    "L2-05": ("Shared effector arms: cytokine (CASP1/IL1B/IL18) and lysis (GSDMD/GSDME/NINJ1)", "shared effector output"),
}
# (candidate, component, new value, reason)
CHANGES = [
    ("L1-01", "random_control_penalty", None, "set from recomputed expression-matched null (R10a)"),
    ("L1-02", "provisional_annotation_penalty", 0.5,
     "Subset defined by canonical monocyte markers after Harmony integration; not author-curated"),
    ("L1-02", "singlecell_compartment_score", 0.8,
     "MEFV detection highest in classical monocytes; monocyte vs macrophage rate ratio per UMI ~4.9 (R12a, R12e)"),
    ("L1-02", "sparse_expression_penalty", 0.4, "Detection ~12% in classical monocytes, versus ~1% of all cells"),
    ("L1-04", "MEFV_detection_score", 0.2,
     "After correcting the AC/PA labels, myeloid MEFV detection is LOWER in core in 3/3 patients (R13a)"),
    ("L1-04", "singlecell_compartment_score", 0.4, "Regional difference largely explained by monocyte composition (R13b-c)"),
    ("L2-05", "myeloid_confounding_penalty", 0.8,
     "Cytokine arm tracks the myeloid-marker score (r=0.96) and is ~95% attenuated after residualization (R22)"),
    ("L2-01", "ODE_leverage_score", 0.5,
     "Unchanged value; v2 model: pyrin gate governs output timing, priming governs IL-1b magnitude (R31-R33)"),
]


# revised narrative fields (supporting/counter-evidence etc.) where the revision changed the evidence
TEXT = {
    "L2-01": {"strongest_counterevidence": "MEFV sparse in scRNA-seq (about 1% of cells); no genotype data; in the revised two-signal model the threshold mainly shifts output timing while priming scales IL-1beta magnitude, so neither lever dominates by construction."},
    "L1-01": {"strongest_supporting_evidence": "PYRIN_BACKBONE highest in the myeloid compartment in all 3 patients and all 6 samples; robust to doublet removal, SoupX/DecontX ambient correction and depth adjustment; survives shared-downstream exclusion (myeloid first for PYRIN_FULL and PYRIN_BACKBONE).",
              "strongest_counterevidence": "Expression-matched random gene sets: 86th percentile (suggestive, not exceptional); specificity delta first in only 2 of 3 patients; compartment is heterogeneous (monocytes, macrophages, dendritic cells); bulk support confounded by myeloid abundance; n = 3."},
    "L1-02": {"strongest_supporting_evidence": "Classical monocytes show the highest MEFV detection of all myeloid subsets (11.7%; first in each patient) with monocyte-marker co-expression (FCN1, S100A8, VCAN, SELL); monocyte vs macrophage rate ratio per UMI 4.9; robust to ambient correction.",
              "strongest_counterevidence": "Subsets defined computationally (no author-curated labels, no protein data); detection varies across patients (5-21%); neutrophils poorly captured; three patients.",
              "recommended_next_validation": "RNA in situ hybridization or spatial transcriptomics of MEFV with CD14/FCN1 and macrophage markers; replication in independent plaque single-cell datasets; pyrin protein in plaque monocytes.",
              "safe_claim": "Within the plaque myeloid compartment, MEFV transcripts are concentrated in classical monocytes, a candidate focus for validating pyrin-axis regulation.",
              "claims_to_avoid": "pyrin-active cell type; macrophage-specific pyrin activity; causal; clinically actionable subset"},
    "L2-03": {"strongest_supporting_evidence": "PKN1/PKN2 phosphorylate pyrin; 14-3-3 (YWHAB/YWHAZ) binding keeps pyrin inhibited, a pyrin-specific switch represented by the gate and inhibition parameters of the two-signal model.",
              "strongest_counterevidence": "Low MEFV context; in the model the gate governs timing more than magnitude; no direct plaque perturbation data."},
    "L1-03": {"strongest_supporting_evidence": "NLRP3_BACKBONE highest in T/NK cells; NLRP3 in atherosclerosis is well established (CANTOS, colchicine trials)."},
    "L1-05": {"strongest_supporting_evidence": "PYRIN_BACKBONE higher in the unstable region in 4/4 plaques (paired d_z 1.44).",
              "strongest_counterevidence": "Residualization on the myeloid-marker score reduces the difference from +0.72 to +0.14 (about 80% attenuation); every module, including NLRP3 and both effector arms, is attenuated by at least 80%; bulk effect sizes track module myeloid specificity."},
    "L2-05": {"strongest_supporting_evidence": "Well-established, tractable shared effectors (IL-1beta: canakinumab/CANTOS); the cytokine arm shows the largest bulk shift.",
              "strongest_counterevidence": "Not pyrin-specific by construction (shared with NLRP3 and other inflammasomes); cytokine arm tracks the myeloid-marker score (r = 0.96) and is about 95% attenuated after adjustment; lysis arm inconsistent.",
              "safe_claim": "Shared cytokine and lysis arms are readouts of pathway output; they are not pyrin-specific and cannot attribute signal to pyrin.",
              "claims_to_avoid": "pyrin-specific target; IL-1beta inhibition proves a pyrin mechanism; evidence of pyroptosis"},
    "L1-04": {"strongest_supporting_evidence": "None after correction of the sample labels; the original signal (higher MEFV detection in core) reflected interchanged AC/PA labels.",
              "strongest_counterevidence": "With verified labels, myeloid MEFV detection is lower in core in 3/3 patients (rate ratio 0.38), largely paralleling fewer monocytes and more TREM2+ lipid-associated macrophages in core.",
              "recommended_next_validation": "Spatial transcriptomics of core versus adjacent tissue with MEFV and lineage probes in a larger paired cohort.",
              "safe_claim": "Regional differences in MEFV detection follow myeloid composition; there is no evidence of pyrin-program enrichment in plaque core.",
              "claims_to_avoid": "pyrin program enriched in plaque core; core drives pyrin activation"},
}


def rebuild_v1():
    comp = pd.read_csv("results/tables/priority_score_components.csv")
    tab = pd.read_csv("results/tables/pyrin_plaque_priority_table.csv").set_index("candidate_id")
    cands = []
    for cid, g in comp.groupby("candidate_id", sort=False):
        t = tab.loc[cid]
        cands.append({"candidate_id": cid, "candidate_type": t.candidate_type, "candidate_name": t.candidate_name,
                      "biological_axis": t.biological_axis,
                      "positives": dict(zip(g[g.kind == "positive"].component, g[g.kind == "positive"].value_0_1)),
                      "penalties": dict(zip(g[g.kind == "penalty"].component, g[g.kind == "penalty"].value_0_1)),
                      **{k: t[k] for k in ("strongest_supporting_evidence", "strongest_counterevidence",
                                           "recommended_next_validation", "safe_claim", "claims_to_avoid")}})
    pdf, _ = P.build_table(cands, W)
    chk = pdf.set_index("candidate_id")["final_priority_score"].sub(tab["final_priority_score"]).abs().max()
    assert chk < 1e-6, f"v1 scores not reproduced (max diff {chk})"
    return cands


def main():
    cands = rebuild_v1()
    null = pd.read_csv(f"{TAB}/R10a_random_geneset_null.csv", index_col=0)
    pct = float(null.loc["Macrophage/Myeloid", "empirical_percentile"])
    # random-control penalty scales with how unexceptional the observed score is (v1: 82nd pct -> 0.8)
    rc = round(min(1.0, max(0.2, (100 - pct) / 100 * 4.4)), 2)
    log = []
    byid = {c["candidate_id"]: c for c in cands}
    for cid, (name, axis) in RENAME.items():
        log.append({"candidate_id": cid, "component": "candidate_name", "v1": byid[cid]["candidate_name"],
                    "revised": name, "reason": "Renamed to match the revised evidence"})
        byid[cid]["candidate_name"], byid[cid]["biological_axis"] = name, axis
    for cid, compn, val, why in CHANGES:
        c = byid[cid]
        kind = "positives" if compn in c["positives"] else "penalties"
        if val is None:
            val = rc
            why = f"{why}: myeloid PYRIN_BACKBONE at the {pct:.0f}th percentile of 200 matched sets"
        log.append({"candidate_id": cid, "component": compn, "v1": c[kind][compn], "revised": val, "reason": why})
        c[kind][compn] = val
    import re as _re
    for c in cands:
        for fld, val in TEXT.get(c["candidate_id"], {}).items():
            log.append({"candidate_id": c["candidate_id"], "component": fld, "v1": c[fld], "revised": val,
                        "reason": "Narrative updated to the revised evidence"})
            c[fld] = val
        for fld in ("strongest_supporting_evidence", "strongest_counterevidence", "recommended_next_validation", "safe_claim"):
            c[fld] = _re.sub(r"Day-[0-9.]+:?\s*", "", str(c[fld])).strip()
    pdf, cdf = P.build_table(cands, W)
    pd.DataFrame(log).to_csv(f"{TAB}/R40_priority_change_log.csv", index=False)
    pdf.to_csv(f"{TAB}/R41_priority_table_revised.csv", index=False)
    cdf.to_csv(f"{TAB}/R42_priority_components_revised.csv", index=False)
    with open("config/priority_candidates_revised.json", "w") as f:
        json.dump(cands, f, indent=1, default=str)
    pd.set_option("display.width", 220)
    print(pdf[["candidate_id", "candidate_name", "raw_positive_score", "total_penalty", "final_priority_score",
               "confidence_grade"]].to_string(index=False))


if __name__ == "__main__":
    main()
