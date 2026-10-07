#!/usr/bin/env python
"""Fill the revised manuscript template from the revision tables, number citations by first
appearance, add tables and figure legends, and write (1) a clean .docx, (2) a tracked-changes
.docx against the originally submitted manuscript and (3) a values.json audit of every number.

Usage (repo root = revision/code):
  python -I tools/build_manuscript.py <template.md> <original.docx> <out_dir>
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "tools")
import docbuild as DB  # noqa: E402

T = "results/revision/tables"
CT = ["Macrophage/Myeloid", "Endothelial", "T/NK", "SMC/Fibroblast", "Mast", "B/Plasma"]
CT_LAB = {"Macrophage/Myeloid": "Myeloid", "Endothelial": "Endothelial", "T/NK": "T/NK",
          "SMC/Fibroblast": "Smooth muscle/fibroblast", "Mast": "Mast", "B/Plasma": "B/plasma"}


def r(name, **kw):
    return pd.read_csv(os.path.join(T, name), **kw)


def f1(x):
    return f"{x:.1f}"


def f2(x):
    return f"{x:.2f}"


def f3(x):
    return f"{x:.3f}"


def sgn(x, d=3):
    s = f"{x:+.{d}f}"
    return s.replace("-", "−")


def neg(x, d=2):
    return f"{x:.{d}f}".replace("-", "−")


def ci(lo, hi, d=2):
    return f"{neg(lo, d)} to {neg(hi, d)}"


def values():
    import anndata as ad
    V = {}
    rep = ad.read_h5ad("data/processed/gse159677_v1repro.h5ad", backed="r")
    V["N_LOADED"] = f"{int(rep.uns['v1repro']['n_loaded']):,}"
    V["N_QC"] = f"{rep.n_obs:,}"
    V["N_GENES"] = f"{rep.n_vars:,}"
    rep.file.close()
    s02 = r("R02_sample_barcode_verification_and_soup.csv")
    V["AGGR_RATIO_RANGE"] = f"{s02.aggregate_to_full_UMI_ratio_median.min():.2f}–{s02.aggregate_to_full_UMI_ratio_median.max():.2f}"
    V["SOUP_MEFV_RANGE"] = f"{1e6 * s02.MEFV_soup_fraction.min():.0f}–{1e6 * s02.MEFV_soup_fraction.max():.0f}"
    d3 = r("R03a_doublets_by_sample.csv", index_col=0)
    nd = int(d3.scDblFinder_doublets.sum())
    V["N_DBL"] = f"{nd:,}"
    V["PCT_DBL"] = f1(100 * nd / d3.n_QC.sum())
    V["DBL_RANGE"] = f"{d3.scDblFinder_rate_pct.min():.1f}–{d3.scDblFinder_rate_pct.max():.1f}%"
    V["N_SCRUB"] = f"{int(d3.scrublet_doublets.sum()):,}"
    ag = r("R03c_doublet_method_agreement.csv", index_col=0)
    V["SCRUB_OVERLAP_PCT"] = f"{100 * ag.loc[True, 'True'] / ag['True'].sum():.0f}"
    comp = r("R05_compartment_summary_revised.csv", index_col=0)
    V["N_SINGLETS"] = f"{int(comp.n_cells.sum()):,}"
    lab_ = {"Macrophage/Myeloid": "myeloid", "Endothelial": "endothelial", "T/NK": "T/NK", "SMC/Fibroblast": "smooth muscle/fibroblast",
            "Mast": "mast", "B/Plasma": "B/plasma"}
    V["CT_COUNTS"] = ", ".join(f"{int(comp.loc[c, 'n_cells']):,} {lab_[c]}" for c in comp.n_cells.sort_values(ascending=False).index) + " cells"
    amb = r("R09b_ambient_by_sample.csv")
    V["SOUPX_RANGE"] = f"{amb.soupx_rho_used.min():.3f}–{amb.soupx_rho_used.max():.3f}"
    V["DECONTX_RANGE"] = f"{100 * amb.decontx_median_contamination.min():.1f}–{100 * amb.decontx_median_contamination.max():.1f}%"
    pos = r("R07a_nUMI_MEFVpos_vs_neg.csv").set_index("subset")
    tot_pos = int(pos.loc["all compartments", "n_MEFV_pos"])
    V["MEFV_ALL_PCT"] = f2(100 * tot_pos / comp.n_cells.sum())
    for k, c in (("MYE", "Macrophage/Myeloid"), ("ENDO", "Endothelial"), ("SMC", "SMC/Fibroblast"), ("TNK", "T/NK")):
        V[f"MEFV_{k}_PCT"] = f2(comp.loc[c, "MEFV_pct"])
    d3b = r("R03b_doublets_by_compartment.csv", index_col=0)
    V["MEFV_MYE_QC_PCT"] = f2(d3b.loc["Macrophage/Myeloid", "MEFV_pct_all"])
    a9 = r("R09a_ambient_MEFV.csv", index_col=0)
    V["MEFV_MYE_SOUPX"] = f2(a9.loc["Macrophage/Myeloid", "MEFV_pct_SoupX"])
    V["MEFV_MYE_DECONTX"] = f2(a9.loc["Macrophage/Myeloid", "MEFV_pct_DecontX"])
    V["MEFV_MYE_AMBIENT"] = f"{a9.loc['Macrophage/Myeloid', 'expected_ambient_MEFV_pct']:.2f}"
    V["NUMI_MYE"] = f"{comp.loc['Macrophage/Myeloid', 'median_nUMI']:,.0f}"
    V["NUMI_TNK"] = f"{comp.loc['T/NK', 'median_nUMI']:,.0f}"
    g = r("R07b_cloglog_glm_depth_offset.csv")
    fe = g[g.model == "compartment + sample FE + offset log(nUMI)"].set_index("contrast")

    def rr(row):
        return f"{row.rate_ratio_per_UMI:.1f}; 95% CI {row.ci95_lo:.1f}–{row.ci95_hi:.1f}"
    V["RR_MYE_TNK"] = rr(fe.loc["Macrophage/Myeloid vs T/NK"])
    V["RR_ENDO_TNK"] = rr(fe.loc["Endothelial vs T/NK"])
    V["RR_MYE_NON"] = rr(g[g.contrast == "myeloid vs non-myeloid"].iloc[0])
    ac = g[g.contrast == "AC vs PA"].iloc[0]
    V["RR_AC_PA"] = f"{ac.rate_ratio_per_UMI:.2f} (95% CI {ac.ci95_lo:.2f}–{ac.ci95_hi:.2f})"
    ds = r("R08_fixed_depth_2000UMI.csv", index_col=0)
    V["DS_MYE"] = f2(ds.loc["Macrophage/Myeloid", "MEFV_pct_at_2000UMI"])
    V["DS_OTHER_MAX"] = f2(ds.drop(index="Macrophage/Myeloid").MEFV_pct_at_2000UMI.max())
    for k, c in (("MYE", "Macrophage/Myeloid"), ("ENDO", "Endothelial"), ("TNK", "T/NK"), ("SMC", "SMC/Fibroblast"),
                 ("MAST", "Mast"), ("B", "B/Plasma")):
        V[f"PB_{k}"] = f3(comp.loc[c, "PYRIN_BACKBONE_mean"])
        V[f"DELTA_{k}"] = sgn(comp.loc[c, "delta_cell_level_z"], 3)
    V["NB_TNK"] = f3(comp.loc["T/NK", "NLRP3_BACKBONE_mean"])
    V["NB_MYE"] = f3(comp.loc["Macrophage/Myeloid", "NLRP3_BACKBONE_mean"])
    V["DELTA_MYE_CM"] = sgn(comp.loc["Macrophage/Myeloid", "delta_compartment_mean_z"], 2)
    V["DELTA_MYE_PAT"] = sgn(comp.loc["Macrophage/Myeloid", "delta_patient_mean"], 2)
    V["DELTA_MYE_PAT_RANGE"] = f"{neg(comp.loc['Macrophage/Myeloid', 'delta_patient_min'])} to {neg(comp.loc['Macrophage/Myeloid', 'delta_patient_max'])}"
    pv = r("R06b_patient_compartment_values.csv").set_index(["patient", "cell_type"])
    V["DELTA_P2_SMC"] = sgn(pv.loc[("Patient 2", "SMC/Fibroblast"), "delta_cell_level_z"], 2)
    V["DELTA_P2_MYE"] = sgn(pv.loc[("Patient 2", "Macrophage/Myeloid"), "delta_cell_level_z"], 2)
    ols = r("R08b_depth_adjusted_module_scores.csv")
    o = ols[(ols.module == "PYRIN_BACKBONE_score") & ols.contrast.str.contains("minus")]
    V["OLS_RANGE"] = f"{abs(o.estimate).min():.2f}–{abs(o.estimate).max():.2f}"
    nl = r("R10a_random_geneset_null.csv", index_col=0)
    V["NULL_PCT"] = f"{nl.loc['Macrophage/Myeloid', 'empirical_percentile']:.0f}th"
    V["NULL_EXCESS"] = f"{nl.loc['Macrophage/Myeloid', 'excess_over_null']:.2f}, {nl.loc['Macrophage/Myeloid', 'excess_over_null_sd_units']:.2f} null SD"
    sub = r("R12a_myeloid_subsets.csv")
    agg = sub.groupby("subset").apply(lambda x: pd.Series({
        "n": x.n_cells.sum(), "MEFV_pct": np.average(x.MEFV_pct, weights=x.n_cells)}), include_groups=False)
    V["N_SUBCLUSTERS"] = str(sub.subcluster.nunique())
    V["N_SUBSETS"] = str(sub.subset.nunique())
    V["MEFV_CMONO_PCT"] = f1(agg.loc["Classical monocyte", "MEFV_pct"])
    V["N_CMONO"] = f"{int(agg.loc['Classical monocyte', 'n']):,}"
    V["MEFV_TREM2_PCT"] = f1(agg.loc["TREM2+ lipid-associated macrophage", "MEFV_pct"])
    V["N_TREM2"] = f"{int(agg.loc['TREM2+ lipid-associated macrophage', 'n']):,}"
    V["MEFV_CDC2_PCT"] = f1(agg.loc["cDC2", "MEFV_pct"])
    V["MEFV_LYVE1_PCT"] = f1(agg.loc["LYVE1+ resident-like macrophage", "MEFV_pct"])
    V["MEFV_NCMONO_PCT"] = f1(agg.loc["Non-classical monocyte", "MEFV_pct"])
    my = pd.read_csv(os.path.join(T, "R12c_sample_by_lineage_MEFV.csv"))
    import anndata as ad2
    mobs = ad2.read_h5ad("data/processed/gse159677_myeloid.h5ad", backed="r").obs
    cm = mobs[mobs.subset == "Classical monocyte"].groupby("patient", observed=True)["MEFV_detected"].mean() * 100
    V["CMONO_BY_PATIENT"] = "; ".join(f"{p.replace('Patient ', 'patient ')} {v:.1f}%" for p, v in cm.items())
    V["CMONO_RANGE"] = f"{cm.min():.0f}–{cm.max():.0f}%"
    V["N_NEUT"] = f"{int((mobs.subset == 'Neutrophil').sum())}"
    lin = r("R12b_myeloid_lineages.csv")
    l1 = lin[lin.definition == "lineage"].set_index("lineage")
    l2 = lin[lin.definition == "lineage_cell_level"].set_index("lineage")
    V["LIN_MONO"], V["LIN_MAC"], V["LIN_DC"] = f1(l1.loc["Monocyte", "MEFV_pct"]), f1(l1.loc["Macrophage", "MEFV_pct"]), f1(l1.loc["Dendritic cell", "MEFV_pct"])
    V["LIN_MONO_CELL"], V["LIN_MAC_CELL"] = f1(l2.loc["Monocyte", "MEFV_pct"]), f1(l2.loc["Macrophage", "MEFV_pct"])
    lg = r("R12e_lineage_cloglog_glm.csv").set_index("contrast")
    x = lg.loc["Monocyte vs Macrophage"]
    V["RR_MONO_MAC"] = f"{x.rate_ratio_per_UMI:.1f} (95% CI {x.ci95_lo:.1f}–{x.ci95_hi:.1f})"
    co = r("R12d_MEFV_coexpression_myeloid.csv").set_index("marker")
    V["COEXP_TEXT"] = "; ".join(f"*{m}* {co.loc[m, 'pct_in_MEFV_pos']:.0f}% versus {co.loc[m, 'pct_in_MEFV_neg']:.0f}%"
                                for m in ("FCN1", "S100A8", "VCAN", "SELL"))
    V["COEXP_TREM2"] = f"{co.loc['TREM2', 'pct_in_MEFV_pos']:.0f}% versus {co.loc['TREM2', 'pct_in_MEFV_neg']:.0f}%"
    v1 = r("R12f_v1_MEFVhigh_cluster_mapping.csv")
    V["V1CLUSTER_TEXT"] = f"{100 * v1.loc[v1.subset == 'Classical monocyte', 'fraction'].sum():.0f}% classical monocytes"
    V["V1CLUSTER_PCT"] = f"{100 * v1.loc[v1.subset == 'Classical monocyte', 'fraction'].sum():.0f}"
    ra = r("R13a_myeloid_AC_vs_PA_corrected.csv").set_index("metric")
    V["REG_MONO_PA"] = f1(ra.loc["Monocyte-lineage fraction (%)", "PA_mean_of_patients"])
    V["REG_MONO_AC"] = f1(ra.loc["Monocyte-lineage fraction (%)", "AC_mean_of_patients"])
    t2 = "TREM2+ lipid-associated macrophage fraction (%)"
    V["REG_TREM2_PA"], V["REG_TREM2_AC"] = f1(ra.loc[t2, "PA_mean_of_patients"]), f1(ra.loc[t2, "AC_mean_of_patients"])
    V["REG_MEFV_AC"] = f1(ra.loc["MEFV detection (%)", "AC_mean_of_patients"])
    V["REG_MEFV_PA"] = f1(ra.loc["MEFV detection (%)", "PA_mean_of_patients"])
    V["REG_MEFV_DIFF"] = sgn(ra.loc["MEFV detection (%)", "mean_AC_minus_PA"], 1)
    sd = r("R13b_region_MEFV_lineage_standardized.csv")
    crude = sd.pivot_table(index="patient", columns="region", values="crude_MEFV_pct")
    std = sd.pivot_table(index="patient", columns="region", values="lineage_standardized_MEFV_pct")
    V["REG_STD_TEXT"] = (f"mean AC − PA difference {sgn((crude.AC - crude.PA).mean(), 1)} percentage points crude and "
                         f"{sgn((std.AC - std.PA).mean(), 1)} after standardization, still lower in AC in "
                         + {3: "all three patients", 2: "two of three patients", 1: "one of three patients", 0: "no patient"}[int(((std.AC - std.PA) < 0).sum())])
    rg = r("R13c_region_glm_lineage_adjusted.csv")
    x = rg[rg.model.str.contains("lineage")].iloc[0]
    V["RR_AC_PA_LIN"] = f"{x.rate_ratio_AC_vs_PA:.2f} (95% CI {x.ci95_lo:.2f}–{x.ci95_hi:.2f})"
    words = {0: "none of the three patients", 1: "one of three patients", 2: "two of three patients", 3: "all three patients"}
    V["REG_PB_N"] = words[int(ra.loc["PYRIN_BACKBONE", "n_patients_AC_higher"])]
    V["REG_PB_DIFF"] = sgn(ra.loc["PYRIN_BACKBONE", "mean_AC_minus_PA"], 3)
    V["REG_NB_NLOW"] = words[3 - int(ra.loc["NLRP3_BACKBONE", "n_patients_AC_higher"])]
    V["REG_NB_DIFF"] = sgn(ra.loc["NLRP3_BACKBONE", "mean_AC_minus_PA"], 3)
    dl = "Specificity delta (cell-level z)"
    V["REG_DELTA_N"] = words[int(ra.loc[dl, "n_patients_AC_higher"])]
    V["REG_DELTA_DIFF"] = sgn(ra.loc[dl, "mean_AC_minus_PA"], 2)
    V["REG_CYT_DIFF"] = "mean difference " + sgn(ra.loc["Cytokine arm (CASP1/IL1B/IL18)", "mean_AC_minus_PA"], 3)
    V["REG_LYS_DIFF"] = "mean difference " + sgn(ra.loc["Lysis arm (GSDMD/GSDME/NINJ1)", "mean_AC_minus_PA"], 3)
    b = r("R22_bulk_paired_results.csv").set_index("module")
    for k, m in (("PB", "PYRIN_BACKBONE"), ("NB", "NLRP3_BACKBONE"), ("CYT", "EFFECTOR_CYTOKINE_ARM"),
                 ("LYS", "EFFECTOR_LYSIS_ARM"), ("MYE", "MYELOID_MARKER")):
        V[f"BULK_{k}"] = sgn(b.loc[m, "mean_paired_diff"], 3)
        V[f"BULK_{k}_DZ"] = f2(b.loc[m, "cohens_dz"])
        V[f"BULK_{k}_CI"] = ci(b.loc[m, "t_ci95_lo"], b.loc[m, "t_ci95_hi"])
    V["BULK_PB_P"] = f"{b.loc['PYRIN_BACKBONE', 'paired_t_p']:.3f}"
    V["BULK_PB_R"] = f2(b.loc["PYRIN_BACKBONE", "r_with_myeloid"])
    V["BULK_CYT_R"] = f2(b.loc["EFFECTOR_CYTOKINE_ARM", "r_with_myeloid"])
    V["BULK_PB_RES"] = sgn(b.loc["PYRIN_BACKBONE", "resid_mean_paired_diff"], 3)
    V["BULK_CYT_RES"] = sgn(b.loc["EFFECTOR_CYTOKINE_ARM", "resid_mean_paired_diff"], 3)
    V["BULK_PB_ATT"] = f"{b.loc['PYRIN_BACKBONE', 'attenuation_pct']:.0f}"
    V["BULK_PB_ATT_ROUND"] = f"{round(b.loc['PYRIN_BACKBONE', 'attenuation_pct'] / 5) * 5:.0f}"
    att = b.drop(index=["MYELOID_MARKER", "LEGACY_GENERIC_8GENE"]).attenuation_pct
    V["BULK_MIN_ATT"] = f"{np.floor(att.min()):.0f}"
    rec = r("R23_bulk_effect_vs_myeloid_enrichment.csv").set_index("module").drop(index=["LEGACY_GENERIC_8GENE"])
    from scipy import stats
    V["BULK_RHO"] = f2(stats.spearmanr(rec.myeloid_log2_enrichment_mean, rec.mean_paired_diff).statistic)
    s4b = r("R04b_pooled_vs_mean_of_samples.csv", index_col=0)
    V["T3_AC_POOLED"], V["T3_PA_POOLED"], V["T3_ALL_POOLED"] = (f2(s4b.loc[k, "pooled_pct"]) for k in ("AC", "PA", "ALL"))
    V["T3_AC_MEAN"], V["T3_PA_MEAN"] = f2(s4b.loc["AC", "mean_of_patient_pct"]), f2(s4b.loc["PA", "mean_of_patient_pct"])
    sc_ = r("R30_ode2_scenario_readouts.csv").set_index("scenario")
    bl, th, pr, ih = (sc_.loc[k] for k in ("baseline", "pyrin_sensitized_low_threshold", "high_inflammatory_priming",
                                          "inhibited_or_delayed_activation"))
    V["ODE_BASE_PA"], V["ODE_BASE_IL1B"], V["ODE_BASE_T50"] = f3(bl.peak_P_active), f3(bl.peak_IL1B_out), f2(bl.time_to_50pct_peak_IL1B)
    V["ODE_THR_T50"] = f2(th.time_to_50pct_peak_IL1B)
    V["ODE_THR_AUC"] = sgn(100 * (th.IL1B_AUC / bl.IL1B_AUC - 1), 1)
    V["ODE_PRIME_AUC"] = f"{100 * (pr.IL1B_AUC / bl.IL1B_AUC - 1):.0f}"
    V["ODE_PRIME_IL18"] = sgn(100 * (pr.IL18_AUC / bl.IL18_AUC - 1), 1)
    V["ODE_PRIME_RATIO"] = f2(pr.IL1B_to_IL18_AUC_ratio)
    V["ODE_INH_PA"], V["ODE_INH_T50"] = f3(ih.peak_P_active), f2(ih.time_to_50pct_peak_IL1B)
    mp = r("R31_ode2_matched_perturbations.csv").set_index("perturbation")
    V["ODE_M_PRIME"] = sgn(mp.loc["priming_plus20pct", "IL1B_AUC_pct_change"], 1)
    V["ODE_M_THR"] = sgn(mp.loc["threshold_minus20pct", "IL1B_AUC_pct_change"], 1)
    V["ODE_M_THR_T50"] = sgn(mp.loc["threshold_minus20pct", "time_to_50pct_peak_IL1B_pct_change"], 1)
    pc = r("R33_ode2_global_prcc.csv").set_index(["output", "parameter"])
    V["PRCC_RATIO_PRIME"] = sgn(pc.loc[("IL1B_to_IL18_AUC_ratio", "priming_input"), "PRCC"], 2)
    V["PRCC_T50_TEXT"] = "; ".join(f"{lab} {sgn(pc.loc[('time_to_50pct_peak_IL1B', p), 'PRCC'], 2)}" for p, lab in (
        ("signal2_input", "signal 2"), ("pyrin_activation_threshold", "threshold"),
        ("CASP1_activation_rate", "caspase-1 activation"), ("GSDMD_cleavage_rate", "GSDMD cleavage"),
        ("priming_input", "priming")))
    pt = r("R41_priority_table_revised.csv")
    top = pt.iloc[0]
    V["PRIO_1_TEXT"] = (f"raw positive score {top.raw_positive_score:.1f}, penalty {top.total_penalty:.1f}, final score "
                        f"{top.final_priority_score:.1f}; {top.confidence_grade} band")
    gs = pt.set_index("candidate_id").final_priority_score
    V["PRIO_REST_TEXT"] = (f"The myeloid-compartment pyrin-backbone state ({gs['L1-01']:.1f}) and the NLRP3 comparator "
                           f"({gs['L2-06']:.1f}) followed, and the classical-monocyte *MEFV*-high subset scored "
                           f"{gs['L1-02']:.1f}. The unstable-plaque bulk signal ({gs['L1-05']:.1f}) and the shared effector "
                           f"arms ({gs['L2-05']:.1f}) were penalized for myeloid composition, and the plaque-core *MEFV* signal "
                           f"ranked last ({gs['L1-04']:.1f}) after the corrected regional analysis.")
    return V


# ------------------------------------------------------------------------------- citations
def number_citations(text, refdb):
    order = []
    for grp in re.findall(r"\[(@[^\]]+)\]", text):
        for k in re.findall(r"@([\w]+)", grp):
            if k not in order:
                order.append(k)
    num = {k: i + 1 for i, k in enumerate(order)}

    def fmt(ns):
        ns = sorted(set(ns))
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(f"{ns[i]}–{ns[j]}" if j - i >= 2 else ",".join(str(x) for x in ns[i:j + 1]))
            i = j + 1
        return "[" + ",".join(out) + "]"
    text = re.sub(r"\[(@[^\]]+)\]", lambda m: fmt([num[k] for k in re.findall(r"@([\w]+)", m.group(1))]), text)
    missing = [k for k in order if k not in refdb]
    assert not missing, missing
    refs = [f"[{num[k]}] {refdb[k]}" for k in order]
    return text, refs


# ------------------------------------------------------------------------------- tables
def tables():
    comp = r("R05_compartment_summary_revised.csv", index_col=0)
    s04 = r("R04_per_sample_myeloid_MEFV.csv")
    s04b = r("R04b_pooled_vs_mean_of_samples.csv", index_col=0)
    sub = r("R12a_myeloid_subsets.csv")
    ra = r("R13a_myeloid_AC_vs_PA_corrected.csv")
    rg = r("R13c_region_glm_lineage_adjusted.csv")
    b = r("R22_bulk_paired_results.csv").set_index("module")
    pt = r("R41_priority_table_revised.csv")
    W = {}

    def t1(d):
        DB.add_table(d, ["Dataset", "Assay", "Analytic sample", "Comparison", "Role", "Principal limitation"],
                     [["GSE159677", "Single-cell RNA-seq (10x)", f"3 patients; 6 samples; {int(comp.n_cells.sum()):,} singlets",
                       "Calcified carotid AC vs patient-matched PA", "Compartment and myeloid-subset mapping; paired regional comparison",
                       "Three patients; marker-based labels; sparse *MEFV*; no genotype"],
                      ["GSE120521", "Bulk RNA-seq (FPKM)", "4 patients; 8 paired regions", "Stable vs unstable regions of the same plaque",
                       "Paired tissue-level triangulation; composition assessment", "Four pairs; FPKM only; no cell resolution"]],
                     [0.9, 1.0, 1.25, 1.2, 1.35, 1.3],
                     caption="**Table 1.** Public datasets and their analytic roles",
                     note="Sample labels for GSE159677 were verified by barcode matching against the depositors' per-sample molecule_info files; GSE120521 pairing and status were taken from GEO sample characteristics. AC, atherosclerotic core; PA, proximal adjacent.")
    W["t1"] = t1

    def t2(d):
        rows = []
        for c in CT:
            x = comp.loc[c]
            rows.append([CT_LAB[c], f"{int(x.n_cells):,}", f"{x.MEFV_pct:.2f} / {x.MEFV_pct_SoupX:.2f} / {x.MEFV_pct_DecontX:.2f}",
                         f"{x.PYRIN_BACKBONE_mean:.3f}", f"{x.NLRP3_BACKBONE_mean:.3f}", sgn(x.delta_cell_level_z, 3),
                         sgn(x.delta_compartment_mean_z, 2),
                         f"{sgn(x.delta_patient_mean, 2)} ({neg(x.delta_patient_min)} to {neg(x.delta_patient_max)})"])
        DB.add_table(d, ["Compartment", "Cells, n", "*MEFV* detected, % (observed / SoupX / DecontX)", "Mean PYRIN_BACKBONE",
                         "Mean NLRP3_BACKBONE", "Specificity delta, cell-level z", "Specificity delta, compartment-mean z",
                         "Specificity delta, patient mean (range)"], rows, [1.15, 0.6, 1.15, 0.8, 0.8, 0.85, 0.85, 1.0],
                     caption="**Table 2.** Single-cell compartment summary (scDblFinder singlets)",
                     note=("Cell-level delta: mean over the compartment's cells of z(PYRIN_BACKBONE) − z(NLRP3_BACKBONE), with z "
                           "computed across all singlets; the cell-count-weighted sum across compartments is zero. Compartment-mean "
                           "delta: z computed across the six compartment means (the standardization used in the original Figure 2B). "
                           "Patient-level values use the cell-level definition within each patient. Scores are relative "
                           "transcriptional summaries and do not indicate inflammasome activation."))
    W["t2"] = t2

    def t3(d):
        rows = []
        for _, x in s04.iterrows():
            rows.append([f"{x.patient}, {x.region} ({x.gsm})", f"{int(x.n_QC_cells):,}", f"{int(x.n_myeloid_QC):,}",
                         f"{int(x.MEFV_pos_myeloid_QC)} ({x.MEFV_pct_myeloid_QC:.2f})", f"{int(x.n_myeloid_singlets):,}",
                         f"{int(x.MEFV_pos_myeloid_singlets)} ({x.MEFV_pct_myeloid_singlets:.2f})",
                         f"{x.MEFV_pct_myeloid_singlets_SoupX:.2f} / {x.MEFV_pct_myeloid_singlets_DecontX:.2f}",
                         f"{x.median_nUMI_myeloid_singlets:,.0f}"])
        for lab, key in (("Pooled, all AC samples", "AC"), ("Pooled, all PA samples", "PA"), ("Pooled, all samples", "ALL")):
            sel = s04 if key == "ALL" else s04[s04.region == key]
            rows.append([lab, f"{int(sel.n_QC_cells.sum()):,}", f"{int(sel.n_myeloid_QC.sum()):,}",
                         f"{int(sel.MEFV_pos_myeloid_QC.sum())} ({100 * sel.MEFV_pos_myeloid_QC.sum() / sel.n_myeloid_QC.sum():.2f})",
                         f"{int(sel.n_myeloid_singlets.sum()):,}",
                         f"{int(sel.MEFV_pos_myeloid_singlets.sum())} ({s04b.loc[key, 'pooled_pct']:.2f})", "", ""])
        DB.add_table(d, ["Sample", "QC cells", "Myeloid cells (QC)", "*MEFV*+ myeloid, n (%)", "Myeloid singlets",
                         "*MEFV*+ myeloid singlets, n (%)", "SoupX / DecontX, %", "Median UMIs (myeloid singlets)"], rows,
                     [1.55, 0.6, 0.7, 0.9, 0.75, 1.0, 0.85, 0.8],
                     caption="**Table 3.** Myeloid *MEFV* detection in each of the six samples",
                     note=("Pooled rates are total *MEFV*-positive cells divided by total myeloid cells; the unweighted means of the "
                           f"three per-patient rates are {s04b.loc['AC', 'mean_of_patient_pct']:.2f}% (AC) and "
                           f"{s04b.loc['PA', 'mean_of_patient_pct']:.2f}% (PA). The two estimands differ because Patient 1 AC "
                           "contributes almost half of all myeloid cells. Regions use the verified sample labels."))
    W["t3"] = t3

    def t4(d):
        rows = []
        mobs = __import__("anndata").read_h5ad("data/processed/gse159677_myeloid.h5ad", backed="r").obs
        for subset, x in sub.groupby("subset"):
            n = x.n_cells.sum()
            g = mobs[mobs.subset == subset]
            pp = g.groupby("patient", observed=True)["MEFV_detected"].agg(["mean", "size"])
            rows.append((np.average(x.MEFV_pct, weights=x.n_cells), [subset, x.lineage.iloc[0], f"{int(n):,}",
                         f"{100 * (g.region == 'AC').mean():.0f}", f"{np.average(x.MEFV_pct, weights=x.n_cells):.1f}",
                         " / ".join(f"{100 * pp.loc[p, 'mean']:.1f}" if p in pp.index and pp.loc[p, 'size'] >= 20 else "–"
                                    for p in ("Patient 1", "Patient 2", "Patient 3")),
                         f"{g.PYRIN_BACKBONE_score.mean():.3f}", f"{g.NLRP3_BACKBONE_score.mean():.3f}"]))
        rows = [x for _, x in sorted(rows, key=lambda t: -t[0])]
        DB.add_table(d, ["Subset", "Lineage", "Cells, n", "From AC, %", "*MEFV* detected, %", "*MEFV* by patient, % (P1 / P2 / P3)",
                         "Mean PYRIN_BACKBONE", "Mean NLRP3_BACKBONE"], rows, [1.7, 0.85, 0.6, 0.6, 0.7, 1.15, 0.8, 0.8],
                     caption="**Table 4.** Myeloid sub-clusters after Harmony integration (singlets)",
                     note=("Sub-clusters were labelled by the highest standardized signature score; sub-clusters sharing a label are "
                           "pooled. Per-patient rates are shown when the patient contributed ≥20 cells. Mixed sub-clusters "
                           "carry T-cell or smooth-muscle transcripts and are reported for completeness."))
    W["t4"] = t4

    def t5(d):
        rows = []
        names = {"MEFV detection (%)": "*MEFV* detection, %", "Monocyte-lineage fraction (%)": "Monocyte-lineage cells, % of myeloid",
                 "TREM2+ lipid-associated macrophage fraction (%)": "TREM2^+^ lipid-associated macrophages, % of myeloid",
                 "PYRIN_BACKBONE": "PYRIN_BACKBONE", "NLRP3_BACKBONE": "NLRP3_BACKBONE",
                 "Specificity delta (cell-level z)": "Specificity delta (cell-level z)",
                 "Cytokine arm (CASP1/IL1B/IL18)": "Cytokine arm (*CASP1, IL1B, IL18*)",
                 "Lysis arm (GSDMD/GSDME/NINJ1)": "Lysis arm (*GSDMD, GSDME, NINJ1*)"}
        for _, x in ra.set_index("metric").loc[list(names)].reset_index().iterrows():
            dd = 1 if "%" in x.metric else 3
            rows.append([names[x.metric], f"{x.PA_mean_of_patients:.{dd}f}", f"{x.AC_mean_of_patients:.{dd}f}",
                         sgn(x.mean_AC_minus_PA, dd), f"{int(x.n_patients_AC_higher)}/3",
                         "; ".join(sgn(float(v), dd) for v in x.per_patient_AC_minus_PA.split(";"))])
        glm = "; ".join(f"{m.split('|')[1].strip()}: rate ratio {x.rate_ratio_AC_vs_PA:.2f} ({x.ci95_lo:.2f}–{x.ci95_hi:.2f})"
                        for m, x in rg.set_index("model").iterrows())
        DB.add_table(d, ["Metric (myeloid compartment)", "PA mean", "AC mean", "AC − PA", "AC higher", "Per-patient differences (P1; P2; P3)"],
                     rows, [2.3, 0.7, 0.7, 0.75, 0.65, 1.9],
                     caption="**Table 5.** Paired core (AC) versus adjacent (PA) comparison in the myeloid compartment, with verified sample labels",
                     note=("Means are unweighted means of the three patient-level values; no formal hypothesis test is emphasized "
                           f"(n = 3). Complementary log–log models of *MEFV* detection (AC vs PA): {glm}."))
    W["t5"] = t5

    def t6(d):
        names = [("PYRIN_BACKBONE", "PYRIN_BACKBONE"), ("PYRIN_FULL", "PYRIN_FULL"), ("NLRP3_BACKBONE", "NLRP3_BACKBONE"),
                 ("NLRP3_FULL", "NLRP3_FULL"), ("EFFECTOR_CYTOKINE_ARM", "Cytokine arm (*CASP1, IL1B, IL18*)"),
                 ("EFFECTOR_LYSIS_ARM", "Lysis arm (*GSDMD, GSDME, NINJ1*)"), ("NONCANONICAL_CASP4_5", "Caspase-4/5"),
                 ("MYELOID_MARKER", "Myeloid marker")]
        rows = []
        for k, lab in names:
            x = b.loc[k]
            res = "" if k == "MYELOID_MARKER" else f"{sgn(x.resid_mean_paired_diff, 3)} ({int(x.resid_n_pairs_positive)}/4)"
            attn = "" if k == "MYELOID_MARKER" else f"{x.attenuation_pct:.0f}"
            rr_ = "" if k == "MYELOID_MARKER" else f"{x.r_with_myeloid:.2f}"
            rows.append([lab, str(int(x.n_genes_present)), sgn(x.mean_paired_diff, 3), f"{int(x.n_pairs_positive)}/4",
                         f"{x.cohens_dz:.2f}", ci(x.t_ci95_lo, x.t_ci95_hi), f"{x.paired_t_p:.3f}", rr_, res, attn])
        DB.add_table(d, ["Module", "Genes, n", "Unstable − stable", "Pairs higher", "d~z~", "95% CI", "P", "r with myeloid score",
                         "Residualized difference (pairs higher)", "Attenuation, %"], rows,
                     [1.6, 0.45, 0.65, 0.5, 0.45, 0.95, 0.45, 0.6, 0.95, 0.6],
                     caption="**Table 6.** Paired bulk stable versus unstable plaque-region analyses (GSE120521, 4 plaques)",
                     note=("Module scores are means of gene-level z scores of log2(FPKM+1). Confidence intervals use the t "
                           "distribution with 3 degrees of freedom; P values are secondary descriptive measures. *GSDME* was mapped "
                           "from its former symbol *DFNA5*. Residualization is a linear adjustment for the myeloid-marker score across "
                           "the eight samples, not formal deconvolution; attenuation >100% indicates a reversed residual difference."))
    W["t6"] = t6

    def t7(d):
        rows = []
        for i, x in pt.reset_index(drop=True).iterrows():
            rows.append([str(i + 1), x.candidate_name, x.biological_axis, f"{x.raw_positive_score:.2f}", f"{x.total_penalty:.2f}",
                         f"{x.final_priority_score:.2f} ({x.confidence_grade})"])
        DB.add_table(d, ["Rank", "Candidate", "Biological axis", "Raw positive score", "Penalty", "Final score (grade)"], rows,
                     [0.45, 2.6, 1.6, 0.8, 0.7, 1.0],
                     caption="**Table 7.** Caveat-aware prioritization of candidate mechanisms and cell states (revised)",
                     note=("Scores are expert-defined research-prioritization outputs, not validated target or clinical rankings. "
                           "Weights and grade thresholds are unchanged from the original submission; every changed component is "
                           "listed in Supplementary Table S23."))
    W["t7"] = t7
    return W


LEGENDS = [
    "**Figure 1. Sensor-proximal programs, separable effector arms and analytic workflow.** (A) The pyrin regulatory backbone and the NLRP3 comparator backbone converge on ASC and caspase-1. Shared downstream genes (*PYCARD, CASP1, GSDMD, IL1B, IL18*) were excluded from the discriminative backbone scores because they cannot attribute activity to one sensor. Downstream of caspase-1, the cytokine arm (*CASP1, IL1B, IL18*) and the lysis arm (*GSDMD, GSDME, NINJ1*) are drawn separately because cytokine release and lytic cell death are separable; non-canonical caspase-4/5 signalling is reported separately. (B) Analytic sequence and the inferential role of each layer. The lower statement lists constraints that applied at every stage. Drawn programmatically from the analysis repository.",
    "**Figure 2. Single-cell localization of *MEFV* and the pyrin regulatory backbone, with sequencing-depth controls (GSE159677 singlets).** (A) Percentage of cells with detectable *MEFV* by compartment (bars: pooled singlets; open symbols: individual patients; vertical ticks: SoupX- and DecontX-corrected values). (B) Compartment means of cell-level z scores (standardized across all singlets) for PYRIN_BACKBONE, NLRP3_BACKBONE, the specificity delta and the two shared-effector arms; this is the definition used in Table 2. (C) Cell-level distributions of PYRIN_BACKBONE and NLRP3_BACKBONE scores; horizontal lines indicate medians. (D) *MEFV* detection within global quintiles of UMI count; the T/NK value in the top quintile rests on 149 cells. (E) Rate ratios per UMI from a complementary log–log binomial model with a log(UMI) offset and sample fixed effects (95% CIs). Module scores are relative transcriptional summaries and do not measure inflammasome activation.",
    "**Figure 3. *MEFV* detection across myeloid subsets.** (A) UMAP of myeloid-compartment singlets after Harmony integration across samples, coloured by lineage and labelled by subset. LAM, lipid-associated macrophage. (B) Percentage of cells with detectable *MEFV* by subset (bars: pooled; symbols: patients contributing ≥20 cells). (C) Detection (dot size) and scaled mean expression (colour) of *MEFV* and lineage markers by subset. (D) Percentage of *MEFV*-positive and *MEFV*-negative myeloid cells in which each marker was detected, with odds ratios.",
    "**Figure 4. Regional and bulk tissue triangulation.** (A) Myeloid compartment of adjacent (PA) and core (AC) tissue in the three patients, using verified sample labels: *MEFV* detection, monocyte-lineage fraction, PYRIN_BACKBONE and the specificity delta. (B) Paired stable and unstable plaque regions in GSE120521 (four plaques) for PYRIN_BACKBONE, NLRP3_BACKBONE and the cytokine and lysis arms. (C) Mean paired unstable − stable differences before and after residualization on the myeloid-marker score for every module. (D) Bulk paired differences plotted against each module's myeloid enrichment in the single-cell data. Bulk results are tissue-level associations and cannot distinguish cell abundance from cell-intrinsic expression.",
    "**Figure 5. Two-signal dimensionless model.** (A) Trajectories of active pyrin, extracellular IL-1β, extracellular IL-18 and the lysed fraction under baseline, lower pyrin activation threshold, higher priming (signal 1) and inhibited/delayed gating. (B) Percentage change from baseline after equal-size (20%) perturbations of the threshold, signal 2 and priming and an inhibition factor of 0.2. (C) Partial rank correlation coefficients from 1,000 Latin hypercube parameter sets (all positive parameters sampled log-uniformly from 0.5- to 2-fold of baseline). Time and states are dimensionless; the model was not calibrated to patients, drugs or *MEFV* variants.",
    "**Figure 6. Caveat-aware hypothesis priority score (revised).** Blue bars show weighted positive evidence, red bars weighted penalties (left of zero) and diamonds the final scores. Dashed lines mark the expert-defined moderate (≥35) and high (≥55) bands. Only the *MEFV*/pyrin threshold axis reached the moderate band; no candidate reached the high band. This is a research-prioritization aid, not a validated target, drug or clinical ranking.",
]


def main():
    template, original, outdir = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    V = values()
    src = open(template, encoding="utf-8").read()
    miss = sorted(set(re.findall(r"\{\{(\w+)\}\}", src)) - set(V))
    assert not miss, f"unfilled placeholders: {miss}"
    for k, v in V.items():
        src = src.replace("{{" + k + "}}", v)
    refdb = json.load(open("tools/refs.json", encoding="utf-8"))
    src, refs = number_citations(src, refdb)
    json.dump(V, open(os.path.join(outdir, "values_used.json"), "w"), indent=1, ensure_ascii=False)
    body = src.replace("[[REFERENCES]]", "\n\n".join(refs))
    body += "\n\n# Tables\n\n" + "\n\n".join(f"[[TABLE:t{i}]]" for i in range(1, 8))
    body += "\n\n# Figure legends\n\n" + "\n\n".join(LEGENDS)
    open(os.path.join(outdir, "manuscript_revised_filled.md"), "w", encoding="utf-8").write(body)
    blocks = DB.parse_markup(body)
    W = tables()

    # clean copy
    d = DB.new_document()
    for blk in blocks:
        p = DB.add_block(d, blk, W)
        if p is not None and blk is blocks[0]:
            p.alignment = 1
    d.save(os.path.join(outdir, "Manuscript_revised_clean.docx"))

    # tracked copy: align body blocks (before Tables) and legend blocks with the original
    old_all = DB.blocks_from_docx(original)
    oi = next(i for i, b in enumerate(old_all) if DB._plain(b["runs"]).strip() in ("Tables", "**Tables**"))
    li = next(i for i, b in enumerate(old_all) if DB._plain(b["runs"]).strip().startswith("Figures and Legends"))
    old_body, old_leg = old_all[:oi], old_all[li + 1:]
    ni = next(i for i, b in enumerate(blocks) if b["kind"] == "h1" and DB._plain(b["runs"]).strip() == "Tables")
    nl = next(i for i, b in enumerate(blocks) if b["kind"] == "h1" and DB._plain(b["runs"]).strip() == "Figure legends")
    new_body, new_tabs, new_leg = blocks[:ni], blocks[ni + 1:nl], blocks[nl + 1:]
    td = DB.new_document()
    sty = {"h1": "Heading 1", "h2": "Heading 2", "h3": "Heading 3", "p": "Normal"}

    def emit(old, new):
        for i, j in DB.align_blocks(old, new):
            if i is not None and j is not None:
                DB.tracked_paragraph(td, sty.get(new[j]["kind"], "Normal"), old[i]["runs"], new[j]["runs"])
            elif i is not None:
                DB.deleted_paragraph(td, sty.get(old[i]["kind"], "Normal"), old[i]["runs"])
            else:
                DB.inserted_paragraph(td, sty.get(new[j]["kind"], "Normal"), new[j]["runs"])
    emit(old_body, new_body)
    DB.tracked_paragraph(td, "Heading 1", [("Tables", "")], [("Tables", "")])
    import docx as _docx
    od = _docx.Document(original)
    for t in od.tables:
        for row in t.rows:
            txt = " | ".join(dict.fromkeys(c.text.strip() for c in row.cells))
            DB.deleted_paragraph(td, "Normal", [(txt, "")])
    n_before = len(td.element.body)
    for blk in new_tabs:
        DB.add_block(td, blk, W)
    # mark every run and row of the newly added tables/captions as inserted
    from docx.oxml.ns import qn
    body_el = td.element.body
    for el in list(body_el)[n_before - 1:]:
        for rr in el.iter(qn("w:r")):
            parent = rr.getparent()
            if parent.tag == qn("w:ins"):
                continue
            w = DB.OxmlElement("w:ins")
            for a, v in (("w:id", DB._next_id()), ("w:author", DB.AUTHOR), ("w:date", DB.DATE)):
                w.set(qn(a), v)
            parent.replace(rr, w)
            w.append(rr)
        for tr in el.iter(qn("w:tr")):
            trpr = tr.find(qn("w:trPr"))
            if trpr is None:
                trpr = DB.OxmlElement("w:trPr")
                tr.insert(0 if tr.find(qn("w:tblPrEx")) is None else 1, trpr)
            m = DB.OxmlElement("w:ins")
            for a, v in (("w:id", DB._next_id()), ("w:author", DB.AUTHOR), ("w:date", DB.DATE)):
                m.set(qn(a), v)
            trpr.append(m)
    DB.tracked_paragraph(td, "Heading 1", [("Figures and Legends", "")], [("Figure legends", "")])
    emit(old_leg, new_leg)
    td.save(os.path.join(outdir, "Manuscript_revised_tracked_changes.docx"))
    words = len(re.sub(r"\[[0-9,–]+\]", "", " ".join(DB._plain(b["runs"]) for b in new_body if b["kind"] == "p")).split())
    abstract = " ".join(DB._plain(b["runs"]) for b in blocks[blocks.index(next(b for b in blocks if DB._plain(b.get("runs", [])) == "Abstract")) + 1:][:5])
    print(f"placeholders filled: {len(V)}; references: {len(refs)}; body words (paragraphs incl. front matter): {words}; "
          f"abstract words: {len(abstract.split())}")


if __name__ == "__main__":
    main()
