"""Per-cell module scoring for pyrin/NLRP3/generic programs. Config-driven; no hard-coded genes."""
import numpy as np
import scanpy as sc
from scipy.stats import rankdata

def detected_genes(adata, genes):
    return [g for g in genes if g in adata.var_names]

def score_modules(adata, scoring_sets, ctrl_size=50, n_bins=25, seed=0):
    """Add <SET>_score (scanpy score_genes) and <SET>_rank (rank-based) to adata.obs."""
    np.random.seed(seed)
    for name, genes in scoring_sets.items():
        g = detected_genes(adata, genes)
        sc.tl.score_genes(adata, g, score_name=f"{name}_score",
                          ctrl_size=ctrl_size, n_bins=n_bins, random_state=seed, use_raw=False)
        adata.obs[f"{name}_rank"] = _rank_score(adata, g)
    return adata

def _rank_score(adata, genes):
    genes = [g for g in genes if g in adata.var_names]
    X = adata[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    ranks = np.apply_along_axis(lambda c: rankdata(c) / len(c), 0, X)
    return ranks.mean(1)

def zscore(series):
    return (series - series.mean()) / series.std()

def specificity_delta(obs, group="cell_type", pyrin="PYRIN_BACKBONE", nlrp3="NLRP3_BACKBONE"):
    """z(PYRIN_BACKBONE) - z(NLRP3_BACKBONE) aggregated by group. Not causal activity."""
    o = obs.copy()
    o[f"{pyrin}_z"] = zscore(o[f"{pyrin}_score"])
    o[f"{nlrp3}_z"] = zscore(o[f"{nlrp3}_score"])
    g = o.groupby(group)
    return (g[f"{pyrin}_z"].mean() - g[f"{nlrp3}_z"].mean()).sort_values(ascending=False)
