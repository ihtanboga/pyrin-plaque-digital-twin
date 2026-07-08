"""Bulk plaque module scoring on log-FPKM. Transparent z-score aggregation; config-driven."""
import numpy as np
import pandas as pd

def load_fpkm(path, gene_col="name"):
    """Load GSE120521 FPKM xlsx -> (expr DataFrame genes x samples, gene index)."""
    df = pd.read_excel(path)
    sample_cols = [c for c in df.columns if "FPKM" in c]
    expr = df[[gene_col] + sample_cols].copy()
    expr = expr.dropna(subset=[gene_col])
    # collapse duplicate gene symbols by max FPKM
    expr = expr.groupby(gene_col)[sample_cols].max()
    return expr, sample_cols

def log_transform(expr):
    """log2(FPKM+1)."""
    return np.log2(expr + 1.0)

def module_score(logexpr, genes):
    """Transparent module score: mean of per-gene z-scores (across samples) for detected genes.

    z per gene across samples, then average across the module's detected genes -> one score per sample.
    Genes absent from the matrix are dropped (reported separately).
    """
    present = [g for g in genes if g in logexpr.index]
    if not present:
        return pd.Series(np.nan, index=logexpr.columns), present
    sub = logexpr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present

def score_all(logexpr, scoring_sets):
    """Return (scores_df samples x sets, detected dict)."""
    scores, detected = {}, {}
    for name, genes in scoring_sets.items():
        s, present = module_score(logexpr, genes)
        scores[name] = s
        detected[name] = present
    return pd.DataFrame(scores), detected

def cohens_d_paired(diffs):
    """Paired Cohen's d = mean(diff)/sd(diff)."""
    diffs = np.asarray(diffs, float)
    sd = diffs.std(ddof=1)
    return float(diffs.mean() / sd) if sd > 0 else np.nan
