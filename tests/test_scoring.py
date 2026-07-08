import sys; sys.path.insert(0, "src")
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from pyrinplaque import scoring as S

def _toy():
    rng = np.random.default_rng(0)
    X = csr_matrix(rng.poisson(0.5, size=(200, 8)).astype("float32"))
    a = ad.AnnData(X)
    a.var_names = [f"G{i}" for i in range(8)]
    return a

def test_rank_score_range():
    a = _toy()
    r = S._rank_score(a, ["G0","G1","G2"])
    assert r.shape[0] == a.n_obs
    assert (r >= 0).all() and (r <= 1).all()

def test_zscore_zero_mean():
    import pandas as pd
    z = S.zscore(pd.Series([1.0, 2, 3, 4, 5]))
    assert abs(z.mean()) < 1e-9

def test_detected_genes_filters_missing():
    a = _toy()
    assert S.detected_genes(a, ["G0","NOPE"]) == ["G0"]
