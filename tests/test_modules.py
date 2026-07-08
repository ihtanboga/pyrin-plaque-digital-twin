import sys; sys.path.insert(0, "src")
from pyrinplaque import modules as M

SHARED = {"PYCARD", "CASP1", "GSDMD", "IL1B", "IL18"}

def test_scoring_sets_present():
    s = M.scoring_sets("config/gene_modules.yaml")
    for k in ["PYRIN_FULL","PYRIN_BACKBONE","NLRP3_FULL","NLRP3_BACKBONE","GENERIC_PYROPTOSIS"]:
        assert k in s and len(s[k]) > 0

def test_backbones_exclude_shared_downstream():
    s = M.scoring_sets("config/gene_modules.yaml")
    assert not (set(s["PYRIN_BACKBONE"]) & SHARED), "pyrin backbone must exclude shared downstream"
    assert not (set(s["NLRP3_BACKBONE"]) & SHARED), "nlrp3 backbone must exclude shared downstream"

def test_mefv_in_pyrin_backbone():
    s = M.scoring_sets("config/gene_modules.yaml")
    assert "MEFV" in s["PYRIN_BACKBONE"]
    assert "NLRP3" in s["NLRP3_BACKBONE"]

def test_pyrin_backbone_is_pyrin_specific_only():
    cfg = M.load_modules("config/gene_modules.yaml")
    bb = M.genes_by_role(cfg, "PYRIN_MODULE", "pyrin_specific")
    assert set(bb) == set(M.scoring_sets("config/gene_modules.yaml")["PYRIN_BACKBONE"])
