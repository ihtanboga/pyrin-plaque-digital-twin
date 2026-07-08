"""Load frozen gene modules and derive scoring sets. No hard-coded gene lists."""
import yaml

SHARED_DOWNSTREAM = {"PYCARD", "CASP1", "GSDMD", "IL1B", "IL18"}

def load_modules(path="config/gene_modules.yaml"):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return {k: v for k, v in cfg.items() if k != "_meta"}

def module_genes(cfg, module):
    """All gene symbols in a module."""
    return list(cfg[module]["genes"].keys())

def genes_by_role(cfg, module, role):
    """Symbols in `module` whose axis_role == role."""
    return [g for g, gi in cfg[module]["genes"].items() if gi.get("axis_role") == role]

def scoring_sets(path="config/gene_modules.yaml"):
    """Derive the five scoring sets required for Day-3 analysis, straight from config.

    PYRIN_BACKBONE / NLRP3_BACKBONE exclude shared downstream genes so the main
    pyrin-vs-NLRP3 comparison is not driven by shared effectors (Day-1 caveat 1).
    """
    cfg = load_modules(path)
    pyrin_full = module_genes(cfg, "PYRIN_MODULE")
    nlrp3_full = module_genes(cfg, "NLRP3_COMPARISON_MODULE")
    generic = module_genes(cfg, "GENERIC_PYROPTOSIS_MODULE")
    pyrin_backbone = genes_by_role(cfg, "PYRIN_MODULE", "pyrin_specific")
    nlrp3_backbone = [g for g in nlrp3_full if g not in SHARED_DOWNSTREAM]
    return {
        "PYRIN_FULL": pyrin_full,
        "PYRIN_BACKBONE": pyrin_backbone,
        "NLRP3_FULL": nlrp3_full,
        "NLRP3_BACKBONE": nlrp3_backbone,
        "GENERIC_PYROPTOSIS": generic,
    }
