"""Provenance capture for reproducible single-cell scoring runs."""
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def run_provenance(inputs, params, outputs, extra=None):
    """Build a provenance dict: input checksums, params, package versions, timestamp."""
    try:
        import scanpy, anndata, numpy, scipy
        versions = {"scanpy": scanpy.__version__, "anndata": anndata.__version__,
                    "numpy": numpy.__version__, "scipy": scipy.__version__}
    except Exception:
        versions = {}
    prov = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "packages": versions,
        "inputs": {p: file_sha256(p) for p in inputs if p},
        "params": params,
        "outputs": outputs,
    }
    if extra:
        prov.update(extra)
    return prov

def write_provenance(prov, path):
    with open(path, "w") as f:
        json.dump(prov, f, indent=2)
    return path
