# Final Package Manifest

**Package:** `pyrin_plaque_digital_twin_submission.zip`
**Size:** ~3.8 MB (exact byte size in `PACKAGE_CHECKSUM.txt`) · **Entries:** 144
**sha256:** see `PACKAGE_CHECKSUM.txt` (external sidecar)

*(A manifest bundled inside the zip cannot record its own container's hash. The authoritative package sha256 is in `PACKAGE_CHECKSUM.txt`, saved alongside the zip — outside the bundle. Verify against that file, or against the saved zip artifact's own stored checksum.)*

## Included (by category, computed from actual zip listing)
| Category | Count |
|----------|-------|
| app | 2 |
| audit | 17 |
| code | 17 |
| config | 5 |
| evidence | 7 |
| figure | 15 |
| metadata | 0 |
| report | 43 |
| table | 23 |
| target_card | 6 |
| top-level | 9 |

## Top-level files
README.md · LICENSE · DATA_LICENSE.md · CITATION.cff · ARTIFACT_INDEX.md · Makefile · environment.yml · requirements-freeze.txt · .gitignore

## Included trees
- `config/` — frozen gene modules, ODE params, scenarios, scoring weights
- `manifests/` · `evidence/` — data manifest, evidence graph, claim cards, novelty audit
- `results/figures/` (5 main + supporting) · `results/tables/` (small derived CSVs only) · `results/target_cards/` (6)
- `reports/` + `reports/figure_cards/` — daily reviews, integrated results, demo, provenance, submission summaries
- `src/pyrinplaque/` · `scripts/` · `tests/` — code + 13 passing tests + smoke test
- `audit/` — consistency, provenance, reviewer flags, inventories
- `app/` — optional static Streamlit demo

## Deliberately EXCLUDED
| Item | Reason |
|------|--------|
| `data/raw/` | Raw public GEO data — not redistributed (download per DATA_LICENSE.md) |
| `data/processed/*.h5ad` | ~1.9 GB scored object — regenerate via `make all` |
| `results/tables/*.parquet` | Large per-cell table (cell_state_scores.parquet) — CSV summaries included instead |
| `__pycache__/`, `.pytest_cache/`, `.git/`, cache dirs | Build/cache artifacts |
| `data/metadata/` | Not in the submission include set (sample metadata; available in the working repo) |

**Verified:** zip listing contains no `data/raw`, no `.h5ad`, no `__pycache__`, no `.git`, no large parquet (checked by filename).

## Reproduce from package
```bash
conda env create -f environment.yml      # or pip install -r requirements-freeze.txt
make smoke    # fast integrity (PASS)
make test     # 13 passed
make all      # full path — requires public GEO downloads (see DATA_LICENSE.md)
```

## Checksum verification
```bash
sha256sum -c PACKAGE_CHECKSUM.txt   # or: shasum -a 256 -c PACKAGE_CHECKSUM.txt
```
The authoritative package sha256 is in `PACKAGE_CHECKSUM.txt`, saved alongside the zip outside the bundle.
