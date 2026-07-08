# Data Licenses & Attribution

The **code** in this repository is MIT-licensed (see LICENSE). The **data** are public third-party resources. **No raw public data are redistributed in this repository** — only derived tables, figures, and metadata. Users must download raw data from the official sources below and comply with each source's terms.

## Datasets (NCBI GEO)
GEO data are subject to NCBI's data usage policies (https://www.ncbi.nlm.nih.gov/home/about/policies/).

| Accession | Type | Use here | Download |
|-----------|------|----------|----------|
| **GSE159677** | scRNA-seq, human carotid plaque | Primary single-cell | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677 |
| **GSE120521** | Bulk RNA-seq (FPKM), stable/unstable plaque | Bulk triangulation | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120521 |
| **GSE253902** | CITE-seq, human carotid | Listed in manifest; NOT analyzed in MVP | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253902 |

The integrated single-cell atlas (Nat Commun 2025, PMID 40931012) is cited as **literature precedent only** and was **not used as a data source**.

## Literature & annotation resources
- **PubMed / PMC** — literature metadata (NCBI; https://www.ncbi.nlm.nih.gov/home/about/policies/). Abstracts/metadata used for the evidence graph; no full-text redistribution.
- **ChEMBL** — EMBL-EBI, CC BY-SA 3.0 (https://www.ebi.ac.uk/chembl/).
- **Open Targets Platform** — target tractability; CC0/CC BY per Open Targets terms (https://platform.opentargets.org/).
- **ClinicalTrials.gov** — trial records (public domain; https://clinicaltrials.gov/).
- **InterPro** (protein-annotation) — EMBL-EBI, CC0 (https://www.ebi.ac.uk/interpro/).
- **MyGene.info** (genes-ontologies) — gene ID mapping (https://mygene.info/).

## Redistribution statement
This repository redistributes **no raw public data**. Derived products (module scores, effect sizes, priority tables, figures) are original computed outputs released under MIT alongside the code. To reproduce, download the raw datasets from the official links above; input checksums are recorded for verification (GSE159677 tar `8cd0b7f8…`, GSE120521 FPKM `27e608a5…`).
