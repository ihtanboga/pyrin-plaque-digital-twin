# How Claude Was Used — PYRIN-PLAQUE Digital Twin

**How did you use Claude? Which products (Claude Science, Claude Code, etc.)? Where did they matter most in your workflow?**

Claude Science was the core scientific workbench. I used it to coordinate literature discovery, public dataset verification, single-cell and bulk omics analysis, figure generation, provenance capture, ODE modeling, and adversarial reviewer-agent audits. Its MCP connectors were essential: GEO/omics-archives verified GSE159677, GSE120521, and GSE253902; PubMed supported the 57-claim evidence graph; MyGene validated all module genes; ChEMBL, Open Targets, ClinicalTrials.gov, and InterPro supported druggability, clinical-alignment, and protein-domain context.

Claude Science mattered most in three places. First, it helped prevent overclaiming by separating pyrin-specific biology from generic NLRP3/shared inflammasome biology. Second, reviewer agents repeatedly found issues that changed the project: shared downstream circularity, unverified atlas usage, provisional cell labels, myeloid confounding, mixed AC-vs-PA signal, and the ODE result that generic priming dominated pyrin threshold effects. Third, every figure and number was linked to code, environment, and provenance.

Claude Code was used for final repository preparation: packaging the Claude Science outputs into a clean open-source GitHub repo, checking figure links, excluding raw data, validating the manifest, and preparing the submission-ready structure.

---

See also: [reports/claude_science_usage_statement.md](reports/claude_science_usage_statement.md) (detailed connector-by-connector breakdown) and [reports/reviewer_audit.md](reports/reviewer_audit.md) (full reviewer-agent gate history).
