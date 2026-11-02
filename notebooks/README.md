# V4 notebooks

Notebook filenames follow one continuous numeric sequence. Run the main workflow in order; Step 8 is an optional audit for the V3-to-V4 migration.

| Step | Notebook | Purpose |
|---:|---|---|
| 00 | `00_configure_run` | Choose pathogen groups, hosts, run mode, and PubMed sample size; write the shared run config |
| 01 | `01_review_taxonomy_and_hosts` | Inspect pinned pathogen and host taxonomy and check references |
| 02 | `02_preview_prompts` | Preview the versioned pathogen and pair-feature prompts without provider calls |
| 03 | `03_acquire_pair_pubmed_records` | Search pathogen–host pairs and acquire PMID results and screening abstracts |
| 04 | `04_screen_pathogen_host_pairs` | Screen sampled abstracts and decide whether each pathogen causes disease in the host |
| 05 | `05_estimate_pathogen_publications` | Estimate pathogen-supporting publications and calculate the pathogenicity index |
| 06 | `06_handoff_disease_evidence` | Validate screening coverage and materialize pair decisions and evidence spans |
| 07 | `07_resolve_feature_profiles` | Resolve the applicable feature set for each pathogen–host context |
| 08 | `08_audit_migration` | Audit the V3 snapshot, conversions, and staging links; optional for native V4 runs |
| 09 | `09_validate_release` | Validate pair tables, provenance, and release readiness |
| 10 | `10_inspect_embeddings` | Inspect coordinate coverage and export an embeddings plot |
| 11 | `11_resolve_wikipedia_references` | Resolve and checkpoint Wikipedia article links |
| 12 | `12_resolve_gbif_references` | Resolve and checkpoint GBIF species references |

## Run sequence

Start with `00_configure_run.ipynb`. Select pathogen groups, hosts, smoke-test or full mode, and the per-pair abstract sample cap (default 200, including the 30 lowest PMID IDs). This writes `configs/pilot_v1.yaml`. Available hosts are Human, Cattle, Sheep, Chicken, Pig, Maize, Rice, Wheat, and Soybean. Start Jupyter from the V4 repository.

Continue through Steps 1–7 in order. The taxonomy review checks shared, pinned inputs. Prompt previews make no model calls. Pair acquisition retains the complete PMID set; abstracts are fetched for the stable screening sample. Screening accepts a pair only when an abstract supplies an exact supporting quote. The publication estimate uses positive sampled records divided by the selected sample, scaled by the pair-search record count. Count-only PubMed searches run for pathogens accepted in at least one selected host.

The disease-evidence handoff consumes completed screening decisions and does not repeat model calls. Feature profiles are resolved separately from evidence scoring. Run the migration audit only when reviewing the V3-to-V4 conversion; then validate release tables and inspect embedding and reference outputs as needed.

The shared helper in `_shared.py` writes logs under `data/v4/logs/<run-id>/`. PubMed notebooks resolve the NCBI key from `NCBI_API_KEY`, `EUTILS_API_KEY`, or the macOS Keychain item used by local V2/V3 workflows. They use `NCBI_EMAIL`, then the configured Git email, then prompt for a contact address. The Wikipedia and GBIF steps checkpoint completed lookups and reuse them on rerun.
