# V4 notebooks

Start with `A_generate_feature_questions` to prepare the versioned feature-question profiles, then run notebooks 00 through 12 in numeric order. Notebook 00 records the exact staging snapshot; downstream notebooks reuse it. Notebook 10 validates and builds the additive Explorer data payload and profile shards. Notebooks 11 and 12 add optional reference-link metadata.

| Step | Notebook | Purpose |
|---:|---|---|
| A | `A_generate_feature_questions` | Generate and maintain the preliminary-active pathogenicity feature profiles across configured host × pathogen-class combinations; export question caches before evidence runs |
| 00 | `00_configure_run` | Choose pathogen groups, hosts, run mode, and PubMed sample size |
| 01 | `01_review_taxonomy_and_hosts` | Review the pinned pathogen and host taxonomy |
| 02 | `02_preview_prompts` | Preview versioned Jev questions without provider calls |
| 03 | `03_acquire_pair_pubmed_records` | Search pathogen–host pairs and acquire sampled abstracts |
| 04 | `04_screen_pathogen_host_pairs` | Screen sampled abstracts and select supported pathogen–host pairs |
| 05 | `05_estimate_pathogen_publications` | Estimate pair-specific publication volume and pathogenicity fraction |
| 06 | `06_handoff_disease_evidence` | Materialize decisions linked to screened abstracts |
| 07 | `07_resolve_feature_profiles` | Resolve applicable feature sets for pathogen–host contexts |
| 08 | `08_score_feature_profiles` | Score every applicable question against accepted pathogen–host abstracts with Jev |
| 09 | `09_validate_release` | Validate pair tables, provenance, and release readiness |
| 10 | `10_inspect_embeddings` | Build the additive Explorer database, fit eligible host × pathogen-group IMSTE maps, and inspect a configured cohort |
| 11 | `11_resolve_wikipedia_references` | Resolve and checkpoint Wikipedia links |
| 12 | `12_resolve_gbif_references` | Resolve and checkpoint GBIF references |

Start Jupyter from the V4 repository. The screening cutoff maps Jev probabilities of yes at 0.7 or above to yes. Each decision remains linked to its PMID and screened-abstract hash. The publication estimate uses the positive sample rate and pair-search total; organism-wide PubMed counts are queried only for pathogens accepted in at least one selected host.

The shared helper in `_shared.py` writes logs under `data/v4/logs/<run-id>/`. PubMed notebooks resolve the NCBI key from `NCBI_API_KEY`, `EUTILS_API_KEY`, or the macOS Keychain and use `NCBI_EMAIL`, then the configured Git email, for contact details. Wikipedia and GBIF lookups checkpoint completed rows and reuse them on rerun.
