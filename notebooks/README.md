# V4 notebooks

The notebooks use a phase-based workflow similar in spirit to the earlier notebooks, with V4 boundaries between taxonomy, pair selection, literature acquisition, disease decisions, feature profiles, release checks, and embeddings. Each phase includes input review, a preflight or validation step, an explicit run switch where calls or promotion are consequential, summary tables, and durable logs.

| Notebook | Phase | Purpose |
|---|---:|---|
| `00_workflow` | 0 | Choose pathogen groups, hosts, run mode, and the per-pair PubMed abstract-screen sample size; author the shared config |
| `01_taxonomy_hosts` | 1 | Inspect pinned host/pathogen taxonomy and check references |
| `01_prompt_preview` | 1a | Print versioned pathogen yes/no and host-specific feature prompts without provider calls |
| `02_runs_and_pubmed` | 2 | Acquire complete pathogen–host pair PubMed matches and article titles |
| `02b_pathogen_screening` | 2b | Screen sampled abstracts for pathogen status in each selected host |
| `02c_pathogenic_interest_index` | 2c | Count total PubMed publications for accepted pathogens and calculate PII |
| `03_disease_gate` | 3 | Estimate task volume, validate evidence spans, run/resume a provider-backed ledger |
| `04_feature_profiles` | 4 | Resolve profiles by pathogen group and host context; save profile hashes |
| `05_migration_audit` | M | Verify the snapshot, inspect conversion reports, check staging links |
| `06_release_validation` | 5 | Validate pair tables, summarize provenance, review release readiness |
| `07_embeddings_explorer` | 6 | Check coordinate coverage, plot one group, and export an SVG |
| `08_wikipedia_references` | 7 | Resolve and checkpoint Wikipedia article links for selected pathogen records |
| `09_gbif_references` | 8 | Resolve and checkpoint exact GBIF species references for selected pathogen records |

## Run logs

The shared helper in `_shared.py` writes console and file logs. File logs are placed under `data/v4/logs/<run-id>/`; run configs, resolved profiles, and plot exports are stored under the ignored `data/v4/` working area.

Start with `00_workflow.ipynb`: its first controls select pathogen groups, common host names, and smoke-test or full mode; it authors and writes `configs/pilot_v1.yaml` from those selections and versions. The later phase notebooks read that shared configuration. Start Jupyter from the V4 repository. Notebook `00_workflow` sets the per-pair abstract-screen sample cap (default 200, including the 30 lowest PMID IDs), oldest-ID count, and random seed. Notebook `01_prompt_preview` prints the pathogen yes/no prompt and the full feature list resolved for a pathogen–host pair, plus one pair-feature prompt. Notebook `02_runs_and_pubmed` retains every PMID and title from pair searches. Notebook `02b_pathogen_screening` fetches abstracts only for the stable per-pair sample and accepts a pathogen–host pair when at least one abstract supplies a supporting quote. Increasing the cap preserves all previous sample IDs and fetches only newly added abstracts. Notebook `02c_pathogenic_interest_index` then runs count-only, organism-name PubMed searches for pathogens accepted in at least one host; this total is not extrapolated from the screen sample. Set `NCBI_EMAIL` before running; `NCBI_API_KEY` is optional and increases the allowed request rate. The Wikipedia and GBIF reference notebooks run after a resolved pathogen selection is available. They checkpoint completed lookups and reuse them on rerun. The migration checksum pass defaults to off because it re-reads the full snapshot; enable it when auditing the source copy.
