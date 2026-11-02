# V4 notebooks

The notebooks use a phase-based workflow similar in spirit to the earlier notebooks, with V4 boundaries between taxonomy, pair selection, literature acquisition, disease decisions, feature profiles, release checks, and embeddings. Each phase includes input review, a preflight or validation step, an explicit run switch where calls or promotion are consequential, summary tables, and durable logs.

| Notebook | Phase | Purpose |
|---|---:|---|
| `00_workflow` | 0 | Choose pathogen groups, hosts, run mode, and the per-pathogen PubMed sample size; author the shared config |
| `01_taxonomy_hosts` | 1 | Inspect pinned host/pathogen taxonomy and check references |
| `02_runs_and_pubmed` | 2 | Acquire complete pair PubMed matches and a resumable pathogen-only count/sample corpus |
| `03_disease_gate` | 3 | Estimate task volume, validate evidence spans, run/resume a provider-backed ledger |
| `04_feature_profiles` | 4 | Resolve profiles by pathogen group and host context; save profile hashes |
| `05_migration_audit` | M | Verify the snapshot, inspect conversion reports, check staging links |
| `06_release_validation` | 5 | Validate pair tables, summarize provenance, review release readiness |
| `07_embeddings_explorer` | 6 | Check coordinate coverage, plot one group, and export an SVG |
| `08_wikipedia_references` | 7 | Resolve and checkpoint Wikipedia article links for selected pathogen records |
| `09_gbif_references` | 8 | Resolve and checkpoint exact GBIF species references for selected pathogen records |

## Run logs

The shared helper in `_shared.py` writes console and file logs. File logs are placed under `data/v4/logs/<run-id>/`; run configs, resolved profiles, and plot exports are stored under the ignored `data/v4/` working area.

Start with `00_workflow.ipynb`: its first controls select pathogen groups, common host names, and smoke-test or full mode; it authors and writes `configs/pilot_v1.yaml` from those selections and versions. The later phase notebooks read that shared configuration. Start Jupyter from the V4 repository. Notebook `00_workflow` sets the host-independent PubMed sample cap (default 200 per pathogen, including 30 oldest PMID IDs), oldest-ID count, and random seed. Notebook `02_runs_and_pubmed` retains every PMID and title from pair searches. After pathogen screening, the host-independent literature count and sample are restricted to accepted pathogen IDs; the configured stable SHA-256 sample extends from 200 to 400 without changing the first 200 or repeating the search. Set `NCBI_EMAIL` before running; `NCBI_API_KEY` is optional and increases the allowed request rate. The Wikipedia and GBIF reference notebooks run after a resolved pathogen selection is available. They checkpoint completed lookups and reuse them on rerun. The migration checksum pass defaults to off because it re-reads the full snapshot; enable it when auditing the source copy.
