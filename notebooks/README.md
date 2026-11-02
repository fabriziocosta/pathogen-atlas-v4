# V4 notebooks

The notebooks follow the modular run lifecycle. They call V4 package contracts and keep taxonomy, run resolution, PubMed query construction, disease attribution, feature profiles, migration auditing, and embedding inspection separate.

| Notebook | Purpose |
|---|---|
| `00_workflow` | Workflow map and repository context |
| `01_taxonomy_hosts` | Inspect pinned host and pathogen taxonomy |
| `02_runs_and_pubmed` | Resolve pairs and preview versioned PubMed queries |
| `03_disease_gate` | Validate pair-gate records and evidence-span rules |
| `04_feature_profiles` | Resolve ordered feature modules |
| `05_migration_audit` | Verify source snapshot and migration accounting |
| `06_explorer_embeddings` | Inspect two-dimensional profile coordinates |

Start Jupyter from the V4 repository root. Run provider calls only through an explicitly configured adapter and resumable decision ledger. The notebooks use small examples for schemas and query inspection; they do not fabricate an acquisition result.
