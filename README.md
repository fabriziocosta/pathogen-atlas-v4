<div align="center">
  <img src="assets/pathogen-atlas-v4-logo.png" alt="Pathogen Atlas V4 logo" width="170">
</div>

# Pathogen Atlas V4

The Pathogen Atlas brings together published evidence about disease relationships between pathogens and their hosts. V4 represents hosts explicitly across humans, livestock, and crops, while preserving the taxonomy, literature sources, and assessment history behind each record.

Pathogen selection covers bacteria, fungi, viruses, parasites, and oomycetes. Parasites and oomycetes are sourced from the pinned NCBI Taxonomy release; their selected taxonomic scopes are described in [pathogen groups](docs/pathogen_groups.md).

## What the Atlas shows

For each pathogen–host pair, the Atlas connects taxonomy, literature matches, assessed evidence, and relevant pathogen or host features. Records retain their PubMed identifiers, supporting text, confidence, uncertainty, and the versions of questions and profiles used to assess them.

The Explorer connects disease evidence to its host and pathogen context. Feature profiles describe what the literature reports and are resolved for each applicable host and pathogen group. Host-independent features can be shared across pairs; host-specific and pair-specific features are assessed in their resolved context.

## Explore

The [Pathogen Atlas Explorer](https://fabriziocosta.github.io/pathogen-atlas-v4/) helps users find pathogens and examine the disease evidence associated with each host. Search and filter by pathogen group and host, review linked PubMed records and their titles, and explore feature profiles in their pathogen–host context.

## Project structure

- `pathogen_atlas_v4/` contains the versioned data contracts, host taxonomy adapter, run configuration, evidence decisions, and release checks.
- `configs/` contains the run defaults, feature catalogue, and bounded pilot cohort.
- `scripts/` contains reproducible source and pilot preparation helpers.
- `docs/` records the data model and workflow.
- `notebooks/` provides modular guides for taxonomy, run resolution, pair evidence, feature profiles, release validation, and reference resolution.
- `data/` holds local source snapshots, caches, staging tables, and releases; data and credentials are not committed to Git.

See the [V4 workflow](docs/workflow.md) for the acquisition sequence, the [V4 data model](docs/data_model.md) for the host, pathogen, pair, and feature contracts, and the [feature policy](docs/feature_policy.md) for question scope and review.

## Local pipeline dashboard

The local V4 dashboard coordinates the notebook pipeline, stores run state under ignored `data/v4/dashboard/runs/`, and can resume a paused run from saved phase checkpoints. It refreshes PubMed on a fresh run while reusing cached abstracts and semantically matching decisions.

Install the dashboard and embedding extras into the V4 environment, then start the server:

```bash
python -m pip install -e '.[dashboard,embedding]'
./scripts/start_dashboard.sh
```

The launcher starts two local servers, opens both pages, and stops them when you press Ctrl+C:

- Smoke test: `http://127.0.0.1:8765` — fixed smoke mode; choose one host, pathogen group, and species.
- Full cohort: `http://127.0.0.1:8766` — fixed full mode; choose hosts and pathogen groups from the YAML-backed menus.

Both pages show the independently maintained dashboard version (`DASHBOARD_VERSION` in `pathogen_atlas_v4/dashboard.py`) and bind to loopback only. Provider credentials stay in the local environment or existing system keychain. A completed full run exposes **Publish to Explorer** for an explicit staged-file review and release action. Smoke runs cannot be published; the dashboard never publishes automatically.

The dashboard offers **Smoke test** and **Full cohort** modes. A smoke test selects one host, one pathogen group, and one pathogen species, screens at most 10 PMIDs, and runs the disease screen, feature scoring, validation, reference lookups, and an isolated Explorer build. It uses a fresh temporary workspace instead of shared run caches or production Explorer files; after successful completion, smoke data and logs are removed and the run cannot be published. Full cohort mode covers every pathogen in the selected groups, with up to 200 PMIDs per pair.
