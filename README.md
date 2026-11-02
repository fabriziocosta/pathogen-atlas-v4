<div align="center">
  <img src="assets/pathogen-atlas-v4-logo.png" alt="Pathogen Atlas V4 logo" width="170">
</div>

# Pathogen Atlas V4

The Pathogen Atlas brings together published evidence about disease relationships between pathogens and their hosts. V4 represents hosts explicitly across humans, livestock, and crops, while preserving the taxonomy, literature sources, and assessment history behind each record.

## What the Atlas shows

For each pathogen–host pair, the Atlas connects taxonomy, literature matches, assessed evidence, and relevant pathogen or host features. Records retain their PubMed identifiers, supporting text, confidence, uncertainty, and the versions of questions and profiles used to assess them.

The Explorer connects disease evidence to its host and pathogen context. Feature profiles describe what the literature reports and are resolved for each applicable host and pathogen group. Host-independent features can be shared across pairs; host-specific and pair-specific features are assessed in their resolved context.

## Explore

The [Pathogen Atlas Explorer](https://fabriziocosta.github.io/pathogen-atlas-v4/) helps users find pathogens and examine the disease evidence associated with each host. Search and filter by pathogen group and host, review linked PubMed records and their titles, and explore feature profiles in their pathogen–host context.

## Project structure

- `pathogen_atlas_v4/` contains the versioned data contracts, host taxonomy adapter, run configuration, evidence decisions, and release checks.
- `configs/` contains the run defaults, feature catalogue, and bounded pilot cohort.
- `scripts/` contains reproducible source and pilot preparation helpers.
- `docs/` records the data model, workflow, and migration methods.
- `data/` holds local source snapshots, caches, staging tables, and releases; data and credentials are not committed to Git.

See the [V4 workflow](docs/workflow.md) for the acquisition sequence, the [V4 data model](docs/data_model.md) for the host, pathogen, pair, and feature contracts, and the [feature policy](docs/feature_policy.md) for question scope and review.
