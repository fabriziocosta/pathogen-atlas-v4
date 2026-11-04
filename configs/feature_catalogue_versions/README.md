# Feature catalogue snapshots

`configs/features.yaml` is the active catalogue used by the V4 notebooks. Each released catalogue version is also saved here as an immutable YAML snapshot. The baseline `v4-feature-catalogue-1.0.8.yaml` is preserved; `v4-feature-catalogue-2.0.2.yaml` is the preliminary-active version.

The 2.0.2 version includes all questions currently listed in the configured modules. Later versions can remove questions from module membership, revise questions with a new `question_version`, or add new feature IDs. Bump both the catalogue and profile versions for a new snapshot, then update the run configuration profile version so decisions and caches remain version-aware. A `curated-draft` question is active in this preliminary version but retains its review provenance until domain experts review it.
