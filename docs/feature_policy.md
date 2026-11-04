# V4 feature policy

The full question design, traits, hierarchy, and migration contract are documented in [the V4 feature-question system](feature_question_system.md). The canonical catalogue is `configs/features.yaml`.

## What a feature states

Each feature describes an observable property of the literature evidence. Its dimension and category describe what is measured; host or pathogen applicability is recorded separately. A feature answer is not automatically a claim about biological truth or disease prevalence.

## Scope and applicability

Every feature has one scope:

- `pathogen`: a pathogen property that does not depend on host context.
- `host`: a host property reusable across pathogen pairs.
- `pathogen_host_pair`: a property of the specific disease relationship.

Profiles are assembled in order from core, pathogen-group, host-group, optional host-taxon, and pair-context modules. The resolved feature list and its version are saved before scoring. Excluded features are recorded as `not_applicable` with the selector that excluded them.

## Wording and evidence

Questions specify the entity scope, the claim being assessed, the permitted evidence, and the interpretation of yes/no or other answer types. Host mentions alone do not support a disease relationship. Pair-level positives cite a PMID and a source passage. Unsupported, uncertain, insufficient-evidence, missing, and not-applicable outcomes remain distinct.

## Versioning and review

Every question has a stable feature ID and question version. Wording, criteria, scope, or interpretation changes require a new question version. The resolved profile has its own version and stable hash. Catalogue version `v4-feature-catalogue-2.0.2` is marked `preliminary-active`: all questions in its configured modules are included for use, including questions that retain `review_status: curated-draft`. This is a preliminary operational promotion, not a claim of domain-expert validation. Later catalogue versions may remove or revise questions; prior snapshots are preserved under `configs/feature_catalogue_versions/`.
