# Pathogen Atlas V4 — V3 migration and host-aware design TODO

## Purpose

Build V4 as a host-aware successor to V3. Preserve V3’s auditable, versioned evidence workflow while making the host an explicit part of run configuration, evidence identity, and feature applicability. V4 should support pathogen–host pairs across human, livestock, and crop hosts without hard-coding host assumptions into shared pipeline stages.

V3 references: [repository](https://github.com/fabriziocosta/pathogen-atlas-v3), [workflow](https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/docs/workflow.md), [data model](https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/docs/data_model.md), [feature taxonomy decision](https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/docs/feature_taxonomy_decision.md).

## Design decisions to preserve

- Keep pathogen and host as separate taxonomic entities, each with stable identifiers, accepted names, aliases, taxonomy source, and release provenance.
- Represent evidence about a pathogen–host relationship explicitly. Do not infer infection merely because an abstract mentions a host.
- Keep feature dimensions and categories about what a feature measures. Store host and pathogen applicability as separate selectors.
- Distinguish pathogen-scoped evidence from host-scoped and pathogen–host-pair evidence so host-independent questions need not be rescored for every pair.
- Preserve question, profile, input, model, and run versions in decision records and cache keys.
- Keep unsupported, missing, and not-applicable outcomes distinct; do not treat lack of retrieved evidence as evidence of absence.

## Ordered implementation checklist

### 1. Establish the V4 baseline

- [ ] Create the V4 repository from a pinned V3 commit; record the source commit and migration date.
- [ ] Copy reusable source code, configuration, notebooks, schemas, documentation, and tests; do not copy generated data, secrets, or local caches.
- [ ] Record V3 components retained, replaced, or deferred in a migration log.
- [ ] Follow the [detailed V3 data migration runbook](docs/V3_DATA_MIGRATION.md) to snapshot, checksum, transfer, validate, and stage the multi-gigabyte V3 data before transforming it.
- [ ] Keep V3 unchanged; make V4 changes in V4 only.

### 2. Define host and interaction records

- [ ] Add a normalized host catalogue with host ID, accepted name, aliases, host group, taxonomy source, release, and mapping status.
- [ ] Add a host taxonomy adapter interface so host data sources can vary without branching the shared workflow.
- [ ] Define a pathogen–host pair record with pathogen ID, host ID, and an explicit relationship context or role.
- [ ] Define pair-level evidence records with PMID, title/abstract hashes, evidence passage or span, decision, confidence, and provenance.
- [ ] Specify how ambiguous host mentions, broad host groups, vectors, reservoirs, food contamination, and experimental hosts are represented.

### 3. Add run configuration and selection

- [ ] Extend the configuration shown in Notebook 0 to select host cohorts and pathogen groups or taxa.
- [ ] Support both cohort runs and explicit pathogen–host pair runs.
- [ ] Define pair-generation rules so a run does not silently expand to every possible combination.
- [ ] Validate selected host IDs, pathogen IDs, pair rules, and taxonomy releases before acquisition starts.
- [ ] Store the resolved run configuration and a stable configuration version with every run.

### 4. Generalize retrieval and the disease gate

- [ ] Refactor PubMed query construction to accept pathogen names, host names and aliases, and a versioned query strategy.
- [ ] Preserve query text/hash, all search results, sampling decisions, and fetched abstract status for each eligible pair.
- [ ] Make the disease-attribution gate ask whether the abstract supports the specified pathogen–host relationship.
- [ ] Return structured outcomes for supported, unsupported, uncertain, and insufficient evidence, with evidence spans.
- [ ] Ensure a host mention alone cannot pass the relationship gate.

### 5. Make features modular and host-aware

- [ ] Keep the existing dimension → category → feature organization for the feature catalogue.
- [ ] Add feature scope: pathogen, host, or pathogen–host pair.
- [ ] Add applicability selectors for pathogen group, host group, optional host taxon, and pair context.
- [ ] Compose profiles from ordered modules: core, pathogen-group, host-group, optional host-taxon, and pair-context.
- [ ] Resolve and save the exact ordered feature profile and its version before scoring.
- [ ] Replace human-specific shared question wording with explicit templates where a feature truly generalizes; mark genuinely human-only features as human-scoped.
- [ ] Preserve question versioning and require a new version whenever wording, criteria, scope, or interpretation changes.
- [ ] Keep not_applicable separate from no, and document why a feature is excluded by the resolved profile.

### 6. Update decision storage and resumability

- [ ] Include host ID and pair context in cache identity for pair-scoped decisions.
- [ ] Keep organism-only decision keys for features that do not depend on host context.
- [ ] Include query, prompt, feature, profile, model, and input versions in the decision ledger.
- [ ] Preserve V3 resume, retry, incomplete-run, and release-blocking behavior where applicable.

### 7. Rebuild release tables and metrics

- [ ] Publish separate host, pathogen, pair, literature-match, gate-decision, and feature-decision tables with documented keys.
- [ ] Make counts and indexes explicit about whether they refer to a pathogen, host, pair, or literature search.
- [ ] Reassess normalization and comparability for host cohorts and pathogen groups before exposing cross-group comparisons.
- [ ] Keep reported literature evidence distinct from biological prevalence, risk, or proof of absence.

### 8. Adapt embeddings and the Explorer

- [ ] Decide whether maps are pathogen-only, pair-level, or separate views; document the interpretation of each.
- [ ] Ensure feature vectors use the resolved profile and handle not-applicable features explicitly.
- [ ] Add host cohort and host filters, and show the host context and evidence behind each pair.
- [ ] Retain source links, PMIDs, question/profile versions, and uncertainty in the interface.

### 9. Validate in bounded pilots

- [ ] Start with a small host set, such as cattle, sheep, maize, and rice, and a bounded pathogen cohort.
- [ ] Create labelled examples for host identification, relationship type, disease attribution, and evidence-span quality.
- [ ] Compare results across pathogen groups and host groups; inspect false positives caused by background mentions or ambiguous host names.
- [ ] Check profile resolution, cache separation, resumability, and release completeness before scaling acquisition.
- [ ] Record expert review, known limitations, and changes made after each pilot.

### 10. Document and release

- [ ] Update the data model, configuration reference, workflow, feature policy, migration notes, and Explorer documentation.
- [ ] Record source licences, taxonomy releases, query versions, prompt versions, and model versions.
- [ ] Define a V4 release format and migration path without overwriting V3 releases.
- [ ] Add a V3-to-V4 comparison note explaining which outputs remain comparable and which do not.

## Initial acceptance criteria

- A run can select one or more hosts and pathogen taxa/groups through configuration.
- The same acquisition, decision, and release stages operate through shared interfaces across host groups.
- Every pair-level positive has a PMID and a supporting passage tied to a normalized pathogen and host.
- Feature profiles are resolved from applicability rules, are reproducible, and distinguish not-applicable from negative.
- Re-running a partial V4 run resumes only decisions whose pair, inputs, and relevant versions match.
- A small pilot can be audited end to end before broad data acquisition.

## Initial pilot hosts

Use cattle, sheep, maize, and rice as candidate starting hosts, subject to checking accepted names, aliases, and taxonomy identifiers in the selected host catalogue.

