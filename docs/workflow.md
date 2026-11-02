# V4 workflow

## Run sequence

1. **Pin inputs.** Load the host taxonomy source and pathogen catalogue, recording source releases and checksums.
2. **Resolve configuration.** Validate selected host and pathogen IDs, taxonomy versions, and explicit pair-generation rules. Save the resolved run configuration before acquisition.
3. **Acquire literature.** Build a versioned PubMed query from pathogen and host accepted names and aliases. Store query text/hash, results, sampling choices, and abstract-fetch status.
4. **Assess disease attribution.** Evaluate the named pathogen–host relationship. A host mention alone cannot pass the gate. Save one structured outcome and the supporting passage and source hashes for each completed assessment.
5. **Resolve features.** Compose ordered core, pathogen-group, host-group, optional host-taxon, and pair-context modules. Persist the exact resolved profile and reasons for not-applicable features.
6. **Score and resume.** Reuse a cached decision only when its scoped entity, inputs, question, feature, profile, model, and relevant configuration versions match. Record failed and incomplete work distinctly for retry.
7. **Validate and release.** Check keys and references across taxonomy, pair, literature, gate, and feature tables. Publish only complete runs; label all counts by whether they describe a pathogen, host, pair, search, or decision.
8. **Explore.** Generate profile-aware feature views and embeddings. Keep host context, PMID links, uncertainty, and feature applicability visible in the Explorer.

## Pair selection

Runs accept host cohorts, pathogen cohorts, or explicit pairs. A Cartesian product is allowed only when the configuration names that rule; otherwise, the run must enumerate the intended pairs. The pilot configuration contains four host IDs and 30 fixed pathogen IDs and enumerates all 120 candidate pairs.

## Evidence and interpretation

The pair gate returns `supported`, `unsupported`, `uncertain`, or `insufficient_evidence`. A supported native pair must cite a PMID and nonempty text passage tied to the normalized host and pathogen. Lack of retrieved evidence remains distinct from an unsupported assessment. Vector, reservoir, experimental-host, and food-contamination mentions may be recorded as annotations but do not count as disease pairs.

## Resumability

Each run stores the resolved configuration and feature profile before scoring. Pathogen-scoped decisions do not vary by host; host-scoped decisions key on the host; pair-scoped decisions key on pathogen, host, and relationship context. All scopes include the question, profile, model, and input versions needed to reproduce a result.

## Current implementation status

The repository currently provides the V4 contracts, explicit cohort/pair validation, resolved-run persistence, NCBI host-taxonomy adapter, versioned query builder, profile resolver, scoped cache identity, a provider-agnostic resumable decision runner, evidence-span validation, release checks, pilot configuration, deterministic data importer, and a copied Explorer shell. A live PubMed acquisition client, provider integration, feature scoring pipeline, embeddings, and V4 Explorer data wiring are the next implementation stages; the contracts and acceptance checks are in place for them.
