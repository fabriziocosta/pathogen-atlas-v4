# V4 workflow

## Run sequence

1. **Pin inputs.** Load the host taxonomy source and pathogen catalogue, recording source releases and checksums.
2. **Resolve configuration.** Validate selected host and pathogen IDs, taxonomy versions, and explicit pair-generation rules. Save the resolved run configuration before acquisition.
3. **Acquire literature.** Build a versioned PubMed query for each pathogen–host pair from accepted names and aliases. Retain the complete PMID search result list, then select a deterministic abstract-screen sample per pair (30 lowest PMID IDs plus a stable hash-ranked remainder, 200 total by default). Increasing the cap preserves the previous sample and fetches only newly selected abstracts.
4. **Screen pathogen status by host.** Assess each sampled abstract for direct evidence that the named pathogen causes disease in the selected host. One yes with an exact supporting span from the abstract accepts the pair. A title alone cannot establish a positive. Pairs with no positive are rejected only after all sampled records with abstracts have been assessed; missing abstracts remain incomplete evidence.
5. **Resolve features.** Compose ordered core, pathogen-group, host-group, optional host-taxon, and pair-context modules. Persist the exact resolved profile and reasons for not-applicable features.
6. **Score and resume.** Reuse a cached decision only when its scoped entity, inputs, question, feature, profile, model, and relevant configuration versions match. Record failed and incomplete work distinctly for retry.
7. **Validate and release.** Check keys and references across taxonomy, pair, literature, gate, and feature tables. Publish only complete runs; label all counts by whether they describe a pathogen, host, pair, search, or decision.
8. **Explore.** Generate profile-aware feature views and embeddings. Keep host context, PMID links, uncertainty, and feature applicability visible in the Explorer.

## Prompt preview

Run `notebooks/01_prompt_preview.ipynb` after Phase 0 and taxonomy loading to print the exact versioned yes/no pathogen-screen request and one pair-scoped feature request for a selected pathogen–host context. The preview accepts an optional title and abstract and does not call a model provider.

## Pair selection

Runs accept host cohorts, pathogen cohorts, or explicit pairs. A Cartesian product is allowed only when the configuration names that rule; otherwise, the run must enumerate the intended pairs. The pilot configuration contains four host IDs and 30 fixed pathogen IDs and enumerates all 120 candidate pairs.

## Evidence and interpretation

The pair gate returns `supported`, `unsupported`, `uncertain`, or `insufficient_evidence`. A supported native pair must cite a PMID and nonempty text passage tied to the normalized host and pathogen. Lack of retrieved evidence remains distinct from an unsupported assessment. Vector, reservoir, experimental-host, and food-contamination mentions may be recorded as annotations but do not count as disease pairs.

## Resumability

Each run stores the resolved configuration and feature profile before scoring. Pathogen-scoped decisions do not vary by host; host-scoped decisions key on the host; pair-scoped decisions key on pathogen, host, and relationship context. All scopes include the question, profile, model, and input versions needed to reproduce a result.

## Current implementation status

The repository provides V4 contracts, explicit cohort/pair validation, resolved-run persistence, NCBI host-taxonomy adapter, versioned PubMed acquisition, deterministic per-pair abstract samples, resumable host-specific pathogen screening, accepted-pathogen publication counts and PII calculation, profile resolution, evidence-span validation, release checks, the deterministic data importer, and a host-filtered Explorer. The Explorer includes human bacteria, fungi, and virus profiles with two-dimensional coordinates and title-linked PMID records. Model-provider handlers, native feature scoring, and the V4 embedding pipeline remain modular backend work.
## Pathogenic Interest Index population

The PII uses the V3 sampling estimator for each accepted pathogen–host pair. Query the total organism-name PubMed record count only for unique pathogens accepted in at least one selected host. Then estimate host-specific pathogen publications as `(positive sampled abstracts / abstracts screened) × total organism PubMed records`. The PII is the log of that estimate, normalized within pathogen group. The saved row includes the sample numerator and denominator, total count, and estimate.
