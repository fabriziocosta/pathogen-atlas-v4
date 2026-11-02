# V4 workflow

## Run sequence

0. **Configure the run.** Select pathogen groups, hosts, run mode, and PMID sample settings in `notebooks/00_configure_run.ipynb`.
1. **Review taxonomy.** Inspect the pinned host and pathogen taxonomy in `notebooks/01_review_taxonomy_and_hosts.ipynb`.
2. **Preview questions.** See the exact versioned Jev questions in `notebooks/02_preview_prompts.ipynb`; this makes no provider calls.
3. **Acquire literature.** `notebooks/03_acquire_pair_pubmed_records.ipynb` searches each pathogen–host pair. Queries use accepted pathogen names and mapped, nomenclatural synonyms that are unique across taxa; ICTV strain/virus aliases, unresolved names, and cross-taxon homonyms remain in taxonomy provenance but are excluded from search. The notebook retains complete PMID result sets and selects a deterministic abstract sample per pair (30 lowest PMID IDs plus a stable hash-ranked remainder, 200 total by default). Increasing the cap preserves the prior sample and fetches only newly selected abstracts.
4. **Screen pathogen status.** `notebooks/04_screen_pathogen_host_pairs.ipynb` assesses sampled abstracts for direct evidence that the pathogen causes disease in the selected host. Jev answers a typed yes/no question and returns the probability of yes; code accepts a sample record at 0.7 or higher. A title alone cannot establish a positive. The ledger keeps each PMID and screened-abstract hash; Jev does not generate text spans. Missing abstracts remain incomplete evidence.
5. **Estimate publication volume.** `notebooks/05_estimate_pathogen_publications.ipynb` estimates disease-supporting publications from the positive sample rate and pair-search total, and calculates the pair-specific pathogenicity fraction.
6. **Handoff disease evidence.** `notebooks/06_handoff_disease_evidence.ipynb` validates screening coverage and materializes pair decisions linked to their screened abstracts.
7. **Resolve feature profiles.** `notebooks/07_resolve_feature_profiles.ipynb` composes ordered core, pathogen-group, host-group, optional host-taxon, and pair-context modules and records the reasons features do not apply.
9. **Validate release.** `notebooks/09_validate_release.ipynb` checks keys and references across taxonomy, pair, literature, gate, and feature tables.
10. **Inspect embeddings.** `notebooks/10_inspect_embeddings.ipynb` generates profile-aware feature views and embeddings while keeping host context, PMID links, uncertainty, and feature applicability visible.
11–12. **Resolve references.** `notebooks/11_resolve_wikipedia_references.ipynb` and `notebooks/12_resolve_gbif_references.ipynb` checkpoint exact external references.

## Prompt preview

Run `notebooks/02_preview_prompts.ipynb` after Step 0 and taxonomy loading to print the exact versioned yes/no pathogen-screen request and one pair-scoped feature request for a selected pathogen–host context. The preview accepts an optional title and abstract and does not call a model provider.

## Pair selection

Runs accept host cohorts, pathogen cohorts, or explicit pairs. A Cartesian product is allowed only when the configuration names that rule; otherwise, the run must enumerate the intended pairs. The pilot configuration contains four host IDs and 30 fixed pathogen IDs and enumerates all 120 candidate pairs.

## Evidence and interpretation

The pair gate returns `supported`, `unsupported`, `uncertain`, or `insufficient_evidence`. Jev System 1 screens return a typed yes/no probability rather than generated text, so a positive decision cites its PMID and screened-abstract hash instead of a generated passage. Text-generating providers must still supply a nonempty exact passage. Lack of retrieved evidence remains distinct from an unsupported assessment. Vector, reservoir, experimental-host, and food-contamination mentions may be recorded as annotations but do not count as disease pairs.

## Resumability

Each run stores the resolved configuration and feature profile before scoring. Pathogen-scoped decisions do not vary by host; host-scoped decisions key on the host; pair-scoped decisions key on pathogen, host, and relationship context. All scopes include the question, profile, model, and input versions needed to reproduce a result. Shared literature caches are keyed by taxonomy and query versions, while shared decision caches also require matching input, prompt, model, and configuration versions. A full run can therefore reuse overlapping smoke-run records when those identities match; increasing the sample cap adds only the newly selected records.

## Current implementation status

The repository provides V4 contracts, explicit cohort/pair validation, resolved-run persistence, NCBI host-taxonomy adapter, versioned PubMed acquisition, deterministic per-pair abstract samples, resumable host-specific pathogen screening, accepted-pathogen publication counts and PII calculation, profile resolution, evidence-span validation, release checks, the deterministic data importer, and a host-filtered Explorer. The Explorer includes human bacteria, fungi, and virus profiles with two-dimensional coordinates and title-linked PMID records. Model-provider handlers, native feature scoring, and the V4 embedding pipeline remain modular backend work.
## Pathogenic Interest Index population

For each accepted pathogen–host pair, estimate publications that assert disease as `(positive sampled abstracts / selected sample records) × pair-search records`. The sample denominator is the selected PMID sample (200 by default); every available abstract in that sample must be screened. Divide the estimate by the organism-only PubMed total to obtain the PII fraction. Query this total only for unique pathogens accepted in at least one selected host. Save both estimates and their component counts.
