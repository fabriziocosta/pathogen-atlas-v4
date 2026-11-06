# V4 workflow

## Run sequence

0. **Configure the run.** Select pathogen groups, hosts, run mode, and PMID sample settings in `notebooks/00_configure_run.ipynb`. The notebook defaults to GPT-6.1 Terra for host-mention expansion, previews the returned cues, and saves their model, prompt, terms, MeSH terms, and effective-term hash in the run config.
1. **Review taxonomy.** Inspect the pinned host and pathogen taxonomy in `notebooks/01_review_taxonomy_and_hosts.ipynb`.
2. **Preview questions.** See the exact versioned Jev questions in `notebooks/02_preview_prompts.ipynb`; this makes no provider calls.
3. **Acquire literature.** `notebooks/03_acquire_pair_pubmed_records.ipynb` searches each pathogen–host pair. Queries use accepted pathogen names, mapped nomenclatural synonyms, and taxon-checked ICTV Virus name(s) and Virus name abbreviation(s) fields that resolve uniquely across taxa; raw isolate/strain labels, unresolved names, and cross-taxon homonyms remain in taxonomy provenance but are excluded from pathogen-name search. The host clause combines NCBI host names, LLM-proposed host-mention cues, and configured MeSH headings. All pathogen alternatives are ORed together; all host alternatives are ORed together; the pathogen clause is ANDed with the host clause. The notebook retains complete PMID result sets and selects a deterministic abstract sample per pair (30 lowest PMID IDs plus a stable hash-ranked remainder, 200 total by default). Abstract fetches deduplicate sampled PMIDs across pairs and request up to 500 PubMed records per EFetch call. Increasing the sample cap preserves the prior sample and fetches only newly selected abstracts.
4. **Screen pathogen status.** `notebooks/04_screen_pathogen_host_pairs.ipynb` assesses sampled abstracts for direct evidence that the pathogen causes disease in the selected host. Jev answers a typed yes/no question and returns the probability of yes; code accepts a sample record at 0.7 or higher. A title alone cannot establish a positive. The ledger keeps each PMID and screened-abstract hash; Jev does not generate text spans. Missing abstracts remain incomplete evidence.
5. **Estimate publication volume.** `notebooks/05_estimate_pathogen_publications.ipynb` estimates disease-supporting publications from the positive sample rate and pair-search total, and calculates the pair-specific pathogenicity fraction.
6. **Handoff disease evidence.** `notebooks/06_handoff_disease_evidence.ipynb` validates screening coverage and materializes pair decisions linked to their screened abstracts.
7. **Resolve feature profiles.** `notebooks/07_resolve_feature_profiles.ipynb` composes ordered core, pathogen-group, host-group, optional host-taxon, and pair-context modules and records the reasons features do not apply.
8. **Score features.** `notebooks/08_score_feature_profiles.ipynb` sends every trait-applicable feature question to Jev for each accepted pathogen–host PMID. Scope-aware task identities allow pathogen and host questions to be reused when their entity and abstract match. The checkpoint ledger resumes interrupted scoring; release data generation rejects accepted evidence with unresolved applicable questions.
9. **Validate release.** `notebooks/09_validate_release.ipynb` checks keys and references across taxonomy, pair, literature, gate, and feature tables.
10. **Build and inspect the Explorer.** `notebooks/10_inspect_embeddings.ipynb` builds the additive Explorer database after release validation, fits and inspects eligible host × pathogen-group IMSTE maps, and stores coordinates with fit provenance. Step 0 records the staging snapshot and later steps reuse that exact snapshot. Maps are never shared across unrelated hosts or pathogen groups.
11–12. **Resolve references.** `notebooks/11_resolve_wikipedia_references.ipynb` and `notebooks/12_resolve_gbif_references.ipynb` checkpoint exact external references.

## Prompt preview

Run `notebooks/02_preview_prompts.ipynb` after Step 0 and taxonomy loading to print the exact versioned yes/no pathogen-screen request and one pair-scoped feature request for a selected pathogen–host context. The preview accepts an optional title and abstract and does not call a model provider.

## Pair selection

Runs accept host cohorts, pathogen cohorts, or explicit pairs. A Cartesian product is allowed only when the configuration names that rule; otherwise, the run must enumerate the intended pairs. The current Rice-virus pilot selects Rice and the virus group, then explicitly forms the configured host–pathogen Cartesian product.

## Host-mention expansion and PubMed use

The host-expansion prompt is generic: it receives a host ID, accepted scientific name, known taxonomy aliases, and a common label, then proposes short phrases that authors may use to indicate that host in biomedical writing. GPT-6.1 Terra is the notebook default (`HOST_SEARCH_EXPANSION_MODEL`); the model can be changed in the first code cell. Running the expansion requires `OPENAI_API_KEY`. The prompt explicitly asks for common names and person, patient, age, and sex descriptors where appropriate. Its returned terms are search cues, not taxonomic synonyms and not evidence that an article concerns the host.

Each returned cue is normalized and checked for length, duplicate aliases, control characters, and PubMed field-tag injection. The exact candidate list, category, rationale, precision-risk label, model, prompt version, and prompt hash are written into the run configuration for review and reproducibility. Terms are quoted and searched in PubMed `[Title/Abstract]`; this field also covers author keywords. Configured controlled-vocabulary headings are searched separately as `[MeSH Terms]`. For Human, the default MeSH heading is `Humans`. The query structure is `(pathogen name OR pathogen synonym ...) AND (host name OR host alias OR host cue ... OR host MeSH heading ...)`.

The expansion intentionally favors recall. Broad cues such as “patient” can retrieve irrelevant records, so query hits remain candidates. The existing abstract-level disease screen must still establish the exact pathogen–host relationship. MeSH adds wording-independent retrieval for indexed MEDLINE records, while text cues remain available for records without MeSH indexing. The effective host-term hash is part of the shared literature-cache key, and pair checkpoints are rerun when a pair's generated query hash changes.

The expansion is generic across hosts, but its vocabulary must be reviewed for the selected host and language. Keep semantic classes in the output (`common_name`, `person_or_patient_descriptor`, `age_or_sex_descriptor`, `host_context_phrase`) so demographic expressions are not mislabeled as taxonomic aliases. Review missed known papers and false-positive samples before broad production runs.

## Evidence and interpretation

The pair gate returns `supported`, `unsupported`, `uncertain`, or `insufficient_evidence`. Jev System 1 screens return a typed yes/no probability rather than generated text, so a positive decision cites its PMID and screened-abstract hash instead of a generated passage. Text-generating providers must still supply a nonempty exact passage. Lack of retrieved evidence remains distinct from an unsupported assessment. Vector, reservoir, experimental-host, and food-contamination mentions may be recorded as annotations but do not count as disease pairs.

## Resumability

Each run stores the resolved configuration and feature profile before scoring. Pathogen-scoped decisions do not vary by host; host-scoped decisions key on the host; pair-scoped decisions key on pathogen, host, and relationship context. All scopes include the question, profile, model, and input versions needed to reproduce a result. Feature scoring covers every configured question for each accepted pair abstract. Shared literature caches are keyed by taxonomy and query versions, while shared decision caches also require matching input, prompt, model, and configuration versions. A full run can therefore reuse overlapping smoke-run records when those identities match; increasing the sample cap adds only the newly selected records.

## Current implementation status

The repository provides V4 contracts, explicit cohort/pair validation, resolved-run persistence, NCBI host-taxonomy adapter, versioned PubMed acquisition, deterministic per-pair abstract samples, resumable host-specific pathogen screening, accepted-pathogen publication counts and PII calculation, trait-aware profile resolution, resumable Jev feature scoring, evidence-span validation, release checks, the deterministic data importer, and a host-filtered Explorer. The Explorer preserves prior host and pathogen profiles, includes literature-curated rice-virus traits plus pair-specific Jev feature scores, and stores host × pathogen-group IMSTE coordinates with PMID links.
## Pathogenic Interest Index population

For each accepted pathogen–host pair, estimate publications that assert disease as `(positive sampled abstracts / selected sample records) × full pair-search results`. The sample denominator is the selected PMID sample (200 by default); when the query returns 200 or fewer records, every result is sampled. Every available abstract in that sample must be screened. Divide the estimate by the full pair-search result count to obtain the pair-specific pathogenicity fraction. Keep the organism-only PubMed total separate; it is not the denominator for this pair-specific fraction. Save both estimates and their component counts.
