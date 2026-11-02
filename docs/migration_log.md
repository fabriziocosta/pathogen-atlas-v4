# V4 migration log

## Initial data import

- Migration date: 2026-11-02.
- Source repository commit: `59409a9136bb1ce8477bc80179c853a9132188a8`.
- Source data snapshot: 147 files, 4,585,858,645 bytes; SHA-256 manifest verified after copy and before conversion.
- Host taxonomy: NCBI Taxonomy taxdump retrieved 2026-11-02; archive SHA-256 `f278697a8ba610bd415eab3bd409afb9bd5dda6326b14d5430ce03580e072a0a`.
- Conversion output: `data/v4/staging/v3-to-v4-20261102-ncbi-f278697a8ba6-v3/`; conversion ID `8b66c7e44bfc97b3a20b`.

## Imported records

- 142,201 pathogen taxonomy rows, 233,550 aliases, and 19,285 virus variants.
- 1,733,750 article identities: 163,687 abstract records plus 1,570,063 explicit `not_fetched_in_v3` placeholders.
- 3,450,866 pathogen–article search matches, all linked to an existing pathogen and article identity.
- Five normalized NCBI host records: human, cattle, sheep, maize, and rice.
- 2,344,486 historical decision records and 510,641 V3 feature assessments retained as legacy decisions.
- 2,674 accepted legacy human disease-evidence rows link to valid pathogens, hosts, and PMIDs. Their evidence passage is the full source abstract and is marked `full_abstract_legacy`; it is not a model-extracted span.
- Additional release, embedding, cache, raw-response, and source files remain byte-preserved in the immutable snapshot. Provenance indexes all 147 source files plus 989,045 JSONL records and 8,146,144 PMID references by query hash and PMID where available.

## Validation and release state

The local snapshot passed its manifest check. The staging tables passed pathogen, host, article, and pair-evidence reference checks; every supported legacy pair has a nonempty source passage. The pilot config resolves 30 pathogen IDs (10 per group) into exactly 120 explicit pairs across the four pilot hosts.

The imported tables remain in staging. They are not a public V4 release: native pair acquisition and scoring, expert-labelled pilot review, feature-complete releases, embeddings, and the Explorer still need implementation and review.
