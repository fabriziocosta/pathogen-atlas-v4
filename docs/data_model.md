# V4 data model

## Taxonomy entities

`taxonomy/pathogens.parquet` and `taxonomy/hosts.parquet` are separate catalogues. Each record carries a stable V4 ID, accepted name, source taxon ID, rank, parent, taxonomy source and release, and mapping status. The pathogen catalogue also retains its organism group. Alias tables link names to an entity ID where the mapping is resolved and preserve ambiguous or unresolved status otherwise. Host groups identify cohorts such as human, livestock, crop, or other. The pinned host set includes Human, Cattle, Sheep, Chicken, Pig, Maize, Rice, Wheat, and Soybean.

## Pairs and literature

A pair is identified by `pathogen_id`, `host_id`, and `relationship_context`. The current relationship context is `disease`; vector, reservoir, experimental-host, and food-contamination mentions remain annotations and cannot create disease positives.

Literature records are PMID-keyed. A search-match row records the pathogen, host, query text hash and version, sampling decision, and fetch status. An article row may contain fetched metadata or an explicit placeholder when a PMID was found but its abstract was not fetched. A missing abstract is not a negative disease decision.

Pair-evidence rows link an article to normalized pathogen and host IDs and include the outcome, confidence, source passage, input hashes, and provenance. Outcomes are `supported`, `unsupported`, `uncertain`, and `insufficient_evidence`. Native supported evidence requires a nonempty passage. Legacy evidence is clearly labelled and retains its original source record and scoring versions.

## Features and decisions

The feature catalogue keeps the dimension → category → feature hierarchy. Each feature declares its scope (`pathogen`, `host`, or `pathogen_host_pair`), question version, applicability selectors, and review state. A resolved profile records the exact ordered feature list, its selectors, exclusions, profile version, and a stable hash.

Decision identity includes the feature, input hash, query, prompt, feature, question, profile, model, run, and configuration versions. Host identity is included for host-scoped and pair-scoped decisions. Pair scope also includes relationship context. Pathogen-scoped cache identities remain reusable across hosts.

## Release table keys

| Table | Key |
| --- | --- |
| Pathogens | `pathogen_id` |
| Hosts | `host_id` |
| Pairs | `pathogen_id`, `host_id`, `relationship_context` |
| Articles | `pmid` |
| Literature matches | `pathogen_id`, `host_id`, `pmid`, `query_hash` |
| Pair evidence | Pair key, `pmid`, `question_version`, `input_hash` |
| Feature decisions | Scoped entity key, `pmid`, `feature_id`, `input_hash`, question/profile/model versions |

Counts are labelled by the entity or acquisition stage they measure. They describe retrieved and assessed literature, not disease prevalence, risk, or absence.
## Pathogenic Interest Index population

The pair-screen sample is 200 PMIDs per pathogen–host pair by default: 30 lowest PMID identifiers followed by a reproducible hash-ranked random sample from the remainder. The random ordering is prefix-stable, so increasing the cap extends the previous sample and reuses cached abstracts. A positive pathogen–host pair requires at least one exact supporting passage from an abstract; a title alone is insufficient.

The PII population is the set of unique pathogens accepted in at least one selected host. Query PubMed for total organism-name publications only for that set, using count-only searches. Do not infer the total from the pair-screen sample. The PII is `log1p(total PubMed records)`, normalized to 0–1 within pathogen group.
