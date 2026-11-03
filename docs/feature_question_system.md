# V4 feature-question system

## Purpose

V4 feature questions turn information that can be explicitly read in an abstract into structured, traceable observations. They are not a general biological knowledge quiz. A positive answer means the supplied abstract reports the claim for the entity or relationship in scope. A negative answer means only that this abstract did not support the claim. Missing evidence must never be interpreted as biological absence.

The feature-question catalogue is a core V4 distinction from V2/V3: each scored statement has a stable identity, a declared scope, a place in a shared hierarchy, explicit applicability rules, and evidence provenance. That makes the same question reusable where it is appropriate, while excluding it from pairs where the question does not make sense.

## Unit of observation and scope

Questions have one of three scopes:

| Scope | What one answer describes | Example |
|---|---|---|
| `pathogen` | A property reported about the named pathogen, independent of a particular host | A bacterial species is reported to form a biofilm |
| `host` | A property reported about the named host, independent of a pathogen | A host is reported to be a reservoir for a pathogen group |
| `pathogen_host_pair` | A property reported for the particular pathogen in the particular host | A virus causes mosaic symptoms in rice |

Use the narrowest scope that matches the claim. Do not turn pathogen-level facts into pair-level facts, or carry a pair-level answer over to another host. Host and pathogen intrinsic answers may be reused across pair profiles only when their scope and provenance permit it.

Questions are intended to characterize pathogenic aspects that abstracts can reasonably establish: pathogen properties relevant to disease, pathogen–host interaction, routes and persistence, tissue involvement, disease signs, and host impact. Each question must be answerable from the title and abstract record available to the pipeline. For a positive feature answer, the abstract must explicitly support the claim and the named entity or pair must be identifiable. The title may help identify an article but is not sufficient evidence on its own. Full-text-only results, plausible biological mechanisms, taxonomic expectations, and unstated implications do not count as abstract support.

## Shared question hierarchy

Questions are organized as **dimension → category → question**. Dimensions group broad parts of the pathogen–host evidence space; categories provide useful subtopics; individual questions are the atomic scoring unit. The hierarchy is for navigation and analysis, not a claim that every question is applicable to every organism.

The catalogue retains existing human question IDs and wording in `legacy_human_profile_membership`, and lists the human questions under the human module. The resolver applies host and pathogen traits. The 12 `rice_*` questions resolve only for rice (`host:rice`, taxon 4530); they describe rice symptoms and rice virus vectors. Shared human disease questions about clinical presentation, exposures, and transmission now include oomycetes alongside bacteria, fungi, and viruses. Fungal-, viral-, and bacterial-specific questions remain restricted to their respective groups. This expands the applicable human profile for oomycetes without treating oomycetes as fungi. Human questions should be revised in a new question version if their meaning changes; do not silently repurpose their IDs.

## Traits and applicability

Host acquisition groups (`human`, `livestock`, `crop`) are useful selectors, but they are not the biological host ontology. V4 therefore adds reusable traits. For example, livestock and human hosts both receive `host:animal`, `host:vertebrate`, and `host:mammal`; crop hosts receive `host:plant` and `host:crop`. Individual taxa may receive narrower curated traits such as `host:cereal_crop`.

Pathogen groups similarly imply traits such as `pathogen:virus`, `pathogen:bacterium`, and `pathogen:fungus`. A question’s `applicability` can constrain pathogen, host, and pair traits independently, using positive `all_of` and `any_of` predicates. For example, a plant symptom question requires `host:plant`; an animal clinical-system question requires `host:animal`; a bacterium-specific question requires the bacterial pathogen group. A question may combine these constraints.

An absent trait means **unknown**, not “does not have this trait.” The resolver uses positive predicates only, and exclusions record whether a question was omitted because its selector did not match or because its required traits were unavailable. `applies_to_label` is for people reading the resolved table; trait predicates, not the label, control selection. This allows one plant question to resolve for rice, wheat, maize, and soybean, while a rice-virus question remains limited to the rice taxon and virus group. Add taxon- or pair-specific traits only from a recorded taxonomy or curated source; do not infer them from the abstract being scored.

The profile is composed in a stable order: core questions, pathogen group, host group, host taxon, then pair context. The resolver then applies trait predicates, rejects duplicate question IDs, saves a trait snapshot, and hashes the resolved profile. Thus a shared animal question can be authored once and apply to both humans and livestock, while an additional human-only question can remain restricted to `host:human`.

## Question record and evidence contract

Each question record contains:

- `feature_id`: stable machine identity.
- `label` and `definition`: concise UI and analyst descriptions.
- `question` and `question_claim`: the complete prompt and the claim that a positive answer asserts.
- `scope`: `pathogen`, `host`, or `pathogen_host_pair`.
- `applies_to_label`: a readable reusability class such as `All viruses`, `All plants · any pathogen`, or `Rice viruses`; the trait predicates remain authoritative.
- `dimension_id`, `category_id`, and display labels: position in the hierarchy.
- `applicability`: positive trait requirements.
- `question_version`: version of wording, scope, and decision meaning.
- `answer_type`: currently `noul` (yes/no plus uncertain/insufficient outcomes supported by the decision workflow).
- `inclusion_examples` and `exclusion_examples`: schematic guidance; examples do not replace evidence.
- `evidence_policy`: abstract source, positive evidence requirement, and interpretation of a negative.
- `review_status`: for example `inherited-human-catalogue`, `curated-draft`, or `reviewed`.

For every positive answer, preserve the PMID and a claim-linked evidence passage. “Not stated,” “uncertain,” “insufficient abstract,” and “not applicable” are distinct from a supported positive. Do not convert an unscored or inapplicable item into a negative.

## How to write questions people can answer from abstracts

1. Ask about a report in the supplied abstract, not whether a biological property is universally true.
2. Name the entity scope in the claim; for pair scope, identify both pathogen and host.
3. Keep one observable claim per question. Split bundled claims such as “transmitted by seed or insects.”
4. Prefer outcomes, contexts, routes, and properties that papers commonly state in abstracts. Avoid claims requiring full text, specialist inference, or an unstated causal mechanism.
5. Use a negative only for “not supported by this abstract.” It is not proof of absence.
6. Make applicability stricter than wording alone: do not ask about plant symptoms of an animal host, or bacterial spores of a virus.
7. Treat the catalogue as versioned scientific measurement. A meaning change requires a new question version and a refreshed resolved-profile hash.

Candidate questions are deliberately marked `curated-draft` and `abstract_observability: review_candidate` until review confirms that abstracts commonly provide enough information and that the claim can be judged consistently. Before broad scoring, review a diverse sample of abstracts per question, refine low-observability or ambiguous items, and promote only questions with stable interpretation. The human catalogue is inherited as the reference style, not presumed to be equally suitable for every host type. Prefer a shorter set of explicit, answerable pathogenic-aspect questions over a large set of speculative or generic questions.

## Configuration and operating notes

`configs/features.yaml` is the canonical editable catalogue. `trait_assignments` maps acquisition groups and known host taxa to reusable traits. `modules` list question IDs, and `questions` stores each question once. Keep the explorer’s separate human display catalogue synchronized only when the public feature browser is intentionally updated; the V4 acquisition catalogue remains the source used to resolve V4 scoring profiles.

The current profiles cover humans, livestock, and four crop taxa represented in the V4 taxonomy. Plant, animal, and pathogen-specific questions are present as candidate questions even if some pathogen-group data or a particular host–pathogen pair has not yet been acquired. A question catalogue does not itself establish a disease relationship or create evidence records; pair-specific questions should be used in the disease-screening workflow and retained alongside that pair’s selected profile.
