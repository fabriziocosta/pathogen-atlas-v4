# Changelog

## [4.0.30] - 2026-11-03

- Render author names in black while keeping the label and link glyphs gray.

## [4.0.29] - 2026-11-03

- Add linked author credits to the Explorer title row.

## [4.0.28] - 2026-11-03

- Place the version label inline with the title text for baseline alignment.


## [4.0.27] - 2026-11-03

- Center Explorer control icons in equal-width slots with balanced button padding.


## [4.0.26] - 2026-11-03

- Align the Explorer version label with the title text baseline.


## [4.0.25] - 2026-11-03

- Place the taxonomy count last and right-align taxonomy source metadata.


## [4.0.24] - 2026-11-03

- Highlight the selected pathogen group in the coverage heading.


## [4.0.23] - 2026-11-03

- Show the complete selected pathogen-group taxonomy count beside the feature count.


## [4.0.22] - 2026-11-03

- Show the selected pathogen taxonomy size, source resource, and source date in the Explorer coverage summary.


## [4.0.21] - 2026-11-03

- Match Explorer control icon size and spacing to V2.

## [4.0.20] - 2026-11-03

- Match Explorer control icons to the Bootstrap Icons used in V2.

## [4.0.19] - 2026-11-03

- Align the host-group CSV download to the far right of the Explorer header.

## [4.0.18] - 2026-11-03

- Export all pathogens in the selected group for the selected host from the main Explorer header.

## [4.0.17] - 2026-11-03

- Match the lower 2D plot colorbar label placement and capitalization to the primary plot.

## [4.0.16] - 2026-11-03

- Add a selected host–pathogen CSV download with taxonomy, pair evidence counts, and feature values and counts.

## [4.0.15] - 2026-11-03

- Display estimated pathogen-supporting publication totals as whole numbers.

## [4.0.14] - 2026-11-03

- Add Atlantic cod and salmon to the host catalogue and pilot selection.
- Replace the chicken in the Pathogen Atlas logo with a fish and refresh the Explorer logo asset.

## [4.0.13] - 2026-11-03

- Make every feature definition a yes/no question in the canonical YAML catalogue.
- Load Rice-virus feature definitions from that YAML catalogue instead of maintaining a duplicate Python feature catalogue.
- Add regression tests requiring every feature definition and screening question to be interrogative.

## [4.0.12] - 2026-11-03

- Display the interrogative feature question in Explorer details instead of the descriptive definition.
- Rewrite spore formation and Rice soil-borne vector questions as direct, abstract-answerable questions.
- Version the revised questions so prior Jev decisions are not reused for changed wording.

## [4.0.11] - 2026-11-03

- Add resumable Jev scoring for every trait-applicable feature question on accepted pathogen–host abstracts.
- Require complete feature decisions before an accepted run can be materialized for the Explorer.
- Merge current completed Rice screens over historical statuses while preserving unsampled historical Rice evidence.
- Keep literature-curated Rice traits and Jev-scored Rice features distinct in the profile and feature picker.
- Preserve existing human bacteria, fungi, and virus profiles while refreshing Rice virus features.

## [4.0.10] - 2026-11-02

- Search rice-virus pairs with ICTV VMR virus names and abbreviations, including their literature names in Jev screening prompts.
- Add 20 literature-listed rice viruses with rice-specific symptom/vector traits and a deterministic rice trait embedding; keep human profiles unchanged.
- Keep three literature names without current ICTV species mappings labeled as literature-associated records.
- Add regression coverage for the 20 rice profiles, host-specific coordinates, and preservation of human profiles.

## [4.0.9] - 2026-11-02

- Use the established confidence label on feature evidence cards as well as pair evidence cards.

## [4.0.8] - 2026-11-02

- Match the pair PubMed cards to the established Explorer design: linked underlined PMIDs, confidence labels, and larger separated cards.


This project follows semantic versioning. Release tags use the form `vMAJOR.MINOR.PATCH`.

## [4.0.7] - 2026-11-02

- Align the pair-evidence sort control with the styled evidence toolbar.

## [4.0.6] - 2026-11-02

- Restore clear card layout for pair-level PubMed evidence, with PMID and confidence in the header and title below.

## [4.0.5] - 2026-11-02

- Keep rice selectable for its saved virus coordinates when it has no feature-value profiles.
- Gray host and pathogen options only when neither feature data nor plotted positions are available.

## [4.0.4] - 2026-11-02

- Allow Pages deployment to proceed after the optional test job is skipped.

## [4.0.3] - 2026-11-02

- Make the fast version check optional for manual deployments.
- Remove browser plot tests and Chromium installation from the deployment workflow.

## [4.0.2] - 2026-11-02

- Display the current semantic version next to the Explorer title.

## [4.0.1] - 2026-11-02

- Restore the V4 species profile's publication-volume and pathogenicity metrics.
- Keep sample-based publication estimates distinct from host-pair search counts.
- Preserve full host and pathogen catalogs while muting options without feature profiles.
- Add data-driven Explorer regression coverage for these behaviors and the 2D plots.
