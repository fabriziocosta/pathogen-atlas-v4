# Changelog

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
