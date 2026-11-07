# Changelog

## [4.0.118] - 2026-11-07

- Show per-host pathogen counts in the host comparison legend.


## [4.0.117] - 2026-11-07

- Updated Explorer data from dashboard run `dashboard-20261107-001230-338233`.


## [4.0.116] - 2026-11-06

- Updated Explorer data from dashboard run `dashboard-20261106-083003-5ae335`.


## [4.0.115] - 2026-11-06

- Limit the dashboard log panel to the latest two panel-heights of entries.
- Bump dashboard UI version to 1.0.31.

## [4.0.114] - 2026-11-06

- Batch dashboard refreshes driven by streamed logs and append log rows incrementally to keep the health indicator responsive during long runs.
- Handle browser disconnects without printing server tracebacks; bump dashboard UI version to 1.0.30.

## [4.0.113] - 2026-11-06

- Start the PubMed health check immediately and time out checks that do not return, so the dashboard cannot remain stuck on “Checking…”.
- Bump the dashboard UI version to 1.0.29.

## [4.0.112] - 2026-11-06

- Calculate the estimated pathogenicity fraction against full host–pathogen PubMed query results, with the 200-record sample remaining a separate count.
- Restore full pair-query totals for historical human profiles and label them correctly in the Explorer.

## [4.0.111] - 2026-11-06

- Show host–pathogen PubMed search results separately from the pathogen's organism-wide publication total.
- Clarify how sampled records, pair-search results, and the organism-wide pathogenicity denominator relate.

## [4.0.110] - 2026-11-06

- Preserve the selected pathogen's publication estimate when switching primary and comparison hosts.
- Improve comparison plot labels and controls, and refine dashboard host-group selection, phase timing, and run-header layout.
- Bump the dashboard UI version to 1.0.28.

## [4.0.109] - 2026-11-06

- Rebuilt the Pages bundle from canonical Explorer files before validation and deployment.
- Made dashboard retries create a fresh, versioned release instead of rerunning the same failed commit.
- Published the completed Rice cohort with matching source and Pages data.

## [4.0.108] - 2026-11-06

- Updated Explorer data from dashboard run `dashboard-20261105-181402-14df69`.


## [4.0.107] - 2026-11-05

- Restore Rice virus points when switching from an embedding that is unavailable for that host and group.
- Keep the map explanation synchronized when the current selection has no visible points.

## [4.0.106] - 2026-11-05

- Fix dashboard publishing so commits include the exact GitHub Pages bundle and source payloads.
- Validate source-to-dist synchronization and the release bundle before confirmation and before push.
- Remove stale taxonomy shards during dashboard builds.
- Publish completed Soybean × Fungi data.

## [4.0.105] - 2026-11-05

- Updated Explorer data from dashboard run `dashboard-20261105-074510-81ad48`.


## [4.0.104] - 2026-11-05

- Align the Custom questions button with the secondary axis selectors and place status text beneath it.

## [4.0.103] - 2026-11-05

- Count only nonzero Color by values in the highlighted pathogen coverage metric.

## [4.0.102] - 2026-11-05

- Omit unavailable or irrelevant features from the secondary axis trait selectors.

## [4.0.101] - 2026-11-05

- Align and add breathing room to the secondary comparison controls, including the Custom questions action.

## [4.0.100] - 2026-11-05

- Attach resolved Wikipedia and GBIF URLs to every pathogen species profile so links appear for any host–pathogen pair.

## [4.0.99] - 2026-11-05

- Add release checks for PubMed documents, PMIDs, feature references, taxonomy shards, and 2D coordinates before and after Pages deployment.
- Fetch evidence with a release-specific cache key so updated feature summaries cannot pair with stale document data.
- Wait for an active Pages deployment before publishing another release.

## [4.0.98] - 2026-11-05

- Updated Explorer data from dashboard run `dashboard-20261104-214740-d2ff2f`.


## [4.0.94] - 2026-11-04

- Updated Explorer data from dashboard run `dashboard-20261104-181846-b4fdf7`.


## [4.0.85] - 2026-11-04

- Updated Explorer data from dashboard run `dashboard-20261104-171147-23ffe2`.


## [4.0.84] - 2026-11-04

- Keep a stable embedding available when another embedding method collapses for a small host–pathogen cohort.
- Publish Soybean × Bacteria with a UMAP view built from the 14 complete scored pair profiles.

## [4.0.83] - 2026-11-04

- Updated Explorer data from dashboard run `dashboard-20261104-142311-bc7692`.


## [4.0.82] - 2026-11-04

- Report per-phase elapsed time, split phase metrics across two lines, and add resumable stop controls to the dashboard.

## [4.0.81] - 2026-11-04

- Tighten host-menu name columns so the longest common name starts two characters after its group heading.

## [4.0.80] - 2026-11-04

- Move the fixed separator closer to the checkmark by tightening the common-name column.

## [4.0.79] - 2026-11-04

- Remove redundant abstract and pair preambles from displayed Explorer feature questions.

## [4.0.78] - 2026-11-04

- Simplify the displayed HOST ENTRY ROUTE question while retaining the versioned screening prompt.

## [4.0.77] - 2026-11-04

- Default to Toxoplasma gondii for the Human × Parasites Explorer view.

## [4.0.76] - 2026-11-04

- Keep host menu columns anchored at the checkmark gutter and fit scientific names within the visible menu.

## [4.0.75] - 2026-11-04

- Show pathogen evidence rows when Explorer evidence records use the pathogen_id schema.

## [4.0.74] - 2026-11-04

- Add the Human × Parasites cohort, evidence, feature scores, and coordinates to the Explorer.
- Speed up Explorer assembly by indexing decisions by PMID and reusing exact-abstract cache matches.


## [4.0.73] - 2026-11-04

- Keep dashboard screening concurrency bounded and make the oldest PMID subset adapt to small sample sizes.

## [4.0.72] - 2026-11-04

- Add a local V4 pipeline dashboard with per-run configuration, phase progress, persistent logs, pause/resume checkpoints, cross-run semantic decision reuse, and explicit Explorer publishing.

## [4.0.71] - 2026-11-04

- Add gray inline species-name completion with Tab acceptance and Enter search; tighten the taxonomy-source spacing.

## [4.0.70] - 2026-11-04

- Align coverage metric values at the top of their cells and vertically center taxonomy provenance beneath the metrics.

## [4.0.69] - 2026-11-04

- Bottom-align the Wikipedia and GBIF links with the species name.

## [4.0.68] - 2026-11-04

- Show the selected pathogen and dotted crosshair in the secondary 2D plot, darken feature-score columns for standard-deviation mark contrast, close color settings on outside click, and use Pythium insidiosum as the human Oomycetes default.

## [4.0.67] - 2026-11-04

- Match the feature-score chart size and colors to V2, keep its legend centered, and color Atlas in the header dark green.

## [4.0.66] - 2026-11-04

- Strengthen only the bacterium silhouette outline in the Explorer masthead.

## [4.0.65] - 2026-11-04

- Add the full publication reference under the Explorer authors and increase the space above the Explorer panels.

## [4.0.64] - 2026-11-04

- Inset the mean plus one standard deviation dash to reveal its rounded ends against the score bars.

## [4.0.63] - 2026-11-04

- Make the mean plus one standard deviation dash thinner with rounded ends.

## [4.0.62] - 2026-11-04

- Thin the mean plus one standard deviation dash and gray the mean marker outline.

## [4.0.61] - 2026-11-04

- Render mean plus one standard deviation as a white dash and add striped space between feature columns.

## [4.0.60] - 2026-11-04

- Show feature means and mean plus one standard deviation as white circular markers.

## [4.0.59] - 2026-11-04

- Use the accepted-profile mean as the feature reference and recolor score bars blue-cyan.

## [4.0.58] - 2026-11-04

- Use a pale yellow for the feature cohort median reference.

## [4.0.57] - 2026-11-04

- Make feature score bars translucent so the host-group median remains visible behind them.

## [4.0.56] - 2026-11-04

- Brighten and square the feature median reference bars behind the pathogen scores.

## [4.0.55] - 2026-11-04

- Show host × pathogen-group feature medians behind pathogen feature scores.

## [4.0.54] - 2026-11-04

- Add selectable IMSTE and UMAP pathogen maps and recompute human bacteria, virus, and fungi embeddings.
- Rename Explorer Color settings to Settings.

## [4.0.53] - 2026-11-03

- Default the human virus profile to Dengue virus instead of rice dwarf polerovirus.

## [4.0.52] - 2026-11-03

- Set YlGnBu as the default feature colormap.

## [4.0.51] - 2026-11-03

- Give the foreground bacterium opaque white interiors and enlarge it to extend beyond the sphere.

## [4.0.50] - 2026-11-03

- Overlay the bacterium as a separate image over the SVG globe and thicken the outer circle.

## [4.0.49] - 2026-11-03

- Build the globe from precise SVG circle geometry and layer the bacterium illustration over it.

## [4.0.48] - 2026-11-03

- Align the globe center vertically with the lowercase “o” in the title.

## [4.0.47] - 2026-11-03

- Refine the globe to a perfectly circular, symmetric sphere behind the bacterium.

## [4.0.46] - 2026-11-03

- Use the supplied hand-drawn bacterium illustration inside the geometric globe logo.

## [4.0.45] - 2026-11-03

- Align the logo center with the visual center of the Pathogen Atlas wordmark.

## [4.0.44] - 2026-11-03

- Use a single ink color and contrast the geometric globe with a hand-drawn bacterium mark.

## [4.0.43] - 2026-11-03

- Replace the illustrated Explorer logo with a minimal two-color atlas and pathogen mark.

## [4.0.42] - 2026-11-03

- Use the shared unboxed disclosure row styling for Feature scores, matching its title and expand icon alignment to the other sections.

## [4.0.41] - 2026-11-03

- Match the Feature scores disclosure title size to the other profile section titles.

## [4.0.40] - 2026-11-03

- Place the Feature scores chart between Supporting PubMed records and Features.

## [4.0.39] - 2026-11-03

- Make the feature score chart collapsible.

## [4.0.38] - 2026-11-03

- Increase feature score bar width and reduce spacing between features.

## [4.0.37] - 2026-11-03

- Reduce and center the square feature score chart.

## [4.0.36] - 2026-11-03

- Add a square feature score overview to each species profile.

## [4.0.35] - 2026-11-03

- Match the taxonomy-source line to coverage labels on portrait phones and prevent mobile text autosizing from changing their relative size.

## [4.0.34] - 2026-11-03

- Improve portrait phone layout, align the coverage metrics in a 3 × 2 grid, and refine header typography and value colors.

## [4.0.33] - 2026-11-03

- Report observed Jev costs and forecast remaining OpenRouter credits during feature screening.

## [4.0.32] - 2026-11-03

- Show cached and current-run decision totals in screening progress output.

## [4.0.31] - 2026-11-03

- Report feature-scoring progress every 5,000 completed decisions with throughput and approximate ETA.

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
