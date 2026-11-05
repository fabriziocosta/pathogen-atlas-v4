# Pathogen Atlas V4: Custom Question Explorer

## Purpose

Custom Question Explorer lets researchers evaluate three questions of their own against the PubMed abstracts that already support disease relationships in Pathogen Atlas V4. It plots the number of abstracts supporting each question for each pathogen–host pair. It is an exploratory analysis of literature evidence, not a measure of disease prevalence, risk, or biological absence.

The app is a static, browser-only companion to the Atlas Explorer. The browser calls PubMed and OpenRouter directly. The researcher supplies their own OpenRouter API key.

## User journey

1. In the Atlas Explorer, the researcher sets the host, pathogen group, search, minimum-support filter, and **Color by** selection.
2. The researcher opens **Custom questions**. The selection is snapshotted from the currently visible, non-gray map points. Hidden points and gray points are excluded. In the current feature-color view, gray means zero or unscored; a selected point must also have an accepted disease relationship and confirmed pathogenicity PMIDs.
3. The app shows the included pathogen–host pairs, confirmed PMID count, and excluded gray/hidden counts. The researcher reviews the summary and continues or returns to the Explorer to adjust it.
4. The researcher enters one natural-language question for each of **X axis**, **Y axis**, and **Colour**, supplies an OpenRouter key, and selects a generative model for reformulation.
5. The app reformulates the questions into the V4 Jev-compatible shape. The researcher reviews and may edit the structured versions, then approves all three before evaluation.
6. The app retrieves the confirmed PMIDs’ abstracts through the V4 PubMed EFetch convention and evaluates each approved question against each available abstract using Jev.
7. The app reports progress and partial results, then plots one point per pathogen–host pair. The researcher can inspect PMID-level results and download CSV and JSON exports.

## Explorer selection and handoff

The handoff represents the rendered set, not a second selection UI. At the moment the researcher opens Custom Questions, snapshot the currently displayed non-gray points after the active Explorer controls have been applied. The selection is the intersection of:

- the currently selected host and pathogen group(s),
- the active search and minimum-support filter,
- points rendered with a non-gray color under the current **Color by** setting, and
- pathogen–host pairs whose disease relationship is accepted and has confirmed positive pathogen PMIDs.

Do not include gray points, hidden points, unaccepted or unscreened pairs, search-result candidates without confirmed disease evidence, or PubMed search-result totals. Do not infer pathogen PMIDs from a point’s color value. Use the accepted PMID-level disease evidence linked to that pair. Show the researcher which points and PMIDs will be included before processing starts.

Pass a versioned JSON selection snapshot. Keep the payload out of the URL and browser history. For a normal handoff, write it to same-origin `sessionStorage`, navigate the same tab to the app, read it once, and remove it from storage. If it exceeds the supported handoff size, offer a JSON download from the Explorer and a matching JSON import in the app. Never silently truncate a selection. The app validates imported data and explains empty selections, invalid IDs, missing PMIDs, and malformed or unsupported schema versions.

Conceptual payload:

```json
{
  "schema_version": "1",
  "snapshot": {
    "host_ids": ["ncbi:9606"],
    "pathogen_groups": ["bacteria"],
    "search": "",
    "minimum_supporting_pmids": 0,
    "color_by": "feature-id"
  },
  "selection": [
    {
      "record_id": "stable-atlas-pair-id",
      "pathogen": {"id": "ncbi:12345", "name": "Pathogen species"},
      "host": {"id": "ncbi:9606", "name": "Homo sapiens"},
      "pmids": ["12345678", "23456789"]
    }
  ]
}
```

Preserve pathogen–host–PMID relationships. Deduplicate PMID retrieval globally, then map each abstract back to every selected association. A PMID may contribute once to each distinct selected pathogen–host pair it supports. Evaluate once per distinct abstract, question, and applicable entity scope within a run; pair-scoped evaluations remain distinct for different pairs. Do not impose a fixed PMID cap or silently sample the confirmed set.

## Question formulation and Jev evaluation

Use `docs/feature_question_system.md`, `docs/data_model.md`, `pathogen_atlas_v4/prompts.py`, and `pathogen_atlas_v4/jev_pathogen_screen.py` as the source of truth for question semantics and Jev request shape. A custom question is run-local and must not be inserted into the released feature catalogue.

Each question must retain a run-local identity, original wording, claim, scope (`pathogen`, `host`, or `pathogen_host_pair`), version, applicability, and evidence policy. Ask whether the supplied title and abstract explicitly report the claim for the named entity or pair. Use the narrowest correct scope. Require the researcher to approve the reformulation before any abstract evaluation begins; show original and structured wording together.

For each question and abstract, construct a scope-aware Jev state containing the question identity/version, structured question, scope, article title and abstract, and the applicable pathogen and/or host identity. Pair-scoped questions include both entities and disease relationship context. Send the Jev typed question with the existing `type`, `instructions`, and `criteria` fields. Keep X, Y, and Colour independent, even when they share an abstract.

Use `typesafe/jev-1.14` as the default Jev evaluation model. The researcher may select/edit the OpenRouter generative model used for reformulation. Preserve both model IDs in the run provenance. For Jev’s yes-probability output, classify `p_yes >= 0.7` as **supported** and lower values as **unsupported**, matching the V4 threshold. Retain `p_yes`; unsupported means only that this abstract did not support the claim. This Jev mode does not generate evidence text: link a result to its PMID and screened-abstract hash rather than manufacturing a quotation. Missing abstracts, not-yet-processed items, provider errors, and failed evaluations remain separate incomplete statuses, never negatives.

Before evaluation, explain that question text, abstracts, and associated pathogen/host names will be sent to OpenRouter for model processing. Display the selected rewrite and Jev models. Send the OpenRouter key only in the Authorization header of OpenRouter requests; never include it in the handoff, URL, PubMed requests, logs, telemetry, or exports. Keep it in memory for the active session only, provide a **Forget key** action, and clear it when the session ends. Do not store run results beyond the active browser session; users can save them through export.

## PubMed retrieval and run behavior

Retrieve the supplied confirmed PMIDs using PubMed EFetch from the existing V4 retrieval convention, not ESearch. Request records in batches of up to 500 PMIDs, using POST for large ID lists, and retain returned citation metadata and abstract text. Use EFetch’s PMID mapping to reconnect each result to its selected associations. Preserve the user’s confirmed PMID list if a record is missing or has no abstract; report those cases separately.

Keep provider adapters, validation, counting, and charting separate. Process all selected PMIDs with bounded concurrency; use configurable retry and exponential backoff for transient network errors, HTTP 429, and HTTP 5xx. Do not retry permanent request or validation errors automatically. Support cancellation and retrying failed items without discarding successful results. Keep completed work in memory for the active session. Show retrieved, evaluated, supported, unsupported, missing-abstract, pending, and failed totals, including partial completion after cancellation.

Run metadata must include a unique run ID and timestamp, selection snapshot, question wording and approved Jev structures, question version/identity, model IDs, prompt versions, retrieval and evaluation statuses, PMID, title/abstract hashes, `p_yes`, threshold, and app/schema version. Do not include the API key.

## Counting and plot semantics

Render one point per pathogen–host pair. For each of X, Y, and Colour, count distinct associated PMIDs classified supported for that question. Do not count a duplicate PMID more than once within a pair. The third count uses a continuous colour gradient over raw counts; it is not a normalized rate or a thresholded category. Host aggregation is out of scope for the initial version.

For every point and dimension, show the supported count, evaluated denominator, and incomplete count. The denominator is the distinct associated PMIDs with an abstract successfully evaluated for that question. Report unsupported separately from missing, pending, and failed items. A zero supported count is meaningful only when the denominator is greater than zero; display **Not evaluated** when it is zero. Failures and unavailable abstracts never count as supported or unsupported.

Tooltips and the point detail view show pathogen and host names, all three question wordings, supported counts, denominators, incomplete totals, and PMID-level outcomes. The detail view links to PubMed records and exposes the abstract/hash used for each evaluation. Keep the counting rule visible beside the plot and in both exports.

## Interface and exports

The app provides:

- Selection summary with included pairs, unique PMIDs, and excluded gray/hidden points.
- Three labelled question editors, a model selector/edit field, API-key field, and **Forget key** control.
- A review/approval step displaying natural-language and Jev-structured questions together.
- Run, progress, cancel, and failed-item retry controls.
- A continuous-colour scatter plot, point details, and a PMID-level results table.
- Clear states for empty or invalid selection, missing abstracts, no evaluable records, no supported answers, partial completion, cancellation, and provider errors.
- A grouped-results CSV and a JSON provenance bundle. CSV contains one row per pathogen–host pair with the X/Y/Colour supported counts, denominators, incomplete totals, question labels, and run ID. JSON contains the selection snapshot, questions, every PMID-level retrieval/evaluation outcome, hashes, metadata, and count definitions.

## Architecture and operational defaults

Keep the app static and browser-only. Use same-origin `sessionStorage` only for the transient Explorer-to-app handoff, then remove the payload after import. Keep the active run and key in memory; there is no backend, account, remote run store, or automatic persistence.

Separate the implementation into an Explorer handoff, selection parser, PubMed EFetch adapter, question reformulation adapter, Jev evaluation adapter, in-memory run store, aggregation layer, and plot/export layer. Use the V4 EFetch batch size (up to 500 IDs), a configurable small Jev concurrency default of 2, and bounded exponential retry for transient errors. Ensure browser-origin requests to both external services work in the deployed site; if browser CORS or service policy prevents a required direct call, do not silently proxy through another service or expose the key—resolve the architecture constraint before release.

The only required user secret is their own OpenRouter API key. No Atlas backend credentials are used. Model IDs and concurrency/retry settings are configurable; Jev defaults to `typesafe/jev-1.14`, while the rewrite model is selected/edited by the researcher.

## Acceptance criteria

1. Opening Custom Questions snapshots exactly the current visible, non-gray Explorer points that have accepted disease relationships and confirmed positive pathogen PMIDs.
2. The handoff preserves pair-to-PMID associations, excludes hidden/gray/unaccepted records, deduplicates retrieval IDs, and offers JSON download/import without truncation.
3. The app retrieves every confirmed PMID through EFetch and reports missing/abstract-less records separately.
4. Researchers can write, review, edit, and approve three independent Jev-compatible questions before evaluation.
5. Jev evaluation uses the existing `noul` typed question shape and 0.7 threshold, retaining each probability and the PMID/abstract-hash provenance.
6. The plot has one point per pathogen–host pair and displays distinct supported PMID counts for X/Y/Colour, with denominators and incomplete totals.
7. Duplicate PMIDs count at most once within each pair; the same PMID can count in each selected pair it supports.
8. Cancellation and failures preserve completed in-session results; retries do not convert unavailable or failed records into negatives.
9. CSV and JSON exports reproduce the grouped plot data and preserve question and PMID-level provenance.
10. The OpenRouter key never enters the URL, Explorer handoff, persistent browser storage, PubMed requests, logs, telemetry, or exports, and can be cleared by the researcher.
11. Browser-origin requests to PubMed and OpenRouter work from the deployed static site, or the release remains blocked pending resolution of that architecture constraint.
