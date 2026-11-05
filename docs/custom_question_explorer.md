# Pathogen Atlas V4: Custom Question Explorer

## Purpose

Build a web app that lets researchers ask their own questions about the literature represented in a user-selected subset of Pathogen Atlas V4. For each selected pathogen, the app evaluates three user-defined questions against the associated PubMed abstracts, counts the positive answers, and displays the counts in an interactive two-dimensional plot with a third question controlling point colour.

The existing Atlas Explorer remains the entry point for selecting records. The new app complements the precomputed views by letting researchers define their own dimensions.

## User journey

1. In the Atlas Explorer, the researcher selects records using the existing search, class, and filtering controls. The current colour selection may help define or preview the selection, but the payload must explicitly identify the selected records; do not infer selection from colour alone.
2. The researcher chooses **Custom questions** (working label) to open the new app.
3. The Explorer passes the selected records and their pathogen, host, and PubMed ID associations to the app.
4. The researcher writes three questions in plain language: one for the X axis, one for the Y axis, and one for colour.
5. A general-purpose LLM reformulates each question into the established structured question format used by the project's Jev workflow. The researcher reviews and confirms each reformulation before running it.
6. The app retrieves abstracts for the supplied PubMed IDs, applies each confirmed question to each abstract through Jev, and records a yes/no result with evidence and processing status.
7. The app counts positive answers for each pathogen on each dimension and plots the results.
8. The researcher can inspect the underlying abstracts, answers, evidence, and failures for any plotted point, then export the results and provenance.

## Data passed from the Explorer

Pass only the selected subset, not the full Atlas dataset. The payload should preserve the association between each pathogen, its host, and its PubMed IDs. Include stable record identifiers where available, plus display names needed by the app.

Suggested conceptual shape (adapt to the Atlas data model):

```json
{
  "schema_version": "1",
  "selection": [
    {
      "record_id": "stable-atlas-record-id",
      "pathogen": {"id": "optional-taxonomy-id", "name": "pathogen name"},
      "host": {"id": "optional-taxonomy-id", "name": "host name"},
      "pmids": ["12345678", "23456789"]
    }
  ]
}
```

Requirements:

- Preserve pathogen–host–PMID relationships; do not flatten the selected PMIDs into one unlabelled list.
- Deduplicate repeated PMIDs for retrieval and abstract evaluation, while retaining all pathogen/host associations that refer to each PMID.
- Validate and report an empty selection, missing PMIDs, invalid identifiers, and payloads that exceed supported limits.
- Do not put the user's OpenRouter key in the Explorer-to-app payload, URL, browser history, or logs.
- Decide whether the handoff uses a short-lived encoded payload, browser storage, or a backend session. Avoid placing a large selection directly in a URL. The app must explain if the selection cannot be transferred and offer a supported import route.

## Question formulation

The user writes one natural-language question for each plot dimension. A general-purpose LLM (for example, the user's preferred model through OpenRouter) converts each into the project's existing Jev-compatible structured question format.

The V4 feature-question system in `docs/feature_question_system.md` and the data model in `docs/data_model.md` are the authoritative conventions. The formulation step must:

- Preserve the researcher's intended meaning and scope.
- Represent each temporary user question using the project's question-record conventions, including a stable run-local identity, wording/claim, scope, question version, applicability, and evidence policy. Keep the Atlas catalogue as the source of existing reusable questions; do not silently add user-authored questions to the released catalogue.
- Ask whether the supplied title and abstract explicitly report the claim for the named pathogen, host, or pathogen–host pair. A positive must be supported by the abstract and the named entity or pair must be identifiable. Do not infer biological facts, unstated mechanisms, or full-text-only results.
- Distinguish supported, unsupported, uncertain, and insufficient-evidence outcomes using the existing decision workflow. An unsupported answer means only that this abstract did not support the claim; it is not evidence of biological absence. Keep missing abstracts and inapplicable questions separate from all scored outcomes.
- Respect the V4 scopes `pathogen`, `host`, and `pathogen_host_pair`. Use the narrowest scope that matches the question and do not reuse a pair-scoped answer for another host.
- Keep the unit of assessment explicit: one question applied to one abstract.
- Ask the LLM to return the structured question and a short readable interpretation for the researcher to review.
- Let the researcher edit the structured question or revise the natural-language version before execution.
- Keep the three questions independent, even when they use the same abstracts.

The researcher confirms all three questions before processing begins. Show the natural-language wording beside the structured Jev version so the mapping can be checked.

## Abstract retrieval and Jev processing

Retrieve the abstracts for the supplied PubMed IDs using the PubMed retrieval facility already used or approved by the project. The transcript refers to this as the “e-search facility”; during implementation, confirm the exact service and endpoint used to fetch records by PMID, and follow the repository's current retrieval convention.

For each distinct PMID:

1. Retrieve and store the abstract text and relevant citation metadata returned by the service.
2. Record missing, retracted, unavailable, or abstract-less records distinctly; do not treat them as negative answers.
3. Apply each confirmed question to the abstract using the established Jev workflow and V4 decision conventions.
4. Store the supported/unsupported/uncertain/insufficient-evidence result, confidence or probability where returned, evidence passage where generated, and model, prompt, question, input hash, and run/version metadata. If the selected Jev mode returns a probability without text, record that mode and link the decision to the PMID and screened-abstract hash, following the V4 data model.
5. Reuse an evaluation for the same abstract and same structured question within a run; do not accidentally count a single PMID multiple times for one pathogen because of duplicate input rows.

Support progress reporting, retry of failed records, and cancellation. Keep successful results if a run is interrupted. Show the user how many abstracts were retrieved, processed, unclear, and failed.

## Counting and plot semantics

For each pathogen–host group, calculate one count per question: the number of distinct associated abstracts classified **supported/Yes** for that question. Use those counts as X, Y, and colour values. The plot is therefore a two-dimensional count plot with a third count encoded by colour. Label these as counts of supporting abstracts; they do not estimate disease prevalence, risk, or biological absence.

Make the counting rule visible in the interface and in exports. In particular:

- Count distinct PMIDs within each pathogen–host group, not duplicate rows.
- Do not count **unsupported**, **uncertain**, **insufficient evidence**, missing abstracts, inapplicable questions, or failed evaluations as supported. Report each separately.
- Display the denominator and the number of unknown/unprocessed records alongside each point, so a low count can be distinguished from incomplete evidence.
- If a PMID is associated with multiple pathogens or hosts, it may contribute to each relevant group, but only once within each group.
- Clarify whether the point represents a pathogen across all hosts or a pathogen–host pair. The initial implementation should use pathogen–host pairs because the input preserves host associations; offer aggregation across hosts only if its counting rule is explicit.
- Provide tooltips and a detail view showing names, X/Y/colour counts, denominators, and the underlying PMID-level results.

Allow the user to change the colour scale and inspect the colour-question wording. Preserve the selection and question set when the plot is refreshed.

## OpenRouter key and privacy

The app uses OpenRouter so the researcher can choose among supported models, including the models already used by the project. Require the researcher to supply their own API key before any LLM processing.

- Keep the key in memory for the active session only by default; do not persist it in local storage, cookies, the URL, source control, telemetry, or application logs.
- Send the key only to the intended OpenRouter request path. Do not send it with PubMed requests or to unrelated services.
- Provide a clear **forget key** action and clear the key when the session ends.
- Explain that abstracts and question text are sent to the selected model provider for processing, and show the selected model before the researcher starts.
- Keep provider-specific request code behind a small adapter so the LLM used for question reformulation or Jev can be changed without changing plot or data logic.
- Never expose the API key in error messages or diagnostic exports.

Confirm the deployment model before implementation. If the app is purely static/client-side, document that the user's browser makes the relevant requests and that the user supplies the key. If any backend is introduced, define secret handling, retention, access control, and logging before sending user keys through it.

## Interface requirements

The app should provide:

- A clear summary of the imported selection, including counts of pathogen–host groups and distinct PMIDs.
- Three labelled question editors: **X axis**, **Y axis**, and **Colour**.
- A review step for the natural-language to Jev-format reformulation.
- OpenRouter key entry, model selection, and a visible key-clearing control.
- A run button with progress, pause/cancel where practical, and retry for failed items.
- A plot with hover details and a selectable point detail view.
- A results table and export that retain the question wording, count definitions, PMID-level outcomes, evidence, and run metadata.
- Clear states for loading, no abstracts, partial completion, errors, and no positive answers.

## Suggested architecture

Keep the implementation modular:

- **Explorer handoff:** serializes the selected records and opens the custom-question app.
- **Selection parser:** validates the payload and builds pathogen–host groups and a unique PMID set.
- **PubMed retrieval adapter:** fetches records for PMID batches and maps results back to all source associations.
- **Question reformulation adapter:** calls the chosen general LLM and returns the established Jev schema.
- **Jev evaluation adapter:** applies a confirmed structured question to one abstract and returns a validated result with evidence.
- **Run store:** tracks inputs, statuses, outcomes, and provenance for the current run.
- **Aggregator:** calculates counts and denominators using the declared rules.
- **Plot and export:** displays and exports the same aggregated results and traceable underlying records.

Keep API calls, limits, retries, concurrency, and error handling configurable. Keep data transformations and count logic independent of any particular LLM provider or chart library.

## Acceptance criteria

The feature is ready for review when:

1. A researcher can select a subset in the Explorer and open the new app with the correct pathogen–host–PMID associations.
2. The app validates the payload and reports missing or invalid data clearly.
3. The researcher can enter three questions, review their Jev-compatible reformulations, and approve them before processing.
4. The app retrieves abstracts by PMID and shows retrieval failures separately from negative answers.
5. The app applies each approved question to each available abstract and retains auditable PMID-level answers and evidence.
6. X and Y show counts of distinct positive abstracts; colour encodes the third question's positive count.
7. Each plotted point shows its denominator and unknown/failed totals, and its detail view reveals the contributing PMID-level results.
8. The researcher can export the plot data and provenance in a documented format.
9. The OpenRouter key is not included in the handoff payload, URL, persistent browser storage, or logs, and the user can clear it.
10. Partial failures and interrupted runs do not silently become No answers or corrupt completed results.

## Decisions to resolve before implementation

- What is the exact Jev structured-question schema and which existing prompt/examples define it?
- Which OpenRouter model should be the default for question reformulation, and which model or models run Jev evaluations?
- What precisely does “colour” mean in the plot: a continuous count gradient, discrete bins, or a thresholded category?
- Should results be grouped by pathogen–host pair only, or should users also be able to aggregate across hosts?
- What are the supported selection size, abstract batch size, concurrency, and rate limits?
- Should results persist beyond the current browser session, and if so, where and for how long?
- What import fallback is needed if a selected payload is too large for the handoff?
- Which export formats are required?
- Does “e-search facility” refer to an existing project wrapper or to PubMed's retrieval interface? Confirm the exact endpoint and existing code path before implementation.
