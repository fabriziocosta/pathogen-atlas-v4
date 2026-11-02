"""Pair-scoped yes/no screening for pathogen disease evidence."""
from __future__ import annotations

from typing import Any, Iterable

import pandas as pd

from pathogen_atlas_v4.decision_identity import decision_cache_key, stable_hash
from pathogen_atlas_v4.prompts import build_pathogen_screen_prompt


def build_pathogen_screen_tasks(sample_matches: pd.DataFrame, articles: pd.DataFrame,
                                pathogens: pd.DataFrame, hosts: pd.DataFrame,
                                prompt_config: dict[str, Any], *, run_id: str,
                                run_config_version: str, model_version: str,
                                profile_version: str = "pathogen-screen-profile-v1") -> list[dict[str, Any]]:
    """Build one auditable yes/no task per sampled PMID and pathogen–host pair."""
    required_matches = {"pathogen_id", "host_id", "pmid", "query_hash", "query_version", "sample_selection"}
    required_articles = {"pmid", "title", "abstract"}
    if not required_matches.issubset(sample_matches.columns):
        raise ValueError(f"Sample matches missing columns: {sorted(required_matches - set(sample_matches.columns))}")
    if not required_articles.issubset(articles.columns):
        raise ValueError(f"Article table missing columns: {sorted(required_articles - set(articles.columns))}")
    pathogen_by_id = pathogens.set_index("pathogen_id")
    host_by_id = hosts.set_index("host_id")
    article_frame = articles.copy()
    article_frame["pmid"] = article_frame.pmid.astype(str)
    article_by_pmid = article_frame.drop_duplicates("pmid", keep="last").set_index("pmid")
    tasks = []
    for match in sample_matches.sort_values(["pathogen_id", "host_id", "pmid"]).to_dict("records"):
        pathogen_id, host_id, pmid = str(match["pathogen_id"]), str(match["host_id"]), str(match["pmid"])
        if pathogen_id not in pathogen_by_id.index or host_id not in host_by_id.index:
            raise ValueError(f"Sample row references unknown pathogen–host pair: {pathogen_id}, {host_id}")
        if pmid not in article_by_pmid.index:
            raise ValueError(f"Sample PMID {pmid} has no fetched article metadata")
        article = article_by_pmid.loc[pmid]
        if isinstance(article, pd.DataFrame):
            article = article.iloc[-1]
        title, abstract = str(article.get("title") or ""), str(article.get("abstract") or "")
        # An abstractless record cannot establish the requested criterion. Leave
        # it unscreened so a pair without another positive remains incomplete.
        if not abstract.strip():
            continue
        pathogen, host = pathogen_by_id.loc[pathogen_id], host_by_id.loc[host_id]
        prompt = build_pathogen_screen_prompt(pathogen_name=str(pathogen.accepted_name), pathogen_id=pathogen_id,
            host_name=str(host.accepted_name), host_id=host_id, article_title=title,
            article_abstract=abstract, prompt_config=prompt_config)
        input_hash = stable_hash({"pathogen_id": pathogen_id, "host_id": host_id,
                                  "query_hash": str(match["query_hash"]), "title": title, "abstract": abstract})
        tasks.append({"scope": "pathogen_host_pair", "decision_kind": "pathogen_screen",
            "pathogen_id": pathogen_id, "host_id": host_id, "relationship_context": "disease",
            "pmid": pmid, "feature_id": "pathogen_screen", "input_hash": input_hash,
            "query_version": str(match["query_version"]), "prompt_version": prompt["prompt_version"],
            "feature_version": "pathogen-screen-system1-v1", "question_version": prompt["question_version"],
            "profile_version": profile_version, "model_version": model_version,
            "run_id": run_id, "run_config_version": run_config_version,
            "query_hash": str(match["query_hash"]), "sample_selection": str(match["sample_selection"]),
            "title_text": title, "abstract_text": abstract, "prompt": prompt})
    # Exact duplicate text under distinct PMIDs can otherwise share one decision
    # key, causing the runner to omit a PMID from the audit ledger. Keep the first
    # stable PMID key reusable and disambiguate later duplicates.
    by_identity = {}
    for task in tasks:
        cache_key = decision_cache_key(scope=task["scope"], pathogen_id=task["pathogen_id"],
            host_id=task["host_id"], relationship_context=task["relationship_context"],
            input_hash=task["input_hash"], feature_id=task["feature_id"],
            query_version=task["query_version"], prompt_version=task["prompt_version"],
            feature_version=task["feature_version"], question_version=task["question_version"],
            profile_version=task["profile_version"], model_version=task["model_version"],
            run_id=task["run_id"], run_config_version=task["run_config_version"])
        by_identity.setdefault(cache_key, []).append(task)
    for duplicates in by_identity.values():
        if len(duplicates) > 1:
            duplicates.sort(key=lambda task: str(task["pmid"]))
            for task in duplicates[1:]:
                task["input_hash"] = stable_hash({"duplicate_source_input_hash": task["input_hash"],
                                                  "pmid": str(task["pmid"])})
    return tasks


def select_pathogen_host_pairs(candidate_pairs: Iterable[tuple[str, str]], sample_matches: pd.DataFrame,
                               ledger: pd.DataFrame) -> pd.DataFrame:
    """Accept a pair when at least one completed sampled abstract supports disease in that host."""
    candidate_keys = list(dict.fromkeys((str(p), str(h)) for p, h in candidate_pairs))
    pair_rows = []
    screen = ledger.copy() if not ledger.empty else pd.DataFrame()
    if not screen.empty and "decision_kind" in screen:
        screen = screen.loc[screen.decision_kind.astype(str).eq("pathogen_screen")]
    samples_by_pair: dict[tuple[str, str], pd.DataFrame] = {}
    if not sample_matches.empty:
        sample_frame = sample_matches.copy()
        sample_frame["pathogen_id"] = sample_frame.pathogen_id.astype(str)
        sample_frame["host_id"] = sample_frame.host_id.astype(str)
        samples_by_pair = {key: group for key, group in sample_frame.groupby(["pathogen_id", "host_id"], sort=False)}
    decisions_by_pair: dict[tuple[str, str], pd.DataFrame] = {}
    if not screen.empty:
        screen["pathogen_id"] = screen.pathogen_id.astype(str)
        screen["host_id"] = screen.host_id.astype(str)
        screen["pmid"] = screen.pmid.astype(str)
        # Preserve the former last-row-wins behavior for retries and resumed ledgers.
        screen = screen.drop_duplicates(["pathogen_id", "host_id", "pmid"], keep="last")
        decisions_by_pair = {key: group for key, group in screen.groupby(["pathogen_id", "host_id"], sort=False)}
    empty_sample = sample_matches.iloc[0:0]
    empty_decisions = screen.iloc[0:0]
    for pathogen_id, host_id in candidate_keys:
        sample = samples_by_pair.get((pathogen_id, host_id), empty_sample)
        decisions = decisions_by_pair.get((pathogen_id, host_id), empty_decisions)
        completed = decisions.loc[decisions.status.astype(str).eq("completed")] if not decisions.empty else decisions
        if not completed.empty and not sample.empty:
            completed = completed.loc[completed.pmid.astype(str).isin(set(sample.pmid.astype(str)))]
        yes = completed.loc[
            completed.get("answer", pd.Series(index=completed.index, dtype="string")).astype(str).str.casefold().eq("yes")
        ] if not completed.empty else completed
        positive_pmids = sorted(set(yes.pmid.astype(str))) if not yes.empty else []
        positive_evidence = ([{"pmid": str(row.pmid), "p_yes": getattr(row, "p_yes", None),
                               "input_hash": str(row.input_hash)}
                             for row in yes.itertuples(index=False)] if not yes.empty else [])
        n_sampled, n_completed = sample.pmid.astype(str).nunique(), completed.pmid.astype(str).nunique() if not completed.empty else 0
        n_with_abstract = int(sample.abstract_available.fillna(False).astype(bool).sum()) if "abstract_available" in sample else n_sampled
        pair_counts = (pd.to_numeric(sample.pair_search_count, errors="coerce").dropna().unique()
                       if "pair_search_count" in sample else [])
        if len(pair_counts) > 1:
            raise ValueError(f"Pair search count is inconsistent for {pathogen_id}, {host_id}")
        n_pair_search_records = int(pair_counts[0]) if len(pair_counts) else 0
        if positive_pmids:
            status = "accepted"
        elif n_sampled == 0:
            status = "no_sampled_records"
        elif n_completed == n_sampled:
            status = "rejected"
        else:
            status = "incomplete"
        pair_rows.append({"pathogen_id": pathogen_id, "host_id": host_id, "relationship_context": "disease",
                          "n_sampled_records": n_sampled, "n_sampled_records_with_abstract": n_with_abstract,
                          "n_completed_screenings": n_completed, "n_pair_search_records": n_pair_search_records,
                          "n_positive_records": len(positive_pmids), "positive_pmids": positive_pmids,
                          "positive_evidence": positive_evidence,
                          "selected_as_pathogen_for_host": bool(positive_pmids), "selection_status": status})
    return pd.DataFrame(pair_rows).sort_values(["pathogen_id", "host_id"]).reset_index(drop=True)
