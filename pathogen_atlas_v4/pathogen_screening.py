"""Pair-scoped yes/no screening for pathogen disease evidence."""
from __future__ import annotations

from typing import Any, Iterable

import pandas as pd

from pathogen_atlas_v4.decision_identity import stable_hash
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
            "feature_version": "pathogen-screen-v2", "question_version": prompt["question_version"],
            "profile_version": profile_version, "model_version": model_version,
            "run_id": run_id, "run_config_version": run_config_version,
            "query_hash": str(match["query_hash"]), "sample_selection": str(match["sample_selection"]),
            "title_text": title, "abstract_text": abstract, "prompt": prompt})
    return tasks


def select_pathogen_host_pairs(candidate_pairs: Iterable[tuple[str, str]], sample_matches: pd.DataFrame,
                               ledger: pd.DataFrame) -> pd.DataFrame:
    """Accept a pair when at least one completed sampled abstract supports disease in that host."""
    pair_rows = []
    screen = ledger.copy() if not ledger.empty else pd.DataFrame()
    if not screen.empty and "decision_kind" in screen:
        screen = screen.loc[screen.decision_kind.astype(str).eq("pathogen_screen")]
    for pathogen_id, host_id in dict.fromkeys((str(p), str(h)) for p, h in candidate_pairs):
        sample = sample_matches.loc[
            sample_matches.pathogen_id.astype(str).eq(pathogen_id)
            & sample_matches.host_id.astype(str).eq(host_id)
        ]
        decisions = screen.loc[
            screen.pathogen_id.astype(str).eq(pathogen_id)
            & screen.host_id.astype(str).eq(host_id)
        ] if not screen.empty else pd.DataFrame()
        if not decisions.empty:
            decisions = decisions.drop_duplicates("pmid", keep="last")
        completed = decisions.loc[decisions.status.astype(str).eq("completed")] if not decisions.empty else decisions
        if not completed.empty and not sample.empty:
            completed = completed.loc[completed.pmid.astype(str).isin(set(sample.pmid.astype(str)))]
        yes = completed.loc[
            completed.get("answer", pd.Series(index=completed.index, dtype="string")).astype(str).str.casefold().eq("yes")
            & completed.get("evidence_span", pd.Series(index=completed.index, dtype="string")).fillna("").astype(str).str.strip().ne("")
        ] if not completed.empty else completed
        positive_pmids = sorted(set(yes.pmid.astype(str))) if not yes.empty else []
        positive_evidence = ([{"pmid": str(row.pmid), "evidence_span": str(row.evidence_span),
                               "input_hash": str(row.input_hash)}
                             for row in yes.itertuples(index=False)] if not yes.empty else [])
        n_sampled, n_completed = sample.pmid.astype(str).nunique(), completed.pmid.astype(str).nunique() if not completed.empty else 0
        n_with_abstract = int(sample.abstract_available.fillna(False).astype(bool).sum()) if "abstract_available" in sample else n_sampled
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
                          "n_completed_screenings": n_completed,
                          "n_positive_records": len(positive_pmids), "positive_pmids": positive_pmids,
                          "positive_evidence": positive_evidence,
                          "selected_as_pathogen_for_host": bool(positive_pmids), "selection_status": status})
    return pd.DataFrame(pair_rows).sort_values(["pathogen_id", "host_id"]).reset_index(drop=True)
