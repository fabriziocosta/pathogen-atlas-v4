"""Trait-resolved feature screening tasks for V4 pathogen–host runs."""
from __future__ import annotations

from typing import Any

import pandas as pd

from pathogen_atlas_v4.decision_identity import stable_hash
from pathogen_atlas_v4.decision_runner import decision_cache_key_for_task
from pathogen_atlas_v4.prompts import build_feature_question


def build_feature_screen_tasks(sample_matches: pd.DataFrame, articles: pd.DataFrame,
                               pair_selection: pd.DataFrame, pathogens: pd.DataFrame,
                               hosts: pd.DataFrame, profiles: list[dict[str, Any]],
                               prompt_config: dict[str, Any], *, run_id: str,
                               run_config_version: str, model_version: str,
                               profile_version: str,
                               feature_version: str = "v4-feature-screen-noul-v1") -> list[dict[str, Any]]:
    """Create one task for each applicable question and accepted pair abstract.

    Profiles are selected by the pair's host/pathogen traits. Abstracts come from
    the positive pathogen–host evidence already accepted by V4. Questions with pathogen or host scope omit the
    unrelated entity from both prompt state and cache identity; pair-scoped
    questions retain both entities and the disease relationship context.
    """
    if sample_matches.empty or articles.empty or pair_selection.empty:
        return []
    required_matches = {"pathogen_id", "host_id", "pmid"}
    required_articles = {"pmid", "title", "abstract"}
    required_selection = {"pathogen_id", "host_id", "selection_status"}
    for frame, required, label in ((sample_matches, required_matches, "sample matches"),
                                   (articles, required_articles, "articles"),
                                   (pair_selection, required_selection, "pair selection")):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} missing columns: {sorted(missing)}")

    pathogen_by_id = pathogens.assign(pathogen_id=pathogens.pathogen_id.astype(str)).set_index("pathogen_id")
    host_by_id = hosts.assign(host_id=hosts.host_id.astype(str)).set_index("host_id")
    article_by_pmid = articles.assign(pmid=articles.pmid.astype(str)).drop_duplicates("pmid", keep="last").set_index("pmid")
    profile_by_selector = {}
    for item in profiles:
        profile = item.get("resolved_profile", item)
        selectors = profile.get("selectors", {})
        key = (str(selectors.get("pathogen_group", item.get("pathogen_group", ""))),
               str(selectors.get("host_group", item.get("host_group", ""))),
               str(selectors.get("host_taxon_id", item.get("host_taxon_id", ""))))
        profile_by_selector[key] = profile

    accepted = pair_selection.loc[pair_selection.selection_status.astype(str).eq("accepted")]
    accepted_pmids = {
        (str(row.pathogen_id), str(row.host_id), str(pmid))
        for row in accepted.itertuples(index=False)
        for pmid in (getattr(row, "positive_pmids", None)
                     if getattr(row, "positive_pmids", None) is not None else [])
    }
    matches = sample_matches.copy()
    matches["pathogen_id"] = matches.pathogen_id.astype(str)
    matches["host_id"] = matches.host_id.astype(str)
    matches["pmid"] = matches.pmid.astype(str)
    matches = matches.loc[[key in accepted_pmids for key in zip(matches.pathogen_id, matches.host_id, matches.pmid)]]

    tasks: list[dict[str, Any]] = []
    for match in matches.sort_values(["pathogen_id", "host_id", "pmid"]).to_dict("records"):
        pathogen_id, host_id, pmid = match["pathogen_id"], match["host_id"], match["pmid"]
        if pathogen_id not in pathogen_by_id.index or host_id not in host_by_id.index:
            raise ValueError(f"Sample row references unknown pair: {pathogen_id}, {host_id}")
        if pmid not in article_by_pmid.index:
            raise ValueError(f"Sample PMID {pmid} has no fetched article metadata")
        article = article_by_pmid.loc[pmid]
        if isinstance(article, pd.DataFrame):
            article = article.iloc[-1]
        title, abstract = str(article.get("title") or ""), str(article.get("abstract") or "")
        if not abstract.strip():
            continue
        pathogen, host = pathogen_by_id.loc[pathogen_id], host_by_id.loc[host_id]
        taxon_id = str(host.get("taxon_id", ""))
        taxon_id = taxon_id.removeprefix("ncbi:")
        group = str(pathogen.pathogen_group)
        host_group = str(host.host_group)
        profile = profile_by_selector.get((group, host_group, taxon_id))
        if profile is None:
            # Profiles can preserve a namespaced host taxon identifier.
            profile = profile_by_selector.get((group, host_group, str(host.get("taxon_id", ""))))
        if profile is None:
            raise ValueError(f"No resolved feature profile for {group} × {host_group} × {taxon_id}")
        for feature in profile.get("ordered_features", []):
            scope = str(feature["scope"])
            prompt = build_feature_question(pathogen_name=str(pathogen.accepted_name), pathogen_id=pathogen_id,
                host_name=str(host.accepted_name), host_id=host_id, feature=feature,
                article_title=title, article_abstract=abstract, prompt_config=prompt_config)
            context = {"feature_id": str(feature["feature_id"]), "question_version": str(feature.get("question_version", "")),
                       "scope": scope, "article_title": title, "article_abstract": abstract,
                       "pathogen": pathogen_id if scope in {"pathogen", "pathogen_host_pair"} else None,
                       "host": host_id if scope in {"host", "pathogen_host_pair"} else None,
                       "relationship_context": "disease" if scope == "pathogen_host_pair" else None}
            input_hash = stable_hash(context)
            task = {"scope": scope, "decision_kind": "feature_screen", "pathogen_id": pathogen_id,
                    "host_id": host_id if scope in {"host", "pathogen_host_pair"} else None,
                    "relationship_context": "disease" if scope == "pathogen_host_pair" else None,
                    "pmid": pmid, "feature_id": str(feature["feature_id"]), "input_hash": input_hash,
                    "query_version": str(match.get("query_version", "")),
                    "prompt_version": prompt["prompt_version"], "feature_version": feature_version,
                    "question_version": prompt["question_version"],
                    "profile_version": str(profile.get("profile_version") or profile_version),
                    "model_version": model_version, "run_id": run_id,
                    "run_config_version": run_config_version,
                    "title_hash": stable_hash(title), "abstract_hash": stable_hash(abstract),
                    "title_text": title, "abstract_text": abstract,
                    "feature_question": feature["question"], "feature_scope": scope,
                    "feature_label": feature.get("label", feature["feature_id"]),
                    "feature_module": feature.get("module", ""),
                    "feature_category": feature.get("category", ""),
                    "feature_dimension": feature.get("dimension", ""),
                    "applicability_label": feature.get("applies_to_label", ""), "prompt": prompt}
            task["cache_key"] = decision_cache_key_for_task(task)
            tasks.append(task)
    return tasks
