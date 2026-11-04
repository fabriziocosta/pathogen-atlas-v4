"""Build the full atlas from taxonomy and completed host-pair releases."""

from __future__ import annotations

import argparse
import hashlib
import gzip
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

HOST_LABELS = {
    "ncbi:9606": "Human · Homo sapiens",
    "ncbi:9913": "Cattle · Bos taurus",
    "ncbi:9940": "Sheep · Ovis aries",
    "ncbi:9031": "Chicken · Gallus gallus",
    "ncbi:9823": "Pig · Sus scrofa",
    "ncbi:4565": "Wheat · Triticum aestivum",
    "ncbi:3847": "Soybean · Glycine max",
    "ncbi:4577": "Maize · Zea mays",
    "ncbi:4530": "Rice · Oryza sativa",
    "ncbi:8049": "Cod · Gadus morhua",
    "ncbi:8030": "Salmon · Salmo salar",
}
HOST_GROUPS = {
    "ncbi:9606": "Animals", "ncbi:9913": "Animals", "ncbi:9940": "Animals",
    "ncbi:9031": "Animals", "ncbi:9823": "Animals", "ncbi:4565": "Plants",
    "ncbi:3847": "Plants", "ncbi:4577": "Plants", "ncbi:4530": "Plants",
    "ncbi:8049": "Animals", "ncbi:8030": "Animals",
}
TAXONOMY_REGISTRY_PATH = Path(__file__).resolve().parents[3] / "configs" / "pathogen_taxonomies.yaml"
with TAXONOMY_REGISTRY_PATH.open(encoding="utf-8") as registry_file:
    _TAXONOMY_REGISTRY = yaml.safe_load(registry_file) or {}
TAXONOMY_SOURCES = {
    group: {
        "name": entry["source"]["name"],
        "url": entry["source"]["url"],
        "date": entry["source"]["date"],
        "date_label": entry["source"]["date_label"],
    }
    for group, entry in _TAXONOMY_REGISTRY.get("pathogen_groups", {}).items()
}

def _json_safe(value: Any) -> Any:
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if hasattr(value, "tolist"):
        try:
            return _json_safe(value.tolist())
        except (ValueError, AttributeError):
            pass
    if pd.api.types.is_scalar(value) and pd.isna(value):
        return None
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, AttributeError):
            pass
    return value


def _taxonomy(path: Any) -> dict[str, str]:
    return {part.split("=", 1)[0]: part.split("=", 1)[1]
            for part in str(path or "").split(";") if "=" in part}


def _taxonomy_provenance(groups: list[str], species_rows: list[dict[str, Any]],
                         ncbi_retrieved_date: str) -> dict[str, dict[str, Any]]:
    provenance = {}
    for group in groups:
        source = TAXONOMY_SOURCES.get(group)
        if not source:
            continue
        item = dict(source)
        if item["date"] == "ncbi_manifest":
            item["date"] = ncbi_retrieved_date
        item["taxa_count"] = sum(row["organism_group"] == group for row in species_rows)
        provenance[group] = item
    return provenance


def _read_csv_index(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    frame = pd.read_csv(path, keep_default_na=False, dtype=str)
    return {str(row["pathogen_id"]): row for row in frame.to_dict("records")}


def _copy_static_site(site_root: Path, dist_root: Path) -> None:
    dist_root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(site_root / "index.html", dist_root / "index.html")
    for name in ("explorer-data.json.gz", "evidence.json.gz", "feature-catalog.json"):
        shutil.copy2(site_root / name, dist_root / name)
    if (site_root / "assets").exists():
        shutil.copytree(site_root / "assets", dist_root / "assets", dirs_exist_ok=True)
    if (site_root / "taxa").exists():
        shutil.copytree(site_root / "taxa", dist_root / "taxa", dirs_exist_ok=True)


def build_site_data(run_dir: Path, staging_root: Path, site_root: Path,
                    *, dist_root: Path | None = None) -> dict[str, Any]:
    """Materialize the full taxonomy atlas and overlay completed host-pair releases."""
    run = Path(run_dir)
    staging = Path(staging_root)
    site = Path(site_root)
    previous_payload_path = (Path(dist_root) / "explorer-data.json.gz") if dist_root is not None else (site / "explorer-data.json.gz")
    try:
        with gzip.open(previous_payload_path, "rt", encoding="utf-8") as previous_file:
            previous_payload = json.load(previous_file)
    except (OSError, json.JSONDecodeError):
        previous_payload = {}
    previous_evidence_path = (Path(dist_root) / "evidence.json.gz") if dist_root is not None else (site / "evidence.json.gz")
    try:
        with gzip.open(previous_evidence_path, "rt", encoding="utf-8") as previous_file:
            previous_evidence = json.load(previous_file)
    except (OSError, json.JSONDecodeError):
        previous_evidence = []
    prior_embedding_runs = previous_payload.get("embedding_runs", {})
    project_root = Path(__file__).parents[3]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from pathogen_atlas_v4.pathogen_taxonomy import load_ncbi_pathogen_additions

    resolved = json.loads((run / "resolved_config.json").read_text(encoding="utf-8"))
    pathogens = pd.read_parquet(staging / "taxonomy/pathogens.parquet")
    pathogens["pathogen_id"] = pathogens.pathogen_id.astype(str)
    pathogens = pathogens.loc[pathogens["rank"].astype(str).str.casefold().eq("species")].copy()
    ncbi_manifest = json.loads((project_root / "data/source_taxonomy/ncbi/source_manifest.json").read_text())
    release = f"ncbi-taxdump-{ncbi_manifest['retrieved_utc'][:10]}-sha256-{ncbi_manifest['archive_sha256'][:12]}"
    ncbi_retrieved_date = str(ncbi_manifest["retrieved_utc"])[:10]
    additions, added_aliases, _ = load_ncbi_pathogen_additions(
        project_root / "data/source_taxonomy/ncbi",
        ["oomycetes", "parasites"], release=release, source_version=release,
    )
    additions = additions.loc[additions["rank"].astype(str).str.casefold().eq("species")]
    pathogens = pd.concat([pathogens, additions], ignore_index=True, sort=False)
    # The resolved run taxonomy is authoritative for any taxa it includes.
    # In particular, host/group cohort selectors can reclassify taxa from the
    # staged catalogue (for example, LPSN bacterial species in the Parasites
    # cohort). Keeping the staged group here makes the Explorer builder look
    # for a feature profile under the wrong selector and drops the new cohort.
    run_taxonomy_path = run / "taxonomy/pathogens.parquet"
    if run_taxonomy_path.exists():
        run_pathogens = pd.read_parquet(run_taxonomy_path)
        run_pathogens["pathogen_id"] = run_pathogens.pathogen_id.astype(str)
        run_pathogens = run_pathogens.loc[run_pathogens["rank"].astype(str).str.casefold().eq("species")]
        replaced_ids = set(run_pathogens.pathogen_id)
        pathogens = pd.concat([
            pathogens.loc[~pathogens.pathogen_id.astype(str).isin(replaced_ids)],
            run_pathogens,
        ], ignore_index=True, sort=False)
    if pathogens.pathogen_id.duplicated().any():
        duplicate_ids = pathogens.loc[pathogens.pathogen_id.duplicated(False), "pathogen_id"].head().tolist()
        raise ValueError(f"Atlas taxonomy contains duplicate pathogen IDs: {duplicate_ids}")
    pairs = pd.read_parquet(run / "pathogen_host_selection.parquet")
    pairs["pathogen_id"] = pairs.pathogen_id.astype(str)
    pairs["host_id"] = pairs.host_id.astype(str)
    if pairs.duplicated(["pathogen_id", "host_id"]).any():
        raise ValueError("Pair selection must contain one row per pathogen and host")
    selected_pathogen_ids = set(pathogens.pathogen_id.astype(str))
    selected_host_ids = set(pd.read_parquet(staging / "taxonomy/hosts.parquet").host_id.astype(str))
    pairs = pairs.loc[pairs.pathogen_id.isin(selected_pathogen_ids)
                      & pairs.host_id.isin(selected_host_ids)].copy()

    pii_path = run / "release/pathogenic_interest_index.parquet"
    pii = pd.read_parquet(pii_path) if pii_path.exists() else pd.DataFrame()
    if not pii.empty:
        pii["pathogen_id"] = pii.pathogen_id.astype(str)
        pii["host_id"] = pii.host_id.astype(str)
    articles = pd.read_parquet(staging / "literature/articles.parquet")
    articles["pmid"] = articles.pmid.astype(str)
    # Add metadata fetched for this run so newly assessed PMIDs retain titles
    # and publication details in Explorer evidence cards.
    source_config = resolved.get("source_config", {})
    acquisition = source_config.get("acquisition", {})
    from pathogen_atlas_v4.pubmed_client import literature_cache_key
    cache_key = literature_cache_key(
        query_version=acquisition.get("query_strategy_version", ""),
        pathogen_literature_query_version=acquisition.get("pathogen_literature_query_version", ""),
        taxonomy_snapshot=staging.name,
        host_expansion_hash=acquisition.get("host_expansion_hash", ""),
    )
    run_articles_path = project_root / "data/v4/cache/literature" / cache_key / "tables/articles.parquet"
    if run_articles_path.exists():
        run_articles = pd.read_parquet(run_articles_path)
        if not run_articles.empty and "pmid" in run_articles:
            run_articles["pmid"] = run_articles.pmid.astype(str)
            articles = pd.concat([articles, run_articles], ignore_index=True, sort=False)
    article_by_pmid = articles.drop_duplicates("pmid", keep="last").set_index("pmid")
    evidence = pd.read_parquet(run / "evidence/pathogen_host_evidence.parquet")
    evidence["pathogen_id"] = evidence.pathogen_id.astype(str)
    evidence["host_id"] = evidence.host_id.astype(str)
    evidence["pmid"] = evidence.pmid.astype(str)
    evidence = evidence.loc[evidence.decision.astype(str).eq("supported")].copy()
    evidence_by_pair = {key: group for key, group in evidence.groupby(["pathogen_id", "host_id"], sort=False)}
    host_taxonomy = pd.read_parquet(staging / "taxonomy/hosts.parquet")
    host_taxonomy["host_id"] = host_taxonomy.host_id.astype(str)
    host_taxonomy_by_id = host_taxonomy.drop_duplicates("host_id", keep="last").set_index("host_id")
    resolved_feature_profiles = {}
    for profile_path in sorted((run / "feature_profiles").glob("*.json")):
        profile = json.loads(profile_path.read_text(encoding="utf-8"))
        selector = profile.get("selectors", {})
        resolved_feature_profiles[(str(selector.get("pathogen_group", "")),
            str(selector.get("host_group", "")), str(selector.get("host_taxon_id", "")).removeprefix("ncbi:"))] = profile
    feature_ledger_path = run / "feature_screen_ledger.parquet"
    feature_ledger = pd.read_parquet(feature_ledger_path) if feature_ledger_path.exists() else pd.DataFrame()
    if not feature_ledger.empty:
        feature_ledger["pathogen_id"] = feature_ledger.pathogen_id.astype(str)
        feature_ledger["pmid"] = feature_ledger.pmid.astype(str)
        if "host_id" in feature_ledger:
            feature_ledger["host_id"] = feature_ledger.host_id.fillna("").astype(str)
        feature_ledger = feature_ledger.loc[feature_ledger.status.astype(str).eq("completed")].copy()
    feature_decisions = feature_ledger.to_dict("records") if not feature_ledger.empty else []

    aliases_by_id: dict[str, list[str]] = {}
    alias_path = staging / "taxonomy/pathogen_aliases.parquet"
    if alias_path.exists():
        aliases = pd.read_parquet(alias_path, columns=["pathogen_id", "name", "name_type"])
        aliases = aliases.loc[aliases.name_type.astype(str).str.casefold().eq("synonym")]
        aliases_by_id = {str(key): sorted(set(group.name.dropna().astype(str)))
                         for key, group in aliases.groupby("pathogen_id", sort=False)}
    added_aliases = added_aliases.loc[added_aliases.name_type.astype(str).str.casefold().eq("synonym")]
    for key, group in added_aliases.groupby("pathogen_id", sort=False):
        aliases_by_id[str(key)] = sorted(set(group.name.dropna().astype(str)))

    reference_dir = run.parents[1] / "references/release"
    wikipedia = _read_csv_index(reference_dir / "11_resolve_wikipedia_references.csv")
    gbif = _read_csv_index(reference_dir / "12_resolve_gbif_references.csv")
    taxonomy_by_id = pathogens.set_index("pathogen_id").to_dict("index")
    pair_by_key = {(str(row.pathogen_id), str(row.host_id)): row._asdict()
                   for row in pairs.itertuples(index=False)}
    pii_by_key = {(str(row.pathogen_id), str(row.host_id)): row._asdict()
                  for row in pii.itertuples(index=False)} if not pii.empty else {}
    host_table = pd.read_parquet(staging / "taxonomy/hosts.parquet")
    hosts = []
    for row in host_table.to_dict("records"):
        host_id = str(row["host_id"])
        hosts.append({**{key: _json_safe(value) for key, value in row.items()},
                      "label": HOST_LABELS.get(host_id, str(row["accepted_name"])),
                      "host_group": HOST_GROUPS.get(host_id, str(row.get("host_group") or "Animals"))})
    # Keep the host picker useful across releases, including hosts whose pair
    # profiles have not been acquired yet.
    from pathogen_atlas_v4.host_taxonomy.ncbi import load_ncbi_hosts
    fish_hosts = load_ncbi_hosts(
        project_root / "data/source_taxonomy/ncbi", ["8049", "8030"], release=release,
    )
    known_host_ids = {host["host_id"] for host in hosts}
    for row in fish_hosts.to_dict("records"):
        if row["host_id"] not in known_host_ids:
            host_id = str(row["host_id"])
            hosts.append({**{key: _json_safe(value) for key, value in row.items()},
                          "label": HOST_LABELS[host_id], "host_group": HOST_GROUPS[host_id]})

    human_source_path = site / "source/human_pair_data.json.gz"
    with gzip.open(human_source_path, "rt", encoding="utf-8") as source_file:
        human_source = json.load(source_file)
    feature_catalogue = yaml.safe_load((project_root / "configs/features.yaml").read_text(encoding="utf-8"))
    current_question_versions = {str(feature_id): str(question.get("question_version", ""))
                                for feature_id, question in feature_catalogue.get("questions", {}).items()}
    feature_decisions = [decision for decision in feature_decisions
        if str(decision.get("question_version", "")) == current_question_versions.get(str(decision.get("feature_id")), "")]
    feature_decisions_by_pmid: dict[str, list[dict[str, Any]]] = {}
    feature_decisions_by_abstract: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for decision in feature_decisions:
        feature_decisions_by_pmid.setdefault(str(decision.get("pmid", "")), []).append(decision)
        abstract_key = (str(decision.get("title_hash") or ""), str(decision.get("abstract_hash") or ""))
        feature_decisions_by_abstract.setdefault(abstract_key, []).append(decision)

    def decisions_for_evidence(record: dict[str, Any], pathogen_id: str, host_id: str) -> list[dict[str, Any]]:
        """Resolve decisions by PMID, falling back to exact abstract identity.

        The decision cache intentionally omits PMID from its semantic key. If two
        PMIDs have byte-identical title/abstract and the same scoped inputs, the
        second task reuses the first result and is not written as a duplicate row.
        """
        def applies(decision: dict[str, Any]) -> bool:
            scope = decision.get("scope")
            return ((scope == "pathogen" and decision.get("pathogen_id") == pathogen_id)
                or (scope == "host" and decision.get("host_id") == host_id)
                or (scope == "pathogen_host_pair" and decision.get("pathogen_id") == pathogen_id
                    and decision.get("host_id") == host_id))

        direct = [decision for decision in feature_decisions_by_pmid.get(str(record.get("pmid", "")), [])
                  if decision.get("decision_kind") == "feature_screen" and applies(decision)]
        if direct:
            return direct
        key = (str(record.get("title_hash") or ""), str(record.get("abstract_hash") or ""))
        if not all(key):
            return []
        return [decision for decision in feature_decisions_by_abstract.get(key, [])
                if decision.get("decision_kind") == "feature_screen" and applies(decision)]
    human_host_id = str(human_source["host_id"])
    human_profiles = {str(row["pathogen_id"]): row for row in human_source["profiles"]}
    human_profile_ids = set(human_profiles)
    human_evidence = human_source["evidence"]
    if not human_profile_ids.issubset(selected_pathogen_ids):
        missing = sorted(human_profile_ids - selected_pathogen_ids)[:5]
        raise ValueError(f"Human evidence profiles are absent from atlas taxonomy: {missing}")

    species_rows: list[dict[str, Any]] = []
    evidence_rows: list[dict[str, Any]] = []
    pair_profiles_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for pathogen_id, profile in human_profiles.items():
        pair_profiles_by_key[(pathogen_id, human_host_id)] = {
            "pair_selection_status": "supported",
            "relationship_context": "disease",
            "n_pathogen_pmids": int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0),
            "n_jev_accepted_pmids": int(profile.get("n_jev_accepted_pmids") or profile.get("n_pathogen_pmids") or 0),
            "n_jev_assessments": int(profile.get("n_jev_assessments") or profile.get("n_abstracts_fetched") or 0),
            "n_positive_records": int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0),
            "n_sampled_records": int(profile.get("n_abstracts_fetched") or 0),
            "n_abstracts_fetched": int(profile.get("n_abstracts_fetched") or 0),
            "n_pair_search_records": 0,
            "n_pubmed_search_records": _json_safe(profile.get("n_pubmed_search_records")),
            "estimated_pathogen_publications": _json_safe(profile.get("estimated_pathogen_publications")),
            "estimated_pathogen_publications_index": _json_safe(profile.get("estimated_pathogen_publications_index")),
            "pathogenicity_fraction": _json_safe(profile.get("pathogenicity_fraction")),
            "pathogenicity_fraction_index": _json_safe(profile.get("pathogenicity_fraction_index")),
            "pathogenic_interest_index": _json_safe(profile.get("pathogenic_interest_index")),
            "index_version": _json_safe(profile.get("index_version")),
            "feature_values": profile.get("feature_values", {}),
            "feature_counts": profile.get("feature_counts", {}),
            "feature_profile_version": "human-evidence-summary",
            "evidence_file": "evidence.json.gz",
            "wikipedia_url": profile.get("wikipedia_url"),
            "gbif_url": profile.get("gbif_url"),
            "x": _json_safe(profile.get("x")),
            "y": _json_safe(profile.get("y")),
        }
    for record in human_evidence:
        record = {**record, "host_id": human_host_id, "pathogen_id": str(record.get("pathogen_id") or record["organism_id"])}
        evidence_rows.append(record)

    taxonomy_by_id = pathogens.set_index("pathogen_id").to_dict("index")
    for key, pair in pairs.set_index(["pathogen_id", "host_id"]).to_dict("index").items():
        pair_profiles_by_key[key] = {
            "pair_selection_status": str(pair.get("selection_status") or "incomplete"),
            "relationship_context": str(pair.get("relationship_context") or "disease"),
            "n_pathogen_pmids": int(pair.get("n_positive_records") or 0),
            "n_jev_accepted_pmids": int(pair.get("n_positive_records") or 0),
            "n_positive_records": int(pair.get("n_positive_records") or 0),
            "n_sampled_records": int(pair.get("n_sampled_records") or 0),
            "n_abstracts_fetched": int(pair.get("n_sampled_records_with_abstract") or 0),
            "n_abstracts_screened": int(pair.get("n_completed_screenings") or 0),
            "n_pair_search_records": int(pair.get("n_pair_search_records") or 0),
            "n_pubmed_search_records": None,
            "estimated_pathogen_publications": None,
            "estimated_pathogen_publications_index": None,
            "pathogenicity_fraction": None,
            "pathogenicity_fraction_index": None,
            "pathogenic_interest_index": None,
            "feature_values": {},
            "feature_counts": {},
            "evidence_file": "evidence.json.gz",
        }
    for key, pair in pair_by_key.items():
        pathogen_id, host_id = key
        taxonomy = taxonomy_by_id.get(pathogen_id)
        if taxonomy is None:
            raise ValueError(f"Pair selection references pathogen absent from atlas taxonomy: {pathogen_id}")
        host = next((item for item in hosts if str(item["host_id"]) == host_id), None)
        if host is None:
            raise ValueError(f"Pair selection references host absent from staging taxonomy: {host_id}")
        metric = pii_by_key.get(key, {})
        matched_evidence = evidence_by_pair.get(key, pd.DataFrame())
        matched_evidence_by_pmid = {str(row["pmid"]): row for row in matched_evidence.to_dict("records")}
        positive_count = int(pair.get("n_positive_records") or 0)
        if len(matched_evidence) != positive_count:
            raise ValueError(f"Pair evidence count mismatch for {pathogen_id}, {host_id}")
        if positive_count:
            host_taxon = host_taxonomy_by_id.loc[host_id]
            selector = (str(taxonomy_by_id[pathogen_id]["pathogen_group"]), str(host_taxon["host_group"]),
                        str(host_taxon["taxon_id"]).removeprefix("ncbi:"))
            resolved_profile = resolved_feature_profiles.get(selector)
            if resolved_profile is None:
                raise ValueError(f"No resolved feature profile for accepted pair {pathogen_id}, {host_id}")
            expected_feature_ids = {str(feature["feature_id"]) for feature in resolved_profile.get("ordered_features", [])}
            if expected_feature_ids and feature_ledger.empty:
                raise ValueError(f"Accepted pair {pathogen_id}, {host_id} has no feature-screen ledger")
            for accepted_pmid in matched_evidence.pmid.astype(str):
                scored_ids = set()
                evidence_record = matched_evidence_by_pmid[str(accepted_pmid)]
                for decision in decisions_for_evidence(evidence_record, pathogen_id, host_id):
                    scored_ids.add(str(decision.get("feature_id")))
                if expected_feature_ids - scored_ids:
                    missing_features = sorted(expected_feature_ids - scored_ids)
                    raise ValueError(f"Feature screening is incomplete for {pathogen_id}, {host_id}, PMID {accepted_pmid}: "
                                     f"{len(missing_features)} applicable questions lack completed Jev decisions")
        pair_profiles_by_key[key].update({
            "n_pubmed_search_records": _json_safe(metric.get("n_pubmed_search_records")),
            "estimated_pathogen_publications": _json_safe(metric.get("estimated_pathogen_publications")),
            "estimated_pathogen_publications_index": _json_safe(metric.get("estimated_pathogen_publications_index")),
            "pathogenicity_fraction": _json_safe(metric.get("pathogenicity_fraction")),
            "pathogenicity_fraction_index": _json_safe(metric.get("pathogenicity_fraction_index")),
            "pathogenic_interest_index": _json_safe(metric.get("pathogenic_interest_index")),
            "index_version": metric.get("index_version"),
        })
        for record in matched_evidence.to_dict("records"):
            article = article_by_pmid.loc[str(record["pmid"])] if str(record["pmid"]) in article_by_pmid.index else {}
            if isinstance(article, pd.DataFrame):
                article = article.iloc[-1]
            accepted_record = {
                "organism_id": pathogen_id,
                "organism_group": str(taxonomy["pathogen_group"]),
                "accepted_name": str(taxonomy["accepted_name"]),
                "pathogen_id": pathogen_id,
                "host_id": host_id,
                "relationship_context": str(pair.get("relationship_context") or "disease"),
                "pmid": str(record["pmid"]),
                "title": str(article.get("title") or ""),
                "publication_date": _json_safe(article.get("publication_date")),
                "publication_year": _json_safe(article.get("publication_year")),
                "doi": str(article.get("doi") or ""),
                "p_yes": _json_safe(record.get("p_yes")),
                "decision": str(record["decision"]),
                "accepted_pathogen_evidence": True,
                "abstract_hash": str(record.get("abstract_hash") or ""),
                "feature_ids_yes": [],
                "feature_confidences": {},
            }
            matching_features = decisions_for_evidence(record, pathogen_id, host_id)
            for decision in matching_features:
                feature_id = str(decision["feature_id"])
                accepted_record["feature_confidences"][feature_id] = _json_safe(decision.get("p_yes"))
                if str(decision.get("answer", "")).casefold() == "yes":
                    accepted_record["feature_ids_yes"].append(feature_id)
            evidence_rows.append(accepted_record)

        if feature_decisions:
            scored_pmids_by_feature: dict[str, set[str]] = {}
            yes_pmids_by_feature: dict[str, set[str]] = {}
            accepted_pmids = set(matched_evidence.pmid.astype(str)) if not matched_evidence.empty else set()
            for evidence_record in matched_evidence.to_dict("records"):
                pmid = str(evidence_record["pmid"])
                for decision in decisions_for_evidence(evidence_record, pathogen_id, host_id):
                    feature_id = str(decision["feature_id"])
                    scored_pmids_by_feature.setdefault(feature_id, set()).add(pmid)
                    if str(decision.get("answer", "")).casefold() == "yes":
                        yes_pmids_by_feature.setdefault(feature_id, set()).add(pmid)
            feature_ids = set(scored_pmids_by_feature)
            feature_values, feature_counts = {}, {}
            for feature_id in feature_ids:
                scored_pmids = scored_pmids_by_feature[feature_id]
                yes_pmids = yes_pmids_by_feature.get(feature_id, set())
                n_scored, n_yes = len(scored_pmids), len(yes_pmids)
                feature_values[feature_id] = n_yes / n_scored if n_scored else None
                feature_counts[feature_id] = {"n_yes": n_yes, "n_scored": n_scored,
                    "n_expected": len(accepted_pmids), "n_failed": 0, "n_pending": max(0, len(accepted_pmids) - n_scored)}
            pair_profiles_by_key[key]["feature_values"].update(feature_values)
            pair_profiles_by_key[key]["feature_counts"].update(feature_counts)
            sources = pair_profiles_by_key[key].setdefault("feature_sources", {})
            sources.update({feature_id: "jev" for feature_id in feature_ids})
            pair_profiles_by_key[key]["feature_profile_version"] = "v4-trait-resolved-jev-v1"

    # Preserve pair profiles from the previously built Explorer. A new run
    # usually covers one host and pathogen cohort, so absent pairs are still
    # useful historical database records (especially the full Rice screen).
    refreshed_pair_keys = set(pair_profiles_by_key)
    prior_taxa_root = (Path(dist_root) / "taxa") if dist_root is not None else (site / "taxa")
    valid_host_ids = {str(host["host_id"]) for host in hosts}
    for shard_path in sorted(prior_taxa_root.glob("*.json.gz")):
        try:
            with gzip.open(shard_path, "rt", encoding="utf-8") as shard_file:
                prior_rows = json.load(shard_file)
        except (OSError, json.JSONDecodeError):
            continue
        for prior_row in prior_rows:
            if not isinstance(prior_row, list) or len(prior_row) < 11:
                continue
            pathogen_id = str(prior_row[0])
            if pathogen_id not in taxonomy_by_id:
                continue
            for host_id, profile in (prior_row[-1] or {}).items():
                key = (pathogen_id, str(host_id))
                if str(host_id) in valid_host_ids and key not in pair_profiles_by_key:
                    profile.setdefault("n_jev_assessments", int(
                        profile.get("n_abstracts_screened")
                        or profile.get("n_abstracts_fetched")
                        or profile.get("n_sampled_records")
                        or 0
                    ))
                    pair_profiles_by_key[key] = profile

    # Historical pair profiles and their scores are carried forward from the
    # previous taxonomy shards. Keep their detailed evidence in sync as well.
    # Older Explorer builds could preserve a profile while dropping its rows
    # from evidence.json.gz; recover those rows from the previous payload first
    # and, when necessary, from the immutable run ledgers on disk.
    recovered_historical_pair_keys: set[tuple[str, str]] = set()
    present_evidence_keys = {
        (str(row.get("organism_id") or row.get("pathogen_id") or ""),
         str(row.get("host_id") or ""), str(row.get("pmid") or ""))
        for row in evidence_rows
    }
    retained_pair_keys = set(pair_profiles_by_key)
    for record in previous_evidence:
        pathogen_id = str(record.get("organism_id") or record.get("pathogen_id") or "")
        host_id = str(record.get("host_id") or "")
        pmid = str(record.get("pmid") or "")
        if (pathogen_id, host_id) not in retained_pair_keys or (pathogen_id, host_id) in refreshed_pair_keys or not pmid:
            continue
        key = (pathogen_id, host_id, pmid)
        if key not in present_evidence_keys:
            evidence_rows.append(record)
            present_evidence_keys.add(key)

    # Backfill older builds whose evidence payload was already incomplete.
    # Run ledgers are the source of truth for supported Jev assessments and
    # completed feature decisions. Restrict recovery to retained pairs and
    # records missing from the current/previous payload.
    runs_root = project_root / "data/v4/runs"
    for historical_run in sorted(runs_root.iterdir()) if runs_root.exists() else []:
        evidence_path = historical_run / "evidence/pathogen_host_evidence.parquet"
        if not evidence_path.exists():
            continue
        try:
            historical = pd.read_parquet(evidence_path)
            historical["pathogen_id"] = historical.pathogen_id.astype(str)
            historical["host_id"] = historical.host_id.astype(str)
            historical["pmid"] = historical.pmid.astype(str)
            historical = historical.loc[historical.decision.astype(str).eq("supported")]
            historical = historical.loc[
                historical.apply(lambda row: (row.pathogen_id, row.host_id) in retained_pair_keys
                                 and (row.pathogen_id, row.host_id) not in refreshed_pair_keys, axis=1)
            ]
            if historical.empty:
                continue
            ledger_path = historical_run / "feature_screen_ledger.parquet"
            historical_ledger = pd.read_parquet(ledger_path) if ledger_path.exists() else pd.DataFrame()
            if not historical_ledger.empty:
                historical_ledger["pathogen_id"] = historical_ledger.pathogen_id.astype(str)
                historical_ledger["host_id"] = historical_ledger.host_id.fillna("").astype(str)
                historical_ledger["pmid"] = historical_ledger.pmid.astype(str)
                historical_ledger = historical_ledger.loc[
                    historical_ledger.status.astype(str).eq("completed")
                    & historical_ledger.decision_kind.astype(str).eq("feature_screen")
                ]
            feature_index = {}
            if not historical_ledger.empty:
                for (pathogen_id, host_id, pmid), decisions in historical_ledger.groupby(
                    ["pathogen_id", "host_id", "pmid"], sort=False
                ):
                    feature_index[(str(pathogen_id), str(host_id), str(pmid))] = decisions.to_dict("records")
            for record in historical.to_dict("records"):
                key = (str(record["pathogen_id"]), str(record["host_id"]), str(record["pmid"]))
                if not key[2] or key in present_evidence_keys:
                    continue
                taxonomy = taxonomy_by_id.get(key[0], {})
                article = article_by_pmid.loc[key[2]] if key[2] in article_by_pmid.index else {}
                if isinstance(article, pd.DataFrame):
                    article = article.iloc[-1]
                feature_decisions = feature_index.get(key, [])
                title = str(article.get("title") or "")
                if not title and feature_decisions:
                    title = str(feature_decisions[0].get("title_text") or "")
                recovered = {
                    "organism_id": key[0],
                    "organism_group": str(taxonomy.get("pathogen_group") or ""),
                    "accepted_name": str(taxonomy.get("accepted_name") or ""),
                    "pathogen_id": key[0], "host_id": key[1],
                    "relationship_context": str(record.get("relationship_context") or "disease"),
                    "pmid": key[2], "title": title,
                    "publication_date": _json_safe(article.get("publication_date")),
                    "publication_year": _json_safe(article.get("publication_year")),
                    "doi": str(article.get("doi") or ""),
                    "p_yes": _json_safe(record.get("p_yes")),
                    "decision": str(record.get("decision") or "supported"),
                    "accepted_pathogen_evidence": True,
                    "abstract_hash": str(record.get("abstract_hash") or ""),
                    "feature_ids_yes": [], "feature_confidences": {},
                    "run_id": str(record.get("run_id") or historical_run.name),
                }
                for decision in feature_decisions:
                    feature_id = str(decision.get("feature_id") or "")
                    if not feature_id:
                        continue
                    recovered["feature_confidences"][feature_id] = _json_safe(decision.get("p_yes"))
                    if str(decision.get("answer") or "").casefold() == "yes":
                        recovered["feature_ids_yes"].append(feature_id)
                evidence_rows.append(recovered)
                present_evidence_keys.add(key)
                recovered_historical_pair_keys.add((key[0], key[1]))
        except (OSError, ValueError, KeyError) as exc:
            print(f"Skipping historical evidence recovery from {historical_run.name}: {exc}")

    shared_virus_coordinates = 0
    rice_host_id = "ncbi:4530"
    for pathogen_id, profile in human_profiles.items():
        taxonomy = taxonomy_by_id.get(pathogen_id, {})
        rice_profile = pair_profiles_by_key.get((pathogen_id, rice_host_id))
        if (taxonomy.get("pathogen_group") == "viruses" and rice_profile is not None
                and profile.get("x") is not None and profile.get("y") is not None):
            rice_profile["x"] = _json_safe(profile["x"])
            rice_profile["y"] = _json_safe(profile["y"])
            rice_profile["coordinate_source_host_id"] = human_host_id
            shared_virus_coordinates += 1

    # Refit both maps for every host × pathogen-group cohort with complete,
    # variable feature data. Incomplete profiles keep their evidence and catalog
    # data; eligible profiles receive fresh, method-specific coordinates.
    from pathogen_atlas_v4.feature_embeddings import fit_feature_embeddings
    embedding_runs: dict[str, Any] = dict(prior_embedding_runs)
    embedding_host_ids: set[str] = set()
    grouped_profiles: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for (pathogen_id, host_id), profile in pair_profiles_by_key.items():
        group = str(taxonomy_by_id.get(pathogen_id, {}).get("pathogen_group") or "")
        if group and any(value is not None for value in (profile.get("feature_values") or {}).values()):
            grouped_profiles.setdefault((host_id, group), []).append(
                {**profile, "pathogen_id": pathogen_id}
            )
    for (host_id, group), profiles in sorted(grouped_profiles.items()):
        embedding_key = f"{host_id}::{group}"
        try:
            fitted_embeddings = fit_feature_embeddings(
                profiles, host_id=host_id, pathogen_group=group,
            )
        except ValueError as exc:
            # A small or incomplete cohort stays available in the catalog.
            print(f"Skipping 2D maps for {host_id} × {group}: {exc}")
            continue
        for method, (coordinates, metadata) in fitted_embeddings.items():
            for pathogen_id, (x, y) in coordinates.items():
                profile = pair_profiles_by_key[(pathogen_id, host_id)]
                profile.setdefault("embedding_coordinates", {})[method] = {"x": x, "y": y}
                if method == "IMSTE":
                    # Retain the original top-level fields for older clients and exports.
                    profile["x"], profile["y"] = x, y
                    profile["coordinate_source_host_id"] = host_id
                    profile["coordinate_embedding_key"] = embedding_key
            embedding_runs[f"{embedding_key}::{method}"] = metadata
            if method == "IMSTE":
                # Preserve the existing metadata key while exposing method-specific keys.
                embedding_runs[embedding_key] = metadata
        embedding_host_ids.add(host_id)

    pair_profiles_by_taxon: dict[str, dict[str, dict[str, Any]]] = {}
    for (pathogen_id, host_id), profile in pair_profiles_by_key.items():
        pair_profiles_by_taxon.setdefault(pathogen_id, {})[host_id] = profile

    for pathogen_id, taxonomy in taxonomy_by_id.items():
        group = str(taxonomy["pathogen_group"])
        item = {
            "organism_id": pathogen_id,
            "pathogen_id": pathogen_id,
            "organism_group": group,
            "accepted_name": str(taxonomy["accepted_name"]),
            "source_taxon_id": str(taxonomy.get("taxon_id") or ""),
            "taxonomy_source": str(taxonomy.get("taxonomy_source") or ""),
            "taxonomy_release": str(taxonomy.get("taxonomy_release") or ""),
            "rank": str(taxonomy.get("rank") or ""),
            "parent_taxon_id": _json_safe(taxonomy.get("parent_taxon_id")),
            "taxonomy_path": str(taxonomy.get("taxonomy_path") or ""),
            "mapping_status": str(taxonomy.get("mapping_status") or ""),
            "synonyms": aliases_by_id.get(pathogen_id, []),
            "pair_profiles": pair_profiles_by_taxon.get(pathogen_id),
        }
        species_rows.append(item)

    species_rows.sort(key=lambda row: (row["organism_group"], row["accepted_name"].casefold(), row["organism_id"]))
    evidence_rows.sort(key=lambda row: (row["pathogen_id"], row["host_id"], str(row["pmid"])))
    groups = sorted({row["organism_group"] for row in species_rows})
    feature_profile_counts = {
        str(host["host_id"]): {group: 0 for group in groups}
        for host in hosts
    }
    profile_availability_counts = {
        str(host["host_id"]): {group: 0 for group in groups}
        for host in hosts
    }
    for (pathogen_id, host_id), profile in pair_profiles_by_key.items():
        values = profile.get("feature_values") or {}
        taxonomy = taxonomy_by_id.get(pathogen_id, {})
        group = str(taxonomy.get("pathogen_group") or "")
        has_features = any(value is not None for value in values.values())
        has_coordinates = all(
            value is not None and math.isfinite(float(value))
            for value in (profile.get("x"), profile.get("y"))
        )
        if group in groups and has_features:
            feature_profile_counts.setdefault(host_id, {}).setdefault(group, 0)
            feature_profile_counts[host_id][group] += 1
        if group in groups and (has_features or has_coordinates):
            profile_availability_counts.setdefault(host_id, {}).setdefault(group, 0)
            profile_availability_counts[host_id][group] += 1
    decisions = yaml.safe_load((Path(__file__).parents[3] / "configs/prompts.yaml").read_text(encoding="utf-8"))
    screen_version = str(decisions["pathogen_screen"]["question_version"])
    supported_pairs = sum(profile.get("pair_selection_status") in {"accepted", "supported"}
                          for profile in pair_profiles_by_key.values())
    feature_profiles = {group: list(human_source.get("feature_profiles", {}).get(group, [])) for group in groups}
    feature_data = {
        "catalog_version": feature_catalogue.get("catalogue_version", "atlas-feature-catalog-v1"),
        "profile_version": feature_catalogue.get("profile_version", "host-aware-profile-v1"),
        "features": human_source.get("features", {}),
        "feature_hierarchy": human_source.get("feature_hierarchy", []),
        "feature_profiles": feature_profiles,
    }
    # Add the features assessed for this run to the Explorer catalogue while
    # retaining inherited profiles and Rice literature features.
    for profile in resolved_feature_profiles.values():
        group = str(profile.get("selectors", {}).get("pathogen_group", ""))
        current = feature_data["feature_profiles"].setdefault(group, [])
        current.extend(str(item["feature_id"]) for item in profile.get("ordered_features", [])
                       if item.get("feature_id"))
        feature_data["feature_profiles"][group] = list(dict.fromkeys(current))
    active_feature_ids = {str(feature_id) for ids in feature_data["feature_profiles"].values() for feature_id in ids}
    for feature_id, question in feature_catalogue.get("questions", {}).items():
        if feature_id not in active_feature_ids:
            continue
        existing = feature_data["features"].get(feature_id, {})
        fields = {key: question[key] for key in (
            "label", "definition", "question", "question_claim", "question_version", "scope",
            "applicability", "applies_to_label", "dimension_id", "dimension", "category_id", "category",
        ) if key in question}
        feature_data["features"][feature_id] = {**existing, **fields}
        dimension_id, category_id = question.get("dimension_id"), question.get("category_id")
        if not dimension_id or not category_id:
            continue
        dimension = next((item for item in feature_data["feature_hierarchy"]
                          if item.get("id") == dimension_id), None)
        if dimension is None:
            dimension = {"id": dimension_id, "label": question.get("dimension", dimension_id), "categories": []}
            feature_data["feature_hierarchy"].append(dimension)
        category = next((item for item in dimension.setdefault("categories", [])
                         if item.get("id") == category_id), None)
        if category is None:
            category = {"id": category_id, "label": question.get("category", category_id), "feature_ids": []}
            dimension["categories"].append(category)
        category["feature_ids"] = list(dict.fromkeys([*category.get("feature_ids", []), str(feature_id)]))
    group_rows = {
        group: [[row["organism_id"], row["accepted_name"], row["source_taxon_id"],
                 row["taxonomy_source"], row["taxonomy_release"], row["rank"],
                 row["parent_taxon_id"], row["taxonomy_path"], row["mapping_status"],
                 row["synonyms"], row.get("pair_profiles")]
                for row in species_rows if row["organism_group"] == group]
        for group in groups
    }
    taxa_shard_versions = {
        group: hashlib.sha256(json.dumps(rows, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()).hexdigest()[:16]
        for group, rows in group_rows.items()
    }
    stats = {
        "run_id": str(resolved.get("run_id", run.name)),
        "screening_question_version": screen_version,
        "screening_yes_probability_threshold": 0.7,
        "pathogen_species": len(species_rows),
        "pathogen_species_by_group": {
            group: sum(row["organism_group"] == group for row in species_rows)
            for group in groups
        },
        "hosts": len(hosts),
        "possible_host_pathogen_pairs": len(species_rows) * len(hosts),
        "pair_profiles_available": len(pair_profiles_by_key),
        "supported_pairs": supported_pairs,
        "current_run_candidate_pairs": len(pairs),
        "current_run_pubmed_hit_pairs": int(pd.to_numeric(pairs.n_pair_search_records, errors="coerce").fillna(0).gt(0).sum()),
        "current_run_sampled_records": int(pairs.n_sampled_records.sum()),
        "human_supported_profiles": len(human_profiles),
        "human_supporting_pmids": len(human_evidence),
        "shared_virus_coordinates_in_rice_view": shared_virus_coordinates,
        "supported_pmids": len(evidence_rows),
        "feature_screen_decisions": len(feature_decisions),
    }
    species_columns = [
        "organism_id", "accepted_name", "source_taxon_id", "taxonomy_source",
        "taxonomy_release", "rank", "parent_taxon_id", "taxonomy_path",
        "mapping_status", "synonyms", "pair_profiles",
    ]
    embedding_methods_by_host: dict[str, set[str]] = {}
    for metadata in embedding_runs.values():
        if not isinstance(metadata, dict):
            continue
        host_id = str(metadata.get("host_id") or "")
        method = str(metadata.get("embedding") or "")
        if host_id and method in {"IMSTE", "UMAP"}:
            embedding_methods_by_host.setdefault(host_id, set()).add(method)

    def embedding_context(host_id: str) -> str:
        methods = [method for method in ("IMSTE", "UMAP")
                   if method in embedding_methods_by_host.get(host_id, set())]
        if not methods:
            return "No stable 2D embedding is available for the selected host and pathogen group."
        label = " and ".join(methods)
        noun = "map" if len(methods) == 1 else "maps"
        verb = "uses" if len(methods) == 1 else "use"
        return f"Host- and pathogen-group-specific {label} {noun} {verb} complete scored feature profiles."

    data = {
        "schema_version": "v4-host-pathogen-atlas-v1",
        "groups": groups,
        "taxonomy_provenance": _taxonomy_provenance(groups, species_rows, ncbi_retrieved_date),
        "hosts": hosts,
        "feature_profile_counts": feature_profile_counts,
        "profile_availability_counts": profile_availability_counts,
        "taxa_shard_versions": taxa_shard_versions,
        "species_columns": species_columns,
        "taxonomy_ranks": sorted({part.split("=", 1)[0] for row in species_rows
                                   for part in row["taxonomy_path"].split(";") if "=" in part}),
        "coordinates_ready": any(profile.get("x") is not None and profile.get("y") is not None
                                  for values in pair_profiles_by_taxon.values() for profile in values.values()),
        "coordinate_host_ids": sorted({
            host_id for (_, host_id), profile in pair_profiles_by_key.items()
            if all(value is not None and math.isfinite(float(value))
                   for value in (profile.get("x"), profile.get("y")))
        }),
        "coordinate_contexts": {
            **{host_id: embedding_context(host_id)
               for host_id in embedding_host_ids},
            human_host_id: embedding_context(human_host_id),
            rice_host_id: "Shared saved virus positions from the human-disease embedding; pair metrics and PMID evidence remain rice-specific.",
        },
        "stats": stats,
        "embedding_runs": embedding_runs,
    }
    from pathogen_atlas_v4.rice_virus_profiles import add_rice_virus_profiles, merge_current_rice_pair_status
    data, species_rows, rice_features, rice_evidence = add_rice_virus_profiles(
        data, species_rows, pathogens,
        site / "source/rice_virus_screen.json",
    )
    species_rows_by_id = {str(row.get("organism_id")): row for row in species_rows}
    # A completed current-run decision supersedes the older Rice screen for
    # that pair. Pairs without a usable current sample keep their historical
    # Rice status and evidence, so a narrow pilot cannot erase atlas coverage.
    current_rice_decisions = {}
    for pair in pairs.loc[pairs.host_id.astype(str).eq(rice_host_id)
                          & pd.to_numeric(pairs.n_sampled_records, errors="coerce").fillna(0).gt(0)].to_dict("records"):
        status = str(pair.get("selection_status") or "incomplete")
        if status not in {"accepted", "rejected"}:
            continue
        key = (str(pair["pathogen_id"]), rice_host_id)
        current_rice_decisions[key] = pair
        taxon = species_rows_by_id.get(key[0])
        if taxon is None:
            continue
        profile = (taxon.get("pair_profiles") or {}).setdefault(rice_host_id, {})
        current_profile = pair_profiles_by_key.get(key, {})
        profile.update(merge_current_rice_pair_status(profile, pair))
        for field in ("n_pubmed_search_records", "estimated_pathogen_publications",
                      "estimated_pathogen_publications_index", "pathogenicity_fraction",
                      "pathogenicity_fraction_index", "pathogenic_interest_index", "index_version"):
            if field in current_profile:
                profile[field] = current_profile[field]
    rejected_current_rice = {key for key, pair in current_rice_decisions.items()
                             if str(pair.get("selection_status")) == "rejected"}
    rice_evidence = [record for record in rice_evidence
        if (str(record.get("pathogen_id")), str(record.get("host_id"))) not in rejected_current_rice]
    # Merge the rice-specific traits into the virus feature catalogue. Existing
    # human profiles and their coordinates/evidence are retained byte-for-byte
    # in their separate source shard.
    feature_data["features"] = {**feature_data.get("features", {}), **rice_features["features"]}
    feature_data["feature_profiles"]["viruses"] = list(dict.fromkeys(
        [*feature_data["feature_profiles"].get("viruses", []), *rice_features["feature_profiles"]["viruses"]]
    ))
    # Rice features need hierarchy entries as well as profile membership;
    # otherwise the Explorer counts them but has no category in which to render them.
    dimensions = feature_data.setdefault("feature_hierarchy", [])
    for rice_dimension in rice_features.get("feature_hierarchy", []):
        dimension = next((item for item in dimensions if item.get("id") == rice_dimension.get("id")), None)
        if dimension is None:
            dimension = {"id": rice_dimension["id"], "label": rice_dimension["label"], "categories": []}
            dimensions.append(dimension)
        for rice_category in rice_dimension.get("categories", []):
            category = next((item for item in dimension.setdefault("categories", [])
                             if item.get("id") == rice_category.get("id")), None)
            if category is None:
                category = {"id": rice_category["id"], "label": rice_category["label"], "feature_ids": []}
                dimension["categories"].append(category)
            category["feature_ids"] = list(dict.fromkeys([
                *category.get("feature_ids", []), *rice_category.get("feature_ids", [])
            ]))
    evidence_rows.extend(rice_evidence)
    # Preserve one record per pair and PMID, preferring the current run's
    # complete assessment when both source ledgers contain the same record.
    evidence_rows.sort(key=lambda row: (
        str(row.get("organism_id", "")), str(row.get("host_id", "")), str(row.get("pmid", "")),
        0 if row.get("run_id") == str(resolved.get("run_id", run.name)) else 1))
    evidence_rows = list({(str(row.get("organism_id", "")), str(row.get("host_id", "")), str(row.get("pmid", ""))): row
                          for row in reversed(evidence_rows)}.values())
    evidence_counts: dict[tuple[str, str], int] = {}
    for row in evidence_rows:
        if row.get("accepted_pathogen_evidence"):
            key = (str(row.get("organism_id") or row.get("pathogen_id") or ""), str(row.get("host_id") or ""))
            evidence_counts[key] = evidence_counts.get(key, 0) + 1
    for key in refreshed_pair_keys | recovered_historical_pair_keys:
        profile = pair_profiles_by_key[key]
        expected = int(profile.get("n_jev_accepted_pmids") or profile.get("n_positive_records") or 0)
        if expected and evidence_counts.get(key, 0) != expected:
            raise ValueError(
                f"Explorer evidence mismatch for {key[0]}, {key[1]}: "
                f"profile reports {expected} accepted Jev PMIDs but evidence payload contains "
                f"{evidence_counts.get(key, 0)}"
            )
    groups = sorted({row["organism_group"] for row in species_rows})
    data["groups"] = groups
    data["stats"]["pathogen_species"] = len(species_rows)
    data["stats"]["rice_supported_virus_taxa"] = sum(
        str(((row.get("pair_profiles") or {}).get(rice_host_id) or {}).get("pair_selection_status"))
        in {"accepted", "supported"}
        for row in species_rows if row.get("organism_group") == "viruses")
    data["stats"]["pathogen_species_by_group"] = {
        group: sum(row["organism_group"] == group for row in species_rows) for group in groups
    }
    data["taxonomy_provenance"] = _taxonomy_provenance(groups, species_rows, ncbi_retrieved_date)
    group_rows = {
        group: [[row["organism_id"], row["accepted_name"], row["source_taxon_id"],
                 row["taxonomy_source"], row["taxonomy_release"], row["rank"],
                 row["parent_taxon_id"], row["taxonomy_path"], row["mapping_status"],
                 row["synonyms"], row.get("pair_profiles")]
                for row in species_rows if row["organism_group"] == group]
        for group in groups
    }
    taxa_shard_versions = {
        group: hashlib.sha256(json.dumps(rows, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()).hexdigest()[:16]
        for group, rows in group_rows.items()
    }
    data["taxa_shard_versions"] = taxa_shard_versions
    data["species_columns"] = species_columns
    data["taxonomy_ranks"] = sorted({part.split("=", 1)[0] for row in species_rows
                                      for part in row["taxonomy_path"].split(";") if "=" in part})
    evidence_rows.sort(key=lambda row: (row["pathogen_id"], row["host_id"], str(row["pmid"])))
    site.mkdir(parents=True, exist_ok=True)
    for name, value in (("explorer-data.json.gz", data), ("evidence.json.gz", evidence_rows)):
        with (site / name).open("wb") as output:
            with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as compressed:
                compressed.write(json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode())
    taxa_root = site / "taxa"
    taxa_root.mkdir(exist_ok=True)
    for group, rows in group_rows.items():
        payload = json.dumps(rows, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
        with (taxa_root / f"{group}.json.gz").open("wb") as output:
            with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as compressed:
                compressed.write(payload)
    (site / "feature-catalog.json").write_text(json.dumps(feature_data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    if dist_root is not None:
        _copy_static_site(site, Path(dist_root))
    return {**stats, "species_profiles": len(species_rows), "evidence_records": len(evidence_rows),
            "hosts": [host["host_id"] for host in hosts], "payload": str(site / "explorer-data.json.gz")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True, help="Completed V4 run directory containing host-pair screening")
    parser.add_argument("--staging-root", type=Path, required=True, help="Pinned atlas taxonomy/literature staging directory")
    parser.add_argument("--site-root", type=Path, default=Path("sites/pathogen-atlas-explorer"))
    parser.add_argument("--dist-root", type=Path, default=Path("sites/pathogen-atlas-explorer/dist"))
    args = parser.parse_args()
    print(json.dumps(build_site_data(args.run_dir, args.staging_root, args.site_root, dist_root=args.dist_root), indent=2))


if __name__ == "__main__":
    main()
