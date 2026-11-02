"""Build the small public Explorer payload from validated V4 staging tables."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

HOST_LABELS = {
    "ncbi:9606": "Human · Homo sapiens",
    "ncbi:9913": "Cattle · Bos taurus",
    "ncbi:9940": "Sheep · Ovis aries",
    "ncbi:4577": "Maize · Zea mays",
    "ncbi:4530": "Rice · Oryza sativa",
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


def build_site_data(staging: Path, site_root: Path) -> dict[str, Any]:
    pathogens = pd.read_parquet(staging / "taxonomy/pathogens.parquet")
    aliases = pd.read_parquet(staging / "taxonomy/pathogen_aliases.parquet")
    hosts = pd.read_parquet(staging / "taxonomy/hosts.parquet")
    evidence = pd.read_parquet(staging / "evidence/pathogen_host_evidence.parquet")
    if "decision" in evidence:
        evidence = evidence.loc[evidence["decision"].astype(str).eq("supported")].copy()
    evidence = evidence.loc[evidence["evidence_span"].fillna("").astype(str).str.strip().ne("")].copy()
    pathogen_ids = set(evidence["pathogen_id"].astype(str))
    pathogens = pathogens.loc[pathogens["pathogen_id"].astype(str).isin(pathogen_ids)].copy()
    pmids = evidence["pmid"].astype(str).unique().tolist()

    # Read only article metadata needed for the source list and lazy profile detail.
    articles_path = staging / "literature/articles.parquet"
    articles_schema = pq.ParquetFile(articles_path).schema.names
    article_columns = [name for name in ("pmid", "title", "publication_date", "publication_year") if name in articles_schema]
    articles = pd.read_parquet(articles_path, columns=article_columns, filters=[("pmid", "in", pmids)])
    articles["pmid"] = articles["pmid"].astype(str)
    article_by_pmid = {str(row["pmid"]): row for row in articles.to_dict("records")}

    alias_map: dict[str, list[str]] = {}
    accepted_names = {str(row.pathogen_id): str(row.accepted_name) for row in pathogens.itertuples(index=False)}
    for row in aliases.loc[aliases["pathogen_id"].astype(str).isin(pathogen_ids)].to_dict("records"):
        pathogen_id = str(row["pathogen_id"])
        name = row.get("name")
        if name and str(name).casefold() != accepted_names.get(pathogen_id, "").casefold():
            alias_map.setdefault(pathogen_id, []).append(str(name))

    match_path = staging / "literature/pathogen_article_matches.parquet"
    match_cols = pq.ParquetFile(match_path).schema.names
    if {"pathogen_id", "pmid", "selected_for_fetch"}.issubset(match_cols):
        matches = pd.read_parquet(match_path, columns=["pathogen_id", "pmid", "selected_for_fetch"])
        matches = matches.loc[matches["pathogen_id"].astype(str).isin(pathogen_ids)].copy()
        search_counts = matches.groupby("pathogen_id")["pmid"].nunique().to_dict()
        fetched = matches.loc[matches["selected_for_fetch"].fillna(False).astype(bool)]
        fetched_counts = fetched.groupby("pathogen_id")["pmid"].nunique().to_dict()
    else:
        search_counts, fetched_counts = {}, {}

    evidence_rows = []
    evidence_counts: dict[str, int] = {}
    for row in evidence.to_dict("records"):
        pathogen_id = str(row["pathogen_id"])
        pmid = str(row["pmid"])
        article = article_by_pmid.get(pmid, {})
        evidence_counts[pathogen_id] = evidence_counts.get(pathogen_id, 0) + 1
        evidence_rows.append({
            "organism_id": pathogen_id,
            "pathogen_id": pathogen_id,
            "host_id": str(row["host_id"]),
            "relationship_context": str(row.get("relationship_context", "disease")),
            "pmid": pmid,
            "title": _json_safe(article.get("title")),
            "publication_date": _json_safe(article.get("publication_date")),
            "p_yes": _json_safe(row.get("confidence")),
            "evidence_span": _json_safe(row.get("evidence_span")),
            "evidence_span_type": str(row.get("evidence_span_type", "full_abstract_legacy")),
            "provenance_status": str(row.get("provenance_status", "legacy")),
            "source_question_version": _json_safe(row.get("source_question_version")),
        })
    evidence_rows.sort(key=lambda row: (row["organism_id"], row["pmid"]))
    species_rows = []
    for row in pathogens.to_dict("records"):
        pathogen_id = str(row["pathogen_id"])
        group = str(row["pathogen_group"])
        species_rows.append({
            "organism_id": pathogen_id, "pathogen_id": pathogen_id,
            "organism_group": group, "accepted_name": str(row["accepted_name"]),
            "source_taxon_id": _json_safe(row.get("taxon_id")),
            "taxonomy_source": _json_safe(row.get("taxonomy_source")),
            "taxonomy_release": _json_safe(row.get("taxonomy_release")),
            "rank": _json_safe(row.get("rank")), "parent_taxon_id": _json_safe(row.get("parent_taxon_id")),
            "taxonomy_path": _json_safe(row.get("taxonomy_path")) or "",
            "synonyms": sorted(set(alias_map.get(pathogen_id, [])), key=str.casefold),
            "host_id": "ncbi:9606", "host_ids": ["ncbi:9606"],
            "feature_values": {}, "feature_counts": {},
            "n_pathogen_pmids": evidence_counts.get(pathogen_id, 0),
            "n_jev_accepted_pmids": evidence_counts.get(pathogen_id, 0),
            "n_pubmed_search_records": int(search_counts.get(pathogen_id, 0)),
            "n_abstracts_fetched": int(fetched_counts.get(pathogen_id, 0)),
            "x": None, "y": None, "pathogenic_interest_index": None,
            "evidence_file": "evidence.json", "provenance_status": "legacy_import",
        })
    species_rows.sort(key=lambda row: (row["organism_group"], row["accepted_name"].casefold()))

    host_rows = []
    for row in hosts.to_dict("records"):
        host_id = str(row["host_id"])
        host_rows.append({**{key: _json_safe(value) for key, value in row.items()},
                          "label": HOST_LABELS.get(host_id, str(row["accepted_name"]))})
    host_rows.sort(key=lambda row: (row["host_id"] != "ncbi:9606", row["accepted_name"]))

    data = {
        "schema_version": "v4-explorer-preview-v1",
        "groups": ["bacteria", "viruses", "fungi"],
        "hosts": host_rows,
        "species": species_rows,
        "features": {}, "feature_hierarchy": [],
        "feature_profiles": {group: [] for group in ("bacteria", "viruses", "fungi")},
        "coordinates_ready": False,
        "preview_status": "historical_human_evidence",
        "note": "Positive records are imported historical human disease evidence. Full abstracts are retained as coarse legacy passages. Other host cohorts have no assessed V4 evidence in this preview.",
    }
    site_root.mkdir(parents=True, exist_ok=True)
    with (site_root / "explorer-data.json.gz").open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as compressed:
            compressed.write(json.dumps(data, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    (site_root / "evidence.json").write_text(json.dumps(evidence_rows, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    catalog = {"catalog_version": "v4-empty-catalog-v1", "profile_version": "host-aware-profile-v1",
               "features": {}, "feature_hierarchy": [], "feature_profiles": data["feature_profiles"]}
    (site_root / "feature-catalog.json").write_text(json.dumps(catalog, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    return {"species": len(species_rows), "hosts": len(host_rows), "evidence_rows": len(evidence_rows),
            "supported_hosts": sorted({row["host_id"] for row in evidence_rows})}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging", type=Path, required=True)
    parser.add_argument("--site-root", type=Path, default=Path("sites/pathogen-atlas-explorer"))
    args = parser.parse_args()
    print(json.dumps(build_site_data(args.staging, args.site_root), indent=2))


if __name__ == "__main__":
    main()
