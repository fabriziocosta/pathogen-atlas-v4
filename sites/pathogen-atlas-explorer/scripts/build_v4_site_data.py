"""Build the Explorer from native V4 run outputs."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
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
}
HOST_GROUPS = {
    "ncbi:9606": "Animals", "ncbi:9913": "Animals", "ncbi:9940": "Animals",
    "ncbi:9031": "Animals", "ncbi:9823": "Animals", "ncbi:4565": "Plants",
    "ncbi:3847": "Plants", "ncbi:4577": "Plants", "ncbi:4530": "Plants",
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


def build_site_data(run_dir: Path, staging_root: Path, site_root: Path,
                    *, dist_root: Path | None = None) -> dict[str, Any]:
    """Materialize host-specific accepted profiles and evidence for one completed run."""
    run = Path(run_dir)
    staging = Path(staging_root)
    site = Path(site_root)
    resolved = json.loads((run / "resolved_config.json").read_text(encoding="utf-8"))
    selected_pathogen_ids = set(map(str, resolved["pathogen_ids"]))
    selected_host_ids = set(map(str, resolved["host_ids"]))
    pathogens = pd.read_parquet(run / "taxonomy/pathogens.parquet")
    pathogens["pathogen_id"] = pathogens.pathogen_id.astype(str)
    pathogens = pathogens.loc[
        pathogens.pathogen_id.isin(selected_pathogen_ids)
        & pathogens["rank"].astype(str).str.casefold().eq("species")
    ].copy()
    pairs = pd.read_parquet(run / "pathogen_host_selection.parquet")
    pairs["pathogen_id"] = pairs.pathogen_id.astype(str)
    pairs["host_id"] = pairs.host_id.astype(str)
    pairs = pairs.loc[pairs.pathogen_id.isin(selected_pathogen_ids)
                      & pairs.host_id.isin(selected_host_ids)].copy()
    if pairs.duplicated(["pathogen_id", "host_id"]).any():
        raise ValueError("Pair selection must contain one row per pathogen and host")
    if len(pairs) != len(selected_pathogen_ids) * len(selected_host_ids):
        raise ValueError("Pair selection does not cover the complete configured host/pathogen scope")

    pii_path = run / "release/pathogenic_interest_index.parquet"
    pii = pd.read_parquet(pii_path) if pii_path.exists() else pd.DataFrame()
    if not pii.empty:
        pii["pathogen_id"] = pii.pathogen_id.astype(str)
        pii["host_id"] = pii.host_id.astype(str)
    articles = pd.read_parquet(staging / "literature/articles.parquet")
    articles["pmid"] = articles.pmid.astype(str)
    article_by_pmid = articles.drop_duplicates("pmid", keep="last").set_index("pmid")
    evidence = pd.read_parquet(run / "evidence/pathogen_host_evidence.parquet")
    evidence["pathogen_id"] = evidence.pathogen_id.astype(str)
    evidence["host_id"] = evidence.host_id.astype(str)
    evidence["pmid"] = evidence.pmid.astype(str)
    evidence = evidence.loc[evidence.decision.astype(str).eq("supported")].copy()
    evidence_by_pair = {key: group for key, group in evidence.groupby(["pathogen_id", "host_id"], sort=False)}

    alias_path = run / "taxonomy/pathogen_aliases.parquet"
    aliases_by_id: dict[str, list[str]] = {}
    if alias_path.exists():
        aliases = pd.read_parquet(alias_path, columns=["pathogen_id", "name", "name_type"])
        aliases = aliases.loc[aliases.name_type.astype(str).str.casefold().eq("synonym")]
        aliases_by_id = {str(key): sorted(set(group.name.dropna().astype(str)))
                         for key, group in aliases.groupby("pathogen_id", sort=False)}

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
    for row in host_table.loc[host_table.host_id.astype(str).isin(selected_host_ids)].to_dict("records"):
        host_id = str(row["host_id"])
        hosts.append({**{key: _json_safe(value) for key, value in row.items()},
                      "label": HOST_LABELS.get(host_id, str(row["accepted_name"])),
                      "host_group": HOST_GROUPS.get(host_id, str(row.get("host_group") or "Animals"))})

    species_rows: list[dict[str, Any]] = []
    evidence_rows: list[dict[str, Any]] = []
    for key, pair in pair_by_key.items():
        pathogen_id, host_id = key
        taxonomy = taxonomy_by_id.get(pathogen_id)
        if taxonomy is None:
            raise ValueError(f"Pair selection references pathogen absent from run taxonomy: {pathogen_id}")
        host = next((item for item in hosts if str(item["host_id"]) == host_id), None)
        if host is None:
            raise ValueError(f"Pair selection references host absent from staging taxonomy: {host_id}")
        metric = pii_by_key.get(key, {})
        matched_evidence = evidence_by_pair.get(key, pd.DataFrame())
        positive_count = int(pair.get("n_positive_records") or 0)
        if len(matched_evidence) != positive_count:
            raise ValueError(f"Pair evidence count mismatch for {pathogen_id}, {host_id}")
        item = {
            "organism_id": pathogen_id,
            "pathogen_id": pathogen_id,
            "organism_group": str(taxonomy["pathogen_group"]),
            "accepted_name": str(taxonomy["accepted_name"]),
            "source_taxon_id": str(taxonomy.get("taxon_id") or ""),
            "taxonomy_source": str(taxonomy.get("taxonomy_source") or ""),
            "taxonomy_release": str(taxonomy.get("taxonomy_release") or ""),
            "rank": str(taxonomy.get("rank") or ""),
            "parent_taxon_id": _json_safe(taxonomy.get("parent_taxon_id")),
            "taxonomy_path": str(taxonomy.get("taxonomy_path") or ""),
            "taxonomy": _taxonomy(taxonomy.get("taxonomy_path")),
            "mapping_status": str(taxonomy.get("mapping_status") or ""),
            "host_id": host_id,
            "host_ids": [host_id],
            "host_group": host["host_group"],
            "relationship_context": str(pair.get("relationship_context") or "disease"),
            "pair_selection_status": str(pair.get("selection_status") or "incomplete"),
            "n_pathogen_pmids": positive_count,
            "n_jev_accepted_pmids": positive_count,
            "n_positive_records": positive_count,
            "n_sampled_records": int(pair.get("n_sampled_records") or 0),
            "n_abstracts_fetched": int(pair.get("n_sampled_records_with_abstract") or 0),
            "n_abstracts_screened": int(pair.get("n_completed_screenings") or 0),
            "n_pair_search_records": int(pair.get("n_pair_search_records") or 0),
            "n_pubmed_search_records": _json_safe(metric.get("n_pubmed_search_records")),
            "estimated_pathogen_publications": _json_safe(metric.get("estimated_pathogen_publications")),
            "estimated_pathogen_publications_index": _json_safe(metric.get("estimated_pathogen_publications_index")),
            "pathogenicity_fraction": _json_safe(metric.get("pathogenicity_fraction")),
            "pathogenicity_fraction_index": _json_safe(metric.get("pathogenicity_fraction_index")),
            "pathogenic_interest_index": _json_safe(metric.get("pathogenic_interest_index")),
            "index_version": metric.get("index_version"),
            "evidence_file": "evidence.json.gz",
            "feature_values": {},
            "feature_counts": {},
            "feature_profile_version": "host-aware-profile-v1",
            "profile_provenance_status": "native_v4",
            "provenance_status": "native_v4",
            "synonyms": aliases_by_id.get(pathogen_id, []),
            "wikipedia_url": wikipedia.get(pathogen_id, {}).get("wikipedia_url") or None,
            "gbif_url": gbif.get(pathogen_id, {}).get("gbif_url") or None,
            "x": None,
            "y": None,
        }
        species_rows.append(item)
        for record in matched_evidence.to_dict("records"):
            article = article_by_pmid.loc[str(record["pmid"])] if str(record["pmid"]) in article_by_pmid.index else {}
            if isinstance(article, pd.DataFrame):
                article = article.iloc[-1]
            evidence_rows.append({
                "organism_id": pathogen_id,
                "organism_group": item["organism_group"],
                "accepted_name": item["accepted_name"],
                "pathogen_id": pathogen_id,
                "host_id": host_id,
                "relationship_context": item["relationship_context"],
                "pmid": str(record["pmid"]),
                "title": str(article.get("title") or ""),
                "publication_date": _json_safe(article.get("publication_date")),
                "publication_year": _json_safe(article.get("publication_year")),
                "doi": str(article.get("doi") or ""),
                "p_yes": _json_safe(record.get("p_yes")),
                "decision": str(record["decision"]),
                "accepted_pathogen_evidence": True,
                "provenance_status": str(record.get("provenance_status") or "native_v4"),
                "abstract_hash": str(record.get("abstract_hash") or ""),
            })

    species_rows.sort(key=lambda row: (row["organism_group"], row["accepted_name"].casefold(), row["host_id"]))
    evidence_rows.sort(key=lambda row: (row["pathogen_id"], row["host_id"], row["pmid"]))
    groups = sorted({row["organism_group"] for row in species_rows})
    decisions = yaml.safe_load((Path(__file__).parents[3] / "configs/prompts.yaml").read_text(encoding="utf-8"))
    screen_version = str(decisions["pathogen_screen"]["question_version"])
    supported_pairs = int(pairs.selection_status.astype(str).eq("accepted").sum())
    stats = {
        "run_id": str(resolved.get("run_id", run.name)),
        "question_version": screen_version,
        "yes_probability_threshold": 0.7,
        "candidate_pairs": len(pairs),
        "pubmed_hit_pairs": int(pd.to_numeric(pairs.n_pair_search_records, errors="coerce").fillna(0).gt(0).sum()),
        "accepted_pairs": supported_pairs,
        "sampled_records": int(pairs.n_sampled_records.sum()),
        "supported_pmids": len(evidence_rows),
    }
    feature_data = {"catalog_version": "native-v4-feature-catalog-v1", "profile_version": "host-aware-profile-v1",
                    "features": {}, "feature_hierarchy": [], "feature_profiles": {group: [] for group in groups}}
    data = {
        "schema_version": "v4-native-pair-screen-v1",
        "groups": groups,
        "hosts": hosts,
        "species": species_rows,
        "features": feature_data["features"],
        "feature_hierarchy": feature_data["feature_hierarchy"],
        "feature_profiles": feature_data["feature_profiles"],
        "coordinates_ready": False,
        "stats": stats,
    }
    site.mkdir(parents=True, exist_ok=True)
    for name, value in (("explorer-data.json.gz", data), ("evidence.json.gz", evidence_rows)):
        with (site / name).open("wb") as output:
            with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as compressed:
                compressed.write(json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode())
    (site / "feature-catalog.json").write_text(json.dumps(feature_data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    if dist_root is not None:
        _copy_static_site(site, Path(dist_root))
    return {**stats, "species_profiles": len(species_rows), "evidence_records": len(evidence_rows),
            "hosts": [host["host_id"] for host in hosts], "payload": str(site / "explorer-data.json.gz")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True, help="Completed V4 run directory")
    parser.add_argument("--staging-root", type=Path, required=True, help="Pinned shared taxonomy/literature staging directory")
    parser.add_argument("--site-root", type=Path, default=Path("sites/pathogen-atlas-explorer"))
    parser.add_argument("--dist-root", type=Path, default=Path("sites/pathogen-atlas-explorer/dist"))
    args = parser.parse_args()
    print(json.dumps(build_site_data(args.run_dir, args.staging_root, args.site_root, dist_root=args.dist_root), indent=2))


if __name__ == "__main__":
    main()
