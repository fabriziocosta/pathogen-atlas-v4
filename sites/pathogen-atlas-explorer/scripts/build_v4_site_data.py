"""Build the full atlas from taxonomy and completed host-pair releases."""

from __future__ import annotations

import argparse
import gzip
import json
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
    if (site_root / "taxa").exists():
        shutil.copytree(site_root / "taxa", dist_root / "taxa", dirs_exist_ok=True)


def build_site_data(run_dir: Path, staging_root: Path, site_root: Path,
                    *, dist_root: Path | None = None) -> dict[str, Any]:
    """Materialize the full taxonomy atlas and overlay completed host-pair releases."""
    run = Path(run_dir)
    staging = Path(staging_root)
    site = Path(site_root)
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
    additions, added_aliases, _ = load_ncbi_pathogen_additions(
        project_root / "data/source_taxonomy/ncbi",
        ["oomycetes", "parasites"], release=release, source_version=release,
    )
    additions = additions.loc[additions["rank"].astype(str).str.casefold().eq("species")]
    pathogens = pd.concat([pathogens, additions], ignore_index=True, sort=False)
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
    article_by_pmid = articles.drop_duplicates("pmid", keep="last").set_index("pmid")
    evidence = pd.read_parquet(run / "evidence/pathogen_host_evidence.parquet")
    evidence["pathogen_id"] = evidence.pathogen_id.astype(str)
    evidence["host_id"] = evidence.host_id.astype(str)
    evidence["pmid"] = evidence.pmid.astype(str)
    evidence = evidence.loc[evidence.decision.astype(str).eq("supported")].copy()
    evidence_by_pair = {key: group for key, group in evidence.groupby(["pathogen_id", "host_id"], sort=False)}

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

    human_source_path = site / "source/human_pair_data.json.gz"
    with gzip.open(human_source_path, "rt", encoding="utf-8") as source_file:
        human_source = json.load(source_file)
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
        positive_count = int(pair.get("n_positive_records") or 0)
        if len(matched_evidence) != positive_count:
            raise ValueError(f"Pair evidence count mismatch for {pathogen_id}, {host_id}")
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
            evidence_rows.append({
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
            })

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
    decisions = yaml.safe_load((Path(__file__).parents[3] / "configs/prompts.yaml").read_text(encoding="utf-8"))
    screen_version = str(decisions["pathogen_screen"]["question_version"])
    supported_pairs = sum(profile.get("pair_selection_status") in {"accepted", "supported"}
                          for profile in pair_profiles_by_key.values())
    feature_profiles = {group: list(human_source.get("feature_profiles", {}).get(group, [])) for group in groups}
    feature_data = {
        "catalog_version": "atlas-feature-catalog-v1",
        "profile_version": "host-aware-profile-v1",
        "features": human_source.get("features", {}),
        "feature_hierarchy": human_source.get("feature_hierarchy", []),
        "feature_profiles": feature_profiles,
    }
    stats = {
        "run_id": str(resolved.get("run_id", run.name)),
        "screening_question_version": screen_version,
        "screening_yes_probability_threshold": 0.7,
        "pathogen_species": len(species_rows),
        "hosts": len(hosts),
        "possible_host_pathogen_pairs": len(species_rows) * len(hosts),
        "pair_profiles_available": len(pair_profiles_by_key),
        "supported_pairs": supported_pairs,
        "current_run_candidate_pairs": len(pairs),
        "current_run_pubmed_hit_pairs": int(pd.to_numeric(pairs.n_pair_search_records, errors="coerce").fillna(0).gt(0).sum()),
        "current_run_sampled_records": int(pairs.n_sampled_records.sum()),
        "human_supported_profiles": len(human_profiles),
        "human_supporting_pmids": len(human_evidence),
        "supported_pmids": len(evidence_rows),
    }
    species_columns = [
        "organism_id", "accepted_name", "source_taxon_id", "taxonomy_source",
        "taxonomy_release", "rank", "parent_taxon_id", "taxonomy_path",
        "mapping_status", "synonyms", "pair_profiles",
    ]
    data = {
        "schema_version": "v4-host-pathogen-atlas-v1",
        "groups": groups,
        "hosts": hosts,
        "species_columns": species_columns,
        "taxonomy_ranks": sorted({part.split("=", 1)[0] for row in species_rows
                                   for part in row["taxonomy_path"].split(";") if "=" in part}),
        "coordinates_ready": False,
        "stats": {**stats, "pathogen_species_by_group": {
            group: sum(row["organism_group"] == group for row in species_rows) for group in groups
        }},
    }
    group_rows = {
        group: [[row["organism_id"], row["accepted_name"], row["source_taxon_id"],
                 row["taxonomy_source"], row["taxonomy_release"], row["rank"],
                 row["parent_taxon_id"], row["taxonomy_path"], row["mapping_status"],
                 row["synonyms"], row.get("pair_profiles")]
                for row in species_rows if row["organism_group"] == group]
        for group in groups
    }
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
