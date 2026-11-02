"""Build V4 Explorer assets with a host-tagged, provenance-preserving human import."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Any

import pandas as pd

HUMAN_ID = "ncbi:9606"
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
    "ncbi:9606": "Animals",
    "ncbi:9913": "Animals",
    "ncbi:9940": "Animals",
    "ncbi:9031": "Animals",
    "ncbi:9823": "Animals",
    "ncbi:4565": "Plants",
    "ncbi:3847": "Plants",
    "ncbi:4577": "Plants",
    "ncbi:4530": "Plants",
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


def _load_legacy_profiles(export_root: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    payload_path = export_root / "explorer-data.json.gz"
    with gzip.open(payload_path, "rt", encoding="utf-8") as stream:
        payload = json.load(stream)
    if not str(payload.get("release", "")).endswith("V3"):
        raise ValueError(f"Expected a pinned legacy human Explorer export at {payload_path}")

    species = []
    evidence = []
    for source_row in payload["species"]:
        # The prior Explorer's released profiles all describe human disease.
        row = {key: _json_safe(value) for key, value in source_row.items()}
        pathogen_id = str(row["organism_id"])
        row.update({
            "pathogen_id": pathogen_id,
            "host_id": HUMAN_ID,
            "host_ids": [HUMAN_ID],
            "host_group": "Animals",
            "profile_provenance_status": "legacy_import",
            "provenance_status": "legacy_import",
            "feature_profile_version": payload.get("stats", {}).get("feature_question_set_version"),
            "evidence_file": "evidence.json.gz",
        })
        evidence_file = row.get("evidence_file_source") or source_row.get("evidence_file")
        if not evidence_file:
            raise ValueError(f"Legacy profile has no evidence shard: {pathogen_id}")
        shard_path = (export_root / evidence_file).resolve()
        if export_root.resolve() not in shard_path.parents:
            raise ValueError(f"Evidence shard escapes the export root: {evidence_file}")
        records = json.loads(shard_path.read_text(encoding="utf-8"))
        expected = int(row.get("n_pathogen_pmids") or 0)
        if len(records) != expected:
            raise ValueError(f"Legacy profile {pathogen_id} expects {expected} evidence records, found {len(records)}")
        species.append(row)
        evidence.extend(records)
    return payload, species, evidence


def build_site_data(staging: Path, site_root: Path, legacy_export: Path) -> dict[str, Any]:
    pathogens = pd.read_parquet(staging / "taxonomy/pathogens.parquet")
    pathogen_ids = set(pathogens["pathogen_id"].astype(str))
    legacy_payload, species_rows, legacy_evidence = _load_legacy_profiles(legacy_export)
    missing_pathogens = sorted({row["pathogen_id"] for row in species_rows} - pathogen_ids)
    if missing_pathogens:
        raise ValueError(f"Legacy human profiles have unresolved V4 taxonomy IDs: {missing_pathogens[:10]}")

    evidence_rows = []
    unresolved_titles = []
    for record in legacy_evidence:
        pathogen_id, pmid = str(record["organism_id"]), str(record["pmid"])
        title = str(record.get("title") or "").strip()
        if not title:
            unresolved_titles.append({"pathogen_id": pathogen_id, "pmid": pmid})
            continue
        record_metadata = {key: _json_safe(value) for key, value in record.items()
                           if key not in {"abstract", "abstract_text", "evidence_span"}}
        evidence_rows.append({
            **record_metadata,
            "pathogen_id": pathogen_id,
            "host_id": HUMAN_ID,
            "relationship_context": "disease",
            "title": title,
            "publication_date": _json_safe(record.get("publication_date")),
            "publication_year": _json_safe(record.get("publication_year")),
            "doi": _json_safe(record.get("doi")),
            "provenance_status": "legacy_import",
            "accepted_pathogen_evidence": True,
        })
    if unresolved_titles:
        raise ValueError(f"{len(unresolved_titles)} legacy human records have no title; first rows: {unresolved_titles[:10]}")

    expected_evidence = int(legacy_payload.get("stats", {}).get("accepted_species_pmid_records", 0))
    if len(evidence_rows) != expected_evidence:
        raise ValueError(f"Expected {expected_evidence} legacy human evidence records, built {len(evidence_rows)}")
    evidence_rows.sort(key=lambda row: (row["pathogen_id"], row["pmid"]))
    species_rows.sort(key=lambda row: (row["organism_group"], row["accepted_name"].casefold()))

    hosts = pd.read_parquet(staging / "taxonomy/hosts.parquet")
    host_rows = []
    for row in hosts.to_dict("records"):
        host_id = str(row["host_id"])
        host_rows.append({**{key: _json_safe(value) for key, value in row.items()},
                          "label": HOST_LABELS.get(host_id, str(row["accepted_name"])),
                          "host_group": HOST_GROUPS.get(host_id, "Animals")})
    host_rows.sort(key=lambda row: (row["host_group"] != "Animals", row["accepted_name"]))

    data = {
        "schema_version": "v4-explorer-human-parity-v2-host-selection",
        "groups": sorted(set(legacy_payload.get("groups", [])) | {"parasites", "oomycetes"}),
        "hosts": host_rows,
        "species": species_rows,
        "features": legacy_payload.get("features", {}),
        "feature_hierarchy": legacy_payload.get("feature_hierarchy", []),
        "feature_profiles": legacy_payload.get("feature_profiles", {}),
        "coordinates_ready": bool(legacy_payload.get("coordinates_ready")),
        "stats": legacy_payload.get("stats", {}),
        "legacy_human_parity": {
            "source_release": legacy_payload.get("release"),
            "feature_profile_version": legacy_payload.get("stats", {}).get("feature_question_set_version"),
            "profile_species": len(species_rows),
            "accepted_species_pmid_records": len(evidence_rows),
            "provenance_status": "legacy_import",
        },
    }
    site_root.mkdir(parents=True, exist_ok=True)
    with (site_root / "explorer-data.json.gz").open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as compressed:
            compressed.write(json.dumps(data, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    with (site_root / "evidence.json.gz").open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as compressed:
            compressed.write(json.dumps(evidence_rows, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    (site_root / "evidence.json").unlink(missing_ok=True)
    catalog = {
        "catalog_version": "legacy-human-feature-catalog-v3.2",
        "profile_version": data["legacy_human_parity"]["feature_profile_version"],
        "features": data["features"],
        "feature_hierarchy": data["feature_hierarchy"],
        "feature_profiles": data["feature_profiles"],
    }
    (site_root / "feature-catalog.json").write_text(
        json.dumps(catalog, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    return {
        "profile_species": len(species_rows),
        "hosts": len(host_rows),
        "accepted_species_pmid_records": len(evidence_rows),
        "supported_hosts": sorted({row["host_id"] for row in evidence_rows}),
        "feature_count": len(data["features"]),
        "coordinates_ready": data["coordinates_ready"],
    }


def update_host_catalog(staging: Path, site_root: Path) -> dict[str, Any]:
    """Add pinned host choices while preserving the released human evidence payload."""
    payload_path = site_root / "explorer-data.json.gz"
    with gzip.open(payload_path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    if not str(data.get("schema_version", "")).startswith("v4-explorer-human-parity-"):
        raise ValueError(f"Expected an existing V4 human-parity Explorer payload at {payload_path}")
    hosts = pd.read_parquet(staging / "taxonomy/hosts.parquet")
    host_rows = []
    for row in hosts.to_dict("records"):
        host_id = str(row["host_id"])
        host_rows.append({**{key: _json_safe(value) for key, value in row.items()},
                          "label": HOST_LABELS.get(host_id, str(row["accepted_name"])),
                          "host_group": HOST_GROUPS.get(host_id, "Animals")})
    host_rows.sort(key=lambda row: (row["host_group"] != "Animals", row["accepted_name"].casefold()))
    data["hosts"] = host_rows
    data["schema_version"] = "v4-explorer-human-parity-v2-host-selection"
    with payload_path.open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as compressed:
            compressed.write(json.dumps(data, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    return {"hosts": len(host_rows), "host_ids": [row["host_id"] for row in host_rows],
            "profile_species_preserved": len(data.get("species", [])),
            "accepted_species_pmid_records_preserved": int(data.get("stats", {}).get("accepted_species_pmid_records", 0)),
            "payload": str(payload_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging", type=Path, required=True)
    parser.add_argument("--legacy-human-export", type=Path,
                        help="Pinned legacy human Explorer bundle (explorer-data.json.gz plus evidence shards)")
    parser.add_argument("--site-root", type=Path, default=Path("sites/pathogen-atlas-explorer"))
    parser.add_argument("--hosts-only", action="store_true",
                        help="Update host choices while preserving the existing human profile/evidence payload")
    args = parser.parse_args()
    if args.hosts_only:
        result = update_host_catalog(args.staging, args.site_root)
    else:
        if args.legacy_human_export is None:
            parser.error("--legacy-human-export is required unless --hosts-only is set")
        result = build_site_data(args.staging, args.site_root, args.legacy_human_export)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
