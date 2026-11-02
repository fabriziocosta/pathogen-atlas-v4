"""Merge title-only human virus profiles from a pinned Explorer data export."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import io
import json
import urllib.request
from pathlib import Path
from typing import Any

HUMAN_ID = "ncbi:9606"
SOURCE_RELEASE = "Pathogen Atlas V3"


def read_json_gzip(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return json.load(stream)


def write_json_gzip(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
            stream.write(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))


def fetch_shard(base_url: str, species: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    shard = str(species["evidence_file"])
    url = base_url.rstrip("/") + "/" + shard.lstrip("/")
    with urllib.request.urlopen(url, timeout=60) as response:
        content = response.read()
    records = json.loads(content)
    expected = int(species.get("n_pathogen_pmids") or 0)
    if len(records) != expected:
        raise ValueError(f"{species['accepted_name']} expected {expected} records, found {len(records)}")
    digest = hashlib.sha256(content).hexdigest()
    normalized = []
    for source in records:
        record = {key: value for key, value in source.items() if key not in {"abstract", "abstract_text"}}
        record.update(pathogen_id=str(record["organism_id"]), host_id=HUMAN_ID,
                      relationship_context="disease", provenance_status="legacy_import")
        if not str(record.get("title") or "").strip():
            raise ValueError(f"Legacy virus PMID {record.get('pmid')} has no title")
        normalized.append(record)
    return normalized, {"species_id": species["organism_id"], "source_shard": shard,
                        "sha256": digest, "records": len(normalized)}


def merge(site_root: Path, source_payload_path: Path, source_base_url: str) -> dict[str, Any]:
    payload_bytes = source_payload_path.read_bytes()
    with gzip.GzipFile(fileobj=io.BytesIO(payload_bytes)) as stream:
        source = json.load(stream)
    if source.get("release") != SOURCE_RELEASE:
        raise ValueError(f"Expected {SOURCE_RELEASE}, got {source.get('release')!r}")
    virus_profiles = [row for row in source["species"] if row.get("organism_group") == "viruses"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        fetched = list(pool.map(lambda row: fetch_shard(source_base_url, row), virus_profiles))
    virus_evidence = [record for records, _ in fetched for record in records]
    shard_manifest = [manifest for _, manifest in fetched]

    payload_path = site_root / "explorer-data.json.gz"
    evidence_path = site_root / "evidence.json.gz"
    current = read_json_gzip(payload_path)
    base_species = [row for row in current["species"] if row.get("organism_group") != "viruses"]
    known_ids = {str(row["organism_id"]) for row in base_species}
    if any(str(row["organism_id"]) in known_ids for row in virus_profiles):
        raise ValueError("Virus profile ID collides with the existing human baseline")
    for species in virus_profiles:
        species.update(pathogen_id=str(species["organism_id"]), host_id=HUMAN_ID,
                       host_ids=[HUMAN_ID], host_group="Animals",
                       profile_provenance_status="legacy_import", provenance_status="legacy_import",
                       feature_profile_version=source.get("stats", {}).get("feature_question_set_version"),
                       evidence_file="evidence.json.gz")
    current["species"] = sorted(base_species + virus_profiles,
                                key=lambda row: (row["organism_group"], row["accepted_name"].casefold()))
    current["groups"] = source["groups"]
    current["coordinates_ready"] = bool(source.get("coordinates_ready"))
    current["stats"] = source.get("stats", {})
    current["legacy_human_parity"].update(
        source_release=SOURCE_RELEASE,
        feature_profile_version=source.get("stats", {}).get("feature_question_set_version"),
        profile_species=len(current["species"]),
        accepted_species_pmid_records=sum(int(row.get("n_pathogen_pmids") or 0) for row in current["species"]),
        virus_pathogen_species=len(virus_profiles),
        virus_accepted_species_pmid_records=len(virus_evidence),
    )
    base_evidence = [row for row in read_json_gzip(evidence_path) if row.get("organism_group") != "viruses"]
    combined = sorted(base_evidence + virus_evidence,
                      key=lambda row: (str(row.get("organism_id", "")), str(row.get("pmid", ""))))
    write_json_gzip(payload_path, current)
    write_json_gzip(evidence_path, combined)
    manifest = {
        "source_release": SOURCE_RELEASE,
        "source_payload_sha256": hashlib.sha256(payload_bytes).hexdigest(),
        "source_payload_path": source_payload_path.name,
        "host_id": HUMAN_ID,
        "profile_count": len(virus_profiles),
        "accepted_pmid_records": len(virus_evidence),
        "titles_only": True,
        "evidence_shards": sorted(shard_manifest, key=lambda row: row["species_id"]),
    }
    (site_root / "legacy-virus-import-manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {"profiles": len(virus_profiles), "evidence_records": len(virus_evidence),
            "species_total": len(current["species"]), "evidence_total": len(combined)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-payload", type=Path, required=True,
                        help="Compressed Explorer data payload from the source release")
    parser.add_argument("--source-base-url", default="https://fabriziocosta.github.io/pathogen-atlas-v3/")
    parser.add_argument("--site-root", type=Path, default=Path("sites/pathogen-atlas-explorer"))
    args = parser.parse_args()
    print(json.dumps(merge(args.site_root, args.source_payload, args.source_base_url), indent=2))


if __name__ == "__main__":
    main()
