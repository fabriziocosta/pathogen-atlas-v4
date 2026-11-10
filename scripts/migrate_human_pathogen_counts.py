#!/usr/bin/env python3
"""Attach validated pathogen-wide PubMed counts to legacy Human profiles."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "sites/pathogen-atlas-explorer/source/human_pair_data.json.gz"


def migrate(source_path: Path, counts_path: Path) -> int:
    with gzip.open(source_path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    counts = pd.read_parquet(counts_path)
    required = {"pathogen_id", "query_hash", "query_version", "n_pubmed_search_records", "status"}
    if not required.issubset(counts.columns):
        raise ValueError(f"Count table is missing columns: {sorted(required - set(counts.columns))}")
    counts["pathogen_id"] = counts.pathogen_id.astype(str)
    if counts.pathogen_id.duplicated().any():
        raise ValueError("Count table must contain one row per pathogen ID")
    if not counts.status.astype(str).eq("completed").all():
        raise ValueError("Every migrated pathogen count must have completed status")
    values = pd.to_numeric(counts.n_pubmed_search_records, errors="coerce")
    if values.isna().any() or values.lt(0).any() or not (values == values.astype("int64")).all():
        raise ValueError("Pathogen-wide counts must be nonnegative integers")
    counts["n_pubmed_search_records"] = values.astype("int64")

    profiles = data.get("profiles", [])
    profile_ids = [str(row.get("pathogen_id") or row.get("organism_id") or "") for row in profiles]
    count_ids = set(counts.pathogen_id)
    missing = sorted(set(profile_ids) - count_ids)
    extra = sorted(count_ids - set(profile_ids))
    if missing or extra:
        raise ValueError(f"Count coverage mismatch: missing={missing[:10]}, extra={extra[:10]}")
    if len(profile_ids) != len(set(profile_ids)):
        raise ValueError("Human source must contain one profile per pathogen ID")

    count_by_id = counts.set_index("pathogen_id").to_dict("index")
    for profile, pathogen_id in zip(profiles, profile_ids):
        count = count_by_id[pathogen_id]
        profile["n_pathogen_pubmed_records"] = int(count["n_pubmed_search_records"])
        profile["n_pathogen_pubmed_query_hash"] = str(count["query_hash"])
        profile["n_pathogen_pubmed_query_version"] = str(count["query_version"])
        # The legacy Human corpus stores the host–pathogen query total in
        # n_pubmed_search_records. The separate n_pathogen_pubmed_records
        # field above is the host-independent organism count.
        for field in (
            "estimated_pathogen_publications", "estimated_pathogen_publications_index",
            "pathogenicity_fraction", "pathogenicity_fraction_index", "pathogenic_interest_index",
        ):
            profile.setdefault(f"legacy_{field}", profile.get(field))
            profile[field] = None
        sampled = int(profile.get("n_abstracts_fetched") or profile.get("n_sampled_records") or 0)
        positive = int(profile.get("n_pathogen_pmids") or profile.get("n_jev_accepted_pmids") or 0)
        if sampled > 0:
            if positive > sampled:
                raise ValueError(f"Human sample has more positives than screened records for {pathogen_id}")
            profile["pathogenicity_fraction"] = positive / sampled
            pair_total = profile.get("n_pair_search_records")
            if pair_total is None:
                pair_total = profile.get("n_pubmed_search_records")
            if pair_total is not None:
                pair_total = int(pair_total)
                profile["n_pair_search_records"] = pair_total
                profile["estimated_pathogen_publications"] = positive / sampled * pair_total

    source_path = Path(source_path)
    source_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{source_path.name}.", suffix=".tmp", dir=source_path.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with gzip.open(temp_path, "wt", encoding="utf-8", compresslevel=9) as stream:
            json.dump(data, stream, ensure_ascii=False, separators=(",", ":"))
        temp_path.replace(source_path)
    finally:
        temp_path.unlink(missing_ok=True)
    return len(profiles)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--counts", type=Path, required=True,
                        help="Complete accepted_pathogen_publication_counts.parquet")
    args = parser.parse_args()
    count = migrate(args.source, args.counts)
    print(f"Migrated pathogen-wide counts for {count} Human profiles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
