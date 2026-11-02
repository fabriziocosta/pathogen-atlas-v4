import gzip
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "sites/pathogen-atlas-explorer"


def _read_gzip_json(path: Path):
    return json.loads(gzip.decompress(path.read_bytes()))


def _profile_availability_counts(metadata: dict) -> tuple[dict, int]:
    columns = metadata["species_columns"]
    pair_index = columns.index("pair_profiles")
    counts = {
        host["host_id"]: {group: 0 for group in metadata["groups"]}
        for host in metadata["hosts"]
    }
    rice_virus_coordinates = 0

    for group in metadata["groups"]:
        rows = _read_gzip_json(SITE / "dist/taxa" / f"{group}.json.gz")
        for row in rows:
            for host_id, profile in (row[pair_index] or {}).items():
                values = profile.get("feature_values") or {}
                has_features = any(value is not None for value in values.values())
                xy = (profile.get("x"), profile.get("y"))
                has_coordinates = all(value is not None and math.isfinite(float(value)) for value in xy)
                if has_features or has_coordinates:
                    counts.setdefault(host_id, {}).setdefault(group, 0)
                    counts[host_id][group] += 1
                if host_id == "ncbi:4530" and group == "viruses" and has_coordinates:
                    rice_virus_coordinates += 1

    return counts, rice_virus_coordinates


def test_profile_availability_includes_rice_virus_coordinates():
    metadata = _read_gzip_json(SITE / "explorer-data.json.gz")
    deployed_metadata = _read_gzip_json(SITE / "dist/explorer-data.json.gz")
    expected, rice_virus_coordinates = _profile_availability_counts(metadata)

    assert metadata["profile_availability_counts"] == expected
    assert deployed_metadata["profile_availability_counts"] == expected
    assert rice_virus_coordinates > 0
    assert metadata["profile_availability_counts"]["ncbi:4530"]["viruses"] >= rice_virus_coordinates

    explorer = (SITE / "index.html").read_text(encoding="utf-8")
    assert "DATA.profile_availability_counts||DATA.feature_profile_counts" in explorer
