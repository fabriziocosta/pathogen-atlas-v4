import gzip
import hashlib
import json
from pathlib import Path

import pytest

from pathogen_atlas_v4.explorer_release import validate_explorer_release


def _write_gzip(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(json.dumps(value, separators=(",", ":")).encode()))


def _site(tmp_path: Path, *, feature_ids_yes=None, coordinates=(0.2, -0.3), n_yes=1, pmid="12345", include_document=True):
    site = tmp_path / "site"
    site.mkdir()
    (site / "index.html").write_text(
        '<span class="title-version">v4.0.99</span> explorer-data.json.gz feature-catalog.json '
        "const separator=path.includes('?')?'&':'?' DATA.evidence_version fetch(url,{cache:'no-store'}) fetch('explorer-data.json.gz',{cache:'no-store'})",
        encoding="utf-8",
    )
    metadata = {
        "groups": ["fungi"], "hosts": [{"host_id": "host:1"}],
        "species_columns": ["organism_id", "accepted_name", "pair_profiles"],
        "stats": {"pathogen_species_by_group": {"fungi": 1}},
        "coordinate_host_ids": ["host:1"], "coordinates_ready": True,
    }
    _write_gzip(site / "explorer-data.json.gz", metadata)
    (site / "feature-catalog.json").write_text(json.dumps({
        "features": {"disease": {"label": "Disease"}},
        "feature_hierarchy": [{"categories": [{"feature_ids": ["disease"]}]}],
    }))
    docs = [{
        "organism_id": "fungi:1", "pathogen_id": "fungi:1", "host_id": "host:1", "pmid": pmid,
        "title": "Evidence title", "accepted_pathogen_evidence": True,
        "feature_ids_yes": feature_ids_yes if feature_ids_yes is not None else ["disease"],
    }] if include_document else []
    _write_gzip(site / "evidence.json.gz", docs)
    x, y = coordinates
    map_y = y if y is not None else 0.0
    profile = {
        "n_jev_accepted_pmids": 1, "n_abstracts_screened": 1,
        "feature_counts": {"disease": {"n_yes": n_yes}},
        "evidence_file": "evidence.json.gz", "x": x, "y": y,
        "embedding_coordinates": {"IMSTE": {"x": x, "y": map_y}, "UMAP": {"x": x + 1, "y": map_y + 1}},
    }
    _write_gzip(site / "taxa/fungi.json.gz", [["fungi:1", "Fungus", {"host:1": profile}]])
    return site


def test_release_validator_accepts_documents_and_2d_coordinates(tmp_path):
    summary = validate_explorer_release(_site(tmp_path), expected_version="4.0.99")
    assert summary == {
        "species": 1, "pair_profiles": 1, "evidence_documents": 1,
        "coordinate_hosts": 1, "features": 1,
    }


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"include_document": False}, "reports 1 accepted PMIDs but has no documents"),
        ({"coordinates": (0.2, None)}, "only one 2D coordinate"),
        ({"pmid": "not-a-pmid"}, "no numeric PMID"),
    ],
)
def test_release_validator_rejects_incomplete_payloads(tmp_path, kwargs, message):
    with pytest.raises(ValueError, match=message):
        validate_explorer_release(_site(tmp_path, **kwargs), expected_version="4.0.99")


def test_release_validator_rejects_stale_displayed_version(tmp_path):
    with pytest.raises(ValueError, match="does not display expected release version"):
        validate_explorer_release(_site(tmp_path), expected_version="4.1.0")


@pytest.mark.parametrize(
    ("metadata_key", "metadata_value", "message"),
    [
        ("evidence_version", "stale", "does not match explorer-data.json.gz evidence_version"),
        ("taxa_shard_versions", {"fungi": "stale"}, "does not match its declared version"),
    ],
)
def test_release_validator_rejects_mixed_payload_versions(tmp_path, metadata_key, metadata_value, message):
    site = _site(tmp_path)
    metadata_path = site / "explorer-data.json.gz"
    metadata = json.loads(gzip.decompress(metadata_path.read_bytes()))
    metadata[metadata_key] = metadata_value
    metadata_path.write_bytes(gzip.compress(json.dumps(metadata, separators=(",", ":")).encode()))
    with pytest.raises(ValueError, match=message):
        validate_explorer_release(site, expected_version="4.0.99")


def test_release_validator_accepts_matching_evidence_digest(tmp_path):
    site = _site(tmp_path)
    evidence = json.loads(gzip.decompress((site / "evidence.json.gz").read_bytes()))
    digest = hashlib.sha256(
        json.dumps(evidence, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()[:16]
    metadata_path = site / "explorer-data.json.gz"
    metadata = json.loads(gzip.decompress(metadata_path.read_bytes()))
    metadata["evidence_version"] = digest
    metadata_path.write_bytes(gzip.compress(json.dumps(metadata, separators=(",", ":")).encode()))
    assert validate_explorer_release(site, expected_version="4.0.99")["evidence_documents"] == 1
