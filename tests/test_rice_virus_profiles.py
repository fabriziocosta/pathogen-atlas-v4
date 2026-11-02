import math
import gzip
import json
from pathlib import Path

import pandas as pd

from pathogen_atlas_v4.rice_virus_profiles import (
    FEATURES,
    LITERATURE_REVIEW_PMID,
    LITERATURE_VIRUSES,
    RICE_HOST_ID,
    add_rice_virus_profiles,
)


def test_rice_virus_literature_profiles_have_traits_and_fitted_coordinates_without_touching_human():
    pathogen_rows = [
        {"pathogen_id": f"virus:{index}", "pathogen_group": "viruses", "accepted_name": item["match"]}
        for index, item in enumerate(LITERATURE_VIRUSES)
        if item.get("match")
    ]
    pathogens = pd.DataFrame(pathogen_rows)
    human_pair = {"feature_values": {"human_feature": 0.75}, "x": 12.0, "y": -4.0}
    taxa = [{
        "organism_id": "virus:0", "pathogen_id": "virus:0", "organism_group": "viruses",
        "accepted_name": "Phytoreovirus alphaoryzae", "synonyms": [],
        "pair_profiles": {"ncbi:9606": dict(human_pair)},
    }]
    data = {
        "groups": ["viruses"], "feature_profile_counts": {}, "profile_availability_counts": {},
        "coordinate_host_ids": ["ncbi:9606"], "coordinate_contexts": {}, "stats": {},
    }

    updated, updated_taxa, feature_data, evidence = add_rice_virus_profiles(data, taxa, pathogens)

    rice_profiles = [row for row in updated_taxa if RICE_HOST_ID in (row.get("pair_profiles") or {})]
    assert len(LITERATURE_VIRUSES) == 20
    assert len(rice_profiles) == 20
    assert set(feature_data["features"]) == set(FEATURES)
    assert "viruses" in updated["groups"]
    assert updated["profile_availability_counts"][RICE_HOST_ID]["viruses"] >= 20
    assert {row["pmid"] for row in evidence if row.get("evidence_kind") == "literature_feature"} == {
        LITERATURE_REVIEW_PMID, "28695276", "34721366", "41649689"
    }
    for taxon in rice_profiles:
        pair = taxon["pair_profiles"][RICE_HOST_ID]
        assert taxon["synonyms"]
        assert math.isfinite(pair["x"]) and math.isfinite(pair["y"])
        assert pair["coordinate_source_host_id"] == RICE_HOST_ID
        assert pair["feature_profile_version"] == "rice-virus-literature-v1"
        assert any(value == 1.0 for value in pair["feature_values"].values())

    # Adding the rice layer must leave the pre-existing human profile untouched.
    human_after = next(row for row in updated_taxa if row["organism_id"] == "virus:0")["pair_profiles"]["ncbi:9606"]
    assert human_after == human_pair
    assert "ncbi:9606" in updated["coordinate_host_ids"]
    assert "ncbi:4530" in updated["coordinate_host_ids"]


def test_built_rice_layer_retains_every_human_profile_and_coordinate():
    root = Path(__file__).resolve().parents[1] / "sites/pathogen-atlas-explorer"
    with gzip.open(root / "source/human_pair_data.json.gz", "rt", encoding="utf-8") as stream:
        source = json.load(stream)
    columns = json.loads(gzip.decompress((root / "explorer-data.json.gz").read_bytes()))["species_columns"]
    pair_column = columns.index("pair_profiles")
    built = {}
    for group in ("bacteria", "fungi", "viruses"):
        with gzip.open(root / f"taxa/{group}.json.gz", "rt", encoding="utf-8") as stream:
            for row in json.load(stream):
                built[str(row[0])] = (row[pair_column] or {}).get("ncbi:9606", {})

    assert len(source["profiles"]) == sum(profile.get("x") is not None and profile.get("y") is not None
                                            for profile in source["profiles"])
    for profile in source["profiles"]:
        restored = built[str(profile["pathogen_id"])]
        assert restored["x"] == profile["x"]
        assert restored["y"] == profile["y"]
        assert restored["feature_values"] == profile["feature_values"]
