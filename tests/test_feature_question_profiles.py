from pathlib import Path

import yaml

from pathogen_atlas_v4.features import resolve_profile


ROOT = Path(__file__).resolve().parents[1]


def _profile(catalogue, host_group, taxon_id, pathogen_group):
    return resolve_profile(
        catalogue,
        pathogen_group=pathogen_group,
        host_group=host_group,
        host_taxon_id=taxon_id,
        relationship_context="disease",
        profile_version=catalogue["profile_version"],
    )


def test_trait_filtering_works_across_every_configured_host_and_pathogen_group():
    catalogue = yaml.safe_load((ROOT / "configs/features.yaml").read_text(encoding="utf-8"))
    host_taxa = {"human": "9606", "livestock": "9913", "crop": "4530"}
    pathogen_groups = [key.split(":", 1)[1] for key in catalogue["modules"]
                       if key.startswith("pathogen_group:")]

    for host_group, taxon_id in host_taxa.items():
        for pathogen_group in pathogen_groups:
            profile = _profile(catalogue, host_group, taxon_id, pathogen_group)
            snapshot = profile["trait_snapshot"]
            feature_ids = [feature["feature_id"] for feature in profile["ordered_features"]]
            assert len(feature_ids) == len(set(feature_ids))
            for feature in profile["ordered_features"]:
                applicability = feature.get("applicability", {})
                for axis in ("pathogen", "host", "pair"):
                    rules = applicability.get(f"{axis}_traits", {})
                    available = set(snapshot[axis])
                    assert set(rules.get("all_of", [])) <= available
                    assert not rules.get("any_of") or set(rules["any_of"]) & available


def test_shared_host_type_questions_and_rice_specific_questions_resolve_by_traits():
    catalogue = yaml.safe_load((ROOT / "configs/features.yaml").read_text(encoding="utf-8"))
    rice = {feature["feature_id"] for feature in _profile(catalogue, "crop", "4530", "viruses")["ordered_features"]}
    wheat = {feature["feature_id"] for feature in _profile(catalogue, "crop", "4565", "viruses")["ordered_features"]}
    human = {feature["feature_id"] for feature in _profile(catalogue, "human", "9606", "viruses")["ordered_features"]}
    livestock = {feature["feature_id"] for feature in _profile(catalogue, "livestock", "9913", "viruses")["ordered_features"]}

    assert "plant_mosaic" in rice & wheat
    assert "rice_stunting" in rice and "rice_stunting" not in wheat | human | livestock
    assert "animal_respiratory_disease" in human & livestock
    assert "animal_respiratory_disease" not in rice | wheat
    assert "virus_genome_type" in rice & wheat & human & livestock


def test_applicability_column_describes_combined_host_and_pathogen_traits():
    catalogue = yaml.safe_load((ROOT / "configs/features.yaml").read_text(encoding="utf-8"))
    cattle_fungi = _profile(catalogue, "livestock", "9913", "fungi")
    selected = {feature["feature_id"]: feature for feature in cattle_fungi["ordered_features"]}
    assert selected["fungal_growth_form"]["applies_to_label"] == "All fungi"
    assert selected["animal_respiratory_disease"]["applies_to_label"] == "All animal pathogens"
    assert selected["pair_device_association"]["applies_to_label"] == "All animal pathogens"
    cattle_viruses = _profile(catalogue, "livestock", "9913", "viruses")
    virus_features = {feature["feature_id"]: feature for feature in cattle_viruses["ordered_features"]}
    assert virus_features["virus_persistence"]["applies_to_label"] == "Animal viruses"
