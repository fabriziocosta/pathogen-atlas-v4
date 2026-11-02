"""Rice-virus literature traits and a rice-specific feature-space embedding."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

RICE_HOST_ID = "ncbi:4530"
LITERATURE_REVIEW_PMID = "36298813"
FEATURES: dict[str, dict[str, str]] = {
    "rice_stunting": {"label": "Stunting or dwarfing", "definition": "Stunting or dwarfing of infected rice plants is reported."},
    "rice_increased_tillering": {"label": "Increased tillering", "definition": "Increased tiller number is reported in infected rice plants."},
    "rice_reduced_tillering": {"label": "Reduced tillering", "definition": "Reduced tiller number is reported in infected rice plants."},
    "rice_leaf_yellowing": {"label": "Leaf yellowing or chlorosis", "definition": "Yellowing or chlorosis of rice leaves is reported."},
    "rice_leaf_striping_mosaic": {"label": "Leaf striping or mosaic", "definition": "Leaf striping or mosaic symptoms are reported."},
    "rice_leaf_distortion_galls": {"label": "Leaf distortion or galls", "definition": "Leaf curling, twisting, ragged margins, galls, or enations are reported."},
    "rice_planthopper_vector": {"label": "Planthopper vector", "definition": "Transmission by a planthopper vector is reported."},
    "rice_leafhopper_vector": {"label": "Leafhopper vector", "definition": "Transmission by a leafhopper vector is reported."},
    "rice_whitefly_aphid_vector": {"label": "Whitefly or aphid vector", "definition": "Transmission by a whitefly or aphid vector is reported."},
    "rice_soilborne_vector": {"label": "Soil-borne fungal vector", "definition": "Soil-borne transmission by Polymyxa graminis is reported."},
    "rice_beetle_vector": {"label": "Beetle vector", "definition": "Transmission by a chrysomelid beetle vector is reported."},
    "rice_symptomatic_plant_association": {"label": "Associated with symptomatic rice", "definition": "The virus was identified in rice plants reported to have disease symptoms."},
}

# Current ICTV species names are resolved from the pinned taxonomy. The three
# unclassified literature names use stable local IDs and retain their tentative
# status rather than being presented as ICTV species.
LITERATURE_VIRUSES: list[dict[str, Any]] = [
    {"name": "Rice dwarf virus", "abbrev": "RDV", "match": "Phytoreovirus alphaoryzae", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_yellowing", "rice_leafhopper_vector"]},
    {"name": "Rice stripe virus", "abbrev": "RSV", "match": "Tenuivirus oryzaclavatae", "traits": ["rice_stunting", "rice_leaf_striping_mosaic", "rice_leafhopper_vector"]},
    {"name": "Rice black-streaked dwarf virus", "abbrev": "RBSDV", "match": "Fijivirus alporyzae", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_distortion_galls", "rice_planthopper_vector"]},
    {"name": "Rice grassy stunt virus", "abbrev": "RGSV", "match": "Tenuivirus oryzabrevis", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_yellowing", "rice_planthopper_vector"]},
    {"name": "Rice tungro bacilliform virus", "abbrev": "RTBV", "match": "Tungrovirus oryzae", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_leaf_yellowing", "rice_leafhopper_vector"]},
    {"name": "Rice tungro spherical virus", "abbrev": "RTSV", "match": "Waikavirus oryzae", "traits": ["rice_leafhopper_vector"]},
    {"name": "Rice yellow stunt virus", "abbrev": "RYSV", "match": "Alphanucleorhabdovirus oryzae", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_leaf_yellowing", "rice_leafhopper_vector"]},
    {"name": "Rice bunchy stunt virus", "abbrev": "RBSV", "match": None, "literature_id": "rice-virus:rbsv", "taxonomic_status": "literature-listed; not mapped to a current ICTV species", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leafhopper_vector"]},
    {"name": "Rice ragged stunt virus", "abbrev": "RRSV", "match": "Oryzavirus oryzae", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_distortion_galls", "rice_planthopper_vector"]},
    {"name": "Rice gall dwarf virus", "abbrev": "RGDV", "match": "Phytoreovirus betaoryzae", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_leaf_distortion_galls", "rice_leafhopper_vector"]},
    {"name": "Rice necrosis mosaic virus", "abbrev": "RNMV", "match": "Bymovirus oryzae", "traits": ["rice_stunting", "rice_leaf_yellowing", "rice_soilborne_vector"]},
    {"name": "Southern rice black-streaked dwarf virus", "abbrev": "SRBSDV", "match": "Fijivirus boryzae", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_distortion_galls", "rice_planthopper_vector"]},
    {"name": "Rice stripe mosaic virus", "abbrev": "RSMV", "match": "Betacytorhabdovirus oryzae", "traits": ["rice_stunting", "rice_increased_tillering", "rice_leaf_striping_mosaic", "rice_leafhopper_vector"]},
    {"name": "Rice virus A", "abbrev": "RVA", "match": None, "literature_id": "rice-virus:rva", "taxonomic_status": "tentative literature name; unassigned to an ICTV species", "source_pmid": "28695276", "source_title": "Complete genome sequence of rice virus A, a new member of the family Tombusviridae", "source_doi": "10.1007/s00705-017-3472-4", "traits": ["rice_symptomatic_plant_association"]},
    {"name": "Rice curl dwarf-associated virus", "abbrev": "RCDaV", "match": None, "literature_id": "rice-virus:rcdav", "taxonomic_status": "tentative literature name; unassigned to an ICTV species", "source_pmid": "34721366", "source_title": "A Novel Rice Curl Dwarf-Associated Picornavirus Encodes a 3C Serine Protease Recognizing Uncommon EPT/S Cleavage Sites", "source_doi": "10.3389/fmicb.2021.757451", "traits": ["rice_stunting", "rice_leaf_distortion_galls"]},
    {"name": "Rice dwarf-associated bunya-like virus", "abbrev": "RDaBV", "match": "Orthodiscovirus oryzae", "traits": ["rice_stunting", "rice_reduced_tillering"]},
    {"name": "Rice yellow mottle virus", "abbrev": "RYMV", "match": "Sobemovirus RYMV", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_leaf_yellowing", "rice_leaf_striping_mosaic", "rice_beetle_vector"]},
    {"name": "Rice stripe necrosis virus", "abbrev": "RSNV", "match": "Benyvirus oryzae", "traits": ["rice_stunting", "rice_leaf_yellowing", "rice_leaf_striping_mosaic", "rice_soilborne_vector"]},
    {"name": "Rice hoja blanca virus", "abbrev": "RHBV", "match": "Tenuivirus oryzalbae", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_leaf_yellowing", "rice_leaf_striping_mosaic", "rice_planthopper_vector"]},
    {"name": "Rice dwarf polerovirus", "abbrev": "RDPV", "match": "Polerovirus RDPV", "source_pmid": "41649689", "source_title": "Discovery of a novel whitefly- and aphid-transmitted polerovirus on rice plants with dwarfing and fewer tillering symptoms", "source_doi": "10.1007/s44297-024-00033-0", "traits": ["rice_stunting", "rice_reduced_tillering", "rice_whitefly_aphid_vector"]},
]


def _coords(records: list[dict[str, Any]]) -> None:
    """Fit a deterministic PCA map from the presence-only rice trait matrix."""
    columns = list(FEATURES)
    matrix = np.array([[1.0 if feature in item["traits"] else 0.0 for feature in columns]
                       for item in records], dtype=float)
    centered = matrix - matrix.mean(axis=0)
    if centered.size and np.linalg.matrix_rank(centered) >= 1:
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        coords = u[:, :2] * singular[:2]
        for dim in range(coords.shape[1]):
            pivot = int(np.argmax(np.abs(coords[:, dim])))
            if coords[pivot, dim] < 0:
                coords[:, dim] *= -1
    else:
        coords = np.zeros((len(records), 2))
    if coords.shape[1] < 2:
        coords = np.pad(coords, ((0, 0), (0, 2 - coords.shape[1])))
    for item, (x, y) in zip(records, coords):
        item["x"], item["y"] = float(x), float(y)
        item["coordinate_source_host_id"] = RICE_HOST_ID


def export_rice_virus_screen(refresh_dir: Path, destination: Path) -> dict[str, Any]:
    """Write a compact, versioned source snapshot of the rice refresh results."""
    selection = pd.read_parquet(refresh_dir / "pathogen_host_selection.parquet")
    searches = pd.read_parquet(refresh_dir / "pair_cache/tables/pair_searches.parquet")
    articles = pd.read_parquet(refresh_dir / "pair_cache/tables/articles.parquet")
    ledger = pd.read_parquet(refresh_dir / "screening/pathogen_screen_ledger_common_names.parquet")
    articles["pmid"] = articles.pmid.astype(str)
    article_by_pmid = articles.drop_duplicates("pmid", keep="last").set_index("pmid").to_dict("index")
    positives = ledger.loc[ledger.status.astype(str).eq("completed") & ledger.answer.astype(str).str.casefold().eq("yes")]
    positives = positives.drop_duplicates(["pathogen_id", "host_id", "pmid"], keep="last")
    evidence = []
    for decision in positives.to_dict("records"):
        article = article_by_pmid.get(str(decision["pmid"]), {})
        evidence.append({"pathogen_id": str(decision["pathogen_id"]), "pmid": str(decision["pmid"]),
            "title": str(article.get("title") or ""), "publication_year": article.get("publication_year"),
            "doi": str(article.get("doi") or ""), "p_yes": float(decision.get("p_yes") or 0),
            "abstract_hash": str(decision.get("abstract_hash") or "")})
    selection_by_id = {str(row["pathogen_id"]): row for row in selection.to_dict("records")}
    pairs = {}
    for search in searches.to_dict("records"):
        pid = str(search["pathogen_id"])
        row = selection_by_id.get(pid, {})
        pairs[pid] = {"pair_selection_status": str(row.get("selection_status") or "no_sampled_records"),
            "n_positive_records": int(row.get("n_positive_records") or 0),
            "n_sampled_records": int(row.get("n_sampled_records") or 0),
            "n_abstracts_fetched": int(row.get("n_sampled_records_with_abstract") or 0),
            "n_jev_assessments": int(row.get("n_completed_screenings") or 0),
            "n_pair_search_records": int(search.get("search_count") or 0)}
    payload = {"source_version": "rice-virus-refresh-v1", "query_version": "pathogen-host-ictv-vmr-virus-names-title-abstract-v3",
        "screening_model": "typesafe/jev-1.14", "screening_threshold": 0.7, "sample_maximum": 20,
        "sample_oldest_count": 10, "host_id": RICE_HOST_ID, "pair_profiles": pairs, "positive_evidence": evidence}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return payload


def build_rice_virus_profiles(pathogens: pd.DataFrame, screen_source_path: Path | None = None) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Resolve 20 literature-recognized rice viruses into profiles and screened evidence."""
    by_name = {str(row.accepted_name).casefold(): row._asdict() for row in pathogens.itertuples(index=False)}
    profiles: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    for virus in LITERATURE_VIRUSES:
        row = by_name.get(str(virus.get("match") or "").casefold()) if virus.get("match") else None
        pid = str(row["pathogen_id"]) if row else str(virus.get("literature_id"))
        group = str(row.get("pathogen_group") or "viruses") if row else "viruses"
        source_pmid = str(virus.get("source_pmid") or LITERATURE_REVIEW_PMID)
        source_title = str(virus.get("source_title") or "A Review of Vector-Borne Rice Viruses")
        source_doi = str(virus.get("source_doi") or "10.3390/v14102258")
        profiles.append({"pathogen_id": pid, "organism_group": group, "accepted_name": str(row.get("accepted_name") or virus["name"]) if row else virus["name"],
            "literature_name": virus["name"], "abbreviation": virus["abbrev"], "taxonomic_status": virus.get("taxonomic_status", "current ICTV species"),
            "traits": list(virus["traits"]), "source_pmid": source_pmid,
            "feature_values": {feature: (1.0 if feature in virus["traits"] else None) for feature in FEATURES},
            "feature_counts": {}, "feature_profile_version": "rice-virus-literature-v1",
            "evidence_file": "evidence.json.gz", "relationship_context": "disease"})
        evidence.append({"organism_id": pid, "organism_group": "viruses",
            "accepted_name": str(row.get("accepted_name") or virus["name"]) if row else virus["name"],
            "pathogen_id": pid, "host_id": RICE_HOST_ID, "relationship_context": "disease",
            "pmid": source_pmid, "title": source_title, "publication_date": None,
            "publication_year": None, "doi": source_doi, "decision": "literature_feature_source",
            "evidence_kind": "literature_feature", "accepted_pathogen_evidence": True,
            "feature_ids_yes": list(virus["traits"])})
    _coords(profiles)
    pair_profiles: dict[str, dict[str, Any]] = {}
    if screen_source_path and screen_source_path.exists():
        source = json.loads(screen_source_path.read_text(encoding="utf-8"))
        pair_profiles = source.get("pair_profiles", {})
        for record in source.get("positive_evidence", []):
            pid = str(record["pathogen_id"])
            evidence.append({"organism_id": pid, "organism_group": "viruses",
                "accepted_name": str(by_name.get(pid, {}).get("accepted_name") or ""),
                "pathogen_id": pid, "host_id": RICE_HOST_ID, "relationship_context": "disease",
                "pmid": str(record["pmid"]), "title": str(record.get("title") or ""),
                "publication_date": None, "publication_year": record.get("publication_year"),
                "doi": str(record.get("doi") or ""), "p_yes": float(record.get("p_yes") or 0),
                "decision": "supported", "accepted_pathogen_evidence": True,
                "abstract_hash": str(record.get("abstract_hash") or "")})
    return profiles, evidence, pair_profiles


def add_rice_virus_profiles(data: dict[str, Any], taxa: list[dict[str, Any]], pathogens: pd.DataFrame,
                            screen_source_path: Path | None = None) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    """Add rice-only traits and coordinates; preserve all other host profiles."""
    profiles, new_evidence, screen_pairs = build_rice_virus_profiles(pathogens, screen_source_path)
    by_id = {str(row["organism_id"]): row for row in taxa}
    previously_available: set[str] = set()
    for taxon in taxa:
        pair = (taxon.get("pair_profiles") or {}).get(RICE_HOST_ID, {})
        values = pair.get("feature_values") or {}
        has_features = any(value is not None for value in values.values())
        xy = (pair.get("x"), pair.get("y"))
        has_coordinates = all(value is not None and np.isfinite(float(value)) for value in xy)
        if has_features or has_coordinates:
            previously_available.add(str(taxon["organism_id"]))
    new_available_count = 0
    custom = []
    for profile in profiles:
        pid = profile["pathogen_id"]
        taxon = by_id.get(pid)
        if taxon is None:
            custom.append({"organism_id": pid, "pathogen_id": pid, "organism_group": "viruses",
                "accepted_name": profile["accepted_name"], "source_taxon_id": "", "taxonomy_source": "rice-virus-literature",
                "taxonomy_release": "rice-virus-literature-v1", "rank": "literature-associated virus",
                "parent_taxon_id": None, "taxonomy_path": "group=viruses", "mapping_status": "literature-associated",
                "synonyms": [profile["abbreviation"]], "pair_profiles": {}})
            taxon = custom[-1]
            by_id[pid] = taxon
        taxon["synonyms"] = sorted(set([*(taxon.get("synonyms") or []), profile["literature_name"], profile["abbreviation"]]),
                                   key=str.casefold)
        pair_map = taxon.get("pair_profiles")
        if not isinstance(pair_map, dict):
            pair_map = {}
            taxon["pair_profiles"] = pair_map
        pair = pair_map.setdefault(RICE_HOST_ID, {"pair_selection_status": "literature profile", "relationship_context": "disease",
            "n_pathogen_pmids": 0, "n_jev_accepted_pmids": 0, "n_jev_assessments": 0, "n_positive_records": 0,
            "n_sampled_records": 0, "n_abstracts_fetched": 0, "n_pair_search_records": 0, "evidence_file": "evidence.json.gz"})
        # Screening is additive: only use refreshed pair decisions when present.
        refreshed = screen_pairs.get(pid)
        if refreshed:
            pair.update({"pair_selection_status": refreshed["pair_selection_status"],
                "n_pathogen_pmids": refreshed["n_positive_records"],
                "n_jev_accepted_pmids": refreshed["n_positive_records"],
                "n_positive_records": refreshed["n_positive_records"],
                "n_sampled_records": refreshed["n_sampled_records"],
                "n_abstracts_fetched": refreshed["n_abstracts_fetched"],
                "n_jev_assessments": refreshed["n_jev_assessments"],
                "n_pair_search_records": refreshed["n_pair_search_records"]})
        pair.update({"feature_values": profile["feature_values"], "feature_counts": profile["feature_counts"],
            "feature_profile_version": profile["feature_profile_version"], "x": profile["x"], "y": profile["y"],
            "feature_source_pmid": profile["source_pmid"],
            "coordinate_source_host_id": RICE_HOST_ID, "evidence_file": "evidence.json.gz"})
        if pid not in previously_available:
            new_available_count += 1
    taxa.extend(custom)
    data["feature_profile_counts"].setdefault(RICE_HOST_ID, {}).setdefault("viruses", 0)
    data["profile_availability_counts"].setdefault(RICE_HOST_ID, {}).setdefault("viruses", 0)
    data["feature_profile_counts"][RICE_HOST_ID]["viruses"] = sum(any(v is not None for v in row["feature_values"].values()) for row in profiles)
    data["profile_availability_counts"][RICE_HOST_ID]["viruses"] = (
        int(data["profile_availability_counts"][RICE_HOST_ID]["viruses"] or 0) + new_available_count)
    data["coordinate_host_ids"] = sorted(set(data.get("coordinate_host_ids", [])) | {RICE_HOST_ID})
    data.setdefault("coordinate_contexts", {})[RICE_HOST_ID] = "20 literature-listed rice-virus positions are fitted from curated rice traits; other saved rice-view positions remain from the human-disease embedding. Blank traits mean not reported, not absence."
    data["coordinates_ready"] = True
    data.setdefault("stats", {})["rice_literature_virus_profiles"] = len(profiles)
    data["stats"]["rice_refreshed_pair_profiles"] = len(screen_pairs)
    feature_data = {"features": FEATURES, "feature_profiles": {"viruses": list(FEATURES)}}
    return data, taxa, feature_data, new_evidence

