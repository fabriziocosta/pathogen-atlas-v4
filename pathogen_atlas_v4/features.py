"""Applicability-based feature profile resolution."""

from __future__ import annotations

from typing import Any

from pathogen_atlas_v4.decision_identity import stable_hash
from pathogen_atlas_v4.schemas import validate_scope

MODULE_ORDER = ("core", "pathogen_group", "host_group", "host_taxon", "pair_context")


def resolve_profile(catalogue: dict[str, Any], *, pathogen_group: str, host_group: str,
                    host_taxon_id: str, relationship_context: str, profile_version: str) -> dict[str, Any]:
    """Resolve ordered modules and list excluded features with their applicability reason."""
    modules = catalogue.get("modules", {})
    selected: list[dict[str, Any]] = []
    excluded: list[dict[str, str]] = []
    keys = [("core", "core"), ("pathogen_group", pathogen_group), ("host_group", host_group),
            ("host_taxon", host_taxon_id), ("pair_context", relationship_context)]
    for module_type, selector in keys:
        module_key = f"{module_type}:{selector}" if module_type != "core" else "core"
        for feature in modules.get(module_key, []):
            scope = validate_scope(feature["scope"])
            selected.append({**feature, "scope": scope, "module": module_key})
    seen: set[str] = set()
    ordered = []
    for feature in selected:
        feature_id = str(feature["feature_id"])
        if feature_id in seen:
            raise ValueError(f"Duplicate feature_id in resolved profile: {feature_id}")
        seen.add(feature_id)
        ordered.append(feature)
    not_applicable = []
    selected_ids = {item["feature_id"] for item in ordered}
    for module_key, module_features in modules.items():
        if module_key == "core" or module_key in {f"pathogen_group:{pathogen_group}", f"host_group:{host_group}",
                                                   f"host_taxon:{host_taxon_id}", f"pair_context:{relationship_context}"}:
            continue
        if ":" not in module_key:
            continue
        selector_type, selector_value = module_key.split(":", 1)
        selector_field = {"pathogen_group": "pathogen_group", "host_group": "host_group",
                          "host_taxon": "host_taxon_id", "pair_context": "relationship_context"}.get(selector_type)
        if not selector_field:
            continue
        for feature in module_features:
            if feature.get("feature_id") not in selected_ids:
                not_applicable.append({"feature_id": str(feature.get("feature_id")),
                                       "reason": f"requires {selector_field}={selector_value}",
                                       "selector": selector_field, "required_value": selector_value})
    profile_id = catalogue.get("profile_id", "host-aware-default")
    resolved = {"profile_id": profile_id, "profile_version": profile_version,
            "selectors": {"pathogen_group": pathogen_group, "host_group": host_group,
                          "host_taxon_id": host_taxon_id, "relationship_context": relationship_context},
            "ordered_features": ordered, "not_applicable": not_applicable}
    resolved["resolved_profile_hash"] = stable_hash(resolved)
    return resolved
