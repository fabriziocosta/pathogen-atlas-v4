"""Trait-aware, versioned feature question profile resolution."""

from __future__ import annotations

from typing import Any, Iterable

from pathogen_atlas_v4.decision_identity import stable_hash
from pathogen_atlas_v4.schemas import validate_scope

MODULE_ORDER = ("core", "pathogen_group", "host_group", "host_taxon", "pair_context")


def _matches(applicability: dict[str, Any], traits: dict[str, set[str]]) -> tuple[bool, str | None]:
    """Evaluate positive trait predicates. Unknown traits never count as false facts."""
    for axis in ("pathogen", "host", "pair"):
        rule = applicability.get(f"{axis}_traits", {}) or {}
        available = traits[axis]
        all_of = set(map(str, rule.get("all_of", [])))
        any_of = set(map(str, rule.get("any_of", [])))
        missing = sorted(all_of - available)
        if missing:
            return False, f"requires {axis} traits: {', '.join(missing)}"
        if any_of and not any_of.intersection(available):
            return False, f"requires any {axis} trait: {', '.join(sorted(any_of))}"
    return True, None


def resolve_profile(catalogue: dict[str, Any], *, pathogen_group: str, host_group: str,
                    host_taxon_id: str, relationship_context: str, profile_version: str,
                    pathogen_traits: Iterable[str] | None = None,
                    host_traits: Iterable[str] | None = None,
                    pair_traits: Iterable[str] | None = None) -> dict[str, Any]:
    """Resolve ordered modules, then apply trait predicates and report exclusions.

    New catalogues store question records once under ``questions`` and refer to them
    by ID in modules. Legacy inline feature records remain readable for migration.
    Selector traits are always included; explicitly supplied traits add curated
    taxon- or pair-level knowledge without changing selector semantics.
    """
    modules = catalogue.get("modules", {})
    questions = catalogue.get("questions", {})
    assignments = catalogue.get("trait_assignments", {})
    traits = {
        "pathogen": {f"pathogen_group:{pathogen_group}",
                     *map(str, assignments.get("pathogen_group", {}).get(pathogen_group, [])),
                     *map(str, pathogen_traits or [])},
        "host": {f"host_group:{host_group}", f"host_taxon:{host_taxon_id}",
                 *map(str, assignments.get("host_group", {}).get(host_group, [])),
                 *map(str, assignments.get("host_taxon", {}).get(str(host_taxon_id), [])),
                 *map(str, host_traits or [])},
        "pair": {f"pair_context:{relationship_context}",
                 *map(str, assignments.get("pair_context", {}).get(relationship_context, [])),
                 *map(str, pair_traits or [])},
    }
    keys = [("core", "core"), ("pathogen_group", pathogen_group), ("host_group", host_group),
            ("host_taxon", host_taxon_id), ("pair_context", relationship_context)]
    selected: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for module_type, selector in keys:
        module_key = f"{module_type}:{selector}" if module_type != "core" else "core"
        for item in modules.get(module_key, []):
            if isinstance(item, str):
                if item not in questions:
                    raise ValueError(f"Unknown feature_id {item!r} referenced by {module_key}")
                feature = {**questions[item], "feature_id": item}
            else:
                feature = dict(item)
            feature["scope"] = validate_scope(feature["scope"])
            feature["module"] = module_key
            applies, reason = _matches(feature.get("applicability", {}), traits)
            if applies:
                selected.append(feature)
            else:
                excluded.append({"feature_id": str(feature.get("feature_id")), "reason": reason or "trait mismatch",
                                 "module": module_key})

    seen: set[str] = set()
    ordered: list[dict[str, Any]] = []
    for feature in selected:
        feature_id = str(feature["feature_id"])
        if feature_id in seen:
            raise ValueError(f"Duplicate feature_id in resolved profile: {feature_id}")
        seen.add(feature_id)
        ordered.append(feature)

    # Modules outside this profile remain visible as selector-level exclusions;
    # trait-level exclusions above distinguish an in-scope module from an
    # inapplicable question inside that module.
    applicable_modules = {"core", f"pathogen_group:{pathogen_group}", f"host_group:{host_group}",
                          f"host_taxon:{host_taxon_id}", f"pair_context:{relationship_context}"}
    for module_key, items in modules.items():
        if module_key in applicable_modules or ":" not in module_key:
            continue
        selector_type, required = module_key.split(":", 1)
        selector_field = {"pathogen_group": "pathogen_group", "host_group": "host_group",
                          "host_taxon": "host_taxon_id", "pair_context": "relationship_context"}.get(selector_type)
        if not selector_field:
            continue
        for item in items:
            feature_id = item if isinstance(item, str) else item.get("feature_id")
            excluded.append({"feature_id": str(feature_id), "reason": f"requires {selector_field}={required}",
                             "selector": selector_field, "required_value": required})

    resolved = {"profile_id": catalogue.get("profile_id", "host-aware-default"),
                "catalogue_version": catalogue.get("catalogue_version", ""),
                "profile_version": profile_version,
                "selectors": {"pathogen_group": pathogen_group, "host_group": host_group,
                              "host_taxon_id": str(host_taxon_id), "relationship_context": relationship_context},
                "trait_snapshot": {axis: sorted(values) for axis, values in traits.items()},
                "ordered_features": ordered, "not_applicable": excluded}
    resolved["resolved_profile_hash"] = stable_hash(resolved)
    return resolved
