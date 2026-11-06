"""Host feature-set similarity matrices and shared-feature pathogen maps."""

from __future__ import annotations

import gzip
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def _load_registry(project_root: Path) -> dict[str, Any]:
    """Load canonical host, pathogen-group, and feature question registries."""
    from pathogen_atlas_v4.host_taxonomy.ncbi import PILOT_HOSTS

    root = Path(project_root)
    with (root / "configs/host_choices.yaml").open(encoding="utf-8") as stream:
        host_config = yaml.safe_load(stream)
    with (root / "configs/pathogen_taxonomies.yaml").open(encoding="utf-8") as stream:
        pathogen_config = yaml.safe_load(stream)
    with (root / "configs/features.yaml").open(encoding="utf-8") as stream:
        feature_catalogue = yaml.safe_load(stream)

    aliases_by_id: dict[str, list[str]] = {}
    for name, taxon_id in host_config["host_name_to_id"].items():
        aliases_by_id.setdefault(str(taxon_id), []).append(str(name))
    hosts = [
        {
            "host_id": host_id,
            "taxon_id": host_id.removeprefix("ncbi:"),
            "label": aliases[0],
            "aliases": aliases[1:],
            "host_group": PILOT_HOSTS.get(host_id.removeprefix("ncbi:"), "other"),
        }
        for host_id, aliases in aliases_by_id.items()
    ]
    # Keep the comparison matrix in biological host-group order: human,
    # livestock, fish, then crops. This places the fish profiles after Pig and
    # before the plant hosts while preserving their configured order within a
    # group.
    host_group_order = {"human": 0, "livestock": 1, "fish": 2, "crop": 3}
    hosts.sort(key=lambda host: host_group_order.get(host["host_group"], 4))
    pathogen_groups = [
        {"pathogen_group": str(group_id), "label": str(group.get("label", group_id))}
        for group_id, group in pathogen_config["pathogen_groups"].items()
    ]
    return {
        "root": root,
        "hosts": hosts,
        "pathogen_groups": pathogen_groups,
        "feature_catalogue": feature_catalogue,
        "relationship_context": "disease",
        "profile_version": str(feature_catalogue.get("profile_version", "")),
    }


def _resolved_feature_signatures(
    registry: dict[str, Any], host: dict[str, str], pathogen_group: str,
    scope_filter: set[str] | None = None,
) -> set[tuple[str, str, str]]:
    from pathogen_atlas_v4.features import resolve_profile

    resolved = resolve_profile(
        registry["feature_catalogue"],
        pathogen_group=pathogen_group,
        host_group=host["host_group"],
        host_taxon_id=host["taxon_id"],
        relationship_context=registry["relationship_context"],
        profile_version=registry["profile_version"],
    )
    records = resolved["ordered_features"]
    if scope_filter is not None:
        records = [record for record in records if record.get("scope") in scope_filter]
    return {
        (
            str(record["feature_id"]),
            str(record.get("question_version", "")),
            str(record["scope"]),
        )
        for record in records
    }


def _positive_feature_support_counts(
    project_root: Path, groups: list[dict[str, str]],
) -> dict[tuple[str, str, str], int]:
    """Count distinct pathogen records with positive values per group/host/feature."""
    site_dist = Path(project_root) / "sites/pathogen-atlas-explorer/dist"
    metadata_path = site_dist / "explorer-data.json.gz"
    if not metadata_path.is_file():
        raise FileNotFoundError(
            f"Explorer release data not found at {metadata_path}; build the Explorer data before comparing profiles."
        )
    with gzip.open(metadata_path, "rt", encoding="utf-8") as stream:
        metadata = json.load(stream)

    supporting_pathogens: dict[tuple[str, str, str], set[str]] = {}
    for group in groups:
        group_id = group["pathogen_group"]
        for taxon in _load_group_rows(site_dist, metadata, group_id):
            pathogen_id = str(taxon["organism_id"])
            for host_id, profile in (taxon.get("pair_profiles") or {}).items():
                for feature_id, value in (profile.get("feature_values") or {}).items():
                    try:
                        score = float(value)
                        has_value = math.isfinite(score) and score > 0
                    except (TypeError, ValueError):
                        has_value = False
                    if has_value:
                        key = (group_id, str(host_id), str(feature_id))
                        supporting_pathogens.setdefault(key, set()).add(pathogen_id)
    return {key: len(pathogens) for key, pathogens in supporting_pathogens.items()}


def build_jaccard_analysis(
    project_root: Path, *, scope_filter: Iterable[str] | None = None,
    min_pathogens_per_host: int = 2,
) -> dict[str, Any]:
    """Build host-by-host Jaccard matrices after a per-host pathogen-support cutoff."""
    registry = _load_registry(Path(project_root))
    if (
        isinstance(min_pathogens_per_host, bool)
        or not isinstance(min_pathogens_per_host, int)
        or min_pathogens_per_host < 1
    ):
        raise ValueError("min_pathogens_per_host must be a positive integer")
    allowed_scopes = set(scope_filter) if scope_filter is not None else None
    hosts = registry["hosts"]
    groups = registry["pathogen_groups"]
    feature_support_counts = _positive_feature_support_counts(registry["root"], groups)
    feature_sets: dict[tuple[str, str], set[tuple[str, str, str]]] = {}
    profile_rows = []

    for group in groups:
        group_id = group["pathogen_group"]
        for host in hosts:
            configured_identities = _resolved_feature_signatures(registry, host, group_id, allowed_scopes)
            positive_identities = {
                identity for identity in configured_identities
                if feature_support_counts.get((group_id, host["host_id"], identity[0]), 0) >= 1
            }
            identities = {
                identity for identity in positive_identities
                if feature_support_counts.get((group_id, host["host_id"], identity[0]), 0)
                >= min_pathogens_per_host
            }
            feature_sets[(group_id, host["host_id"])] = identities
            profile_rows.append({
                "pathogen_group": group_id,
                "host": host["label"],
                "host_id": host["host_id"],
                "host_group": host["host_group"],
                "configured_feature_count": len(configured_identities),
                "positive_feature_count": len(positive_identities),
                "feature_count": len(identities),
            })

    host_labels = [host["label"] for host in hosts]
    matrices: dict[str, pd.DataFrame] = {}
    shared_feature_matrices: dict[str, pd.DataFrame] = {}
    pairwise_rows = []
    for group in groups:
        group_id = group["pathogen_group"]
        matrix = pd.DataFrame(np.nan, index=host_labels, columns=host_labels, dtype=float)
        shared_matrix = pd.DataFrame(0, index=host_labels, columns=host_labels, dtype=int)
        for left, right in combinations(hosts, 2):
            left_set = feature_sets[(group_id, left["host_id"])]
            right_set = feature_sets[(group_id, right["host_id"])]
            intersection = left_set & right_set
            union = left_set | right_set
            score = len(intersection) / len(union) if union else np.nan
            matrix.loc[left["label"], right["label"]] = score
            matrix.loc[right["label"], left["label"]] = score
            shared_matrix.loc[left["label"], right["label"]] = len(intersection)
            shared_matrix.loc[right["label"], left["label"]] = len(intersection)
            pairwise_rows.append({
                "pathogen_group": group_id,
                "host_a": left["label"],
                "host_b": right["label"],
                "features_a": len(left_set),
                "features_b": len(right_set),
                "intersection": len(intersection),
                "union": len(union),
                "jaccard": score,
                "shared_feature_ids": "; ".join(sorted(item[0] for item in intersection)),
            })
        for host in hosts:
            own_set = feature_sets[(group_id, host["host_id"])]
            matrix.loc[host["label"], host["label"]] = 1.0 if own_set else np.nan
            shared_matrix.loc[host["label"], host["label"]] = len(own_set)
        matrices[group_id] = matrix
        shared_feature_matrices[group_id] = shared_matrix

    profile_counts = pd.DataFrame(profile_rows)
    pairwise_comparisons = pd.DataFrame(pairwise_rows)
    return {
        "registry": registry,
        "profile_counts": profile_counts,
        "matrices": matrices,
        "shared_feature_matrices": shared_feature_matrices,
        "pairwise_comparisons": pairwise_comparisons,
        "scope_filter": allowed_scopes,
        "min_pathogens_per_host": min_pathogens_per_host,
    }


def plot_jaccard_matrices(analysis: dict[str, Any]) -> None:
    """Render relative Jaccard similarity beside absolute shared-feature counts."""
    hosts = analysis["registry"]["hosts"]
    host_labels = [host["label"] for host in hosts]
    for group in analysis["registry"]["pathogen_groups"]:
        group_id = group["pathogen_group"]
        matrix = analysis["matrices"][group_id]
        shared_matrix = analysis["shared_feature_matrices"][group_id]
        max_shared = int(shared_matrix.to_numpy().max()) if shared_matrix.size else 0
        fig, axes = plt.subplots(1, 2, figsize=(18, 9), constrained_layout=True)
        panels = [
            (axes[0], matrix, 1, "Relative feature overlap (Jaccard)", "Jaccard similarity", ".2f"),
            (axes[1], shared_matrix, max(max_shared, 1), "Absolute shared feature count", "Shared features", "d"),
        ]
        for ax, values, vmax, title, colorbar_label, number_format in panels:
            image = ax.imshow(values.to_numpy(dtype=float), vmin=0, vmax=vmax, cmap="YlGnBu")
            ax.set_xticks(range(len(host_labels)), host_labels, rotation=45, ha="right")
            ax.set_yticks(range(len(host_labels)), host_labels)
            ax.set_title(title, pad=12)
            ax.set_xlabel("Host")
            ax.set_ylabel("Host")
            ax.set_aspect("equal")
            for row in range(len(host_labels)):
                for col in range(len(host_labels)):
                    value = values.iat[row, col]
                    if pd.notna(value):
                        ax.text(
                            col, row, format(value, number_format), ha="center", va="center",
                            color="white" if value >= 0.58 * vmax else "#183247", fontsize=8,
                        )
            fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label=colorbar_label)
        fig.suptitle(f"{group['label']} · host feature-profile comparison", fontsize=14)
        plt.show()


def _load_group_rows(site_dist: Path, metadata: dict[str, Any], group_id: str) -> list[dict[str, Any]]:
    with gzip.open(site_dist / "taxa" / f"{group_id}.json.gz", "rt", encoding="utf-8") as stream:
        packed_rows = json.load(stream)
    columns = metadata["species_columns"]
    return [dict(zip(columns, values)) for values in packed_rows]


def _is_supported(pair: dict[str, Any]) -> bool:
    return str(pair.get("pair_selection_status", "")).casefold() in {"accepted", "supported"}


def _support_count(pair: dict[str, Any]) -> int:
    value = pair.get("n_jev_accepted_pmids")
    if value is None:
        value = pair.get("n_pathogen_pmids", 0)
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _fit_shared_host_map(
    registry: dict[str, Any], site_dist: Path, metadata: dict[str, Any],
    host_a_id: str, host_b_id: str, group_id: str,
    scope_filter: set[str] | None,
) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    from imste import IteratedMinimumSpanningTreeEmbedder
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if host_a_id == host_b_id:
        return None, {"message": "Choose two different hosts."}
    host_by_id = {host["host_id"]: host for host in registry["hosts"]}
    host_a, host_b = host_by_id[host_a_id], host_by_id[host_b_id]
    label_a, label_b = host_a["label"], host_b["label"]
    signatures_a = _resolved_feature_signatures(registry, host_a, group_id, scope_filter)
    signatures_b = _resolved_feature_signatures(registry, host_b, group_id, scope_filter)
    shared_signatures = signatures_a & signatures_b
    shared_feature_ids = sorted({signature[0] for signature in shared_signatures})
    if not shared_feature_ids:
        return None, {"message": "These hosts have no applicable shared features for this pathogen group."}

    pathogen_vectors = []
    excluded_no_shared_values = 0
    group_rows = _load_group_rows(site_dist, metadata, group_id)
    confirmed_for_either_host = 0
    for taxon in group_rows:
        profiles = taxon.get("pair_profiles") or {}
        pair_a, pair_b = profiles.get(host_a_id) or {}, profiles.get(host_b_id) or {}
        affects_a, affects_b = _is_supported(pair_a), _is_supported(pair_b)
        if not (affects_a or affects_b):
            continue
        confirmed_for_either_host += 1
        values_a = (pair_a.get("feature_values") or {}) if affects_a else {}
        values_b = (pair_b.get("feature_values") or {}) if affects_b else {}
        feature_values = {}
        for feature_id in shared_feature_ids:
            observed = []
            for source in (values_a, values_b):
                try:
                    value = float(source.get(feature_id))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    observed.append(value)
            if observed:
                feature_values[feature_id] = float(sum(observed) / len(observed))
        if not feature_values:
            excluded_no_shared_values += 1
            continue

        if affects_a and affects_b:
            class_name = "Both"
        elif affects_a:
            class_name = "First"
        else:
            class_name = "Second"
        pathogen_vectors.append({
            "pathogen_id": str(taxon["organism_id"]),
            "pathogen": str(taxon["accepted_name"]),
            "class": class_name,
            "feature_values": feature_values,
            "host_1_status": pair_a.get("pair_selection_status") or "not screened",
            "host_2_status": pair_b.get("pair_selection_status") or "not screened",
            "host_1_supporting_pmids": _support_count(pair_a),
            "host_2_supporting_pmids": _support_count(pair_b),
        })

    if len(pathogen_vectors) < 2:
        return None, {
            "message": f"Only {len(pathogen_vectors)} supported pathogen(s) have scored shared-feature values; a 2D comparison needs at least 2.",
            "shared_features": len(shared_feature_ids),
            "excluded_no_shared_values": excluded_no_shared_values,
        }

    # Require an observed value for every plotted pathogen; missing values stay unknown, never zero.
    complete_feature_ids = [
        feature_id for feature_id in shared_feature_ids
        if all(feature_id in row["feature_values"] for row in pathogen_vectors)
    ]
    if not complete_feature_ids:
        return None, {
            "message": "No shared feature has observed values for every supported pathogen; no values were imputed.",
            "shared_features": len(shared_feature_ids),
            "pathogens_with_any_shared_values": len(pathogen_vectors),
            "excluded_no_shared_values": excluded_no_shared_values,
        }

    matrix = np.asarray([
        [row["feature_values"][feature_id] for feature_id in complete_feature_ids]
        for row in pathogen_vectors
    ], dtype=float)
    variable_mask = np.nanstd(matrix, axis=0) > 0
    used_feature_ids = complete_feature_ids
    if not variable_mask.any():
        return None, {
            "message": "The complete shared feature values are constant; IMSTE has no variable dimensions to fit.",
            "shared_features": len(shared_feature_ids),
            "excluded_no_shared_values": excluded_no_shared_values,
        }

    matrix = StandardScaler().fit_transform(matrix)
    n_species = len(pathogen_vectors)
    if n_species == 2:
        # With two samples, the neighbor graph has no non-neighbor pair for IMSTE
        # repulsion. PCA still displays both observations using all shared features.
        coordinates = PCA(n_components=min(2, matrix.shape[1])).fit_transform(matrix)
        if coordinates.shape[1] == 1:
            coordinates = np.column_stack([coordinates[:, 0], np.zeros(n_species)])
        embedding_method = "PCA (two-pathogen fallback)"
        n_msts = 0
    else:
        n_msts = min(10, max(2, n_species - 3))
        mapper = IteratedMinimumSpanningTreeEmbedder(
            n_msts=n_msts, n_components=2, n_epochs=300, batch_size=32768,
            learning_rate=0.05, negative_ratio=5, lambda_rep=0.5,
            random_state=42, epsilon=0.0001, device="cpu", logistic_margin=1.0,
            logistic_temperature=0.5, mst_method="famst",
            mst_neighbors=min(15, n_species - 1), mst_inter_component_edges=5,
            mst_max_neighbors=None, mst_representatives_per_component=10,
        )
        try:
            coordinates = np.asarray(mapper.fit_transform(matrix), dtype=float)
        except (ValueError, RuntimeError, FloatingPointError) as error:
            return None, {"message": f"IMSTE could not fit this selection: {error}"}
        embedding_method = "IMSTE"
    if coordinates.shape != (n_species, 2) or not np.isfinite(coordinates).all():
        return None, {"message": "IMSTE returned invalid coordinates for this selection."}
    if float(np.max(coordinates.std(axis=0))) < 1e-6:
        return None, {"message": "IMSTE returned a collapsed map for this selection."}

    for row, (x, y) in zip(pathogen_vectors, coordinates):
        row["IMSTE 1"], row["IMSTE 2"] = float(x), float(y)
    report = {
        "message": "",
        "group": group_id,
        "host_a": label_a,
        "host_b": label_b,
        "shared_features": len(shared_feature_ids),
        "used_features": used_feature_ids,
        "embedding_method": embedding_method,
        "n_pathogens": len(pathogen_vectors),
        "n_group_pathogens": len(group_rows),
        "n_confirmed_for_either_host": confirmed_for_either_host,
        "excluded_no_shared_values": excluded_no_shared_values,
        "n_msts": n_msts,
    }
    return pd.DataFrame(pathogen_vectors), report


def create_interactive_host_comparison(
    project_root: Path, *, scope_filter: Iterable[str] | None = None,
) -> Any:
    """Create host/group dropdowns and a Plotly shared-feature IMSTE view."""
    import ipywidgets as widgets
    import plotly.express as px
    from IPython.display import clear_output, display

    root = Path(project_root)
    registry = _load_registry(root)
    site_dist = root / "sites/pathogen-atlas-explorer/dist"
    with gzip.open(site_dist / "explorer-data.json.gz", "rt", encoding="utf-8") as stream:
        metadata = json.load(stream)
    allowed_scopes = set(scope_filter) if scope_filter is not None else None

    host_options = [(host["label"], host["host_id"]) for host in registry["hosts"]]
    default_a = next((host["host_id"] for host in registry["hosts"] if host["host_id"] == "ncbi:9606"), host_options[0][1])
    default_b = next((host["host_id"] for host in registry["hosts"] if host["host_id"] == "ncbi:4530"), None)
    if default_b is None or default_b == default_a:
        default_b = next(host_id for _, host_id in host_options if host_id != default_a)
    group_options = [(group["label"], group["pathogen_group"]) for group in registry["pathogen_groups"]]
    default_group = next((group_id for _, group_id in group_options if group_id == "viruses"), group_options[0][1])

    host_a = widgets.Dropdown(options=host_options, value=default_a, description="Host 1:", layout=widgets.Layout(width="28%"))
    host_b = widgets.Dropdown(options=host_options, value=default_b, description="Host 2:", layout=widgets.Layout(width="28%"))
    pathogen_group = widgets.Dropdown(options=group_options, value=default_group, description="Group:", layout=widgets.Layout(width="24%"))
    run_button = widgets.Button(
        description="Fit shared-feature IMSTE", button_style="success", icon="area-chart",
        layout=widgets.Layout(width="220px", flex="0 0 220px"),
    )
    output = widgets.Output()

    def run(_):
        with output:
            clear_output(wait=True)
            points, report = _fit_shared_host_map(
                registry, site_dist, metadata, host_a.value, host_b.value,
                pathogen_group.value, allowed_scopes,
            )
            if points is None:
                print(report["message"])
                for key in ("shared_features", "pathogens_with_any_shared_values", "excluded_no_shared_values"):
                    if key in report:
                        print(f"{key.replace('_', ' ')}: {report[key]}")
                return

            label_a, label_b = report["host_a"], report["host_b"]
            symbols = {"Both": "circle", "First": "square", "Second": "diamond"}
            colors = {"Both": "#16877d", "First": "#315f8c", "Second": "#d28732"}
            figure = px.scatter(
                points, x="IMSTE 1", y="IMSTE 2", color="class", symbol="class",
                symbol_map=symbols, color_discrete_map=colors,
                hover_name="pathogen",
                hover_data={
                    "pathogen_id": True,
                    "host_1_status": True,
                    "host_2_status": True,
                    "host_1_supporting_pmids": True,
                    "host_2_supporting_pmids": True,
                    "class": True,
                    "IMSTE 1": ":.3f",
                    "IMSTE 2": ":.3f",
                },
                title=f"{report['group'].title()} · shared-feature {report['embedding_method']} · {label_a} vs {label_b}",
                labels={"class": "Confirmed disease association", "IMSTE 1": "IMSTE 1", "IMSTE 2": "IMSTE 2"},
            )
            figure.update_traces(marker={"size": 11, "line": {"width": 0.7, "color": "white"}})
            figure.update_layout(
                width=950, height=760, legend_title_text="Confirmed association",
                xaxis_title="IMSTE 1", yaxis_title="IMSTE 2",
                margin={"l": 55, "r": 30, "t": 80, "b": 55},
            )
            display(figure)
            print(f"Applicable shared features: {report['shared_features']}")
            print(f"Complete, variable features used: {len(report['used_features'])} — {', '.join(report['used_features'])}")
            print(f"Embedding: {report['embedding_method']} ({report['n_pathogens']} pathogen rows × {len(report['used_features'])} feature columns)")
            print(
                f"Pathogens plotted: {report['n_pathogens']} of {report['n_confirmed_for_either_host']} "
                f"with confirmed evidence for either host ({report['n_group_pathogens']} in the group); "
                f"{report['excluded_no_shared_values']} confirmed pathogen(s) lacked shared-feature values."
            )
            if report["n_msts"]:
                print(f"n_msts: {report['n_msts']}")

    run_button.on_click(run)
    controls = widgets.VBox([
        widgets.HBox(
            [host_a, host_b, pathogen_group, run_button],
            layout=widgets.Layout(width="100%", flex_flow="row nowrap", align_items="flex-end", gap="10px"),
        ),
        output,
    ])
    return controls
