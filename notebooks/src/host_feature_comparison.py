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
            numeric_values = values.to_numpy(dtype=float)
            zero_or_missing = ~np.isfinite(numeric_values) | (numeric_values == 0)
            masked_values = np.ma.masked_array(numeric_values, mask=zero_or_missing)
            cmap = plt.get_cmap("YlGnBu").copy()
            cmap.set_bad("white")
            image = ax.imshow(masked_values, vmin=0, vmax=vmax, cmap=cmap)
            ax.set_xticks(range(len(host_labels)), host_labels, rotation=45, ha="right")
            ax.set_yticks(range(len(host_labels)), host_labels)
            ax.set_title(title, pad=12)
            ax.set_xlabel("Host")
            ax.set_ylabel("Host")
            ax.set_aspect("equal")
            for row in range(len(host_labels)):
                for col in range(len(host_labels)):
                    value = values.iat[row, col]
                    if pd.notna(value) and value != 0:
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
        pathogen_points = []
        for class_name, host_label, pair, is_supported in (
            ("First", label_a, pair_a, affects_a),
            ("Second", label_b, pair_b, affects_b),
        ):
            if not is_supported:
                continue
            source_values = pair.get("feature_values") or {}
            feature_values = {}
            for feature_id in shared_feature_ids:
                try:
                    value = float(source_values.get(feature_id))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    feature_values[feature_id] = value
            if not feature_values:
                continue
            pathogen_points.append({
                "pathogen_id": str(taxon["organism_id"]),
                "pathogen": str(taxon["accepted_name"]),
                "class": class_name,
                "host": host_label,
                "pair_status": pair.get("pair_selection_status") or "not screened",
                "supporting_pmids": _support_count(pair),
                "estimated_pathogen_publications_index": pair.get("estimated_pathogen_publications_index"),
                "pathogenicity_fraction": pair.get("pathogenicity_fraction"),
                "feature_values": feature_values,
            })
        if not pathogen_points:
            excluded_no_shared_values += 1
            continue
        pathogen_vectors.extend(pathogen_points)

    if len(pathogen_vectors) < 2:
        return None, {
            "message": f"Only {len(pathogen_vectors)} host-specific pathogen point(s) have scored shared-feature values; a 2D comparison needs at least 2.",
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
        "n_pathogens": len({row["pathogen_id"] for row in pathogen_vectors}),
        "n_points": len(pathogen_vectors),
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
    import asyncio
    from html import escape

    import ipywidgets as widgets
    import plotly.graph_objects as go
    from IPython.display import HTML, display

    root = Path(project_root)
    registry = _load_registry(root)
    site_dist = root / "sites/pathogen-atlas-explorer/dist"
    with gzip.open(site_dist / "explorer-data.json.gz", "rt", encoding="utf-8") as stream:
        metadata = json.load(stream)
    allowed_scopes = set(scope_filter) if scope_filter is not None else None

    host_options = [(host["label"], host["host_id"]) for host in registry["hosts"]]
    default_a = next((host["host_id"] for host in registry["hosts"] if host["host_id"] == "ncbi:3847"), host_options[0][1])
    default_b = next((host["host_id"] for host in registry["hosts"] if host["host_id"] == "ncbi:4530"), None)
    if default_b is None or default_b == default_a:
        default_b = next(host_id for _, host_id in host_options if host_id != default_a)
    group_options = [(group["label"], group["pathogen_group"]) for group in registry["pathogen_groups"]]
    default_group = next((group_id for _, group_id in group_options if group_id == "bacteria"), group_options[0][1])

    host_a = widgets.Dropdown(
        options=host_options, value=default_a,
        layout=widgets.Layout(width="100%", height="38px"),
    )
    host_b = widgets.Dropdown(
        options=host_options, value=default_b,
        layout=widgets.Layout(width="100%", height="38px"),
    )
    pathogen_group = widgets.Dropdown(
        options=group_options, value=default_group,
        layout=widgets.Layout(width="100%", height="38px"),
    )
    def selector_shell(dropdown: Any) -> Any:
        dropdown.layout = widgets.Layout(width="100%", height="100%")
        dropdown.add_class("host-feature-dropdown")
        shell = widgets.Box(
            [dropdown],
            layout=widgets.Layout(width="100%", height="38px", padding="0"),
        )
        shell.add_class("host-feature-select-shell")
        return shell

    def selector_field(label: str, control: Any, *, last: bool = False) -> Any:
        return widgets.VBox(
            [
                widgets.HTML(
                    value=f'<span style="color:#475d70;font-size:11px;font-weight:700">{label}</span>'
                ),
                control,
            ],
            layout=widgets.Layout(
                flex="1 1 0%", min_width="190px", gap="4px",
                margin="0" if last else "0 14px 0 0",
            ),
        )

    selector_fields = [
        selector_field("Host 1", selector_shell(host_a)),
        selector_field("Host 2", selector_shell(host_b)),
        selector_field("Group", selector_shell(pathogen_group), last=True),
    ]
    connect_pairs = widgets.Checkbox(
        value=True,
        description="Connect same pathogen",
        indent=False,
        layout=widgets.Layout(width="auto"),
    )
    run_button = widgets.Button(
        description="Fit 2D map", button_style="", icon="area-chart",
        tooltip="Fit a 2D map from the selected hosts’ shared feature values",
        layout=widgets.Layout(width="180px", height="38px", flex="0 0 180px"),
    )
    run_button.style.button_color = "#116b73"
    run_button.style.text_color = "#ffffff"
    run_button.style.font_weight = "600"
    run_button.add_class("host-feature-fit-button")
    display(HTML(
        '<style>.host-feature-fit-button{border-radius:10px!important;'
        'padding:8px 16px!important;box-shadow:0 1px 2px rgba(20,65,75,.16)}'
        '.host-feature-fit-button:focus-visible{outline:3px solid rgba(17,107,115,.22);outline-offset:2px}'
        '.host-feature-select-shell{box-sizing:border-box!important;width:100%!important;min-width:0!important;'
        'height:38px!important;border:1px solid #aebdc9!important;border-radius:8px!important;'
        'overflow:hidden!important;padding:0!important;background:#fff!important;box-shadow:none!important}'
        '.host-feature-dropdown{box-sizing:border-box!important;width:100%!important;height:100%!important;'
        'border:0!important;padding:0!important;background:transparent!important;box-shadow:none!important}'
        '.host-feature-dropdown select{box-sizing:border-box!important;width:100%!important;height:100%!important;'
        'border:0!important;border-radius:0!important;outline:0!important;'
        'padding:0 10px!important;background:transparent!important;color:#25384a!important}</style>'
    ))
    output = widgets.Output()
    status = widgets.HTML(
        value="",
        layout=widgets.Layout(width="260px", min_width="240px", height="32px"),
    )
    last_result: dict[str, Any] = {}

    def html_output(content: str) -> dict[str, Any]:
        bundle = HTML(content)._repr_mimebundle_()
        if isinstance(bundle, tuple):
            data, output_metadata = bundle
        else:
            data, output_metadata = bundle, {}
        return {
            "output_type": "display_data",
            "data": data,
            "metadata": output_metadata,
        }

    intro = widgets.HTML(
        value=(
            '<div style="box-sizing:border-box;width:100%;padding:14px 18px;'
            'background:#f1f6f9;border-bottom:1px solid #e1e9ef;'
            'border-radius:11px 11px 0 0;font-family:Inter,-apple-system,BlinkMacSystemFont,Segoe UI,Arial,sans-serif">'
            '<div style="display:flex;align-items:center;gap:14px;flex-wrap:wrap">'
            '<strong style="font-size:16px;color:#203b52">Compare host profiles</strong>'
            '<span style="width:1px;height:20px;background:#cbd6df"></span>'
            '<span style="font-size:13px;color:#62798d">Shared-feature map</span>'
            '</div></div>'
        ),
    )
    selection_row = widgets.HBox(
        selector_fields,
        layout=widgets.Layout(
            width="100%", flex_flow="row wrap", align_items="flex-end", gap="0px",
        ),
    )
    action_group = widgets.HBox(
        [status, run_button],
        layout=widgets.Layout(
            flex="0 1 auto", flex_flow="row wrap", align_items="center",
            justify_content="flex-end", gap="10px",
        ),
    )
    action_row = widgets.HBox(
        [connect_pairs, action_group],
        layout=widgets.Layout(
            width="100%", flex_flow="row wrap", align_items="center",
            justify_content="space-between", gap="12px",
        ),
    )
    panel_body = widgets.VBox(
        [selection_row, action_row],
        layout=widgets.Layout(
            width="calc(100% - 40px)", margin="0 20px", padding="14px 0 18px", gap="16px",
        ),
    )
    control_card = widgets.VBox(
        [intro, panel_body],
        layout=widgets.Layout(
            width="100%", max_width="950px", padding="0", gap="0",
            border_radius="12px",
            background_color="#ffffff",
        ),
    )

    def render_result(points: pd.DataFrame, report: dict[str, Any]) -> None:
        label_a, label_b = report["host_a"], report["host_b"]
        host_symbols = {label_a: "circle", label_b: "diamond"}
        connector_x, connector_y = [], []
        if connect_pairs.value:
            points_by_pathogen: dict[str, dict[str, dict[str, Any]]] = {}
            for row in points.to_dict("records"):
                points_by_pathogen.setdefault(row["pathogen_id"], {})[row["class"]] = row
            for host_points in points_by_pathogen.values():
                first, second = host_points.get("First"), host_points.get("Second")
                if first is None or second is None:
                    continue
                connector_x.extend([first["IMSTE 1"], second["IMSTE 1"], None])
                connector_y.extend([first["IMSTE 2"], second["IMSTE 2"], None])

        figure_data = []
        if connector_x:
            figure_data.insert(0, go.Scatter(
                x=connector_x,
                y=connector_y,
                mode="lines",
                line={"color": "rgba(75, 91, 108, 0.55)", "width": 1, "dash": "dot"},
                hoverinfo="skip",
                showlegend=False,
                name="Same pathogen",
            ))
        colorbar_shown = False
        for host_label in (label_a, label_b):
            host_points = points.loc[points["host"] == host_label]
            if host_points.empty:
                continue
            positive_rows, gray_rows = [], []
            for row in host_points.to_dict("records"):
                try:
                    fraction = float(row["pathogenicity_fraction"])
                except (TypeError, ValueError):
                    fraction = float("nan")
                if math.isfinite(fraction) and fraction > 0:
                    positive_rows.append((row, min(1.0, max(0.0, fraction))))
                else:
                    gray_rows.append(row)

            def add_host_trace(rows: list[dict[str, Any]], *, gray: bool) -> None:
                nonlocal colorbar_shown
                if not rows:
                    return
                sizes = []
                customdata = []
                for row in rows:
                    try:
                        volume_score = float(row["estimated_pathogen_publications_index"])
                    except (TypeError, ValueError):
                        volume_score = 0.0
                    sizes.append(min(1.0, max(0.0, volume_score)) if math.isfinite(volume_score) else 0.0)
                    customdata.append([
                        row["host"], row["pair_status"], row["supporting_pmids"],
                        row["estimated_pathogen_publications_index"], row["pathogenicity_fraction"],
                    ])
                marker: dict[str, Any] = {
                    "symbol": host_symbols[host_label],
                    "size": sizes,
                    "sizemode": "area",
                    "sizeref": 1 / (24 * 24),
                    "sizemin": 4,
                    "line": {"color": "#111827", "width": 0.45},
                    "opacity": 0.16 if gray else 0.84,
                }
                if gray:
                    marker["color"] = "#aeb9c7"
                else:
                    marker.update({
                        "color": [fraction for _, fraction in positive_rows],
                        # Match Explorer's reversed YlGnBu pathogenicity scale.
                        "colorscale": "YlGnBu",
                        "reversescale": True,
                        "cmin": 0,
                        "cmax": 1,
                        "showscale": not colorbar_shown,
                    })
                    if not colorbar_shown:
                        marker["colorbar"] = {
                            "x": 1.07,
                            "thickness": 13,
                            "len": 0.72,
                            "tickformat": ".0%",
                            "title": {"text": "Pathogenicity fraction", "side": "right"},
                        }
                        colorbar_shown = True
                figure_data.append(go.Scattergl(
                    x=[row["IMSTE 1"] for row in rows],
                    y=[row["IMSTE 2"] for row in rows],
                    mode="markers",
                    name=host_label,
                    text=[row["pathogen"] for row in rows],
                    customdata=customdata,
                    hovertemplate=(
                        "<b>%{text}</b><br>%{customdata[0]}<br>"
                        "Pair status: %{customdata[1]}<br>Supporting PMIDs: %{customdata[2]}<br>"
                        "Publication volume score: %{customdata[3]:.3f}<br>"
                        "Estimated pathogenicity fraction: %{customdata[4]:.1%}<extra></extra>"
                    ),
                    marker=marker,
                    showlegend=False,
                ))

            add_host_trace([row for row, _ in positive_rows], gray=False)
            add_host_trace(gray_rows, gray=True)

        figure = go.Figure(data=figure_data)
        figure.update_layout(
            title={
                "text": f"{report['group'].title()} · shared-feature {report['embedding_method']} · {label_a} vs {label_b}",
                "x": 0.5, "y": 1.0, "xanchor": "center", "yanchor": "top",
            },
            autosize=True, height=650,
            paper_bgcolor="#ffffff", plot_bgcolor="#ffffff",
            font={"family": "Inter,-apple-system,BlinkMacSystemFont,Segoe UI,Arial,sans-serif", "color": "#29425f", "size": 11},
            xaxis={"title": "", "showticklabels": False, "ticks": "", "showline": False, "showgrid": True, "gridcolor": "#edf1f6", "zeroline": False},
            yaxis={"title": "", "showticklabels": False, "ticks": "", "showline": False, "showgrid": True, "gridcolor": "#edf1f6", "zeroline": False, "scaleanchor": "x", "scaleratio": 1},
            hovermode="closest",
            showlegend=False,
            margin={"l": 15, "r": 112, "t": 65, "b": 45},
        )
        shape_key = (
            '<div style="display:flex;align-items:center;justify-content:center;flex-wrap:wrap;gap:18px;'
            'margin:6px 0 12px;padding:7px 0;background:#fff;'
            'font-family:Inter,-apple-system,BlinkMacSystemFont,Segoe UI,Arial,sans-serif">'
            '<span style="margin-right:3px;color:#687b8d;font-size:10px;font-weight:750;letter-spacing:.09em;text-transform:uppercase">Host key</span>'
            f'<span style="display:inline-flex;align-items:center;gap:7px;color:#29425f;font-size:12px">'
            f'<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><circle cx="7" cy="7" r="4.5" fill="#dbeafe" stroke="#29425f" stroke-width="1.5"/></svg>'
            f'<strong>{escape(label_a)}</strong><span style="color:#8190a0;font-size:10px">Circle</span></span>'
            f'<span style="display:inline-flex;align-items:center;gap:7px;color:#29425f;font-size:12px">'
            f'<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><path d="M7 1.5 12.5 7 7 12.5 1.5 7Z" fill="#dbeafe" stroke="#29425f" stroke-width="1.5"/></svg>'
            f'<strong>{escape(label_b)}</strong><span style="color:#8190a0;font-size:10px">Diamond</span></span>'
            '</div>'
        )
        outputs = [html_output(shape_key)]
        figure_bundle = figure._repr_mimebundle_()
        if isinstance(figure_bundle, tuple):
            figure_data_bundle, figure_metadata = figure_bundle
        else:
            figure_data_bundle, figure_metadata = figure_bundle, {}
        plotly_payload = figure_data_bundle.get("application/vnd.plotly.v1+json")
        if isinstance(plotly_payload, dict):
            plotly_payload.setdefault("config", {})["responsive"] = True
        outputs.append({
            "output_type": "display_data",
            "data": figure_data_bundle,
            "metadata": figure_metadata,
        })
        summary_lines = [
            f"<div><strong>Applicable shared features:</strong> {report['shared_features']}</div>",
            f"<div><strong>Complete, variable features used:</strong> {len(report['used_features'])} — {escape(', '.join(report['used_features']))}</div>",
            (
                f"<div><strong>Embedding:</strong> {escape(report['embedding_method'])} "
                f"({report['n_points']} host-specific points from {report['n_pathogens']} pathogens × "
                f"{len(report['used_features'])} feature columns)</div>"
            ),
            (
                f"<div><strong>Pathogens plotted:</strong> {report['n_pathogens']} of "
                f"{report['n_confirmed_for_either_host']} with confirmed evidence for either host "
                f"({report['n_group_pathogens']} in the group); {report['excluded_no_shared_values']} "
                "confirmed pathogen(s) lacked shared-feature values.</div>"
            ),
        ]
        points_by_host = points.groupby("host")["pathogen_id"].nunique().to_dict()
        for host_label in (label_a, label_b):
            count = points_by_host.get(host_label, 0)
            suffix = " (no shared-feature values in this host's pair records)" if count == 0 else ""
            summary_lines.append(
                f"<div><strong>{escape(host_label)}:</strong> {count} plotted pathogen(s){escape(suffix)}</div>"
            )
        if report["n_msts"]:
            summary_lines.append(f"<div><strong>n_msts:</strong> {report['n_msts']}</div>")
        details_html = (
            '<details style="margin:10px 0 0;font:12px Inter,-apple-system,BlinkMacSystemFont,Segoe UI,Arial,sans-serif;color:#52677b">'
            '<summary style="cursor:pointer;font-weight:650;color:#52677b">Fit details</summary>'
            '<div style="display:grid;gap:5px;margin:8px 0 0;padding:10px 12px;'
            'border-left:2px solid #cbd9e1;background:#f8fafc;line-height:1.45">'
            + "".join(summary_lines)
            + '</div></details>'
        )
        outputs.append(html_output(details_html))
        # Assigning outputs replaces the previous plot in one widget state update.
        output.outputs = tuple(outputs)

    def set_busy(busy: bool) -> None:
        for control in (host_a, host_b, pathogen_group, connect_pairs):
            control.disabled = busy
        run_button.disabled = busy
        run_button.description = "Working…" if busy else "Fit 2D map"
        run_button.icon = "" if busy else "area-chart"
        if busy:
            status.value = (
                '<style>@keyframes host-feature-spin{to{transform:rotate(360deg)}}</style>'
                '<span style="display:inline-flex;align-items:center;gap:8px;color:#426276;font-size:12px;white-space:nowrap">'
                '<span style="width:16px;height:16px;box-sizing:border-box;border:2px solid #cbd5e1;'
                'border-top-color:#14877e;border-radius:50%;animation:host-feature-spin .75s linear infinite"></span>'
                '<strong>Fitting 2D map…</strong></span>'
            )

    async def fit_selection() -> None:
        try:
            points, report = await asyncio.to_thread(
                _fit_shared_host_map,
                registry, site_dist, metadata, host_a.value, host_b.value,
                pathogen_group.value, allowed_scopes,
            )
            if points is None:
                lines = [f"<strong>{escape(report['message'])}</strong>"]
                for key in ("shared_features", "pathogens_with_any_shared_values", "excluded_no_shared_values"):
                    if key in report:
                        lines.append(f"<div><strong>{escape(key.replace('_', ' '))}:</strong> {report[key]}</div>")
                output.outputs = (html_output(
                    '<div style="padding:10px 12px;border-left:3px solid #d5a642;background:#fffaf0;color:#755b24">'
                    + "".join(lines) + '</div>'
                ),)
                last_result.clear()
                status.value = '<span style="color:#9a5c18;font-size:12px">No map for this selection.</span>'
                return
            last_result.update(points=points, report=report)
            render_result(points, report)
            status.value = (
                f'<span style="color:#2f6f62;font-size:12px">Map ready · '
                f'{report["n_points"]} host-specific points</span>'
            )
        except Exception as error:
            last_result.clear()
            output.outputs = (html_output(
                '<div style="padding:10px 12px;border-left:3px solid #c2413b;background:#fff7f6;color:#8c2722">'
                f'<strong>Could not fit this comparison ({escape(type(error).__name__)}):</strong> '
                f'{escape(str(error))}</div>'
            ),)
            status.value = '<span style="color:#b42318;font-size:12px">Fit failed. Review the message below.</span>'
        finally:
            set_busy(False)

    def run(_):
        if run_button.disabled:
            return
        set_busy(True)
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            asyncio.run(fit_selection())
        else:
            loop.create_task(fit_selection())

    def update_connectors(change: dict[str, Any]) -> None:
        if change.get("name") != "value" or not last_result:
            return
        render_result(last_result["points"], last_result["report"])

    run_button.on_click(run)
    connect_pairs.observe(update_connectors, names="value")
    controls = widgets.VBox(
        [control_card, output],
        layout=widgets.Layout(width="100%", max_width="950px", gap="10px"),
    )
    return controls
