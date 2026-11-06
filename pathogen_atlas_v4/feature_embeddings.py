"""Fit reproducible IMSTE coordinates for host-specific pathogen feature profiles."""

from __future__ import annotations

from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


def fit_host_comparison_embeddings(
    species_rows: list[dict[str, Any]], *, host_a_id: str, host_b_id: str,
    pathogen_group: str, shared_feature_ids: list[str],
) -> dict[str, dict[str, Any]]:
    """Fit IMSTE and UMAP maps for a shared-feature comparison of two hosts.

    Only accepted/supported pathogen–host profiles with at least one observed
    shared feature participate. One set of complete feature columns is used for
    every plotted point; missing values are never imputed.
    """
    if host_a_id == host_b_id:
        return {}

    pathogen_points: list[dict[str, Any]] = []
    for taxon in species_rows:
        profiles = taxon.get("pair_profiles") or {}
        for host_id in (host_a_id, host_b_id):
            pair = profiles.get(host_id) or {}
            if str(pair.get("pair_selection_status", "")).casefold() not in {"accepted", "supported"}:
                continue
            values: dict[str, float] = {}
            for feature_id in shared_feature_ids:
                try:
                    value = float((pair.get("feature_values") or {}).get(feature_id))
                except (TypeError, ValueError):
                    continue
                if np.isfinite(value):
                    values[feature_id] = value
            if values:
                pathogen_points.append({
                    "pathogen_id": str(taxon["organism_id"]),
                    "host_id": host_id,
                    "feature_values": values,
                })

    if len(pathogen_points) < 2:
        return {}
    complete_feature_ids = [
        feature_id for feature_id in shared_feature_ids
        if all(feature_id in point["feature_values"] for point in pathogen_points)
    ]
    if not complete_feature_ids:
        return {}
    matrix = np.asarray([
        [point["feature_values"][feature_id] for feature_id in complete_feature_ids]
        for point in pathogen_points
    ], dtype=float)
    variable_mask = np.nanstd(matrix, axis=0) > 0
    used_feature_ids = [feature_id for feature_id, keep in zip(complete_feature_ids, variable_mask) if keep]
    matrix = matrix[:, variable_mask]
    if not used_feature_ids:
        return {}
    matrix = StandardScaler().fit_transform(matrix)

    fitted: dict[str, tuple[np.ndarray, int]] = {}
    if len(pathogen_points) == 2:
        from sklearn.decomposition import PCA

        coordinates = PCA(n_components=min(2, matrix.shape[1])).fit_transform(matrix)
        if coordinates.shape[1] == 1:
            coordinates = np.column_stack([coordinates[:, 0], np.zeros(len(pathogen_points))])
        fitted["IMSTE"] = (coordinates, 0)
    else:
        from imste import IteratedMinimumSpanningTreeEmbedder

        n_msts = min(10, max(2, len(pathogen_points) - 3))
        estimator = IteratedMinimumSpanningTreeEmbedder(
            n_msts=n_msts, n_components=2, n_epochs=300, batch_size=32768,
            learning_rate=0.05, negative_ratio=5, lambda_rep=0.5,
            random_state=42, epsilon=0.0001, device="cpu", logistic_margin=1.0,
            logistic_temperature=0.5, mst_method="famst",
            mst_neighbors=min(15, len(pathogen_points) - 1), mst_inter_component_edges=5,
            mst_max_neighbors=None, mst_representatives_per_component=10,
        )
        try:
            fitted["IMSTE"] = (np.asarray(estimator.fit_transform(matrix), dtype=float), n_msts)
        except (ValueError, RuntimeError, FloatingPointError):
            pass

        from umap import UMAP

        try:
            umap = UMAP(
                n_components=2, n_neighbors=min(15, len(pathogen_points) - 1),
                min_dist=0.1, metric="euclidean", random_state=42,
                transform_seed=42, n_jobs=1,
            )
            fitted["UMAP"] = (np.asarray(umap.fit_transform(matrix), dtype=float), 0)
        except (ValueError, RuntimeError, FloatingPointError):
            pass

    results: dict[str, dict[str, Any]] = {}
    for method, (coordinates, n_msts) in fitted.items():
        if coordinates.shape != (len(pathogen_points), 2) or not np.isfinite(coordinates).all():
            continue
        coordinate_std = coordinates.std(axis=0)
        if float(np.max(coordinate_std)) < 1e-6:
            continue
        points = [
            {"pathogen_id": point["pathogen_id"], "host_id": point["host_id"],
             "x": float(x), "y": float(y)}
            for point, (x, y) in zip(pathogen_points, coordinates)
        ]
        results[method] = {
            "points": points,
            "metadata": {
                "host_ids": [host_a_id, host_b_id], "pathogen_group": pathogen_group,
                "embedding": method, "n_points": len(points),
                "n_pathogens": len({point["pathogen_id"] for point in points}),
                "shared_feature_ids": shared_feature_ids, "used_feature_ids": used_feature_ids,
                "n_msts": n_msts,
            },
        }
    return results

def fit_feature_embeddings(
    profiles: list[dict[str, Any]], *, host_id: str, pathogen_group: str,
) -> dict[str, tuple[dict[str, tuple[float, float]], dict[str, Any]]]:
    """Fit IMSTE and UMAP maps to accepted profiles with complete feature values."""
    selected = [profile for profile in profiles
                if str(profile.get("pair_selection_status", "")).casefold() in {"accepted", "supported"}]
    observed_sets = [
        {str(feature_id) for feature_id, value in (profile.get("feature_values") or {}).items()
         if value is not None}
        for profile in selected
    ]
    feature_ids = sorted(set.intersection(*observed_sets)) if observed_sets else []
    complete = [profile for profile, observed in zip(selected, observed_sets)
                if feature_ids and set(feature_ids).issubset(observed)]
    if len(complete) < 3:
        raise ValueError(f"IMSTE needs at least three complete {host_id} × {pathogen_group} profiles; found {len(complete)}")

    matrix = np.asarray([[float((profile.get("feature_values") or {})[feature_id])
                          for feature_id in feature_ids] for profile in complete], dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError(f"IMSTE input contains non-finite feature values for {host_id} × {pathogen_group}")
    variable_mask = np.nanstd(matrix, axis=0) > 0
    matrix = matrix[:, variable_mask]
    feature_ids = [feature_id for feature_id, keep in zip(feature_ids, variable_mask) if keep]
    if not feature_ids:
        raise ValueError(f"IMSTE has no variable feature values for {host_id} × {pathogen_group}")
    matrix = StandardScaler().fit_transform(matrix)

    from imste import IteratedMinimumSpanningTreeEmbedder
    from umap import UMAP

    # Very small cohorts can collapse with the rice default (n_species - 1)
    # trees. Keep at least two trees and leave a small margin below cohort size.
    n_msts = min(10, max(2, len(complete) - 3))
    estimators = {
        "IMSTE": (
            IteratedMinimumSpanningTreeEmbedder(
                n_msts=n_msts, n_components=2, n_epochs=300, batch_size=32768,
                learning_rate=0.05, negative_ratio=5, lambda_rep=0.5,
                random_state=42, epsilon=0.0001, device="cpu", logistic_margin=1.0,
                logistic_temperature=0.5, mst_method="famst",
                mst_neighbors=min(15, len(complete) - 1), mst_inter_component_edges=5,
                mst_max_neighbors=None, mst_representatives_per_component=10,
            ),
            {"estimator": "IteratedMinimumSpanningTreeEmbedder",
             "imste_version": version("imste"), "n_msts": n_msts,
             "n_epochs": 300, "mst_method": "famst",
             "mst_neighbors": min(15, len(complete) - 1), "device": "cpu"},
        ),
        "UMAP": (
            UMAP(n_components=2, n_neighbors=min(15, len(complete) - 1),
                 min_dist=0.1, metric="euclidean", random_state=42,
                 transform_seed=42, n_jobs=1),
            {"estimator": "umap-learn.UMAP", "umap_version": version("umap-learn"),
             "n_neighbors": min(15, len(complete) - 1), "min_dist": 0.1,
             "metric": "euclidean", "n_jobs": 1},
        ),
    }
    results = {}
    failures = []
    for method, (estimator, parameters) in estimators.items():
        try:
            coordinates_array = np.asarray(estimator.fit_transform(matrix), dtype=float)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            failures.append(f"{method}: {exc}")
            continue
        if coordinates_array.shape != (len(complete), 2) or not np.isfinite(coordinates_array).all():
            failures.append(f"{method}: returned invalid coordinates")
            continue
        coordinate_std = coordinates_array.std(axis=0)
        if float(np.max(coordinate_std)) < 1e-6:
            failures.append(f"{method}: returned collapsed coordinates")
            continue
        coordinates = {
            str(profile["pathogen_id"]): (float(point[0]), float(point[1]))
            for profile, point in zip(complete, coordinates_array)
        }
        metadata = {
            "embedding": method, "host_id": host_id,
            "pathogen_group": pathogen_group, "n_species": len(complete),
            "n_features": len(feature_ids), "feature_ids": feature_ids,
            "input": "accepted pair feature yes fractions on features scored for every fitted profile",
            "scaling": "StandardScaler", "random_state": 42,
            "coordinate_std": coordinate_std.tolist(),
            "coordinate_range": [[float(coordinates_array[:, i].min()),
                                  float(coordinates_array[:, i].max())] for i in range(2)],
            **parameters,
        }
        results[method] = (coordinates, metadata)
    if not results:
        details = "; ".join(failures) or "no embedding method completed"
        raise ValueError(f"No stable 2D map for {host_id} × {pathogen_group}: {details}")
    return results


def fit_feature_embedding(
    profiles: list[dict[str, Any]], *, host_id: str, pathogen_group: str,
) -> tuple[dict[str, tuple[float, float]], dict[str, Any]]:
    """Compatibility wrapper returning the IMSTE map."""
    return fit_feature_embeddings(profiles, host_id=host_id, pathogen_group=pathogen_group)["IMSTE"]
