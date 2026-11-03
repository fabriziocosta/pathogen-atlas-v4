"""Fit reproducible IMSTE coordinates for host-specific pathogen feature profiles."""

from __future__ import annotations

from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


def fit_feature_embedding(
    profiles: list[dict[str, Any]], *, host_id: str, pathogen_group: str,
) -> tuple[dict[str, tuple[float, float]], dict[str, Any]]:
    """Fit IMSTE to accepted profiles with complete, variable feature values."""
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

    # Very small cohorts can collapse with the rice default (n_species - 1)
    # trees. Keep at least two trees and leave a small margin below cohort size.
    n_msts = min(10, max(2, len(complete) - 3))
    embedder = IteratedMinimumSpanningTreeEmbedder(
        n_msts=n_msts, n_components=2, n_epochs=300, batch_size=32768,
        learning_rate=0.05, negative_ratio=5, lambda_rep=0.5,
        random_state=42, epsilon=0.0001, device="cpu", logistic_margin=1.0,
        logistic_temperature=0.5, mst_method="famst",
        mst_neighbors=min(15, len(complete) - 1), mst_inter_component_edges=5,
        mst_max_neighbors=None, mst_representatives_per_component=10,
    )
    coordinates_array = np.asarray(embedder.fit_transform(matrix), dtype=float)
    if coordinates_array.shape != (len(complete), 2) or not np.isfinite(coordinates_array).all():
        raise ValueError(f"IMSTE returned invalid coordinates for {host_id} × {pathogen_group}")
    coordinate_std = coordinates_array.std(axis=0)
    if float(np.max(coordinate_std)) < 1e-6:
        raise ValueError(f"IMSTE returned collapsed coordinates for {host_id} × {pathogen_group}")
    coordinates = {
        str(profile["pathogen_id"]): (float(point[0]), float(point[1]))
        for profile, point in zip(complete, coordinates_array)
    }
    metadata = {
        "embedding": "IMSTE", "estimator": "IteratedMinimumSpanningTreeEmbedder",
        "imste_version": version("imste"), "host_id": host_id,
        "pathogen_group": pathogen_group, "n_species": len(complete),
        "n_features": len(feature_ids), "feature_ids": feature_ids,
        "input": "accepted pair feature yes fractions on features scored for every fitted profile",
        "scaling": "StandardScaler", "n_msts": n_msts, "n_epochs": 300,
        "mst_method": "famst", "mst_neighbors": min(15, len(complete) - 1),
        "random_state": 42, "device": "cpu",
        "coordinate_std": coordinate_std.tolist(),
        "coordinate_range": [[float(coordinates_array[:, i].min()),
                              float(coordinates_array[:, i].max())] for i in range(2)],
    }
    return coordinates, metadata
