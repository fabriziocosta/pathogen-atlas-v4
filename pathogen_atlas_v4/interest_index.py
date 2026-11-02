"""Host-qualified Pathogenic Interest Index from organism publication counts."""
from __future__ import annotations

import numpy as np
import pandas as pd

INDEX_VERSION = "pathogenic-interest-index-v4.1-log1p-publication-count"


def build_pathogenic_interest_index(publication_counts: pd.DataFrame,
                                    accepted_pairs: pd.DataFrame,
                                    pathogen_taxonomy: pd.DataFrame) -> pd.DataFrame:
    """Normalize PubMed publication volume for pathogens accepted in at least one host.

    The index is log1p(total distinct PubMed results) normalized to [0, 1]
    within pathogen group. Only organisms represented by an accepted pair are
    allowed into the output table.
    """
    required_counts = {"pathogen_id", "pathogen_group", "n_pubmed_search_records", "query_hash", "query_version"}
    required_pairs = {"pathogen_id", "host_id", "selected_as_pathogen_for_host", "selection_status"}
    if not required_counts.issubset(publication_counts.columns):
        raise ValueError(f"Publication counts missing columns: {sorted(required_counts - set(publication_counts.columns))}")
    if not required_pairs.issubset(accepted_pairs.columns):
        raise ValueError(f"Accepted pair table missing columns: {sorted(required_pairs - set(accepted_pairs.columns))}")
    accepted = accepted_pairs.loc[
        accepted_pairs.selected_as_pathogen_for_host.fillna(False).astype(bool)
        & accepted_pairs.selection_status.astype(str).eq("accepted"), "pathogen_id"
    ].astype(str).unique()
    if not len(accepted):
        return pd.DataFrame(columns=["pathogen_id", "pathogen_group", "n_pubmed_search_records", "log_publication_count", "pathogenic_interest_index", "index_version", "query_hash", "query_version"])
    counts = publication_counts.loc[publication_counts.pathogen_id.astype(str).isin(set(accepted))].copy()
    if counts.pathogen_id.astype(str).duplicated().any():
        raise ValueError("Publication counts must contain exactly one row per accepted pathogen ID")
    missing = sorted(set(accepted) - set(counts.pathogen_id.astype(str)))
    extra = sorted(set(counts.pathogen_id.astype(str)) - set(accepted))
    if missing or extra:
        raise ValueError(f"Accepted pathogen count coverage mismatch: missing={missing[:10]}, extra={extra[:10]}")
    known = pathogen_taxonomy[["pathogen_id", "accepted_name", "pathogen_group"]].drop_duplicates("pathogen_id")
    result = counts.merge(known, on=["pathogen_id", "pathogen_group"], how="left", validate="one_to_one")
    if result.accepted_name.isna().any():
        raise ValueError("Publication count refers to a pathogen absent from taxonomy")
    n = pd.to_numeric(result.n_pubmed_search_records, errors="coerce")
    if n.isna().any() or n.lt(0).any():
        raise ValueError("PubMed publication counts must be nonnegative integers")
    result["n_pubmed_search_records"] = n.astype("int64")
    result["log_publication_count"] = np.log1p(result.n_pubmed_search_records.astype(float))
    result["pathogenic_interest_index"] = 0.0
    for _, indices in result.groupby("pathogen_group", sort=False).groups.items():
        maximum = float(result.loc[indices, "log_publication_count"].max())
        if maximum > 0:
            result.loc[indices, "pathogenic_interest_index"] = result.loc[indices, "log_publication_count"] / maximum
    result["index_version"] = INDEX_VERSION
    return result[["pathogen_id", "accepted_name", "pathogen_group", "n_pubmed_search_records",
                   "log_publication_count", "pathogenic_interest_index", "index_version",
                   "query_hash", "query_version"]].sort_values(["pathogen_group", "pathogen_id"]).reset_index(drop=True)
