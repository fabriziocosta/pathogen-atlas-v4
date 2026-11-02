"""Host-qualified Pathogenic Interest Index from sampled pair evidence."""
from __future__ import annotations

import numpy as np
import pandas as pd

INDEX_VERSION = "pathogenic-interest-index-v4.2-sampled-extrapolated"

OUTPUT_COLUMNS = [
    "pathogen_id", "host_id", "accepted_name", "pathogen_group",
    "n_positive_records", "n_abstracts_screened", "n_pubmed_search_records",
    "estimated_pathogen_pmids", "log_estimated_pathogen_pmids",
    "pathogenic_interest_index", "index_version", "query_hash", "query_version",
]


def build_pathogenic_interest_index(publication_counts: pd.DataFrame,
                                    accepted_pairs: pd.DataFrame,
                                    pathogen_taxonomy: pd.DataFrame) -> pd.DataFrame:
    """Extrapolate the sampled positive rate to each accepted pair's PubMed total.

    As in V3, estimated pathogen publications are the sampled positive rate
    multiplied by the organism-only PubMed search count. The log estimate is
    normalized within pathogen group. Total-count searches must cover only the
    unique pathogens represented by at least one accepted host pair.
    """
    required_counts = {"pathogen_id", "pathogen_group", "n_pubmed_search_records", "query_hash", "query_version"}
    required_pairs = {"pathogen_id", "host_id", "selected_as_pathogen_for_host", "selection_status",
                      "n_sampled_records_with_abstract", "n_completed_screenings", "n_positive_records"}
    if not required_counts.issubset(publication_counts.columns):
        raise ValueError(f"Publication counts missing columns: {sorted(required_counts - set(publication_counts.columns))}")
    if not required_pairs.issubset(accepted_pairs.columns):
        raise ValueError(f"Accepted pair table missing columns: {sorted(required_pairs - set(accepted_pairs.columns))}")
    accepted = accepted_pairs.loc[
        accepted_pairs.selected_as_pathogen_for_host.fillna(False).astype(bool)
        & accepted_pairs.selection_status.astype(str).eq("accepted")
    ].copy()
    if accepted.duplicated(["pathogen_id", "host_id"]).any():
        raise ValueError("Accepted pair table must contain one row per pathogen_id and host_id")
    if accepted.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    accepted["pathogen_id"] = accepted.pathogen_id.astype(str)
    accepted["host_id"] = accepted.host_id.astype(str)
    positive = pd.to_numeric(accepted.n_positive_records, errors="coerce")
    screened = pd.to_numeric(accepted.n_completed_screenings, errors="coerce")
    available = pd.to_numeric(accepted.n_sampled_records_with_abstract, errors="coerce")
    if positive.isna().any() or screened.isna().any() or available.isna().any():
        raise ValueError("Accepted pairs require numeric positive, screened, and abstract-available counts")
    if screened.ne(available).any():
        raise ValueError("Finish screening every sampled abstract before estimating the PII")
    if positive.le(0).any() or positive.gt(screened).any() or screened.le(0).any():
        raise ValueError("Accepted pairs require 1..n screened positive abstracts")
    accepted["n_positive_records"] = positive.astype("int64")
    accepted["n_abstracts_screened"] = screened.astype("int64")
    pathogen_ids = set(accepted.pathogen_id)
    count_ids = set(publication_counts.pathogen_id.astype(str))
    missing = sorted(pathogen_ids - count_ids)
    extra = sorted(count_ids - pathogen_ids)
    if missing or extra:
        raise ValueError(f"Accepted pathogen count coverage mismatch: missing={missing[:10]}, extra={extra[:10]}")
    counts = publication_counts.copy()
    if counts.pathogen_id.astype(str).duplicated().any():
        raise ValueError("Publication counts must contain exactly one row per accepted pathogen ID")
    known = pathogen_taxonomy[["pathogen_id", "accepted_name", "pathogen_group"]].drop_duplicates("pathogen_id")
    result = counts.merge(known, on=["pathogen_id", "pathogen_group"], how="left", validate="one_to_one")
    if result.accepted_name.isna().any():
        raise ValueError("Publication count refers to a pathogen absent from taxonomy")
    n = pd.to_numeric(result.n_pubmed_search_records, errors="coerce")
    if n.isna().any() or n.lt(0).any():
        raise ValueError("PubMed publication counts must be nonnegative integers")
    if not np.equal(n, np.floor(n)).all():
        raise ValueError("PubMed publication counts must be integers")
    counts["n_pubmed_search_records"] = n.astype("int64")
    result = accepted.merge(counts, on="pathogen_id", how="left", validate="many_to_one")
    result = result.merge(known, on="pathogen_id", how="left", validate="many_to_one", suffixes=("", "_taxonomy"))
    if result.accepted_name.isna().any() or result.pathogen_group.isna().any():
        raise ValueError("Accepted pair or publication count refers to a pathogen absent from taxonomy")
    if not result.pathogen_group.eq(result.pathogen_group_taxonomy).all():
        raise ValueError("Publication count pathogen group does not match taxonomy")
    result["estimated_pathogen_pmids"] = (
        result.n_positive_records / result.n_abstracts_screened * result.n_pubmed_search_records
    )
    result["log_estimated_pathogen_pmids"] = 0.0
    positive_estimates = result.estimated_pathogen_pmids.gt(0)
    result.loc[positive_estimates, "log_estimated_pathogen_pmids"] = np.log(
        result.loc[positive_estimates, "estimated_pathogen_pmids"]
    )
    result["pathogenic_interest_index"] = 0.0
    for _, indices in result.groupby("pathogen_group", sort=False).groups.items():
        maximum = float(result.loc[indices, "log_estimated_pathogen_pmids"].max())
        if maximum > 0:
            result.loc[indices, "pathogenic_interest_index"] = result.loc[
                indices, "log_estimated_pathogen_pmids"
            ] / maximum
    result["index_version"] = INDEX_VERSION
    return result[OUTPUT_COLUMNS].sort_values(["pathogen_group", "pathogen_id", "host_id"]).reset_index(drop=True)
