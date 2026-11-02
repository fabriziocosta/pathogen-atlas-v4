"""Host-pair Pathogenic Interest Index from sampled evidence and PubMed counts."""
from __future__ import annotations

import pandas as pd

INDEX_VERSION = "pathogenic-interest-index-v4.3-host-pair-fraction"

OUTPUT_COLUMNS = [
    "pathogen_id", "host_id", "accepted_name", "pathogen_group",
    "n_positive_records", "n_sampled_records", "n_abstracts_screened", "n_pair_search_records",
    "estimated_pathogen_publications", "n_pubmed_search_records",
    "pathogenicity_fraction", "pathogenic_interest_index", "index_version",
    "query_hash", "query_version",
]


def build_pathogenic_interest_index(publication_counts: pd.DataFrame,
                                    accepted_pairs: pd.DataFrame,
                                    pathogen_taxonomy: pd.DataFrame) -> pd.DataFrame:
    """Estimate the pathogenicity fraction for each accepted pathogen–host pair.

    Extrapolate the observed positive rate across all selected sample records
    to all PubMed records returned by the pair query, then divide that estimate
    by the count-only organism-wide PubMed total. Organism-wide queries must
    cover exactly the unique pathogens accepted in at least one selected host.
    """
    required_counts = {"pathogen_id", "pathogen_group", "n_pubmed_search_records", "query_hash", "query_version"}
    required_pairs = {"pathogen_id", "host_id", "selected_as_pathogen_for_host", "selection_status",
                      "n_sampled_records", "n_sampled_records_with_abstract", "n_completed_screenings", "n_positive_records",
                      "n_pair_search_records"}
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
    positives = pd.to_numeric(accepted.n_positive_records, errors="coerce")
    sampled = pd.to_numeric(accepted.n_sampled_records, errors="coerce")
    screened = pd.to_numeric(accepted.n_completed_screenings, errors="coerce")
    abstracts_available = pd.to_numeric(accepted.n_sampled_records_with_abstract, errors="coerce")
    pair_total = pd.to_numeric(accepted.n_pair_search_records, errors="coerce")
    if positives.isna().any() or sampled.isna().any() or screened.isna().any() or abstracts_available.isna().any() or pair_total.isna().any():
        raise ValueError("Accepted pairs require numeric sample, screened, positive, and pair-search counts")
    if sampled.lt(1).any() or abstracts_available.gt(sampled).any() or screened.gt(sampled).any():
        raise ValueError("Sample and screened-abstract counts must fit within the selected PMID sample")
    if screened.ne(abstracts_available).any():
        raise ValueError("Finish screening every sampled abstract before estimating the PII")
    if positives.le(0).any() or positives.gt(sampled).any() or screened.le(0).any():
        raise ValueError("Accepted pairs require positive evidence within a nonempty PMID sample")
    if pair_total.lt(positives).any() or pair_total.le(0).any():
        raise ValueError("Pair-search records must cover every positive sampled abstract")
    accepted["n_positive_records"] = positives.astype("int64")
    accepted["n_sampled_records"] = sampled.astype("int64")
    accepted["n_abstracts_screened"] = screened.astype("int64")
    accepted["n_pair_search_records"] = pair_total.astype("int64")

    pathogen_ids = set(accepted.pathogen_id)
    count_ids = set(publication_counts.pathogen_id.astype(str))
    missing, extra = sorted(pathogen_ids - count_ids), sorted(count_ids - pathogen_ids)
    if missing or extra:
        raise ValueError(f"Accepted pathogen count coverage mismatch: missing={missing[:10]}, extra={extra[:10]}")
    counts = publication_counts.copy()
    counts["pathogen_id"] = counts.pathogen_id.astype(str)
    if counts.pathogen_id.duplicated().any():
        raise ValueError("Publication counts must contain exactly one row per accepted pathogen ID")
    counts["n_pubmed_search_records"] = pd.to_numeric(counts.n_pubmed_search_records, errors="coerce")
    count_values = counts.n_pubmed_search_records
    if count_values.isna().any() or count_values.lt(0).any() or not (count_values == count_values.astype("int64")).all():
        raise ValueError("PubMed publication counts must be nonnegative integers")
    counts["n_pubmed_search_records"] = count_values.astype("int64")
    if counts.n_pubmed_search_records.eq(0).any():
        raise ValueError("An accepted pathogen must have a nonzero organism-wide PubMed count")

    known = pathogen_taxonomy[["pathogen_id", "accepted_name", "pathogen_group"]].drop_duplicates("pathogen_id")
    known["pathogen_id"] = known.pathogen_id.astype(str)
    result = accepted.merge(counts, on="pathogen_id", how="left", validate="many_to_one")
    result = result.merge(known, on="pathogen_id", how="left", validate="many_to_one", suffixes=("", "_taxonomy"))
    if result.accepted_name.isna().any() or result.pathogen_group.isna().any():
        raise ValueError("Accepted pair or publication count refers to a pathogen absent from taxonomy")
    if not result.pathogen_group.eq(result.pathogen_group_taxonomy).all():
        raise ValueError("Publication count pathogen group does not match taxonomy")

    result["estimated_pathogen_publications"] = (
        result.n_positive_records / result.n_sampled_records * result.n_pair_search_records
    )
    result["pathogenicity_fraction"] = (
        result.estimated_pathogen_publications / result.n_pubmed_search_records
    )
    result["pathogenic_interest_index"] = result.pathogenicity_fraction
    result["index_version"] = INDEX_VERSION
    return result[OUTPUT_COLUMNS].sort_values(["pathogen_group", "pathogen_id", "host_id"]).reset_index(drop=True)
