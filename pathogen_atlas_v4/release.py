"""Pair-aware release validation and literature-evidence metrics."""

from __future__ import annotations

import pandas as pd


def validate_pair_release(pathogens: pd.DataFrame, hosts: pd.DataFrame, articles: pd.DataFrame,
                          evidence: pd.DataFrame) -> None:
    required = {"pathogen_id", "host_id", "pmid", "decision", "provenance_status", "evidence_span"}
    if not required.issubset(evidence.columns):
        raise ValueError(f"Pair evidence missing columns: {sorted(required - set(evidence.columns))}")
    pathogen_ids = set(pathogens["pathogen_id"].astype(str))
    host_ids = set(hosts["host_id"].astype(str))
    pmids = set(articles["pmid"].astype(str))
    if not set(evidence["pathogen_id"].astype(str)) <= pathogen_ids:
        raise ValueError("Pair evidence references unknown pathogen IDs")
    if not set(evidence["host_id"].astype(str)) <= host_ids:
        raise ValueError("Pair evidence references unknown host IDs")
    if not set(evidence["pmid"].astype(str)) <= pmids:
        raise ValueError("Pair evidence references unknown PMIDs")
    positives = evidence["decision"].astype(str).eq("supported")
    missing_span = evidence.loc[positives, "evidence_span"].fillna("").astype(str).str.strip().eq("")
    system1 = evidence["provenance_status"].astype(str).eq("system1_jev")
    if (missing_span & ~system1.loc[positives]).any():
        raise ValueError("Every supported text-generating decision must carry a supporting evidence span")
    if (positives & system1).any():
        if "p_yes" not in evidence or evidence.loc[positives & system1, "p_yes"].isna().any():
            raise ValueError("Every supported Jev decision must retain its yes probability")
        p_yes = pd.to_numeric(evidence.loc[positives & system1, "p_yes"], errors="coerce")
        if p_yes.isna().any() or p_yes.lt(0.5).any():
            raise ValueError("Supported Jev decisions must meet the 0.5 yes-probability threshold")
        if evidence.loc[positives & system1, "abstract_hash"].fillna("").astype(str).str.strip().eq("").any():
            raise ValueError("Every supported Jev decision must identify its screened abstract")
    if evidence.loc[positives, "provenance_status"].astype(str).isin({"v4_native", "system1_jev"}).any():
        if "relationship_context" not in evidence or not evidence.loc[positives, "relationship_context"].eq("disease").all():
            raise ValueError("Native supported pairs must use the disease relationship context")


def evidence_metrics(matches: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Report literature search/acceptance counts by pair without prevalence claims."""
    if matches.empty:
        return pd.DataFrame(columns=["pathogen_id", "host_id", "n_search_matches", "n_fetched", "n_supported"])
    grouped = matches.groupby(["pathogen_id", "host_id"], dropna=False)
    result = grouped.agg(n_search_matches=("pmid", "nunique"), n_fetched=("selected_for_fetch", "sum")).reset_index()
    positives = evidence.loc[evidence["decision"].astype(str).eq("supported")]
    counts = positives.groupby(["pathogen_id", "host_id"])["pmid"].nunique().rename("n_supported").reset_index()
    return result.merge(counts, on=["pathogen_id", "host_id"], how="left").fillna({"n_supported": 0})
