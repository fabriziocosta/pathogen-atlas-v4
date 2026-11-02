"""Resumable, pair-scoped PubMed acquisition through NCBI E-utilities."""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from pathogen_atlas_v4.pubmed import build_pair_query, build_pathogen_query

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200


def _xml(response: requests.Response) -> ET.Element:
    response.raise_for_status()
    return ET.fromstring(response.content)


def _search(session: requests.Session, query: str, *, email: str, api_key: str | None, logger: logging.Logger):
    """Return a complete (at most 10,000) PubMed UID page and result count."""
    params = {"db": "pubmed", "term": query, "usehistory": "y", "retmax": 10000,
              "email": email, "tool": "pathogen_atlas_v4"}
    if api_key:
        params["api_key"] = api_key
    for attempt in range(7):
        try:
            root = _xml(session.post(f"{EUTILS}/esearch.fcgi", data=params, timeout=90))
            error = root.findtext("Error") or root.findtext("ERROR")
            if error:
                raise RuntimeError(f"NCBI ESearch error: {error}")
            count = int(root.findtext("Count", "0"))
            pmids = [node.text.strip() for node in root.findall(".//IdList/Id") if node.text]
            time.sleep(0.11 if api_key else 0.36)
            return count, pmids
        except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
            if attempt == 6:
                raise
            delay = min(2 ** attempt, 60)
            logger.warning("ESearch retry %d/6 after %s; sleeping %ds", attempt + 1, type(exc).__name__, delay)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _fetch_pmids(session: requests.Session, pmids: list[str], *, email: str, api_key: str | None, logger: logging.Logger):
    """Fetch metadata only for selected sample PMIDs, reusing the PMID cache."""
    articles = []
    params = {"db": "pubmed", "retmode": "xml", "email": email, "tool": "pathogen_atlas_v4"}
    if api_key:
        params["api_key"] = api_key
    for start in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[start:start + BATCH_SIZE]
        params["id"] = ",".join(batch)
        for attempt in range(7):
            try:
                root = _xml(session.post(f"{EUTILS}/efetch.fcgi", data=params, timeout=120))
                rows = _parse_articles(root)
                if len(rows) != len(batch):
                    raise RuntimeError(f"EFetch returned {len(rows)} of {len(batch)} requested PubMed records")
                articles.extend(rows)
                break
            except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
                if attempt == 6:
                    raise
                delay = min(2 ** attempt, 60)
                logger.warning("EFetch retry at sample offset %d after %s; sleeping %ds", start, type(exc).__name__, delay)
                time.sleep(delay)
        time.sleep(0.11 if api_key else 0.36)
    return articles


def _parse_articles(root: ET.Element) -> list[dict[str, str]]:
    rows = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID", "").strip()
        medline = article.find(".//MedlineCitation")
        art = medline.find(".//Article") if medline is not None else None
        if not pmid or art is None:
            continue
        title_node = art.find("ArticleTitle")
        title = "".join(title_node.itertext()).strip() if title_node is not None else ""
        journal = art.findtext(".//Journal/Title", "")
        year = art.findtext(".//JournalIssue/PubDate/Year", "") or art.findtext(".//JournalIssue/PubDate/MedlineDate", "")
        doi = ""
        for ident in art.findall(".//ELocationID"):
            if ident.attrib.get("EIdType", "").casefold() == "doi":
                doi = "".join(ident.itertext()).strip()
                break
        rows.append({"pmid": pmid, "title": title, "journal": journal, "publication_year": year, "doi": doi})
    return rows


def _date_range_query(query: str, start: str, end: str) -> str:
    return f'({query}) AND ("{start}"[DP] : "{end}"[DP])'


def _all_pmids(session, query: str, *, email, api_key, logger):
    """Collect all UIDs, segmenting large result sets before ESearch's 10,000 UID limit."""
    count, pmids = _search(session, query, email=email, api_key=api_key, logger=logger)
    if count <= 10000:
        if len(pmids) != count:
            raise RuntimeError(f"ESearch returned {len(pmids)} of {count} PubMed UIDs")
        return pmids, []
    logger.warning("Search returned %d records; segmenting by publication year to avoid PubMed's 10,000 UID retrieval ceiling", count)
    from datetime import date
    import calendar
    segmented_pmids, segment_counts = [], []
    for year in range(1800, date.today().year + 1):
        year_query = _date_range_query(query, f"{year}/01/01", f"{year}/12/31")
        year_count, year_pmids = _search(session, year_query, email=email, api_key=api_key, logger=logger)
        if year_count <= 10000:
            if len(year_pmids) != year_count:
                raise RuntimeError(f"ESearch returned {len(year_pmids)} of {year_count} UIDs for {year}")
            segmented_pmids.extend(year_pmids)
            segment_counts.append({"period": str(year), "count": year_count})
            continue
        for month in range(1, 13):
            last_day = calendar.monthrange(year, month)[1]
            month_query = _date_range_query(query, f"{year}/{month:02d}/01", f"{year}/{month:02d}/{last_day:02d}")
            month_count, month_pmids = _search(session, month_query, email=email, api_key=api_key, logger=logger)
            if month_count <= 10000:
                if len(month_pmids) != month_count:
                    raise RuntimeError(f"ESearch returned {len(month_pmids)} of {month_count} UIDs for {year}-{month:02d}")
                segmented_pmids.extend(month_pmids)
                segment_counts.append({"period": f"{year}-{month:02d}", "count": month_count})
                continue
            for day in range(1, last_day + 1):
                daily_query = _date_range_query(query, f"{year}/{month:02d}/{day:02d}", f"{year}/{month:02d}/{day:02d}")
                daily_count, daily_pmids = _search(session, daily_query, email=email, api_key=api_key, logger=logger)
                if daily_count > 10000:
                    raise RuntimeError(f"A single publication date still has {daily_count} PubMed records; cannot guarantee complete retrieval")
                if len(daily_pmids) != daily_count:
                    raise RuntimeError(f"ESearch returned {len(daily_pmids)} of {daily_count} UIDs for {year}-{month:02d}-{day:02d}")
                segmented_pmids.extend(daily_pmids)
                segment_counts.append({"period": f"{year}-{month:02d}-{day:02d}", "count": daily_count})
    unique = sorted(set(segmented_pmids), key=lambda x: (0, int(x)) if x.isdigit() else (1, x))
    if len(unique) != count:
        raise RuntimeError(f"Date-segmented ESearch recovered {len(unique)} of {count} PubMed UIDs; stopping to avoid incomplete results")
    return unique, segment_counts


def select_stable_sample(pmids: Iterable[str], *, maximum: int | None = 200,
                         oldest_count: int = 30, seed: int = 42) -> dict[str, str]:
    """Select oldest PMID identifiers plus a nested, deterministic hash-ranked sample.

    Raising ``maximum`` preserves every previously selected PMID. The random
    priority is a stable SHA-256 rank, so cap 400 extends cap 200 by 200.
    """
    values = sorted({str(pmid) for pmid in pmids}, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))
    if maximum is not None and int(maximum) < 1:
        raise ValueError("maximum must be positive or None")
    if int(oldest_count) < 0:
        raise ValueError("oldest_count cannot be negative")
    if maximum is None or len(values) <= int(maximum):
        return {pmid: "all" for pmid in values}
    cap = int(maximum)
    oldest = values[:min(int(oldest_count), cap)]
    oldest_set = set(oldest)
    remainder = [pmid for pmid in values if pmid not in oldest_set]
    random_budget = cap - len(oldest)
    random_order = sorted(remainder, key=lambda pmid: hashlib.sha256(f"{seed}|{pmid}".encode()).hexdigest())
    return {**{pmid: "oldest" for pmid in oldest}, **{pmid: "random" for pmid in random_order[:random_budget]}}


def acquire_pathogen_literature(pathogen_ids: Iterable[str], pathogens: pd.DataFrame, aliases: pd.DataFrame,
                                output_dir: str | Path, *, query_version: str, email: str,
                                sample_maximum: int | None = 200, oldest_count: int = 30,
                                seed: int = 42, sample_policy_version: str = "stable-prefix-hash-v1",
                                api_key: str | None = None,
                                logger: logging.Logger | None = None) -> dict[str, object]:
    """Count the host-independent PubMed corpus and fetch a stable sample per pathogen.

    Full PMID lists are cached independently from the sample cap. Increasing the
    cap reuses both the search and previously fetched article metadata.
    """
    if not email or "@" not in email:
        raise ValueError("Set a valid NCBI contact email in NCBI_EMAIL before acquisition.")
    logger = logger or logging.getLogger("pathogen_atlas_v4.pubmed.pathogen_literature")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    checkpoint, article_file = root / "pathogen_search_checkpoints.jsonl", root / "articles.jsonl"
    done = {}
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                record = json.loads(line)
                done[str(record["pathogen_id"])] = record
    article_cache = {}
    if article_file.exists():
        with article_file.open(encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    row = json.loads(line)
                    article_cache[str(row["pmid"])] = row
    pidx = pathogens.set_index("pathogen_id")
    aliases_by_pathogen = {}
    if not aliases.empty and {"pathogen_id", "name"}.issubset(aliases.columns):
        for pid, group in aliases.groupby(aliases.pathogen_id.astype(str)):
            aliases_by_pathogen[str(pid)] = group.name.dropna().astype(str).tolist()
    session = requests.Session()
    session.headers.update({"User-Agent": "PathogenAtlasV4/1.0 (contact: " + email + ")"})
    selected_ids = list(dict.fromkeys(map(str, pathogen_ids)))
    started, updated = time.monotonic(), 0
    with checkpoint.open("a", encoding="utf-8") as out, article_file.open("a", encoding="utf-8") as article_out:
        for index, pathogen_id in enumerate(selected_ids, 1):
            record = done.get(pathogen_id)
            if not record or record.get("status") not in {"searched", "completed"}:
                pathogen = pidx.loc[pathogen_id]
                query_record = build_pathogen_query([pathogen.accepted_name, *aliases_by_pathogen.get(pathogen_id, [])], strategy_version=query_version)
                pmids, segments = _all_pmids(session, query_record["query"], email=email, api_key=api_key, logger=logger)
                record = {"pathogen_id": pathogen_id, "pathogen_group": str(pathogen.pathogen_group), "query": query_record["query"],
                          "query_hash": query_record["query_hash"], "query_version": query_record["query_version"], "pmids": pmids,
                          "search_count": len(pmids), "segmented_search": bool(segments), "segments": segments, "status": "searched"}
                out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                out.flush()
                done[pathogen_id] = record
            sample = select_stable_sample(record["pmids"], maximum=sample_maximum, oldest_count=oldest_count, seed=seed)
            missing_pmids = [pmid for pmid in sample if pmid not in article_cache]
            same_settings = (record.get("sample_maximum") == sample_maximum
                             and record.get("sample_oldest_count") == oldest_count
                             and record.get("sample_seed") == seed
                             and record.get("sample_policy_version") == sample_policy_version)
            if record.get("status") == "completed" and same_settings and not missing_pmids:
                continue
            for row in _fetch_pmids(session, missing_pmids, email=email, api_key=api_key, logger=logger):
                if row["pmid"] not in article_cache:
                    article_cache[row["pmid"]] = row
                    article_out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            article_out.flush()
            record.update({"sample_maximum": sample_maximum, "sample_oldest_count": oldest_count,
                           "sample_seed": seed, "sample_policy_version": sample_policy_version,
                           "selected_sample_count": len(sample), "status": "completed"})
            out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            out.flush()
            done[pathogen_id] = record
            updated += 1
            if updated % 25 == 0 or index == len(selected_ids):
                elapsed = max(time.monotonic() - started, 0.001)
                pending = sum(1 for item in selected_ids[index:] if done.get(item, {}).get("status") != "completed")
                eta = pending / (updated / elapsed) / 60 if updated else float("inf")
                logger.info("Organism literature %d/%d (%.1f%%); updated=%d; ETA=%.1f min; articles cached=%d; sample cap=%s",
                            index, len(selected_ids), 100 * index / max(len(selected_ids), 1), updated, eta, len(article_cache), sample_maximum)
            time.sleep(0.11 if api_key else 0.36)
    count_rows, match_rows = [], []
    for pathogen_id in selected_ids:
        record = done.get(pathogen_id)
        if not record:
            continue
        sample = select_stable_sample(record["pmids"], maximum=sample_maximum, oldest_count=oldest_count, seed=seed)
        count_rows.append({"pathogen_id": pathogen_id, "pathogen_group": record["pathogen_group"], "query_hash": record["query_hash"],
                           "query_version": record["query_version"], "n_pubmed_search_records": record["search_count"],
                           "n_selected_for_review": len(sample), "sample_maximum": sample_maximum,
                           "sample_oldest_count": oldest_count, "sample_seed": seed,
                           "sample_policy_version": sample_policy_version, "status": record["status"]})
        for pmid in record["pmids"]:
            match_rows.append({"pathogen_id": pathogen_id, "pathogen_group": record["pathogen_group"], "pmid": str(pmid),
                               "query_hash": record["query_hash"], "query_version": record["query_version"],
                               "selected_for_sampling": str(pmid) in sample, "sample_selection": sample.get(str(pmid), "not_sampled"),
                               "metadata_available": str(pmid) in article_cache})
    tables = root / "tables"; tables.mkdir(exist_ok=True)
    pd.DataFrame(count_rows).sort_values("pathogen_id").to_parquet(tables / "pathogen_literature_counts.parquet", index=False)
    pd.DataFrame(match_rows, columns=["pathogen_id", "pathogen_group", "pmid", "query_hash", "query_version", "selected_for_sampling", "sample_selection", "metadata_available"]).sort_values(["pathogen_id", "pmid"]).to_parquet(tables / "pathogen_literature_matches.parquet", index=False)
    pd.DataFrame([article_cache[key] for key in sorted(article_cache, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))]).to_parquet(tables / "articles.parquet", index=False)
    return {"pathogens_total": len(selected_ids), "pathogens_complete": len(count_rows), "pathogens_updated_this_session": updated,
            "total_literature_records": sum(row["n_pubmed_search_records"] for row in count_rows),
            "sampled_pathogen_pmids": sum(row["n_selected_for_review"] for row in count_rows),
            "sample_cap_per_pathogen": sample_maximum, "oldest_per_pathogen": oldest_count, "sample_seed": seed,
            "sample_policy_version": sample_policy_version,
            "unique_article_metadata_cached": len(article_cache), "output_dir": str(tables), "checkpoint": str(checkpoint)}


def acquire_pubmed_pairs(pairs: Iterable[tuple[str, str]], pathogens: pd.DataFrame, hosts: pd.DataFrame,
                         aliases: pd.DataFrame, output_dir: str | Path, *, query_version: str,
                         email: str, api_key: str | None = None, logger: logging.Logger | None = None,
                         sample_maximum: int | None = None, sample_oldest_count: int = 30,
                         sample_seed: int = 42) -> dict[str, object]:
    """Search every pair and fetch every matching PMID, resuming by query and PMID cache."""
    if not email or "@" not in email:
        raise ValueError("Set a valid NCBI contact email in NCBI_EMAIL before acquisition.")
    logger = logger or logging.getLogger("pathogen_atlas_v4.pubmed")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    checkpoint = root / "pair_search_checkpoints.jsonl"
    article_file = root / "articles.jsonl"
    done = {}
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                record = json.loads(line)
                done[(record["pathogen_id"], record["host_id"])] = record
    article_cache = {}
    if article_file.exists():
        with article_file.open(encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    row = json.loads(line)
                    article_cache[row["pmid"]] = row
    pidx, hidx = pathogens.set_index("pathogen_id"), hosts.set_index("host_id")
    aliases_by_pathogen = {}
    if not aliases.empty and {"pathogen_id", "name"}.issubset(aliases.columns):
        for pid, group in aliases.groupby(aliases.pathogen_id.astype(str)):
            aliases_by_pathogen[str(pid)] = group.name.dropna().astype(str).tolist()
    session = requests.Session()
    session.headers.update({"User-Agent": "PathogenAtlasV4/1.0 (contact: " + email + ")"})
    selected_pairs = list(pairs)
    total = len(selected_pairs)
    started = time.monotonic()
    completed_this_run = 0
    with checkpoint.open("a", encoding="utf-8") as out, article_file.open("a", encoding="utf-8") as article_out:
        for index, (pid, hid) in enumerate(selected_pairs, 1):
            record = done.get((pid, hid))
            if not record or record.get("status") not in {"searched", "completed"}:
                pathogen, host = pidx.loc[pid], hidx.loc[hid]
                pnames = [pathogen.accepted_name, *aliases_by_pathogen.get(pid, [])]
                haliases = host.get("aliases", [])
                if hasattr(haliases, "tolist"):
                    haliases = haliases.tolist()
                query_record = build_pair_query(pnames, [host.accepted_name, *(haliases or [])], strategy_version=query_version)
                pmids, segments = _all_pmids(session, query_record["query"], email=email, api_key=api_key, logger=logger)
                record = {"pathogen_id": str(pid), "host_id": str(hid), "relationship_context": "disease", "query": query_record["query"],
                          "query_hash": query_record["query_hash"], "query_version": query_record["query_version"], "pmids": pmids,
                          "search_count": len(pmids), "segmented_search": bool(segments), "segments": segments, "status": "searched"}
                out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                out.flush()
                done[(pid, hid)] = record
            sample = select_stable_sample(record["pmids"], maximum=sample_maximum,
                                          oldest_count=sample_oldest_count, seed=sample_seed)
            missing_pmids = [pmid for pmid in sample if pmid not in article_cache]
            settings_match = (record.get("sample_maximum") == sample_maximum
                              and record.get("sample_oldest_count") == sample_oldest_count
                              and record.get("sample_seed") == sample_seed)
            if record.get("status") == "completed" and settings_match and not missing_pmids:
                continue
            for row in _fetch_pmids(session, missing_pmids, email=email, api_key=api_key, logger=logger):
                if row["pmid"] not in article_cache:
                    article_cache[row["pmid"]] = row
                    article_out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            article_out.flush()
            record.update({"sample_maximum": sample_maximum, "sample_oldest_count": sample_oldest_count,
                           "sample_seed": sample_seed, "selected_sample_count": len(sample), "status": "completed"})
            out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            out.flush()
            done[(pid, hid)] = record
            completed_this_run += 1
            if completed_this_run % 25 == 0 or index == total:
                elapsed = max(time.monotonic() - started, 0.001)
                rate = completed_this_run / elapsed
                remaining = sum(1 for pair in selected_pairs[index:] if done.get(pair, {}).get("status") != "completed")
                eta_minutes = remaining / rate / 60 if rate else float("inf")
                logger.info("PubMed pair checkpoints %d/%d (%.1f%%); updated this session=%d; ETA %.1f min; cached article PMIDs=%d; sample cap=%s",
                            index, total, index * 100 / max(total, 1), completed_this_run, eta_minutes, len(article_cache), sample_maximum)
            time.sleep(0.11 if api_key else 0.36)
    pair_rows = []
    for pair in selected_pairs:
        record = done.get(pair)
        if record:
            pair_rows.append({k: record.get(k) for k in ("pathogen_id", "host_id", "relationship_context", "query", "query_hash", "query_version", "search_count", "segmented_search", "sample_maximum", "sample_oldest_count", "sample_seed", "selected_sample_count", "status")})
    match_rows = []
    for pair in selected_pairs:
        record = done.get(pair)
        if record:
            current_sample = select_stable_sample(record["pmids"], maximum=sample_maximum,
                                                  oldest_count=sample_oldest_count, seed=sample_seed)
            for pmid in record["pmids"]:
                match_rows.append({"pathogen_id": pair[0], "host_id": pair[1], "relationship_context": "disease", "pmid": str(pmid),
                                   "query_hash": record["query_hash"], "query_version": record["query_version"],
                                   "selected_for_sampling": str(pmid) in current_sample,
                                   "selected_for_fetch": str(pmid) in current_sample,
                                   "sample_selection": current_sample.get(str(pmid), "not_sampled"),
                                   "metadata_available": str(pmid) in article_cache})
    out_root = root / "tables"
    out_root.mkdir(exist_ok=True)
    pd.DataFrame(pair_rows).sort_values(["pathogen_id", "host_id"]).to_parquet(out_root / "pair_searches.parquet", index=False)
    pd.DataFrame(match_rows, columns=["pathogen_id", "host_id", "relationship_context", "pmid", "query_hash", "query_version", "selected_for_sampling", "selected_for_fetch", "sample_selection", "metadata_available"]).sort_values(["pathogen_id", "host_id", "pmid"]).to_parquet(out_root / "pathogen_host_article_matches.parquet", index=False)
    pd.DataFrame([article_cache[key] for key in sorted(article_cache, key=lambda x: int(x))]).to_parquet(out_root / "articles.parquet", index=False)
    return {"pairs_total": total, "pairs_complete": len(pair_rows), "pairs_updated_this_session": completed_this_run,
            "pair_matches": len(match_rows), "sample_selected_matches": sum(bool(row["selected_for_sampling"]) for row in match_rows),
            "sample_cap_per_pair": sample_maximum, "oldest_per_pair": sample_oldest_count, "sample_seed": sample_seed,
            "unique_sample_articles_cached": len(article_cache), "output_dir": str(out_root), "checkpoint": str(checkpoint)}
