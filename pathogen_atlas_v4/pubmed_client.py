"""Resumable, pair-scoped PubMed acquisition through NCBI E-utilities."""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from pathogen_atlas_v4.pubmed import build_pair_query, build_pathogen_query

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200
_worker_local = threading.local()


def literature_cache_key(*, query_version: str, pathogen_literature_query_version: str,
                         taxonomy_snapshot: str, host_expansion_hash: str = "") -> str:
    """Stable shared literature-cache identity used by acquisition and downstream notebooks."""
    identity = {"query_version": str(query_version),
                "pathogen_literature_query_version": str(pathogen_literature_query_version),
                "taxonomy_snapshot": str(taxonomy_snapshot)}
    # Empty expansion hashes preserve the cache identity used before host
    # expansions were added to the key. Non-empty hashes isolate expanded runs.
    if host_expansion_hash:
        identity["host_expansion_hash"] = str(host_expansion_hash)
    return hashlib.sha256(json.dumps(identity, sort_keys=True, default=list).encode()).hexdigest()


def ncbi_api_key() -> str | None:
    """Resolve NCBI credentials from environment or the local macOS Keychain."""
    value = (os.environ.get("NCBI_API_KEY") or os.environ.get("EUTILS_API_KEY") or "").strip()
    if value:
        return value
    if sys.platform != "darwin" or shutil.which("security") is None:
        return None
    try:
        result = subprocess.run(
            ["security", "find-generic-password", "-s", "NCBI", "-a", "NCBI_API_KEY", "-w"],
            capture_output=True, text=True, check=False, timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    key = result.stdout.strip()
    return key if result.returncode == 0 and key else None


class _RateLimiter:
    """Pace shared request starts below NCBI's keyed and unkeyed ceilings."""
    def __init__(self, api_key: str | None):
        self.interval = 0.125 if api_key else 0.34
        self.lock = threading.Lock()
        self.next_time = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            slot = max(now, self.next_time)
            self.next_time = slot + self.interval
        delay = slot - now
        if delay > 0:
            time.sleep(delay)


def _format_eta(seconds: float) -> str:
    """Format a remaining duration with days, hours, and minutes as needed."""
    if not math.isfinite(seconds):
        return "unknown"
    total_minutes = max(0, math.ceil(seconds / 60))
    if total_minutes == 0:
        return "0m"
    days, remainder = divmod(total_minutes, 24 * 60)
    hours, minutes = divmod(remainder, 60)
    parts = []
    if days:
        parts.append(f"{days}d")
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    return " ".join(parts)


def pubmed_aliases_by_pathogen(pathogens: pd.DataFrame, aliases: pd.DataFrame
                               ) -> tuple[dict[str, list[str]], dict[str, int]]:
    """Return unique, taxon-specific names suitable for PubMed queries.

    Alongside accepted names and nomenclatural synonyms, extract ICTV VMR's
    explicit virus-name and virus-name-abbreviation fields. VMR isolate names
    and strain labels are not used. Every candidate must map to one pathogen;
    homonyms shared across taxa are excluded. The accepted taxon name itself
    is added separately by each query builder.
    """
    counts = {
        "aliases_total": int(len(aliases)),
        "aliases_missing_name": 0,
        "aliases_unmapped": 0,
        "aliases_non_nomenclatural": 0,
        "aliases_shared_across_pathogens": 0,
        "aliases_vmr_records_checked": 0,
        "aliases_vmr_terms_extracted": 0,
        "aliases_searchable": 0,
    }
    if aliases.empty or not {"pathogen_id", "name"}.issubset(aliases.columns):
        return {}, counts

    rows = aliases.copy()
    rows["_pid"] = rows["pathogen_id"].astype("string").str.strip()
    rows["_name"] = rows["name"].astype("string").str.strip()
    rows["_normalized"] = rows["_name"].str.casefold().str.replace(r"\s+", " ", regex=True)
    missing_name = rows["_name"].isna() | rows["_name"].eq("")
    pathogen_ids = set(pathogens["pathogen_id"].astype(str)) if "pathogen_id" in pathogens else set()
    unmapped = rows["_pid"].isna() | ~rows["_pid"].isin(pathogen_ids)
    name_types = rows.get("name_type", pd.Series("", index=rows.index)).astype("string").str.casefold().str.strip()
    mapping_status = rows.get("mapping_status", pd.Series("", index=rows.index)).astype("string").str.casefold().str.strip()
    sources = rows.get("source", pd.Series("", index=rows.index)).astype("string").str.casefold().str.strip()
    nomenclatural = name_types.isin({"accepted", "synonym"}) & mapping_status.isin({"accepted", "synonym"})
    valid_nomenclatural = ~(missing_name | unmapped | ~nomenclatural)
    names_by_pathogen: dict[str, set[str]] = {}
    accepted_by_pathogen: dict[str, set[str]] = {}
    taxonomy_ids_by_pathogen: dict[str, set[str]] = {}
    if {"pathogen_id", "accepted_name"}.issubset(pathogens.columns):
        for row in pathogens.to_dict("records"):
            pid = str(row["pathogen_id"])
            accepted_value = row.get("accepted_name")
            name = " ".join(str(accepted_value).split()).casefold() if pd.notna(accepted_value) else ""
            if name:
                accepted_by_pathogen.setdefault(pid, set()).add(name)
                names_by_pathogen.setdefault(pid, set()).add(name)
            taxon_value = row.get("taxon_id")
            taxon_id = str(taxon_value).strip() if pd.notna(taxon_value) else ""
            if taxon_id:
                taxonomy_ids_by_pathogen.setdefault(pid, set()).add(taxon_id.casefold())

    candidate_aliases: list[tuple[str, str]] = []
    for pid, name in rows.loc[valid_nomenclatural, ["_pid", "_name"]].itertuples(index=False, name=None):
        candidate_aliases.append((str(pid), " ".join(str(name).split())))

    def split_vmr_names(value: object) -> list[str]:
        if value is None or pd.isna(value):
            return []
        return [" ".join(part.split()) for part in re.split(r"[;\n]+", str(value)) if len(part.strip()) >= 3]

    safe_vmr_rows = pd.Series(False, index=rows.index)
    if "source_record_json" in rows:
        for index, row in rows.loc[(sources == "ictv_vmr") & ~unmapped & ~missing_name].iterrows():
            raw_record = row.get("source_record_json")
            try:
                source_record = raw_record if isinstance(raw_record, dict) else json.loads(str(raw_record))
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if not isinstance(source_record, dict):
                continue
            pid = str(row["_pid"])
            source_species = " ".join(str(source_record.get("Species") or "").split()).casefold()
            source_ictv_id = str(source_record.get("ICTV_ID") or "").strip().casefold()
            if (source_species not in accepted_by_pathogen.get(pid, set())
                    and source_ictv_id not in taxonomy_ids_by_pathogen.get(pid, set())):
                continue
            safe_vmr_rows.loc[index] = True
            # VMR has distinct fields for actual virus names and abbreviations.
            # Do not search isolate designations or raw source-alias strings.
            for field in ("Virus name(s)", "Virus name abbreviation(s)"):
                terms = split_vmr_names(source_record.get(field))
                candidate_aliases.extend((pid, value) for value in terms)
                counts["aliases_vmr_terms_extracted"] += len(terms)
            counts["aliases_vmr_records_checked"] += 1

    non_nomenclatural = ~nomenclatural & ~safe_vmr_rows
    counts["aliases_missing_name"] = int(missing_name.sum())
    counts["aliases_unmapped"] = int((~missing_name & unmapped).sum())
    counts["aliases_non_nomenclatural"] = int((~missing_name & ~unmapped & non_nomenclatural).sum())

    pathogen_ids_by_name: dict[str, set[str]] = {}
    for pid, names in names_by_pathogen.items():
        for name in names:
            pathogen_ids_by_name.setdefault(name, set()).add(pid)
    normalized_candidates = []
    for pid, display_name in candidate_aliases:
        name = " ".join(display_name.split()).casefold()
        if len(name) < 3:
            continue
        normalized_candidates.append((pid, name, display_name))
        pathogen_ids_by_name.setdefault(name, set()).add(pid)
    shared_names = {name for name, ids in pathogen_ids_by_name.items() if len(ids) > 1}
    searchable: dict[str, list[str]] = {}
    seen_by_pathogen: dict[str, set[str]] = {pid: set(names) for pid, names in accepted_by_pathogen.items()}
    counts["aliases_shared_across_pathogens"] = sum(name in shared_names for _, name, _ in normalized_candidates)
    for pid, name, display_name in normalized_candidates:
        if name in shared_names or name in seen_by_pathogen.setdefault(pid, set()):
            continue
        searchable.setdefault(pid, []).append(display_name)
        seen_by_pathogen[pid].add(name)
    counts["aliases_searchable"] = sum(map(len, searchable.values()))
    return searchable, counts


def _thread_session(email: str) -> requests.Session:
    if not hasattr(_worker_local, "session"):
        _worker_local.session = requests.Session()
        _worker_local.session.headers.update({"User-Agent": f"PathogenAtlasV4/1.0 (contact: {email})"})
    return _worker_local.session


def _xml(response: requests.Response) -> ET.Element:
    response.raise_for_status()
    return ET.fromstring(response.content)


def _search(session: requests.Session, query: str, *, email: str, api_key: str | None, logger: logging.Logger, limiter=None):
    """Return a complete (at most 10,000) PubMed UID page and result count."""
    params = {"db": "pubmed", "term": query, "usehistory": "y", "retmax": 10000,
              "email": email, "tool": "pathogen_atlas_v4"}
    if api_key:
        params["api_key"] = api_key
    for attempt in range(7):
        try:
            if limiter:
                limiter.wait()
            root = _xml(session.post(f"{EUTILS}/esearch.fcgi", data=params, timeout=90))
            error = root.findtext("Error") or root.findtext("ERROR")
            if error:
                raise RuntimeError(f"NCBI ESearch error: {error}")
            count = int(root.findtext("Count", "0"))
            pmids = [node.text.strip() for node in root.findall(".//IdList/Id") if node.text]
            if not limiter:
                time.sleep(0.11 if api_key else 0.36)
            return count, pmids
        except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
            if attempt == 6:
                raise
            delay = min(2 ** attempt, 60)
            logger.warning("ESearch retry %d/6 after %s; sleeping %ds", attempt + 1, type(exc).__name__, delay)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _fetch_pmids(session: requests.Session, pmids: list[str], *, email: str, api_key: str | None,
                 logger: logging.Logger, include_abstract: bool = False):
    """Fetch journal and book records, recovering partial EFetch batches by PMID."""
    articles = []
    params = {"db": "pubmed", "retmode": "xml", "email": email, "tool": "pathogen_atlas_v4"}
    if api_key:
        params["api_key"] = api_key

    def fetch(batch: list[str], offset: int, attempts: int = 7) -> list[dict[str, str]]:
        params["id"] = ",".join(batch)
        for attempt in range(attempts):
            try:
                root = _xml(session.post(f"{EUTILS}/efetch.fcgi", data=params, timeout=120))
                error = root.findtext("Error") or root.findtext("ERROR")
                if error:
                    raise RuntimeError(f"NCBI EFetch error: {error}")
                return _parse_articles(root, include_abstract=include_abstract)
            except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
                if attempt == attempts - 1:
                    raise
                delay = min(2 ** attempt, 60)
                logger.warning("EFetch retry at sample offset %d after %s: %s; sleeping %ds",
                               offset, type(exc).__name__, exc, delay)
                time.sleep(delay)
        raise AssertionError("unreachable")

    for start in range(0, len(pmids), BATCH_SIZE):
        batch = [str(pmid) for pmid in pmids[start:start + BATCH_SIZE]]
        rows_by_pmid = {str(row["pmid"]): row for row in fetch(batch, start)}
        missing = [pmid for pmid in batch if pmid not in rows_by_pmid]
        if missing:
            logger.warning("EFetch batch returned %d of %d PMIDs; retrying %d unmatched PMIDs individually",
                           len(batch) - len(missing), len(batch), len(missing))
            for pmid in missing:
                time.sleep(0.11 if api_key else 0.36)
                try:
                    individual = fetch([pmid], start, attempts=3)
                except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
                    logger.warning("EFetch could not retrieve PMID %s after individual retries: %s", pmid, exc)
                    individual = []
                match = next((row for row in individual if str(row["pmid"]) == pmid), None)
                if match is not None:
                    rows_by_pmid[pmid] = match
                else:
                    logger.warning("EFetch returned no article for PMID %s; retaining its pair match with metadata unavailable", pmid)
                    rows_by_pmid[pmid] = {
                        "pmid": pmid, "title": "", "abstract": "",
                        "abstract_fetched": include_abstract,
                        "journal": "", "publication_year": "", "doi": "",
                        "metadata_status": "unavailable", "record_type": "unavailable",
                    }
        articles.extend(rows_by_pmid[pmid] for pmid in batch)
        time.sleep(0.11 if api_key else 0.36)
    return articles


def _parse_articles(root: ET.Element, *, include_abstract: bool = False) -> list[dict[str, str]]:
    """Parse journal articles and book articles returned by PubMed EFetch."""
    rows = []
    records = [(item, "journal_article") for item in root.findall(".//PubmedArticle")]
    records.extend((item, "book_article") for item in root.findall(".//PubmedBookArticle"))
    for article, record_type in records:
        pmid = (article.findtext(".//PMID", "") or "").strip()
        if not pmid:
            continue
        medline = article.find(".//MedlineCitation")
        art = medline.find(".//Article") if medline is not None else None
        book = article.find(".//Book") if record_type == "book_article" else None
        title_node = art.find("ArticleTitle") if art is not None else None
        title = "".join(title_node.itertext()).strip() if title_node is not None else ""
        if not title and book is not None:
            title = (book.findtext("BookTitle", "") or "").strip()
        abstract_nodes = article.findall(".//Abstract/AbstractText")
        abstract = " ".join("".join(node.itertext()).strip() for node in abstract_nodes
                            if "".join(node.itertext()).strip()) if include_abstract else ""
        journal = art.findtext(".//Journal/Title", "") if art is not None else ""
        if not journal and book is not None:
            journal = book.findtext("CollectionTitle", "") or ""
        year = ""
        if art is not None:
            year = art.findtext(".//JournalIssue/PubDate/Year", "") or art.findtext(".//JournalIssue/PubDate/MedlineDate", "")
        if not year and book is not None:
            year = book.findtext(".//PubDate/Year", "") or ""
        doi = ""
        if art is not None:
            for ident in art.findall(".//ELocationID"):
                if ident.attrib.get("EIdType", "").casefold() == "doi":
                    doi = "".join(ident.itertext()).strip()
                    break
        rows.append({"pmid": pmid, "title": title, "abstract": abstract, "abstract_fetched": include_abstract,
                     "journal": journal, "publication_year": year, "doi": doi,
                     "metadata_status": "available", "record_type": record_type})
    return rows

def _date_range_query(query: str, start: str, end: str) -> str:
    return f'({query}) AND ("{start}"[DP] : "{end}"[DP])'


def _all_pmids(session, query: str, *, email, api_key, logger, limiter=None):
    """Collect all UIDs, segmenting large result sets before ESearch's 10,000 UID limit."""
    count, pmids = _search(session, query, email=email, api_key=api_key, logger=logger, limiter=limiter)
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
        year_count, year_pmids = _search(session, year_query, email=email, api_key=api_key, logger=logger, limiter=limiter)
        if year_count <= 10000:
            if len(year_pmids) != year_count:
                raise RuntimeError(f"ESearch returned {len(year_pmids)} of {year_count} UIDs for {year}")
            segmented_pmids.extend(year_pmids)
            segment_counts.append({"period": str(year), "count": year_count})
            continue
        for month in range(1, 13):
            last_day = calendar.monthrange(year, month)[1]
            month_query = _date_range_query(query, f"{year}/{month:02d}/01", f"{year}/{month:02d}/{last_day:02d}")
            month_count, month_pmids = _search(session, month_query, email=email, api_key=api_key, logger=logger, limiter=limiter)
            if month_count <= 10000:
                if len(month_pmids) != month_count:
                    raise RuntimeError(f"ESearch returned {len(month_pmids)} of {month_count} UIDs for {year}-{month:02d}")
                segmented_pmids.extend(month_pmids)
                segment_counts.append({"period": f"{year}-{month:02d}", "count": month_count})
                continue
            for day in range(1, last_day + 1):
                daily_query = _date_range_query(query, f"{year}/{month:02d}/{day:02d}", f"{year}/{month:02d}/{day:02d}")
                daily_count, daily_pmids = _search(session, daily_query, email=email, api_key=api_key, logger=logger, limiter=limiter)
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


def prepare_pair_screen_samples(pairs: Iterable[tuple[str, str]], output_dir: str | Path, *,
                                sample_maximum: int | None = 200, oldest_count: int = 30,
                                seed: int = 42, sample_policy_version: str = "stable-prefix-hash-v1",
                                email: str, api_key: str | None = None,
                                logger: logging.Logger | None = None) -> dict[str, object]:
    """Sample pair-search results and fetch abstracts needed for the yes/no screen.

    PMID search lists are reused from the complete pair-search checkpoint. The
    deterministic sample extends when the cap grows; abstract EFetches are
    limited to newly selected PMIDs missing from the shared cache.
    """
    if not email or "@" not in email:
        raise ValueError("Set a valid NCBI contact email in NCBI_EMAIL before acquisition.")
    api_key = (api_key or ncbi_api_key() or "").strip() or None
    logger = logger or logging.getLogger("pathogen_atlas_v4.pubmed.pair_screen")
    root = Path(output_dir)
    checkpoint, article_file = root / "pair_search_checkpoints.jsonl", root / "articles.jsonl"
    if not checkpoint.exists():
        raise FileNotFoundError(f"Run Phase 2 pair acquisition first: {checkpoint}")
    records = {}
    for line in checkpoint.read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            records[(str(record["pathogen_id"]), str(record["host_id"]))] = record
    article_cache = {}
    if article_file.exists():
        with article_file.open(encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    row = json.loads(line)
                    row.setdefault("metadata_status", "available")
                    row.setdefault("record_type", "journal_article")
                    article_cache[str(row["pmid"])] = row
    selected_pairs = list(dict.fromkeys((str(pathogen), str(host)) for pathogen, host in pairs))
    missing_pairs = [pair for pair in selected_pairs if pair not in records or records[pair].get("status") != "completed"]
    if missing_pairs:
        raise RuntimeError(f"Phase 2 is incomplete for {len(missing_pairs)} selected pathogen–host pairs")
    session = requests.Session()
    session.headers.update({"User-Agent": "PathogenAtlasV4/1.0 (contact: " + email + ")"})
    sample_rows, missing_abstracts = [], set()
    for pathogen_id, host_id in selected_pairs:
        record = records[(pathogen_id, host_id)]
        sample = select_stable_sample(record["pmids"], maximum=sample_maximum, oldest_count=oldest_count, seed=seed)
        for pmid, selection in sample.items():
            cached = article_cache.get(str(pmid), {})
            if (not bool(cached.get("abstract_fetched", False))
                    or cached.get("metadata_status") == "unavailable"):
                missing_abstracts.add(str(pmid))
            sample_rows.append({"pathogen_id": pathogen_id, "host_id": host_id, "relationship_context": "disease",
                                "pmid": str(pmid), "query_hash": record["query_hash"], "query_version": record["query_version"],
                                "pair_search_count": int(record.get("search_count", len(record["pmids"]))),
                                "sample_selection": selection, "abstract_fetched": bool(cached.get("abstract_fetched", False)),
                                "abstract_available": bool(str(cached.get("abstract") or "").strip())})
    updated = 0
    with article_file.open("a", encoding="utf-8") as article_out:
        missing_list = sorted(missing_abstracts, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))
        for row in _fetch_pmids(session, missing_list, email=email, api_key=api_key, logger=logger, include_abstract=True):
            article_cache[row["pmid"]] = row
            article_out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            updated += 1
        article_out.flush()
    for row in sample_rows:
        cached = article_cache.get(row["pmid"], {})
        row["abstract_fetched"] = bool(cached.get("abstract_fetched", False))
        row["abstract_available"] = bool(str(cached.get("abstract") or "").strip())
    tables = root / "tables"; tables.mkdir(exist_ok=True)
    matches = pd.DataFrame(sample_rows, columns=["pathogen_id", "host_id", "relationship_context", "pmid", "query_hash", "query_version", "pair_search_count", "sample_selection", "abstract_fetched", "abstract_available"])
    matches.sort_values(["pathogen_id", "host_id", "pmid"]).to_parquet(tables / "pathogen_screen_matches.parquet", index=False)
    article_columns = ["pmid", "title", "abstract", "abstract_fetched", "journal", "publication_year", "doi",
                       "metadata_status", "record_type"]
    pd.DataFrame([article_cache[key] for key in sorted(article_cache, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))],
                 columns=article_columns).to_parquet(tables / "articles.parquet", index=False)
    manifest = matches.groupby(["pathogen_id", "host_id"], as_index=False).agg(sample_count=("pmid", "nunique"), abstracts_available=("abstract_available", "sum"))
    manifest["sample_maximum"] = sample_maximum
    manifest["oldest_count"] = oldest_count
    manifest["sample_seed"] = seed
    manifest["sample_policy_version"] = sample_policy_version
    manifest.to_parquet(tables / "pathogen_screen_sample_manifest.parquet", index=False)
    logger.info("Prepared pathogen screen samples: pairs=%d sample_rows=%d new_abstract_fetches=%d cap=%s",
                len(selected_pairs), len(matches), updated, sample_maximum)
    return {"pairs": len(selected_pairs), "sample_rows": len(matches), "new_abstract_fetches": updated,
            "sample_cap_per_pair": sample_maximum, "oldest_per_pair": oldest_count,
            "sample_seed": seed, "sample_policy_version": sample_policy_version,
            "outputs": str(tables)}


def _count_search(session: requests.Session, query: str, *, email: str, api_key: str | None,
                  logger: logging.Logger) -> int:
    """Count PubMed matches without returning PMIDs; called only for selected pathogens."""
    params = {"db": "pubmed", "term": query, "usehistory": "y", "retmax": 0,
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
            time.sleep(0.11 if api_key else 0.36)
            return count
        except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
            if attempt == 6:
                raise
            delay = min(2 ** attempt, 60)
            logger.warning("PubMed count retry %d/6 after %s; sleeping %ds", attempt + 1, type(exc).__name__, delay)
            time.sleep(delay)
    raise AssertionError("unreachable")


def count_accepted_pathogen_publications(accepted_pathogen_ids: Iterable[str], pathogens: pd.DataFrame,
                                         aliases: pd.DataFrame, output_dir: str | Path, *,
                                         query_version: str, email: str, api_key: str | None = None,
                                         logger: logging.Logger | None = None) -> pd.DataFrame:
    """Count host-independent PubMed records only for pathogen IDs accepted in a host."""
    if not email or "@" not in email:
        raise ValueError("Set a valid NCBI contact email in NCBI_EMAIL before acquisition.")
    api_key = (api_key or ncbi_api_key() or "").strip() or None
    accepted_ids = sorted(set(map(str, accepted_pathogen_ids)))
    known_ids = set(pathogens.pathogen_id.astype(str))
    unknown = sorted(set(accepted_ids) - known_ids)
    if unknown:
        raise ValueError(f"Accepted pathogen IDs are not in the resolved taxonomy: {unknown[:10]}")
    logger = logger or logging.getLogger("pathogen_atlas_v4.pubmed.accepted_counts")
    root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
    checkpoint = root / "accepted_pathogen_publication_counts.jsonl"
    done = {}
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                record = json.loads(line)
                done[str(record["pathogen_id"])] = record
    aliases_by_pathogen, alias_counts = pubmed_aliases_by_pathogen(pathogens, aliases)
    logger.info("PubMed alias safeguard: %s", json.dumps(alias_counts, sort_keys=True))
    pidx = pathogens.set_index("pathogen_id")
    session = requests.Session()
    session.headers.update({"User-Agent": "PathogenAtlasV4/1.0 (contact: " + email + ")"})
    with checkpoint.open("a", encoding="utf-8") as out:
        for index, pathogen_id in enumerate(accepted_ids, 1):
            pathogen = pidx.loc[pathogen_id]
            query_record = build_pathogen_query([pathogen.accepted_name, *aliases_by_pathogen.get(pathogen_id, [])], strategy_version=query_version)
            cached = done.get(pathogen_id)
            if cached and cached.get("query_hash") == query_record["query_hash"] and cached.get("status") == "completed":
                continue
            count = _count_search(session, query_record["query"], email=email, api_key=api_key, logger=logger)
            record = {"pathogen_id": pathogen_id, "pathogen_group": str(pathogen.pathogen_group),
                      "query_hash": query_record["query_hash"], "query_version": query_version,
                      "n_pubmed_search_records": count, "status": "completed"}
            out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            out.flush(); done[pathogen_id] = record
            if index % 50 == 0 or index == len(accepted_ids):
                logger.info("Accepted pathogen publication counts %d/%d", index, len(accepted_ids))
    rows = [done[pathogen_id] for pathogen_id in accepted_ids if pathogen_id in done]
    result = pd.DataFrame(rows, columns=["pathogen_id", "pathogen_group", "query_hash", "query_version", "n_pubmed_search_records", "status"])
    tables = root / "tables"; tables.mkdir(exist_ok=True)
    result.sort_values("pathogen_id").to_parquet(tables / "accepted_pathogen_publication_counts.parquet", index=False)
    return result


def acquire_pubmed_pairs(pairs: Iterable[tuple[str, str]], pathogens: pd.DataFrame, hosts: pd.DataFrame,
                         aliases: pd.DataFrame, output_dir: str | Path, *, query_version: str,
                         host_search_expansions: dict[str, dict[str, object]] | None = None,
                         email: str, api_key: str | None = None, logger: logging.Logger | None = None,
                         sample_maximum: int | None = 200, sample_oldest_count: int = 30,
                         sample_seed: int = 42, max_workers: int = 8) -> dict[str, object]:
    """Search every pair concurrently, retaining all PMIDs and resuming checkpoints.

    Article metadata is fetched later by the pair-screening sample workflow.
    """
    if not email or "@" not in email:
        raise ValueError("Set a valid NCBI contact email in NCBI_EMAIL before acquisition.")
    api_key = (api_key or ncbi_api_key() or "").strip() or None
    if max_workers < 1:
        raise ValueError("max_workers must be positive.")
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
                    row.setdefault("metadata_status", "available")
                    row.setdefault("record_type", "journal_article")
                    article_cache[row["pmid"]] = row
    pidx, hidx = pathogens.set_index("pathogen_id"), hosts.set_index("host_id")
    aliases_by_pathogen, alias_counts = pubmed_aliases_by_pathogen(pathogens, aliases)
    host_search_expansions = {str(host_id): record for host_id, record in (host_search_expansions or {}).items()}
    logger.info("PubMed pair-search alias safeguard: %s", json.dumps(alias_counts, sort_keys=True))
    selected_pairs = list(pairs)
    total = len(selected_pairs)
    started = time.monotonic()
    completed_this_run = 0
    limiter = _RateLimiter(api_key)

    def pair_query(pid, hid):
        pathogen, host = pidx.loc[pid], hidx.loc[hid]
        pnames = [pathogen.accepted_name, *aliases_by_pathogen.get(pid, [])]
        haliases = host.get("aliases", [])
        if hasattr(haliases, "tolist"):
            haliases = haliases.tolist()
        expansion = host_search_expansions.get(str(hid), {})
        text_terms = [str(item.get("term", "")) for item in expansion.get("terms", [])
                      if isinstance(item, dict) and item.get("term")]
        mesh_terms = [str(term) for term in expansion.get("mesh_terms", []) if str(term).strip()]
        query_record = build_pair_query(
            pnames, [host.accepted_name, *(haliases or [])], host_mention_terms=text_terms,
            host_mesh_terms=mesh_terms, strategy_version=query_version)
        return query_record

    def search_pair(pid, hid):
        pathogen = pidx.loc[pid]
        query_record = pair_query(pid, hid)
        pmids, segments = _all_pmids(_thread_session(email), query_record["query"], email=email,
                                     api_key=api_key, logger=logger, limiter=limiter)
        record = {"pathogen_id": str(pid), "host_id": str(hid), "relationship_context": "disease", "query": query_record["query"],
                  "query_hash": query_record["query_hash"], "query_version": query_record["query_version"], "pmids": pmids,
                  "search_count": len(pmids), "segmented_search": bool(segments), "segments": segments, "status": "searched"}
        return pid, hid, record

    with checkpoint.open("a", encoding="utf-8") as out:
        expected_hashes = {pair: pair_query(*pair)["query_hash"] for pair in selected_pairs}
        pending = [(pid, hid) for pid, hid in selected_pairs
                   if done.get((pid, hid), {}).get("status") not in {"searched", "completed"}
                   or done.get((pid, hid), {}).get("query_hash") != expected_hashes[(pid, hid)]]
        already_complete = sum(1 for pair in selected_pairs
                               if done.get(pair, {}).get("status") == "completed"
                               and done.get(pair, {}).get("query_hash") == expected_hashes[pair])
        logger.info("PubMed pair acquisition settings: pairs=%d already_checkpointed=%d searches_to_run=%d; sample cap=%s oldest=%d seed=%d; workers=%d; API key=%s; cached article PMIDs=%d",
                    total, already_complete, len(pending), sample_maximum, sample_oldest_count,
                    sample_seed, max_workers, bool(api_key), len(article_cache))
        for pair in selected_pairs:
            record = done.get(pair)
            if (record and record.get("status") == "searched"
                    and record.get("query_hash") == expected_hashes[pair]):
                record.update({"sample_maximum": sample_maximum, "sample_oldest_count": sample_oldest_count,
                               "sample_seed": sample_seed,
                               "selected_sample_count": len(select_stable_sample(record["pmids"], maximum=sample_maximum,
                                                                                  oldest_count=sample_oldest_count, seed=sample_seed)),
                               "status": "completed"})
                out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                done[pair] = record
        out.flush()
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {pool.submit(search_pair, pid, hid): (pid, hid) for pid, hid in pending}
            failures = []
            for future in as_completed(futures):
                try:
                    pid, hid, record = future.result()
                except Exception as exc:
                    pid, hid = futures[future]
                    failures.append((pid, hid, exc))
                    logger.error("PubMed pair search failed pathogen=%s host=%s: %s", pid, hid, exc)
                    continue
                record.update({"sample_maximum": sample_maximum, "sample_oldest_count": sample_oldest_count,
                               "sample_seed": sample_seed,
                               "selected_sample_count": len(select_stable_sample(record["pmids"], maximum=sample_maximum,
                                                                                  oldest_count=sample_oldest_count, seed=sample_seed)),
                               "status": "completed"})
                out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                out.flush()
                done[(pid, hid)] = record
                completed_this_run += 1
                if completed_this_run % 25 == 0 or completed_this_run == len(pending):
                    out.flush()
                    finished = sum(1 for pair in selected_pairs if done.get(pair, {}).get("status") == "completed")
                    elapsed = max(time.monotonic() - started, 0.001)
                    rate = completed_this_run / elapsed
                    remaining = max(0, total - finished)
                    eta = _format_eta(remaining / rate if rate else float("inf"))
                    logger.info("PubMed pair checkpoints %d/%d (%.1f%%); updated this session=%d; ETA %s",
                                finished, total, finished * 100 / max(total, 1), completed_this_run, eta)
            if failures:
                examples = ", ".join(f"{pid}/{hid}" for pid, hid, _ in failures[:5])
                raise RuntimeError(f"{len(failures)} PubMed pair searches failed and were left resumable; examples: {examples}") from failures[0][2]
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
                                   "metadata_available": (str(pmid) in article_cache
                                                          and article_cache[str(pmid)].get("metadata_status", "available") == "available")})
    out_root = root / "tables"
    out_root.mkdir(exist_ok=True)
    pd.DataFrame(pair_rows).sort_values(["pathogen_id", "host_id"]).to_parquet(out_root / "pair_searches.parquet", index=False)
    pd.DataFrame(match_rows, columns=["pathogen_id", "host_id", "relationship_context", "pmid", "query_hash", "query_version", "selected_for_sampling", "selected_for_fetch", "sample_selection", "metadata_available"]).sort_values(["pathogen_id", "host_id", "pmid"]).to_parquet(out_root / "pathogen_host_article_matches.parquet", index=False)
    article_columns = ["pmid", "title", "abstract", "abstract_fetched", "journal", "publication_year", "doi",
                       "metadata_status", "record_type"]
    pd.DataFrame([article_cache[key] for key in sorted(article_cache, key=lambda x: int(x))],
                 columns=article_columns).to_parquet(out_root / "articles.parquet", index=False)
    return {"pairs_total": total, "pairs_complete": len(pair_rows), "pairs_updated_this_session": completed_this_run,
            "pair_matches": len(match_rows), "sample_selected_matches": sum(bool(row["selected_for_sampling"]) for row in match_rows),
            "sample_cap_per_pair": sample_maximum, "oldest_per_pair": sample_oldest_count, "sample_seed": sample_seed,
            "unique_sample_articles_cached": len(article_cache), "output_dir": str(out_root), "checkpoint": str(checkpoint)}
