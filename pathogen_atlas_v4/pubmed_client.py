"""Resumable, pair-scoped PubMed acquisition through NCBI E-utilities."""
from __future__ import annotations

import json
import logging
import os
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from pathogen_atlas_v4.pubmed import build_pair_query

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200


def _xml(response: requests.Response) -> ET.Element:
    response.raise_for_status()
    return ET.fromstring(response.content)


def _search(session: requests.Session, query: str, *, email: str, api_key: str | None, logger: logging.Logger) -> tuple[int, str | None, str | None]:
    params = {"db": "pubmed", "term": query, "usehistory": "y", "retmax": 0, "email": email, "tool": "pathogen_atlas_v4"}
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
            return count, root.findtext("QueryKey"), root.findtext("WebEnv")
        except (requests.RequestException, ET.ParseError, RuntimeError) as exc:
            if attempt == 6:
                raise
            delay = min(2 ** attempt, 60)
            logger.warning("ESearch retry %d/6 after %s; sleeping %ds", attempt + 1, type(exc).__name__, delay)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _fetch_history(session: requests.Session, query_key: str, webenv: str, count: int, *, email: str, api_key: str | None, logger: logging.Logger):
    params = {"db": "pubmed", "query_key": query_key, "WebEnv": webenv, "retmode": "xml", "email": email, "tool": "pathogen_atlas_v4"}
    if api_key:
        params["api_key"] = api_key
    articles = []
    for start in range(0, count, BATCH_SIZE):
        params.update({"retstart": start, "retmax": min(BATCH_SIZE, count - start)})
        for attempt in range(7):
            try:
                root = _xml(session.post(f"{EUTILS}/efetch.fcgi", data=params, timeout=120))
                articles.extend(_parse_articles(root))
                break
            except (requests.RequestException, ET.ParseError) as exc:
                if attempt == 6:
                    raise
                delay = min(2 ** attempt, 60)
                logger.warning("EFetch retry at record %d after %s; sleeping %ds", start, type(exc).__name__, delay)
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
    """Retrieve complete history sets; recursively split overlarge result sets by publication date."""
    count, key, webenv = _search(session, query, email=email, api_key=api_key, logger=logger)
    if count <= 10000:
        if not count:
            return [], []
        return _fetch_history(session, key, webenv, count, email=email, api_key=api_key, logger=logger), []
    logger.warning("Search returned %d records; segmenting by publication year to avoid PubMed's 10,000 UID retrieval ceiling", count)
    from datetime import date
    year = date.today().year
    all_articles = []
    segment_counts = []
    for y in range(1800, year + 1):
        subquery = _date_range_query(query, f"{y}/01/01", f"{y}/12/31")
        n, subkey, subwebenv = _search(session, subquery, email=email, api_key=api_key, logger=logger)
        if n > 10000:
            for month in range(1, 13):
                import calendar
                last = calendar.monthrange(y, month)[1]
                monthly = _date_range_query(query, f"{y}/{month:02d}/01", f"{y}/{month:02d}/{last:02d}")
                mn, mkey, mweb = _search(session, monthly, email=email, api_key=api_key, logger=logger)
                if mn > 10000:
                    for day in range(1, last + 1):
                        daily = _date_range_query(query, f"{y}/{month:02d}/{day:02d}", f"{y}/{month:02d}/{day:02d}")
                        dn, dkey, dweb = _search(session, daily, email=email, api_key=api_key, logger=logger)
                        if dn > 10000:
                            raise RuntimeError(f"A single publication date still has {dn} PubMed records; cannot guarantee complete retrieval for query {query!r}")
                        if dn:
                            all_articles.extend(_fetch_history(session, dkey, dweb, dn, email=email, api_key=api_key, logger=logger))
                        segment_counts.append({"period": f"{y}-{month:02d}-{day:02d}", "count": dn})
                elif mn:
                    all_articles.extend(_fetch_history(session, mkey, mweb, mn, email=email, api_key=api_key, logger=logger))
                    segment_counts.append({"period": f"{y}-{month:02d}", "count": mn})
        elif n:
            all_articles.extend(_fetch_history(session, subkey, subwebenv, n, email=email, api_key=api_key, logger=logger))
            segment_counts.append({"period": str(y), "count": n})
    unique = {row["pmid"]: row for row in all_articles}
    return [unique[k] for k in sorted(unique, key=lambda x: int(x))], segment_counts


def acquire_pubmed_pairs(pairs: Iterable[tuple[str, str]], pathogens: pd.DataFrame, hosts: pd.DataFrame,
                         aliases: pd.DataFrame, output_dir: str | Path, *, query_version: str,
                         email: str, api_key: str | None = None, logger: logging.Logger | None = None) -> dict[str, object]:
    """Search/fetch every pair, checkpoint each completed pair, and materialize normalized outputs."""
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
            if (pid, hid) in done and done[(pid, hid)].get("status") == "completed":
                continue
            pathogen, host = pidx.loc[pid], hidx.loc[hid]
            pnames = [pathogen.accepted_name, *aliases_by_pathogen.get(pid, [])]
            haliases = host.get("aliases", [])
            if hasattr(haliases, "tolist"):
                haliases = haliases.tolist()
            query_record = build_pair_query(pnames, [host.accepted_name, *(haliases or [])], strategy_version=query_version)
            articles, segments = _all_pmids(session, query_record["query"], email=email, api_key=api_key, logger=logger)
            pmids = sorted({str(row["pmid"]) for row in articles}, key=lambda x: int(x))
            for row in articles:
                if row["pmid"] not in article_cache:
                    article_cache[row["pmid"]] = row
                    article_out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            article_out.flush()
            record = {"pathogen_id": str(pid), "host_id": str(hid), "relationship_context": "disease", "query": query_record["query"],
                      "query_hash": query_record["query_hash"], "query_version": query_record["query_version"], "pmids": pmids,
                      "search_count": len(pmids), "segmented_search": bool(segments), "segments": segments, "status": "completed"}
            out.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            out.flush()
            done[(pid, hid)] = record
            completed_this_run += 1
            if completed_this_run % 25 == 0 or index == total:
                elapsed = max(time.monotonic() - started, 0.001)
                rate = completed_this_run / elapsed
                remaining = max(total - index, 0)
                eta_minutes = remaining / rate / 60 if rate else float("inf")
                logger.info("PubMed pairs %d/%d (%.1f%%); this session %.2f pairs/s; ETA %.1f min; checkpointed=%d PMIDs=%d",
                            index, total, index * 100 / max(total, 1), rate, eta_minutes, len(done), len(article_cache))
            time.sleep(0.11 if api_key else 0.36)
    pair_rows = []
    for pair in selected_pairs:
        record = done.get(pair)
        if record:
            pair_rows.append({k: record.get(k) for k in ("pathogen_id", "host_id", "relationship_context", "query", "query_hash", "query_version", "search_count", "segmented_search", "status")})
    match_rows = []
    for pair in selected_pairs:
        record = done.get(pair)
        if record:
            for pmid in record["pmids"]:
                match_rows.append({"pathogen_id": pair[0], "host_id": pair[1], "relationship_context": "disease", "pmid": str(pmid), "query_hash": record["query_hash"], "query_version": record["query_version"]})
    out_root = root / "tables"
    out_root.mkdir(exist_ok=True)
    pd.DataFrame(pair_rows).sort_values(["pathogen_id", "host_id"]).to_parquet(out_root / "pair_searches.parquet", index=False)
    pd.DataFrame(match_rows, columns=["pathogen_id", "host_id", "relationship_context", "pmid", "query_hash", "query_version"]).sort_values(["pathogen_id", "host_id", "pmid"]).to_parquet(out_root / "pathogen_host_article_matches.parquet", index=False)
    pd.DataFrame([article_cache[key] for key in sorted(article_cache, key=lambda x: int(x))]).to_parquet(out_root / "articles.parquet", index=False)
    return {"pairs_total": total, "pairs_complete": len(pair_rows), "pairs_run_this_session": completed_this_run,
            "pair_matches": len(match_rows), "unique_articles": len(article_cache), "output_dir": str(out_root), "checkpoint": str(checkpoint)}
