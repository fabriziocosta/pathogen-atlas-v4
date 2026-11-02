"""Resumable Wikipedia and GBIF references for V4 pathogen taxonomy records."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import pandas as pd

GROUP_TO_KINGDOM = {"bacteria": "Bacteria", "fungi": "Fungi", "viruses": "Viruses"}
WIKIPEDIA_FIELDS = ["pathogen_id", "pathogen_group", "accepted_name", "wikipedia_title", "wikipedia_url", "match_status", "match_type"]
GBIF_FIELDS = ["pathogen_id", "pathogen_group", "accepted_name", "gbif_taxon_key", "gbif_scientific_name", "gbif_url", "match_status", "match_type", "confidence"]


def _request_json(request: Request, *, timeout: int, retries: int = 7):
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"Reference lookup request failed: {error}") from error
            retry_after = error.headers.get("Retry-After") if isinstance(error, HTTPError) else None
            time.sleep(float(retry_after) if retry_after and retry_after.isdigit() else min(2 * (2**attempt), 60))
    raise AssertionError("unreachable")


def _wikipedia_batch(names: list[str]) -> dict[str, dict[str, str]]:
    query = urlencode({"action": "query", "format": "json", "formatversion": "2", "redirects": "1", "maxlag": "5", "titles": "|".join(names)})
    request = Request("https://en.wikipedia.org/w/api.php?" + query, headers={"User-Agent": "PathogenAtlasV4/1.0 (taxonomy reference link resolver)"})
    payload = _request_json(request, timeout=30)
    if not isinstance(payload, dict) or payload.get("error"):
        raise RuntimeError(f"Wikipedia API returned an error: {payload.get('error') if isinstance(payload, dict) else payload}")
    query_data = payload.get("query", {})
    title_map = {mapping["from"].casefold(): mapping["to"] for kind in ("normalized", "redirects") for mapping in query_data.get(kind, [])}
    found_by_title = {page["title"].casefold(): page["title"] for page in query_data.get("pages", []) if "missing" not in page}
    resolved = {}
    for requested in names:
        title, visited = requested, set()
        while title.casefold() in title_map and title.casefold() not in visited:
            visited.add(title.casefold())
            title = title_map[title.casefold()]
        canonical = found_by_title.get(title.casefold())
        resolved[requested] = ({"wikipedia_title": canonical,
                                "wikipedia_url": "https://en.wikipedia.org/wiki/" + quote(canonical.replace(" ", "_"), safe="()_,-."),
                                "match_status": "found", "match_type": "exact_or_redirect"}
                               if canonical else {"wikipedia_title": "", "wikipedia_url": "", "match_status": "not_found", "match_type": "exact_title_query"})
    return resolved


def _gbif_batch(records: list[dict[str, str]]) -> dict[str, dict[str, str | int | float]]:
    body = json.dumps([{"scientificName": row["accepted_name"], "kingdom": GROUP_TO_KINGDOM[row["pathogen_group"]]} for row in records]).encode()
    request = Request("https://api.gbif.org/v2/species/match", data=body,
                      headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "PathogenAtlasV4/1.0 (GBIF reference link resolver)"}, method="POST")
    payload = _request_json(request, timeout=60)
    if not isinstance(payload, list) or len(payload) != len(records):
        raise RuntimeError(f"GBIF returned {len(payload) if isinstance(payload, list) else 'an invalid number of'} results for {len(records)} names.")
    output = {}
    for row, result in zip(records, payload):
        usage, diagnostics = result.get("usage") or {}, result.get("diagnostics") or {}
        match_type = str(diagnostics.get("matchType") or "NONE")
        canonical, key = str(usage.get("canonicalName") or ""), str(usage.get("key") or "")
        exact = bool(key) and str(usage.get("rank") or "").upper() == "SPECIES" and match_type == "EXACT" and canonical.casefold() == row["accepted_name"].casefold()
        output[row["pathogen_id"]] = {"gbif_taxon_key": key if exact else "", "gbif_scientific_name": str(usage.get("name") or "") if exact else "",
            "gbif_url": f"https://www.gbif.org/species/{quote(key, safe='')}" if exact else "", "match_status": "found" if exact else "not_found",
            "match_type": match_type, "confidence": diagnostics.get("confidence", "")}
    return output


def _selected_records(pathogen_table: str | Path, groups: list[str], pathogen_ids: set[str] | None) -> list[dict[str, str]]:
    unknown = sorted(set(groups) - set(GROUP_TO_KINGDOM))
    if unknown:
        raise ValueError(f"Unsupported pathogen groups: {unknown}")
    frame = pd.read_parquet(pathogen_table, columns=["pathogen_id", "pathogen_group", "accepted_name", "rank"])
    frame = frame.loc[frame.pathogen_group.isin(groups) & frame["rank"].astype(str).str.casefold().eq("species")]
    if pathogen_ids is not None:
        frame = frame.loc[frame.pathogen_id.astype(str).isin(pathogen_ids)]
    frame = frame.loc[frame.accepted_name.notna() & frame.pathogen_id.notna()].drop_duplicates("pathogen_id")
    return [{"pathogen_id": str(r.pathogen_id), "pathogen_group": str(r.pathogen_group), "accepted_name": str(r.accepted_name).strip()}
            for r in frame.itertuples(index=False) if str(r.accepted_name).strip()]


def _write(rows: dict[str, dict], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame((rows[key] for key in sorted(rows)), columns=fields).to_csv(path, index=False)


def _resolve(pathogen_table: str | Path, output_path: str | Path, groups: list[str], *, force_refresh: bool,
             fields: list[str], batch_size: int, delay: float, lookup: Callable, pathogen_ids: set[str] | None) -> pd.DataFrame:
    records = _selected_records(pathogen_table, groups, pathogen_ids)
    if not records:
        raise ValueError(f"No species-rank pathogen records found for {groups} in {pathogen_table}")
    path, rows = Path(output_path), {}
    selected_ids, selected_groups = {r["pathogen_id"] for r in records}, set(groups)
    if path.exists():
        previous = pd.read_csv(path, keep_default_na=False, dtype=str)
        if set(fields).issubset(previous.columns):
            for item in previous.to_dict("records"):
                pid = str(item["pathogen_id"])
                if item.get("match_status") in {"found", "not_found"} and (item.get("pathogen_group") not in selected_groups or (pid in selected_ids and not force_refresh)):
                    rows[pid] = {field: item[field] for field in fields}
    pending = [r for r in records if r["pathogen_id"] not in rows]
    print(f"Resolving {len(pending):,} pathogen records; reusing {len(rows):,} completed rows.", flush=True)
    for start in range(0, len(pending), batch_size):
        batch = pending[start:start + batch_size]
        if fields is WIKIPEDIA_FIELDS:
            names = list(dict.fromkeys(r["accepted_name"] for r in batch))
            results = lookup(names)
            for record in batch:
                rows[record["pathogen_id"]] = {**record, **results[record["accepted_name"]]}
        else:
            results = lookup(batch)
            for record in batch:
                rows[record["pathogen_id"]] = {**record, **results[record["pathogen_id"]]}
        _write(rows, path, fields)
        print(f"Resolved {min(start + len(batch), len(pending)):,}/{len(pending):,} pending records.", flush=True)
        if delay:
            time.sleep(delay)
    _write(rows, path, fields)
    print(f"Saved {len(rows):,} associations to {path}", flush=True)
    return pd.read_csv(path, keep_default_na=False)


def resolve_wikipedia_links(pathogen_table: str | Path, output_path: str | Path, pathogen_groups: list[str], *, force_refresh: bool = False, pathogen_ids: set[str] | None = None) -> pd.DataFrame:
    return _resolve(pathogen_table, output_path, pathogen_groups, force_refresh=force_refresh, fields=WIKIPEDIA_FIELDS, batch_size=25, delay=1.1, lookup=_wikipedia_batch, pathogen_ids=pathogen_ids)


def resolve_gbif_links(pathogen_table: str | Path, output_path: str | Path, pathogen_groups: list[str], *, force_refresh: bool = False, pathogen_ids: set[str] | None = None) -> pd.DataFrame:
    return _resolve(pathogen_table, output_path, pathogen_groups, force_refresh=force_refresh, fields=GBIF_FIELDS, batch_size=100, delay=0.25, lookup=_gbif_batch, pathogen_ids=pathogen_ids)
