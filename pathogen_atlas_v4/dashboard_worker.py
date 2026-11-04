"""Restartable, checkpoint-oriented runner for the local V4 dashboard."""

from __future__ import annotations

import argparse
import asyncio
import json
from datetime import datetime
from statistics import median
import os
import re
import sys
import time
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import nbformat
import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
STATE_ROOT = ROOT / "data/v4/dashboard/runs"
PHASES = [
    ("setup", "00_configure_run.ipynb"), ("taxonomy", "01_review_taxonomy_and_hosts.ipynb"),
    ("prompts", "02_preview_prompts.ipynb"), ("search", "03_acquire_pair_pubmed_records.ipynb"),
    ("screen", "04_screen_pathogen_host_pairs.ipynb"), ("estimates", "05_estimate_pathogen_publications.ipynb"),
    ("handoff", "06_handoff_disease_evidence.ipynb"), ("profiles", "07_resolve_feature_profiles.ipynb"),
    ("scoring", "08_score_feature_profiles.ipynb"), ("validation", "09_validate_release.ipynb"),
    ("wikipedia", "11_resolve_wikipedia_references.ipynb"), ("gbif", "12_resolve_gbif_references.ipynb"),
    ("build", "10_inspect_embeddings.ipynb"),
]


def read_json(path, default=None):
    try: return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError): return default


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    temp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    temp.replace(path)


def emit(run_dir, kind, **data):
    event = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "kind": kind, **data}
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(event, ensure_ascii=False, default=str) + "\n")
    manifest = read_json(run_dir / "manifest.json", {})
    manifest["updated_at"] = event["time"]
    save_json(run_dir / "manifest.json", manifest)


def status(run_dir, run_status=None, phase=None, **updates):
    manifest = read_json(run_dir / "manifest.json", {})
    if run_status: manifest["status"] = run_status
    if phase:
        manifest["phases"].setdefault(phase, {}).update(updates)
    save_json(run_dir / "manifest.json", manifest)


def should_pause(run_dir):
    return (run_dir / "pause.request").exists()


def _patch_notebook(notebook, key, config):
    env_config = "__import__('pathlib').Path(__import__('os').environ['PATHOGEN_ATLAS_RUN_CONFIG'])"
    smoke = config.get("mode") == "smoke_test"
    smoke_root = f"ROOT / 'data/v4/dashboard/smoke-work' / {config['run_id']!r}"
    host_lookup = yaml.safe_load((ROOT / "configs/host_choices.yaml").read_text(encoding="utf-8"))
    host_lookup = host_lookup.get("host_name_to_id", host_lookup)
    configured_host_names = config.get("host_names") or []
    host_names = []
    for index, host_id_value in enumerate(config["hosts"]):
        configured_name = configured_host_names[index] if index < len(configured_host_names) else None
        if configured_name and str(host_lookup.get(configured_name)) == str(host_id_value):
            host_names.append(configured_name)
        else:
            host_names.append(next(name for name, host_id in host_lookup.items()
                                  if str(host_id) == str(host_id_value)))
    group_names = config["groups"]
    for cell in notebook.cells:
        if cell.cell_type != "code": continue
        src = cell.source
        # Every downstream notebook reads this path; keep run configuration entirely under ignored data/v4.
        src = re.sub(r"\b(?:PROJECT_ROOT|ROOT)\s*/\s*['\"]configs/pilot_v1\.yaml['\"]", env_config, src)
        if smoke:
            # Isolate all mutable smoke outputs from full-run caches and release files.
            for shared_path, isolated_path in (
                ("ROOT / 'data/v4/runs'", f"{smoke_root} / 'runs'"),
                ('ROOT / "data/v4/runs"', f"{smoke_root} / 'runs'"),
                ("ROOT / 'data/v4/cache/literature'", f"{smoke_root} / 'cache/literature'"),
                ('ROOT / "data/v4/cache/literature"', f"{smoke_root} / 'cache/literature'"),
                ("ROOT / 'data/v4/cache/decisions'", f"{smoke_root} / 'cache/decisions'"),
                ('ROOT / "data/v4/cache/decisions"', f"{smoke_root} / 'cache/decisions'"),
                ("ROOT / 'data/v4/cache/decisions", f"{smoke_root} / 'cache/decisions"),
                ('ROOT / "data/v4/cache/decisions', f"{smoke_root} / 'cache/decisions"),
                ("ROOT / 'data/v4/references'", f"{smoke_root} / 'references'"),
                ('ROOT / "data/v4/references"', f"{smoke_root} / 'references'"),
                ("ROOT / 'data/v4/references", f"{smoke_root} / 'references"),
                ('ROOT / "data/v4/references', f"{smoke_root} / 'references"),
                ("ROOT / 'data/v4/plots'", f"{smoke_root} / 'plots'"),
                ('ROOT / "data/v4/plots"', f"{smoke_root} / 'plots'"),
            ):
                src = src.replace(shared_path, isolated_path)
            if key == "setup":
                src = src.replace("'selection_mode': RUN_MODE,",
                                  "'selection_mode': ('full' if RUN_MODE == 'smoke_test' else RUN_MODE),")
                src = src.replace("'screen_pathogen_pairs': RUN_MODE == 'full',",
                                  "'screen_pathogen_pairs': True,")
            if key == "build":
                production_site = "SITE_ROOT = ROOT / 'sites/pathogen-atlas-explorer'"
                isolated_site = (
                    f"SITE_ROOT = {smoke_root} / 'site'\n"
                    "import shutil\n"
                    "if not SITE_ROOT.exists():\n"
                    "    shutil.copytree(ROOT / 'sites/pathogen-atlas-explorer', SITE_ROOT, "
                    "ignore=shutil.ignore_patterns('dist', 'taxa', 'explorer-data.json.gz', "
                    "'evidence.json.gz', 'feature-catalog.json'))"
                )
                src = src.replace(production_site, isolated_site)
        if key == "setup" and config.get("taxonomy_snapshot"):
            src = src.replace("TAXONOMY_ROOT = sorted((ROOT / 'data/v4/staging').glob('*/taxonomy'))[-1]",
                              f"TAXONOMY_ROOT = ROOT / 'data/v4/staging' / {config['taxonomy_snapshot']!r} / 'taxonomy'")
        if key == "setup":
            src = src.replace("'decision_workers': 1024", "'decision_workers': 16")
        if key == "setup" and "RUN_ID = 'human-parasites" in src:
            smoke_ids = [config["smoke_pathogen_id"]] if config.get("mode") == "smoke_test" and config.get("smoke_pathogen_id") else []
            substitutions = {
                r"^RUN_ID = .*?$": f"RUN_ID = {config['run_id']!r}",
                r"^RUN_MODE = .*?$": f"RUN_MODE = {config.get('mode', 'full')!r}",
                r"^SELECTED_HOST_NAMES = .*?$": f"SELECTED_HOST_NAMES = {host_names!r}",
                r"^SELECTED_PATHOGEN_GROUPS = .*?$": f"SELECTED_PATHOGEN_GROUPS = {group_names!r}",
                r"^SMOKE_PATHOGEN_IDS = .*?$": f"SMOKE_PATHOGEN_IDS = {smoke_ids!r}",
                r"^MAX_PMIDS_PER_PAIR = .*?$": f"MAX_PMIDS_PER_PAIR = {config['sample_per_pair']!r}",
                r"^OLDEST_PMIDS_PER_PAIR = .*?$": f"OLDEST_PMIDS_PER_PAIR = min({1 if smoke else 30}, {config['sample_per_pair']!r})",
                r"^SELECTION_SEED = .*?$": f"SELECTION_SEED = {config['run_id']!r}",
                r"^HOST_SEARCH_EXPANSION_ENABLED = .*?$": "HOST_SEARCH_EXPANSION_ENABLED = True",
            }
            for pattern, replacement in substitutions.items(): src = re.sub(pattern, replacement, src, flags=re.M)
            src = src.replace("1 <= MAX_PMIDS_PER_PAIR <= 200", "MAX_PMIDS_PER_PAIR >= 1")
            src = src.replace("MAX_PMIDS_PER_PAIR must be an integer from a positive integer; all-hit abstract screening is disabled",
                              "MAX_PMIDS_PER_PAIR must be a positive integer; all-hit abstract screening is disabled")
            src = src.replace("(maximum 200)", "")
            src = src.replace("(maximum 200) |", "|")
            src = src.replace("(maximum 200)", "")
            src = src.replace("hard maximum abstract-screen sample", "abstract-screen sample")
        src = src.replace("1 <= SAMPLE_MAXIMUM <= 200", "SAMPLE_MAXIMUM >= 1")
        src = src.replace("Full workflow screens at most 200 abstracts per pathogen–host pair; set max_pmids_per_pair to an integer from 1 to 200.",
                          "Full workflow requires a positive per-pair PMID sample maximum.")
        src = src.replace("1 to 200", "a positive integer")
        src = src.replace("<= 200", ">= 1") if "MAX_PMIDS_PER_PAIR" in src or "SAMPLE_MAXIMUM" in src else src
        if key == "setup" and "CONFIG_PATH = ROOT / 'configs/pilot_v1.yaml'" in src:
            src = src.replace("CONFIG_PATH = ROOT / 'configs/pilot_v1.yaml'", f"CONFIG_PATH = {env_config}")
        if key == "setup" and "MAX_PMIDS_PER_PAIR =" in src:
            src = "import os\nfrom pathlib import Path\n" + src
        # Inject a cooperative safe-checkpoint stop in provider-backed decision loops.
        if key == "screen" and "LEDGER = run_decisions(" in src:
            src = src.replace("progress_callback=LOG.info)",
                              "progress_callback=LOG.info, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists())")
            if "import os" not in src: src = "import os\nfrom pathlib import Path\n" + src
        if key == "screen" and "prepare_pair_screen_samples(" in src:
            src = src.replace("api_key=NCBI_API_KEY, logger=LOG)",
                              "api_key=NCBI_API_KEY, logger=LOG, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists())")
            if "import os" not in src: src = "import os\nfrom pathlib import Path\n" + src
        if key == "search" and "acquire_pubmed_pairs(" in src:
            src = src.replace("max_workers=CONFIG['acquisition'].get('max_concurrent_requests', 8),",
                              "max_workers=CONFIG['acquisition'].get('max_concurrent_requests', 8),\n    refresh_id=RUN_ID, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists(),")
            if "import os" not in src: src = "import os\nfrom pathlib import Path\n" + src
        if key == "estimates" and "count_accepted_pathogen_publications(" in src:
            src = src.replace("email=NCBI_EMAIL, api_key=NCBI_API_KEY, logger=LOG)",
                              "email=NCBI_EMAIL, api_key=NCBI_API_KEY, logger=LOG, refresh_id=RUN_ID, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists())")
            if "import os" not in src: src = "import os\nfrom pathlib import Path\n" + src
        if key == "scoring" and "FEATURE_LEDGER = decision_runner.run_decisions(" in src:
            src = src.replace("progress_callback=LOG.info)",
                              "progress_callback=LOG.info, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists())")
            if "import os" not in src: src = "import os\nfrom pathlib import Path\n" + src
            anchor = "FEATURE_RAW_PATH = RUN_DIR / 'feature_screen_raw_responses.jsonl'"
            shared = """FEATURE_RAW_PATH = RUN_DIR / 'feature_screen_raw_responses.jsonl'
GLOBAL_FEATURE_LEDGER_PATH = ROOT / 'data/v4/cache/decisions/feature_screen_ledger.parquet'
if GLOBAL_FEATURE_LEDGER_PATH.exists() and FEATURE_TASKS:
    global_feature_ledger = pd.read_parquet(GLOBAL_FEATURE_LEDGER_PATH)
    requested_feature_keys = {task['cache_key'] for task in FEATURE_TASKS}
    reusable = global_feature_ledger.loc[
        global_feature_ledger.cache_key.astype(str).isin(requested_feature_keys)
        & global_feature_ledger.status.astype(str).eq('completed')
    ].copy()
    if not reusable.empty:
        reusable['run_id'] = RUN_ID
        prior_run_ledger = pd.read_parquet(FEATURE_LEDGER_PATH) if FEATURE_LEDGER_PATH.exists() else pd.DataFrame()
        FEATURE_LEDGER_PATH.parent.mkdir(parents=True, exist_ok=True)
        pd.concat([prior_run_ledger, reusable], ignore_index=True).drop_duplicates('cache_key', keep='last').to_parquet(FEATURE_LEDGER_PATH, index=False)
        LOG.info('Reused %d semantically matching feature decisions from the shared cache', len(reusable))"""
            if not smoke:
                src = src.replace(anchor, shared)
            merge_anchor = "task_keys = {task['cache_key'] for task in FEATURE_TASKS}"
            merge_code = """# Extend the cross-run cache only with semantically keyed completed decisions.
    if not FEATURE_LEDGER.empty:
        GLOBAL_FEATURE_LEDGER_PATH.parent.mkdir(parents=True, exist_ok=True)
        prior_global = pd.read_parquet(GLOBAL_FEATURE_LEDGER_PATH) if GLOBAL_FEATURE_LEDGER_PATH.exists() else pd.DataFrame()
        combined_global = pd.concat([prior_global, FEATURE_LEDGER], ignore_index=True).drop_duplicates('cache_key', keep='last')
        global_temp_path = GLOBAL_FEATURE_LEDGER_PATH.with_suffix('.tmp.parquet')
        combined_global.to_parquet(global_temp_path, index=False)
        global_temp_path.replace(GLOBAL_FEATURE_LEDGER_PATH)
    task_keys = {task['cache_key'] for task in FEATURE_TASKS}"""
            if not smoke:
                src = src.replace(merge_anchor, merge_code)
        cell.source = src
    return notebook


def _output_text(cell):
    chunks = []
    for output in cell.get("outputs", []):
        if output.get("output_type") == "stream": chunks.append(output.get("text", ""))
        elif output.get("output_type") == "error": chunks.append("\n".join(output.get("traceback", [])))
        elif "text/plain" in output.get("data", {}): chunks.append(output["data"]["text/plain"])
    return "\n".join(chunks)


def _capture_metrics(run_dir, phase, line):
    manifest = read_json(run_dir / "manifest.json", {})
    metric = manifest.get("phases", {}).get(phase, {})
    search = re.search(r"PubMed pair checkpoints\s+(\d+)/(\d+).*?ETA\s+([^;.]*)", line)
    decision_start = re.search(r"Decision screening:\s+[\d,]+ pending;\s*([\d,]+) cached;", line)
    decision = re.search(r"Decision progress:\s*([\d,]+)/([\d,]+) this run.*?cached\s+([\d,]+).*?ETA ~([^|.]+)", line)
    decision_complete = re.search(r"Decision screening complete:\s*([\d,]+) completed this run", line)
    reported_cost = re.search(r"reported cost \$([0-9.]+)", line)
    cost_coverage = re.search(r"cost coverage\s+([\d,]+)/([\d,]+)", line)
    cost_unavailable = "reported cost unavailable" in line
    projected_cost = re.search(r"projected run cost \$([0-9.]+)", line)
    search_plan = re.search(r"pairs=(\d+) already_checkpointed=(\d+) searches_to_run=(\d+)", line)
    query_preparation = re.search(r"PubMed query preparation (\d+)/(\d+)", line)
    screen_tasks = re.search(r"Built pathogen screen tasks=(\d[\d,]*)", line)
    feature_tasks = re.search(r"Feature scoring queued tasks=(\d[\d,]*)", line)
    setup = re.search(r"candidate pairs.?\s*(\d[\d,]*)", line, re.I)
    if not any((search, decision_start, decision, decision_complete, reported_cost, search_plan, query_preparation,
                screen_tasks, feature_tasks, setup)):
        return
    if decision_start:
        metric["cost_segment_base_usd"] = float(metric.get("cost_usd") or 0)
        metric["cost_segment_base_known"] = int(metric.get("cost_tasks_known", 0))
        metric["cost_segment_base_total"] = int(metric.get("cost_tasks_total", 0))
    if search:
        done, total = map(int, search.group(1, 2))
        metric.update({"completed": done, "total": total, "progress": round(done / max(total, 1) * 100, 1),
                       "eta_seconds": _parse_eta(search.group(3))})
    elif decision:
        done, total = (int(value.replace(",", "")) for value in decision.group(1, 2))
        cached = int(decision.group(3).replace(",", ""))
        metric.update({"completed": done + cached, "total": total + cached,
                       "progress": round((done + cached) / max(total + cached, 1) * 100, 1),
                       "count_label": "scored features" if phase == "scoring" else "decisions",
                       "cache_hits": cached, "eta_seconds": _parse_eta(decision.group(4))})
    if decision_complete:
        metric.update({"completed": int(decision_complete.group(1).replace(",", "")),
                       "total": int(decision_complete.group(1).replace(",", "")),
                       "count_label": "decisions", "progress": 100})
    if reported_cost:
        metric["cost_usd"] = float(metric.get("cost_segment_base_usd", 0)) + float(reported_cost.group(1))
    elif cost_unavailable:
        if not metric.get("cost_tasks_known", 0):
            metric["cost_usd"] = None
    if cost_coverage:
        metric["cost_tasks_known"] = (int(metric.get("cost_segment_base_known", 0))
                                      + int(cost_coverage.group(1).replace(",", "")))
        metric["cost_tasks_total"] = (int(metric.get("cost_segment_base_total", 0))
                                     + int(cost_coverage.group(2).replace(",", "")))
    if projected_cost:
        metric["estimate_usd"] = float(projected_cost.group(1))
    elif search_plan:
        total, cached, remaining = map(int, search_plan.groups())
        metric.update({"completed": cached, "total": total, "count_label": "pair searches",
                       "cache_hits": cached, "progress": round(cached / max(total, 1) * 100, 1)})
    elif query_preparation:
        completed, total = map(int, query_preparation.groups())
        metric.update({"completed": completed, "total": total, "count_label": "queries prepared",
                       "progress": round(completed / max(total, 1) * 100, 1)})
    elif screen_tasks:
        total = int(screen_tasks.group(1).replace(",", ""))
        metric.update({"completed": 0, "total": total, "count_label": "screen tasks", "progress": 0})
    elif feature_tasks:
        total = int(feature_tasks.group(1).replace(",", ""))
        metric.update({"completed": 0, "total": total, "count_label": "scored features", "progress": 0})
    elif setup:
        metric.update({"total": int(setup.group(1).replace(",", "")), "completed": 0})
    manifest["phases"][phase] = metric
    save_json(run_dir / "manifest.json", manifest)


def _parse_eta(value):
    value = value.strip().replace("~", "")
    seconds = 0
    for number, unit in re.findall(r"(\d+)\s*(h|m|s)", value):
        seconds += int(number) * {"h": 3600, "m": 60, "s": 1}[unit]
    return seconds if seconds else None


def _usd_to_gbp_rate():
    """Return a dated USD-to-GBP ECB reference rate, or None when unavailable."""
    try:
        response = requests.get("https://www.ecb.europa.eu/stats/eurofxref/eurofxref-daily.xml", timeout=5)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        dated_rates = root.find(".//{*}Cube[@time]")
        if dated_rates is None:
            return None
        rates = {node.attrib.get("currency"): float(node.attrib["rate"])
                 for node in dated_rates if node.attrib.get("currency") in {"USD", "GBP"}}
        if not rates.get("USD") or not rates.get("GBP"):
            return None
        return {"usd_to_gbp": rates["GBP"] / rates["USD"],
                "date": dated_rates.attrib["time"], "source": "ECB daily reference rates"}
    except (requests.RequestException, ET.ParseError, KeyError, ValueError):
        return None


def _finalize_costs(manifest):
    cost_phases = [manifest.get("phases", {}).get(key, {}) for key in ("screen", "scoring")]
    known_costs = [float(phase["cost_usd"]) for phase in cost_phases
                   if phase.get("cost_usd") is not None and phase.get("cost_tasks_known", 0) > 0]
    known_tasks = sum(int(phase.get("cost_tasks_known", 0)) for phase in cost_phases)
    total_tasks = sum(int(phase.get("cost_tasks_total", 0)) for phase in cost_phases)
    manifest["total_cost_usd"] = sum(known_costs) if known_costs else None
    manifest["total_cost_tasks_known"] = known_tasks
    manifest["total_cost_tasks_total"] = total_tasks
    fx = manifest.get("fx_rate") or {}
    rate = fx.get("usd_to_gbp")
    manifest["total_cost_gbp"] = (manifest["total_cost_usd"] * rate
                                   if manifest["total_cost_usd"] is not None and rate else None)


def _phase_counts(phase, run_dir, run_id):
    import pandas as pd
    manifest = read_json(run_dir / "manifest.json", {})
    smoke = manifest.get("config", {}).get("mode") == "smoke_test"
    if smoke:
        smoke_root = ROOT / "data/v4/dashboard/smoke-work" / run_id
        run = smoke_root / "runs" / run_id
        literature_root = smoke_root / "cache/literature"
        references_root = smoke_root / "references"
        site_root = smoke_root / "site"
    else:
        run = ROOT / "data/v4/runs" / run_id
        literature_root = ROOT / "data/v4/cache/literature"
        references_root = ROOT / "data/v4/references"
        site_root = ROOT / "sites/pathogen-atlas-explorer"
    def rows(path):
        return len(pd.read_parquet(path)) if path.exists() else None
    if phase == "setup":
        resolved = read_json(run / "resolved_config.json", {})
        if resolved.get("pairs"):
            return ("pairs", len(resolved["pairs"]))
        config_path = run_dir / "dashboard-run.yaml"
        run_config = yaml.safe_load(config_path.read_text(encoding="utf-8")) if config_path.exists() else {}
        run_settings = run_config.get("run", {})
        pair_count = len(run_settings.get("pairs", []))
        if not pair_count:
            pair_count = len(run_settings.get("host_ids", [])) * len(run_settings.get("pathogen_ids", []))
        return ("pairs", pair_count)
    if phase == "taxonomy":
        count = rows(run / "taxonomy/pathogens.parquet")
        return ("taxa", count)
    if phase == "prompts":
        return ("prompt previews", 1)
    if phase == "search":
        resolved = read_json(run / "resolved_config.json", {})
        source = resolved.get("source_config", {})
        acquisition = source.get("acquisition", {})
        from pathogen_atlas_v4.pubmed_client import literature_cache_key
        cache_key = literature_cache_key(
            query_version=acquisition.get("query_strategy_version", ""),
            pathogen_literature_query_version=acquisition.get("pathogen_literature_query_version", ""),
            taxonomy_snapshot=source.get("taxonomy", {}).get("staging_snapshot", ""),
            host_expansion_hash=acquisition.get("host_expansion_hash", ""))
        path = literature_root / cache_key / "tables/pair_searches.parquet"
        return ("pair searches", rows(path))
    if phase == "screen":
        return ("screened pairs", rows(run / "pathogen_host_selection.parquet"))
    if phase == "estimates":
        return ("accepted profiles", rows(run / "release/pathogenic_interest_index.parquet"))
    if phase == "handoff":
        return ("evidence records", rows(run / "evidence/pathogen_host_evidence.parquet"))
    if phase == "profiles":
        folder = run / "feature_profiles"
        return ("profiles", len(list(folder.glob("*.json"))) if folder.exists() else 0)
    if phase == "scoring":
        path = run / "feature_screen_ledger.parquet"
        return ("scored features", rows(path))
    if phase == "validation":
        return ("validation checks", 1)
    if phase == "wikipedia":
        path = references_root / "release/11_resolve_wikipedia_references.csv"
        return ("references", len(pd.read_csv(path)) if path.exists() else 0)
    if phase == "gbif":
        path = references_root / "release/12_resolve_gbif_references.csv"
        return ("references", len(pd.read_csv(path)) if path.exists() else 0)
    if phase == "build":
        payload = site_root / "explorer-data.json.gz"
        return ("Explorer payloads", 1 if payload.exists() else 0)
    return ("items", None)


def _historical_phase_duration(run_dir, phase):
    """Return a median completed duration for this phase and run mode."""
    current = read_json(run_dir / "manifest.json", {})
    mode = current.get("config", {}).get("mode")
    durations = []
    manifests = sorted(STATE_ROOT.glob("*/manifest.json"), key=lambda path: path.stat().st_mtime, reverse=True)
    for path in manifests:
        if path.parent == run_dir:
            continue
        previous = read_json(path, {})
        if previous.get("status") != "completed":
            continue
        if mode and previous.get("config", {}).get("mode") != mode:
            continue
        item = previous.get("phases", {}).get(phase, {})
        if item.get("status") != "completed" or not item.get("started_at") or not item.get("completed_at"):
            continue
        try:
            start = datetime.strptime(item["started_at"], "%Y-%m-%d %H:%M:%S")
            end = datetime.strptime(item["completed_at"], "%Y-%m-%d %H:%M:%S")
        except ValueError:
            continue
        seconds = (end - start).total_seconds()
        if 1 <= seconds <= 7 * 24 * 60 * 60:
            durations.append(seconds)
        if len(durations) == 20:
            break
    return round(median(durations)) if durations else None


async def _execute_notebook(nb, phase, run_dir):
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError
    client = NotebookClient(nb, timeout=None, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}})
    manifest = read_json(run_dir / "manifest.json", {})
    manifest["status"] = "running"
    manifest["phases"][phase].pop("error", None)
    phase_metrics = {"status": "running", "started_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    if phase == "build":
        phase_metrics.update({"count_label": "Explorer payloads", "completed": 0, "total": 1,
                              "estimated_duration_seconds": _historical_phase_duration(run_dir, phase)})
    manifest["phases"][phase].update(phase_metrics)
    save_json(run_dir / "manifest.json", manifest)
    emit(run_dir, "phase", phase=phase, status="running", message=f"Started {manifest['phases'][phase]['label']}.")
    try:
        async with client.async_setup_kernel():
            for index, cell in enumerate(nb.cells):
                if cell.cell_type != "code": continue
                if should_pause(run_dir) and index == 0:
                    raise InterruptedError("Pause requested at phase boundary")
                await client.async_execute_cell(cell, index)
                text = _output_text(cell).strip()
                if text:
                    for line in text.splitlines():
                        if line.strip():
                            _capture_metrics(run_dir, phase, line)
                            emit(run_dir, "log", phase=phase, message=line[:3000])
                cell_error = next((o for o in cell.get("outputs", []) if o.get("output_type") == "error"), None)
                if cell_error: raise CellExecutionError.from_cell_and_msg(cell, cell_error)
                if should_pause(run_dir):
                    raise InterruptedError("Pause requested; phase checkpoint completed")
    except InterruptedError:
        manifest = read_json(run_dir / "manifest.json", {})
        manifest["phases"][phase]["status"] = "paused"
        manifest["status"] = "paused"
        save_json(run_dir / "manifest.json", manifest)
        emit(run_dir, "phase", phase=phase, status="paused", message="Paused at a saved checkpoint.")
        return False
    except Exception as exc:
        if should_pause(run_dir):
            manifest = read_json(run_dir / "manifest.json", {})
            manifest["phases"][phase]["status"] = "paused"
            manifest["status"] = "paused"
            save_json(run_dir / "manifest.json", manifest)
            emit(run_dir, "phase", phase=phase, status="paused",
                 message="Provider work checkpointed; resuming will continue remaining tasks.")
            return False
        manifest = read_json(run_dir / "manifest.json", {})
        manifest["phases"][phase].update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        manifest["status"] = "failed"
        save_json(run_dir / "manifest.json", manifest)
        emit(run_dir, "error", phase=phase, message=f"{type(exc).__name__}: {exc}")
        raise
    manifest = read_json(run_dir / "manifest.json", {})
    # Smoke notebooks are temporary artifacts and are removed with their workspace.
    if manifest.get("config", {}).get("mode") == "smoke_test":
        output_path = ROOT / "data/v4/dashboard/smoke-work" / manifest["run_id"] / "executed" / f"{phase}.executed.ipynb"
    else:
        output_path = run_dir / f"{phase}.executed.ipynb"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, output_path)
    count_label, count = _phase_counts(phase, run_dir, manifest["run_id"])
    manifest["phases"][phase].update({"status": "completed", "progress": 100,
                                      "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                                      "count_label": count_label})
    if count is not None:
        manifest["phases"][phase].update({"completed": count, "total": count})
    save_json(run_dir / "manifest.json", manifest)
    emit(run_dir, "phase", phase=phase, status="completed", message=f"Completed {manifest['phases'][phase]['label']}.")
    return True


def run(run_id):
    run_dir = STATE_ROOT / run_id
    manifest = read_json(run_dir / "manifest.json")
    if not manifest: raise RuntimeError(f"Unknown run: {run_id}")
    config = {**manifest["config"], "run_id": run_id}
    smoke = config.get("mode") == "smoke_test"
    smoke_root = ROOT / "data/v4/dashboard/smoke-work" / run_id
    if smoke:
        smoke_root.mkdir(parents=True, exist_ok=True)
    # Pause requests are removed by the explicit resume endpoint; new runs start without one.
    manifest["status"] = "running"
    manifest.pop("error", None)
    if "fx_rate" not in manifest:
        manifest["fx_rate"] = _usd_to_gbp_rate()
    save_json(run_dir / "manifest.json", manifest)
    emit(run_dir, "run", message="Pipeline worker started.")
    try:
        for phase, filename in PHASES:
            if phase == "wikipedia" and not config.get("wikipedia", True):
                status(run_dir, phase=phase, status="skipped", progress=100)
                emit(run_dir, "phase", phase=phase, status="skipped", message="Skipped by run configuration.")
                continue
            if phase == "gbif" and not config.get("gbif", True):
                status(run_dir, phase=phase, status="skipped", progress=100)
                emit(run_dir, "phase", phase=phase, status="skipped", message="Skipped by run configuration.")
                continue
            current = read_json(run_dir / "manifest.json", {})
            if current.get("phases", {}).get(phase, {}).get("status") == "completed": continue
            if should_pause(run_dir):
                status(run_dir, "paused")
                emit(run_dir, "run", message="Paused before the next phase.")
                return
            nb = nbformat.read(ROOT / "notebooks" / filename, as_version=4)
            _patch_notebook(nb, phase, config)
            env = os.environ
            env["PATHOGEN_ATLAS_RUN_CONFIG"] = str(smoke_root / "dashboard-run.yaml" if smoke else run_dir / "dashboard-run.yaml")
            env["PATHOGEN_ATLAS_PAUSE_FILE"] = str(run_dir / "pause.request")
            # Notebook kernels inherit process environment; subprocess worker is isolated per run.
            complete = asyncio.run(_execute_notebook(nb, phase, run_dir))
            if not complete: return
        manifest = read_json(run_dir / "manifest.json", {})
        manifest["status"] = "completed"
        manifest["publish_ready"] = not smoke
        manifest["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
        _finalize_costs(manifest)
        save_json(run_dir / "manifest.json", manifest)
        if smoke:
            import shutil
            shutil.rmtree(smoke_root, ignore_errors=True)
            shutil.rmtree(ROOT / "data/v4/logs" / run_id, ignore_errors=True)
            (run_dir / "worker.json").unlink(missing_ok=True)
            (run_dir / "pause.request").unlink(missing_ok=True)
            emit(run_dir, "run", message="Smoke pipeline completed end to end; isolated run files and checkpoints were removed.")
        else:
            emit(run_dir, "run", message="Explorer payload and embedding phase completed; release is ready for review.")
    except Exception as exc:
        manifest = read_json(run_dir / "manifest.json", {})
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        save_json(run_dir / "manifest.json", manifest)
        emit(run_dir, "error", message=manifest["error"])
        return


def bump_release(run_id):
    # Reserved for the explicit publish action; update patch version and prepend changelog entry.
    import re
    version_path = ROOT / "pyproject.toml"
    text = version_path.read_text(encoding="utf-8")
    match = re.search(r'^version = "(\d+)\.(\d+)\.(\d+)"', text, re.M)
    if not match: raise RuntimeError("Could not parse V4 version")
    version = f"{match[1]}.{match[2]}.{int(match[3]) + 1}"
    version_path.write_text(text[:match.start(1)] + version + text[match.end(3):], encoding="utf-8")
    init = ROOT / "pathogen_atlas_v4/__init__.py"
    init_text = init.read_text(encoding="utf-8") if init.exists() else ""
    if "__version__" in init_text:
        init_text = re.sub(r'__version__\s*=\s*["\'][^"\']+["\']', f'__version__ = "{version}"', init_text)
        init.write_text(init_text, encoding="utf-8")
    for page in (ROOT / "sites/pathogen-atlas-explorer/index.html",
                 ROOT / "sites/pathogen-atlas-explorer/dist/index.html"):
        if page.exists():
            page_text = page.read_text(encoding="utf-8")
            page_text = re.sub(r"v\d+\.\d+\.\d+", f"v{version}", page_text, count=1)
            page.write_text(page_text, encoding="utf-8")
    changelog = ROOT / "CHANGELOG.md"
    old = changelog.read_text(encoding="utf-8") if changelog.exists() else "# Changelog\n"
    entry = f"## [{version}] - {time.strftime('%Y-%m-%d')}\n\n- Updated Explorer data from dashboard run `{run_id}`.\n\n"
    changelog.write_text(old.replace("# Changelog\n", "# Changelog\n\n" + entry, 1), encoding="utf-8")
    print(version)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id")
    parser.add_argument("--bump-release")
    args = parser.parse_args()
    if args.bump_release: bump_release(args.bump_release)
    elif args.run_id: run(args.run_id)
    else: parser.error("--run-id is required")


if __name__ == "__main__": main()
