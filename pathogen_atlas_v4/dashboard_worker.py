"""Restartable, checkpoint-oriented runner for the local V4 dashboard."""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import re
import sys
import time
from pathlib import Path

import nbformat
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
    temp = path.with_suffix(path.suffix + ".tmp")
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
    host_lookup = yaml.safe_load((ROOT / "configs/host_choices.yaml").read_text(encoding="utf-8"))
    host_lookup = host_lookup.get("host_name_to_id", host_lookup)
    host_names = config.get("host_names") or [
        next(name for name, host_id in host_lookup.items() if str(host_id) == host_id_value)
        for host_id_value in config["hosts"]]
    group_names = config["groups"]
    for cell in notebook.cells:
        if cell.cell_type != "code": continue
        src = cell.source
        # Every downstream notebook reads this path; keep run configuration entirely under ignored data/v4.
        src = re.sub(r"\b(?:PROJECT_ROOT|ROOT)\s*/\s*['\"]configs/pilot_v1\.yaml['\"]", env_config, src)
        if key == "setup" and config.get("taxonomy_snapshot"):
            src = src.replace("TAXONOMY_ROOT = sorted((ROOT / 'data/v4/staging').glob('*/taxonomy'))[-1]",
                              f"TAXONOMY_ROOT = ROOT / 'data/v4/staging' / {config['taxonomy_snapshot']!r} / 'taxonomy'")
        if key == "setup":
            src = src.replace("'decision_workers': 1024", "'decision_workers': 16")
        if key == "setup" and "RUN_ID = 'human-parasites" in src:
            substitutions = {
                r"^RUN_ID = .*?$": f"RUN_ID = {config['run_id']!r}",
                r"^SELECTED_HOST_NAMES = .*?$": f"SELECTED_HOST_NAMES = {host_names!r}",
                r"^SELECTED_PATHOGEN_GROUPS = .*?$": f"SELECTED_PATHOGEN_GROUPS = {group_names!r}",
                r"^MAX_PMIDS_PER_PAIR = .*?$": f"MAX_PMIDS_PER_PAIR = {config['sample_per_pair']!r}",
                r"^OLDEST_PMIDS_PER_PAIR = .*?$": f"OLDEST_PMIDS_PER_PAIR = min(30, {config['sample_per_pair']!r})",
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
            src = src.replace("checkpoint_every=SCREEN_CHECKPOINT_EVERY)",
                              "checkpoint_every=SCREEN_CHECKPOINT_EVERY, stop_requested=lambda: Path(os.environ['PATHOGEN_ATLAS_PAUSE_FILE']).exists())")
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
    decision = re.search(r"Decision progress:\s*(\d+)/(\d+) this run.*?cached\s+([\d,]+).*?ETA ~([^|.]+)(?:\s*\| reported cost \$([0-9.]+); projected run cost \$([0-9.]+))?", line)
    search_plan = re.search(r"pairs=(\d+) already_checkpointed=(\d+) searches_to_run=(\d+)", line)
    setup = re.search(r"candidate pairs.?\s*(\d[\d,]*)", line, re.I)
    if not any((search, decision, search_plan, setup)):
        return
    if search:
        done, total = map(int, search.group(1, 2))
        metric.update({"completed": done, "total": total, "progress": round(done / max(total, 1) * 100, 1),
                       "eta_seconds": _parse_eta(search.group(3))})
    elif decision:
        done, total = map(int, decision.group(1, 2))
        cached = int(decision.group(3).replace(",", ""))
        metric.update({"completed": done + cached, "total": total + cached,
                       "progress": round((done + cached) / max(total + cached, 1) * 100, 1),
                       "cache_hits": cached, "eta_seconds": _parse_eta(decision.group(4))})
        if decision.group(5): metric["cost_usd"] = float(decision.group(5))
        if decision.group(6): metric["estimate_usd"] = float(decision.group(6))
    elif search_plan:
        total, cached, remaining = map(int, search_plan.groups())
        metric.update({"completed": cached, "total": total, "cache_hits": cached,
                       "progress": round(cached / max(total, 1) * 100, 1)})
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


def _phase_counts(phase, run_dir, run_id):
    import pandas as pd
    run = ROOT / "data/v4/runs" / run_id
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
        path = ROOT / "data/v4/cache/literature" / cache_key / "tables/pair_searches.parquet"
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
        path = ROOT / "data/v4/references/release/11_resolve_wikipedia_references.csv"
        return ("references", len(pd.read_csv(path)) if path.exists() else 0)
    if phase == "gbif":
        path = ROOT / "data/v4/references/release/12_resolve_gbif_references.csv"
        return ("references", len(pd.read_csv(path)) if path.exists() else 0)
    if phase == "build":
        payload = ROOT / "sites/pathogen-atlas-explorer/explorer-data.json.gz"
        return ("Explorer payloads", 1 if payload.exists() else 0)
    return ("items", None)


async def _execute_notebook(nb, phase, run_dir):
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError
    client = NotebookClient(nb, timeout=None, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}})
    manifest = read_json(run_dir / "manifest.json", {})
    manifest["status"] = "running"
    manifest["phases"][phase].update({"status": "running", "started_at": time.strftime("%Y-%m-%d %H:%M:%S")})
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
    # Save rendered outputs for audit/debugging, separate from source notebooks.
    nbformat.write(nb, run_dir / f"{phase}.executed.ipynb")
    manifest = read_json(run_dir / "manifest.json", {})
    count_label, count = _phase_counts(phase, run_dir, manifest["run_id"])
    manifest["phases"][phase].update({"status": "completed", "progress": 100,
                                      "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                                      "count_label": count_label})
    if count is not None:
        manifest["phases"][phase].update({"completed": count, "total": count})
    if phase in {"taxonomy", "prompts", "search", "estimates", "handoff", "profiles", "validation", "wikipedia", "gbif", "build"}:
        manifest["phases"][phase]["cost_usd"] = 0.0
    save_json(run_dir / "manifest.json", manifest)
    emit(run_dir, "phase", phase=phase, status="completed", message=f"Completed {manifest['phases'][phase]['label']}.")
    return True


def run(run_id):
    run_dir = STATE_ROOT / run_id
    manifest = read_json(run_dir / "manifest.json")
    if not manifest: raise RuntimeError(f"Unknown run: {run_id}")
    config = {**manifest["config"], "run_id": run_id}
    # Pause requests are removed by the explicit resume endpoint; new runs start without one.
    manifest["status"] = "running"
    save_json(run_dir / "manifest.json", manifest)
    emit(run_dir, "run", message="Pipeline worker started.")
    try:
        for phase, filename in PHASES:
            if phase == "wikipedia" and not config.get("wikipedia", True):
                status(run_dir, phase=phase, status="skipped", progress=100, cost_usd=0)
                emit(run_dir, "phase", phase=phase, status="skipped", message="Skipped by run configuration.")
                continue
            if phase == "gbif" and not config.get("gbif", True):
                status(run_dir, phase=phase, status="skipped", progress=100, cost_usd=0)
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
            env["PATHOGEN_ATLAS_RUN_CONFIG"] = str(run_dir / "dashboard-run.yaml")
            env["PATHOGEN_ATLAS_PAUSE_FILE"] = str(run_dir / "pause.request")
            # Notebook kernels inherit process environment; subprocess worker is isolated per run.
            complete = asyncio.run(_execute_notebook(nb, phase, run_dir))
            if not complete: return
        manifest = read_json(run_dir / "manifest.json", {})
        manifest["status"] = "completed"
        manifest["publish_ready"] = True
        manifest["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
        save_json(run_dir / "manifest.json", manifest)
        emit(run_dir, "run", message="Explorer payload and embedding phase completed; release is ready for review.")
    except Exception:
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
