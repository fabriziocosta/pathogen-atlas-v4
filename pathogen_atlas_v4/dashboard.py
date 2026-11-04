"""Local-only V4 pipeline dashboard server."""

from __future__ import annotations

import json
import gzip
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
STATE_ROOT = ROOT / "data/v4/dashboard/runs"
UI = Path(__file__).with_name("dashboard.html")
DASHBOARD_MODE = os.environ.get("PATHOGEN_ATLAS_DASHBOARD_MODE", "").strip()
SMOKE_DASHBOARD_PORT = os.environ.get("PATHOGEN_ATLAS_SMOKE_PORT", "8765")
DASHBOARD_VERSION = "1.0.13"
SERVER_INSTANCE_ID = uuid.uuid4().hex
PHASE_LOGS = {
    "setup": "00_configure_run", "taxonomy": "01_review_taxonomy_and_hosts",
    "prompts": "02_preview_prompts", "search": "03_acquire_pair_pubmed_records",
    "screen": "04_screen_pathogen_host_pairs", "estimates": "05_estimate_pathogen_publications",
    "handoff": "06_handoff_disease_evidence", "profiles": "07_resolve_feature_profiles",
    "scoring": "08_score_feature_profiles", "validation": "09_validate_release",
    "wikipedia": "11_resolve_wikipedia_references", "gbif": "12_resolve_gbif_references",
    "build": "10_inspect_embeddings",
}
_LIVE_LOG_STATE = {}
_SCREEN_LEDGER_STATE = {}
_LIVE_LOG_LOCK = threading.Lock()
_PUBLISH_MONITORS = set()
_PUBLISH_MONITOR_LOCK = threading.Lock()
PHASES = [
    ("setup", "Taxonomy, setup & host query expansion"), ("taxonomy", "Taxonomy snapshot review"),
    ("prompts", "Prompt preview"), ("search", "PubMed search & sampling"),
    ("screen", "Pathogen–host screening"), ("estimates", "Publication estimates"),
    ("handoff", "Evidence handoff"), ("profiles", "Feature profiles"),
    ("scoring", "Feature scoring"), ("validation", "Release validation"),
    ("wikipedia", "Wikipedia resolution"), ("gbif", "GBIF resolution"),
    ("build", "Explorer build & embeddings"),
]


def _json(path: Path, default):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def _save(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(
        f".{path.name}.{os.getpid()}.{threading.get_ident()}.{uuid.uuid4().hex}.tmp"
    )
    temp.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temp.replace(path)


@lru_cache(maxsize=1)
def _explorer_metadata():
    path = ROOT / "sites/pathogen-atlas-explorer/explorer-data.json.gz"
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, ValueError):
        return {}


@lru_cache(maxsize=5)
def _explorer_pathogens(group: str):
    path = ROOT / "sites/pathogen-atlas-explorer/taxa" / f"{group}.json.gz"
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            rows = json.load(stream)
        columns = _explorer_metadata().get("species_columns", [])
        id_col, name_col, rank_col = (columns.index(key) for key in ("organism_id", "accepted_name", "rank"))
        return tuple((str(row[id_col]), str(row[name_col])) for row in rows
                     if len(row) > max(id_col, name_col, rank_col) and str(row[rank_col]).casefold() == "species")
    except (OSError, ValueError, KeyError, StopIteration):
        return ()


def _taxonomy_options():
    snapshots = sorted((ROOT / "data/v4/staging").glob("*/taxonomy"))
    hosts, groups, counts = [], [], {}
    explorer_hosts = _explorer_metadata().get("hosts", [])
    explorer_labels = {str(row.get("host_id")): str(row.get("label")) for row in explorer_hosts}
    explorer_order = {str(row.get("host_id")): index for index, row in enumerate(explorer_hosts)}
    host_config = yaml.safe_load((ROOT / "configs/host_choices.yaml").read_text()) or {}
    host_mapping = host_config.get("host_name_to_id", host_config) if isinstance(host_config, dict) else {}
    seen_host_ids = set()
    for label, host_id in host_mapping.items():
        host_id = str(host_id)
        if host_id not in seen_host_ids:
            hosts.append({"id": host_id, "label": explorer_labels.get(host_id, str(label)),
                          "config_name": str(label)})
            seen_host_ids.add(host_id)
    hosts.sort(key=lambda host: explorer_order.get(host["id"], len(explorer_order)))

    group_config = yaml.safe_load((ROOT / "configs/pathogen_taxonomies.yaml").read_text()) or {}
    group_registry = group_config.get("pathogen_groups", {})
    if snapshots:
        tax = snapshots[-1]
        path_frame = pd.read_parquet(tax / "pathogens.parquet")
        counts = path_frame.groupby("pathogen_group").pathogen_id.nunique().to_dict()
    explorer_counts = _explorer_metadata().get("stats", {}).get("pathogen_species_by_group", {})
    if explorer_counts:
        counts = explorer_counts
    groups = [
        {"id": str(group), "label": (f"{config.get('label', str(group).title())} · {int(counts[group]):,} taxa"
                                      if group in counts else config.get("label", str(group).title())),
         "count": int(counts.get(group, 0))}
        for group, config in group_registry.items()
    ]
    return {"hosts": hosts, "groups": groups, "snapshot": snapshots[-1].parent.name if snapshots else None,
            "pathogen_counts": {str(k): int(v) for k, v in counts.items()}}


def _event(run_dir: Path, kind: str, **data):
    row = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "kind": kind, **data}
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(row, ensure_ascii=False) + "\n")
    manifest = _json(run_dir / "manifest.json", {})
    manifest["updated_at"] = row["time"]
    _save(run_dir / "manifest.json", manifest)


def _publish_state(run_id: str, status: str, *, message: str | None = None, **fields) -> None:
    run_dir = STATE_ROOT / run_id
    manifest_path = run_dir / "manifest.json"
    manifest = _json(manifest_path, None)
    if not manifest:
        return
    manifest.update(fields)
    manifest["publish_status"] = status
    manifest["publish_updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    _save(manifest_path, manifest)
    if message:
        _event(run_dir, "publish", message=message, publish_status=status)


def _query_pages_run(commit: str, workflow_run_id: str | None = None) -> dict | None:
    if not shutil.which("gh"):
        return None
    if workflow_run_id:
        command = ["gh", "run", "view", str(workflow_run_id), "--json",
                   "status,conclusion,url,headSha"]
    else:
        command = ["gh", "run", "list", "--workflow", "deploy-pages.yml", "--branch", "main",
                   "--commit", commit, "--limit", "5", "--json",
                   "databaseId,status,conclusion,url,headSha"]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode:
        return None
    try:
        rows = json.loads(result.stdout)
    except (TypeError, ValueError):
        return None
    if workflow_run_id:
        if not isinstance(rows, dict) or rows.get("headSha") != commit:
            return None
        return {**rows, "databaseId": workflow_run_id}
    if not isinstance(rows, list):
        return None
    return next((row for row in rows if row.get("headSha") == commit), None)


def _monitor_pages(run_id: str, commit: str, workflow_run_id: str | None = None) -> None:
    run_dir = STATE_ROOT / run_id
    last_status = None
    try:
        if not shutil.which("gh"):
            _publish_state(run_id, "unknown", message="Push succeeded, but GitHub CLI is unavailable to monitor Pages.",
                           publish_commit=commit)
            return
        while True:
            run = _query_pages_run(commit, workflow_run_id)
            if run:
                status = str(run.get("status") or "queued")
                conclusion = run.get("conclusion")
                url = str(run.get("url") or "")
                pages_status = ("deployed" if conclusion == "success" else "failed") if status == "completed" else (
                    "deploying" if status == "in_progress" else "queued")
                fields = {"publish_commit": commit,
                          "publish_workflow_run_id": str(run.get("databaseId") or workflow_run_id or ""),
                          "publish_workflow_url": url,
                          "publish_workflow_status": status,
                          "publish_workflow_conclusion": conclusion}
                if pages_status != last_status:
                    if pages_status == "deployed":
                        message = f"GitHub Pages deployment for commit {commit[:7]} completed successfully."
                        fields["published_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                    elif pages_status == "failed":
                        message = f"GitHub Pages deployment failed for commit {commit[:7]}. Use Retry deployment to rerun it."
                    elif pages_status == "deploying":
                        message = f"GitHub Pages deployment is running for commit {commit[:7]}."
                    else:
                        message = f"GitHub Pages deployment is queued for commit {commit[:7]}."
                    _publish_state(run_id, pages_status, message=message, **fields)
                    last_status = pages_status
                if status == "completed":
                    return
            time.sleep(10)
    finally:
        with _PUBLISH_MONITOR_LOCK:
            _PUBLISH_MONITORS.discard(run_id)


def _start_pages_monitor(run_id: str, commit: str, workflow_run_id: str | None = None) -> None:
    with _PUBLISH_MONITOR_LOCK:
        if run_id in _PUBLISH_MONITORS:
            return
        _PUBLISH_MONITORS.add(run_id)
    thread = threading.Thread(target=_monitor_pages, args=(run_id, commit, workflow_run_id),
                              name=f"pages-monitor-{run_id}", daemon=True)
    thread.start()


def _resume_pages_monitors() -> None:
    for folder in STATE_ROOT.iterdir():
        if not folder.is_dir():
            continue
        manifest = _json(folder / "manifest.json", {})
        if manifest.get("publish_status") in {"pushed", "queued", "deploying", "unknown"} and manifest.get("publish_commit"):
            _start_pages_monitor(folder.name, str(manifest["publish_commit"]),
                                 manifest.get("publish_workflow_run_id") or None)


def _sync_active_notebook_logs(run_dir: Path, manifest: dict) -> None:
    """Forward new notebook log lines into run events while a phase is active."""
    phase = next((key for key, value in (manifest.get("phases") or {}).items()
                  if value.get("status") == "running" and key in PHASE_LOGS), None)
    if not phase:
        return
    log_path = ROOT / "data/v4/logs" / run_dir.name / f"{PHASE_LOGS[phase]}.log"
    if not log_path.is_file():
        return
    key = (run_dir.name, phase)
    with _LIVE_LOG_LOCK:
        state = _LIVE_LOG_STATE.get(key)
        if state is None:
            seen = set()
            phase_metrics = (manifest.get("phases") or {}).get(phase, {})
            for field in ("cost_usd", "estimate_usd", "cost_tasks_known", "cost_tasks_total",
                          "cost_segment_base_usd", "cost_segment_base_known", "cost_segment_base_total"):
                phase_metrics.pop(field, None)
            _save(run_dir / "manifest.json", manifest)
            event_path = run_dir / "events.jsonl"
            if event_path.exists():
                from .dashboard_worker import _capture_metrics
                for raw in event_path.read_text(encoding="utf-8").splitlines():
                    try:
                        event = json.loads(raw)
                    except ValueError:
                        continue
                    if event.get("phase") == phase and event.get("kind") == "log":
                        message = event.get("message", "")
                        seen.add(message)
                        _capture_metrics(run_dir, phase, message)
            state = {"offset": 0, "seen": seen}
            _LIVE_LOG_STATE[key] = state
        size = log_path.stat().st_size
        if size < state["offset"]:
            state["offset"] = 0
        with log_path.open(encoding="utf-8") as stream:
            stream.seek(state["offset"])
            chunk = stream.read()
            state["offset"] = stream.tell()
        if not chunk:
            _sync_screen_ledger_progress(run_dir)
            return
        from .dashboard_worker import _capture_metrics
        for line in chunk.splitlines():
            message = line.strip()
            if not message or message in state["seen"]:
                continue
            state["seen"].add(message)
            _capture_metrics(run_dir, phase, message)
            timestamp = re.match(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}(?::\d{2})?)", message)
            _event(run_dir, "log", phase=phase, message=message[:3000],
                   time=timestamp.group(1) if timestamp else time.strftime("%Y-%m-%d %H:%M:%S"))
        _sync_screen_ledger_progress(run_dir)


def _sync_screen_ledger_progress(run_dir: Path) -> None:
    """Reflect checkpointed disease-screen decisions while the notebook cell runs."""
    phase_data = (_json(run_dir / "manifest.json", {}).get("phases") or {}).get("screen", {})
    if phase_data.get("status") != "running":
        return
    total = phase_data.get("total")
    if not isinstance(total, int) or total < 1:
        return
    ledger_path = ROOT / "data/v4/cache/decisions/pathogen_screen_ledger.parquet"
    if not ledger_path.is_file():
        return
    try:
        modified = ledger_path.stat().st_mtime_ns
    except OSError:
        return
    state = _SCREEN_LEDGER_STATE.setdefault(run_dir.name, {"modified": None, "completed": 0})
    if modified != state["modified"]:
        try:
            ledger = pd.read_parquet(ledger_path, columns=["run_id", "status"],
                                     filters=[("run_id", "==", run_dir.name)])
            state["completed"] = int(ledger.status.astype(str).eq("completed").sum())
            state["modified"] = modified
        except Exception:
            return
    completed = min(int(state["completed"]), total)
    manifest = _json(run_dir / "manifest.json", {})
    phase_data = (manifest.get("phases") or {}).get("screen", {})
    completed = max(completed, min(int(phase_data.get("completed") or 0), total))
    if completed <= int(phase_data.get("completed") or 0):
        return
    phase_data.update({"completed": completed, "total": total, "count_label": "screen tasks",
                       "progress": round(completed / total * 100, 1)})
    if phase_data.get("eta_seconds") is None and completed > 0 and phase_data.get("started_at"):
        try:
            started = time.mktime(time.strptime(phase_data["started_at"], "%Y-%m-%d %H:%M:%S"))
            elapsed = max(time.time() - started, 1.0)
            rate = completed / elapsed
            phase_data["eta_seconds"] = round(max(0, total - completed) / rate) if rate else None
        except (ValueError, OverflowError):
            pass
    _save(run_dir / "manifest.json", manifest)



def _worker_alive(run_dir: Path):
    pid = _json(run_dir / "worker.json", {}).get("pid")
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def _other_worker_active(except_run_id: str | None = None):
    for folder in STATE_ROOT.iterdir() if STATE_ROOT.exists() else []:
        if folder.is_dir() and folder.name != except_run_id and _worker_alive(folder):
            return True
    return False


class Handler(BaseHTTPRequestHandler):
    server_version = "PathogenAtlasDashboard/1"

    def log_message(self, fmt, *args):
        return

    def _body(self):
        length = int(self.headers.get("content-length", "0"))
        return json.loads(self.rfile.read(length) or b"{}")

    def _reply(self, data, status=200, content_type="application/json"):
        body = data if isinstance(data, bytes) else json.dumps(data).encode()
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/":
            html = UI.read_text(encoding="utf-8")
            mode_labels = {"smoke_test": "Smoke test", "full": "Full cohort"}
            html = html.replace("__DASHBOARD_VERSION__", DASHBOARD_VERSION)
            html = html.replace("__DASHBOARD_MODE__", DASHBOARD_MODE)
            html = html.replace("__DASHBOARD_MODE_LABEL__", mode_labels.get(DASHBOARD_MODE, "Choose run mode"))
            html = html.replace('href="http://127.0.0.1:8765/"',
                                f'href="http://127.0.0.1:{SMOKE_DASHBOARD_PORT}/"')
            return self._reply(html.encode("utf-8"), content_type="text/html; charset=utf-8")
        if path == "/logo.png":
            logo = ROOT / "sites/pathogen-atlas-explorer/assets/pathogen-atlas-v4-logo.png"
            return self._reply(logo.read_bytes(), content_type="image/png")
        if path == "/api/options":
            return self._reply(_taxonomy_options())
        if path == "/api/status":
            return self._reply({"instance_id": SERVER_INSTANCE_ID, "version": DASHBOARD_VERSION})
        if path == "/api/pathogens":
            query = parse_qs(urlparse(self.path).query)
            group = str((query.get("group") or [""])[0])
            search = str((query.get("q") or [""])[0]).strip()[:100]
            options = _taxonomy_options()
            if group not in options["pathogen_counts"]:
                return self._reply({"pathogens": []})
            pathogens = [(pathogen_id, name) for pathogen_id, name in _explorer_pathogens(group)
                         if not search or search.casefold() in name.casefold()]
            pathogens = sorted(pathogens, key=lambda item: item[1].casefold())[:40]
            return self._reply({"pathogens": [
                {"id": pathogen_id, "label": f"{name} · {pathogen_id}"}
                for pathogen_id, name in pathogens
            ]})
        if path == "/api/runs":
            items = []
            for folder in sorted(STATE_ROOT.glob("*"), reverse=True):
                if folder.is_dir():
                    items.append(_json(folder / "manifest.json", {}))
            return self._reply(items)
        if path.startswith("/api/run/"):
            run_id = path.rsplit("/", 1)[-1]
            run_dir = STATE_ROOT / run_id
            manifest = _json(run_dir / "manifest.json", None)
            if manifest is None:
                return self._reply({"error": "Run not found"}, 404)
            _sync_active_notebook_logs(run_dir, manifest)
            manifest = _json(run_dir / "manifest.json", manifest)
            events = []
            event_path = run_dir / "events.jsonl"
            if event_path.exists():
                for line in event_path.read_text(encoding="utf-8").splitlines()[-600:]:
                    try: events.append(json.loads(line))
                    except ValueError: pass
            manifest["events"] = events
            manifest["worker_alive"] = _worker_alive(run_dir)
            return self._reply(manifest)
        if path.startswith("/api/stream/"):
            run_id = path.rsplit("/", 1)[-1]
            run_dir = STATE_ROOT / run_id
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Connection", "keep-alive")
            self.end_headers()
            offset = 0
            try:
                while True:
                    live_manifest = _json(run_dir / "manifest.json", {})
                    _sync_active_notebook_logs(run_dir, live_manifest)
                    file = run_dir / "events.jsonl"
                    if file.exists():
                        with file.open("r", encoding="utf-8") as stream:
                            stream.seek(offset)
                            for line in stream:
                                self.wfile.write(f"data: {line.strip()}\n\n".encode())
                                self.wfile.flush()
                            offset = stream.tell()
                    time.sleep(1)
            except (BrokenPipeError, ConnectionResetError):
                return
        return self._reply({"error": "Not found"}, 404)

    def do_POST(self):
        path = urlparse(self.path).path
        if path == "/api/restart":
            self._reply({"ok": True, "instance_id": SERVER_INSTANCE_ID}, 202)
            self.wfile.flush()
            timer = threading.Timer(
                0.75,
                lambda: os.execv(
                    sys.executable,
                    [sys.executable, "-m", "pathogen_atlas_v4.dashboard"],
                ),
            )
            timer.daemon = True
            timer.start()
            return
        if path == "/api/runs":
            if _other_worker_active():
                return self._reply({"error": "Another V4 pipeline run is active. Pause or finish it before starting a new run."}, 409)
            data = self._body()
            hosts = list(dict.fromkeys(map(str, data.get("hosts", []))))
            groups = list(dict.fromkeys(map(str, data.get("groups", []))))
            mode = str(data.get("mode", "smoke_test"))
            smoke_pathogen_id = str(data.get("smoke_pathogen_id", "")).strip()
            sample = data.get("sample", 10 if mode == "smoke_test" else 200)
            if isinstance(sample, str) and sample.isdecimal():
                sample = int(sample)
            if mode not in {"smoke_test", "full"}:
                return self._reply({"error": "Mode must be 'smoke_test' or 'full'."}, 400)
            if DASHBOARD_MODE and mode != DASHBOARD_MODE:
                return self._reply({"error": f"This address is reserved for {DASHBOARD_MODE.replace('_', ' ')} runs."}, 400)
            if not hosts or not groups or isinstance(sample, bool) or not isinstance(sample, int) or sample < 1:
                return self._reply({"error": "Choose at least one host and pathogen group; sample must be a positive integer."}, 400)
            maximum_sample = 10 if mode == "smoke_test" else 200
            if sample > maximum_sample:
                return self._reply({"error": f"{mode} allows at most {maximum_sample} PMIDs per pair."}, 400)
            if mode == "smoke_test" and (len(hosts) != 1 or len(groups) != 1):
                return self._reply({"error": "Smoke test requires exactly one host and one pathogen group."}, 400)
            options = _taxonomy_options()
            valid_hosts = {h["id"] for h in options["hosts"]}
            valid_groups = {group["id"] for group in options["groups"]}
            if not set(hosts) <= valid_hosts or not set(groups) <= valid_groups:
                return self._reply({"error": "Selection contains an unknown host or pathogen group."}, 400)
            if mode == "smoke_test":
                if not smoke_pathogen_id:
                    return self._reply({"error": "Choose one pathogen species for the smoke test."}, 400)
                if smoke_pathogen_id not in {pathogen_id for pathogen_id, _ in _explorer_pathogens(groups[0])}:
                    return self._reply({"error": "The smoke pathogen must be a species in the selected pathogen group."}, 400)
            run_id = f"dashboard-{time.strftime('%Y%m%d-%H%M%S')}-{uuid.uuid4().hex[:6]}"
            run_dir = STATE_ROOT / run_id
            run_dir.mkdir(parents=True, exist_ok=False)
            manifest = {"run_id": run_id, "status": "queued", "config": {"hosts": hosts, "groups": groups,
                "host_names": [h.get("config_name", h["label"]) for h in options["hosts"] if h["id"] in hosts],
                "mode": mode, "smoke_pathogen_id": smoke_pathogen_id,
                "sample_per_pair": sample, "wikipedia": bool(data.get("wikipedia", True)),
                "gbif": bool(data.get("gbif", True)), "taxonomy_snapshot": options["snapshot"]},
                "phases": {key: {"label": label, "status": "pending", "progress": 0, "completed": 0,
                                  "total": None, "eta_seconds": None, "cache_hits": 0, "cost_usd": None,
                                  "estimate_usd": None} for key, label in PHASES}, "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "publish_ready": False}
            _save(run_dir / "manifest.json", manifest)
            _event(run_dir, "run", message="Run configuration saved.")
            self._launch(run_id)
            return self._reply({"run_id": run_id}, 201)
        pieces = path.strip("/").split("/")
        if len(pieces) >= 4 and pieces[:2] == ["api", "run"]:
            run_id, action = pieces[2], pieces[3]
            run_dir = STATE_ROOT / run_id
            manifest = _json(run_dir / "manifest.json", None)
            if not manifest:
                return self._reply({"error": "Run not found"}, 404)
            if action == "pause":
                (run_dir / "pause.request").touch()
                manifest["status"] = "pausing"
                _save(run_dir / "manifest.json", manifest)
                _event(run_dir, "control", message="Pause requested; current provider batch will checkpoint.")
                return self._reply({"ok": True})
            if action == "resume":
                if _worker_alive(run_dir):
                    return self._reply({"error": "Worker is already active."}, 409)
                if _other_worker_active(run_id):
                    return self._reply({"error": "Another V4 pipeline run is active. Pause or finish it before resuming this run."}, 409)
                (run_dir / "pause.request").unlink(missing_ok=True)
                self._launch(run_id)
                return self._reply({"ok": True})
            if action == "publish-preview":
                return self._reply(self._publish_preview(run_dir, manifest))
            if action == "publish-retry":
                return self._reply(self._retry_publish(run_dir, manifest))
            if action == "publish":
                data = self._body()
                if data.get("confirm") is not True:
                    return self._reply({"error": "Explicit confirmation required."}, 400)
                return self._reply(self._publish(run_dir, manifest))
        return self._reply({"error": "Not found"}, 404)

    def _launch(self, run_id):
        run_dir = STATE_ROOT / run_id
        command = [sys.executable, "-m", "pathogen_atlas_v4.dashboard_worker", "--run-id", run_id]
        env = dict(os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""))
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   start_new_session=True)
        _save(run_dir / "worker.json", {"pid": process.pid})

    def _publish_preview(self, run_dir, manifest):
        if manifest.get("status") != "completed" or not manifest.get("publish_ready"):
            return {"error": "Publishing is enabled after Explorer validation and packaging complete."}
        paths = _release_paths()
        return {"run_id": manifest["run_id"], "version": _next_version(), "files": paths,
                "summary": "Additive V4 Explorer update. Only listed Explorer payloads and release metadata will be staged."}

    def _publish(self, run_dir, manifest):
        preview = self._publish_preview(run_dir, manifest)
        if "error" in preview:
            return preview
        existing_status = manifest.get("publish_status")
        if existing_status in {"queued", "deploying", "deployed", "pushed", "unknown"} and manifest.get("publish_commit"):
            return {"error": f"This release is already {existing_status}."}
        if existing_status == "failed":
            return {"error": "The release was pushed, but its deployment failed. Use Retry deployment."}
        branch = subprocess.run(["git", "branch", "--show-current"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
        if branch != "main":
            return {"error": f"Publishing requires the Pages deployment branch main; current branch is {branch}."}
        upstream = subprocess.run(["git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{upstream}"],
                                  cwd=ROOT, text=True, capture_output=True, check=False)
        remote_branch = upstream.stdout.strip()
        if upstream.returncode != 0 or "/" not in remote_branch:
            return {"error": "The current branch has no configured upstream remote."}
        remote, remote_main = remote_branch.split("/", 1)
        if remote_main != "main":
            return {"error": f"Configured upstream {remote_branch} is not the Pages branch main."}
        fetched = subprocess.run(["git", "fetch", remote], cwd=ROOT, text=True, capture_output=True)
        if fetched.returncode:
            return {"error": fetched.stderr[-3000:] or fetched.stdout[-3000:]}
        counts = subprocess.run(["git", "rev-list", "--left-right", "--count", f"HEAD...{remote}/main"],
                                cwd=ROOT, text=True, capture_output=True, check=True).stdout.split()
        if len(counts) == 2 and int(counts[1]) > 0:
            return {"error": f"{remote}/main has {counts[1]} commit(s) missing locally. Sync the checkout before publishing."}
        # Only after the release target and remote are verified do we bump release metadata.
        result = subprocess.run([sys.executable, "-m", "pathogen_atlas_v4.dashboard_worker", "--bump-release",
                                 manifest["run_id"]], cwd=ROOT, text=True, capture_output=True)
        if result.returncode:
            return {"error": result.stderr[-3000:] or result.stdout[-3000:]}
        paths = _release_paths() + ["CHANGELOG.md", "pyproject.toml", "pathogen_atlas_v4/__init__.py",
                                    "sites/pathogen-atlas-explorer/dist/index.html"]
        # A path-limited commit excludes any existing staged or working-tree changes outside the allowlist.
        version = _version()
        commit = subprocess.run(["git", "commit", "--only", "-m", f"Release Pathogen Atlas V4 {version}", "--", *paths], cwd=ROOT,
                                text=True, capture_output=True)
        if commit.returncode:
            return {"error": commit.stderr[-3000:] or commit.stdout[-3000:]}
        push = subprocess.run(["git", "push", remote, "main"], cwd=ROOT, text=True, capture_output=True)
        if push.returncode:
            return {"error": push.stderr[-3000:] or push.stdout[-3000:], "commit": commit.stdout.strip()}
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
        _publish_state(manifest["run_id"], "queued", message=f"V4 {version} pushed; GitHub Pages deployment is queued.",
                       publish_version=version, publish_commit=sha,
                       publish_workflow_run_id="", publish_workflow_url="",
                       publish_workflow_status="queued", publish_workflow_conclusion=None)
        _start_pages_monitor(manifest["run_id"], sha)
        return {"status": "queued", "version": version, "commit": sha}

    def _retry_publish(self, run_dir, manifest):
        if manifest.get("status") != "completed" or not manifest.get("publish_ready"):
            return {"error": "A completed full run is required before deployment can be retried."}
        if manifest.get("publish_status") != "failed":
            return {"error": "There is no failed Pages deployment to retry."}
        workflow_run_id = str(manifest.get("publish_workflow_run_id") or "")
        if not workflow_run_id or not shutil.which("gh"):
            return {"error": "The failed Pages workflow could not be identified for retry."}
        result = subprocess.run(["gh", "run", "rerun", workflow_run_id, "--failed"],
                                cwd=ROOT, text=True, capture_output=True)
        if result.returncode:
            return {"error": result.stderr[-3000:] or result.stdout[-3000:]}
        _publish_state(manifest["run_id"], "queued", message="Pages deployment retry queued.",
                       publish_workflow_status="queued", publish_workflow_conclusion=None,
                       publish_error="")
        _start_pages_monitor(manifest["run_id"], str(manifest.get("publish_commit") or ""))
        return {"status": "queued", "version": manifest.get("publish_version"),
                "commit": manifest.get("publish_commit"), "workflow": workflow_run_id}


def _release_paths():
    site = ROOT / "sites/pathogen-atlas-explorer"
    candidates = [site / "index.html", site / "explorer-data.json.gz", site / "evidence.json.gz", site / "feature-catalog.json"]
    candidates += sorted((site / "dist/taxa").glob("*.json.gz")) if (site / "dist/taxa").exists() else []
    return [str(path.relative_to(ROOT)) for path in candidates if path.is_file()]


def _version():
    import re
    text = (ROOT / "pyproject.toml").read_text()
    return re.search(r'^version = "([^"]+)"', text, re.M).group(1)


def _next_version():
    import re
    match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", _version())
    if not match: return _version()
    return f"{match[1]}.{match[2]}.{int(match[3]) + 1}"


def main():
    os.chdir(ROOT)
    STATE_ROOT.mkdir(parents=True, exist_ok=True)
    _resume_pages_monitors()
    server = ThreadingHTTPServer(("127.0.0.1", int(os.environ.get("PATHOGEN_ATLAS_DASHBOARD_PORT", "8765"))), Handler)
    mode_label = {"smoke_test": "Smoke", "full": "Full"}.get(DASHBOARD_MODE, "Pipeline")
    print(f"Pathogen Atlas V4 {mode_label} dashboard: http://{server.server_address[0]}:{server.server_address[1]}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
