#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

PYTHON="${PYTHON:-python}"
if ! command -v "${PYTHON}" >/dev/null 2>&1; then
  PYTHON=python3
fi

SMOKE_PORT="${PATHOGEN_ATLAS_SMOKE_PORT:-${PATHOGEN_ATLAS_DASHBOARD_PORT:-8765}}"
FULL_PORT="${PATHOGEN_ATLAS_FULL_PORT:-8766}"
if [[ "${SMOKE_PORT}" == "${FULL_PORT}" ]]; then
  echo "Smoke and full dashboards need different ports." >&2
  exit 1
fi
SMOKE_URL="http://127.0.0.1:${SMOKE_PORT}"
FULL_URL="http://127.0.0.1:${FULL_PORT}"

stop_previous_dashboard() {
  local port="$1" pids pid command remaining
  if ! command -v lsof >/dev/null 2>&1; then
    echo "lsof is required to check for an existing dashboard on port ${port}." >&2
    return 1
  fi
  pids="$(lsof -nP -t -iTCP:"${port}" -sTCP:LISTEN 2>/dev/null || true)"
  for pid in ${pids}; do
    command="$(ps -p "${pid}" -o command= 2>/dev/null || true)"
    if [[ "${command}" != *"pathogen_atlas_v4.dashboard"* ]]; then
      echo "Port ${port} is occupied by another process; leaving it untouched." >&2
      return 1
    fi
    kill "${pid}" 2>/dev/null || true
  done
  for _ in {1..50}; do
    remaining="$(lsof -nP -t -iTCP:"${port}" -sTCP:LISTEN 2>/dev/null || true)"
    [[ -z "${remaining}" ]] && return 0
    sleep 0.1
  done
  echo "Previous dashboard on port ${port} did not stop." >&2
  return 1
}

# Re-running the launcher replaces its old smoke/full servers before binding
# the ports again, so stale dashboard processes cannot accumulate.
stop_previous_dashboard "${SMOKE_PORT}"
stop_previous_dashboard "${FULL_PORT}"

PATHOGEN_ATLAS_DASHBOARD_MODE=smoke_test PATHOGEN_ATLAS_DASHBOARD_PORT="${SMOKE_PORT}" \
  PATHOGEN_ATLAS_SMOKE_PORT="${SMOKE_PORT}" \
  "${PYTHON}" -m pathogen_atlas_v4.dashboard &
SMOKE_PID=$!
PATHOGEN_ATLAS_DASHBOARD_MODE=full PATHOGEN_ATLAS_DASHBOARD_PORT="${FULL_PORT}" \
  PATHOGEN_ATLAS_SMOKE_PORT="${SMOKE_PORT}" \
  "${PYTHON}" -m pathogen_atlas_v4.dashboard &
FULL_PID=$!

cleanup() {
  kill "${SMOKE_PID}" "${FULL_PID}" 2>/dev/null || true
  wait "${SMOKE_PID}" 2>/dev/null || true
  wait "${FULL_PID}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

wait_for_server() {
  local pid="$1" url="$2" label="$3"
  for _ in {1..100}; do
    if ! kill -0 "${pid}" 2>/dev/null; then
      wait "${pid}"
      return 1
    fi
    if curl --silent --fail "${url}/" >/dev/null; then
      echo "${label}: ${url}"
      return 0
    fi
    sleep 0.2
  done
  echo "${label} dashboard did not become available at ${url}" >&2
  return 1
}

wait_for_server "${SMOKE_PID}" "${SMOKE_URL}" "Smoke dashboard"
wait_for_server "${FULL_PID}" "${FULL_URL}" "Full dashboard"

if command -v open >/dev/null 2>&1; then
  open "${FULL_URL}"
elif command -v xdg-open >/dev/null 2>&1; then
  xdg-open "${FULL_URL}"
else
  echo "Open the full dashboard at ${FULL_URL}; its Smoke test link opens ${SMOKE_URL}."
fi

wait "${SMOKE_PID}"
