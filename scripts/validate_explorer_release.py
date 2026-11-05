#!/usr/bin/env python3
"""Validate a local Explorer artifact or the live GitHub Pages release."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pathogen_atlas_v4.explorer_release import REQUIRED_FILES, validate_explorer_release  # noqa: E402


def _fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"Cache-Control": "no-cache", "Pragma": "no-cache"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def _validate_live(base_url: str, expected_version: str | None) -> dict[str, int]:
    base = base_url.rstrip("/") + "/"
    with tempfile.TemporaryDirectory(prefix="pathogen-atlas-live-check-") as folder:
        site = Path(folder)
        cache_bust = urllib.parse.urlencode({"release-check": os.environ.get("GITHUB_SHA", "latest")[:12]})
        def save(relative: str) -> bytes:
            separator = "&" if "?" in relative else "?"
            payload = _fetch(urllib.parse.urljoin(base, relative + separator + cache_bust))
            destination = site / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
            return payload

        save("index.html")
        for name in REQUIRED_FILES[1:]:
            save(name)
        metadata_path = site / "explorer-data.json.gz"
        metadata = json.loads(gzip.decompress(metadata_path.read_bytes()))
        for group in metadata.get("groups", []):
            save(f"taxa/{group}.json.gz")
        return validate_explorer_release(site, expected_version=expected_version)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--site-root", type=Path, help="Directory containing index.html and deployed payloads")
    source.add_argument("--url", help="Deployed site base URL, such as the GitHub Pages URL")
    parser.add_argument("--expected-version", help="Require this V4 version in the Explorer HTML")
    args = parser.parse_args()
    try:
        result = (_validate_live(args.url, args.expected_version) if args.url else
                  validate_explorer_release(args.site_root, expected_version=args.expected_version))
    except (OSError, ValueError, urllib.error.URLError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print("Explorer release checks passed: " + json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
