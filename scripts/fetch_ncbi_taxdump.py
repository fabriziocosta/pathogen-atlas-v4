"""Fetch the current NCBI taxonomy dump and record its content hash."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import requests

URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/source_taxonomy/ncbi"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    archive = args.output / "taxdump.tar.gz"
    digest = hashlib.sha256()
    with requests.get(URL, stream=True, timeout=(30, 180)) as response:
        response.raise_for_status()
        with archive.open("wb") as stream:
            for block in response.iter_content(8 * 1024 * 1024):
                if block:
                    digest.update(block)
                    stream.write(block)
    wanted = {"names.dmp", "nodes.dmp"}
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            name = Path(member.name).name
            if name not in wanted:
                continue
            if not member.isfile():
                raise ValueError(f"Unexpected non-file member {member.name}")
            extracted = bundle.extractfile(member)
            if extracted is None:
                raise ValueError(f"Could not read {member.name}")
            (args.output / name).write_bytes(extracted.read())
    if not all((args.output / name).is_file() for name in wanted):
        raise RuntimeError("NCBI taxdump did not contain names.dmp and nodes.dmp")
    manifest = {"source": "NCBI Taxonomy taxdump", "url": URL,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "archive_sha256": digest.hexdigest(), "archive_bytes": archive.stat().st_size}
    (args.output / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
