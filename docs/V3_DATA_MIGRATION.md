# V3 data migration runbook

## Purpose and migration boundary

This runbook describes how to carry the existing local V3 dataset into V4 without changing V3 or treating old human-disease scores as newly evaluated host-specific results.

V3’s data are local-only and excluded from Git. The GitHub repository contains the configuration and migration description, but not the several-gigabyte working dataset. The transfer must therefore run on a machine or mounted volume that can read the V3 data directory. Keep the data themselves out of Git and GitHub; track the migration code, schema, small migration summary, and this runbook in V4.

Use the V3 configuration and docs as the source of truth for input paths and semantics:

- V3 config: https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/configs/default.yaml
- V3 data model: https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/docs/data_model.md
- V3 workflow: https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/docs/workflow.md
- V3 V2-snapshot provenance: https://github.com/fabriziocosta/pathogen-atlas-v3/blob/main/configs/migrations/v2_snapshot.yaml

The V3 paths below are confirmed by its current configuration and workflow. V4 target tables are the proposed migration contract in this document; implement and version those schemas before running the converter.

## V3 inputs to capture

Copy the complete V3 data root, not just its release tables. Preserve all files under these paths when present:

| V3 path relative to repository root | Role in migration |
| --- | --- |
| data/01_taxonomy/release/ | Normalized pathogen taxonomy, aliases, and virus variants |
| data/02_pubmed/cache/ | Reusable query and PMID caches |
| data/02_pubmed/release/ | Abstract and pathogen-to-PMID release tables |
| data/02_pubmed/domains/ | Per-domain PubMed outputs |
| data/03_decisions/parsed/decision_ledger.parquet | Parsed V3 gate and feature decisions |
| data/03_decisions/raw_responses/jev.jsonl | Original model responses |
| data/03_decisions/legacy_v2/ | V2 decisions retained as legacy inputs, if present |
| data/04_release/ | V3-derived release tables and summaries |
| data/05_embeddings/ | V3-derived coordinates and embedding artifacts |
| data/migration/version2_snapshot/ | V2 data snapshot carried forward by V3 |
| data/source_taxonomy/ | Pinned source taxonomy inputs used by V3 |

Also record the exact V3 Git commit, the V3 configuration and question/profile versions, the taxonomy source releases, and the migration date. Do not copy .env files, API keys, local credentials, virtual environments, caches outside the V3 data root, or unrelated files.

## Step 1: freeze and identify the source

1. Stop V3 acquisition, decision, release, and embedding jobs. Do not migrate while files are being written.
2. Set V3_ROOT to the checked-out V3 repository and V4_ROOT to the V4 checkout on a machine that can read both.
3. Confirm that V3_ROOT/data exists and note the V3 commit. Data are ignored by Git, so a clean Git status does not prove that the data tree is unchanged; the file manifest below is the data snapshot identity.
4. Create the V4 snapshot directories. The snapshot is immutable input; conversions go to a separate staging directory.

~~~bash
V3_ROOT=/absolute/path/to/pathogen-atlas-v3
V4_ROOT=/absolute/path/to/pathogen-atlas-v4
V3_DATA="$V3_ROOT/data"
V4_SNAPSHOT="$V4_ROOT/data/migration/v3_snapshot"

test -d "$V3_DATA"
mkdir -p "$V4_SNAPSHOT/data"
git -C "$V3_ROOT" rev-parse HEAD > "$V4_SNAPSHOT/source_git_commit.txt"
~~~

## Step 2: make a SHA-256 manifest before transfer

Create a manifest outside the V3 data root. It records every regular file’s relative path, byte count, and SHA-256, as well as the source Git commit. Reject symlinks for manual review so the copy cannot silently read outside the declared data root.

~~~bash
python3 - "$V3_ROOT" "$V3_DATA" "$V4_SNAPSHOT/transfer_manifest.json" <<'PY'
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

repo = Path(sys.argv[1]).resolve()
data_root = Path(sys.argv[2]).resolve()
manifest_path = Path(sys.argv[3])

symlinks = [p for p in data_root.rglob("*") if p.is_symlink()]
if symlinks:
    raise SystemExit("Review symlinks before migration: " + ", ".join(map(str, symlinks)))

files = []
for path in sorted(p for p in data_root.rglob("*") if p.is_file()):
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    files.append({
        "path": path.relative_to(data_root).as_posix(),
        "size_bytes": size,
        "sha256": digest.hexdigest(),
    })

commit = subprocess.check_output(
    ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
).strip()

manifest = {
    "migration_id": "v3-to-v4",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "source_repository": "https://github.com/fabriziocosta/pathogen-atlas-v3",
    "source_commit": commit,
    "source_data_root": "data/",
    "file_count": len(files),
    "total_bytes": sum(item["size_bytes"] for item in files),
    "files": files,
}
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
print(f"Manifest: {manifest_path}")
print(f"Files: {manifest['file_count']}; bytes: {manifest['total_bytes']}")
PY
~~~

Retain transfer_manifest.json beside the snapshot. It is migration provenance, not a replacement for the V3 source files. Also record the selected taxonomy, feature-profile, gate-question, and feature-question versions in the V4 migration summary.

## Step 3: copy the snapshot without using Git

For source and destination paths available on the same machine or mounted filesystem:

~~~bash
rsync -a --partial --human-readable --progress \
  "$V3_DATA/" \
  "$V4_SNAPSHOT/data/"
~~~

For a remote V4 data volume, use rsync over the institution’s approved SSH or storage mount, copying the same source tree and manifest. Resume an interrupted transfer with the same command. Do not add --delete. Do not upload this data tree to the GitHub repository.

Ensure V4 ignores local data before continuing. Confirm that a test path under data/ is ignored by Git, and check git status to make sure no data files are staged.

## Step 4: verify the copied snapshot

Run this verifier from any directory. It checks missing files, extra files, size changes, and SHA-256 mismatches. Stop on any mismatch; do not convert a partial snapshot.

~~~bash
python3 - "$V4_SNAPSHOT/transfer_manifest.json" "$V4_SNAPSHOT/data" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
root = Path(sys.argv[2])
expected = {row["path"]: row for row in manifest["files"]}
actual = {
    p.relative_to(root).as_posix(): p
    for p in root.rglob("*")
    if p.is_file()
}

errors = []
for rel, row in expected.items():
    path = actual.get(rel)
    if path is None:
        errors.append(f"missing: {rel}")
        continue
    if path.stat().st_size != row["size_bytes"]:
        errors.append(f"size mismatch: {rel}")
        continue
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != row["sha256"]:
        errors.append(f"SHA-256 mismatch: {rel}")

for rel in sorted(set(actual) - set(expected)):
    errors.append(f"unexpected file: {rel}")

if errors:
    raise SystemExit("\n".join(errors))
print(f"Verified {len(expected)} files, {manifest['total_bytes']} bytes.")
print(f"Source commit: {manifest['source_commit']}")
PY
~~~

Do not edit, rename, or normalize files inside data/migration/v3_snapshot/data after verification. Make the snapshot read-only for the migration job where practical. If the converter needs scratch data, it must write to a separate staging directory.

## Step 5: transform into the V4 staging tables

Implement one deterministic converter, planned at scripts/migration/import_v3_data.py. It must accept the verified snapshot path, a versioned V4 migration configuration, an output staging path, and a dry-run option. It must never modify its source snapshot. Give each conversion a run ID derived from the source commit and migration-config version. Write a conversion report with source and output row counts, unmapped records, duplicate keys, rejected rows, and warnings.

Proposed V4 output layout:

~~~text
data/v4/staging/<migration_run_id>/
  taxonomy/pathogens.parquet
  taxonomy/pathogen_aliases.parquet
  taxonomy/pathogen_variants.parquet
  taxonomy/hosts.parquet
  literature/articles.parquet
  literature/pathogen_article_matches.parquet
  evidence/pathogen_host_evidence.parquet
  decisions/legacy_v3_decisions.parquet
  decisions/legacy_v3_feature_assessments.parquet
  provenance/source_records.parquet
  migration_report.json
~~~

Apply these mappings:

| V3 source | V4 staging destination | Conversion rule |
| --- | --- | --- |
| 01_taxonomy/release/species.parquet | taxonomy/pathogens.parquet | Preserve the V3 organism ID, accepted name, organism group, taxon ID, parent, rank, taxonomy source, and release. Add source_version=v3. Do not create new pathogen IDs from names. |
| 01_taxonomy/release/aliases.parquet | taxonomy/pathogen_aliases.parquet | Preserve every alias, mapping type, source, and unresolved/ambiguous status. Link by the V3 organism ID, not by name matching. |
| 01_taxonomy/release/virus_variants.parquet, if present | taxonomy/pathogen_variants.parquet | Preserve MSL species linkage, isolate/designation/accession fields, host-source fields, and original source row. |
| 02_pubmed/release/abstracts.parquet | literature/articles.parquet | Deduplicate by PMID only after verifying that duplicate rows agree on title/abstract metadata. Preserve the selected source row and report disagreements. |
| 02_pubmed/release/species_abstracts.parquet | literature/pathogen_article_matches.parquet | Preserve pathogen ID, PMID, query hash, sample/fetch flags, selection reason, and source counts. These are pathogen-to-literature matches, not automatically host-pair evidence. |
| 02_pubmed/cache/, 02_pubmed/domains/ | provenance/source_records.parquet plus archived snapshot | Keep the original files unchanged. Index cache records by their query hash and PMID where available; do not turn cache presence into a positive disease decision. |
| 03_decisions/parsed/decision_ledger.parquet and legacy_v2/decisions.parquet, if present | decisions/legacy_v3_decisions.parquet and decisions/legacy_v3_feature_assessments.parquet | Preserve raw probabilities, labels, model, run ID, input hash, question/profile versions, source decision, and threshold fields. Label all imported rows as V3 legacy provenance. |
| 03_decisions/raw_responses/jev.jsonl | archived snapshot and provenance/source_records.parquet | Preserve byte-for-byte. Do not rewrite raw provider responses during migration. |
| 04_release/ | archived snapshot; use only as a comparison source | Rebuild V4 releases from V4 staging decisions; do not publish V3 release tables as V4 outputs. |
| 05_embeddings/ | archived snapshot; use only as a comparison source | Refit V4 embeddings from V4-resolved features. V3 coordinates use V3’s feature profiles and are not V4 coordinates. |
| migration/version2_snapshot/ and source_taxonomy/ | archived snapshot and provenance/source_records.parquet | Preserve source files and source checksums. Use them to explain V3 provenance; do not silently replace them with newer downloads. |

### Human host mapping and legacy decisions

V3’s disease gate and feature wording use a human-disease context. Resolve that context through an explicit V4 host mapping entry for human; do not infer the host from a pathogen name or from a PubMed search match.

For each V3 gate row, retain the original V3 result in legacy_v3_decisions.parquet. Where the V3 release accepted the row as human disease evidence, create a corresponding row in pathogen_host_evidence.parquet with the mapped human host ID, PMID, original V3 probability, original threshold, question version, and provenance status legacy_v3_accepted. This is imported historical evidence, not a V4 model decision. Keep negative, uncertain, missing, and failed V3 outcomes in the legacy table; do not turn them into positive V4 pair records.

V3 did not store model-extracted evidence offsets. The importer therefore retains the full source abstract as `evidence_span` for accepted legacy rows and marks `evidence_span_type=full_abstract_legacy`. This coarse historical passage must not be presented as a V4-extracted quotation. Native V4 supported rows require an explicit nonempty evidence span.

Keep V3 feature scores in legacy_v3_feature_assessments.parquet with their original human-scoped question and profile versions. Reuse a score as a V4 cache seed only after an explicit feature-equivalence mapping confirms that the wording, criteria, evidence context, host scope, and pathogen scope are identical. Otherwise, retain it as historical evidence and rescore with V4 prompts.

V3’s pathogenic-interest estimates and embeddings are derived outputs. Preserve them for comparison, but recompute V4 metrics and coordinates from V4’s pair-aware evidence and resolved feature profiles.

## Step 6: validate before promotion

The converter must be deterministic: the same verified snapshot, migration configuration, and converter version must produce identical sorted tables and the same migration report, apart from explicitly excluded run timestamps.

Before promoting staging to a V4 release, verify:

- Every file in the snapshot manifest passes size and SHA-256 checks.
- Every V3 taxonomy row is either mapped to a V4 pathogen row or listed in an explicit rejection/unmapped report; no row disappears silently.
- V3 organism IDs remain stable in the imported pathogen table; aliases and variants link to valid pathogen IDs.
- Every imported article match references an existing article and pathogen.
- Every migrated human pair evidence row references an explicit human host row and an existing pathogen and article.
- Every legacy decision retains its source question/profile/model/run identifiers and original probability.
- No legacy decision is labelled as a native V4 decision.
- V3 release and embedding outputs are present only in the immutable snapshot or comparison area, not substituted for regenerated V4 releases.
- Counts, duplicate keys, unresolved IDs, and rejected rows are recorded in migration_report.json and reviewed before release.

Run on a small copy first. Promote the complete verified staging directory atomically to the V4 release directory only after the report is accepted. Keep the original snapshot and report so the migration can be reproduced.

## Known count caveat

V3’s configs/migrations/v2_snapshot.yaml records V2-to-V3 source counts, including 22,804 species, 110,213 abstracts, 2,076,387 species–PMID matches, and 2,134,734 legacy decisions. These are V2 snapshot/import counts, not guaranteed counts of the final normalized V3 releases. Use the file-level V3 manifest and actual V3 table counts as the authoritative transfer checks.

## Migration acceptance criteria

- V3’s data tree is preserved as a checksum-verified, immutable input snapshot.
- V4 outputs are written only to a separate staging/release tree.
- Each imported row has a traceable V3 source file and source identifier.
- Human evidence imported from V3 remains explicitly legacy evidence; V4 host-pair decisions are generated only under V4’s versioned questions.
- The conversion report accounts for every input row, including duplicate, ambiguous, unresolved, rejected, and migrated rows.
- The entire conversion can be rerun without modifying V3 or changing the verified snapshot.

## V4 implementation entry points

The V4 package provides versioned contracts in `pathogen_atlas_v4/schemas.py`, run and pair validation in `runs.py`, host-scoped cache identity in `decision_identity.py`, host-aware query construction in `pubmed.py`, NCBI host parsing in `host_taxonomy/ncbi.py`, profile resolution in `features.py`, a resumable provider-agnostic ledger runner in `decision_runner.py`, and positive evidence/release checks in `evidence.py` and `release.py`.

Install the project and run its tests with `uv venv .venv`, `uv pip install --python .venv/bin/python '.[dev]'`, and `.venv/bin/python -m pytest -q`. Run the full import from the V4 root:

~~~bash
.venv/bin/pathogen-atlas-v4-import-v3 \
  --source-data ../version3/data \
  --snapshot-dir data/migration/v3_snapshot \
  --output data/v4/staging/v3-to-v4-20261102-ncbi-f278697a8ba6 \
  --snapshot \
  --host-taxdump data/source_taxonomy/ncbi \
  --host-taxonomy-release ncbi-taxdump-2026-11-02-sha256-f278697a8ba6 \
  --config-version v3-to-v4-import-v3-ncbi-f278697a8ba6
~~~

The import CLI snapshots and verifies the data tree before writing separate staging tables. Staging is not a release until its migration report and referential checks have been reviewed.
