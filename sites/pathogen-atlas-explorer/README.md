# Pathogen Atlas V4 Explorer

The Explorer presents native V4 pathogen–host screening results. Profiles show the selected host, pair decision status, sampled and supporting PMID counts, publication estimates when available, taxonomy, external taxon references, and links to the screened PubMed records.

Build the payload from a completed run:

```sh
python sites/pathogen-atlas-explorer/scripts/build_v4_site_data.py \
  --run-dir data/v4/runs/pilot-v1 \
  --staging-root data/v4/staging/<taxonomy-snapshot>
```

The builder writes the Explorer payloads to the site source and `dist` directories. Publishing to the configured Pages branch deploys the checked-in site assets.
