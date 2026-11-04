# Pathogen Atlas V4 Explorer

The Explorer catalogs all species in the pinned pathogen taxonomy across the available human, livestock, and crop hosts. Host-specific pair profiles are shown where screening data is available. An unassessed pair is labeled “not screened”; it is not treated as a negative result.

Build the atlas payload from the pinned taxonomy and a completed host-pair run:

```sh
python sites/pathogen-atlas-explorer/scripts/build_v4_site_data.py \
  --run-dir data/v4/runs/pilot-v1 \
  --staging-root data/v4/staging/<taxonomy-snapshot>
```

The builder combines the full staged taxonomy, configured host list, saved human evidence profiles in `source/human_pair_data.json.gz`, and pair results from the selected run. It writes the Explorer payloads to the site source and `dist` directories. Publishing to the configured Pages branch deploys the checked-in site assets.

Each build overlays the selected run onto the prior Explorer taxonomy shards. The current run replaces records for pairs it contains; pair profiles and their detailed Jev evidence from other hosts and cohorts remain available. Historical ledgers backfill detailed assessments if an older payload preserved the profile but dropped its evidence rows. The builder checks evidence counts for refreshed profiles and profiles restored from run ledgers. Eligible host × pathogen-group feature cohorts receive reproducible IMSTE and UMAP coordinate maps, with both coordinate sets and fit provenance stored in the payload.
