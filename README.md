# PhenoDigm2

PhenoDigm2 builds phenotype annotation databases and computes PhenoDigm scores
that compare animal models with diseases. A completed build can be distributed
as SQLite data and as a Parquet document bundle for the PhenoDigm Solr index.

## Installation

PhenoDigm2 is a Python package built with
[`uv`](https://docs.astral.sh/uv/) and requires Python 3.12 or later.

```bash
git clone [REPO_URL]
cd PhenoDigm2
uv sync
uv run phenodigm --help
```

`phenodigm` is the supported command-line entry point. The legacy
`python3 phenodigm2.py ...` form is no longer used in the main documentation.

Owltools is no longer supported. Current builds use downloaded
Phenio/Semsimian mappings through the `ontology_mapping` action.

## Build overview

Use a release-specific directory such as `vTODAY`:

```bash
uv run phenodigm --init vTODAY

# Check the latest version-specific Phenio/Semsimian Zenodo record first.
uvx zenodo_get -w - 18474575
# Update the record ID in vTODAY/resources/dependencies.json.

uv run phenodigm download --db vTODAY
uv run phenodigm build --db vTODAY
uv run phenodigm ontology_mapping --db vTODAY
uv run phenodigm score --fast --db vTODAY
uv run phenodigm index --db vTODAY
uv run phenodigm parquet --db vTODAY
```

The Parquet action creates `vTODAY/output/parquet` with the document types
required by EBI to build the PhenoDigm Solr index. A local Solr core remains an
optional compatibility workflow.

See the following guides for complete instructions:

- [BUILD.md](BUILD.md) — package setup and the complete build sequence.
- [IMPC_RELEASE.md](IMPC_RELEASE.md) — IMPC release checklist and validation.
- [PARQUET.md](PARQUET.md) — Parquet schemas and bundle layout.
- [SOLR.md](SOLR.md) — optional local Solr configuration.
- [EXAMPLES.md](EXAMPLES.md) — database query and export examples.

## Useful commands

Inspect a completed database:

```bash
uv run phenodigm status --db vTODAY
```

Use a non-default ontology-mapping IC threshold:

```bash
uv run phenodigm ontology_mapping \
  --ontology_mapping_min_ic <IC> \
  --db vTODAY
```

Prepare and build an optional local Solr core:

```bash
uv run phenodigm solr-prepare --db vTODAY
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml up -d
uv run phenodigm solr \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml down
```

`solr-prepare` is idempotent and removes the previous manual Compose-file copy.
See [SOLR.md](SOLR.md) for lifecycle and customization details.

Run the test suite:

```bash
uv run pytest
```
