# PhenoDigm2

PhenoDigm2 builds phenotype annotation databases and computes PhenoDigm scores
that compare animal models with diseases. A completed build can be distributed
as SQLite data and as a Parquet document bundle for the PhenoDigm Solr index.

PhenoDigm2 is used **primarily as a Python library** (the `PhenoDigm` class); an
**equivalent command-line interface** is also available. Both drive the same
pipeline stages, so pick whichever fits your workflow — a script/orchestrator
(Airflow) or a shell.

## Installation

PhenoDigm2 is a Python package built with
[`uv`](https://docs.astral.sh/uv/) and requires Python 3.12 or later.

```bash
git clone [REPO_URL]
cd PhenoDigm2
uv sync
uv run phenodigm --help
```

### Install as a dependency

The package is not published to PyPI. To use it from another project (for
example an Airflow deployment), reference this repository from a
`requirements.txt`:

```text
# from git (pin a tag/branch/commit):
phenodigm2 @ git+https://github.com/<org>/PhenoDigm2.git@<ref>
# or from a local checkout:
phenodigm2 @ file:///abs/path/to/PhenoDigm2
```

Both entry points are then available: `from phenodigm2 import PhenoDigm` (the
library) and the `phenodigm` / `phenodigm2` console command (the CLI).

Owltools is no longer supported. Current builds use downloaded
Phenio/Semsimian mappings through the `ontology_mapping` action.

## Build overview

A build populates a release-specific directory such as `vTODAY`.

**1. Initialize the build directory.** This materializes `vTODAY/resources/` so
you can review and edit `dependencies.json` before downloading:

- Python: `PhenoDigm("vTODAY").init()`
- CLI: `uv run phenodigm --init vTODAY`

**2. Set the data version.** Check the latest Phenio/Semsimian Zenodo record and
update its ID in `vTODAY/resources/dependencies.json`:

```bash
uvx zenodo_get -w - 18474575
```

**3. Run the pipeline** in either form below.

### Python API (primary)

```python
from phenodigm2 import PhenoDigm

pd = PhenoDigm("vTODAY")           # "vTODAY" is the build dir (CLI: --db)
pd.download()
pd.build()
pd.ontology_mapping()
pd.score(fast=True)                # per-call kwargs map to CLI flags (--fast)
pd.index()
pd.parquet()
```

- Constructor options apply to **every** stage: `PhenoDigm("vTODAY", cores=8, verbose=True)`.
- Per-call keyword arguments **override** them for that call and use the same names as the CLI flags.
- Methods return `self` (so calls chain) and **raise** on error instead of calling `exit()`, so an orchestrator can fail the task.

A runnable Airflow DAG wiring the full pipeline is provided at
[`examples/airflow/phenodigm_dag.py`](examples/airflow/phenodigm_dag.py). In
automation where `dependencies.json` is not hand-edited, step 1 is optional —
`download` creates and seeds the build directory on its first run.

### Command-line interface

The same stages are available as the `phenodigm` command:

```bash
uv run phenodigm download --db vTODAY
uv run phenodigm build --db vTODAY
uv run phenodigm ontology_mapping --db vTODAY
uv run phenodigm score --fast --db vTODAY
uv run phenodigm index --db vTODAY
uv run phenodigm parquet --db vTODAY
```

### Stage reference

| Python method | CLI action | Purpose |
|---|---|---|
| `init()` | `--init` | Create + seed a fresh build directory (fails if it exists) |
| `download()` | `download` | Download declared data inputs |
| `build()` | `build` | Build the SQLite schema and load tables |
| `ontology_mapping()` | `ontology_mapping` | Phenio/Semsimian term-term similarities |
| `score(fast=...)` | `score [--fast]` | Compute disease-model associations |
| `index()` | `index` | Create SQLite indexes |
| `parquet()` | `parquet` | Export the Parquet document bundle |
| `solr_prepare()` | `solr-prepare` | Prepare an optional local Solr bundle |
| `solr()` | `solr` | Build an optional local Solr core |
| `post_process()` | `post-process` | Optional R post-processing; pass `benchmark_only=True` / `--benchmark-only` for Polars-only PheVal outputs |
| `status()` | `status` | Report build status |
| `query()` / `compute()` / `export()` | `query` / `compute` / `export` | Inspect / export a built db |

Method names match the CLI actions with hyphens replaced by underscores;
`pd.run("solr-prepare")` also dispatches by action name (hyphen or underscore).

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

Each has a Python API equivalent (shown first) and a CLI form.

Inspect a completed database — `pd.status()`:

```bash
uv run phenodigm status --db vTODAY
```

Use a non-default ontology-mapping IC threshold —
`pd.ontology_mapping(ontology_mapping_min_ic=<IC>)`:

```bash
uv run phenodigm ontology_mapping \
  --ontology_mapping_min_ic <IC> \
  --db vTODAY
```

Prepare and build an optional local Solr core — `pd.solr_prepare()` then
`pd.solr(solr_url="http://localhost:8984/solr/")`:

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
