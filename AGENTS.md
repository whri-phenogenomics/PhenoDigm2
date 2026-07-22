# Agent Instructions

## Use the Codebase Knowledge Graph First

This repository has a codebase-memory-mcp knowledge graph under
`.codebase-memory/`. Treat it as the first stop for code discovery instead of
rereading or broadly searching the repository.

Preferred workflow:

1. Use `search_graph` to find functions, classes, methods, routes, and important
   symbols.
2. Use `trace_path` to inspect callers, callees, and impact before changing
   behavior.
3. Use `get_code_snippet` with the exact `qualified_name` returned by
   `search_graph` when source for a specific symbol is needed.
4. Use `query_graph` for broader relationship or complexity questions.
5. Use `get_architecture` for a high-level map of packages, dependencies,
   hotspots, and entry points.

The graph project name is `PhenoDigm2`.

The persisted artifact is `.codebase-memory/graph.db.zst`. If the graph is
missing, stale, or unavailable, re-index the repository with `index_repository`
before large-scale code discovery. Do not modify `.codebase-memory/` by hand.

Use normal file reads or text search only when:

- Looking for exact strings, errors, configuration values, or shell snippets.
- Inspecting Markdown, YAML, TOML, shell scripts, lockfiles, or other non-code
  files.
- Verifying the exact current lines before editing.
- Graph results are ambiguous or incomplete.

The graph is a navigation aid, not the source of truth. Always inspect current
file contents before editing and preserve user changes already in the working
tree or index.

## Current Project Model

PhenoDigm2 is a Python package built with `uv`. It builds phenotype annotation
databases, computes PhenoDigm scores, and exports shared Solr/Parquet documents.

Important paths:

- `pyproject.toml`: package metadata, Python requirement, dependencies, and CLI
  entry points.
- `pd2/`: main Python package.
- `pd2/cli.py`: implementation behind the `phenodigm` command.
- `pd2/resources/`: resources bundled into the wheel and seeded into new build
  directories.
- `pd2/resources/dc-solr-7.5.yml`: bundled optional Solr Compose template.
- `pd2/rscripts/`: R scripts bundled for the `post-process` workflow.
- `tests/`: pytest suite.
- `post_process_config.yaml`: optional post-processing overrides.
- `BUILD.md`: canonical database-build workflow.
- `IMPC_RELEASE.md`: operational IMPC release checklist.
- `PARQUET.md` and `SOLR.md`: output-specific documentation.

The supported CLI entry point is `phenodigm`, exposed by the installed package:

```bash
uv run phenodigm --help
```

Do not introduce new documentation or workflows based on
`python3 phenodigm2.py ...`. The `phenodigm2` console alias and root wrapper may
remain for compatibility, but `uv run phenodigm ...` is canonical.

Owltools is no longer supported. Ontology mappings come from downloaded
Phenio/Semsimian files and are processed by the `ontology_mapping` action. Do
not restore the old Owltools action, module, importer resources, or instructions
unless the user explicitly requests that legacy behavior.

## Development Commands

`pyproject.toml` requires Python 3.12 or later. Use `uv` to create and run the
project environment:

```bash
uv sync
uv run phenodigm --help
uv run pytest
uv run ruff check .
uv run ruff format .
```

The repository also provides Tox environments for linting and formatting. If
Tox is not already installed, run it as an isolated tool:

```bash
uvx --from 'tox>=4.24.1' tox -e lint
uvx --from 'tox>=4.24.1' tox -e format
```

Useful non-mutating smoke checks include:

```bash
uv run phenodigm --help
uv run phenodigm status --db PATH
```

When changing behavior, run the narrowest relevant tests first. Run a broader
pytest selection when shared producers, CLI dispatch, resources, database
models, or output adapters are affected. Keep Luigi and dependency deprecation
warnings distinct from actual test failures.

## Canonical Build and Release Workflow

The release examples use `vTODAY` as a placeholder bundle name:

```bash
uv run phenodigm --init vTODAY
uvx zenodo_get -w - 18474575
# Update the Semsimian record ID in vTODAY/resources/dependencies.json.
uv run phenodigm download --db vTODAY
uv run phenodigm build --db vTODAY
uv run phenodigm ontology_mapping --db vTODAY
uv run phenodigm score --fast --db vTODAY
uv run phenodigm index --db vTODAY
uv run phenodigm parquet --db vTODAY
```

Before downloading release resources:

- Resolve Zenodo concept record `18474575` and update the version-specific
  Semsimian record ID in the initialized bundle's `resources/dependencies.json`.
- Review the remaining dependency URLs.
- Ensure `OMIM_API_KEY` is set when OMIM files are required.
- Validate downloads against a previous release before continuing.

The explicit `build` action remains required before ontology mapping and
scoring. Use `--ontology_mapping_min_ic <IC>` only when the release needs a
non-default ontology-mapping threshold.

Download, build, mapping, scoring, indexing, post-processing, Docker, and full
release operations may require network access, substantial memory, many hours,
or external services. Do not run these heavy or externally mutating stages
unless the user explicitly requests them. Prefer unit tests and temporary
fixtures for verification.

## Output Workflows

### Parquet

```bash
uv run phenodigm parquet --db vTODAY
```

The default destination is `vTODAY/output/parquet`. The bundle is published
atomically and an existing target requires `--overwrite`. Solr and Parquet use
the same document definitions and producers; changes to shared producers must
be checked for output parity. Read `PARQUET.md` before changing schemas, bundle
layout, filters, manifests, or publication behavior.

### Optional local Solr

Prepare the self-contained Solr output, start the container explicitly, build
the core, then stop the container:

```bash
uv run phenodigm solr-prepare --db vTODAY
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml up -d
uv run phenodigm solr \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml down
```

`solr-prepare` is idempotent, preserves an existing Compose file, and creates
`vTODAY/output/solr/solrcores7.5`. The `solr` action also checks preparation
defensively. Do not silently start, stop, delete, or reset Docker containers or
Solr core directories. Container lifecycle and destructive cleanup require
explicit user intent. Read `SOLR.md` before changing Solr behavior.

## Working Guidelines

- Prefer focused changes that match the existing Python style.
- Treat `pyproject.toml` and the current documentation as authoritative over
  legacy comments or commands.
- Preserve unrelated staged, unstaged, and untracked user work.
- Keep generated databases, release bundles, Parquet outputs, Solr cores,
  caches, wheels, and temporary build artifacts out of commits unless the user
  explicitly asks otherwise.
- Add package runtime assets under `pd2/resources/`, not only at repository
  root, and verify important assets are included in a built wheel.
- Use temporary directories for output workflow tests. Do not run a real Solr
  server or download release datasets merely to test path preparation.
- For database builds, downloads, IMPC releases, Parquet, Solr, or
  post-processing, read the relevant guide before running commands.
- Update documentation and tests when changing CLI choices, defaults, bundle
  paths, package resources, or release steps.
