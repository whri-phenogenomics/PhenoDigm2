# Building a PhenoDigm2 database

PhenoDigm2 is a Python package built and run with
[`uv`](https://docs.astral.sh/uv/). Run all commands from the repository root.

## Setup

The package requires Python 3.12 or later. Create or update the project
environment and confirm that the command-line interface is available:

```bash
uv sync
uv run phenodigm --help
```

PhenoDigm2 can be driven two equivalent ways: **primarily** as a Python library
(`from phenodigm2 import PhenoDigm`) and via the `phenodigm` / `phenodigm2`
command-line entry point. This guide walks through the CLI steps in detail; the
equivalent Python API calls are summarized under [Release build](#release-build)
and in the [README stage reference](README.md#stage-reference).

Owltools is no longer supported by the build workflow. Ontology mappings are
produced from the downloaded Phenio/Semsimian mapping files with the
`ontology_mapping` action.

## Release build

The examples use `vTODAY` as the release directory. Replace it with the actual
version or date for the build.

### Python API equivalent

The whole sequence below is also available from Python, which is the primary way
to drive PhenoDigm2 from a script or an Airflow DAG:

```python
from phenodigm2 import PhenoDigm

pd = PhenoDigm("vTODAY")
pd.init()                     # step 1 (edit dependencies.json before download)
pd.download()                 # step 3
pd.build()                    # step 4
pd.ontology_mapping()         # step 5  (ontology_mapping_min_ic=<IC> to tune)
pd.score(fast=True, cores=4)  # step 6
pd.index()                    # step 7
pd.parquet()                  # step 8  (overwrite=True to replace)
```

Per-call keyword arguments mirror the CLI flags used in the detailed steps that
follow. See [`examples/airflow/phenodigm_dag.py`](examples/airflow/phenodigm_dag.py)
for the pipeline wired as an Airflow DAG.

### 1. Initialize the release directory

```bash
uv run phenodigm --init vTODAY
```

Initialization creates the directory layout and copies the packaged resources,
including `dependencies.json`, into `vTODAY/resources`. It fails if the target
already exists.

### 2. Check the Phenio Zenodo record

Before every release download, resolve the Semsimian Zenodo concept record to
its latest version-specific record:

```bash
uvx zenodo_get -w - 18474575
```

`-w -` prints the latest record's direct file URLs to standard output, so it is
slightly simpler than creating a temporary file. Read the numeric ID following
`/records/` in those URLs and replace the existing Semsimian record ID in:

```text
vTODAY/resources/dependencies.json
```

To keep the URL list for release records instead, use:

```bash
uvx zenodo_get -w zenodo_urls 18474575
```

Review the other dependency URLs before continuing. If OMIM downloads are
required, provide a valid API key in the environment:

```bash
export OMIM_API_KEY="<your_key_here>"
```

### 3. Download resources

```bash
uv run phenodigm download --db vTODAY
```

The download action uses `vTODAY/resources/dependencies.json` and writes the
release inputs under `vTODAY/data_raw`. Check file sizes and contents against a
previous release before continuing.

### 4. Build the SQLite database

```bash
uv run phenodigm build --db vTODAY
```

This action creates and populates `phenodigm2-vTODAY.sqlite`. It remains a
required step in the current package: ontology mapping loads its results into
this database, and scoring reads the populated annotation tables.

### 5. Load Phenio ontology mappings

```bash
uv run phenodigm ontology_mapping --db vTODAY
```

The action transforms the downloaded Phenio/Semsimian HP-to-HP and HP-to-MP
files, filters them by information content and the source `phenodigm_score`, and
loads the mappings into SQLite. The `phenodigm_score` value is used only for
filtering and is not written to the generated Parquet caches or SQLite table.

The default minimum ontology-ontology score is `1.5`. To select another value
with the Python API:

```python
pd.ontology_mapping(output_min_ontology_ontology_score=2.0)
```

The equivalent CLI invocation is:

```bash
uv run phenodigm ontology_mapping \
  --output_min_ontology_ontology_score 2.0 \
  --db vTODAY
```

The score threshold can be combined with a non-default information-content
threshold:

```bash
uv run phenodigm ontology_mapping \
  --ontology_mapping_min_ic <IC> \
  --output_min_ontology_ontology_score <SCORE> \
  --db vTODAY
```

These filters are applied when each `phenio-cache-*.parquet` file is created.
An existing cache is reused and logged as `Skipping`, so changing either flag
alone does not rebuild or re-filter that cache. Cache-generation logs record the
thresholds that were applied.

`ontology_mapping` is the positional action. The similarly named
`--ontology_mapping` option configures an executable and is not the action
selector.

### 6. Score disease-model associations

```bash
uv run phenodigm score --fast --db vTODAY
```

`--fast` limits the calculation to disease-model associations. Without it, the
package also calculates model-model and disease-disease associations, requiring
considerably more time and disk space. Use `--cores <N>` to select the scoring
worker count.

### 7. Create database indexes

```bash
uv run phenodigm index --db vTODAY
```

### 8. Write the Parquet release output

```bash
uv run phenodigm parquet --db vTODAY
```

This creates `vTODAY/output/parquet`, containing the document datasets required
by EBI to build the PhenoDigm Solr index. The output is published atomically and
is not replaced unless `--overwrite` is supplied. See [PARQUET.md](PARQUET.md)
for the schemas, bundle layout, filtering options, and custom destination flag.

## Optional actions

Inspect the completed database with:

```bash
uv run phenodigm status --db vTODAY
```

A local Solr core can still be built. First prepare the self-contained Solr
output directory:

```bash
uv run phenodigm solr-prepare --db vTODAY
```

This creates `vTODAY/output/solr/dc-solr-7.5.yml` and the matching
`solrcores7.5` volume directory; no manual copy is needed. The action is safe to
rerun and preserves an existing Compose file. Start the bundled server, then
populate its core:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml up -d
uv run phenodigm solr \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
```

Stop the bundled server when validation is complete:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml down
```

See [SOLR.md](SOLR.md) for server, core-name, URL, and output-directory options.
The release-oriented local Solr procedure is also described in
[IMPC_RELEASE.md](IMPC_RELEASE.md).
