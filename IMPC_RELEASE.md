# Procedure for an IMPC release

This document is the operational checklist for producing PhenoDigm2 data for an
IMPC release.

Related documentation:

- [BUILD.md](BUILD.md) explains each database-build stage.
- [PARQUET.md](PARQUET.md) describes the EBI Parquet bundle and its schemas.
- [SOLR.md](SOLR.md) describes optional local Solr configuration.

## Prerequisites

### Package environment

PhenoDigm2 is a Python package built with `uv` and requires Python 3.12 or
later. From the repository root, synchronize the environment and verify the
installed command:

```bash
uv sync
uv run phenodigm --help
```

The examples below assume the repository root is the current directory and use
`vTODAY` as the release directory. Replace `vTODAY` with the real release name
or an absolute output path as appropriate.

Owltools is no longer supported, and its legacy CLI action has been removed.
Do not run the old Protege/Robot conversion procedure. The current workflow obtains ontology
mappings from Phenio/Semsimian data and loads them with `ontology_mapping`.

### OMIM credentials

The OMIM inputs require a valid API key. Export it before the download step:

```bash
export OMIM_API_KEY="<your_key_here>"
```

Without a valid key, the downloaded OMIM files may contain an error response
instead of release data.

## Build the release

The release build is driven **primarily** from Python via the `PhenoDigm` class;
the numbered steps below are the detailed CLI equivalents, kept for their
per-step validation guidance. The whole sequence in Python:

```python
from phenodigm2 import PhenoDigm

pd = PhenoDigm("vTODAY")
pd.init()                     # step 1 (then update dependencies.json, step 2)
pd.download()                 # step 3 (validate data_raw afterwards)
pd.build()                    # step 4
pd.ontology_mapping()         # step 5  (ontology_mapping_min_ic=<IC> to tune)
pd.score(fast=True, cores=4)  # step 6
pd.index()                    # step 7
pd.parquet()                  # step 8  (overwrite=True to replace an existing bundle)
```

Set `OMIM_API_KEY` in the environment before `pd.download()`, exactly as for the
CLI. See [`examples/airflow/phenodigm_dag.py`](examples/airflow/phenodigm_dag.py)
for the same pipeline wired as an Airflow DAG.

### 1. Initialize

```bash
uv run phenodigm --init vTODAY
```

This creates the release directory and seeds `vTODAY/resources` from the
package. Do not create `vTODAY` manually; initialization deliberately fails if
the destination already exists.

### 2. Update the Phenio Zenodo record

Before downloading, resolve the Semsimian concept record to the latest
version-specific record:

```bash
uvx zenodo_get -w - 18474575
```

The output contains direct URLs for the current files. Copy the numeric ID after
`/records/` and replace the Semsimian record ID in:

```text
vTODAY/resources/dependencies.json
```

The equivalent command that preserves the URL list in a file is:

```bash
uvx zenodo_get -w zenodo_urls 18474575
```

Review all entries in the release copy of `dependencies.json` before the
download. Updating the source file under `phenodigm2/resources` is unnecessary for a
one-off release and would not change an already initialized bundle.

### 3. Download and validate resources

```bash
uv run phenodigm download --db vTODAY
```

Validate the resulting files under `vTODAY/data_raw`. Compare their sizes and
basic contents with the previous IMPC release, paying particular attention to:

- the Phenio/Semsimian archives;
- OMIM files, which should contain data rather than an API error;
- `data_raw/annotations/human_mouse_mapping.txt.gz`; and
- the main IMPC, MGI, HGNC, Ensembl, Orphanet, and ontology inputs.

### 4. Build the database

```bash
uv run phenodigm build --db vTODAY
```

The current implementation still requires this explicit action to create and
populate the SQLite tables before ontology mappings and scores can be loaded.

### 5. Load Phenio ontology mappings

```bash
uv run phenodigm ontology_mapping --db vTODAY
```

If the release requires a different information-content threshold:

```bash
uv run phenodigm ontology_mapping \
  --ontology_mapping_min_ic <IC> \
  --db vTODAY
```

### 6. Score disease-model associations

```bash
uv run phenodigm score --fast --db vTODAY
```

The score stage uses four cores by default. Select another count with
`--cores <N>`. Keep `--fast` for the standard release workflow; omitting it also
calculates disease-disease and model-model scores.

### 7. Create database indexes

```bash
uv run phenodigm index --db vTODAY
```

Optionally inspect the completed database:

```bash
uv run phenodigm status --db vTODAY
```

### 8. Write the EBI Parquet bundle

```bash
uv run phenodigm parquet --db vTODAY
```

The release bundle is written to `vTODAY/output/parquet`. It contains the nine
document types and manifest needed by EBI to build the PhenoDigm Solr index.
See [PARQUET.md](PARQUET.md) for the exact layout and validation details.

If the target already exists after a failed or repeated release run, inspect it
before choosing whether to rerun with `--overwrite`.

Once these steps complete, retain the SQLite database, Parquet bundle, release
resources, and logs together in the release archive.

## Optional: build a local Solr core

Parquet is the standard hand-off for EBI, but a local Solr core can still be
built for compatibility testing. The Solr version must match the target
environment. The existing Docker Compose setup uses Solr 7.5.

Prepare the Compose file and its matching core volume inside the release
bundle:

```bash
uv run phenodigm solr-prepare --db vTODAY
```

This idempotent action creates:

```text
vTODAY/output/solr/
├── dc-solr-7.5.yml
└── solrcores7.5/
```

No manual copy from the repository root is required, including for an older
release bundle. An existing Compose file is preserved. The preparation action
does not start Docker or contact Solr. Start the server explicitly and confirm
that the `phenodigm` core does not already exist:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml up -d
curl 'localhost:8984/solr/phenodigm/select?q=*:*'
```

To reset an existing local test server, stop it, remove its core directory only
after confirming the path, and restart it. Alternatively, rename the existing
core in `core.properties` before starting Solr.

Create and populate the core with the package entry point:

```bash
uv run phenodigm solr \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
```

Without `--solr_cores_dir`, the package writes to the same
`vTODAY/output/solr/solrcores7.5` directory mounted by the bundled Compose file.
The default core name is `phenodigm`. See [SOLR.md](SOLR.md) for all relevant
options. Validate the generated document types, for example:

```bash
curl 'localhost:8984/solr/phenodigm/select?q=*:*&rows=0&facet=true&facet.field=type'
curl 'localhost:8984/solr/phenodigm/select?q=type:gene&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:gene_gene&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:ontology&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_gene_summary&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_model_summary&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_search&rows=2'
```

Stop the server before archiving or transferring its core files:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml down
```

## Optional: post-processing analysis

The `post-process` Luigi workflow produces outputs for the disease models
portal and PheVal benchmarking. It requires at least one core with 32 GB RAM and
a compatible R environment.

The OMIM curation file and bundled R scripts are copied into the release at run
time, so a standard run needs no extra configuration. To override a bundled
input, copy `post_process_config.yaml` into `vTODAY`, edit the relevant path,
then run:

```bash
module load R/4.4.1
export R_LIBS_USER=/data/WHRI-Phenogenomics/projects/PhenoDigm2/post_processing_dependencies/r_lib_paths/R/x86_64-pc-linux-gnu-library/4.4.1
uv run phenodigm post-process --db vTODAY
```
