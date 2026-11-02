# Parquet export

PhenoDigm2 can export the documents used by its Solr core as a portable
Parquet bundle. The export reads only the completed SQLite database and does
not require a running Solr server.

Shown **primarily** with the Python API (`PhenoDigm`), with the equivalent
`phenodigm` CLI command alongside:

```python
from phenodigm2 import PhenoDigm

PhenoDigm("/path/to/build").parquet()
```

```bash
phenodigm parquet --db /path/to/build
```

By default, this creates `/path/to/build/output/parquet`. Use `parquet_dir=` (CLI
`--parquet-dir`) to choose another destination:

```python
PhenoDigm("/path/to/build").parquet(parquet_dir="/path/to/phenodigm-parquet")
```

```bash
phenodigm parquet --db /path/to/build \
  --parquet-dir /path/to/phenodigm-parquet
```

An existing destination is never replaced unless `overwrite=True` (CLI
`--overwrite`) is supplied. The replacement is staged beside the destination and
published only after all datasets and the manifest have been written
successfully.

Solr and Parquet share two document-filtering options:
`output_min_ontology_ontology_score` thresholds `sqrt(simJ * ic)` mappings, and
`output_min_disease_model_2d_score` thresholds the combined average/max raw score
used for computed disease-model associations (CLI: the `--output_min_*` flags).
The ontology threshold is also applied to the source `phenodigm_score` when the
earlier `ontology_mapping` action creates its internal mapping caches. Existing
mapping caches are reused rather than re-filtered when this option changes.

## Shared document pipeline

`DATASET_SPECS` is the central catalog connecting each document type to its
ordered logical schema and SQLite producer. The producers contain all joins,
filtering, association expansion, and derived fields; Solr and Parquet consume
the same producer output.

Normalization adds `type`, discards undeclared fields, and preserves schema
order. Lists remain lists. For historical Solr compatibility, non-empty sets
become space-separated strings and empty sets become null. Solr keeps this
document sparse, omitting absent fields. Parquet materializes the same document
onto every schema column, representing those absent fields as null.

Several compatibility rules are intentionally preserved:

- `ontology_ontology` retains unique, known HP/MP pairs whose composite
  `sqrt(simJ * ic)` score is at least the configured minimum.
- `disease_gene_summary` treats non-MGI disease-gene mappings as curated,
  accepts valid HGNC genes, and expands them to known mouse orthologs. A valid
  HGNC gene with no ortholog mapping remains a human-only document.
- `disease_model_summary` requires at least one known model gene, then includes
  a model when its raw 2D score is strictly greater than the configured minimum
  or its gene set overlaps the disease gene set. The latter is the curated
  association path.
- `disease_search` uses the same computed and curated tests for its facet flags.
  Its historical `impc_*` bucket means every model source other than the
  literal `MGI`, including unknown or other sources; it is not a provenance
  assertion.

## Bundle layout

The bundle contains one directory for each Solr document type. Large datasets
are split into numbered parts, each compressed with Zstandard.

```text
output/parquet/
├── gene/part-00000.parquet
├── gene_gene/part-00000.parquet
├── disease/part-00000.parquet
├── mouse_model/part-00000.parquet
├── ontology/part-00000.parquet
├── ontology_ontology/part-00000.parquet
├── disease_gene_summary/part-00000.parquet
├── disease_model_summary/part-00000.parquet
├── disease_search/part-00000.parquet
└── manifest.json
```

Every dataset has an explicit schema. `type` is the first, non-null string
column; other columns are nullable. Multi-valued fields are native Parquet
lists, while marker fields already represented as space-joined strings in
Solr remain strings. Even an empty dataset contains one zero-row part with the
correct schema.

`manifest.json` is written last. It records bundle schema version 1, UTC
creation time, source database path, every dataset schema, row and part counts,
and relative part paths.

## Adding another output writer

Document schemas use the storage-neutral `FieldKind` values defined in
`phenodigm2/document_definitions.py`. Each output writer provides a complete mapping
from those logical kinds to its own native types; for example, the Parquet
adapter maps `FieldKind.STRING` to `polars.String`. A new writer should add its
own adapter rather than putting writer-specific types in the shared schema.
