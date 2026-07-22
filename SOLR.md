# Building a PhenoDigm2 Solr core

A local Solr core is optional. The standard EBI hand-off is the Parquet bundle,
but a Solr 7.5 core can still be built for compatibility testing.

## Prepare the Solr output bundle

Run the preparation action after the SQLite database is complete:

```bash
uv run phenodigm solr-prepare --db vTODAY
```

This creates an idempotent, self-contained layout:

```text
vTODAY/output/solr/
├── dc-solr-7.5.yml
└── solrcores7.5/
```

No manual copy from the repository root is required. The same command can be
used to add the layout to a release bundle created by an older package version.
An existing Compose file is preserved so release-specific edits are not
overwritten. The Compose file mounts its adjacent `solrcores7.5` directory into
the Solr container.

## Start Solr and build the core

Container lifecycle remains explicit; `phenodigm` does not silently start or
stop Docker. Start the bundled server with:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml up -d
```

The bundled Compose configuration exposes Solr at port 8984. Create and fill
the core with:

```bash
uv run phenodigm solr \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
```

The `solr` action defensively runs the same preparation logic before building.
When `--solr_cores_dir` is omitted, core data is written to
`vTODAY/output/solr/solrcores7.5`, matching the Compose volume. An explicit
`--solr_cores_dir PATH` remains available for an externally managed server.

If Solr is already running outside the bundled Compose setup, `solr-prepare` can
be omitted: the `solr` action will still create any missing bundle files before
building the core. The bundled Compose workflow needs the explicit preparation
step first so its configuration exists before Docker starts.

The default Solr core name is `phenodigm`; override it with
`--solr_corename NAME`. The Solr and Parquet writers use the same document
filters:

- `--output_min_ontology_ontology_score`
- `--output_min_disease_model_2d_score`

Stop the bundled server with:

```bash
docker compose -f vTODAY/output/solr/dc-solr-7.5.yml down
```
