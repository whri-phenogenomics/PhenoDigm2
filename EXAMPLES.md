# Examples

Usage examples for inspecting a PhenoDigm2 database: targeted queries, on-demand
score computation, and full-table exports.

These are shown **primarily** with the Python API (`PhenoDigm`); the equivalent
`phenodigm` CLI command follows each one. Set up a handle to the build directory
once — replace `[DB]` with the path to the directory holding the SQLite
database:

```python
from phenodigm2 import PhenoDigm

pd = PhenoDigm("[DB]")
```

The CLI equivalents all take the form `phenodigm query --db [DB] ...`. Both forms
print their results to stdout.

## Queries

`query` extracts bits of information from a PhenoDigm database, including a
breakdown of how stored PhenoDigm scores are computed from the raw data.

### What does a phenotype term mean?

A single term, or several space-separated terms (note the quotes around the
list):

```python
pd.query(term="HP:0009695")
pd.query(term="HP:0009695 HP:0025332")
```

```bash
phenodigm query --db [DB] --term "HP:0009695 HP:0025332"
```

### What does a disease id mean?

```python
pd.query(disease="OMIM:114480 OMIM:134300")
```

```bash
phenodigm query --db [DB] --disease "OMIM:114480 OMIM:134300"
```

### What models contain a given modified gene?

```python
pd.query(gene="Irf Rad9a")
```

```bash
phenodigm query --db [DB] --gene "Irf Rad9a"
```

Search strings can be incomplete gene names: `Irf` matches models with mutated
`Irf5` as well as `Irf8`.

### What phenotypes are associated with a model/disease?

Add `phenotype=True` alongside a `model` or `disease` list (it is not possible to
query models and diseases at once):

```python
pd.query(phenotype=True, model="XXX YYY")
pd.query(phenotype=True, disease="OMIM:114480 OMIM:134300")
```

```bash
phenodigm query --db [DB] --phenotype --model "XXX YYY"
```

### How similar are phenotype terms?

Phenio/Semsimian ontology-mapping scores, via `sim=True`:

```python
pd.query(sim=True, term="HP:0009695 HP:0025332")
```

```bash
phenodigm query --db [DB] --sim --term "HP:0009695 HP:0025332"
```

The output is a table for a symmetric matrix. These scores are loaded from the
Phenio/Semsimian mappings and form the raw inputs for PhenoDigm calculations.

### How similar are a model and a disease?

Stored PhenoDigm scores, for any combination of models and diseases, via
`score=True`:

```python
pd.query(score=True, disease="DECIPHER:14 DECIPHER:18", model="MODEL:10")
```

```bash
phenodigm query --db [DB] --score --disease "DECIPHER:14 DECIPHER:18" --model "MODEL:10"
```

The resulting table holds all model-model, model-disease, disease-model, and
disease-disease entries for the given keys.

## Compute

`compute` computes or re-computes PhenoDigm associations between two models, two
diseases, or a model and a disease. It is more verbose than `query` and reports
details even for associations that are not stored in the database due to
thresholding:

```python
pd.compute(model="MODEL:10", disease="OMIM:3000")
```

```bash
phenodigm compute --db [DB] --model "MODEL:10" --disease "OMIM:3000"
```

Compute commands can take a long time to finish. They cache their results, so
repeat calls are much faster.

## Export

`export` extracts entire tables from the database into tab-separated text on
stdout:

```python
pd.export(table="[TABLENAME]")
pd.export(table="[TABLENAME]", where="id LIKE 'MGI%'")   # export part of a table
```

```bash
phenodigm export --db [DB] --table [TABLENAME]
phenodigm export --db [DB] --table [TABLENAME] --where "id LIKE 'MGI%'"
```

To save to disk, redirect or pipe the CLI output through a compressor — a shell
convenience the CLI form is best suited for:

```bash
phenodigm export --db [DB] --table [TABLENAME] | gzip > [TABLENAME].tsv.gz
phenodigm export --db [DB] --table [TABLENAME] \
  --where "id LIKE 'MGI%'" | gzip > [TABLENAME].MGI.tsv.gz
```
