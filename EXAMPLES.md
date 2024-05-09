# Examples

This file contains a number of usage examples. Its parts
include a section of performing targeted queries and performing
exports of entire database tables.


## Examples of `phenodigm2.py query`

`phenodigm2.py query` is a utility that extract bits of information
from a phenodigm database. It can also provide a breakdown of
how phenodigm scores (stored in the db) are computed from the raw data.

The utility covers a number of usage cases. Some are outlined in
a question/answer format here.

All the commands should start with

```
python3 phenodigm2.py query --db [DB] ...
```

Here, replace the `[DB]` placeholder with the path to the directory 
holding the sqlite database. 


### What does a phenotype term mean?

To obtain information about a single term,

```
... --term HP:0009695
```

To obtain information about multiple terms

```
... --term "HP:0009695 HP:0025332"
```

In this case, note the quotations about the space-separated list.


### What does a disease id mean?

To query multiple diseases,

```
... --disease "OMIM:114480 OMIM:134300"
```


### What models contain a given modified gene?

To query multiple diseases,

```
... --gene "Irf Rad9a"
```

Note that the search strings can contain incomplete gene names. 
In this case, for example, the query `Irf` can display 
data for models with mutated `Irf5` as well as `Irf8`. 



### What phenotypes are associated with a model/disease?

The database holds associations between models/diseases and
their phenotypes. To query these associations,

```
... --phenotype --model "XXX YYY"
... --phenotype --disease "OMIM:114480 OMIM:134300"
```

Note the presence of an option `--phenotype` in addition
to the `--model` or `--disease` lists. Also note that
it is not possible to query models and disease at once.  


### How similar are phenotypes terms?

This question is about owltools similarity scores.

```
... --sim --term "HP:0009695 HP:0025332"
```

Note the presence of an option `--score` in addition to
the list of phenotype terms.

The output to this command consists of a table for a symmetric matrix.
The scores are computed by owltools, not phenodigm; they
form the raw inputs for phenodigm calculations (see below). 


### How similar are a model and a disease?

Comparisons of models and diseases are based on phenodigm scores.
The stored values can be queries for any combination of models and
diseases.

```
... --score --disease "DECIPHER:14 DECIPHER:18" --model "MODEL:10"
```

The resulting table holds all model-model, model-disease, disease-model,
and disease-disease entries for all the given keys.   


## Examples of `phenodigm2.py compute`

The `compute` utility computes or re-computes phenodigm associations 
between two models, two diseases, or a model and a disease. The
output is more verbose that the `query` and can provide details
even for those associations that are not recorded in the database due
to thresholding.

```
python3 phenodigm2.py compute --db [DB] 
        --model "MODEL:10" --diseases "OMIM:3000"
```

Compute commands can take a long time to finish. They save results
in a cache file, so repeat queries will be much faster.


## Examples of `phenodigm2.py export`

`phenodigm2.py export` is a utility to extract entire tables
from the database into a tab-separated format. 

The core command has the following form

```
python3 phenodigm2.py export --db [DB] --table [TABLENAME]
```

This displays information on screen. To save to disk, the output 
can be piped into a compressor and then into a target file.

```
python3 phenodigm2.py export --db [DB] --table [TABLENAME]
		| gzip > [TABLENAME].tsv.gz
```


It is also possible to export a part of a table by supplementing
the base query with a `WHERE` clause.

```
python3 phenodigm2.py export --db [DB] --table [TABLENAME]
        --where "id LIKE 'MGI%'" | gzip > [TABLENAME.MGI].tsv.gz
```
