# Procedure for IMPC release

This document describes how to use the `phenodigm2` program to support an IMPC data release. 

Some related documents:

- `BUILD.md` - describes the general procedure for how to build a phenodigm2 database.
- `SOLR.md` - describes how to create a solr core



## Prerequisites

### File locations

These instructions assume:

- that the `PhenoDigm2` repository is available at a location `/code/PhenoDigm2/`. 
- data files will be stored at a location `/data/PhenoDigm2/`.



### Owltools

Owltools is a java program for computing similarities between ontology terms. It is required during the database build stage. A compiled copy is available on apocrita at `/data/WHRI-Phenogenomics/software/opt/owltools`. 

Test that owltools is available and working. 

```
/path/to/owltools/owltools --help
```

(This should display a long list of command-line arguments.)


### Protege

[Protege](https://protege.stanford.edu/) is a java program and GUI for viewing ontology files. It is required for a manual step during the database build stage.


### Python environment

The `phenodigm2` program runs using python (tested with versions 3.6, 3.7, and 3.8). 

Set up a new python environment (here called `venv`), activate it, and install required packages. At the end, deactivate the environment. 

```
cd /data/PhenoDigm2/
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r /code/PhenoDigm2/requirements.txt
# test that phenodigm2 program works
python3 /code/PhenoDigm2/phenodigm2.py --help
deactivate
```

(The help command should display a summary of command-line arguments.)


### Solr server

A solr server is required to store data in the format required for the IMPC website. It is important that the version of the solr server matches that used by the IMPC. The easiest way to achieve this is using docker.

Launch a solr server using `docker-compose` and the provided compose configuration. Test that it responds to queries. When done, stop the solr server.

```
cd /data/PhenoDigm2/
cp /code/PhenoDigm2/dc-solr-7.5.yml .
docker-compose -f dc-solr-7.5.yml up -d
curl http://localhost:8984/solr/
docker-compose -f dc-solr-7.5.yml down
```

(The curl command should display text that looks like html. Alternatively, 
navigate to that address in a browser; the browser should display the solr dashboard.)


## Building a database

To begin a database build, create a new directory and copy the `resources` from the code repository into that new directory.

```
cd /data/PhenoDigm2/
# create new directory - replace TODAY by a version identifier, e.g. a date 
mkdir vTODAY
cp -r /code/PhenoDigm2/resources vTODAY/
```

Activate the python environment. Then download all the required raw data (IMPC, MGI, obofoundry, HGNC, Ensemble, OMIM, ORPHANET).

**NOTE**: To download OMIM data we require an OMIM API key. This can be passed as a local global variable: `OMIM_API_KEY`. If no key/invalid key is passed, morbidmap and mimTitle and  will contain an error. 
```
cd /data/PhenoDigm2/
source venv/bin/activate
export OMIM_API_KEY="<your_key_here>"
python3 /code/PhenoDigm2/phenodigm2.py download --db vTODAY
```

The download can take a few minutes. When complete, manually check that the files in the `raw_data` subdirectory contain the expected data. Note that some ontology-related files will be empty or show error messages. This is unfortunate, but normal. The most useful strategy to check the downloads is to compare file sizes against a previous release.

After the download, manual intervention is required on two files.

- Find file `data_raw/obo/hp.obo` and open it in Protege. Save the contents
  in owl format to replace `data_raw/obo/hp/hp-edit.owl`.
- Find file `data_raw/obo/mp.obo` and open it in Protege. Save the contents in owl format to replace `data_raw/obo/mp/mp-edit.owl`.

After this adjustment, the rest of the pipeline should run without
the need for further intervention. (It is straightforward to put the next few steps into a single script, but keeping them separate creates breakpoints to check files manually, if desired.)

```
# build a database with a subset of tables (~2 mins)
python3 /code/PhenoDigm2/phenodigm2.py build --db vTODAY
# compute HP-MP similarities (~24 hours)
python3 /code/PhenoDigm2/phenodigm2.py owltools --db vTODAY \
                           -owltools /path/to/owltools/owltools
# score associations between diseases and models (~24 hours)
python3 /code/PhenoDigm2/phenodigm2.py score --fast --db vTODAY
# create indexes on database table (~10 min)
python3 /code/PhenoDigm2/phenodigm2.py index --db vTODAY
# optional - print a summary of the database tables (~1 min)
python3 /code/PhenoDigm2/phenodigm2.py status --db vTODAY
``` 

The **owltools** step runs in single-core mode and requires a lot of memory. The default allocation is 26G; lower values will likely fail. Use the command-line argument `--owltools_mem` to change the allocation.

The **score** steps runs in multi-core mode. The default is to use 4 cores. Use argument `--cores` to change this as needed. Each core will require around 4GB of memory.

Note that in the **score** step uses the flag **--fast**. This limits the calculations to associations between diseases and models. If the flag is omitted, the step also calculates associations between pairs of models and pairs of diseases, taking a lot more time and disk space.

Once complete, the database is ready-to-use, copy, transfer, archive. The primary output of the pipeline - the database - is a single file with extension `sqlite`. However, archives on apocrita contain the full set of files.


## Building a solr core

Once a local database is complete, the `phenodigm2` software can transfer the data into a solr core.

An important point when preparing a solr core for an IMPC data release is that the IMPC requires a core named 'phenodigm' irrespective of the release number. When creating a new core, it is necessary to ensure that a core with that name does not already exist. There are two strategies to ensure this. 

One strategy relies on resetting the solr server. To start from scratch, make sure the solr server is not running (i.e. it is down), delete the folder with the core data, and restart the server. 

```
cd /data/PhenoDigm2/
docker-compose -f dc-solr-7.5.yml down
rm -fr solrcores7.5
docker-compose -f dc-solr-7.5.yml up -d
```

An alternative strategy is to find an existing core named 'phenodigm' and change its name. To do this, make sure the solr server is not running (i.e. it is down), navigate into the folder with the core data, edit file 'core.properties' to assign a new name to the existing core, and restart the server.

To proceed, start the solr server and check that a core named 'phenodigm' does not exist. 

```
cd /data/PhenoDigm2/
docker-compose -f dc-solr-7.5.yml up -d
curl 'localhost:8984/solr/phenodigm/select?q=*:*'
```

(This should report that the core does not exist.)

Transfer data from a local database to a solr core.

```
cd /data/PhenoDigm2/
# create a new solr core (~30 mins)
python3 /code/PhenoDigm2/phenodigm2 solr \
  --solr_cores_dir /data/PhenoDigm2/solrcores7.5/ \
  --solr_url http://localhost:8984/solr/ \
  --db vTODAY
```

(The process may emit warnings about invalid identifiers. Messages about MPATH can be ignored; they are ids that exist in the database but do not participate in phenotype scoring.)

Once the data is transfered into solr, check that it is accessible. Some example queries are below.

```
curl 'localhost:8984/solr/phenodigm/select?q=*:*&rows=0&facet=true&facet.field=type'
curl 'localhost:8984/solr/phenodigm/select?q=*:*&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:gene&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:gene_gene&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:ontology&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_gene_summary&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_model_summary&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_search&rows=2'
curl 'localhost:8984/solr/phenodigm/select?q=type:disease_model_summary%20AND%20marker_id:"MGI:1919819"&rows=1'
```

(Each query should return json documents with data.)

As the queries execute, the solr server takes internal steps to optimize how the data is stored. The optimizations takes place in the background and are useful because they reduce the size of the final data bundle (see below). To help the optimizations run, execute a few more queries, restart the server a few times, wait a few minutes, run a few more queries. 

Finally, stop the solr server and prepare a zip bundle. 

```
docker-compose -f dc-solr-7.5.yml down
cd solrcores7.5
zip -r phenodigm_vTODAY.zip phenodigm_vTODAY
```

The resulting zip file is ready to transfer to IMPC.

## Release end notes

In spring 2021, the IMPC team decided to update how disease-model associations are displayed on the data portal. The new proposal requires displaying information about 'matching phenotypes'. These fields are now part of the main solr core.

## Post processing analysis pipeline
In 2023, the [disease models portal](https://diseasemodels.research.its.qmul.ac.uk) was published and in 2024 the pheval benchmarking was completed. To feed these resources, the pipeline was extended and was called `post_process` using [luigi](https://luigi.readthedocs.io/en/stable/running_luigi.html). This pipeline produces the files needed for these resources. **This is an experimental feature and is considered unstable**.

To execute the pipeline you need at least 1 core of 32GB of RAM:

1. Copy the post_process_config.yaml into the db directory
```
cp -r /code/PhenoDigm2/post_process_config.yaml vTODAY/
```
2. Using a text editor of your choice (e.g nano) edit `post_process_config.yaml` to write the path to the following files :
- "omim_curation.tsv"
- "DR_22_Update_DM_pipeline.R
- "hgnc_symbol_checker.R"

Both R scripts are within this repo, so you could pass the following paths and the omim file is available in the HPC PhenoDigm2 directory, so as an example use:
```
post_process_config.yaml

omim_curation_path: "Path/to/omim/curation/file/omim_curation.tsv"
main_r_script_path: "code/PhenoDigm2/RScripts/DR_22_Update_DM_pipeline.R"
hgnc_symbol_checker_script_path: "code/PhenoDigm2/RScripts/auxiliary/hgnc_symbol_checker.R"
```
This will create a copy of the three files into the bundle, which can be useful for tracking.

3. Load a version of R and run the pipeline
```
module load R/4.2.0
python3 /code/PhenoDigm2/phenodigm2.py post_process --db vTODAY
```


