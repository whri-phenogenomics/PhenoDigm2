# Building a PhenoDigm2 database

## Background

The repo consists of the following important parts:

 - **phenogigm2.py** is the executable program
 - **pd2** is a folder with the source code
 - **pd2tests** is a folder with unit tests
 - **resources** is a folder with input/configuration files 

 
## Building a database

Preparing the phenodigm database is a multi-step process. The preparation
must be done manually; the other steps can be executed separately or
through a single command.
 

#### Prep [~5 min]

To build a new instance of a phenodigm database, you must first prepare 
an output directory. This might be, for example `version-DATE`, where you
can substitute DATE by the current date.

Create such a directory, then copy the entire `resources` folder from the 
repo into this directory. After this step, you would thus have `resources` 
as a sub-directory of `version-DATE`.

If you would like to customize the build, you can edit appropriate files
in the `resources` directory. For example, the paths for data downloads are 
encoded in a configuration file `dependencies.json`. Further information 
is in that directory's README.md.


### Multi-step build

The build process takes several steps. The concise recipe is as follows 
(including manual intervention steps):

```
python3 phenodigm2.py download --db [PATH-TO-OUPUT-DIR]
```

After the download, manual intervention is required to adjust some of
the ontology files. 

 - Find file `hp-edit.owl` in folder `data_raw/obo/hp/`. If the content
 shows an error, this file should be replaced manually. Open file `data_raw/obo/hp.obo`
 in Protege and save it in owl format to replace `hp-edit.owl`.
 - Find file `mp-edit.owl` in folder `data_raw/obo/mp/`. If the content
 shows an error, this file should be replaced manually, similarly as for hp.
 Open file `data_raw/obo/mp.obo` in Protege and save it in owl format to replace
 `mp-edit.owl`.

After the manual adjustment, the rest of the pipeline should run without
the need for further intervention.

```
python3 phenodigm2.py build --db [PATH-TO-OUPUT-DIR]
python3 phenodigm2.py owltools --db [PATH-TO-OUPUT-DIR]
                               --owltools [PATH-TO-OWLTOOLS]
python3 phenodigm2.py score --fast --db [PATH-TO-OUPUT-DIR]
python3 phenodigm2.py index --db [PATH-TO-OUPUT-DIR]
python3 phenodigm2.py status --db [PATH-TO-OUPUT-DIR]
``` 



#### Download [~10 min]

This stage reads a list of dependencies from the `resources` directory
and downloads data files.

The downloaded files are stored under `data_raw` in the 
output directory, which is itself partitioned into subdirectories.


#### Database build [<1 min]

This stage creates an sqlite database within the output directory and
transfers data into the database. 

The database file is stored in the output directory. It's name starts 
with a prefix `phenodigm2` and includes the name of the directory.

This stage also creates some processed data files. These are stored
under `data_processed` in the output directory.


#### Owltools [~5 hours]

This stage launches `owltools` processes that compute similarities
between phenotype ontology terms. This can be a time-consuming
and memory-intensive step. 

Option `--owltools_mem` determines the amount of heap space allocated
to each java process.

Option `--owltools_min_ic` determines the information content threshold
for outputing term-to-term matches. 


#### Scoring [~30 hours]

This stage uses owltools scores for individual phenotype terms to
compute scores comparing sets of terms. Such sets of terms can be 
associated with diseases or mouse models. This stage thus computes
model-model, model-disease, disease-model, and disease-disease scores.  

Option `--fast` instructs the program to compute only disease-model 
associations. The estimated running time is based on this fast option. 
The full calculation will take several times longer.


#### Indexing [~1 hours]

This stage create database indexes on some database tables.


#### Status

This stage is optional. It does not introduce any changes to the database. 
It does create a number of files in the output directory holding descriptive
statistics about the database tables.


