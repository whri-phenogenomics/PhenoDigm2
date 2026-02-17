# PhenoDigm2

The PhenoDigm2 software serves two purposes. 

 - It provides an implementation for [phenodigm scoring](https://academic.oup.com/database/article/333089). 
 This scoring technique uses phenotype annotations to assess similarities
 between animal models and diseases. 

 - It creates a database with phenotype annotations and 
 precomputed phenodigm scores. This database can be distributed
 following the licence conditions [explain].



## Installation

### Hardware prerequisites

PhenoDigm2 requires around 15 GB of disk space to store raw and 
processed data (more if computing scores between models and between diseases).

The database creation requires a minimum 28 GB of RAM for the `owltools` 
phase, and about 4 GB per core for the `score` phase (see below).


### Software prerequisites

PhenoDigm2 requires `python3` (v. 3.11) and packages `sqlite3` and `requests`. Installation
of the python framework and dependent packages is covered in several web resources.
For the present application, you should be able to execute `python3` from the command line
and obtain the interactive python prompt.

Within the python prompt, you should be able to execute `import sqlite3` and `import requests`
without errror.

PhenoDigm2 also requires [`owltools`](https://github.com/owlcollab/owltools) (release 0.3.0). 
Its installation procedure is covered in the owltools github repo README. Once set up, you
should be able to execute `owltools -h` from the command line and observe the 
help menu.


### Local setup

To set up PhenoDigm2 on a local computer, fetch the code from the github repo.

```
git clone [REPO_URL]
```

The repo should have a file `phenodigm2.py` in the root directory. 
This is the main program file and running this with `-h` switch should display 
a summary of the command line interface. 

```
python3 phenodigm2.py -h
```


### Setup in Apocrita
PhenoDigm2 runs in Apocrita using `python3` (v 3.6).

Since May 2023, a parsing error appeared after running `python3 code/PhenoDigm2/phenodigm.py build...`. The working code for python3.6 can be found as a tag through `git checkout python3.6`.


### Testing

A limited number of unittests are available under the `pd2tests` folder. 

```
python3 -m unittest pd2tests/*
```


## Usage

### Build

Building a phenodigm database consists of a setup stage and calculation
stage. The procedure is explained in the [BUILD](BUILD.md) page.


### Status

After a successful database build, you can browse the resuling database
with any sqlite client (e.g. in R, in a browser, etc.). The `phenodigm2.py` 
executable also provides a simple utility to view the status of the 
database from the command line. 

```
python3 phenodigm2.py status --db [PATH]
```


### Examples

There are also some utilities to query the database; this can be
useful during manual inspections and to export data to other formats.

These utilities are explained in the [EXAMPLES.md](EXAMPLES.md) document.

