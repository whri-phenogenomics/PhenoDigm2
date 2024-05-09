## Building a Phenodigm2 Solr core

### Background

A solr core is an alternative manner of storing and indexing data. 
After building an sqlite database, you can create a solr core by running
the following command


```
python3 phenodigm2.py solr --db [PATH-TO-OUPUT-DIR]
                           --solr_cores_dir [PATH-TO-SOLR-CORES]
```

This does not change the database. It does read the sqlite db and create a new 
core at the indicated location. The core will be saved in a directory matching the 
db name; the core identifier is always `phenodigm2`. 

Option `--solr_url` sets the url of an active Solr server. The default is to use 
a localhost server with port 8983 (the default Solr setting).

Option `--solr_corename` sets the solr name of the core. The default is to use
`phenodigm`.

