"""Display the status of the phenodigm db.

Displays number of rows in all the phenodigm db tables.

@author: Tomasz Konopka
"""

import gzip
import json
import os.path
from os.path import abspath
from os.path import join
from datetime import datetime

from . import tools as pd2tools 
from .dbextractors import *


# ############################################################################
# Some helpers

def statusFilename(config, suffix):
    """Create a filepath for a status file"""
                       
    rootdir = abspath(config.db)
    dbname = os.path.basename(rootdir)           
    suffix += "-" + datetime.now().strftime('%Y-%m-%d--%H-%M')
    return join(rootdir, "status-" + dbname + "-" + suffix + ".txt.gz")
    

def warninfo(msg, warning=True):
    """log a warning or an information string."""
    
    if warning:
        msg = "*** " + msg + " ***"
    pd2tools.log(msg, 4)
    
    
# ############################################################################
# Individual status update functions

def checkTableRows(config):
    """Compute the row count on all tables."""
    
    outfile = statusFilename(config, "rowcounts")
    pd2tools.log("checking: rowcounts", 2)
                         
    # define the directories for outputs as ../downloads/
    rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
    
    # get the tablenames from the db schema
    schemafile = join(resourcesdir, "db_schema.json")
    with open(schemafile, "r") as f:
        dbschema = json.load(f)
    tablenames = sorted(list(dbschema.keys()))
        
    # count number of tables with zero rows
    numzeros = 0
    
    # obtain a connection to interact with the db
    conn = pd2tools.getDBconn(config.dbfile)
    cursor = conn.cursor()    
    
    with gzip.open(outfile, "wt") as f:
        f.write("# Number of rows in db tables\n")
        f.write("table\trows\n")
        # query all tables in the db
        for t in tablenames:
            # find the first field in the table
            field0 = dbschema[t]["fields"][0].split(" ")[0]
            # avoid SELECT COUNT(*); instead SELECT COUNT(col1)        
            sql = "SELECT COUNT("+field0+") FROM " + t
            result = cursor.execute(sql).fetchone()[0]
            f.write(t + "\t" + str(result) + "\n")
            if result == 0:
                numzeros += 1

    conn.close()
    
    # log a message on screen
    warninfo(str(numzeros) + " empty tables", numzeros > 0)
    

def checkDiseases(config):
    """Check if all declared diseases have phenotypes,
    and that all diseases phenotypes have a disease description."""
        
    phenfile = statusFilename(config, "disease-phenotypes")
    descfile = statusFilename(config, "disease-descriptions")
    pd2tools.log("checking: diseases", 2)
    
    # define the directories for outputs as ../downloads/
    rootdir, datadir, resourcesdir, dbdir = pd2tools.getPD2dirs(config)
        
    # obtain a connection to interact with the db
    conn = pd2tools.getDBconn(config.dbfile)
    cur = conn.cursor()

    # get a list of all disease_ids in diseases table
    diseases = dict()
    sql1 = "SELECT id from disease"
    cur.execute(sql1)       
    for row in cur:                
        diseases[row["id"]] = 0                    
    
    # find phenotypes associated with diseases
    sql2 = "SELECT id, phenotype FROM disease_phenotype"
    missing = set()
    cur.execute(sql2)       
    for row in cur:
        if row["id"] in diseases:                
            diseases[row["id"]] = 1
        else:    
            missing.add(row["id"])
    conn.close()
    
    # count the number of diseases without a description
    num_nopheno = 0
    for d in diseases.keys():        
        num_nopheno += (diseases[d] == 0)
            
    # print out messages about missing phenotypes
    with gzip.open(phenfile, "wt") as f:
        f.write("# diseases without phenotypes\n")
        f.write("disease_id\n")
        for did in diseases.keys():
            if diseases[did] == 0:
                f.write(did + "\n")
    warninfo(str(num_nopheno) + " diseases without phenotypes",
             num_nopheno > 0)
    
    # print out messages about missing descriptions        
    with gzip.open(descfile, "wt") as f:
        f.write("# diseases without descriptions\n")
        f.write("disease_id\n")
        for did in missing:
            f.write(did + "\n")
    num_missing = len(missing)
    warninfo(str(num_missing) + " diseases without descriptions",
             num_missing > 0)
    

def checkModels(config):
    """Check if all models have phenotypes,
    and that all model phenotypes have a disease description."""
    
    pd2tools.log("checking: models", 2)
     
    # retrieve model definitions, genotypes, phenotypes
    dbfile = config.dbfile    
    models = getDbMouseModelMap(dbfile)
    genotypes = getDbMapSets(ModelModelGenotype(dbfile))    
    phenomodel = ModelIdPhenotype(dbfile, "model_phenotype")
    phenotypes = getDbMapSets(phenomodel)    
    
    fileheader = "model_id\tcomment\n"
    
    # look for models without phenotypes
    num_nopheno = 0
    phenfile = statusFilename(config, "model-phenotypes")
    with gzip.open(phenfile, "wt") as f:
        f.write("# models without phenotypes\n")
        f.write(fileheader)
        for mid in models.keys():
            if not phenotypes.has(mid):
                f.write(mid+"\t"+"missing in model_phenotype\n")
                num_nopheno += 1
    warninfo(str(num_nopheno) + " models without phenotypes",
             num_nopheno > 0)
    
    # look for models without genotypes
    num_nogeno = 0
    genofile = statusFilename(config, "model-genotypes")
    with gzip.open(genofile, "wt") as f:
        f.write("# models without genotypes\n")
        f.write(fileheader)
        for mid in models.keys():
            if not genotypes.has(mid):
                f.write(mid + "\t" + "missing in model_genotype\n")
                num_nogeno += 1
    warninfo(str(num_nogeno) + " models without genotypes",
             num_nogeno > 0)
    
    # look for models with phenotypes but without descriptions
    num_nodesc = 0
    descfile = statusFilename(config, "model-descriptions")
    with gzip.open(descfile, "wt") as f:
        f.write("# models without descriptions\n")
        f.write(fileheader)
        for mid in phenotypes.keys():
            if mid not in models:
                f.write(mid + "\t" + "missing in model\n")
                num_nodesc += 1
    warninfo(str(num_nodesc) + " models without descriptions",
             num_nodesc > 0)


# ############################################################################
# Run this function from outside module

def statusPhenodigm(config):
    """Run a series of status checks on the db."""
   
    checkTableRows(config)
    checkDiseases(config)
    checkModels(config)
    pd2tools.log("done")

