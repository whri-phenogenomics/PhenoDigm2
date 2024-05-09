"""Main interface to the PhenoDigm2.

@author: Tomasz Konopka
"""

import argparse
from sys import exit
from pd2 import load
from pd2 import dbbuild
from pd2 import export
from pd2 import owl
from pd2 import prep
from pd2 import query
from pd2 import score
from pd2 import solr
from pd2 import status
from pd2 import tools


# ############################################################################
# create a parser object for the phenodigm executable

parser = argparse.ArgumentParser(description="PhenoDigm2")

# for outputs
parser.add_argument("--db", action="store", 
                    help="Database directory")                    

# for internal tuning/thresholding
parser.add_argument("--phenodigm_min_perc", action="store",
                    help="minimal perc value for phenodigm score",
                    default=60)
parser.add_argument("--phenodigm_min_raw", action="store",
                    help="minimal raw value for phenodigm score",
                    default=1.3)
parser.add_argument("--impc_pval", action="store",
                    help="p-value threshold for IMPC statistical tests",
                    default=1e-4)
parser.add_argument("--cores", action="store",
                    help="number of cores used in phenodigm scoring",
                    default=4, type=int)
parser.add_argument("--fast", action="store_true",
                    help="only compute disease-model associations",
                    default=False)
parser.add_argument("--verbose", action="store_true",
                    help="Write extra output to stdout",
                    default=False)

# for "export"
parser.add_argument("--table", action="store", 
                    help="For use with --export",
                    default="")    
parser.add_argument("--where", action="store", 
                    help="for use with --export",
                    default="")


# for "explain"
parser.add_argument("--model", action="store",
                    help="used with --query and --compute, a model id", default="")
parser.add_argument("--disease", action="store",
                    help="used with --query and --compute, a disease id", default="")
parser.add_argument("--gene", action="store",
                    help="used with --explain, a gene id", default="")
parser.add_argument("--term", action="store",
                    help="used with --query, an ontology term", default="")
parser.add_argument("--sim", action="store_true",
                    help="used with --query, show score data", 
                    default=False)
parser.add_argument("--phenotype", action="store_true",
                    help="used with --query, show phenotypes", 
                    default=False)
parser.add_argument("--association", action="store_true",
                    help="used with --query, show phenodigm associations", 
                    default=False)
parser.add_argument("--score", action="store_true",
                    help="used --query, show phenodigm calculations", 
                    default=False)


# running owltools to obtain ontology term-term similarities
parser.add_argument("--owltools", action="store", 
                    help="complete command to run owltools",
                    default="owltools")
parser.add_argument("--owltools_mem", action="store",
                    help="heap size allocation for owltools",
                    default="26G")
parser.add_argument("--owltools_min_ic", action="store",
                    help="minimal ic for owltools output",
                    default=2.5)


# creating solr core, server addresses, new sets of thresholds
parser.add_argument("--solr_url", action="store",
                    help="url for communicating with solr",
                    default="http://localhost:8983/solr/")
parser.add_argument("--solr_cores_dir", action="store",
                    help="directory holding localhost solr instance",
                    default="/tmp/phenodigm2/")
parser.add_argument("--solr_corename", action="store",
                    help="name of core set in core.properties",
                    default="phenodigm")
parser.add_argument("--solr_min_mapscore", action="store",
                    help="minimum ontology-ontology score",
                    default=1.5)
parser.add_argument("--solr_min_2dscore", action="store",
                    help="minimum raw phenodigm score",
                    default=2.2)


# determining what part of the calculation to perform
parser.add_argument("action", action="store",
                    help="Type of calculation/action to perform", 
                    choices=["download", "build", "owltools", "score", 
                             "index", "solr", "query", "export", "compute",
                             "status"])


# ############################################################################
# Execute the program if module is used as an executable

if __name__ == "__main__":
    
    config = parser.parse_args()
    try:        
        config.dbfile = tools.getDBfile(config)
    except:
        exit("Error creating connection to db (check --db)")        
    
    # handle requests for manual inspection
    # and one-off tasks that don't require start/done messages
    if config.action == "query":
        query.queryPhenodigm(config)
        exit()
    if config.action == "compute":    
        query.computeScores(config)
        exit()
    if config.action == "export":    
        export.exportTables(config)
        exit()
        
    tools.log("Starting PhenoDigm2 ["+config.action+"]")

    if config.action == "status":
        status.statusPhenodigm(config)
        exit()
        
    prep.runDirPrep(config)
    if config.action == "download":
        prep.runDownloads(config)     
    if config.action == "build":
        dbbuild.runDBBuild(config)
        load.runLoad(config)
    if config.action == "owltools":
        owl.runOwltools(config)
    if config.action == "score":
        score.runScoring(config)
    if config.action == "index":
        dbbuild.runDBIndexing(config)   
    if config.action == "solr":
        solr.runSolrCoreBuild(config)
        
    tools.log("Done")        

