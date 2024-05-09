"""Use contents of PhenoDigm2 db to create a Solr core.

@author: Tomasz Konopka
"""

import os.path
import requests
from shutil import rmtree, copytree

from . import tools as pd2tools
from . import solrdata
from . import solrlinks
from . import solrsearch

 
# ############################################################################
#

def initSolrCore(config):
    """Initialize a solr core (solr version 7)"""
    
    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)
    corename = config.solr_corename
     
    if config.db.endswith("/"):
        config.db = config.db[:-1]
    instance = "phenodigm2_" + os.path.basename(config.db)

    # get paths to directories for this core
    coredir, confdir, datadir = pd2tools.getPD2coredir(config)
    
    # unload and remove the core if necessary
    if os.path.exists(coredir):
        pd2tools.log("Removing existing core", 2)
        unload = config.solr_url + "admin/cores?action=UNLOAD"
        unload += "&core="+corename+"&deleteIndex=true"
        requests.get(url=unload)
        rmtree(coredir)
    
    pd2tools.log("Creating new solr core", 2)
    
    # create the core from scratch and let anyone read/write execute
    for onedir in [coredir, datadir]:
        if not os.path.exists(onedir):
            os.makedirs(onedir)
            os.chmod(onedir, 0o777)
    
    # copy the configuration files into the conf directory
    resdirsolr = os.path.join(resdir, "solr7")
    copytree(resdirsolr, confdir)
    os.chmod(confdir, 0o777)
    
    create = config.solr_url + "admin/cores?action=CREATE"
    create += "&name=" + corename + "&instanceDir=mycores/" + instance
    
    r = requests.get(url=create)
    return r.ok, r.text


# ############################################################################
# run this from outside of module

def runSolrCoreBuild(config):
    """Create and fill a Solr core from a Phenodigm2 db."""

    pd2tools.log("Initializing solr core")
    init_ok, init_result = initSolrCore(config)
    if not init_ok:
        pd2tools.log("Error: "+init_result, 2)
        return
    
    pd2tools.log("Transferring data to solr core")
    solrlinks.runSolrGeneGene(config)
    solrdata.runSolrOntologies(config)
    solrdata.runSolrGenes(config)
    solrdata.runSolrDiseases(config)
    solrdata.runSolrMouseModels(config)
    solrlinks.runSolrOntoOnto(config)
    solrsearch.runSolrDiseaseSearch(config)
    solrlinks.runSolrDiseaseGenes(config)
    solrlinks.runSolrDiseaseModels(config)

