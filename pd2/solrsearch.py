"""Transfer data from an sqlite db into a solr core.

This module prepares documents for solr searches.

Produces documents with type: "disease_search"

@author: Tomasz Konopka
"""

# in solr modules, use imports with * to avoid pd2models, pd2dss, etc.
from .dbextractors import *
from .solrcoremodels import *
from .solrdata import getBaseDiseaseDoc
from .solrlinks import score2d
from . import tools as pd2tools


# ############################################################################
# some helper functions used during generation of solr docs


def isImpc(source):
    """Determine if a model source qualifies as IMPC.
    This is a rather naive implementation, all that is non MGI qualifies.
    """

    return source != "MGI"
    

# ############################################################################
# run these from outside the module

def runSolrDiseaseSearch(config):
    """create docs for disease-search.
    
    The purpose here is to write out basic disease info and fields
    relevant for search faceting."""
    
    solr = SolrDiseaseSearch(config.solr_url, config.solr_corename)    
    pd2tools.log("Indexing disease-search (type:'"+solr.type+"')", 2)                            
            
    dbfile = config.dbfile 
    # load background information about models (including source)
    models = getDbMouseModelMap(dbfile)   
    # load all disease data and curated disease gene associations
    diseases = getDbDiseaseMap(dbfile)
    diseasegenes = getAllDiseaseGenes(dbfile)
    # load all genes associated with models
    modelgenes = getDbMapSets(ModelModelGenotype(dbfile)) 
    
    # thresholds
    min2d = config.solr_min_2dscore   
            
    # prepare a set of documents in memory
    searchdocs = dict()       
    for id in diseases.keys():                                
        # create a base document for the disease                
        doc = getBaseDiseaseDoc(diseases[id])                                        
        # define disease-search specific fields
        searchqf = [doc["disease_id"], doc["disease_term"]]
        for _ in doc["disease_alts"]:
            searchqf.append(_)
        doc["search_qf"] = searchqf                        
        # add booleans for all fields required in the SolrDiseaseSearch
        for _ in solr.fieldnames:
            if _ not in doc:
                doc[_] = False                        
        # store in memory (will complete just below)
        searchdocs[id] = doc
    
    # encode when diseases have a curated associated gene
    # when a disease is in diseasegenes, it must have an associated gene
    for id in diseasegenes.keys():        
        searchdocs[id]["human_curated_gene"] = True
              
    # scan disease_model associations, and complete disease annotations
    dmm = ModelAssociation(dbfile)    
    dmm.tabname = "disease_model_association"        
    generator = PhenodigmSimpleGenerator(dmm)    
    for row in generator.next():        
                         
        diseaseid, modelid = row["query"], row["match"] 
        id = diseaseid       
        mgenes = modelgenes.getset(modelid)
        dgenes = diseasegenes.getset(diseaseid)
    
        prefix = "mgi"
        if isImpc(models[modelid]["source"]):
            prefix = "impc"                            
                                                        
        # check if this row passes scoring criteria                    
        scorepass = score2d(row["score_avg_raw"], row["score_max_raw"]) > min2d            
        if scorepass:
            searchdocs[id][prefix+"_model_with_computed_association"] = True                                            
                
        # check if this disease-gene is a known association
        curated = not mgenes.isdisjoint(dgenes)#
        if curated:
            searchdocs[id][prefix+"_model_with_curated_gene"] = True         
                                        
    # send all the disease docs to the solr
    for id in searchdocs.keys():
        solr.add(searchdocs[id])    
    solr.save()

