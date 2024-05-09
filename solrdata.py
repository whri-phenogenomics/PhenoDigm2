"""Transfer data from an sqlite db into a solr core.

This module prepares annotation documents, i.e. documents 
prepared from data parsed from annotation files. 

Produces documents with type: "disease", "gene", 
"phenotype", "mouse_model"

@author: Tomasz Konopka
"""

import re

# in solr modules, use imports with * to avoid pd2models, pd2dss, etc.
from .dbextractors import *
from .solrcoremodels import *
from . import tools as pd2tools


# ############################################################################
# helper functions, e.g. using custom merging of sql tables

def getDbModelsForGenes(dbfile):
    """Retrieve associations from gene-ids to models, id->set."""
    
    # a set of genes with impc models
    impc_genes = set()
    # associations from genes to models (all models)
    models_map = MapSets() 
        
    mm = ModelModel(dbfile)
    mmtab = mm.tabname
    mmg = ModelModelGenotype(dbfile)
    mmgtab = mmg.tabname
    
    generator = PhenodigmJoinGenerator(mm, mmg, ["id", "id"]) 
    for row in generator.next():
        # for models with impc source, record gene name
        if re.match("IMPC", row[mmtab+"_source"]):
            impc_genes.add(row[mmgtab+"_gene_id"])
        # for all models, recover mapping from gene id to model id
        models_map.add(row[mmgtab+"_gene_id"], 
                       row[mmtab+"_id"])
  
    return models_map, impc_genes


def getBaseDiseaseDoc(disease_data):
    """Create a solr doc with disease info from a db dict."""
                
    doc = dict()                
    # transfer general disease info into a document        
    for key in ["id", "source", "term", "alts", "classes"]:
        doc["disease_"+key] = disease_data[key]
                    
    return doc


# ############################################################################
# run this from outside the module

def runSolrGenes(config):
    """Transfer gene definitions into solr core."""
    
    solr = SolrGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing genes (type:'"+solr.type+"')", 2)                            
                    
    # retrieve all gene definitions
    genes = getDbGeneMap(config.dbfile)
    for id in genes.keys():
        gdata = genes[id]
        if not gdata.isValid():
            continue
        doc = dict()
        mhprefix = ""            
        if gdata.id.startswith("HGNC"):
            mhprefix = "hgnc_"
        doc[mhprefix+"gene_id"] = gdata.id
        doc[mhprefix+"gene_symbol"] = gdata.official_symbol
        doc[mhprefix+"gene_symbols_withdrawn"] = gdata.symbols
        doc[mhprefix+"gene_locus"] = gdata.locus            

        solr.add(doc)
    solr.save()            


def runSolrDiseases(config):
    """Create documents in core pertaining to diseaess."""
        
    solr = SolrDisease(config.solr_url, config.solr_corename)    
    pd2tools.log("Indexing diseases (type:'"+solr.type+"')", 2)                            

    # retrieve all disease definitions, all ontologies
    dbfile = config.dbfile
    ontos = getDbOntologyMap(dbfile) 
    diseases = getDbDiseaseMap(dbfile)
    phenomodel = ModelIdPhenotype(dbfile, "disease_phenotype")    
    phenotypes = getDbMapSets(phenomodel)
    for id in diseases.keys():                            
        # create a base document for the disease             
        doc = getBaseDiseaseDoc(diseases[id])                                        
        # add-in phenotypes
        phenos = set()
        if phenotypes.has(id):
            for p in phenotypes.get(id):                
                phenos.add(p+" "+ontos[p])
        doc["disease_phenotypes"] = [_ for _ in phenos]
        # add to solr
        solr.add(doc)
    solr.save()


def runSolrOntologies(config):
    """Create documents in core pertaining to phenotypes."""
    
    solr = SolrOntology(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing ontologies (type:'"+solr.type+"')", 2)

    # retrieve ontology definitions and synonyms
    ontos = getDbOntologyMap(config.dbfile)
    synonyms = getDbMapSets(ModelOntologySynonym(config.dbfile))
    for id in ontos.keys():
        term = ontos[id]
        doc = dict()
        doc["phenotype_id"] = id
        doc["phenotype_term"] = term
        doc["ontology"] = id[:2]
        # get synonyms
        doc["phenotype_synonym"] = []
        if synonyms.has(id):
            for _ in synonyms.get(id):
                if _ != term:
                    doc["phenotype_synonym"].append(_)
                        
        solr.add(doc)
    solr.save()
    

def runSolrMouseModels(config):
    """Create documents in core pertaining to mouse models."""
    
    solr = SolrMouseModel(config.solr_url, config.solr_corename) 
    pd2tools.log("Indexing models (type:'"+solr.type+"')", 2)                            
         
    dbfile = config.dbfile 
    # retrieve all model definitions, etc
    ontos = getDbOntologyMap(dbfile)
    models = getDbMouseModelMap(dbfile)
    genotypes = getDbMapSets(ModelModelGenotype(dbfile))
    genes = getDbGeneMap(dbfile)    
    phenomodel = ModelIdPhenotype(dbfile, "model_phenotype")
    phenotypes = getDbMapSets(phenomodel)    
                         
    for id in models.keys():
        mdata = models[id]
        doc = dict()
        modelid = mdata["id"]
        # transfer general info into a document
        for _ in ["id", "source", "description", "genetic_background"]:
            doc["model_"+_] = mdata[_]        
                   
        if genotypes.has(modelid):            
            symbols = set()
            accessions = set()            
            for g in genotypes.get(modelid):
                if g in genes:
                    symbols.add(genes[g].official_symbol)
                else:
                    print("Unknown gene id: " + str(g))
                accessions.add(g)                                        
            doc["marker_id"] = " ".join(accessions)            
            doc["marker_symbol"] = " ".join(symbols)
                        
        # transfer phenotypes
        phenos = set()
        if phenotypes.has(modelid):
            for p in phenotypes.get(modelid):
                if p in ontos:
                    phenos.add(p+" "+ontos[p])
                else:
                    print("Unknown phenotype: " + str(p))
                    phenos.add(p)
        doc["model_phenotypes"] = [_ for _ in phenos]

        if doc["marker_symbol"] != '':       
            solr.add(doc)        

    solr.save()

