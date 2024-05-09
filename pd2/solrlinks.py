"""Transfer data from an sqlite db into a solr core.

This module prepares documents that describe links, 
i.e. ontology-ontology mappings or computed disease-model associations.

Produces documents with type: "ontology_ontology", "disease_gene_summary",
"disease_model_summary"

@author: Tomasz Konopka
"""

# in solr modules, use imports with * to avoid pd2models, pd2dss, etc.
from math import sqrt
from .dbextractors import *
from .solrcoremodels import *
from .solrdata import getBaseDiseaseDoc
from . import tools as pd2tools


# ############################################################################
# some helper functions used during generation of solr docs

def getGeneSymbols(dbfile):
    """Get a dict with a mapping from gene id to official symbol."""
    
    genes = getDbGeneMap(dbfile)
    genesymbols = dict()
    for geneid in genes.keys():
        if genes[geneid].isValid():
            genesymbols[geneid] = genes[geneid].official_symbol
    
    return genesymbols
       

def makeDiseaseGeneDocsSymbols(basedoc, geneids, orthologs, genesymbols):
    """Create an array of documents based on basedoc.
    
    The genes (mouse, human) are all using valid (not withdrawn)
    symbols.
    
    geneids - iterable of disease associated genes
    orthologs - map with all possible orthologs
    genesymbols - map of all official symbols (id->symbol)
    genedata - map form id to object that includes official_symbol and others 
    """
    
    # create a set with gene pairs
    genepairs = set()
    for id in geneids:
        # check that the associated gene ids are valid
        if id not in genesymbols:
            continue
        
        # avoid entries that don't have orthologs
        if orthologs.has(id):
            for ortholog in orthologs.get(id):
                if ortholog in genesymbols:
                    if id>ortholog:
                        id, ortholog = ortholog, id
                    genepairs.add(id+" "+ortholog)
        else:
            genepairs.add(id)
                                                                
    def setInDoc(doc, id):
        """Helper to augment a doc with gene id and symbols."""
        if id.startswith("MGI"):
            doc["marker_id"] = id
            doc["marker_symbol"] = genesymbols[id]
        elif id.startswith("HGNC"):
            doc["hgnc_gene_id"] = id
            doc["hgnc_gene_symbol"] = genesymbols[id]
    
    result = []
    for pair in genepairs:
        doc = basedoc.copy()
        ids = pair.split(" ")
        setInDoc(doc, ids[0])
        if len(ids) > 1:
            setInDoc(doc, ids[1])
        result.append(doc)
            
    return result


def makeDiseaseGeneDocs(basedoc, geneids, orthologs, genedata):
    """Create an array of documents based on basedoc.
    
    The genes (mouse, human) are all using valid (not withdrawn)
    symbols.
    
    geneids - iterable of disease associated genes
    orthologs - map with all possible orthologs    
    genedata - map form id to object that includes official_symbol and others 
    """
         
    # create a set with gene pairs
    genepairs = set()
    for id in geneids:
        # check that the associated gene ids are valid
        if id not in genedata:
            continue
        if not genedata[id].isValid():
            continue
        # only consider geneid that are human, i.e. curated by human dbs
        if not id.startswith("HGNC"):
            continue
        
        # construct all possible pairings
        if orthologs.has(id):
            for ortholog in orthologs.get(id):                
                if ortholog in genedata:
                    if id > ortholog:
                        id, ortholog = ortholog, id
                    genepairs.add(id + " " + ortholog)
        else:
            genepairs.add(id)
                                                                               
    def setInDoc(doc, id):
        """Helper to augment a doc with gene id and symbols."""
        hg = "hgnc_gene_"
        mm = "marker_"
        if id.startswith("MGI"):
            doc[mm+"id"] = id
            doc[mm+"symbol"] = genedata[id].official_symbol
            doc[mm+"symbols_withdrawn"] = [_ for _ in genedata[id].symbols]            
        elif id.startswith("HGNC"):
            doc[hg+"id"] = id
            doc[hg+"symbol"] = genedata[id].official_symbol
            doc[hg+"symbols_withdrawn"] = [_ for _ in genedata[id].symbols]
            doc[hg+"locus"] = genedata[id].locus
                
    result = []
    for pair in genepairs:
        doc = basedoc.copy()
        ids = pair.split(" ")
        setInDoc(doc, ids[0])        
        if len(ids) > 1:
            setInDoc(doc, ids[1])
        result.append(doc)
            
    return result


def getMarkerModelCounts(dbfile):
    """Get a map between marker id and (number of models for that gene)"""
     
    # get a mapping from marker gene to a set of models
    mm = ModelModelGenotype(dbfile)
    markermodels = getDbMapSets(mm, field_key=1, field_value=0)
    
    # marker models now provides sets, simplify to just a count
    result = dict()
    for marker in markermodels.keys():
        result[marker] = len(markermodels.getset(marker))
        
    return result
 

def score2d(a, b):
    """Helper function to compute euclidean distance of two scores."""
    return sqrt((a*a)+(b*b))
    

# ############################################################################
# run these from outside the module


def runSolrGeneGene(config):
    """Transfer gene-gene mappings into solr core."""
    
    solr = SolrGeneGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing mappings (type:'"+solr.type+"')", 2)                            
    
    # helper functions to identify HGNC and MGI genes
    def isHGNC(x):
        return x.startswith("HGNC")

    def isMGI(x):
        return x.startswith("MGI")
    
    # scan all gene-gene mapping and transfer to solr 
    ggm = ModelGeneGeneMapping(config.dbfile)
    generator = PhenodigmSimpleGenerator(ggm)
    for row in generator.next():
        if isHGNC(row["query"]) and isMGI(row["match"]):
            doc = dict()
            doc["gene_id"] = row["match"]
            doc["hgnc_gene_id"] = row["query"]
            solr.add(doc)
        
    solr.save()


def runSolrOntoOnto(config):
    """Transfer ontology-ontology mappings into solr core."""
    
    solr = SolrOntoOnto(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing mappings (type:'"+solr.type+"')", 2)

    # retrieve ontology definitions and synonyms
    ontos = getDbOntologyMap(config.dbfile)
    minscore = config.solr_min_mapscore
    
    # retrieve onto-onto mappings from db
    oom = ModelOntologyOntologyMapping(config.dbfile)
    generator = PhenodigmSimpleGenerator(oom)
    for row in generator.next():
        # use symmetry to avoid duplicate documents
        if row["match"] >= row["query"]:
            continue
        # assign ids for cross ontology mappings
        mpid, hpid = row["match"], row["query"]
        if mpid[:2] != "MP" or hpid[:2] != "HP":
            hpid, mpid = mpid, hpid
        if mpid[:2] != "MP" or hpid[:2] != "HP":
            continue  # this indicates a self ontology
        # only consider a link with a minimum score
        if ooscore(row["simJ"], row["ic"]) < minscore:
            continue
        # make sure the ids exist in the onto dictionary
        if mpid not in ontos or hpid not in ontos:
            continue
                                   
        doc = dict()
        doc["mp_id"] = mpid
        doc["hp_id"] = hpid
        doc["mp_term"] = ontos[mpid]
        doc["hp_term"] = ontos[hpid]

        solr.add(doc)

    solr.save()


def runSolrDiseaseGenes(config):
    """create docs for disease-gene associations."""
             
    solr = SolrDiseaseGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing disease-gene (type:'"+solr.type+"')", 2)
    
    dbfile = config.dbfile
    # retrieve just the official symbols
    genedata = getDbGeneMap(dbfile)     
    # retrieve all disease definitions, all phenotypes, all gene data
    ontos = getDbOntologyMap(dbfile)
    orthologs = getDbMapSets(ModelGeneGeneMapping(dbfile))
    diseases = getDbDiseaseMap(dbfile)
    # get mapping from disease -> sets of hgnc gene ids (avoid MGI associations)
    diseasegenemodel = ModelDiseaseGeneMapping(dbfile)
    diseasegenes = MapSets()
    generator = PhenodigmSimpleGenerator(diseasegenemodel)
    for row in generator.next():
        if (row["source"] != "MGI"):
            diseasegenes.add(row["query"], row["match"])
    
    # look at all curated disease/gene relations
    for id in diseases.keys():
        if diseasegenes.has(id):
            # get a basic disease doc
            basedoc = getBaseDiseaseDoc(diseases[id])

            # generate array of docs with disease-gene
            alldocs = makeDiseaseGeneDocs(basedoc, diseasegenes.get(id),
                                          orthologs, genedata)

            # add all the docs
            for doc in alldocs:
                solr.add(doc)
        
        solr.presave()            
        
    solr.save()


def runSolrDiseaseModels(config):
    """create docs for disease-model associations.
    
    The purpose here is to write out hp-mp matching terms."""
    
    solr = SolrDiseaseModel(config.solr_url, config.solr_corename)    
    pd2tools.log("Indexing disease-model (type:'"+solr.type+"')", 2)                            
    
    # thresholds
    min2d = config.solr_min_2dscore
    
    # retrieve background dictionaries about diseases and models
    dbfile = config.dbfile
    genedata = getDbGeneMap(dbfile)
    models = getDbMouseModelMap(dbfile)
    diseases = getDbDiseaseMap(dbfile)
    diseasegenes = getAllDiseaseGenes(dbfile)
    modelgenes = getDbMapSets(ModelModelGenotype(dbfile))
    nmodels = getMarkerModelCounts(dbfile)
    ontos = getDbOntologyMap(dbfile)
                        
    def phen_array(idstring):
        """Turn comma-separated ids into an array,
        e.g. MP:1,MP:2 -> [MP:1 term1, MP:2 term2]"""
        
        if idstring == "":
            return []
        return [_ + " " + ontos[_] for _ in idstring.split(",")]
                                        
    # retrieve disease_model association from db one at a time    
    dmm = ModelAssociation(dbfile)
    dmm.tabname = "disease_model_association"
    generator = PhenodigmSimpleGenerator(dmm)
    for row in generator.next():
        diseaseid, modelid = row["query"], row["match"]        
        mgenes = modelgenes.getset(modelid)
        dgenes = diseasegenes.getset(diseaseid)
        
        mgenes = set([_ for _ in mgenes if _ in genedata])
        if len(mgenes) == 0:
            print("Unknown gene(s) in model: "+modelid)
            continue
        
        # check if this row passes scoring criteria                    
        scorepass = score2d(row["score_avg_raw"], row["score_max_raw"]) > min2d
        # check if this disease-gene is a known association
        curated = not mgenes.isdisjoint(dgenes)

        # to record this association, require manual curation or high scores
        if not scorepass and not curated:
            continue
                                                                                    
        # prepare a doc for this association                                    
        doc = dict()
        doc["disease_id"] = diseaseid
        doc["disease_term"] = diseases[diseaseid]["term"]
        doc["model_id"] = modelid
        for _ in ["source", "description", "genetic_background"]:
            doc["model_" + _] = models[modelid][_]
        doc["marker_id"] = " ".join(mgenes)
        doc["marker_symbol"] = " ".join([genedata[_].official_symbol for _ in mgenes])
        doc["marker_locus"] = " ".join([genedata[_].locus for _ in mgenes])
        doc["marker_num_models"] = sum([nmodels[_] for _ in mgenes])
        doc["association_curated"] = curated
        for _ in ["avg_norm", "avg_raw", "max_norm", "max_raw"]:
            doc["disease_model_" + _] = row["score_" + _]
        doc["disease_matched_phenotypes"] = phen_array(row["query_phenotype"])
        doc["model_matched_phenotypes"] = phen_array(row["match_phenotype"])

        # send the doc to solr                            
        solr.add(doc)
        solr.presave()
        
    solr.save()
    

