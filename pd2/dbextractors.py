"""Tools for extracting info from sqlite into in-memory objects.

@author: Tomasz Konopka
"""

from .dss import *
from .dbmodels import *


# ############################################################################
# Some multi-purpose retrievers

def getDbOntologyMap(dbfile):
    """Retrieve ontology ids and terms."""
    
    result = dict()
    
    om = ModelOntology(dbfile)
    ontoid = om.fieldnames[0]
    ontoterm = om.fieldnames[1]
    generator = PhenodigmSimpleGenerator(om)    
    for row in generator.next():
        result[row[ontoid]] = row[ontoterm]        
    
    return result


def getDbMapSets(model, field_key=0, field_value=1):
    """Retrieve from db into an id->set format.
    
    Suitable for getting orthologs, genotypes, synonyms
    from two-column db tables.
    
    model - PhenodigmTable model
    field_key - int, position of the models "fieldnames" that 
                will provide the key for the output map
    field_value - int, position of the model's "fieldnames" that
                will give the values in the output map
    """
        
    # get the two fieldnames from the model
    f_key = model.fieldnames[field_key]
    f_value = model.fieldnames[field_value]
    
    # build the mapset
    # result = MapSets()
    result = OptimisedMapSets()
    generator = PhenodigmSimpleGenerator(model)
    for row in generator.next():
        result.add(row[f_key], row[f_value])
        
    return result
    

# ############################################################################
# Some phenodigm specific retrievers

def getDbGeneMap(dbfile):
    """Retrieve gene definitions into a id->GeneInfo map."""
    
    gene_map = dict()
    
    mg = ModelGene(dbfile)
    generator = PhenodigmSimpleGenerator(mg)
    for row in generator.next():
        rid = row["id"]        
        if rid not in gene_map:
            gene_map[rid] = GeneInfo(rid, row["organism"])        
        if row["withdrawn"] == 0:
            gene_map[rid].addSymbol(row["symbol"])
            gene_map[rid].type = row["type"]
            gene_map[rid].name = row["name"]
            gene_map[rid].locus = row["locus"]                        
        else:   
            gene_map[rid].addSymbol(row["symbol"], True)
                    
    return gene_map


def getDbDiseaseMap(dbfile):
    """Retrieve summaries of diseases."""
    
    diseases = dict()
    
    dm = ModelDisease(dbfile)
    generator = PhenodigmSimpleGenerator(dm)
    for row in generator.next():        
        did = row["id"]
        
        ddict = dict()
        for _ in ["id", "term"]:
            ddict[_] = row[_]
        ddict["alts"] = []
        if row["alts"] is not None:
            temp = row["alts"].split(";; ")
            ddict["alts"] = [_ for _ in temp if _ != ""]                    
        ddict["classes"] = ["unclassified"]
        if row["class"] is not None and row["class"] != "":
            ddict["classes"] = row["class"].split(",")        
        ddict["source"] = did.split(":")[0]  
        if ddict["source"] == "ORPHA":
            ddict["source"] = "ORPHANET"
                      
        diseases[did] = ddict        
            
    return diseases


def getDbMouseModelMap(dbfile):
    """Retrieve data about mouse models."""
    
    models = dict()
    
    mm = ModelModel(dbfile)
    generator = PhenodigmSimpleGenerator(mm)
    for row in generator.next():            
        mdict = dict()
        for _ in ["id", "source", "genetic_background", "description"]:
            mdict[_] = row[_]        
        if row["species"] == "Mus musculus":
            models[row["id"]] = mdict
            
    return models

    
def getAllDiseaseGenes(dbfile):
    """Get a mapping from disease id to gene ids.
    This includes human curated genes only and their orthologs
    (i.e. avoid MGI disease-mousegene associations)."""
    
    # build a map from disease to human genes
    # diseasegenes = MapSets()
    diseasegenes = OptimisedMapSets()
    diseasegenemodel = ModelDiseaseGeneMapping(dbfile)    
    generator = PhenodigmSimpleGenerator(diseasegenemodel)
    for row in generator.next():
        if row["source"] != "MGI":
            diseasegenes.add(row["query"], row["match"])
            
    # get maps from genes to genes (including human to mouse)
    orthologs = getDbMapSets(ModelGeneGeneMapping(dbfile))
    # get definitions of genes from ids
    genes = getDbGeneMap(dbfile) 
    
    # create a new MapSets objects with disease->mouse
    # raw = MapSets()
    raw = OptimisedMapSets()    
    for did in diseasegenes.keys():        
        for g1 in diseasegenes.get(did):
            if not g1.startswith("HGNC"):
                continue            
            raw.add(did, g1)                            
            for g2 in orthologs.getset(g1):
                raw.add(did, g2)               

    # check that genes inserted in first pass are valid gene ids
    # (i.e. the ids are not withdrawn)
    # result = MapSets()
    result = OptimisedMapSets()
    for did in raw.keys():
        for g in raw.get(did):
            if genes[g].isValid():
                result.add(did, g)
    return result

