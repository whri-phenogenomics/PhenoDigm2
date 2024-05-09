""" Load ontology data from disk files into the sqlite db.

@author: Tomasz Konopka
"""

from os.path import join as join
from . import tools as pd2tools
from . import dbmodels as pd2models
from . import parsers as pd2parsers


# ############################################################################
# Run this function from outside module 

def runLoadOntologies(config):
    """Load ontology definitions into db."""

    pd2tools.log("Loading: ontologies", 2)
    
    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)    
    obodir = join(downdir, "obo")            
    
    # db interaction classes    
    terms = pd2models.ModelOntology(config.dbfile)
    synonyms = pd2models.ModelOntologySynonym(config.dbfile)      
    
    for ont in ["hp", "mp"]:
        
        # parse the obo file
        obofile = join(obodir, ont + ".obo")
        ont_data = pd2parsers.OboSynonymsParser(obofile)                
        
        # transfer the data from the parser into the db
        for x in ont_data.terms:
            terms.addData(x, ont_data.terms[x])
        terms.save()
                    
        for x in ont_data.synonyms:                    
            synonyms.addData(x, ont_data.synonyms[x])                    
        synonyms.save()

