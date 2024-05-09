""" Load data from disk files into the sqlite db.

This module just defers the detailed work to other modules.

@author: Tomasz Konopka
"""

from . import tools as pd2tools
from . import loadgenes
from . import loadontologies
from . import loadmodels
from . import loaddiseases


# ############################################################################
# Run this function from outside module to load all data

def runLoad(config):
    """Runs all the load functions"""
    
    pd2tools.log("Loading data into db tables")

    loadgenes.runLoadGenes(config)
    loadontologies.runLoadOntologies(config)
    loaddiseases.runLoadDiseases(config)
    loadmodels.runLoadModels(config)

