"""Transfer data from an sqlite db into a solr core.

This module prepares documents for solr searches.

Produces documents with type: "disease_search"

@author: Tomasz Konopka
"""

from .solrcoremodels import SolrDiseaseSearch
from .documents import iter_disease_search_documents
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
    pd2tools.log("Indexing disease-search (type:'" + solr.type + "')", 2)

    for doc in iter_disease_search_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()
