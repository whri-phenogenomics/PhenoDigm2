"""Transfer data from an sqlite db into a solr core.

This module prepares annotation documents, i.e. documents
prepared from data parsed from annotation files.

Produces documents with type: "disease", "gene",
"phenotype", "mouse_model"

@author: Tomasz Konopka
"""

import re

from .dbmodels import ModelModel, ModelModelGenotype, PhenodigmJoinGenerator
from .dss import MapSets
from .solrcoremodels import SolrDisease, SolrGene, SolrMouseModel, SolrOntology
from .documents import (
    iter_disease_documents,
    iter_gene_documents,
    iter_mouse_model_documents,
    iter_ontology_documents,
)
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
        if re.match("IMPC", row[mmtab + "_source"]):
            impc_genes.add(row[mmgtab + "_gene_id"])
        # for all models, recover mapping from gene id to model id
        models_map.add(row[mmgtab + "_gene_id"], row[mmtab + "_id"])

    return models_map, impc_genes


def getBaseDiseaseDoc(disease_data):
    """Create a solr doc with disease info from a db dict."""

    doc = dict()
    # transfer general disease info into a document
    for key in ["id", "source", "term", "alts", "classes"]:
        doc["disease_" + key] = disease_data[key]

    return doc


# ############################################################################
# run this from outside the module


def runSolrGenes(config):
    """Transfer gene definitions into solr core."""

    solr = SolrGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing genes (type:'" + solr.type + "')", 2)

    for doc in iter_gene_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrDiseases(config):
    """Create documents in core pertaining to diseaess."""

    solr = SolrDisease(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing diseases (type:'" + solr.type + "')", 2)

    for doc in iter_disease_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrOntologies(config):
    """Create documents in core pertaining to phenotypes."""

    solr = SolrOntology(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing ontologies (type:'" + solr.type + "')", 2)

    for doc in iter_ontology_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrMouseModels(config):
    """Create documents in core pertaining to mouse models."""

    solr = SolrMouseModel(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing models (type:'" + solr.type + "')", 2)

    for doc in iter_mouse_model_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()
