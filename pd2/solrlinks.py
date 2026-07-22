"""Transfer data from an sqlite db into a solr core.

This module prepares documents that describe links,
i.e. ontology-ontology mappings or computed disease-model associations.

Produces documents with type: "ontology_ontology", "disease_gene_summary",
"disease_model_summary"

@author: Tomasz Konopka
"""

from math import sqrt

from .dbextractors import getDbGeneMap, getDbMapSets
from .dbmodels import ModelModelGenotype
from .solrcoremodels import (
    SolrDiseaseGene,
    SolrDiseaseModel,
    SolrGeneGene,
    SolrOntoOnto,
)
from .documents import (
    iter_disease_gene_documents,
    iter_disease_model_documents,
    iter_gene_gene_documents,
    iter_ontology_ontology_documents,
)
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
                    if id > ortholog:
                        id, ortholog = ortholog, id
                    genepairs.add(id + " " + ortholog)
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
            doc[mm + "id"] = id
            doc[mm + "symbol"] = genedata[id].official_symbol
            doc[mm + "symbols_withdrawn"] = [_ for _ in genedata[id].symbols]
        elif id.startswith("HGNC"):
            doc[hg + "id"] = id
            doc[hg + "symbol"] = genedata[id].official_symbol
            doc[hg + "symbols_withdrawn"] = [_ for _ in genedata[id].symbols]
            doc[hg + "locus"] = genedata[id].locus

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
    return sqrt((a * a) + (b * b))


# ############################################################################
# run these from outside the module


def runSolrGeneGene(config):
    """Transfer gene-gene mappings into solr core."""

    solr = SolrGeneGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing mappings (type:'" + solr.type + "')", 2)

    for doc in iter_gene_gene_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrOntoOnto(config):
    """Transfer ontology-ontology mappings into solr core."""

    solr = SolrOntoOnto(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing mappings (type:'" + solr.type + "')", 2)

    for doc in iter_ontology_ontology_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrDiseaseGenes(config):
    """create docs for disease-gene associations."""

    solr = SolrDiseaseGene(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing disease-gene (type:'" + solr.type + "')", 2)

    for doc in iter_disease_gene_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()


def runSolrDiseaseModels(config):
    """create docs for disease-model associations.

    The purpose here is to write out hp-mp matching terms."""

    solr = SolrDiseaseModel(config.solr_url, config.solr_corename)
    pd2tools.log("Indexing disease-model (type:'" + solr.type + "')", 2)

    for doc in iter_disease_model_documents(config):
        solr.add(doc)
        solr.presave()
    solr.save()
