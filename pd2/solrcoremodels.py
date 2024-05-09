"""Classes for interacting with phenodigm solr core

@author: Tomasz Konopka
"""

import json
import requests
from collections import OrderedDict


# ############################################################################
# Generic Interface for Solr document creation

class PhenodigmSolr:
    """A base class to interact with a solr core."""
    
    # url for solr server http
    url = "http://localhost:8983/solr/"
    core = "phenodigm2"
    headers = {'content-type': 'application/json'}
    
    # setsep is a character that converts sets into strings
    setsep = " "
    
    # insertN determines the number of docs that are sent at once
    insertN = 100000
    
    # a set of keys 
    # (classes that extend PhenodigmSolr should replace this)
    fieldnames = ["a", "b"]
        
    # the document type
    # (classes that extend PhenodigmSolr should replace this)
    type = "doc"
    
    def __init__(self, solr_url, solr_corename):
        """Instance creates an empty set of documents."""
        self.url = solr_url
        self.core = solr_corename
        self.data = []

    def add(self, obj):
        """Register a document in the class instance.
        This does not send a request to solr yet, use save().
        """
                        
        # create a dict with the right data
        doc = OrderedDict()
        doc["type"] = self.type
        for f in self.fieldnames:
            if f in obj:
                objf = obj[f]
                if type(objf) is set:
                    if len(objf) > 0:
                        objf = self.setsep.join(objf)
                    else:
                        objf = None
                doc[f] = objf
        # append the new document to self.data (not sent to core yet)        
        self.data.append(doc)
        
    def clear(self):
        """Erase staged data in this object."""
        self.data = []

    def save(self):
        """Sends a set of documents to the solr core."""
                        
        url = self.url+self.core+"/update?commit=true&wt=json"
                                             
        # send documents to solr core in chunks
        for x in range(0, len(self.data), self.insertN):
            x_data = json.dumps(self.data[x:x+self.insertN])
            requests.post(url, data=x_data, headers=self.headers)

        self.clear()

    def presave(self):
        """Similar to save, but does not save unless data is full."""
        if len(self.data) >= self.insertN:
            self.save()


# ############################################################################
# Implementation of document types


class SolrGene(PhenodigmSolr):
    """Docs describing mouse genes."""        
    
    type = "gene"
    fieldnames = ["gene_id", "gene_symbol", "gene_symbols_withdrawn", 
                  "gene_locus", "hgnc_gene_id", "hgnc_gene_symbol", 
                  "hgnc_gene_symbols_withdrawn", "hgnc_gene_locus",
                  "impc_model", "mouse_model"]


class SolrGeneGene(PhenodigmSolr):
    """Docs describing orthology mapping between mouse and human genes."""        
    
    type = "gene_gene"
    fieldnames = ["gene_id", "hgnc_gene_id"]
    
    
class SolrDisease(PhenodigmSolr):
    """Docs describing diseases."""    
    
    type = "disease"
    fieldnames = ["disease_id", "disease_source", "disease_term", 
                  "disease_alts", "disease_classes", "disease_phenotypes"]
    
        
class SolrMouseModel(PhenodigmSolr):
    """Docs summarizing mouse models. """    
    
    type = "mouse_model"
    fieldnames = ["model_id", "model_source",
                  "model_description", "model_genetic_background", 
                  "marker_id", "marker_symbol", "model_phenotypes"]
    

class SolrOntology(PhenodigmSolr):
    """Docs holding ontology definitions/synonyms."""
    
    type = "ontology"
    fieldnames = ["ontology", "phenotype_id", 
                  "phenotype_term", "phenotype_synonym"]


class SolrOntoOnto(PhenodigmSolr):
    """Docs linking ontology terms (based on owlsim)."""
    
    type = "ontology_ontology"
    fieldnames = ["mp_id", "mp_term", "hp_id", "hp_term"]
    
    
class SolrDiseaseGene(PhenodigmSolr):
    """Docs with mappings from diseases to genes."""
    
    type = "disease_gene_summary"
    fieldnames = ["disease_id", "disease_term",                   
                  "marker_id", "marker_symbol", 
                  "marker_symbols_withdrawn",
                  "hgnc_gene_id", "hgnc_gene_symbol", 
                  "hgnc_gene_symbols_withdrawn", "hgnc_gene_locus"]        
    
    
class SolrDiseaseModel(PhenodigmSolr):
    """Docs with computed assocs from diseases to models."""
    
    type = "disease_model_summary"
    fieldnames = ["disease_id", "disease_term",
                  "model_id", "model_source", 
                  "model_description", "model_genetic_background",  
                  "marker_id", "marker_symbol", "marker_locus", 
                  "marker_num_models",
                  "disease_model_avg_raw", "disease_model_avg_norm",
                  "disease_model_max_raw", "disease_model_max_norm",                  
                  "association_curated",
                  "disease_matched_phenotypes", "model_matched_phenotypes"]


class SolrDiseaseSearch(PhenodigmSolr):
    """Docs with summary of associations for diseases.
    These docs are suitable for faceted solr search."""
    
    type = "disease_search"
    fieldnames = ["disease_id",  "disease_term", "disease_source",
                  "disease_alts", 
                  "disease_classes", 
                  "search_qf",                   
                  "human_curated_gene",                  
                  "impc_model_with_curated_gene",
                  "impc_model_with_computed_association",
                  "mgi_model_with_curated_gene",
                  "mgi_model_with_computed_association"]
