"""Classes for interacting with phenodigm solr core

@author: Tomasz Konopka
"""

import json
import requests

from .document_definitions import (
    STRING,
    Document,
    NormalizedDocument,
    define_document,
    normalize_document,
)
from .documents import get_dataset_spec


# ############################################################################
# Generic Interface for Solr document creation


class PhenodigmSolr:
    """A base class to interact with a solr core."""

    # url for solr server http
    url = "http://localhost:8983/solr/"
    core = "phenodigm2"
    headers = {"content-type": "application/json"}

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
    definition = define_document(type, ("a", STRING), ("b", STRING))

    def __init__(self, solr_url: str, solr_corename: str) -> None:
        """Instance creates an empty set of documents."""
        self.url = solr_url
        self.core = solr_corename
        self.data: list[NormalizedDocument] = []

    def add(self, obj: Document) -> None:
        """Normalize and stage one sparse document for Solr.

        Missing declared fields stay absent because Solr treats them like null.
        This does not send a request yet; use :meth:`save` to post the batch.
        """

        definition = self.definition
        if definition.document_type != self.type or definition.fieldnames != tuple(
            self.fieldnames
        ):
            definition = define_document(
                self.type, *((field, STRING) for field in self.fieldnames)
            )
        doc = normalize_document(definition, obj)
        # append the new document to self.data (not sent to core yet)
        self.data.append(doc)

    def clear(self) -> None:
        """Erase staged data in this object."""
        self.data = []

    def save(self) -> None:
        """Sends a set of documents to the solr core."""

        url = self.url + self.core + "/update?commit=true&wt=json"

        # send documents to solr core in chunks
        for x in range(0, len(self.data), self.insertN):
            x_data = json.dumps(self.data[x : x + self.insertN])
            requests.post(url, data=x_data, headers=self.headers)

        self.clear()

    def presave(self) -> None:
        """Similar to save, but does not save unless data is full."""
        if len(self.data) >= self.insertN:
            self.save()


# ############################################################################
# Implementation of document types


class SolrGene(PhenodigmSolr):
    """Docs describing mouse genes."""

    definition = get_dataset_spec("gene").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrGeneGene(PhenodigmSolr):
    """Docs describing orthology mapping between mouse and human genes."""

    definition = get_dataset_spec("gene_gene").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrDisease(PhenodigmSolr):
    """Docs describing diseases."""

    definition = get_dataset_spec("disease").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrMouseModel(PhenodigmSolr):
    """Docs summarizing mouse models."""

    definition = get_dataset_spec("mouse_model").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrOntology(PhenodigmSolr):
    """Docs holding ontology definitions/synonyms."""

    definition = get_dataset_spec("ontology").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrOntoOnto(PhenodigmSolr):
    """Docs linking ontology terms (based on owlsim)."""

    definition = get_dataset_spec("ontology_ontology").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrDiseaseGene(PhenodigmSolr):
    """Docs with mappings from diseases to genes."""

    definition = get_dataset_spec("disease_gene_summary").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrDiseaseModel(PhenodigmSolr):
    """Docs with computed assocs from diseases to models."""

    definition = get_dataset_spec("disease_model_summary").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)


class SolrDiseaseSearch(PhenodigmSolr):
    """Docs with summary of associations for diseases.
    These docs are suitable for faceted solr search."""

    definition = get_dataset_spec("disease_search").definition
    type = definition.document_type
    fieldnames = list(definition.fieldnames)
