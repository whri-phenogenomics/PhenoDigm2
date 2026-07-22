"""End-to-end parity tests across shared producers, Solr, and Parquet."""

import json

import polars as pl

from pd2 import parquet, solrdata, solrlinks, solrsearch
from pd2.document_definitions import materialize_row, normalize_document
from pd2.documents import DATASET_SPECS, get_dataset_spec
from tests.parquet_test_support import canonical_rows


SOLR_WRAPPERS = {
    "gene": solrdata.runSolrGenes,
    "gene_gene": solrlinks.runSolrGeneGene,
    "disease": solrdata.runSolrDiseases,
    "mouse_model": solrdata.runSolrMouseModels,
    "ontology": solrdata.runSolrOntologies,
    "ontology_ontology": solrlinks.runSolrOntoOnto,
    "disease_gene_summary": solrlinks.runSolrDiseaseGenes,
    "disease_model_summary": solrlinks.runSolrDiseaseModels,
    "disease_search": solrsearch.runSolrDiseaseSearch,
}


def test_solr_and_parquet_match_shared_producers(document_database, monkeypatch):
    raw_documents = {
        document_type: list(spec.producer(document_database))
        for document_type, spec in DATASET_SPECS.items()
    }
    posted = []

    def collect_post(url, data, headers):
        posted.extend(json.loads(data))

    monkeypatch.setattr("pd2.solrcoremodels.requests.post", collect_post)
    for document_type, wrapper in SOLR_WRAPPERS.items():
        posted.clear()
        wrapper(document_database)
        definition = get_dataset_spec(document_type).definition
        expected = [
            normalize_document(definition, doc) for doc in raw_documents[document_type]
        ]
        assert canonical_rows(posted) == canonical_rows(expected)

    target = parquet.runParquetBundleBuild(document_database)
    manifest = json.loads((target / "manifest.json").read_text())
    for document_type, spec in DATASET_SPECS.items():
        definition = spec.definition
        paths = manifest["datasets"][document_type]["paths"]
        actual = pl.concat([pl.read_parquet(target / path) for path in paths])
        expected = [
            materialize_row(definition, doc) for doc in raw_documents[document_type]
        ]
        assert canonical_rows(actual.to_dicts()) == canonical_rows(expected)


def test_solr_payload_matches_parquet_bundle(document_database, monkeypatch):
    """Assert direct equality between actual Solr payloads and Parquet rows."""

    target = parquet.runParquetBundleBuild(document_database)
    manifest = json.loads((target / "manifest.json").read_text())
    posted = []

    def collect_post(url, data, headers):
        posted.extend(json.loads(data))

    monkeypatch.setattr("pd2.solrcoremodels.requests.post", collect_post)

    for document_type, wrapper in SOLR_WRAPPERS.items():
        definition = get_dataset_spec(document_type).definition
        posted.clear()
        wrapper(document_database)

        # Solr is sparse and Parquet is dense. Densifying the posted documents
        # models their equivalent Parquet representation because Solr treats an
        # absent field and explicit null identically.
        solr_rows = [
            {name: doc.get(name) for name in definition.schema} for doc in posted
        ]
        paths = manifest["datasets"][document_type]["paths"]
        parquet_rows = pl.concat(
            [pl.read_parquet(target / path) for path in paths]
        ).to_dicts()

        assert canonical_rows(solr_rows) == canonical_rows(parquet_rows)
