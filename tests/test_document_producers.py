"""Schemas, normalization, and domain-rule tests for shared document producers."""

import sqlite3

import polars as pl

from pd2 import documents
from pd2.document_definitions import FieldKind, materialize_row, normalize_document
from pd2.documents import (
    DATASET_SPECS,
    DISEASE_SEARCH_FLAG_DEFAULTS,
    get_dataset_spec,
)
from pd2.parquet import to_polars_schema


def test_document_registry_and_normalization():
    assert list(DATASET_SPECS) == [
        "gene",
        "gene_gene",
        "disease",
        "mouse_model",
        "ontology",
        "ontology_ontology",
        "disease_gene_summary",
        "disease_model_summary",
        "disease_search",
    ]

    definition = get_dataset_spec("gene").definition
    assert definition.schema["gene_id"] is FieldKind.STRING
    assert definition.columns[0].name == "type"
    assert definition.columns[0].nullable is False
    assert to_polars_schema(definition)["gene_id"] == pl.String
    normalized = normalize_document(
        definition,
        {
            "gene_id": "MGI:1",
            "gene_symbols_withdrawn": {"old"},
            "unknown": "discarded",
        },
    )
    assert normalized == {
        "type": "gene",
        "gene_id": "MGI:1",
        "gene_symbols_withdrawn": "old",
    }
    assert "unknown" not in normalized

    empty_set = normalize_document(definition, {"gene_symbols_withdrawn": set()})
    assert empty_set["gene_symbols_withdrawn"] is None
    row = materialize_row(definition, {"gene_id": "MGI:1"})
    assert list(row) == list(definition.schema)
    assert row["gene_symbol"] is None

    search_definition = get_dataset_spec("disease_search").definition
    search_boolean_fields = {
        field.name
        for field in search_definition.fields
        if field.kind is FieldKind.BOOLEAN
    }
    assert set(DISEASE_SEARCH_FLAG_DEFAULTS) == search_boolean_fields
    assert all(value is False for value in DISEASE_SEARCH_FLAG_DEFAULTS.values())


def test_all_document_producers_preserve_domain_rules(document_database):
    raw_documents = {
        document_type: list(spec.producer(document_database))
        for document_type, spec in DATASET_SPECS.items()
    }
    assert {key: len(value) for key, value in raw_documents.items()} == {
        "gene": 2,
        "gene_gene": 1,
        "disease": 2,
        "mouse_model": 3,
        "ontology": 4,
        "ontology_ontology": 1,
        "disease_gene_summary": 1,
        "disease_model_summary": 3,
        "disease_search": 2,
    }
    assert {
        doc.get("gene_id") or doc.get("hgnc_gene_id") for doc in raw_documents["gene"]
    } == {"MGI:1", "HGNC:1"}
    assert raw_documents["ontology_ontology"][0]["mp_id"] == "MP:1"
    assert {
        (doc["disease_id"], doc["model_id"])
        for doc in raw_documents["disease_model_summary"]
    } == {
        ("OMIM:1", "IMPC:1"),
        ("OMIM:1", "MGI:1"),
        ("ORPHA:2", "MGI:2"),
    }
    disease_gene = raw_documents["disease_gene_summary"][0]
    assert disease_gene["marker_symbols_withdrawn"] == ["OLDM1"]
    assert disease_gene["hgnc_gene_symbols_withdrawn"] == ["OLDH1"]

    search = {doc["disease_id"]: doc for doc in raw_documents["disease_search"]}
    assert search["OMIM:1"]["human_curated_gene"] is True
    assert search["OMIM:1"]["impc_model_with_curated_gene"] is True
    assert search["OMIM:1"]["impc_model_with_computed_association"] is False
    assert search["ORPHA:2"]["mgi_model_with_computed_association"] is True
    assert search["ORPHA:2"]["mgi_model_with_curated_gene"] is False


def test_producer_diagnostics_use_project_logger(document_database, monkeypatch):
    messages = []
    monkeypatch.setattr(
        documents,
        "log",
        lambda message, indent=0: messages.append((message, indent)),
    )

    with sqlite3.connect(document_database.dbfile) as connection:
        connection.execute(
            "INSERT INTO model_genotype VALUES (?, ?, ?)",
            ("MGI:1", "MGI:unknown", "unknown"),
        )
        connection.execute(
            "INSERT INTO model_phenotype VALUES (?, ?)",
            ("MGI:1", "MP:unknown"),
        )

    list(documents.iter_mouse_model_documents(document_database))
    assert set(messages) == {
        ("Unknown gene id: MGI:unknown", 4),
        ("Unknown phenotype: MP:unknown", 4),
    }

    messages.clear()
    with sqlite3.connect(document_database.dbfile) as connection:
        connection.execute(
            "INSERT INTO model VALUES (?, ?, ?, ?, ?, ?)",
            (
                "MGI:unknown-model",
                "MGI",
                "Mus musculus",
                "B6",
                "adult",
                "unknown genes",
            ),
        )
        connection.execute(
            "INSERT INTO model_genotype VALUES (?, ?, ?)",
            ("MGI:unknown-model", "MGI:missing", "unknown"),
        )
        connection.execute(
            "INSERT INTO disease_model_association VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ("OMIM:1", "MGI:unknown-model", 1.0, 1.0, 1.0, 1.0, "", ""),
        )

    list(documents.iter_disease_model_documents(document_database))
    assert messages == [("Unknown gene(s) in model: MGI:unknown-model", 4)]


def test_legacy_disease_search_source_buckets():
    assert documents.is_impc("MGI") is False
    assert documents.is_impc("IMPC") is True
    assert documents.is_impc("OTHER") is True
