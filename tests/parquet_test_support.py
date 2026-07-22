"""Shared fixture builders and comparison helpers for document-output tests."""

import sqlite3
from pathlib import Path
from types import SimpleNamespace


def build_document_database(tmp_path: Path) -> SimpleNamespace:
    """Build the minimal SQLite fixture covering all nine document types."""

    database = tmp_path / "source.sqlite"
    with sqlite3.connect(database) as connection:
        connection.executescript(
            """
            CREATE TABLE gene (
                id TEXT, organism TEXT, symbol TEXT, name TEXT,
                altname TEXT, type TEXT, locus TEXT, withdrawn INTEGER
            );
            CREATE TABLE gene_gene_mapping (query TEXT, match TEXT);
            CREATE TABLE ontology (id TEXT, term TEXT);
            CREATE TABLE ontology_synonym (id TEXT, synonym TEXT);
            CREATE TABLE disease (id TEXT, term TEXT, alts TEXT, class TEXT);
            CREATE TABLE disease_gene_mapping (
                query TEXT, match TEXT, locus TEXT, source TEXT
            );
            CREATE TABLE disease_phenotype (id TEXT, phenotype TEXT);
            CREATE TABLE model (
                id TEXT, source TEXT, species TEXT,
                genetic_background TEXT, life_stage TEXT, description TEXT
            );
            CREATE TABLE model_genotype (
                id TEXT, gene_id TEXT, description TEXT
            );
            CREATE TABLE model_phenotype (id TEXT, phenotype TEXT);
            CREATE TABLE ontology_ontology_mapping (
                query TEXT, match TEXT, simJ REAL, ic REAL, lcs TEXT
            );
            CREATE TABLE disease_model_association (
                query TEXT, match TEXT,
                score_avg_norm REAL, score_avg_raw REAL,
                score_max_norm REAL, score_max_raw REAL,
                query_phenotype TEXT, match_phenotype TEXT
            );
            """
        )
        connection.executemany(
            "INSERT INTO gene VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("MGI:1", "mouse", "M1", "mouse", None, "Gene", "1", 0),
                ("MGI:1", "mouse", "OLDM1", None, None, None, None, 1),
                ("HGNC:1", "human", "H1", "human", None, "Gene", "2", 0),
                ("HGNC:1", "human", "OLDH1", None, None, None, None, 1),
                ("MGI:invalid", "mouse", "OLD", None, None, None, None, 1),
            ],
        )
        connection.execute(
            "INSERT INTO gene_gene_mapping VALUES (?, ?)",
            ("HGNC:1", "MGI:1"),
        )
        connection.executemany(
            "INSERT INTO ontology VALUES (?, ?)",
            [
                ("HP:1", "human phenotype"),
                ("HP:2", "other human phenotype"),
                ("MP:1", "mouse phenotype"),
                ("MP:2", "other mouse phenotype"),
            ],
        )
        connection.executemany(
            "INSERT INTO ontology_synonym VALUES (?, ?)",
            [("HP:1", "human phenotype"), ("HP:1", "human synonym")],
        )
        connection.executemany(
            "INSERT INTO disease VALUES (?, ?, ?, ?)",
            [
                ("OMIM:1", "disease one", "ALT:1;; ALT:2", "class-a,class-b"),
                ("ORPHA:2", "disease two", None, None),
            ],
        )
        connection.executemany(
            "INSERT INTO disease_gene_mapping VALUES (?, ?, ?, ?)",
            [
                ("OMIM:1", "HGNC:1", "2", "OMIM"),
                ("ORPHA:2", "MGI:1", "1", "MGI"),
            ],
        )
        connection.executemany(
            "INSERT INTO disease_phenotype VALUES (?, ?)",
            [("OMIM:1", "HP:1"), ("ORPHA:2", "HP:2")],
        )
        connection.executemany(
            "INSERT INTO model VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("IMPC:1", "IMPC", "Mus musculus", "B6", "adult", "impc"),
                ("MGI:1", "MGI", "Mus musculus", "B6", "adult", "mgi one"),
                ("MGI:2", "MGI", "Mus musculus", "B6", "adult", "mgi two"),
                ("OTHER:1", "OTHER", "Rattus norvegicus", "x", "adult", "rat"),
            ],
        )
        connection.executemany(
            "INSERT INTO model_genotype VALUES (?, ?, ?)",
            [
                ("IMPC:1", "MGI:1", "genotype"),
                ("MGI:1", "MGI:1", "genotype"),
                ("MGI:2", "MGI:1", "genotype"),
            ],
        )
        connection.executemany(
            "INSERT INTO model_phenotype VALUES (?, ?)",
            [("IMPC:1", "MP:1"), ("MGI:1", "MP:1"), ("MGI:2", "MP:2")],
        )
        connection.executemany(
            "INSERT INTO ontology_ontology_mapping VALUES (?, ?, ?, ?, ?)",
            [
                ("HP:1", "MP:1", 1.0, 4.0, "root"),
                ("MP:1", "HP:1", 1.0, 4.0, "root"),
                ("HP:2", "MP:2", 0.1, 1.0, "root"),
                ("MP:2", "HP:2", 0.1, 1.0, "root"),
                ("HP:1", "HP:2", 1.0, 4.0, "root"),
            ],
        )
        connection.executemany(
            "INSERT INTO disease_model_association VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("OMIM:1", "IMPC:1", 10.0, 0.1, 20.0, 0.1, "HP:1", "MP:1"),
                ("OMIM:1", "MGI:1", 30.0, 2.0, 40.0, 2.0, "HP:1", "MP:1"),
                ("ORPHA:2", "MGI:2", 50.0, 2.0, 60.0, 2.0, "HP:2", "MP:2"),
                ("ORPHA:2", "IMPC:1", 1.0, 0.1, 1.0, 0.1, "", ""),
            ],
        )

    return SimpleNamespace(
        db=str(tmp_path),
        dbfile=database,
        solr_url="http://unused/",
        solr_corename="unused",
        output_min_ontology_ontology_score=1.5,
        output_min_disease_model_2d_score=2.2,
        parquet_dir=tmp_path / "bundle",
        overwrite=False,
    )


def canonical_rows(rows):
    """Return order-insensitive rows while preserving field/value pairing."""

    def canonical_value(value):
        if isinstance(value, list):
            return tuple(sorted(value))
        return value

    values = [
        tuple((key, canonical_value(value)) for key, value in row.items())
        for row in rows
    ]
    return sorted(values, key=repr)
