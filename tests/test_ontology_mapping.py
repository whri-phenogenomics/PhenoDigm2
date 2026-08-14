from collections import Counter
from io import BytesIO
from pathlib import Path
import sqlite3
import tarfile
from types import SimpleNamespace

import polars as pl
import pytest

from phenodigm2.dbmodels import ModelOntologyOntologyMapping
from phenodigm2.ontology_mapping import (
    load_one_ontology_cache_file,
    run_ontology_mapping_processing,
    transform_one_ontology_mapping_file,
)


SOURCE_COLUMNS = (
    "subject_id",
    "object_id",
    "jaccard_similarity",
    "ancestor_information_content",
    "phenodigm_score",
    "ancestor_id",
)


def make_config(tmp_path, min_ic=2.5, min_phenodigm_score=1.5):
    root = tmp_path / "build"
    raw_ontology_dir = root / "data_raw" / "obo"
    processed_dir = root / "data_processed"
    raw_ontology_dir.mkdir(parents=True)
    processed_dir.mkdir(parents=True)
    db_file = root / "phenodigm2-build.sqlite"
    with sqlite3.connect(db_file) as connection:
        connection.execute(
            """
            CREATE TABLE ontology_ontology_mapping (
                query TEXT,
                match TEXT,
                simJ REAL,
                ic REAL,
                lcs TEXT
            )
            """
        )
    return SimpleNamespace(
        db=str(root),
        dbfile=str(db_file),
        ontology_mapping_min_ic=min_ic,
        output_min_ontology_ontology_score=min_phenodigm_score,
    )


def write_phenio_archive(path, rows, columns=SOURCE_COLUMNS):
    contents = "\n".join(
        ("\t".join(columns), *("\t".join(row) for row in rows))
    ).encode()
    member = tarfile.TarInfo("release/test_phenio_mapping.tsv")
    member.size = len(contents)
    with tarfile.open(path, "w:gz") as archive:
        archive.addfile(member, BytesIO(contents))


def read_mapping_rows(db_file):
    with sqlite3.connect(db_file) as connection:
        return connection.execute(
            "SELECT query, match, simJ, ic, lcs FROM ontology_ontology_mapping"
        ).fetchall()


def test_transform_writes_filtered_db_shaped_parquet(tmp_path, capsys):
    config = make_config(tmp_path)
    archive_path = (
        Path(config.db) / "data_raw" / "obo" / "HP_vs_HP_semsimian_phenio.tsv.tar.gz"
    )
    write_phenio_archive(
        archive_path,
        [
            ("HP:0001", "HP:0002", "0.75", "3.0", "1.5", " HP:0000 "),
            ("HP:0003", "HP:0004", "0.25", "2.49", "2.0", "HP:0005"),
            ("HP:0006", "HP:0007", "0.5", "3.0", "1.49", "HP:0008"),
        ],
    )

    transform_one_ontology_mapping_file(config, "hp-hp")

    output = capsys.readouterr().out
    assert "ontology_mapping_min_ic=2.5" in output
    assert "output_min_ontology_ontology_score=1.5" in output

    output_path = Path(config.db) / "data_processed" / "phenio-cache-hp-hp.parquet"
    assert output_path.is_file()
    result = pl.read_parquet(output_path)
    assert result.schema == pl.Schema(
        {
            "query": pl.String,
            "match": pl.String,
            "simJ": pl.Float64,
            "ic": pl.Float64,
            "lcs": pl.String,
        }
    )
    assert result.rows() == [("HP:0001", "HP:0002", 0.75, 3.0, "HP_0000;")]


def test_loader_normalizes_mirrors_and_preserves_duplicates(tmp_path, monkeypatch):
    config = make_config(tmp_path)
    cache_path = Path(config.db) / "data_processed" / "phenio-cache-hp-mp.parquet"
    pl.DataFrame(
        {
            "query": [" HP:0001 ", "ZZ:0001", " SELF ", " HP:0001 "],
            "match": [" MP:0001 ", "AA:0001", "SELF", " MP:0001 "],
            "simJ": ["0.5", "0.2", "1.0", "0.5"],
            "ic": ["3.0", "4.0", "5.0", "3.0"],
            "lcs": [" HP_0000; ", "DROP;", " SELF; ", " HP_0000; "],
        }
    ).write_parquet(cache_path)
    monkeypatch.setattr(ModelOntologyOntologyMapping, "insertN", 2)

    load_one_ontology_cache_file(config, "hp-mp")

    assert Counter(read_mapping_rows(config.dbfile)) == Counter(
        {
            ("HP:0001", "MP:0001", 0.5, 3.0, "HP_0000;"): 2,
            ("MP:0001", "HP:0001", 0.5, 3.0, "HP_0000;"): 2,
            ("SELF", "SELF", 1.0, 5.0, "SELF;"): 1,
        }
    )


def test_loader_reports_missing_parquet(tmp_path):
    config = make_config(tmp_path)

    with pytest.raises(
        FileNotFoundError,
        match=r"phenio-cache-hp-hp\.parquet",
    ):
        load_one_ontology_cache_file(config, "hp-hp")


def test_loader_reports_missing_required_columns(tmp_path):
    config = make_config(tmp_path)
    cache_path = Path(config.db) / "data_processed" / "phenio-cache-hp-hp.parquet"
    pl.DataFrame(
        {
            "query": ["HP:0001"],
            "match": ["HP:0002"],
            "simJ": [0.5],
            "ic": [3.0],
        }
    ).write_parquet(cache_path)

    with pytest.raises(ValueError, match=r"Missing required columns.*lcs"):
        load_one_ontology_cache_file(config, "hp-hp")


def test_run_processes_archives_through_parquet_into_sqlite(tmp_path):
    config = make_config(tmp_path)
    ontology_dir = Path(config.db) / "data_raw" / "obo"
    write_phenio_archive(
        ontology_dir / "HP_vs_HP_semsimian_phenio.tsv.tar.gz",
        [("HP:0001", "HP:0002", "0.8", "3.0", "2.0", "HP:0000")],
    )
    write_phenio_archive(
        ontology_dir / "HP_vs_MP_semsimian_phenio.tsv.tar.gz",
        [("HP:0003", "MP:0004", "0.6", "4.0", "2.0", "UPHENO:0001")],
    )

    run_ontology_mapping_processing(config)

    processed_dir = Path(config.db) / "data_processed"
    assert (processed_dir / "phenio-cache-hp-hp.parquet").is_file()
    assert (processed_dir / "phenio-cache-hp-mp.parquet").is_file()
    assert Counter(read_mapping_rows(config.dbfile)) == Counter(
        {
            ("HP:0001", "HP:0002", 0.8, 3.0, "HP_0000;"): 1,
            ("HP:0002", "HP:0001", 0.8, 3.0, "HP_0000;"): 1,
            ("HP:0003", "MP:0004", 0.6, 4.0, "UPHENO_0001;"): 1,
            ("MP:0004", "HP:0003", 0.6, 4.0, "UPHENO_0001;"): 1,
        }
    )


def test_run_is_idempotent_across_reruns(tmp_path):
    config = make_config(tmp_path)
    ontology_dir = Path(config.db) / "data_raw" / "obo"
    write_phenio_archive(
        ontology_dir / "HP_vs_HP_semsimian_phenio.tsv.tar.gz",
        [("HP:0001", "HP:0002", "0.8", "3.0", "2.0", "HP:0000")],
    )
    write_phenio_archive(
        ontology_dir / "HP_vs_MP_semsimian_phenio.tsv.tar.gz",
        [("HP:0003", "MP:0004", "0.6", "4.0", "2.0", "UPHENO:0001")],
    )

    run_ontology_mapping_processing(config)
    first_run = Counter(read_mapping_rows(config.dbfile))

    # Second run reuses the cached Parquet (transform skips) and must clear the
    # table before loading, so row counts stay identical instead of doubling.
    run_ontology_mapping_processing(config)
    second_run = Counter(read_mapping_rows(config.dbfile))

    assert second_run == first_run


def test_transform_drops_rows_with_unparseable_scores(tmp_path):
    config = make_config(tmp_path)
    archive_path = (
        Path(config.db) / "data_raw" / "obo" / "HP_vs_HP_semsimian_phenio.tsv.tar.gz"
    )
    write_phenio_archive(
        archive_path,
        [
            ("HP:0001", "HP:0002", "0.75", "3.0", "2.0", "HP:0000"),
            (
                "HP:0003",
                "HP:0004",
                "0.5",
                "not_a_number",
                "2.0",
                "HP:0005",
            ),
            ("HP:0006", "HP:0007", "bad", "3.0", "2.0", "HP:0008"),
            ("HP:0009", "HP:0010", "0.5", "3.0", "bad", "HP:0011"),
        ],
    )

    # A non-numeric score must drop only that row, not abort the whole file.
    transform_one_ontology_mapping_file(config, "hp-hp")

    output_path = Path(config.db) / "data_processed" / "phenio-cache-hp-hp.parquet"
    result = pl.read_parquet(output_path)
    assert result.rows() == [("HP:0001", "HP:0002", 0.75, 3.0, "HP_0000;")]
