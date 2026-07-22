"""CLI parsing and dispatch tests for Parquet export."""

import pytest

from pd2 import cli


def test_parquet_cli_options_and_dispatch_without_solr_or_prep(tmp_path, monkeypatch):
    output = tmp_path / "custom"
    assert "<db>/output/parquet" in cli.build_parser().format_help()
    parsed = cli.build_parser().parse_args(
        [
            "parquet",
            "--db",
            str(tmp_path),
            "--parquet-dir",
            str(output),
            "--output_min_ontology_ontology_score",
            "1.75",
            "--output_min_disease_model_2d_score",
            "2.5",
            "--overwrite",
        ]
    )
    assert parsed.parquet_dir == str(output)
    assert parsed.output_min_ontology_ontology_score == 1.75
    assert parsed.output_min_disease_model_2d_score == 2.5
    assert not hasattr(parsed, "solr_min_mapscore")
    assert not hasattr(parsed, "solr_min_2dscore")
    assert parsed.overwrite is True

    captured = []
    monkeypatch.setattr(
        cli.tools, "getDBfile", lambda config: tmp_path / "source.sqlite"
    )
    monkeypatch.setattr(
        cli.parquet_export,
        "runParquetBundleBuild",
        lambda config: captured.append(config),
    )
    monkeypatch.setattr(
        cli.prep,
        "runDirPrep",
        lambda config: pytest.fail("parquet must not run directory prep"),
    )

    cli.main(
        [
            "parquet",
            "--db",
            str(tmp_path),
            "--parquet-dir",
            str(output),
            "--overwrite",
        ]
    )

    assert len(captured) == 1
    assert captured[0].dbfile == tmp_path / "source.sqlite"
    assert captured[0].parquet_dir == str(output)
    assert captured[0].overwrite is True
