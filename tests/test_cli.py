"""Tests for top-level PhenoDigm2 command dispatch."""

import pytest

from pd2 import cli


def assert_initialized(build_dir):
    """Assert the directory layout produced by ``--init``."""

    assert build_dir.is_dir()
    assert (build_dir / "resources" / "dependencies.json").is_file()
    assert (build_dir / "data_raw").is_dir()
    assert (build_dir / "data_processed").is_dir()


def test_init_creates_custom_database_directory(tmp_path):
    build_dir = tmp_path / "nested" / "phenodigm-build"

    cli.main(["--init", str(build_dir)])

    assert_initialized(build_dir)


def test_init_requires_bundle_name():
    parser = cli.build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["--init"])


def test_init_fails_when_destination_exists(tmp_path):
    build_dir = tmp_path / "existing-build"
    build_dir.mkdir()

    with pytest.raises(SystemExit, match="destination already exists"):
        cli.main(["--init", str(build_dir)])


def test_init_rejects_db_option(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["--init", "new-build", "--db", str(tmp_path)])


def test_init_is_mutually_exclusive_with_pipeline_actions():
    parser = cli.build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args([])
    with pytest.raises(SystemExit):
        parser.parse_args(["--init", "new-build", "download"])
