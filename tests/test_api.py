"""Tests for the programmatic PhenoDigm API (``phenodigm2.api``)."""

import pytest

from phenodigm2 import PhenoDigm


def test_config_defaults_and_db(tmp_path):
    """_config reuses the CLI parser, so it carries the CLI defaults."""
    db = tmp_path / "build"
    config = PhenoDigm(db)._config("build")

    assert config.db == str(db)
    assert config.dbfile.endswith(".sqlite")
    # defaults come straight from cli.build_parser (single source of truth)
    assert config.cores == 4
    assert config.fast is False


def test_config_overrides_and_options(tmp_path):
    """Constructor options apply globally; per-call overrides win."""
    pd = PhenoDigm(tmp_path / "build", verbose=True)
    config = pd._config("score", cores=8)

    assert config.verbose is True  # constructor-level option
    assert config.cores == 8  # per-call override


def test_init_creates_build_layout(tmp_path):
    db = tmp_path / "build"
    result = PhenoDigm(db).init()

    assert (db / "resources").is_dir()
    assert (db / "data_raw").is_dir()
    assert (db / "data_processed").is_dir()
    assert isinstance(result, PhenoDigm)  # methods return self for chaining


def test_init_raises_when_dir_exists(tmp_path):
    """API raises rather than calling exit(), so Airflow can fail the task."""
    db = tmp_path / "build"
    PhenoDigm(db).init()
    with pytest.raises(FileExistsError):
        PhenoDigm(db).init()


def test_run_dispatches_by_action_name(tmp_path):
    db = tmp_path / "build"
    PhenoDigm(db).run("init")  # generic dispatch resolves "init" -> init()
    assert (db / "resources").is_dir()
