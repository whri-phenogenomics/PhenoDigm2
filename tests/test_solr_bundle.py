"""Tests for the self-contained Solr output bundle."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from pd2 import cli, solr


def make_config(build_dir, solr_cores_dir=None):
    return SimpleNamespace(db=str(build_dir), solr_cores_dir=solr_cores_dir)


def test_prepare_solr_bundle_uses_output_default_and_preserves_compose(tmp_path):
    build_dir = tmp_path / "vtest"
    build_dir.mkdir()
    config = make_config(build_dir)

    output_dir, compose_path, cores_dir = solr.prepareSolrBundle(config)

    assert output_dir == build_dir / "output" / "solr"
    assert compose_path == output_dir / "dc-solr-7.5.yml"
    assert "./solrcores7.5" in compose_path.read_text()
    assert cores_dir == output_dir / "solrcores7.5"
    assert cores_dir.is_dir()
    assert Path(config.solr_cores_dir) == cores_dir

    compose_path.write_text("release-specific compose\n")
    solr.prepareSolrBundle(config)
    assert compose_path.read_text() == "release-specific compose\n"


def test_prepare_solr_bundle_preserves_custom_core_directory(tmp_path):
    build_dir = tmp_path / "vtest"
    build_dir.mkdir()
    custom_cores = tmp_path / "custom-cores"
    config = make_config(build_dir, custom_cores)

    _, _, cores_dir = solr.prepareSolrBundle(config)

    assert cores_dir == custom_cores.resolve()
    assert config.solr_cores_dir == custom_cores
    assert cores_dir.is_dir()


def test_solr_cli_documents_bundle_default(tmp_path):
    parser = cli.build_parser()
    parsed = parser.parse_args(["solr", "--db", str(tmp_path)])

    assert parsed.solr_cores_dir is None
    assert "<db>/output/solr/solrcores7.5" in parser.format_help()


def test_solr_prepare_action_materializes_bundle(tmp_path):
    build_dir = tmp_path / "existing-release"

    cli.main(["solr-prepare", "--db", str(build_dir)])

    solr_output = build_dir / "output" / "solr"
    assert (solr_output / "dc-solr-7.5.yml").is_file()
    assert (solr_output / "solrcores7.5").is_dir()


def test_solr_build_prepares_bundle_and_reports_core_creation_failure(
    tmp_path, monkeypatch
):
    build_dir = tmp_path / "existing-release"
    build_dir.mkdir()
    config = make_config(build_dir)
    monkeypatch.setattr(solr, "initSolrCore", lambda config: (False, "not ready"))

    with pytest.raises(RuntimeError, match="not ready"):
        solr.runSolrCoreBuild(config)

    assert (build_dir / "output" / "solr" / "dc-solr-7.5.yml").is_file()
