"""Bundle orchestration, publication, cleanup, and rollback tests."""

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

from phenodigm2 import parquet
from phenodigm2.documents import DATASET_SPECS
from phenodigm2.parquet import to_polars_schema


@pytest.fixture
def mocked_bundle_producers(monkeypatch):
    specs = {
        document_type: replace(
            spec,
            producer=lambda config, current=document_type: iter(
                ({"unexpected": current},)
            ),
        )
        for document_type, spec in DATASET_SPECS.items()
    }
    monkeypatch.setattr(parquet, "DATASET_SPECS", specs)


def test_bundle_manifest_default_collision_and_overwrite(
    tmp_path, mocked_bundle_producers
):
    database = tmp_path / "build" / "phenodigm2-build.sqlite"
    database.parent.mkdir()
    database.touch()
    config = SimpleNamespace(
        db=str(database.parent),
        dbfile=database,
        parquet_dir=None,
        overwrite=False,
    )

    target = parquet.runParquetBundleBuild(config)
    assert target == database.parent / "output" / "parquet"
    manifest = json.loads((target / "manifest.json").read_text())
    assert manifest["bundle_version"] == 1
    assert manifest["created_at"].endswith("Z")
    assert manifest["source_database"] == str(database.resolve())
    assert list(manifest["datasets"]) == list(DATASET_SPECS)

    for document_type, spec in DATASET_SPECS.items():
        definition = spec.definition
        dataset = manifest["datasets"][document_type]
        assert dataset["row_count"] == 1
        assert dataset["part_count"] == 1
        assert dataset["paths"] == [f"{document_type}/part-00000.parquet"]
        frame = pl.read_parquet(target / dataset["paths"][0])
        assert frame.schema == pl.Schema(to_polars_schema(definition))
        assert frame["type"].to_list() == [document_type]

    with pytest.raises(FileExistsError):
        parquet.runParquetBundleBuild(config)

    config.overwrite = True
    assert parquet.runParquetBundleBuild(config) == target
    assert not list(target.parent.glob(".parquet.backup-*"))
    assert not list(target.parent.glob(".parquet.staging-*"))


def test_bundle_failure_cleans_staging(tmp_path, monkeypatch):
    target = tmp_path / "bundle"
    config = SimpleNamespace(
        db=str(tmp_path),
        dbfile=tmp_path / "source.sqlite",
        parquet_dir=target,
        overwrite=False,
    )

    def fail(config):
        raise RuntimeError("producer failed")

    specs = dict(parquet.DATASET_SPECS)
    specs["gene"] = replace(specs["gene"], producer=fail)
    monkeypatch.setattr(parquet, "DATASET_SPECS", specs)
    with pytest.raises(RuntimeError, match="producer failed"):
        parquet.runParquetBundleBuild(config)
    assert not target.exists()
    assert not list(tmp_path.glob(".bundle.staging-*"))


def test_overwrite_publish_rolls_back(tmp_path, monkeypatch):
    target = tmp_path / "bundle"
    target.mkdir()
    (target / "old").write_text("old")
    staging = tmp_path / ".bundle.staging-test"
    staging.mkdir()
    (staging / "new").write_text("new")

    real_replace = parquet.os.replace
    failed = False

    def fail_new_publish(source, destination):
        nonlocal failed
        if Path(source) == staging and not failed:
            failed = True
            raise OSError("publish failed")
        return real_replace(source, destination)

    monkeypatch.setattr(parquet.os, "replace", fail_new_publish)
    with pytest.raises(OSError, match="publish failed"):
        parquet._publish_bundle(staging, target, overwrite=True)

    assert (target / "old").read_text() == "old"
    assert not (target / "new").exists()
    assert not staging.exists()
    assert not list(tmp_path.glob(".bundle.backup-*"))
