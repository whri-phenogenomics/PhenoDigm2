"""Write the nine PhenoDigm document datasets as an atomic Parquet bundle."""

import json
import os
import shutil
import tempfile
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn, Self
from uuid import uuid4

import polars as pl

from . import tools as pd2tools
from .document_definitions import (
    DatasetSpec,
    Document,
    DocumentBuildConfig,
    DocumentDefinition,
    FieldKind,
    NormalizedDocument,
    materialize_row,
)
from .documents import DATASET_SPECS


BUNDLE_VERSION = 1
DEFAULT_BATCH_SIZE = 100_000
DEFAULT_OUTPUT_DIRECTORY = "output/parquet"

POLARS_TYPE_BY_FIELD_KIND = MappingProxyType(
    {
        FieldKind.STRING: pl.String,
        FieldKind.BOOLEAN: pl.Boolean,
        FieldKind.INTEGER: pl.Int64,
        FieldKind.FLOAT: pl.Float64,
        FieldKind.STRING_LIST: pl.List(pl.String),
    }
)


def to_polars_schema(definition: DocumentDefinition) -> dict[str, pl.DataType]:
    """Translate a storage-neutral document definition into a Polars schema."""

    return {
        field.name: POLARS_TYPE_BY_FIELD_KIND[field.kind]
        for field in definition.columns
    }


class ParquetDatasetWriter:
    """Materialize dense rows and stream them into numbered Parquet parts."""

    def __init__(
        self,
        directory: str | Path,
        definition: DocumentDefinition,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> None:
        if batch_size <= 0:
            raise ValueError("batch_size must be greater than zero")
        self.directory = Path(directory)
        self.definition = definition
        self.schema = to_polars_schema(definition)
        self.batch_size = batch_size
        self.data: list[NormalizedDocument] = []
        self.row_count = 0
        self.parts: list[Path] = []
        self._finished = False
        self.directory.mkdir(parents=True, exist_ok=True)

    def add(self, obj: Document) -> None:
        """Normalize, densify, validate, and stage one document row."""

        if self._finished:
            raise RuntimeError(
                f"{self.definition.document_type}: writer is already finished"
            )
        row = materialize_row(self.definition, obj)
        self._validate(row)
        self.data.append(row)
        self.row_count += 1
        if len(self.data) >= self.batch_size:
            self._write_part(self.batch_size)

    def write_all(self, documents: Iterable[Document]) -> Self:
        """Write an iterable of documents and finish the dataset."""

        for document in documents:
            self.add(document)
        return self.finish()

    def finish(self) -> Self:
        """Flush remaining rows, or write one schema-only part if empty."""

        if self._finished:
            return self
        if self.data:
            self._write_part()
        elif not self.parts:
            self._write_part(empty=True)
        self._finished = True
        return self

    def _validate(self, row: NormalizedDocument) -> None:
        for field in self.definition.columns:
            value = row[field.name]
            dtype = self.schema[field.name]
            if value is None:
                if not field.nullable:
                    self._schema_error(field.name, dtype, value)
                continue

            valid = False
            if field.kind is FieldKind.STRING:
                valid = isinstance(value, str)
            elif field.kind is FieldKind.BOOLEAN:
                valid = type(value) is bool
            elif field.kind is FieldKind.INTEGER:
                valid = isinstance(value, Integral) and type(value) is not bool
            elif field.kind is FieldKind.FLOAT:
                valid = isinstance(value, Real) and type(value) is not bool
            elif field.kind is FieldKind.STRING_LIST:
                valid = isinstance(value, list) and all(
                    isinstance(item, str) for item in value
                )
            if not valid:
                self._schema_error(field.name, dtype, value)

    def _schema_error(self, field: str, dtype: pl.DataType, value: object) -> NoReturn:
        raise TypeError(
            f"{self.definition.document_type}.{field}: expected {dtype}, "
            f"got {type(value).__name__}"
        )

    def _write_part(self, row_limit: int | None = None, empty: bool = False) -> None:
        part = self.directory / f"part-{len(self.parts):05d}.parquet"
        rows = self.data if row_limit is None else self.data[:row_limit]
        try:
            if empty:
                frame = pl.DataFrame(schema=self.schema)
            else:
                frame = pl.from_dicts(rows, schema=self.schema, strict=True)
            frame.write_parquet(part, compression="zstd")
        except Exception as error:
            raise RuntimeError(
                f"Failed writing {self.definition.document_type} dataset: {error}"
            ) from error
        self.parts.append(part)
        self.data = self.data[len(rows) :]


def _schema_manifest(definition: DocumentDefinition) -> list[dict[str, object]]:
    return [
        {
            "name": field.name,
            "type": str(POLARS_TYPE_BY_FIELD_KIND[field.kind]),
            "nullable": field.nullable,
        }
        for field in definition.columns
    ]


def _resolve_bundle_target(config: DocumentBuildConfig) -> Path:
    """Resolve the configured bundle location or the default beside the database."""

    configured_target = getattr(config, "parquet_dir", None)
    if configured_target:
        return Path(configured_target).expanduser().resolve()

    database_dir = getattr(config, "db", None)
    if database_dir is None:
        database_dir = Path(config.dbfile).parent
    return (Path(database_dir).expanduser() / DEFAULT_OUTPUT_DIRECTORY).resolve()


def _write_dataset(
    staging: Path, spec: DatasetSpec, config: DocumentBuildConfig
) -> dict[str, object]:
    """Write one registered dataset and return its manifest entry."""

    document_type = spec.document_type
    definition = spec.definition
    pd2tools.log(f"Exporting Parquet dataset (type:'{document_type}')", 2)
    writer = ParquetDatasetWriter(staging / document_type, definition)
    writer.write_all(spec.producer(config))
    return {
        "schema": _schema_manifest(definition),
        "row_count": writer.row_count,
        "part_count": len(writer.parts),
        "paths": [part.relative_to(staging).as_posix() for part in writer.parts],
    }


def _build_manifest(
    config: DocumentBuildConfig, datasets: Mapping[str, dict[str, object]]
) -> dict[str, object]:
    """Build bundle metadata after every dataset has been written."""

    return {
        "bundle_version": BUNDLE_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "source_database": str(Path(config.dbfile).expanduser().resolve()),
        "datasets": datasets,
    }


def _write_manifest(staging: Path, manifest: Mapping[str, object]) -> None:
    """Write the manifest as the final file in the staged bundle."""

    manifest_path = staging / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as output:
        json.dump(manifest, output, indent=2)
        output.write("\n")


def _remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def _publish_bundle(staging: Path, target: Path, overwrite: bool) -> None:
    backup = None
    try:
        if target.exists():
            if not overwrite:
                raise FileExistsError(f"Parquet bundle already exists: {target}")
            backup = target.with_name(f".{target.name}.backup-{uuid4().hex}")
            os.replace(target, backup)

        try:
            os.replace(staging, target)
        except Exception:
            if backup is not None and backup.exists() and not target.exists():
                os.replace(backup, target)
            raise

        if backup is not None:
            _remove_path(backup)
    except Exception:
        if staging.exists():
            _remove_path(staging)
        raise


def runParquetBundleBuild(config: DocumentBuildConfig) -> Path:
    """Build and atomically publish all nine Parquet datasets."""

    target = _resolve_bundle_target(config)
    overwrite = bool(getattr(config, "overwrite", False))

    if target.exists() and not overwrite:
        raise FileExistsError(f"Parquet bundle already exists: {target}")

    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{target.name}.staging-", dir=target.parent)
    )

    try:
        datasets = {
            document_type: _write_dataset(staging, spec, config)
            for document_type, spec in DATASET_SPECS.items()
        }
        _write_manifest(staging, _build_manifest(config, datasets))
        _publish_bundle(staging, target, overwrite)
    except Exception:
        if staging.exists():
            _remove_path(staging)
        raise

    return target


__all__ = [
    "BUNDLE_VERSION",
    "DEFAULT_BATCH_SIZE",
    "DEFAULT_OUTPUT_DIRECTORY",
    "ParquetDatasetWriter",
    "runParquetBundleBuild",
    "to_polars_schema",
]
