"""Unit tests for schema validation, buffering, and Parquet part writing."""

import polars as pl
import pytest

from pd2.documents import get_dataset_spec
from pd2.parquet import ParquetDatasetWriter, to_polars_schema


def test_parquet_writer_multipart_types_and_unknown_fields(tmp_path):
    definition = get_dataset_spec("disease_search").definition
    writer = ParquetDatasetWriter(tmp_path / "disease_search", definition, 2)
    for index in range(5):
        writer.add(
            {
                "disease_id": f"OMIM:{index}",
                "disease_alts": [],
                "search_qf": [f"OMIM:{index}", "term"],
                "human_curated_gene": index % 2 == 0,
                "ignored": "value",
            }
        )
    assert [part.name for part in writer.parts] == [
        "part-00000.parquet",
        "part-00001.parquet",
    ]
    assert writer.finish() is writer
    assert writer.finish() is writer

    assert writer.row_count == 5
    assert [part.name for part in writer.parts] == [
        "part-00000.parquet",
        "part-00001.parquet",
        "part-00002.parquet",
    ]
    frames = [pl.read_parquet(part) for part in writer.parts]
    frame = pl.concat(frames)
    assert frame.schema == pl.Schema(to_polars_schema(definition))
    assert frame["disease_alts"].to_list() == [[], [], [], [], []]
    assert frame["human_curated_gene"].to_list() == [
        True,
        False,
        True,
        False,
        True,
    ]
    assert "ignored" not in frame.columns


def test_parquet_writer_empty_and_schema_errors(tmp_path):
    definition = get_dataset_spec("disease_model_summary").definition
    writer = ParquetDatasetWriter(tmp_path / "empty", definition)
    assert writer.write_all(()) is writer
    with pytest.raises(RuntimeError, match="writer is already finished"):
        writer.add({})
    frame = pl.read_parquet(writer.parts[0])
    assert frame.height == 0
    assert frame.schema == pl.Schema(to_polars_schema(definition))

    bad = ParquetDatasetWriter(tmp_path / "bad", definition)
    with pytest.raises(TypeError, match=r"disease_model_summary\.marker_num_models"):
        bad.add({"marker_num_models": "not an integer"})
    with pytest.raises(TypeError, match=r"disease_model_summary\.association_curated"):
        bad.add({"association_curated": 1})
