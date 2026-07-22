"""Transform and load Phenio ontology-mapping files.

@author: Diego Pava
Based on code by Tomasz Konopka

"""

from pathlib import Path
from typing import NamedTuple
import tarfile

import polars as pl

from . import tools as pd2tools
from . import dbmodels as pd2models


# Semantic handles for the columns the transform/load logic names directly.
QUERY_COLUMN, MATCH_COLUMN = "query", "match"
SIMILARITY_COLUMN, INFORMATION_CONTENT_COLUMN, LCS_COLUMN = "simJ", "ic", "lcs"


class _Column(NamedTuple):
    phenio: str  # header name in the source Phenio TSV
    cache: str  # name in our Parquet cache and the DB table
    dtype: pl.DataType


# Single source of truth — row order matches the DB table column order.
_COLUMNS = (
    _Column("subject_id", QUERY_COLUMN, pl.String),
    _Column("object_id", MATCH_COLUMN, pl.String),
    _Column("jaccard_similarity", SIMILARITY_COLUMN, pl.Float64),
    _Column("ancestor_information_content", INFORMATION_CONTENT_COLUMN, pl.Float64),
    _Column("ancestor_id", LCS_COLUMN, pl.String),
)

PHENIO_COLUMNS = tuple(c.phenio for c in _COLUMNS)
CACHE_COLUMNS = tuple(c.cache for c in _COLUMNS)
PHENIO_TO_CACHE_COLUMNS = {c.phenio: c.cache for c in _COLUMNS}
CACHE_TO_PHENIO_COLUMNS = {c.cache: c.phenio for c in _COLUMNS}
CACHE_SCHEMA = {c.cache: c.dtype for c in _COLUMNS}


def _validate_required_columns(source, required_columns, source_name):
    """Raise a clear error when a lazy data source lacks required columns."""

    available_columns = set(source.collect_schema().names())
    missing_columns = sorted(set(required_columns) - available_columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in {source_name}: {', '.join(missing_columns)}"
        )


# ############################################################################
# Functions to load data from phenodigm-cache into db


def load_one_ontology_cache_file(config, onto_pair):
    """Stream one Parquet ontology-mapping cache into the database."""

    _, _, _, dbdir = pd2tools.getPD2dirs(config)
    cache_file = f"phenio-cache-{onto_pair}.parquet"
    cache_path = Path(dbdir, cache_file)
    pd2tools.log(f"Loading: {cache_file}", 2)

    if not cache_path.is_file():
        raise FileNotFoundError(f"Ontology mapping cache not found: {cache_path}")

    source = pl.scan_parquet(cache_path)
    _validate_required_columns(source, CACHE_COLUMNS, cache_file)

    # Replaces old ModelOntologyOntologyMapping.addData() logic where: Accept insertions of type A, B, data, But do not acccept of type B, A, data.
    normalized_columns = [
        pl.col(column).cast(data_type).str.strip_chars()
        if data_type == pl.String
        else pl.col(column).cast(data_type)
        for column, data_type in CACHE_SCHEMA.items()
    ]
    canonical = source.select(normalized_columns).filter(
        pl.col(MATCH_COLUMN) >= pl.col(QUERY_COLUMN)
    )
    reversed_rows = canonical.filter(
        pl.col(QUERY_COLUMN) != pl.col(MATCH_COLUMN)
    ).select(
        pl.col(MATCH_COLUMN).alias(QUERY_COLUMN),
        pl.col(QUERY_COLUMN).alias(MATCH_COLUMN),
        SIMILARITY_COLUMN,
        INFORMATION_CONTENT_COLUMN,
        LCS_COLUMN,
    )
    # Concat both datasets
    rows = pl.concat((canonical, reversed_rows))

    # Init model and write a generator containing batches of rows from the df
    mapdata = pd2models.ModelOntologyOntologyMapping(config.dbfile)
    batches = (
        batch.iter_rows() for batch in rows.collect_batches(chunk_size=mapdata.insertN)
    )
    mapdata.save_batches(batches)


def load_ontology_mapping_cache_files(config, ontology_pairs):
    """Transfer Parquet ontology-mapping caches into the database."""

    # This action rebuilds ontology_ontology_mapping from scratch. Clear it once
    # up front: transform skips existing caches, but the load below always
    # INSERTs, so without this a re-run appends a second copy of every mapping.
    mapdata = pd2models.ModelOntologyOntologyMapping(config.dbfile)
    mapdata.emptyTable()
    pd2tools.log(f"Cleared existing rows from {mapdata.tabname}", 2)

    for onto_pair in ontology_pairs:
        load_one_ontology_cache_file(config, onto_pair)


def transform_one_ontology_mapping_file(config, onto_pair):
    """Transform one archived Phenio TSV into a filtered Parquet cache."""

    _, datadir, _, dbdir = pd2tools.getPD2dirs(config)

    prefix = "phenio-cache-"
    suffix = ".parquet"

    # Makeshift name matching:
    file_matching = {
        "hp-hp": "HP_vs_HP_semsimian_phenio.tsv.tar.gz",
        "hp-mp": "HP_vs_MP_semsimian_phenio.tsv.tar.gz",
    }

    # Validate before indexing the mapping so callers see the supported pairs
    # instead of an opaque KeyError.
    if onto_pair not in file_matching:
        expected_pairs = ", ".join(sorted(file_matching))
        raise ValueError(
            f"Unknown ontology pair {onto_pair!r}; expected one of: {expected_pairs}"
        )

    input_filename = file_matching[onto_pair]
    input_path = Path(datadir, "obo", input_filename)
    output_file = Path(dbdir, prefix + onto_pair + suffix)
    # Skip if cache file exists
    if output_file.exists():
        pd2tools.log(f"Skipping: {output_file.name}", 2)
        return

    # The source archive is only required when the cache must be regenerated.
    if not input_path.is_file():
        raise FileNotFoundError(f"Ontology mapping archive not found: {input_path}")

    pd2tools.log(f"Generating: {output_file.name}", 2)

    with tarfile.open(input_path, "r:gz") as archive:
        tsv_members = [
            member
            for member in archive.getmembers()
            if member.isfile()
            and "_phenio_" in Path(member.name).name
            and member.name.endswith(".tsv")
        ]
        # A release archive must contain exactly one regular Phenio TSV. Separate
        # errors make unexpected archive layouts straightforward to diagnose.
        if not tsv_members:
            raise ValueError(
                f"No regular _phenio_*.tsv file found in archive: {input_path}"
            )
        if len(tsv_members) > 1:
            member_names = ", ".join(member.name for member in tsv_members)
            raise ValueError(
                f"Multiple _phenio_*.tsv files found in archive {input_path}: "
                f"{member_names}"
            )

        tsv_stream = archive.extractfile(tsv_members[0])
        # Defend against tarfile being unable to expose an otherwise valid member.
        if tsv_stream is None:
            raise ValueError(f"Could not read {tsv_members[0].name} from {input_path}")

        # Apply transformations to the ontology tsv file with Polars
        with tsv_stream:
            # Keep raw fields as strings until the header is validated. This
            # prevents identifier coercion and lets us report all missing columns
            # before Polars attempts to resolve the transformation expressions.
            source = pl.scan_csv(
                tsv_stream,
                separator="\t",
                infer_schema=False,
            )
            _validate_required_columns(
                source,
                PHENIO_COLUMNS,
                tsv_members[0].name,
            )

            similarity_column = CACHE_TO_PHENIO_COLUMNS[SIMILARITY_COLUMN]
            information_content_column = CACHE_TO_PHENIO_COLUMNS[
                INFORMATION_CONTENT_COLUMN
            ]
            lcs_column = CACHE_TO_PHENIO_COLUMNS[LCS_COLUMN]

            transformed = (
                # Only score columns need numeric types for filtering and output.
                # Cast non-strict so a single malformed score becomes null and is
                # dropped below, instead of aborting the whole file.
                source.with_columns(
                    pl.col(similarity_column).cast(
                        CACHE_SCHEMA[SIMILARITY_COLUMN], strict=False
                    ),
                    pl.col(information_content_column).cast(
                        CACHE_SCHEMA[INFORMATION_CONTENT_COLUMN], strict=False
                    ),
                )
                .filter(
                    pl.col(similarity_column).is_not_null()
                    & (
                        pl.col(information_content_column)
                        >= float(config.ontology_mapping_min_ic)
                    )
                )
                .with_columns(
                    (
                        pl.col(lcs_column)
                        .str.replace_all(":", "_", literal=True)
                        .str.strip_chars()
                        + pl.lit(";")
                    ).alias(lcs_column)
                )
                .select(PHENIO_COLUMNS)
                .rename(PHENIO_TO_CACHE_COLUMNS)
            )
            transformed.sink_parquet(output_file)


def transform_ontology_mapping_files(config, ontology_pairs):
    """Transform and filter Phenio ontology-mapping files."""

    for onto_pair in ontology_pairs:
        transform_one_ontology_mapping_file(config, onto_pair)


# ############################################################################
# run this from outside the module


def run_ontology_mapping_processing(config):
    """Transform and load the ontology files to produce ontology-ontology scores."""

    pd2tools.log("Running ontology mapping processing")

    # Define ontology pairs:
    ontology_pairs = (
        "hp-hp",
        "hp-mp",
    )

    # transform ontology_mapping files
    transform_ontology_mapping_files(config, ontology_pairs)
    pd2tools.log("Transferring ontology mapping data into db")
    load_ontology_mapping_cache_files(config, ontology_pairs)
    pd2tools.log("Complete")
