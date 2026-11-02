"""Transform and load Phenio ontology-mapping files.

@author: Diego Pava
Based on code by Tomasz Konopka

"""

from pathlib import Path
import tarfile

import polars as pl

from . import tools as pd2tools
from . import dbmodels as pd2models

# Phenio file handles for columns
SUBJECT_ID, OBJECT_ID = "subject_id", "object_id"
JACCARD_SIMILARITY = "jaccard_similarity"
ANCESTOR_INFORMATION_CONTENT = "ancestor_information_content"
PHENODIGM_SCORE, ANCESTOR_ID = "phenodigm_score", "ancestor_id"

# Cache file handles for the columns the transform/load logic names directly.
QUERY, MATCH = "query", "match"
PD, SIMJ, IC, LCS = "pd", "simJ", "ic", "lcs"

PHENIO_SCHEMA = pl.Schema(
    {
        SUBJECT_ID: pl.String,
        OBJECT_ID: pl.String,
        JACCARD_SIMILARITY: pl.Float64,
        ANCESTOR_INFORMATION_CONTENT: pl.Float64,
        PHENODIGM_SCORE: pl.Float64,
        ANCESTOR_ID: pl.String,
    }
)

CACHE_SCHEMA = pl.Schema(
    {
        QUERY: pl.String,
        MATCH: pl.String,
        SIMJ: pl.Float64,
        IC: pl.Float64,
        PD: pl.Float64,
        LCS: pl.String,
    }
)

EQUIVALENCE_SCHEMA = {
    SUBJECT_ID: QUERY,
    OBJECT_ID: MATCH,
    JACCARD_SIMILARITY: SIMJ,
    ANCESTOR_INFORMATION_CONTENT: IC,
    PHENODIGM_SCORE: PD,
    ANCESTOR_ID: LCS,
}


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

    minimum_ic = float(config.ontology_mapping_min_ic)
    minimum_phenodigm_score = float(config.output_min_ontology_ontology_score)
    pd2tools.log(
        f"Generating: {output_file.name} "
        f"(ontology_mapping_min_ic={minimum_ic:g}, "
        "output_min_ontology_ontology_score="
        f"{minimum_phenodigm_score:g})",
        2,
    )

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
            # Validate the phenio schema
            _validate_required_columns(
                source,
                PHENIO_SCHEMA.names(),
                tsv_members[0].name,
            )

            transformed = (
                # Only score columns need numeric types for filtering and output.
                # Cast non-strict so a single malformed score becomes null and is
                # dropped below, instead of aborting the whole file.
                source.with_columns(
                    pl.col(PHENODIGM_SCORE).cast(
                        PHENIO_SCHEMA[PHENODIGM_SCORE], strict=False
                    ),
                    pl.col(JACCARD_SIMILARITY).cast(
                        PHENIO_SCHEMA[JACCARD_SIMILARITY], strict=False
                    ),
                    pl.col(ANCESTOR_INFORMATION_CONTENT).cast(
                        PHENIO_SCHEMA[ANCESTOR_INFORMATION_CONTENT], strict=False
                    ),
                )
                .filter(
                    pl.col(JACCARD_SIMILARITY).is_not_null()
                    & (pl.col(ANCESTOR_INFORMATION_CONTENT) >= minimum_ic)
                )
                .filter(
                    pl.col(PHENODIGM_SCORE).is_not_null()
                    & (pl.col(PHENODIGM_SCORE) >= minimum_phenodigm_score)
                )
                .with_columns(
                    (
                        pl.col(ANCESTOR_ID)
                        .str.replace_all(":", "_", literal=True)
                        .str.strip_chars()
                        + pl.lit(";")
                    ).alias(ANCESTOR_ID)
                )
                .rename(EQUIVALENCE_SCHEMA)
                .select(CACHE_SCHEMA.names())
            )
            transformed.sink_parquet(output_file)


def load_one_ontology_cache_file(config, onto_pair):
    """Stream one Parquet ontology-mapping cache into the database."""

    _, _, _, dbdir = pd2tools.getPD2dirs(config)
    cache_file = f"phenio-cache-{onto_pair}.parquet"
    cache_path = Path(dbdir, cache_file)
    pd2tools.log(f"Loading: {cache_file}", 2)

    if not cache_path.is_file():
        raise FileNotFoundError(f"Ontology mapping cache not found: {cache_path}")

    source = pl.scan_parquet(cache_path)
    _validate_required_columns(source, CACHE_SCHEMA.names(), cache_file)

    # Replaces old ModelOntologyOntologyMapping.addData() logic where: Accept insertions of type A, B, data, But do not acccept of type B, A, data.
    normalized_columns = [
        pl.col(column).cast(data_type).str.strip_chars()
        if data_type == pl.String
        else pl.col(column).cast(data_type)
        for column, data_type in CACHE_SCHEMA.items()
    ]

    # DROP PD - not needed for db
    canonical = (
        source.select(normalized_columns)
        .drop(PD)
        .filter(pl.col(MATCH) >= pl.col(QUERY))
    )
    reversed_rows = canonical.filter(pl.col(QUERY) != pl.col(MATCH)).select(
        pl.col(MATCH).alias(QUERY),
        pl.col(QUERY).alias(MATCH),
        SIMJ,
        IC,
        LCS,
    )
    # Concat both datasets
    rows = pl.concat((canonical, reversed_rows))

    # Init model and write a generator containing batches of rows from the df.
    # LazyFrame.collect_batches streams the query in row-count chunks, but it was
    # only added in newer polars. On older polars (>=1.32.3) fall back to
    # collecting once and slicing, which materializes this (IC-filtered) frame.
    mapdata = pd2models.ModelOntologyOntologyMapping(config.dbfile)
    if hasattr(rows, "collect_batches"):
        frames = rows.collect_batches(chunk_size=mapdata.insertN)
    else:
        frames = rows.collect().iter_slices(mapdata.insertN)
    batches = (batch.iter_rows() for batch in frames)
    mapdata.save_batches(batches)


def _validate_required_columns(source, required_columns, source_name):
    """Raise a clear error when a lazy data source lacks required columns."""

    available_columns = set(source.collect_schema().names())
    missing_columns = sorted(set(required_columns) - available_columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in {source_name}: {', '.join(missing_columns)}"
        )


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
