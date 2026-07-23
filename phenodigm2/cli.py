"""Command-line interface for PhenoDigm2.

Provides an importable ``main()`` entry point (exposed as the ``phenodigm2``
console script) that parses arguments and dispatches to the pipeline stages.

@author: Tomasz Konopka
"""

import argparse
from sys import exit

from . import load
from . import dbbuild
from . import export
from . import ontology_mapping
from . import parquet as parquet_export
from . import prep
from . import query
from . import score
from . import solr
from . import status
from . import tools
from . import post_process


# ############################################################################
# create a parser object for the phenodigm executable


def build_parser():
    """Build the argparse parser for the phenodigm2 executable."""

    parser = argparse.ArgumentParser(description="PhenoDigm2")

    # for outputs
    parser.add_argument(
        "--db",
        action="store",
        help="Database directory",
    )

    # for internal tuning/thresholding
    parser.add_argument(
        "--phenodigm_min_perc",
        action="store",
        help="minimal perc value for phenodigm score",
        default=60,
    )
    parser.add_argument(
        "--phenodigm_min_raw",
        action="store",
        help="minimal raw value for phenodigm score",
        default=1.3,
    )
    parser.add_argument(
        "--impc_pval",
        action="store",
        help="p-value threshold for IMPC statistical tests",
        default=1e-4,
    )
    parser.add_argument(
        "--cores",
        action="store",
        help="number of cores used in phenodigm scoring",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="only compute disease-model associations",
        default=False,
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write extra output to stdout",
        default=False,
    )

    # for "export"
    parser.add_argument(
        "--table", action="store", help="For use with --export", default=""
    )
    parser.add_argument(
        "--where", action="store", help="for use with --export", default=""
    )

    # for "parquet"
    parser.add_argument(
        "--parquet-dir",
        action="store",
        help=(
            "Parquet bundle directory "
            f"(default: <db>/{parquet_export.DEFAULT_OUTPUT_DIRECTORY})"
        ),
        default=None,
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing Parquet bundle",
        default=False,
    )

    # for "explain"
    parser.add_argument(
        "--model",
        action="store",
        help="used with --query and --compute, a model id",
        default="",
    )
    parser.add_argument(
        "--disease",
        action="store",
        help="used with --query and --compute, a disease id",
        default="",
    )
    parser.add_argument(
        "--gene", action="store", help="used with --explain, a gene id", default=""
    )
    parser.add_argument(
        "--term", action="store", help="used with --query, an ontology term", default=""
    )
    parser.add_argument(
        "--sim",
        action="store_true",
        help="used with --query, show score data",
        default=False,
    )
    parser.add_argument(
        "--phenotype",
        action="store_true",
        help="used with --query, show phenotypes",
        default=False,
    )
    parser.add_argument(
        "--association",
        action="store_true",
        help="used with --query, show phenodigm associations",
        default=False,
    )
    parser.add_argument(
        "--score",
        action="store_true",
        help="used --query, show phenodigm calculations",
        default=False,
    )

    # processing Phenio mappings to obtain ontology term-term similarities
    parser.add_argument(
        "--ontology_mapping",
        action="store",
        help="complete command to TL ontology mapping (phenio files)",
        default="ontology_mapping",
    )
    parser.add_argument(
        "--ontology_mapping_min_ic",
        action="store",
        help="minimal ic for ontology_mapping output",
        default=2.5,
    )

    # creating a Solr core and configuring its server addresses
    parser.add_argument(
        "--solr_url",
        action="store",
        help="url for communicating with solr",
        default="http://localhost:8983/solr/",
    )
    parser.add_argument(
        "--solr_cores_dir",
        action="store",
        help=(
            "directory holding localhost Solr cores "
            "(default: <db>/output/solr/solrcores7.5)"
        ),
        default=None,
    )
    parser.add_argument(
        "--solr_corename",
        action="store",
        help="name of core set in core.properties",
        default="phenodigm",
    )

    # These shared output thresholds replace the legacy Solr-only names:
    # solr_min_mapscore -> output_min_ontology_ontology_score
    # solr_min_2dscore  -> output_min_disease_model_2d_score
    parser.add_argument(
        "--output_min_ontology_ontology_score",
        action="store",
        type=float,
        help="minimum sqrt(simJ * ic) ontology-ontology mapping score",
        default=1.5,
    )
    parser.add_argument(
        "--output_min_disease_model_2d_score",
        action="store",
        type=float,
        help="minimum 2D raw score for computed disease-model associations",
        default=2.2,
    )

    # Initialization is a standalone mode, not a pipeline stage. Keeping it
    # mutually exclusive with the positional action prevents an accidental
    # download or build while a user is only preparing a directory.
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument(
        "--init",
        action="store",
        metavar="BUNDLE_NAME",
        help="Create a new build directory and seed its bundled resources",
    )
    action_group.add_argument(
        "action",
        action="store",
        nargs="?",
        help="Type of calculation/action to perform",
        choices=[
            "download",
            "build",
            "ontology_mapping",
            "score",
            "index",
            "solr-prepare",
            "solr",
            "parquet",
            "query",
            "export",
            "compute",
            "status",
            "post-process",
        ],
    )

    return parser


# ############################################################################
# Execute the program


def main(argv=None):
    """Parse arguments and run the requested PhenoDigm2 action."""

    parser = build_parser()
    config = parser.parse_args(argv)

    if config.init is not None:
        if config.db is not None:
            parser.error(
                "--db cannot be used with --init; pass the bundle path directly"
            )
        config.db = config.init
        tools.log("Starting PhenoDigm2 [init]")
        try:
            prep.runDirInit(config)
        except FileExistsError as error:
            exit(f"Error initializing build directory: {error}")
        tools.log("Done")
        return

    try:
        config.dbfile = tools.getDBfile(config)
    except Exception:
        exit("Error creating connection to db (check --db)")

    # handle requests for manual inspection
    # and one-off tasks that don't require start/done messages
    if config.action == "query":
        query.queryPhenodigm(config)
        exit()
    if config.action == "compute":
        query.computeScores(config)
        exit()
    if config.action == "export":
        export.exportTables(config)
        exit()
    if config.action == "parquet":
        tools.log("Starting PhenoDigm2 [parquet]")
        parquet_export.runParquetBundleBuild(config)
        tools.log("Done")
        return

    tools.log("Starting PhenoDigm2 [" + config.action + "]")

    if config.action == "status":
        status.statusPhenodigm(config)
        exit()

    prep.runDirPrep(config)
    if config.action == "download":
        prep.runDownloads(config)
    if config.action == "build":
        dbbuild.runDBBuild(config)
        load.runLoad(config)
    if config.action == "ontology_mapping":
        ontology_mapping.run_ontology_mapping_processing(config)
    if config.action == "score":
        score.runScoring(config)
    if config.action == "index":
        dbbuild.runDBIndexing(config)
    if config.action == "solr-prepare":
        solr.prepareSolrBundle(config)
    if config.action == "solr":
        solr.runSolrCoreBuild(config)
    if config.action == "post-process":
        post_process.run_post_process(config)

    tools.log("Done")


if __name__ == "__main__":
    main()
