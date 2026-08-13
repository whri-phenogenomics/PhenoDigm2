"""Runs post processing of Phenodigm release.
This involves prodcing necessary files for data analysis of each DR and the Disease Models ShinyApp

By: Diego Pava"""

import json
import shutil
import subprocess
from importlib.resources import as_file
from pathlib import Path
from typing import Dict

import polars as pl
import requests

from . import tools as pd2tools
from . import post_process_config as pd2PostProcConfig


_SQLITE_TO_POLARS_TYPE = {
    "INT": pl.Int64,
    "INTEGER": pl.Int64,
    "REAL": pl.Float64,
    "TEXT": pl.String,
}


# Contants with paths to download data
def post_process_paths(config) -> Dict[str, Path]:
    post_proc_dir = _get_post_proc_dir(config)
    return {
        "data": Path(post_proc_dir, "data"),
        "data_aux": Path(post_proc_dir, "data", "auxiliary"),
        "hpo": Path(post_proc_dir, "data", "hpo"),
        "impc": Path(post_proc_dir, "data", "impc"),
        "intermediate": Path(post_proc_dir, "data", "intermediate"),
        "output": Path(post_proc_dir, "data", "output"),
        "phenodigm": Path(post_proc_dir, "data", "phenodigm"),
        "scripts": Path(post_proc_dir, "scripts"),
        "scripts_aux": Path(post_proc_dir, "scripts", "auxiliary"),
    }


# helper functions here
def _get_post_proc_dir(config) -> Path:
    rootdir = Path(config.db).resolve()
    return Path(rootdir, "post_processing")


def _get_pipeline_status_dir(config) -> Path:
    post_proc_dir = _get_post_proc_dir(config)
    return Path(post_proc_dir, "pipeline_status")


def _run_post_process_step(
    config,
    marker_name,
    action,
    start_message,
    success_message,
    marker_message,
    failure_message,
):
    """Run one idempotent pipeline step and record its successful completion."""
    marker = _get_pipeline_status_dir(config) / marker_name
    if marker.exists():
        pd2tools.log(f"Skipping completed post-processing step: {marker_name}")
        return

    pd2tools.log(start_message)
    try:
        action(config)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(marker_message, encoding="utf-8")
        pd2tools.log(success_message)
    except Exception as error:
        pd2tools.log(f"{failure_message}: {error}")
        raise RuntimeError(f"{failure_message} - check logs for details") from error


def run_post_process(config):
    """Run the post-processing stages sequentially without an external scheduler."""
    pd2tools.log("Running post processing pipeline")
    _run_post_process_step(
        config,
        ".directories_created",
        create_post_process_dirs,
        "Creating post-processing directories...",
        "Created post-processing directories successfully.",
        "Directories created successfully.",
        "Creating directories failed",
    )
    _run_post_process_step(
        config,
        ".files_downloaded",
        download_resources,
        "Downloading post_process resources...",
        "Downloaded post-processing files successfully.",
        "Files downloaded successfully.",
        "Downloading resources failed",
    )
    _run_post_process_step(
        config,
        ".tables_created",
        export_tables,
        "Exporting tables:",
        "Exported tables successfully.",
        "Tables exported successfully.",
        "Exporting tables failed",
    )
    _run_post_process_step(
        config,
        ".resources_relocated",
        relocate_external_resources,
        "Relocating external resources...",
        "Relocated external resources successfully.",
        "Relocated external resources successfully.",
        "External resources relocation failed",
    )
    _run_post_process_step(
        config,
        ".post_processing_analysis_complete",
        run_post_processing_analysis,
        "Running post processing analysis...",
        "Post processing analysis successful.",
        "Post processing analysis successful.",
        "Running post processing analysis failed",
    )


def download_data(url, filename):
    """Generic function to download data passing a URL.

    Args:
        url (str): URL of the ontology/gwas catalogue
        filename (str): Filename to save the data. Can pass a path.
    """
    pd2tools.log(f"Fetching {filename}")
    response = requests.get(url, timeout=10)
    with open(filename, "wb") as f:
        f.write(response.content)


def downloads_dict(impc_data_release: str = "latest"):
    """Generic dictionary containing urls, filenames and paths to download post_processing files

    Args:
        impc_data_release (str, optional): IMPC data release to fetch. For older DRs type in the format: 'release-21.1' Defaults to 'latest'.

    Returns:
        Dict[str]: Returns a dictionary with the urls, filenames and paths to download post_processing files
    """
    downloads = {
        "orthology": {
            "url": "https://www.gentar.org/orthology-api/api/ortholog/one_to_one/impc/write_to_tsv_file",
            "filename": "one_to_one_orthologs.tsv",
            "targetdir": "data_aux",
        },
        # TODO: Instead of downloading Move from vTODAY/data_raw/annotations/to where it's needed OR create a symlink (preferable)
        "impc_geno_pheno_assertions": {
            "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/genotype-phenotype-assertions-ALL.csv.gz",
            "filename": "genotype-phenotype-assertions-ALL.csv.gz",
            "targetdir": "impc",
        },
        # NOTE: In the future, in airflow, this might be available locally.
        "impc_viability": {
            "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/viability.csv.gz",
            "filename": "viability.csv.gz",
            "targetdir": "impc",
        },
        # TODO: Instead of downloading Move from vTODAY/data_raw/annotations/to where it's needed OR create a symlink (preferable)
        "impc_stat_result": {
            "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/statistical-results-ALL.csv.gz",
            "filename": "statistical-results-ALL.csv.gz",
            "targetdir": "impc",
        },
        "hpo_gene_to_pheno": {
            "url": "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-06-23/genes_to_phenotype.txt",
            "filename": "genes_to_phenotype.txt",
            "targetdir": "hpo",
        },
    }
    return downloads


# Function to extract the tables from the phenodigm database
def _read_table(config, table, table_config, where=""):
    """Read one database table into a typed Polars LazyFrame."""
    fields = [field.split() for field in table_config["fields"]]
    columns = [field[0] for field in fields]
    schema = {
        column: _SQLITE_TO_POLARS_TYPE[sql_type.upper()]
        for column, sql_type, *_ in fields
    }
    query = f"SELECT {', '.join(columns)} FROM {table}"
    if where:
        query += f" WHERE {where}"

    connection = pd2tools.getDBconn(config.dbfile)
    try:
        return pl.read_database(
            query=query, connection=connection, schema_overrides=schema
        ).lazy()
    finally:
        connection.close()


def export_tables(config):
    """Export the post-processing database tables as Parquet files."""
    output_path = post_process_paths(config)["phenodigm"]
    _, _, resources_dir, _ = pd2tools.getPD2dirs(config)
    with Path(resources_dir, "db_schema.json").open(encoding="utf-8") as schema_file:
        db_schema = json.load(schema_file)

    tables = [
        "model",
        "model_genotype",
        "disease",
        "disease_gene_mapping",
        "gene_gene_mapping",
    ]

    # Iterate over tables to export
    pd2tools.log("Exporting standard tables...")
    for table in tables:
        output_file_path = output_path / f"{table}.parquet"

        # Skip if the file already exists
        if output_file_path.exists():
            pd2tools.log(f"Skipping table: {table}.parquet (already exists)")
            continue
        pd2tools.log(f"Exporting table: {table}.parquet")
        frame = _read_table(config, table, db_schema[table])
        frame.sink_parquet(output_file_path, compression="zstd")

    # Extract the conditional tables
    pd2tools.log("Extracting conditional tables...")

    conditions = [
        "query LIKE '%OMIM%' AND match LIKE '%#%'",
        "query LIKE '%ORPHA%' AND match LIKE '%#%'",
        "query LIKE '%DECIPHER%' AND match LIKE '%#%'",
        "query LIKE '%OMIM%' AND match NOT LIKE '%#%'",
        "query LIKE '%ORPHA%' AND match NOT LIKE '%#%'",
        "query LIKE '%DECIPHER%' AND match NOT LIKE '%#%'",
        "query LIKE '%OMIM%'",
    ]

    filenames = [
        "disease_model_association_omim_impc",
        "disease_model_association_orphanet_impc",
        "disease_model_association_decipher_impc",
        "disease_model_association_omim_nonimpc",
        "disease_model_association_orphanet_nonimpc",
        "disease_model_association_decipher_nonimpc",
        "disease_model_association_omim",
    ]

    for cond, filename in zip(conditions, filenames):
        output_file_path = output_path / f"{filename}.parquet"
        # Skip if the file already exists
        if output_file_path.exists():
            pd2tools.log(f"Skipping table: {filename}.parquet (already exists)")
            continue
        pd2tools.log(f"Exporting table: {filename}.parquet")
        frame = _read_table(
            config,
            "disease_model_association",
            db_schema["disease_model_association"],
            where=cond,
        )
        frame.sink_parquet(output_file_path, compression="zstd")


def _copy_input(override, bundled, dest_dir):
    """Copy a post-processing input (omim file or R script) into dest_dir.

    Uses the external override path when the yaml provides one, otherwise the
    copy bundled in the package (extracted to a real path via ``as_file``). The
    file is copied into the --db bundle so it can be edited there without
    touching the installed package.
    """
    if override is not None:
        shutil.copy2(override, dest_dir)
    else:
        with as_file(bundled) as src:
            shutil.copy2(src, dest_dir)


def relocate_external_resources(config):
    rootdir = Path(config.db).resolve()

    # Set the target path to extract using post_process_paths
    output_path = post_process_paths(config)

    # Load overrides from the yaml if present; otherwise use bundled defaults.
    config_file = Path(rootdir, "post_process_config.yaml")
    if config_file.exists():
        post_proc_paths = pd2PostProcConfig.load_config(config_file)
    else:
        post_proc_paths = pd2PostProcConfig.PostProcessConfig()

    # The omim curation file and both R scripts default to the copies bundled in
    # the package; the yaml may override any of them with an external path (e.g.
    # a freshly curated omim file or a locally edited script).

    # TODO: These three now live inside the bundle. No need for config file to find their path
    resources = pd2tools.getBundledResourcesDir()
    rscripts = pd2tools.getBundledRScriptsDir()
    _copy_input(
        post_proc_paths.omim_curation_path,
        resources / "omim_curation.tsv",
        output_path["data_aux"],
    )
    _copy_input(
        post_proc_paths.main_r_script_path,
        rscripts / "post_processing_DM_pipeline.R",
        output_path["scripts"],
    )
    _copy_input(
        post_proc_paths.hgnc_symbol_checker_script_path,
        rscripts / "auxiliary" / "hgnc_symbol_checker.R",
        output_path["scripts_aux"],
    )


def run_post_processing_analysis(config):
    # Locate the script
    r_script_path = Path(
        post_process_paths(config)["scripts"], "post_processing_DM_pipeline.R"
    )

    # Post process path for the script to execute:
    post_proc_dir = _get_post_proc_dir(config)

    command = ["Rscript", str(r_script_path), str(post_proc_dir)]
    subprocess.run(command, text=True, check=True)


def create_post_process_dirs(config):
    """Create all directories required by the post-processing workflow."""
    for directory in post_process_paths(config).values():
        directory.mkdir(parents=True, exist_ok=True)
        pd2tools.log(f"Created directory: {directory}")


def download_resources(config):
    """Download the external inputs required by post-processing."""
    downloads = downloads_dict(impc_data_release="latest")
    download_paths = post_process_paths(config)
    for values in downloads.values():
        target_path = download_paths[values["targetdir"]]
        download_data(values["url"], target_path / values["filename"])

    # The input versions intentionally remain flexible and cannot all be
    # determined when the main database build starts.
