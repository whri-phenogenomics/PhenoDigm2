"""Build the PheVal benchmarking outputs without running the R workflow."""

from pathlib import Path

import polars as pl

from . import tools as pd2tools


ALL_MODELS_OUTPUT = "phenodigm_scores_benchmarking_all_models.parquet"
IMPC_MODELS_OUTPUT = "phenodigm_scores_benchmarking_impc_models.parquet"

_MODEL_QUERY = """
    SELECT model.id, model.source, model_genotype.gene_id AS mgi_id
    FROM model
    INNER JOIN model_genotype ON model.id = model_genotype.id
"""

_DISEASE_GENE_QUERY = """
    SELECT disease_gene_mapping.query AS disorder_id,
           disease_gene_mapping.match AS hgnc_id
    FROM disease_gene_mapping
    INNER JOIN disease ON disease_gene_mapping.query = disease.id
    WHERE disease_gene_mapping.source != 'MGI'
"""

_ASSOCIATION_COLUMNS = """
    query, match, score_avg_norm, score_avg_raw, score_max_norm, score_max_raw
"""


def _read_database_inputs(dbfile):
    """Read only the database columns used by the benchmark calculation."""
    connection = pd2tools.getDBconn(dbfile)
    try:
        model = pl.read_database(
            query=_MODEL_QUERY,
            connection=connection,
            schema_overrides={
                "id": pl.String,
                "source": pl.String,
                "mgi_id": pl.String,
            },
        )
        gene_disease = pl.read_database(
            query=_DISEASE_GENE_QUERY,
            connection=connection,
            schema_overrides={"disorder_id": pl.String, "hgnc_id": pl.String},
        )
        all_models = pl.read_database(
            query=(
                f"SELECT {_ASSOCIATION_COLUMNS} "
                "FROM disease_model_association WHERE query LIKE '%OMIM%'"
            ),
            connection=connection,
        )
        impc_models = pl.read_database(
            query=(
                f"SELECT {_ASSOCIATION_COLUMNS} "
                "FROM disease_model_association "
                "WHERE query LIKE '%OMIM%' AND match LIKE '%#%'"
            ),
            connection=connection,
        )
    finally:
        connection.close()

    return model, gene_disease, all_models, impc_models


def _read_orthologs(path):
    """Read and normalize the one-to-one human/mouse ortholog mapping."""
    return pl.read_csv(
        path,
        separator="\t",
        columns=["Human Gene Symbol", "Hgnc Acc Id", "Mgi Gene Acc Id"],
        schema_overrides={
            "Human Gene Symbol": pl.String,
            "Hgnc Acc Id": pl.String,
            "Mgi Gene Acc Id": pl.String,
        },
    ).rename(
        {
            "Human Gene Symbol": "gene_symbol",
            "Hgnc Acc Id": "hgnc_id",
            "Mgi Gene Acc Id": "mgi_id",
        }
    )


def _read_tested_impc_genes(path):
    """Read the distinct IMPC genes present in the statistical results."""
    return (
        pl.scan_csv(
            path,
            schema_overrides={"marker_accession_id": pl.String},
        )
        .select(pl.col("marker_accession_id").alias("mgi_id"))
        .drop_nulls()
        .filter(pl.col("mgi_id") != "-")
        .unique()
        .collect(engine="streaming")
    )


def _top_genes(
    associations,
    model,
    orthologs,
    eligible_diseases,
    *,
    min_percentage,
    min_raw,
):
    """Keep the highest-scoring row for each OMIM disease/HGNC gene pair."""
    score = (pl.col("score_avg_norm") + pl.col("score_max_norm")) / 2
    passes_threshold = (
        (pl.col("score_avg_norm") > min_percentage)
        | (pl.col("score_max_norm") > min_percentage)
        | (pl.col("score_avg_raw") > min_raw)
        | (pl.col("score_max_raw") > min_raw)
    )

    return (
        associations.with_columns(score.alias("score"))
        .filter((pl.col("score") > 0) & passes_threshold)
        .join(model, left_on="match", right_on="id", how="inner")
        .join(orthologs, on="mgi_id", how="inner")
        .join(eligible_diseases, left_on="query", right_on="disorder_id", how="inner")
        .sort("score", descending=True)
        .unique(subset=["query", "hgnc_id"], keep="first", maintain_order=True)
        .select(
            pl.col("query").alias("disorder_id"),
            "gene_symbol",
            "hgnc_id",
            "score",
            "source",
        )
        .unique(maintain_order=True)
    )


def write_benchmark_files(
    *,
    dbfile,
    statistical_results_path,
    orthologs_path,
    output_dir,
    phenodigm_min_perc=60,
    phenodigm_min_raw=1.3,
):
    """Write the two PheVal benchmark Parquet files using Polars."""
    dbfile = Path(dbfile)
    statistical_results_path = Path(statistical_results_path)
    orthologs_path = Path(orthologs_path)
    output_dir = Path(output_dir)

    required_inputs = (dbfile, statistical_results_path, orthologs_path)
    missing_inputs = [path for path in required_inputs if not path.is_file()]
    if missing_inputs:
        missing = ", ".join(str(path) for path in missing_inputs)
        raise FileNotFoundError(f"Missing benchmarking input(s): {missing}")

    model, gene_disease, all_models, impc_models = _read_database_inputs(dbfile)
    orthologs = _read_orthologs(orthologs_path)
    tested_impc_genes = _read_tested_impc_genes(statistical_results_path)

    eligible_diseases = (
        gene_disease.join(orthologs, on="hgnc_id", how="inner")
        .filter(pl.col("disorder_id").str.contains("OMIM"))
        .join(tested_impc_genes, on="mgi_id", how="inner")
        .select("disorder_id")
        .unique()
    )

    options = {
        "model": model,
        "orthologs": orthologs,
        "eligible_diseases": eligible_diseases,
        "min_percentage": float(phenodigm_min_perc),
        "min_raw": float(phenodigm_min_raw),
    }
    all_models_benchmark = _top_genes(all_models, **options)
    impc_models_benchmark = _top_genes(impc_models, **options)

    output_dir.mkdir(parents=True, exist_ok=True)
    all_models_benchmark.write_parquet(
        output_dir / ALL_MODELS_OUTPUT,
        compression="zstd",
        compression_level=5,
    )
    impc_models_benchmark.write_parquet(
        output_dir / IMPC_MODELS_OUTPUT,
        compression="zstd",
        compression_level=5,
    )

    return {
        "all_models": output_dir / ALL_MODELS_OUTPUT,
        "impc_models": output_dir / IMPC_MODELS_OUTPUT,
    }
