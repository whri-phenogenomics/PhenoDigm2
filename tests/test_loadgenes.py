import gzip

from phenodigm2.loadgenes import _write_ensembl_orthologs


def _write_gzip_lines(path, lines):
    with gzip.open(path, "wt", encoding="utf-8") as file:
        file.write("\n".join(lines) + "\n")


def test_write_ensembl_orthologs_writes_gzipped_tsv(tmp_path):
    annotations = tmp_path / "annotations"
    headers = tmp_path / "headers"
    annotations.mkdir()
    headers.mkdir()

    (headers / "human_genes_ensembl.header").write_text(
        "human_ensembl_gene_id\thgnc_id\thgnc_symbol\tchromosome_name\t"
        "start_position\tend_position\tband\n"
    )
    (headers / "human_to_homolog_ensembl.header").write_text(
        "human_ensembl_gene_id\tmouse_ensembl_gene_id\t"
        "mmusculus_homolog_associated_gene_name\t"
        "mmusculus_homolog_orthology_type\n"
    )
    (headers / "mouse_genes_ensembl.header").write_text(
        "mouse_ensembl_gene_id\tmgi_id\tmgi_symbol\n"
    )
    (headers / "human_mouse_mapping.header").write_text(
        "hgnc_id\thgnc_symbol\tchromosome_name\tstart_position\t"
        "end_position\tband\tmgi_symbol\tmgi_id\n"
    )

    _write_gzip_lines(
        annotations / "human_genes_ensembl.txt.gz",
        ["ENSG1\tHGNC:1\tHUMAN1\t1\t100\t200\tq1"],
    )
    _write_gzip_lines(
        annotations / "human_to_homolog_ensembl.txt.gz",
        ["ENSG1\tENSMUSG1\tMOUSE1\tortholog_one2one"],
    )
    _write_gzip_lines(
        annotations / "mouse_genes_ensembl.txt.gz",
        ["ENSMUSG1\tMGI:1\tMOUSE1"],
    )

    _write_ensembl_orthologs(annotations, headers)

    with gzip.open(
        annotations / "human_mouse_mapping.txt.gz", "rt", encoding="utf-8"
    ) as mapping_file:
        assert mapping_file.read() == (
            "HGNC:1\tHUMAN1\t1\t100\t200\tq1\tMOUSE1\tMGI:1\n"
        )
