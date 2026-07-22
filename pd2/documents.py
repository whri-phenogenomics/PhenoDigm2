"""Shared producers for Solr and Parquet documents."""

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from math import sqrt
from types import MappingProxyType
from typing import Self

from .dbextractors import (
    ModelAssociation,
    ModelDiseaseGeneMapping,
    ModelGeneGeneMapping,
    ModelIdPhenotype,
    ModelModelGenotype,
    ModelOntologyOntologyMapping,
    ModelOntologySynonym,
    PhenodigmSimpleGenerator,
    getAllDiseaseGenes,
    getDbDiseaseMap,
    getDbGeneMap,
    getDbMapSets,
    getDbMouseModelMap,
    getDbOntologyMap,
    ooscore,
)
from .document_definitions import (
    BOOLEAN,
    FLOAT,
    INTEGER,
    STRING,
    STRINGS,
    DatasetSpec,
    Document,
    DocumentBuildConfig,
    DocumentValue,
    define_document,
)
from .dss import MapSets
from .tools import log


@dataclass(frozen=True, slots=True)
class DocumentBuildOptions:
    """Output-neutral thresholds used while producing shared documents."""

    min_ontology_ontology_score: float
    min_disease_model_2d_score: float

    @classmethod
    def from_config(cls, config: DocumentBuildConfig) -> Self:
        """Read document thresholds from the CLI/application configuration."""

        # These filters originally lived in the Solr pipeline as
        # solr_min_mapscore and solr_min_2dscore. They now live here so every
        # output writer applies the same document-selection behavior.
        return cls(
            min_ontology_ontology_score=config.output_min_ontology_ontology_score,
            min_disease_model_2d_score=config.output_min_disease_model_2d_score,
        )


DISEASE_SEARCH_FLAG_DEFAULTS = MappingProxyType(
    {
        "human_curated_gene": False,
        "impc_model_with_curated_gene": False,
        "impc_model_with_computed_association": False,
        "mgi_model_with_curated_gene": False,
        "mgi_model_with_computed_association": False,
    }
)


def get_base_disease_document(
    disease_data: Mapping[str, DocumentValue],
) -> dict[str, DocumentValue]:
    """Create a document with the shared disease fields."""

    return {
        "disease_" + key: disease_data[key]
        for key in ("id", "source", "term", "alts", "classes")
    }


def make_disease_gene_documents(
    basedoc: Document, geneids: Iterable[str], orthologs, genedata
) -> list[Document]:
    """Expand valid curated HGNC genes into human/mouse ortholog documents.

    An HGNC gene produces a human-only document when it has no ortholog entry.
    When known mouse orthologs exist, each HGNC/MGI pair becomes one document.
    Non-HGNC, invalid, and unknown genes are ignored.
    """

    gene_pairs: set[str] = set()
    for gene_id in geneids:
        if gene_id not in genedata or not genedata[gene_id].isValid():
            continue
        if not gene_id.startswith("HGNC"):
            continue

        if orthologs.has(gene_id):
            for ortholog in orthologs.get(gene_id):
                if ortholog in genedata:
                    pair = sorted((gene_id, ortholog))
                    gene_pairs.add(" ".join(pair))
        else:
            gene_pairs.add(gene_id)

    def set_in_document(doc: dict[str, DocumentValue], gene_id: str) -> None:
        if gene_id.startswith("MGI"):
            doc["marker_id"] = gene_id
            doc["marker_symbol"] = genedata[gene_id].official_symbol
            doc["marker_symbols_withdrawn"] = list(genedata[gene_id].symbols)
        elif gene_id.startswith("HGNC"):
            doc["hgnc_gene_id"] = gene_id
            doc["hgnc_gene_symbol"] = genedata[gene_id].official_symbol
            doc["hgnc_gene_symbols_withdrawn"] = list(genedata[gene_id].symbols)
            doc["hgnc_gene_locus"] = genedata[gene_id].locus

    result: list[Document] = []
    for pair in gene_pairs:
        doc = basedoc.copy()
        for gene_id in pair.split(" "):
            set_in_document(doc, gene_id)
        result.append(doc)
    return result


def get_marker_model_counts(dbfile) -> dict[str, int]:
    """Get the number of mouse models associated with each marker."""

    marker_models = getDbMapSets(ModelModelGenotype(dbfile), field_key=1, field_value=0)
    return {
        marker: len(marker_models.getset(marker)) for marker in marker_models.keys()
    }


def score_2d(a: float, b: float) -> float:
    """Combine average/max raw scores using their Euclidean magnitude."""

    return sqrt((a * a) + (b * b))


def is_impc(source: str) -> bool:
    """Return the historical disease-search facet bucket for a model source.

    This intentionally treats every source except the literal ``MGI`` as IMPC,
    including other or previously unseen sources. It is a compatibility rule,
    not a claim about the model's true originating organization.
    """

    return source != "MGI"


def iter_gene_documents(config: DocumentBuildConfig) -> Iterator[Document]:
    """Yield gene documents."""

    genes = getDbGeneMap(config.dbfile)
    for gene_id in genes.keys():
        gene = genes[gene_id]
        if not gene.isValid():
            continue
        prefix = "hgnc_" if gene.id.startswith("HGNC") else ""
        yield {
            prefix + "gene_id": gene.id,
            prefix + "gene_symbol": gene.official_symbol,
            prefix + "gene_symbols_withdrawn": gene.symbols,
            prefix + "gene_locus": gene.locus,
        }


def iter_gene_gene_documents(config: DocumentBuildConfig) -> Iterator[Document]:
    """Yield human-to-mouse orthology documents."""

    generator = PhenodigmSimpleGenerator(ModelGeneGeneMapping(config.dbfile))
    for row in generator.next():
        if row["query"].startswith("HGNC") and row["match"].startswith("MGI"):
            yield {"gene_id": row["match"], "hgnc_gene_id": row["query"]}


def iter_disease_documents(config: DocumentBuildConfig) -> Iterator[Document]:
    """Yield disease documents with expanded phenotype labels."""

    ontologies = getDbOntologyMap(config.dbfile)
    diseases = getDbDiseaseMap(config.dbfile)
    phenotypes = getDbMapSets(ModelIdPhenotype(config.dbfile, "disease_phenotype"))
    for disease_id in diseases.keys():
        doc = get_base_disease_document(diseases[disease_id])
        terms = set()
        if phenotypes.has(disease_id):
            for phenotype_id in phenotypes.get(disease_id):
                terms.add(phenotype_id + " " + ontologies[phenotype_id])
        doc["disease_phenotypes"] = list(terms)
        yield doc


def iter_mouse_model_documents(config: DocumentBuildConfig) -> Iterator[Document]:
    """Yield mouse model documents."""

    dbfile = config.dbfile
    ontologies = getDbOntologyMap(dbfile)
    models = getDbMouseModelMap(dbfile)
    genotypes = getDbMapSets(ModelModelGenotype(dbfile))
    genes = getDbGeneMap(dbfile)
    phenotypes = getDbMapSets(ModelIdPhenotype(dbfile, "model_phenotype"))

    for model_key in models.keys():
        model = models[model_key]
        model_id = model["id"]
        doc = {
            "model_" + field: model[field]
            for field in ("id", "source", "description", "genetic_background")
        }

        if genotypes.has(model_id):
            symbols = set()
            accessions = set()
            for gene_id in genotypes.get(model_id):
                if gene_id in genes:
                    symbols.add(genes[gene_id].official_symbol)
                else:
                    log("Unknown gene id: " + str(gene_id), 4)
                accessions.add(gene_id)
            doc["marker_id"] = " ".join(accessions)
            doc["marker_symbol"] = " ".join(symbols)

        terms = set()
        if phenotypes.has(model_id):
            for phenotype_id in phenotypes.get(model_id):
                if phenotype_id in ontologies:
                    terms.add(phenotype_id + " " + ontologies[phenotype_id])
                else:
                    log("Unknown phenotype: " + str(phenotype_id), 4)
                    terms.add(phenotype_id)
        doc["model_phenotypes"] = list(terms)

        if doc["marker_symbol"] != "":
            yield doc


def iter_ontology_documents(config: DocumentBuildConfig) -> Iterator[Document]:
    """Yield ontology term and synonym documents."""

    ontologies = getDbOntologyMap(config.dbfile)
    synonyms = getDbMapSets(ModelOntologySynonym(config.dbfile))
    for phenotype_id in ontologies.keys():
        term = ontologies[phenotype_id]
        doc = {
            "phenotype_id": phenotype_id,
            "phenotype_term": term,
            "ontology": phenotype_id[:2],
            "phenotype_synonym": [],
        }
        if synonyms.has(phenotype_id):
            for synonym in synonyms.get(phenotype_id):
                if synonym != term:
                    doc["phenotype_synonym"].append(synonym)
        yield doc


def iter_ontology_ontology_documents(
    config: DocumentBuildConfig,
) -> Iterator[Document]:
    """Yield unique, known HP-to-MP mappings meeting the composite threshold.

    Reciprocal rows are collapsed by the historical lexical-order check. Only
    HP/MP pairs with known labels and ``sqrt(simJ * ic)`` at or above the
    configured minimum are retained.
    """

    options = DocumentBuildOptions.from_config(config)
    ontologies = getDbOntologyMap(config.dbfile)
    generator = PhenodigmSimpleGenerator(ModelOntologyOntologyMapping(config.dbfile))
    for row in generator.next():
        if row["match"] >= row["query"]:
            continue
        mp_id, hp_id = row["match"], row["query"]
        if mp_id[:2] != "MP" or hp_id[:2] != "HP":
            hp_id, mp_id = mp_id, hp_id
        if mp_id[:2] != "MP" or hp_id[:2] != "HP":
            continue
        if ooscore(row["simJ"], row["ic"]) < options.min_ontology_ontology_score:
            continue
        if mp_id not in ontologies or hp_id not in ontologies:
            continue
        yield {
            "mp_id": mp_id,
            "hp_id": hp_id,
            "mp_term": ontologies[mp_id],
            "hp_term": ontologies[hp_id],
        }


def iter_disease_gene_documents(
    config: DocumentBuildConfig,
) -> Iterator[Document]:
    """Yield non-MGI curated disease/HGNC genes with mouse ortholog expansion."""

    dbfile = config.dbfile
    gene_data = getDbGeneMap(dbfile)
    orthologs = getDbMapSets(ModelGeneGeneMapping(dbfile))
    diseases = getDbDiseaseMap(dbfile)
    disease_genes = MapSets()
    generator = PhenodigmSimpleGenerator(ModelDiseaseGeneMapping(dbfile))
    for row in generator.next():
        if row["source"] != "MGI":
            disease_genes.add(row["query"], row["match"])

    for disease_id in diseases.keys():
        if not disease_genes.has(disease_id):
            continue
        basedoc = get_base_disease_document(diseases[disease_id])
        yield from make_disease_gene_documents(
            basedoc, disease_genes.get(disease_id), orthologs, gene_data
        )


def iter_disease_model_documents(
    config: DocumentBuildConfig,
) -> Iterator[Document]:
    """Yield disease models passing either computed or curated inclusion.

    Computed inclusion requires the average/max raw 2D score to be strictly
    greater than the configured minimum. Curated inclusion means the model and
    disease gene sets overlap and can retain a model below that score threshold.
    """

    options = DocumentBuildOptions.from_config(config)
    dbfile = config.dbfile
    gene_data = getDbGeneMap(dbfile)
    models = getDbMouseModelMap(dbfile)
    diseases = getDbDiseaseMap(dbfile)
    disease_genes = getAllDiseaseGenes(dbfile)
    model_genes = getDbMapSets(ModelModelGenotype(dbfile))
    model_counts = get_marker_model_counts(dbfile)
    ontologies = getDbOntologyMap(dbfile)

    def phenotype_array(id_string: str) -> list[str]:
        if id_string == "":
            return []
        return [
            phenotype_id + " " + ontologies[phenotype_id]
            for phenotype_id in id_string.split(",")
        ]

    association = ModelAssociation(dbfile)
    association.tabname = "disease_model_association"
    generator = PhenodigmSimpleGenerator(association)
    for row in generator.next():
        disease_id, model_id = row["query"], row["match"]
        genes = {
            gene_id for gene_id in model_genes.getset(model_id) if gene_id in gene_data
        }
        if not genes:
            log("Unknown gene(s) in model: " + model_id, 4)
            continue

        score_passes = (
            score_2d(row["score_avg_raw"], row["score_max_raw"])
            > options.min_disease_model_2d_score
        )
        curated = not genes.isdisjoint(disease_genes.getset(disease_id))
        if not score_passes and not curated:
            continue

        doc = {
            "disease_id": disease_id,
            "disease_term": diseases[disease_id]["term"],
            "model_id": model_id,
        }
        for field in ("source", "description", "genetic_background"):
            doc["model_" + field] = models[model_id][field]
        doc["marker_id"] = " ".join(genes)
        doc["marker_symbol"] = " ".join(
            gene_data[gene_id].official_symbol for gene_id in genes
        )
        doc["marker_locus"] = " ".join(gene_data[gene_id].locus for gene_id in genes)
        doc["marker_num_models"] = sum(model_counts[gene_id] for gene_id in genes)
        doc["association_curated"] = curated
        for score in ("avg_norm", "avg_raw", "max_norm", "max_raw"):
            doc["disease_model_" + score] = row["score_" + score]
        doc["disease_matched_phenotypes"] = phenotype_array(row["query_phenotype"])
        doc["model_matched_phenotypes"] = phenotype_array(row["match_phenotype"])
        yield doc


def iter_disease_search_documents(
    config: DocumentBuildConfig,
) -> Iterator[Document]:
    """Yield disease search rows with curated/computed model facet flags.

    Computed flags use the same strict 2D threshold as disease-model summaries;
    curated flags use disease/model gene overlap. Model sources are assigned to
    the legacy MGI-versus-IMPC buckets implemented by :func:`is_impc`.
    """

    options = DocumentBuildOptions.from_config(config)
    dbfile = config.dbfile
    models = getDbMouseModelMap(dbfile)
    diseases = getDbDiseaseMap(dbfile)
    disease_genes = getAllDiseaseGenes(dbfile)
    model_genes = getDbMapSets(ModelModelGenotype(dbfile))

    search_documents = {}
    for disease_id in diseases.keys():
        doc = get_base_disease_document(diseases[disease_id])
        doc["search_qf"] = [doc["disease_id"], doc["disease_term"]]
        doc["search_qf"].extend(doc["disease_alts"])
        doc.update(DISEASE_SEARCH_FLAG_DEFAULTS)
        search_documents[disease_id] = doc

    for disease_id in disease_genes.keys():
        search_documents[disease_id]["human_curated_gene"] = True

    association = ModelAssociation(dbfile)
    association.tabname = "disease_model_association"
    generator = PhenodigmSimpleGenerator(association)
    for row in generator.next():
        disease_id, model_id = row["query"], row["match"]
        prefix = "impc" if is_impc(models[model_id]["source"]) else "mgi"
        if (
            score_2d(row["score_avg_raw"], row["score_max_raw"])
            > options.min_disease_model_2d_score
        ):
            search_documents[disease_id][
                prefix + "_model_with_computed_association"
            ] = True
        if not model_genes.getset(model_id).isdisjoint(
            disease_genes.getset(disease_id)
        ):
            search_documents[disease_id][prefix + "_model_with_curated_gene"] = True

    yield from search_documents.values()


def _index_dataset_specs(
    specs: tuple[DatasetSpec, ...],
) -> Mapping[str, DatasetSpec]:
    """Build a read-only, insertion-ordered index and reject duplicate names."""

    by_type = {spec.document_type: spec for spec in specs}
    if len(by_type) != len(specs):
        raise ValueError("Dataset document types must be unique")
    return MappingProxyType(by_type)


DATASET_SPECS = _index_dataset_specs(
    (
        DatasetSpec(
            define_document(
                "gene",
                ("gene_id", STRING),
                ("gene_symbol", STRING),
                ("gene_symbols_withdrawn", STRING),
                ("gene_locus", STRING),
                ("hgnc_gene_id", STRING),
                ("hgnc_gene_symbol", STRING),
                ("hgnc_gene_symbols_withdrawn", STRING),
                ("hgnc_gene_locus", STRING),
                ("impc_model", BOOLEAN),
                ("mouse_model", BOOLEAN),
            ),
            iter_gene_documents,
        ),
        DatasetSpec(
            define_document(
                "gene_gene",
                ("gene_id", STRING),
                ("hgnc_gene_id", STRING),
            ),
            iter_gene_gene_documents,
        ),
        DatasetSpec(
            define_document(
                "disease",
                ("disease_id", STRING),
                ("disease_source", STRING),
                ("disease_term", STRING),
                ("disease_alts", STRINGS),
                ("disease_classes", STRINGS),
                ("disease_phenotypes", STRINGS),
            ),
            iter_disease_documents,
        ),
        DatasetSpec(
            define_document(
                "mouse_model",
                ("model_id", STRING),
                ("model_source", STRING),
                ("model_description", STRING),
                ("model_genetic_background", STRING),
                ("marker_id", STRING),
                ("marker_symbol", STRING),
                ("model_phenotypes", STRINGS),
            ),
            iter_mouse_model_documents,
        ),
        DatasetSpec(
            define_document(
                "ontology",
                ("ontology", STRING),
                ("phenotype_id", STRING),
                ("phenotype_term", STRING),
                ("phenotype_synonym", STRINGS),
            ),
            iter_ontology_documents,
        ),
        DatasetSpec(
            define_document(
                "ontology_ontology",
                ("mp_id", STRING),
                ("mp_term", STRING),
                ("hp_id", STRING),
                ("hp_term", STRING),
            ),
            iter_ontology_ontology_documents,
        ),
        DatasetSpec(
            define_document(
                "disease_gene_summary",
                ("disease_id", STRING),
                ("disease_term", STRING),
                ("marker_id", STRING),
                ("marker_symbol", STRING),
                ("marker_symbols_withdrawn", STRINGS),
                ("hgnc_gene_id", STRING),
                ("hgnc_gene_symbol", STRING),
                ("hgnc_gene_symbols_withdrawn", STRINGS),
                ("hgnc_gene_locus", STRING),
            ),
            iter_disease_gene_documents,
        ),
        DatasetSpec(
            define_document(
                "disease_model_summary",
                ("disease_id", STRING),
                ("disease_term", STRING),
                ("model_id", STRING),
                ("model_source", STRING),
                ("model_description", STRING),
                ("model_genetic_background", STRING),
                ("marker_id", STRING),
                ("marker_symbol", STRING),
                ("marker_locus", STRING),
                ("marker_num_models", INTEGER),
                ("disease_model_avg_raw", FLOAT),
                ("disease_model_avg_norm", FLOAT),
                ("disease_model_max_raw", FLOAT),
                ("disease_model_max_norm", FLOAT),
                ("association_curated", BOOLEAN),
                ("disease_matched_phenotypes", STRINGS),
                ("model_matched_phenotypes", STRINGS),
            ),
            iter_disease_model_documents,
        ),
        DatasetSpec(
            define_document(
                "disease_search",
                ("disease_id", STRING),
                ("disease_term", STRING),
                ("disease_source", STRING),
                ("disease_alts", STRINGS),
                ("disease_classes", STRINGS),
                ("search_qf", STRINGS),
                ("human_curated_gene", BOOLEAN),
                ("impc_model_with_curated_gene", BOOLEAN),
                ("impc_model_with_computed_association", BOOLEAN),
                ("mgi_model_with_curated_gene", BOOLEAN),
                ("mgi_model_with_computed_association", BOOLEAN),
            ),
            iter_disease_search_documents,
        ),
    )
)


def get_dataset_spec(document_type: str) -> DatasetSpec:
    """Return the registered schema and producer for a document type."""

    return DATASET_SPECS[document_type]
