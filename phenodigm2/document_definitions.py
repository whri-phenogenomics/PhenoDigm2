"""Shared document schema models and normalization rules."""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
from os import PathLike
from typing import Protocol

DocumentValue = str | bool | int | float | list[str] | set[str] | None
NormalizedDocumentValue = str | bool | int | float | list[str] | None
Document = Mapping[str, DocumentValue]
NormalizedDocument = dict[str, NormalizedDocumentValue]


class DocumentBuildConfig(Protocol):
    """Configuration required by the shared document producers."""

    dbfile: str | PathLike[str]
    output_min_ontology_ontology_score: float
    output_min_disease_model_2d_score: float


DocumentProducer = Callable[[DocumentBuildConfig], Iterable[Document]]


class FieldKind(str, Enum):
    """Storage-neutral value types supported by exported documents."""

    STRING = "string"
    BOOLEAN = "boolean"
    INTEGER = "integer"
    FLOAT = "float"
    STRING_LIST = "string_list"


@dataclass(frozen=True, slots=True)
class FieldDefinition:
    """One named field in a document schema."""

    name: str
    kind: FieldKind
    nullable: bool = True


TYPE_FIELD = FieldDefinition("type", FieldKind.STRING, nullable=False)


@dataclass(frozen=True, slots=True)
class DocumentDefinition:
    """The name and ordered fields of one exported document type."""

    document_type: str
    fields: tuple[FieldDefinition, ...]

    @property
    def fieldnames(self) -> tuple[str, ...]:
        """Return declared fields, excluding the logical ``type`` column."""

        return tuple(field.name for field in self.fields)

    @property
    def columns(self) -> tuple[FieldDefinition, ...]:
        """Return every logical column, including the required ``type`` field."""

        return (TYPE_FIELD, *self.fields)

    @property
    def schema(self) -> dict[str, FieldKind]:
        """Return the storage-neutral schema in export column order."""

        return {field.name: field.kind for field in self.columns}


@dataclass(frozen=True, slots=True)
class DatasetSpec:
    """Connect one document definition to the function that produces its rows."""

    definition: DocumentDefinition
    producer: DocumentProducer

    @property
    def document_type(self) -> str:
        """Return the stable name used by Solr and output datasets."""

        return self.definition.document_type


STRING = FieldKind.STRING
BOOLEAN = FieldKind.BOOLEAN
INTEGER = FieldKind.INTEGER
FLOAT = FieldKind.FLOAT
STRINGS = FieldKind.STRING_LIST


def define_document(
    document_type: str, *fields: tuple[str, FieldKind]
) -> DocumentDefinition:
    """Build an immutable document definition while preserving field order."""

    return DocumentDefinition(
        document_type,
        tuple(FieldDefinition(name, kind) for name, kind in fields),
    )


def normalize_document(
    definition: DocumentDefinition, obj: Document
) -> NormalizedDocument:
    """Create an ordered, sparse writer-neutral document.

    Unknown fields are discarded. Declared lists remain lists; sets retain the
    historical Solr representation of a space-joined string, with empty sets
    becoming null. Missing declared fields remain absent at this stage.
    """

    doc = {"type": definition.document_type}
    for field in definition.fieldnames:
        if field not in obj:
            continue
        value = obj[field]
        if type(value) is set:
            value = " ".join(value) if value else None
        doc[field] = value
    return doc


def materialize_row(
    definition: DocumentDefinition, obj: Document
) -> NormalizedDocument:
    """Create a dense row by filling every absent schema column with null."""

    normalized = normalize_document(definition, obj)
    return {name: normalized.get(name) for name in definition.schema}
