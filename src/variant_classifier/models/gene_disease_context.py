"""GeneDiseaseContext — the gene/disease-level configuration a criterion
evaluator needs before it can reason about any specific variant: inheritance,
mechanism, and which rule specification applies. See ACMG Engine Detailed
Design Guide, Section 6, "Applicability and rule precedence".

This is deliberately a small config file per gene-disease pair, not
something embedded in code (per the user's original request).
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, require_bool, require_dict, require_str
from .enums import DiseaseMechanism, Inheritance, SpecificationType


@dataclass(frozen=True)
class Specification:
    type: SpecificationType
    version: str

    def __post_init__(self) -> None:
        if not self.version or not self.version.strip():
            raise SchemaValidationError("specification.version must be a non-empty string")

    @classmethod
    def from_dict(cls, data: dict, context: str) -> "Specification":
        data = require_dict(data, f"{context}.specification")
        spec_type = coerce_enum(SpecificationType, data.get("type"), "type", f"{context}.specification")
        version = require_str(data, "version", f"{context}.specification")
        return cls(type=spec_type, version=version)


@dataclass(frozen=True)
class GeneDiseaseContext:
    gene: str
    disease: str
    inheritance: Inheritance
    mechanism: DiseaseMechanism
    lof_established: bool
    specification: Specification
    gene_disease_validity: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"GeneDiseaseContext[{self.gene}]"
        if self.mechanism == DiseaseMechanism.LOSS_OF_FUNCTION and not self.lof_established and self.lof_established is not False:
            # lof_established must be an explicit bool either way; this branch
            # exists so a future reader knows the check was considered, not
            # merely omitted. (bool type is already enforced in from_dict.)
            pass

    @classmethod
    def from_dict(cls, data: dict) -> "GeneDiseaseContext":
        data = require_dict(data, "GeneDiseaseContext")
        gene = require_str(data, "gene", "GeneDiseaseContext")
        context = f"GeneDiseaseContext[{gene}]"
        disease = require_str(data, "disease", context)
        inheritance = coerce_enum(Inheritance, data.get("inheritance"), "inheritance", context)
        mechanism = coerce_enum(DiseaseMechanism, data.get("mechanism"), "mechanism", context)
        lof_established = require_bool(data, "lof_established", context)
        specification = Specification.from_dict(data.get("specification", {}), context)
        gene_disease_validity = data.get("gene_disease_validity")
        return cls(
            gene=gene,
            disease=disease,
            inheritance=inheritance,
            mechanism=mechanism,
            lof_established=lof_established,
            specification=specification,
            gene_disease_validity=gene_disease_validity,
        )
