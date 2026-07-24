"""CaseInterpretation — the output of Milestone 4's case-level reasoning:
does what was found in this patient explain their disease, given how the
gene's disease is inherited? Deliberately separate from
ProvisionalClassification, which answers a different question ("how
strong is the evidence for this one variant") — CaseInterpretation
consumes one or two ProvisionalClassifications as input rather than
replacing them.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, require_dict, require_str
from .enums import CaseInterpretationStatus


@dataclass(frozen=True)
class CaseInterpretation:
    case_id: str
    gene: str
    status: CaseInterpretationStatus
    rationale: str

    def __post_init__(self) -> None:
        if not self.rationale.strip():
            raise SchemaValidationError(f"CaseInterpretation[{self.case_id}]: rationale must not be empty")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CaseInterpretation":
        data = require_dict(data, context or "CaseInterpretation")
        case_id = require_str(data, "case_id", context or "CaseInterpretation")
        ctx = f"CaseInterpretation[{case_id}]"
        gene = require_str(data, "gene", ctx)
        status = coerce_enum(CaseInterpretationStatus, data.get("status"), "status", ctx)
        rationale = require_str(data, "rationale", ctx)
        return cls(case_id=case_id, gene=gene, status=status, rationale=rationale)
