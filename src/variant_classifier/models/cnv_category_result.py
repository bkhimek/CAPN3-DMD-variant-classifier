"""CnvCategoryResult -- one evaluated Section-2 (dosage-sensitivity)
category from the ACMG/ClinGen CNV rubric (Riggs et al. 2020), as a
structured, explainable record. Deliberately parallel to CriterionResult,
NOT a reuse of it -- CriterionResult.code is validated against
ACMG_CRITERION_CODES (Richards et al. 2015's 28-code vocabulary), which a
CNV category code like "2A" or "2E" is not a member of, and CriterionResult
carries a CriterionStrength (STRONG/MODERATE/...) rather than a raw point
value. Collapsing the two would either force a fake ACMG code onto a CNV
category or force a fake strength tier onto a point value that has no
tier -- both would misrepresent which framework produced the result. See
cnv_deletion_evidence.py's docstring for the full "why a new evidence
family, not an extension" writeup.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, require_dict, require_list, require_str
from .enums import CNV_LOSS_CATEGORY_CODES, CriterionStatus, EvidenceDirection


@dataclass(frozen=True)
class CnvCategoryResult:
    code: str
    status: CriterionStatus
    direction: EvidenceDirection
    points: float
    rule_source: str
    rule_version: str
    rationale: str
    evidence_ids: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        context = f"CnvCategoryResult[{self.code}]"
        if self.code not in CNV_LOSS_CATEGORY_CODES:
            raise SchemaValidationError(
                f"{context}: '{self.code}' is not a recognised CNV loss category code "
                f"(expected one of {sorted(CNV_LOSS_CATEGORY_CODES)})"
            )
        if self.status == CriterionStatus.MET and self.points == 0.0:
            raise SchemaValidationError(f"{context}: status=MET requires a non-zero points value")
        if self.status != CriterionStatus.MET and self.points != 0.0:
            raise SchemaValidationError(
                f"{context}: status={self.status.value} must carry points=0.0 -- points only apply once "
                "a category is MET (mirrors CriterionResult's MET-requires-strength convention)"
            )
        if not self.rationale.strip():
            raise SchemaValidationError(f"{context}: rationale must not be empty")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CnvCategoryResult":
        data = require_dict(data, context or "CnvCategoryResult")
        code = require_str(data, "code", context or "CnvCategoryResult")
        ctx = f"CnvCategoryResult[{code}]"
        status = coerce_enum(CriterionStatus, data.get("status"), "status", ctx)
        direction = coerce_enum(EvidenceDirection, data.get("direction"), "direction", ctx)
        points_raw = data.get("points")
        if not isinstance(points_raw, (int, float)) or isinstance(points_raw, bool):
            raise SchemaValidationError(f"{ctx}: 'points' must be a number, got {points_raw!r}")
        rule_source = require_str(data, "rule_source", ctx)
        rule_version = require_str(data, "rule_version", ctx)
        rationale = require_str(data, "rationale", ctx)
        evidence_ids = require_list(data, "evidence_ids", ctx)
        for eid in evidence_ids:
            if not isinstance(eid, str) or not eid.strip():
                raise SchemaValidationError(f"{ctx}: every entry in evidence_ids must be a non-empty string")
        return cls(
            code=code,
            status=status,
            direction=direction,
            points=float(points_raw),
            rule_source=rule_source,
            rule_version=rule_version,
            rationale=rationale,
            evidence_ids=list(evidence_ids),
        )
