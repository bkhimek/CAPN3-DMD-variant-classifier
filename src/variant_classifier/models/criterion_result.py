"""CriterionResult — one evaluated ACMG/AMP criterion, as a structured,
explainable record. This is the central design constraint carried through
every guide in the set: an evidence code is never collapsed into a bare
label. See ACMG Engine Detailed Design Guide, Section 5, "Criterion result
model".
"""

from dataclasses import dataclass, field
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_dict, require_list, require_str
from .enums import ACMG_CRITERION_CODES, AutomationConfidence, CriterionStatus, CriterionStrength, EvidenceDirection


@dataclass(frozen=True)
class CriterionResult:
    code: str
    status: CriterionStatus
    direction: EvidenceDirection
    rule_source: str
    rule_version: str
    rationale: str
    strength: Optional[CriterionStrength] = None
    evidence_ids: List[str] = field(default_factory=list)
    automation: Optional[AutomationConfidence] = None

    def __post_init__(self) -> None:
        context = f"CriterionResult[{self.code}]"
        if self.code not in ACMG_CRITERION_CODES:
            raise SchemaValidationError(
                f"{context}: '{self.code}' is not a recognised ACMG/AMP criterion code"
            )
        if self.status == CriterionStatus.MET and self.strength is None:
            raise SchemaValidationError(f"{context}: status=MET requires an explicit strength")
        if self.status in (CriterionStatus.NOT_MET, CriterionStatus.NOT_EVALUATED, CriterionStatus.NOT_APPLICABLE) and self.strength is not None:
            raise SchemaValidationError(
                f"{context}: status={self.status.value} must not carry a strength value "
                "(strength only applies once a criterion is MET)"
            )
        if not self.rationale.strip():
            raise SchemaValidationError(f"{context}: rationale must not be empty")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CriterionResult":
        data = require_dict(data, context or "CriterionResult")
        code = require_str(data, "code", context or "CriterionResult")
        ctx = f"CriterionResult[{code}]"
        status = coerce_enum(CriterionStatus, data.get("status"), "status", ctx)
        direction = coerce_enum(EvidenceDirection, data.get("direction"), "direction", ctx)
        rule_source = require_str(data, "rule_source", ctx)
        rule_version = require_str(data, "rule_version", ctx)
        rationale = require_str(data, "rationale", ctx)
        strength_raw = data.get("strength")
        strength = coerce_enum(CriterionStrength, strength_raw, "strength", ctx) if strength_raw is not None else None
        automation_raw = data.get("automation")
        automation = (
            coerce_enum(AutomationConfidence, automation_raw, "automation", ctx) if automation_raw is not None else None
        )
        evidence_ids = require_list(data, "evidence_ids", ctx)
        for eid in evidence_ids:
            if not isinstance(eid, str) or not eid.strip():
                raise SchemaValidationError(f"{ctx}: every entry in evidence_ids must be a non-empty string")
        return cls(
            code=code,
            status=status,
            direction=direction,
            rule_source=rule_source,
            rule_version=rule_version,
            rationale=rationale,
            strength=strength,
            evidence_ids=list(evidence_ids),
            automation=automation,
        )
