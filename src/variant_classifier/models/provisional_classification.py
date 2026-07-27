"""ProvisionalClassification — the automated engine's output. Always
PROVISIONAL_AUTOMATED at the point this model is constructed; only a
Scientist Review & Sign-off action (outside Milestone-1 scope) can move a
case to FINAL. See Workflow Architecture Guide v7, Reasoning layer, and
ACMG Engine Detailed Design Guide, Section 6.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import optional_str, require_dict, require_list, require_str
from .criterion_result import CriterionResult
from .enums import ClassificationStatus, ProvisionalClass


@dataclass(frozen=True)
class ProvisionalClassification:
    provisional_class: ProvisionalClass
    status: ClassificationStatus
    criteria: List[CriterionResult]
    combining_rule_source: str
    combining_rule_version: str
    rationale: str
    conflicting_evidence_flag: bool = False
    manual_review_required: bool = False
    # Net Bayesian point total (Tavtigian et al. 2020), when this result came from
    # bayesian.combine_bayesian() rather than engine.py's classic Table 5 combine().
    # None for Table 5 results -- points aren't a concept Table 5 uses. Added batch 20
    # (Milestone 5) alongside bayesian.py; kept optional rather than added as a new
    # required field so every pre-existing Table 5 result and golden case stays valid
    # with no migration needed.
    points: Optional[int] = None

    def __post_init__(self) -> None:
        if not self.criteria:
            raise SchemaValidationError("ProvisionalClassification: criteria list must not be empty")
        if not self.rationale.strip():
            raise SchemaValidationError("ProvisionalClassification: rationale must not be empty")
        if self.status != ClassificationStatus.PROVISIONAL_AUTOMATED:
            raise SchemaValidationError(
                "ProvisionalClassification: Milestone 1 only produces PROVISIONAL_AUTOMATED results "
                "(FINAL requires Scientist Review & Sign-off, which is out of scope here)"
            )
        codes_seen = [c.code for c in self.criteria]
        if len(codes_seen) != len(set(codes_seen)):
            raise SchemaValidationError(
                "ProvisionalClassification: duplicate criterion codes in the same result "
                "(one CriterionResult per code, per the double-counting-avoidance constraint)"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "ProvisionalClassification":
        data = require_dict(data, context or "ProvisionalClassification")
        ctx = context or "ProvisionalClassification"
        provisional_class_raw = require_str(data, "provisional_class", ctx)
        try:
            provisional_class = ProvisionalClass(provisional_class_raw)
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in ProvisionalClass))
            raise SchemaValidationError(
                f"{ctx}: 'provisional_class' = {provisional_class_raw!r} is invalid; expected one of: {valid}"
            ) from exc
        status_raw = require_str(data, "status", ctx)
        try:
            status = ClassificationStatus(status_raw)
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in ClassificationStatus))
            raise SchemaValidationError(
                f"{ctx}: 'status' = {status_raw!r} is invalid; expected one of: {valid}"
            ) from exc
        raw_criteria = require_list(data, "criteria", ctx)
        criteria = [CriterionResult.from_dict(c, f"{ctx}.criteria[{i}]") for i, c in enumerate(raw_criteria)]
        combining_rule_source = require_str(data, "combining_rule_source", ctx)
        combining_rule_version = require_str(data, "combining_rule_version", ctx)
        rationale = require_str(data, "rationale", ctx)
        conflicting_evidence_flag = bool(data.get("conflicting_evidence_flag", False))
        manual_review_required = bool(data.get("manual_review_required", False))
        points_raw = data.get("points")
        points = int(points_raw) if points_raw is not None else None
        return cls(
            provisional_class=provisional_class,
            status=status,
            criteria=criteria,
            combining_rule_source=combining_rule_source,
            combining_rule_version=combining_rule_version,
            rationale=rationale,
            conflicting_evidence_flag=conflicting_evidence_flag,
            manual_review_required=manual_review_required,
            points=points,
        )
