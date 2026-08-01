"""CnvProvisionalClassification -- the CNV scoring module's output.
Deliberately parallel to ProvisionalClassification, NOT a reuse of it:
ProvisionalClassification.criteria is typed List[CriterionResult] and its
__post_init__ enforces "no duplicate ACMG codes," neither of which fits a
CNV result's List[CnvCategoryResult]. The five-tier output vocabulary
(ProvisionalClass: PATHOGENIC/LIKELY_PATHOGENIC/VUS/LIKELY_BENIGN/BENIGN)
and ClassificationStatus ARE reused as-is -- both frameworks genuinely
produce the same five-tier ACMG/AMP-style output, so sharing that
vocabulary is accurate, not a convenience shortcut.
"""

from dataclasses import dataclass
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import require_dict, require_list, require_str
from .cnv_category_result import CnvCategoryResult
from .enums import ClassificationStatus, ProvisionalClass


@dataclass(frozen=True)
class CnvProvisionalClassification:
    provisional_class: ProvisionalClass
    status: ClassificationStatus
    categories: List[CnvCategoryResult]
    combining_rule_source: str
    combining_rule_version: str
    rationale: str
    points: float
    manual_review_required: bool = False

    def __post_init__(self) -> None:
        if not self.categories:
            raise SchemaValidationError("CnvProvisionalClassification: categories list must not be empty")
        if not self.rationale.strip():
            raise SchemaValidationError("CnvProvisionalClassification: rationale must not be empty")
        if self.status != ClassificationStatus.PROVISIONAL_AUTOMATED:
            raise SchemaValidationError(
                "CnvProvisionalClassification: this prototype only produces PROVISIONAL_AUTOMATED "
                "results (FINAL requires Scientist Review & Sign-off, out of scope here)"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CnvProvisionalClassification":
        data = require_dict(data, context or "CnvProvisionalClassification")
        ctx = context or "CnvProvisionalClassification"
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
        raw_categories = require_list(data, "categories", ctx)
        categories = [CnvCategoryResult.from_dict(c, f"{ctx}.categories[{i}]") for i, c in enumerate(raw_categories)]
        combining_rule_source = require_str(data, "combining_rule_source", ctx)
        combining_rule_version = require_str(data, "combining_rule_version", ctx)
        rationale = require_str(data, "rationale", ctx)
        points_raw = data.get("points")
        if not isinstance(points_raw, (int, float)) or isinstance(points_raw, bool):
            raise SchemaValidationError(f"{ctx}: 'points' must be a number, got {points_raw!r}")
        manual_review_required = bool(data.get("manual_review_required", False))
        return cls(
            provisional_class=provisional_class,
            status=status,
            categories=categories,
            combining_rule_source=combining_rule_source,
            combining_rule_version=combining_rule_version,
            rationale=rationale,
            points=float(points_raw),
            manual_review_required=manual_review_required,
        )
