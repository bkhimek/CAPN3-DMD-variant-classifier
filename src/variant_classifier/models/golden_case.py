"""GoldenCase — an expected result, curated independently of the evidence
bundle it judges. Golden cases live in validation/golden_cases/, deliberately
separate from data/curated/, so the test suite never validates the engine
against its own assumptions (see Validation and Verification Design Guide,
"Golden-case philosophy").
"""

from dataclasses import dataclass
from typing import Dict, Optional

from ..errors import SchemaValidationError
from ._coerce import optional_str, require_dict, require_str
from .enums import ACMG_CRITERION_CODES, CriterionStatus, ProvisionalClass


@dataclass(frozen=True)
class GoldenCase:
    variant_id: str
    expected_provisional_class: ProvisionalClass
    expected_criterion_status: Dict[str, CriterionStatus]
    source: str
    curator_note: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"GoldenCase[{self.variant_id}]"
        if not self.expected_criterion_status:
            raise SchemaValidationError(f"{context}: expected_criterion_status must not be empty")
        for code in self.expected_criterion_status:
            if code not in ACMG_CRITERION_CODES:
                raise SchemaValidationError(f"{context}: '{code}' is not a recognised ACMG/AMP criterion code")
        if not self.source.strip():
            raise SchemaValidationError(f"{context}: source must be a non-empty string (where did this expectation come from?)")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "GoldenCase":
        data = require_dict(data, context or "GoldenCase")
        variant_id = require_str(data, "variant_id", context or "GoldenCase")
        ctx = f"GoldenCase[{variant_id}]"
        expected_class_raw = require_str(data, "expected_provisional_class", ctx)
        try:
            expected_provisional_class = ProvisionalClass(expected_class_raw)
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in ProvisionalClass))
            raise SchemaValidationError(
                f"{ctx}: 'expected_provisional_class' = {expected_class_raw!r} is invalid; expected one of: {valid}"
            ) from exc
        raw_status_map = data.get("expected_criterion_status")
        if not isinstance(raw_status_map, dict) or not raw_status_map:
            raise SchemaValidationError(f"{ctx}: 'expected_criterion_status' must be a non-empty mapping of code -> status")
        expected_criterion_status = {}
        for code, status_raw in raw_status_map.items():
            try:
                expected_criterion_status[code] = CriterionStatus(status_raw)
            except ValueError as exc:
                valid = ", ".join(sorted(v.value for v in CriterionStatus))
                raise SchemaValidationError(
                    f"{ctx}: expected_criterion_status[{code!r}] = {status_raw!r} is invalid; expected one of: {valid}"
                ) from exc
        source = require_str(data, "source", ctx)
        curator_note = optional_str(data, "curator_note")
        return cls(
            variant_id=variant_id,
            expected_provisional_class=expected_provisional_class,
            expected_criterion_status=expected_criterion_status,
            source=source,
            curator_note=curator_note,
        )
