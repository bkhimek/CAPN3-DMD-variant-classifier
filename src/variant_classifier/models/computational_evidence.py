"""ComputationalEvidence — the single calibrated computational-evidence
call that feeds PP3 (pathogenic-supporting) or BP4 (benign-supporting).

Deliberately one record per variant, not one per tool: the ACMG Engine
Detailed Design Guide's "Computational and splicing evidence" section
warns against stacking multiple correlated in-silico tool votes as if they
were independent evidence. A future evaluator is expected to have already
done any tool-combining/calibration upstream and hand this model the
result of that calibration, not raw per-tool scores.

Reuses PopulationRetrievalStatus for retrieval_status: the same
missing-vs-absent-vs-unassessed distinction applies to "was this variant
successfully scored by the calibrated predictor" as it does to population
frequency lookups. ABSENT is not a meaningful state for a predictor (a
tool doesn't observe a variant as "absent" the way a frequency database
does), so it is rejected here even though the enum technically allows it.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_float, require_dict, require_str
from .enums import ComputationalPrediction, PopulationRetrievalStatus


@dataclass(frozen=True)
class ComputationalEvidence:
    tool: str
    tool_version: str
    calibration_source: str
    retrieval_status: PopulationRetrievalStatus
    score: Optional[float] = None
    calibrated_prediction: Optional[ComputationalPrediction] = None

    def __post_init__(self) -> None:
        context = f"ComputationalEvidence[{self.tool} {self.tool_version}]"
        if self.retrieval_status == PopulationRetrievalStatus.ABSENT:
            raise SchemaValidationError(
                f"{context}: retrieval_status=ABSENT is not meaningful for computational evidence "
                "(use NOT_APPLICABLE if the predictor does not apply to this variant class, "
                "or NOT_ASSESSED if it simply was not run)"
            )
        if self.retrieval_status == PopulationRetrievalStatus.OBSERVED and self.calibrated_prediction is None:
            raise SchemaValidationError(
                f"{context}: retrieval_status=OBSERVED requires calibrated_prediction to be populated"
            )
        if self.retrieval_status in (
            PopulationRetrievalStatus.UNAVAILABLE,
            PopulationRetrievalStatus.NOT_ASSESSED,
            PopulationRetrievalStatus.NOT_APPLICABLE,
            PopulationRetrievalStatus.UNKNOWN,
        ) and (self.score is not None or self.calibrated_prediction is not None):
            raise SchemaValidationError(
                f"{context}: retrieval_status={self.retrieval_status.value} but score/calibrated_prediction "
                "are populated — a predictor that was not successfully run cannot also report a result"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "ComputationalEvidence":
        data = require_dict(data, context or "ComputationalEvidence")
        tool = require_str(data, "tool", context or "ComputationalEvidence")
        ctx = f"ComputationalEvidence[{tool}]"
        tool_version = require_str(data, "tool_version", ctx)
        calibration_source = require_str(data, "calibration_source", ctx)
        retrieval_status = coerce_enum(PopulationRetrievalStatus, data.get("retrieval_status"), "retrieval_status", ctx)
        prediction_raw = data.get("calibrated_prediction")
        calibrated_prediction = (
            coerce_enum(ComputationalPrediction, prediction_raw, "calibrated_prediction", ctx)
            if prediction_raw is not None
            else None
        )
        return cls(
            tool=tool,
            tool_version=tool_version,
            calibration_source=calibration_source,
            retrieval_status=retrieval_status,
            score=optional_float(data, "score", ctx, minimum=0.0, maximum=1.0),
            calibrated_prediction=calibrated_prediction,
        )
