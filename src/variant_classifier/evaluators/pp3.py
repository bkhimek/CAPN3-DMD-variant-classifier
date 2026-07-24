"""PP3 evaluator — "computational evidence supports a deleterious effect"
(Richards et al. 2015). Supporting pathogenic evidence.

Reads VariantEvidenceBundle.computational_evidence — the single calibrated
prediction added in Milestone 1 specifically so PP3/BP4 would have
something concrete to evaluate from (see models/computational_evidence.py).
No bundle has computational evidence gathered at all unless that field is
populated; that's a real "we didn't check" state, distinct from "we
checked and it was indeterminate."
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import ComputationalPrediction, CriterionStatus, CriterionStrength, EvidenceDirection, PopulationRetrievalStatus

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_pp3(bundle: VariantEvidenceBundle) -> CriterionResult:
    variant_id = bundle.variant.variant_id
    ce = bundle.computational_evidence

    if ce is None:
        return CriterionResult(
            code="PP3",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"No computational_evidence recorded for {variant_id} — PP3 was never assessed.",
            evidence_ids=[],
        )

    evidence_id = f"{ce.tool}:{ce.tool_version}"

    if ce.retrieval_status != PopulationRetrievalStatus.OBSERVED:
        return CriterionResult(
            code="PP3",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"{evidence_id} has retrieval_status={ce.retrieval_status.value}; PP3 cannot be evaluated without a successful prediction.",
            evidence_ids=[evidence_id],
        )

    if ce.calibrated_prediction == ComputationalPrediction.PATHOGENIC:
        return CriterionResult(
            code="PP3",
            status=CriterionStatus.MET,
            strength=CriterionStrength.SUPPORTING,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"{evidence_id} calibrated prediction is PATHOGENIC (score={ce.score}).",
            evidence_ids=[evidence_id],
        )

    # BENIGN or INDETERMINATE: the pathogenic-supporting condition isn't met either way.
    return CriterionResult(
        code="PP3",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=f"{evidence_id} calibrated prediction is {ce.calibrated_prediction.value} (score={ce.score}), not PATHOGENIC.",
        evidence_ids=[evidence_id],
    )
