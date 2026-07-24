"""BP4 evaluator — "computational evidence supports a benign effect"
(Richards et al. 2015). Supporting benign evidence; PP3's mirror image,
reading the same single calibrated ComputationalEvidence record. Added in
Milestone 1 specifically because BS1 alone could never combine to
LIKELY_BENIGN — see SUPPORTED_CRITERIA_MILESTONE_1 in models/enums.py.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import ComputationalPrediction, CriterionStatus, CriterionStrength, EvidenceDirection, PopulationRetrievalStatus

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_bp4(bundle: VariantEvidenceBundle) -> CriterionResult:
    variant_id = bundle.variant.variant_id
    ce = bundle.computational_evidence

    if ce is None:
        return CriterionResult(
            code="BP4",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"No computational_evidence recorded for {variant_id} — BP4 was never assessed.",
            evidence_ids=[],
        )

    evidence_id = f"{ce.tool}:{ce.tool_version}"

    if ce.retrieval_status != PopulationRetrievalStatus.OBSERVED:
        return CriterionResult(
            code="BP4",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"{evidence_id} has retrieval_status={ce.retrieval_status.value}; BP4 cannot be evaluated without a successful prediction.",
            evidence_ids=[evidence_id],
        )

    if ce.calibrated_prediction == ComputationalPrediction.BENIGN:
        return CriterionResult(
            code="BP4",
            status=CriterionStatus.MET,
            strength=CriterionStrength.SUPPORTING,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"{evidence_id} calibrated prediction is BENIGN (score={ce.score}).",
            evidence_ids=[evidence_id],
        )

    # PATHOGENIC or INDETERMINATE: the benign-supporting condition isn't met either way.
    return CriterionResult(
        code="BP4",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.BENIGN,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=f"{evidence_id} calibrated prediction is {ce.calibrated_prediction.value} (score={ce.score}), not BENIGN.",
        evidence_ids=[evidence_id],
    )
