"""BA1 evaluator — "allele frequency too high for any plausible rare
Mendelian disease" (Richards et al. 2015). Stand-Alone benign evidence:
on its own, sufficient to call a variant Benign.

Mirrors PM2's structure deliberately, in the opposite direction, including
the founder-enrichment ambiguity check: if the overall frequency is below
the BA1 threshold but an ancestry-specific frequency clears it, that's the
same kind of ambiguity PM2 flags — a variant genuinely common only in one
tested population isn't obviously "too common for any disease" the way a
uniformly common variant is. Given the 5% BA1 threshold is so high, this
branch is expected to trigger rarely in practice, but it exists for
consistency with PM2/BS1 rather than treating this evaluator as a special
case with weaker standards.
"""

from ..errors import SchemaValidationError
from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, CriterionStrength, EvidenceDirection, PopulationRetrievalStatus

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_ba1(bundle: VariantEvidenceBundle, thresholds: dict) -> CriterionResult:
    """thresholds is the dict returned by loader.load_frequency_thresholds()."""
    variant_id = bundle.variant.variant_id
    context = f"evaluate_ba1[{variant_id}]"

    if len(bundle.population_evidence) != 1:
        raise SchemaValidationError(
            f"{context}: expected exactly one population_evidence entry, found "
            f"{len(bundle.population_evidence)} — multi-source aggregation is not implemented"
        )
    evidence = bundle.population_evidence[0]
    evidence_id = f"{evidence.source}:{evidence.source_version}"

    if evidence.retrieval_status in (
        PopulationRetrievalStatus.NOT_ASSESSED,
        PopulationRetrievalStatus.UNAVAILABLE,
        PopulationRetrievalStatus.UNKNOWN,
        PopulationRetrievalStatus.NOT_APPLICABLE,
    ):
        return CriterionResult(
            code="BA1",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Population frequency from {evidence_id} has retrieval_status={evidence.retrieval_status.value}; BA1 cannot be evaluated.",
            evidence_ids=[evidence_id],
        )

    if evidence.retrieval_status == PopulationRetrievalStatus.ABSENT:
        return CriterionResult(
            code="BA1",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Absent from {evidence_id} — cannot be 'too common' if unobserved.",
            evidence_ids=[evidence_id],
        )

    # OBSERVED
    ba1_threshold = thresholds["ba1_stand_alone_af"]
    overall_af = evidence.overall_af
    ancestry_af = evidence.ancestry_specific_max_af

    if overall_af >= ba1_threshold:
        return CriterionResult(
            code="BA1",
            status=CriterionStatus.MET,
            strength=CriterionStrength.STAND_ALONE,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Overall AF {overall_af:.6f} in {evidence_id} is at or above the BA1 threshold {ba1_threshold:.6f}.",
            evidence_ids=[evidence_id],
        )

    if ancestry_af is not None and ancestry_af >= ba1_threshold:
        return CriterionResult(
            code="BA1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Overall AF {overall_af:.6f} in {evidence_id} is below the BA1 threshold "
                f"({ba1_threshold:.6f}), but the ancestry-specific maximum AF {ancestry_af:.6f} is not — "
                "a founder-enrichment pattern. Flagged for manual review rather than auto-decided, "
                "consistent with how PM2/BS1 handle the same pattern."
            ),
            evidence_ids=[evidence_id],
        )

    return CriterionResult(
        code="BA1",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.BENIGN,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=f"Overall AF {overall_af:.6f} in {evidence_id} is below the BA1 threshold {ba1_threshold:.6f}.",
        evidence_ids=[evidence_id],
    )
