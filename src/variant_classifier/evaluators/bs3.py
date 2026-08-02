"""BS3 evaluator -- "well-established in vitro or in vivo functional
studies show no damaging effect on protein function or splicing"
(Richards et al. 2015, Table 3). Strong benign evidence in the base
ACMG/AMP framework, added batch 25 alongside PS3 (evaluators/ps3.py) --
see that file and models/functional_evidence.py for the shared design
writeup. Structurally the mirror image of evaluate_ps3(): same three-way
assay_result branch, opposite MET direction.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, EvidenceDirection, FunctionalAssayResult

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015); strength per Brnich et al. 2019 (Genome Medicine 11:98)"
RULE_VERSION = "2015 / 2019"


def evaluate_bs3(bundle: VariantEvidenceBundle) -> CriterionResult:
    fe = bundle.functional_evidence
    evidence_id = f"functional_evidence:{bundle.variant.variant_id}"

    if fe is None:
        return CriterionResult(
            code="BS3",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="No functional_evidence recorded for this variant -- BS3 was never assessed.",
            evidence_ids=[evidence_id],
        )

    if fe.assay_result == FunctionalAssayResult.INDETERMINATE:
        return CriterionResult(
            code="BS3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                "Functional assay was performed but did not clearly discriminate a non-damaging effect "
                f"({fe.notes or 'see curated notes'}) -- an indeterminate result supports neither BS3 "
                "nor PS3."
            ),
            evidence_ids=[evidence_id],
        )

    if fe.assay_result == FunctionalAssayResult.ABNORMAL:
        return CriterionResult(
            code="BS3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="Functional assay showed an abnormal (damaging) result -- this is PS3 evidence, not BS3.",
            evidence_ids=[evidence_id],
        )

    # NORMAL
    return CriterionResult(
        code="BS3",
        status=CriterionStatus.MET,
        strength=fe.validation_strength,
        direction=EvidenceDirection.BENIGN,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Functional assay ({fe.assay_description or 'assay not further described'}) showed no "
            f"damaging effect, at curated validation strength {fe.validation_strength.value} "
            "(Brnich et al. 2019 validation-control tier)."
        ),
        evidence_ids=[evidence_id],
    )
