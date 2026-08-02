"""PS3 evaluator -- "well-established in vitro or in vivo functional
studies supportive of a damaging effect on the gene or gene product"
(Richards et al. 2015, Table 3). Strong pathogenic evidence in the base
ACMG/AMP framework, added batch 25 alongside BS3 (evaluators/bs3.py).

See models/functional_evidence.py for the full design writeup: why
"well-established" is treated as a curated validation_strength tier
(Brnich et al. 2019's Supporting/Moderate/Strong ladder) rather than a
fixed Strong-or-nothing call, why INDETERMINATE is a real third assay
outcome distinct from "no functional evidence at all," and why this
project represents one already-reconciled FunctionalEvidence record per
variant rather than re-deriving a verdict from multiple raw assays.

Scope, disclosed rather than silently assumed: this evaluator applies
the base Richards et al. 2015 / Brnich et al. 2019 framework only. It
does not implement a VCEP-specific PS3 refinement (e.g. the CAPN3 LGMD
VCEP's own functional-assay specification, if one exists beyond the base
framework) -- same "deliberately partial, base-framework-only" treatment
PS1/PM5/PVS1 have had since earlier milestones.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, EvidenceDirection, FunctionalAssayResult

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015); strength per Brnich et al. 2019 (Genome Medicine 11:98)"
RULE_VERSION = "2015 / 2019"


def evaluate_ps3(bundle: VariantEvidenceBundle) -> CriterionResult:
    fe = bundle.functional_evidence
    evidence_id = f"functional_evidence:{bundle.variant.variant_id}"

    if fe is None:
        return CriterionResult(
            code="PS3",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="No functional_evidence recorded for this variant -- PS3 was never assessed.",
            evidence_ids=[evidence_id],
        )

    if fe.assay_result == FunctionalAssayResult.INDETERMINATE:
        return CriterionResult(
            code="PS3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                "Functional assay was performed but did not clearly discriminate a damaging effect "
                f"({fe.notes or 'see curated notes'}) -- an indeterminate result supports neither PS3 "
                "nor BS3."
            ),
            evidence_ids=[evidence_id],
        )

    if fe.assay_result == FunctionalAssayResult.NORMAL:
        return CriterionResult(
            code="PS3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="Functional assay showed a normal (non-damaging) result -- this is BS3 evidence, not PS3.",
            evidence_ids=[evidence_id],
        )

    # ABNORMAL
    return CriterionResult(
        code="PS3",
        status=CriterionStatus.MET,
        strength=fe.validation_strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Functional assay ({fe.assay_description or 'assay not further described'}) showed a "
            f"damaging effect, at curated validation strength {fe.validation_strength.value} "
            "(Brnich et al. 2019 validation-control tier)."
        ),
        evidence_ids=[evidence_id],
    )
