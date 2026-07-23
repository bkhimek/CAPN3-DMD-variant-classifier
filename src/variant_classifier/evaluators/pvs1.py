"""PVS1 evaluator — "null variant in a gene where loss-of-function is a
known disease mechanism" (Richards et al. 2015).

This is a deliberately partial evaluator. The full PVS1 decision tree
(Abou Tayoun et al. 2018) branches on protein-domain criticality and
constitutive-exon-splicing information that this project does not model
yet (see README.md, "PVS1 scope"). What's implemented:

1. Gate on established disease mechanism. If loss-of-function isn't the
   gene's known mechanism, PVS1 doesn't apply at all — NOT_APPLICABLE.
2. For frameshift and nonsense (stop_gained) variants — the two
   consequence types TranscriptConsequence requires an explicit
   nmd_predicted value for — follow the NMD branch:
   - nmd_predicted=True: the transcript is expected to be degraded before
     translation. This is the clean case PVS1 is built for — MET,
     Very Strong.
   - nmd_predicted=False: no NMD, so a truncated protein is made instead.
     Whether that shortened protein still functions depends on protein
     structure this project doesn't model — MANUAL_REVIEW, not a guess.
3. Splice donor/acceptor and start-loss variants are within PVS1's scope
   in principle, but resolving them needs a predicted splicing outcome or
   an alternative-start-codon check that isn't implemented —
   MANUAL_REVIEW, explaining why.
4. Every other consequence type isn't a null-variant class PVS1 applies
   to — NOT_APPLICABLE.

Deliberately conservative: this evaluator only ever returns MET for the
one case it can defend end-to-end (early truncation, NMD predicted, in a
gene with an established loss-of-function mechanism). Everything harder
is MANUAL_REVIEW, never a guessed MET or NOT_MET.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import Consequence, CriterionStatus, CriterionStrength, DiseaseMechanism, EvidenceDirection

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"

_NMD_RELEVANT_CONSEQUENCES = (Consequence.FRAMESHIFT_VARIANT, Consequence.STOP_GAINED)
_SPLICE_OR_START_CONSEQUENCES = (
    Consequence.SPLICE_DONOR_VARIANT,
    Consequence.SPLICE_ACCEPTOR_VARIANT,
    Consequence.START_LOST,
)


def evaluate_pvs1(bundle: VariantEvidenceBundle) -> CriterionResult:
    ctx = bundle.gene_disease_context
    transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)
    evidence_id = f"transcript:{transcript.transcript_id}"

    if ctx.mechanism != DiseaseMechanism.LOSS_OF_FUNCTION or not ctx.lof_established:
        return CriterionResult(
            code="PVS1",
            status=CriterionStatus.NOT_APPLICABLE,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"PVS1 requires an established loss-of-function disease mechanism for {ctx.gene}; "
                f"gene_disease_context.mechanism={ctx.mechanism.value}, lof_established={ctx.lof_established}."
            ),
            evidence_ids=[evidence_id],
        )

    consequence = transcript.consequence

    if consequence in _NMD_RELEVANT_CONSEQUENCES:
        if transcript.nmd_predicted:
            return CriterionResult(
                code="PVS1",
                status=CriterionStatus.MET,
                strength=CriterionStrength.VERY_STRONG,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{consequence.value} in {transcript.transcript_id}, predicted to trigger "
                    f"nonsense-mediated decay, in a gene with an established loss-of-function mechanism."
                ),
                evidence_ids=[evidence_id],
            )
        return CriterionResult(
            code="PVS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} in {transcript.transcript_id} is not predicted to trigger "
                "nonsense-mediated decay (likely last-exon or near the final exon-exon junction), so a "
                "truncated protein is expected instead of no protein. Whether that truncated protein "
                "still functions depends on protein-domain information this evaluator does not model — "
                "flagged for manual review rather than guessed."
            ),
            evidence_ids=[evidence_id],
        )

    if consequence in _SPLICE_OR_START_CONSEQUENCES:
        return CriterionResult(
            code="PVS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} in {transcript.transcript_id} falls within PVS1's scope in "
                "principle, but resolving it requires a predicted splicing outcome or an "
                "alternative-start-codon check that this evaluator does not yet implement — "
                "flagged for manual review rather than guessed."
            ),
            evidence_ids=[evidence_id],
        )

    return CriterionResult(
        code="PVS1",
        status=CriterionStatus.NOT_APPLICABLE,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=f"{consequence.value} is not a null-variant consequence class that PVS1 applies to.",
        evidence_ids=[evidence_id],
    )
