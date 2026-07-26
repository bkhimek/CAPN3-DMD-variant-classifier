"""PM4 evaluator — "protein length changes as a result of in-frame
deletions/insertions in a nonrepeat region or stop-loss variants"
(Richards et al. 2015, Table 3). Moderate pathogenic evidence.

Added in batch 14, this project's first new criterion since Milestone 3.
It was scoped out of Milestones 1-3 by accident of fixture selection, not
by design: `Consequence` has always included INFRAME_DELETION,
INFRAME_INSERTION, and STOP_LOST, but no curated fixture ever used any of
them, so the gap stayed invisible until batch 14 went looking for a real
stop-loss variant to test PVS1's scope and found — after checking the
real ACMG/AMP criterion definitions rather than assuming — that
stop-loss variants aren't a PVS1 concern at all (PVS1 is null/loss-of-
function-specific; a stop-loss variant produces an elongated protein,
not a null one). PM4 is the correct home for them, and for in-frame
indels, per Richards et al. 2015 Table 3 directly.

Design:

1. Only INFRAME_DELETION, INFRAME_INSERTION, and STOP_LOST are in scope
   — everything else is NOT_APPLICABLE. Unlike PVS1, PM4 is not gated on
   an established loss-of-function mechanism: a protein-length change can
   matter for genes with other mechanisms too, and Richards et al. 2015
   does not condition PM4 on mechanism the way it does PVS1.
2. PM4 explicitly excludes repeat/low-complexity regions ("in a
   nonrepeat region") — an indel in a repetitive region is far less
   likely to be functionally significant, and far more likely to be a
   sequencing/alignment artifact. `TranscriptConsequence.repeat_region`
   must be stated explicitly for every PM4-relevant consequence (see its
   __post_init__) — never guessed. If it's True, PM4 is NOT_MET, no
   matter how the protein length changed.
3. Strength: MODERATE by default. The ClinGen SVI's own caution note on
   PM4 (echoed across several VCEP specifications) is that a single
   amino acid in-frame indel should generally be downgraded to
   SUPPORTING unless gene-specific evidence justifies Moderate — this
   evaluator applies that downgrade whenever
   `protein_length_change_aa` is recorded and equal to 1. When the size
   isn't recorded at all, this evaluator defaults to MODERATE rather
   than refusing to answer — a disclosed simplification (like PM2's
   omitted confidence-interval math), not a silent guess: the rationale
   always states plainly whether the downgrade was actually evaluated.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import Consequence, CriterionStatus, CriterionStrength, EvidenceDirection

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"

_PM4_RELEVANT_CONSEQUENCES = (Consequence.INFRAME_DELETION, Consequence.INFRAME_INSERTION, Consequence.STOP_LOST)


def evaluate_pm4(bundle: VariantEvidenceBundle) -> CriterionResult:
    transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)
    evidence_id = f"transcript:{transcript.transcript_id}"
    consequence = transcript.consequence

    if consequence not in _PM4_RELEVANT_CONSEQUENCES:
        return CriterionResult(
            code="PM4",
            status=CriterionStatus.NOT_APPLICABLE,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} is not a protein-length-changing consequence class "
                "(in-frame deletion/insertion, stop-loss) that PM4 applies to."
            ),
            evidence_ids=[evidence_id],
        )

    if transcript.repeat_region:
        return CriterionResult(
            code="PM4",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} in {transcript.transcript_id} falls within a repeat/"
                "low-complexity region — PM4 explicitly excludes these by definition (Richards "
                "et al. 2015), regardless of the resulting protein-length change."
            ),
            evidence_ids=[evidence_id],
        )

    if transcript.protein_length_change_aa == 1:
        strength = CriterionStrength.SUPPORTING
        size_note = (
            "a single amino acid in-frame change — downgraded to Supporting per the ClinGen SVI's "
            "caution against using full Moderate strength for single-residue indels absent "
            "gene-specific evidence otherwise"
        )
    elif transcript.protein_length_change_aa is not None:
        strength = CriterionStrength.MODERATE
        size_note = f"a {transcript.protein_length_change_aa}-amino-acid protein-length change"
    else:
        strength = CriterionStrength.MODERATE
        size_note = (
            "an unrecorded protein-length change (protein_length_change_aa not curated for this "
            "fixture) — defaulting to the standard Moderate strength rather than guessing whether "
            "a single-residue downgrade applies; a disclosed simplification, not a silent guess"
        )

    return CriterionResult(
        code="PM4",
        status=CriterionStatus.MET,
        strength=strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"{consequence.value} in {transcript.transcript_id}, confirmed outside a repeat "
            f"region, resulting in {size_note}."
        ),
        evidence_ids=[evidence_id],
    )
