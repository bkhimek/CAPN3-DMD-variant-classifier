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
3. Splice donor/acceptor variants (batch 26 update — see "PS3 and BS3"-
   adjacent design note "Splice-RNA evidence feeds PVS1 directly (batch
   26)" in README.md for the full writeup): resolving these in general
   needs a predicted splicing outcome (in-frame vs out-of-frame exon
   skip, intron retention) this project does not compute. But when a real
   RNA/splicing assay result is curated (`TranscriptConsequence.splicing_rna_evidence`,
   see models/enums.py's `SplicingRnaEvidence`), it is used directly,
   mirroring the real ClinGen LGMD VCEP CAPN3 specification's own
   instruction that experimental splicing evidence should be scored
   under PVS1 (not PS3), per the ClinGen SVI Splicing Subgroup's decision
   tree (Walker et al. 2023, PMID 37352859):
   - CONFIRMED_NULL_EQUIVALENT: the assay confirms an aberrant transcript
     functionally equivalent to a null allele — treated identically to a
     confirmed-NMD frameshift/nonsense variant — MET, Very Strong.
   - CONFIRMED_IN_FRAME_OR_PARTIAL_FUNCTION: aberrant splicing confirmed,
     but in-frame — same open protein-domain-criticality gap as the
     NMD-escape branch — MANUAL_REVIEW.
   - CONFIRMED_NORMAL_SPLICING: the assay directly contradicts the
     predicted splice disruption — the null-variant mechanism this
     evaluator would otherwise assume is refuted — NOT_MET, not a guess
     in either remaining direction.
   - INCONCLUSIVE: an assay was performed but didn't clearly establish
     any of the above — MANUAL_REVIEW, same "checked, still uncertain"
     treatment as everything else.
   - Not set (the default/common case — no RNA assay at all): falls
     through to the original, unchanged predicted-only MANUAL_REVIEW
     path.
   Only the CONFIRMED_NULL_EQUIVALENT/CONFIRMED_NORMAL_SPLICING branches
   are exact-threshold-free (a real assay directly establishes the
   outcome, no percentage cutoffs needed); everything requiring the
   original decision tree's numeric thresholds (percentage of transcript
   affected, protein-region criticality) remains unimplemented, since the
   primary Walker et al. 2023 paper and the CAPN3-specific PVS1 flowchart
   PDF were both unreachable during this batch's research (see README).
4. Start-loss variants remain within PVS1's scope in principle, but
   resolving them needs an alternative-start-codon check this project
   does not implement — MANUAL_REVIEW, explaining why (unchanged from
   before batch 26; not addressed this round).
5. Every other consequence type isn't a null-variant class PVS1 applies
   to — NOT_APPLICABLE.

Deliberately conservative: this evaluator only ever returns MET for the
cases it can defend end-to-end (early truncation with NMD predicted or
confirmed via RNA evidence, in a gene with an established loss-of-function
mechanism). Everything harder is MANUAL_REVIEW, never a guessed MET or
NOT_MET — except CONFIRMED_NORMAL_SPLICING, where "not met" is itself the
defensible, non-guessed conclusion a direct experimental contradiction
supports.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import (
    Consequence,
    CriterionStatus,
    CriterionStrength,
    DiseaseMechanism,
    EvidenceDirection,
    SplicingRnaEvidence,
)

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015); splice-RNA-evidence branch per the ClinGen LGMD VCEP CAPN3 specification and Walker et al. 2023 (PMID 37352859)"
RULE_VERSION = "2015 / 2023"

_NMD_RELEVANT_CONSEQUENCES = (Consequence.FRAMESHIFT_VARIANT, Consequence.STOP_GAINED)
_SPLICE_CONSEQUENCES = (Consequence.SPLICE_DONOR_VARIANT, Consequence.SPLICE_ACCEPTOR_VARIANT)
_SPLICE_OR_START_CONSEQUENCES = _SPLICE_CONSEQUENCES + (Consequence.START_LOST,)


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

    if consequence in _SPLICE_CONSEQUENCES:
        rna_evidence = transcript.splicing_rna_evidence

        if rna_evidence == SplicingRnaEvidence.CONFIRMED_NULL_EQUIVALENT:
            return CriterionResult(
                code="PVS1",
                status=CriterionStatus.MET,
                strength=CriterionStrength.VERY_STRONG,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{consequence.value} in {transcript.transcript_id}, with a real RNA/splicing "
                    "assay confirming the resulting transcript is functionally equivalent to a null "
                    "allele (out-of-frame exon skip, intron retention producing a frameshift, or a "
                    "confirmed premature termination codon) — treated identically to a confirmed-NMD "
                    "frameshift/nonsense variant, in a gene with an established loss-of-function "
                    "mechanism."
                ),
                evidence_ids=[evidence_id],
            )

        if rna_evidence == SplicingRnaEvidence.CONFIRMED_NORMAL_SPLICING:
            return CriterionResult(
                code="PVS1",
                status=CriterionStatus.NOT_MET,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{consequence.value} in {transcript.transcript_id} was predicted to disrupt "
                    "splicing, but a real RNA/splicing assay directly confirms splicing is NOT "
                    "disrupted — the null-variant mechanism this consequence class assumes is "
                    "refuted by direct evidence, so PVS1 does not apply on that basis (the variant "
                    "may still act through some other mechanism this evaluator does not check)."
                ),
                evidence_ids=[evidence_id],
            )

        if rna_evidence == SplicingRnaEvidence.CONFIRMED_IN_FRAME_OR_PARTIAL_FUNCTION:
            return CriterionResult(
                code="PVS1",
                status=CriterionStatus.MANUAL_REVIEW,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{consequence.value} in {transcript.transcript_id}: a real RNA/splicing assay "
                    "confirms aberrant splicing, but the resulting transcript remains in-frame. "
                    "Whether that altered protein still functions depends on protein-domain "
                    "criticality information this evaluator does not model — flagged for manual "
                    "review rather than guessed, the same open gap as the NMD-escape branch."
                ),
                evidence_ids=[evidence_id],
            )

        if rna_evidence == SplicingRnaEvidence.INCONCLUSIVE:
            return CriterionResult(
                code="PVS1",
                status=CriterionStatus.MANUAL_REVIEW,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{consequence.value} in {transcript.transcript_id}: a real RNA/splicing assay "
                    "was performed but did not clearly establish the resulting transcript's frame or "
                    "NMD status — flagged for manual review rather than guessed either way."
                ),
                evidence_ids=[evidence_id],
            )

        # No RNA/splicing assay curated at all — the original, unchanged predicted-only path.
        return CriterionResult(
            code="PVS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} in {transcript.transcript_id} falls within PVS1's scope in "
                "principle, but resolving it requires a predicted splicing outcome this evaluator "
                "does not compute, and no real RNA/splicing assay evidence is curated for this "
                "variant (transcript_consequences.splicing_rna_evidence is unset) — flagged for "
                "manual review rather than guessed."
            ),
            evidence_ids=[evidence_id],
        )

    if consequence == Consequence.START_LOST:
        return CriterionResult(
            code="PVS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{consequence.value} in {transcript.transcript_id} falls within PVS1's scope in "
                "principle, but resolving it requires an alternative-start-codon check that this "
                "evaluator does not yet implement — flagged for manual review rather than guessed."
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
