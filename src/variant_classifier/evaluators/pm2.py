"""PM2 evaluator — "absent from controls, or at extremely low frequency if
recessive" (Richards et al. 2015).

Design, worked out against the CAPN3_c.550del fixture specifically:
1. If the population source was never successfully queried (NOT_ASSESSED,
   UNAVAILABLE, UNKNOWN), we cannot say anything about frequency — return
   NOT_EVALUATED rather than guessing.
2. If the source doesn't apply to this variant (NOT_APPLICABLE), likewise
   NOT_EVALUATED — there's nothing to evaluate PM2 against.
3. If the variant was queried and genuinely ABSENT from a well-covered
   locus, that is exactly what PM2 asks for — MET, at the gene's
   configured PM2 strength.
4. If the variant was OBSERVED with a frequency:
   a. Compare the overall allele frequency against the gene's configured
      max-credible-AF threshold. At or above it, PM2 is NOT_MET — too
      common to be this disorder's causal allele.
   b. Below it — but if an ancestry-specific frequency is also reported
      and it clears the threshold while the overall frequency doesn't
      (a founder-enrichment pattern, exactly what CAPN3_c.550del shows:
      0.023% overall vs 0.75% in one ancestry group), do NOT silently
      auto-pass PM2. Flag MANUAL_REVIEW: whether "extremely low frequency"
      holds depends on the tested individual's ancestry, which this
      evaluator has no way to know. This is what keeps the automated
      engine from overclaiming certainty a founder mutation doesn't
      support (see ACMG Engine Detailed Design Guide, Section 7).
   c. Below the threshold with no such ancestry-specific enrichment — MET.

Strength: the evaluator no longer hardcodes Moderate. PM2's strength is
read per-gene from config/population_thresholds.yaml ("pm2_strength"),
defaulting to MODERATE (the generic ACMG/AMP convention) when a gene has
no override. CAPN3 is configured at SUPPORTING, per the real ClinGen LGMD
VCEP specification (see that config file's threshold_source note) — the
VCEP only defines a Supporting-strength PM2 for CAPN3, not Moderate.

This intentionally does not implement per-inheritance-pattern max-credible
allele frequency math (Whiffin et al. 2017) beyond what a VCEP spec
directly supplies; config/population_thresholds.yaml documents remaining
simplifications (e.g. no confidence-interval computation).

Batch 31 (BRCA1) addition — indel/delins exclusion: BRCA1's real ENIGMA
spec does not apply PM2 to insertion, deletion, or delins variants at all
(unlike CAPN3/DMD's generic ACMG PM2, which carries no such exclusion).
This evaluator has no consequence-type awareness by default, and its two
BRCA1 Ashkenazi-founder fixtures (185delAG, 5382insC) need real, OBSERVED
population evidence to drive BA1/BS1's founder-frequency handling — the
same evidence PM2 would otherwise see and evaluate. The original batch 31
design assumed this could be avoided by fixture curation alone, the same
way PM4/PM5 are avoided — it could not: population_evidence is a
required, shared field, not something PM2-specific a curator can omit.
Rather than leave PM2 evaluating a real frequency comparison against a
variant type the real spec excludes entirely, this adds a small,
additive, opt-in gate: a gene may set `pm2_excludes_indel_delins: true`
in population_thresholds.yaml, checked before any retrieval-status
branch. Absent or false for every gene except BRCA1 — CAPN3/DMD are
completely unaffected, since their real PM2 has no such exclusion. See
PM2_INDEL_DELINS_CONSEQUENCES below for exactly which Consequence values
this covers, and tests/unit/test_loader.py's BRCA1 fixture guard test for
the regression check that this gate stays configured and firing.
"""
from typing import Dict

from ..errors import SchemaValidationError
from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import (
    Consequence,
    CriterionStatus,
    CriterionStrength,
    EvidenceDirection,
    PopulationRetrievalStatus,
)

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"

# Batch 31 -- the Consequence values this project treats as "indel/delins"
# for the purpose of a gene's pm2_excludes_indel_delins gate. Frameshift
# variants are indels that shift the reading frame; in-frame deletions/
# insertions are indels that don't. Deliberately narrower than
# TranscriptConsequence.PM4_RELEVANT_CONSEQUENCES (which also includes
# STOP_LOST) -- a stop-loss is not itself an indel/delins variant class,
# it's a consequence that can arise from either an SNV or an indel, and
# BRCA1's real PM2 exclusion is about variant type, not this consequence.
PM2_INDEL_DELINS_CONSEQUENCES = (
    Consequence.FRAMESHIFT_VARIANT,
    Consequence.INFRAME_DELETION,
    Consequence.INFRAME_INSERTION,
)


def evaluate_pm2(bundle: VariantEvidenceBundle, thresholds: dict) -> CriterionResult:
    """thresholds is the dict returned by loader.load_frequency_thresholds():
    {"ba1_stand_alone_af": float, "genes": {gene: {"pm2_max_credible_af": ..., "bs1_min_af": ..., ...}}}.
    """
    variant_id = bundle.variant.variant_id
    gene = bundle.variant.gene
    context = f"evaluate_pm2[{variant_id}]"
    if len(bundle.population_evidence) != 1:
        raise SchemaValidationError(
            f"{context}: expected exactly one population_evidence entry, found "
            f"{len(bundle.population_evidence)} — multi-source aggregation is not implemented in Milestone 2"
        )
    evidence = bundle.population_evidence[0]
    evidence_id = f"{evidence.source}:{evidence.source_version}"

    gene_thresholds = thresholds.get("genes", {})
    gene_config = gene_thresholds.get(gene, {})

    if gene_config.get("pm2_excludes_indel_delins", False):
        # Safe without a fallback: VariantEvidenceBundle.__post_init__ already
        # guarantees exactly one clinically_relevant=True entry exists before
        # this bundle can be constructed at all (the same invariant
        # evaluators/pm4.py's identical next(...) call relies on).
        transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)
        if transcript.consequence in PM2_INDEL_DELINS_CONSEQUENCES:
            return CriterionResult(
                code="PM2",
                status=CriterionStatus.NOT_APPLICABLE,
                direction=EvidenceDirection.PATHOGENIC,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"{gene}'s real specification does not apply PM2 to indel/delins variants "
                    f"(pm2_excludes_indel_delins=True in population_thresholds.yaml); "
                    f"{transcript.consequence.value} in {transcript.transcript_id} is such a "
                    "variant, so PM2 does not apply regardless of population frequency data."
                ),
                evidence_ids=[evidence_id],
            )

    if evidence.retrieval_status in (
        PopulationRetrievalStatus.NOT_ASSESSED,
        PopulationRetrievalStatus.UNAVAILABLE,
        PopulationRetrievalStatus.UNKNOWN,
    ):
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Population frequency from {evidence_id} has retrieval_status="
                f"{evidence.retrieval_status.value}; PM2 cannot be evaluated without a successful lookup."
            ),
            evidence_ids=[evidence_id],
        )
    if evidence.retrieval_status == PopulationRetrievalStatus.NOT_APPLICABLE:
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"{evidence_id} does not apply to this variant; PM2 has no frequency evidence to evaluate.",
            evidence_ids=[evidence_id],
        )
    pm2_strength = CriterionStrength[gene_config.get("pm2_strength", "MODERATE")]
    if evidence.retrieval_status == PopulationRetrievalStatus.ABSENT:
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.MET,
            strength=pm2_strength,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Absent from {evidence_id} at a locus with adequate coverage — "
                "meets PM2's 'absent from controls' condition."
            ),
            evidence_ids=[evidence_id],
        )
    # retrieval_status == OBSERVED from here on (the only remaining enum value).
    # Unlike strength (which defaults to MODERATE if unconfigured), the numeric
    # max-credible-AF threshold has no safe default — a gene must be explicitly
    # configured before an OBSERVED frequency can be compared against anything.
    if gene not in gene_thresholds:
        raise SchemaValidationError(
            f"{context}: no PM2 frequency threshold configured for gene {gene!r} "
            "in config/population_thresholds.yaml"
        )
    max_credible_af = gene_thresholds[gene]["pm2_max_credible_af"]
    overall_af = evidence.overall_af
    ancestry_af = evidence.ancestry_specific_max_af
    if overall_af >= max_credible_af:
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Overall AF {overall_af:.6f} in {evidence_id} is at or above the configured "
                f"max-credible-AF threshold {max_credible_af:.6f} for {gene} — too common for PM2."
            ),
            evidence_ids=[evidence_id],
        )
    if ancestry_af is not None and ancestry_af >= max_credible_af:
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Overall AF {overall_af:.6f} in {evidence_id} is below the {gene} threshold "
                f"({max_credible_af:.6f}), but the ancestry-specific maximum AF {ancestry_af:.6f} "
                "is not — a founder-enrichment pattern. Whether 'extremely low frequency' holds "
                "depends on the tested individual's ancestry, which is not available here; "
                "flagged for manual review rather than auto-decided."
            ),
            evidence_ids=[evidence_id],
        )
    return CriterionResult(
        code="PM2",
        status=CriterionStatus.MET,
        strength=pm2_strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Overall AF {overall_af:.6f} in {evidence_id} is below the configured max-credible-AF "
            f"threshold {max_credible_af:.6f} for {gene}, with no ancestry-specific enrichment on record."
        ),
        evidence_ids=[evidence_id],
    )
