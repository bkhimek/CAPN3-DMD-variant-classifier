"""PM2 evaluator — "absent from controls, or at extremely low frequency if
recessive" (Richards et al. 2015).

Design, worked out against the CAPN3_c.550del fixture specifically:

1. If the population source was never successfully queried (NOT_ASSESSED,
   UNAVAILABLE, UNKNOWN), we cannot say anything about frequency — return
   NOT_EVALUATED rather than guessing.
2. If the source doesn't apply to this variant (NOT_APPLICABLE), likewise
   NOT_EVALUATED — there's nothing to evaluate PM2 against.
3. If the variant was queried and genuinely ABSENT from a well-covered
   locus, that is exactly what PM2 asks for — MET, Moderate strength.
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

This intentionally does not implement per-inheritance-pattern max-credible
allele frequency math (Whiffin et al. 2017); config/population_thresholds.yaml
documents that as a known simplification.
"""

from typing import Dict

from ..errors import SchemaValidationError
from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, CriterionStrength, EvidenceDirection, PopulationRetrievalStatus

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_pm2(bundle: VariantEvidenceBundle, thresholds: Dict[str, dict]) -> CriterionResult:
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

    if evidence.retrieval_status == PopulationRetrievalStatus.ABSENT:
        return CriterionResult(
            code="PM2",
            status=CriterionStatus.MET,
            strength=CriterionStrength.MODERATE,
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
    if gene not in thresholds:
        raise SchemaValidationError(
            f"{context}: no PM2 frequency threshold configured for gene {gene!r} "
            "in config/population_thresholds.yaml"
        )
    max_credible_af = thresholds[gene]["pm2_max_credible_af"]
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
        strength=CriterionStrength.MODERATE,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Overall AF {overall_af:.6f} in {evidence_id} is below the configured max-credible-AF "
            f"threshold {max_credible_af:.6f} for {gene}, with no ancestry-specific enrichment on record."
        ),
        evidence_ids=[evidence_id],
    )
