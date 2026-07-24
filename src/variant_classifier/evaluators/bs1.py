"""BS1 evaluator — "allele frequency greater than expected for this
disorder" (Richards et al. 2015). Strong benign evidence.

Structurally the mirror image of PM2: same gene-specific threshold
pattern (bs1_min_af in config/population_thresholds.yaml), same
founder-enrichment ambiguity check. PM2 asks "rare enough to be
plausible"; BS1 asks "common enough to rule out" — treated as separate
config values per gene (not one shared cutoff enforced in code) since a
real VCEP might set them differently, even though they currently coincide
for CAPN3.
"""

from ..errors import SchemaValidationError
from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, CriterionStrength, EvidenceDirection, PopulationRetrievalStatus

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_bs1(bundle: VariantEvidenceBundle, thresholds: dict) -> CriterionResult:
    """thresholds is the dict returned by loader.load_frequency_thresholds()."""
    variant_id = bundle.variant.variant_id
    gene = bundle.variant.gene
    context = f"evaluate_bs1[{variant_id}]"

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
            code="BS1",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Population frequency from {evidence_id} has retrieval_status={evidence.retrieval_status.value}; BS1 cannot be evaluated.",
            evidence_ids=[evidence_id],
        )

    if evidence.retrieval_status == PopulationRetrievalStatus.ABSENT:
        return CriterionResult(
            code="BS1",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Absent from {evidence_id} — cannot be 'more common than expected' if unobserved.",
            evidence_ids=[evidence_id],
        )

    # OBSERVED
    gene_thresholds = thresholds.get("genes", {})
    if gene not in gene_thresholds:
        raise SchemaValidationError(
            f"{context}: no BS1 frequency threshold configured for gene {gene!r} in config/population_thresholds.yaml"
        )
    bs1_threshold = gene_thresholds[gene]["bs1_min_af"]
    overall_af = evidence.overall_af
    ancestry_af = evidence.ancestry_specific_max_af

    if overall_af >= bs1_threshold:
        return CriterionResult(
            code="BS1",
            status=CriterionStatus.MET,
            strength=CriterionStrength.STRONG,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Overall AF {overall_af:.6f} in {evidence_id} is at or above the {gene} BS1 threshold "
                f"{bs1_threshold:.6f} — more common than expected for this disorder."
            ),
            evidence_ids=[evidence_id],
        )

    if ancestry_af is not None and ancestry_af >= bs1_threshold:
        return CriterionResult(
            code="BS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.BENIGN,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Overall AF {overall_af:.6f} in {evidence_id} is below the {gene} BS1 threshold "
                f"({bs1_threshold:.6f}), but the ancestry-specific maximum AF {ancestry_af:.6f} is not — "
                "a founder-enrichment pattern. A real founder pathogenic allele can be locally common "
                "without being benign, so this is not auto-resolved either way — flagged for manual "
                "review, consistent with how PM2 handles the same pattern on this exact variant."
            ),
            evidence_ids=[evidence_id],
        )

    return CriterionResult(
        code="BS1",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.BENIGN,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=f"Overall AF {overall_af:.6f} in {evidence_id} is below the {gene} BS1 threshold {bs1_threshold:.6f}.",
        evidence_ids=[evidence_id],
    )
