"""Tests for the BA1 and BS1 evaluators — PM2's mirror images on the
benign side. Same two-part structure: golden-case cross-check against the
real fixtures, then hand-built edge cases.
"""

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators import evaluate_ba1, evaluate_bs1
from variant_classifier.models import (
    GeneDiseaseContext,
    PopulationEvidence,
    Specification,
    TranscriptConsequence,
    VariantEvidenceBundle,
    VariantIdentity,
)
from variant_classifier.models.enums import (
    Consequence,
    CriterionStatus,
    CriterionStrength,
    DiseaseMechanism,
    GenomeBuild,
    Inheritance,
    PopulationRetrievalStatus,
    SpecificationType,
)


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


# ------------------------------------------------------- against real fixtures

def test_ba1_and_bs1_match_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()
    thresholds = loader.load_frequency_thresholds()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        for code, evaluate in (("BA1", evaluate_ba1), ("BS1", evaluate_bs1)):
            if code not in golden.expected_criterion_status:
                continue
            result = evaluate(bundle, thresholds)
            expected = golden.expected_criterion_status[code]
            assert result.status == expected, (
                f"{bundle.variant.variant_id}: {code} evaluated to {result.status}, "
                f"golden case expects {expected}. Rationale: {result.rationale}"
            )
            checked += 1
    # Two codes (BA1, BS1) checked per bundle -- see test_pm2_evaluator.py
    # for why this is derived from len(bundles) rather than a literal.
    assert checked == 2 * len(bundles) > 0


def test_bs1_founder_case_is_flagged_manual_review_not_silently_not_met():
    # This is the case that motivated correcting the golden case from
    # NOT_MET to MANUAL_REVIEW — see validation/golden_cases/variant_golden_cases.yaml.
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = evaluate_bs1(bundle, thresholds)
    assert result.status == CriterionStatus.MANUAL_REVIEW
    assert "founder" in result.rationale.lower()


# ------------------------------------------------------- hand-built edge cases

def _bundle_with_population_evidence(pe: PopulationEvidence) -> VariantEvidenceBundle:
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene="CAPN3", disease="LGMDR1", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True,
            specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
        ),
        transcript_consequences=[
            TranscriptConsequence(transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT)
        ],
        population_evidence=[pe],
    )


def _thresholds():
    return {
        "ba1_stand_alone_af": 0.05,
        "genes": {"CAPN3": {"pm2_max_credible_af": 0.001, "bs1_min_af": 0.001, "threshold_source": "test"}},
    }


def test_ba1_absent_yields_not_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
            locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
        )
    )
    result = evaluate_ba1(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_MET


def test_ba1_above_threshold_yields_met_stand_alone():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.08)
    )
    result = evaluate_ba1(bundle, _thresholds())
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STAND_ALONE


def test_ba1_below_threshold_yields_not_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.001)
    )
    result = evaluate_ba1(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_MET


def test_ba1_not_assessed_yields_not_evaluated():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED)
    )
    result = evaluate_ba1(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_ba1_gene_specific_override_is_honored():
    # Covers the batch-4 addition: BA1's threshold can be overridden per
    # gene via genes.<GENE>.ba1_af (CAPN3's real ClinGen LGMD VCEP value is
    # 0.003, far stricter than the generic 0.05 default). An AF that would
    # be NOT_MET under the generic default should be MET under the
    # gene-specific override.
    thresholds = {
        "ba1_stand_alone_af": 0.05,
        "genes": {"CAPN3": {"pm2_max_credible_af": 0.0001, "bs1_min_af": 0.001, "ba1_af": 0.003, "threshold_source": "test"}},
    }
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.008)
    )
    result = evaluate_ba1(bundle, thresholds)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STAND_ALONE


def test_ba1_falls_back_to_global_default_when_gene_unconfigured():
    # A gene with no ba1_af override (e.g. DMD, no VCEP spec adopted here)
    # should keep using the global generic-ACMG default (0.05), not error
    # or silently use some other value.
    thresholds = {"ba1_stand_alone_af": 0.05, "genes": {}}
    bundle_pe = PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.008)
    bundle = VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE_DMD", gene="DMD", genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene="DMD", disease="dystrophinopathy", inheritance=Inheritance.X_LINKED_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True,
            specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
        ),
        transcript_consequences=[
            TranscriptConsequence(transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT)
        ],
        population_evidence=[bundle_pe],
    )
    result = evaluate_ba1(bundle, thresholds)
    # 0.008 is below the 5% global default, so NOT_MET -- if the code
    # incorrectly fell back to CAPN3's 0.3% override or errored, this would
    # catch it.
    assert result.status == CriterionStatus.NOT_MET


def test_bs1_absent_yields_not_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
            locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
        )
    )
    result = evaluate_bs1(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_MET


def test_bs1_above_threshold_yields_met_strong():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.01)
    )
    result = evaluate_bs1(bundle, _thresholds())
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STRONG


def test_bs1_below_threshold_yields_not_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.00001)
    )
    result = evaluate_bs1(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_MET


def test_bs1_missing_gene_threshold_raises():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED, overall_af=0.01)
    )
    expect_schema_error(lambda: evaluate_bs1(bundle, {"ba1_stand_alone_af": 0.05, "genes": {}}))
