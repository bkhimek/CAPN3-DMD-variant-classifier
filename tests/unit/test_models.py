"""Schema-validation tests for the seven typed models.

Written in pytest convention (test_*.py filename, test_ function names,
plain `assert`, `pytest.raises`-free style using try/except instead) so
these run unmodified with real pytest once your environment has it
installed, or right now via tests/run_tests.py, which does not depend on
pytest being present. See README.md for why pytest could not be installed
during Milestone 1 development.
"""

from variant_classifier.errors import SchemaValidationError
from variant_classifier.models import (
    ComputationalEvidence,
    CriterionResult,
    GeneDiseaseContext,
    GoldenCase,
    PopulationEvidence,
    ProvisionalClassification,
    Specification,
    TranscriptConsequence,
    VariantEvidenceBundle,
    VariantIdentity,
)
from variant_classifier.models.enums import (
    ClassificationStatus,
    ComputationalPrediction,
    CriterionStatus,
    CriterionStrength,
    DiseaseMechanism,
    EvidenceDirection,
    GenomeBuild,
    Inheritance,
    PopulationRetrievalStatus,
    ProvisionalClass,
    SpecificationType,
)


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


# ---------------------------------------------------------------- VariantIdentity

def test_variant_identity_valid_construction():
    v = VariantIdentity(variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38)
    assert v.variant_id == "v1"
    assert v.coordinate_verified is True


def test_variant_identity_unverified_coordinates_cannot_carry_position():
    expect_schema_error(
        lambda: VariantIdentity(
            variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38,
            position=100, coordinate_verified=False,
        )
    )


def test_variant_identity_unverified_coordinates_ok_without_position():
    v = VariantIdentity(variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38, coordinate_verified=False)
    assert v.coordinate_verified is False


def test_variant_identity_from_dict_rejects_bad_genome_build():
    expect_schema_error(
        lambda: VariantIdentity.from_dict({"variant_id": "v1", "gene": "CAPN3", "genome_build": "hg19"})
    )


# ---------------------------------------------------------------- GeneDiseaseContext

def _valid_spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="Richards et al. 2015")


def test_gene_disease_context_valid_construction():
    ctx = GeneDiseaseContext(
        gene="CAPN3", disease="LGMDR1", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
        mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_valid_spec(),
    )
    assert ctx.gene == "CAPN3"


def test_specification_rejects_empty_version():
    expect_schema_error(lambda: Specification(type=SpecificationType.GENERIC_ACMG, version="   "))


# ---------------------------------------------------------------- TranscriptConsequence

def test_transcript_consequence_frameshift_requires_nmd_predicted():
    from variant_classifier.models.enums import Consequence
    expect_schema_error(
        lambda: TranscriptConsequence(
            transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.FRAMESHIFT_VARIANT,
        )
    )


def test_transcript_consequence_frameshift_with_nmd_predicted_ok():
    from variant_classifier.models.enums import Consequence
    tc = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.FRAMESHIFT_VARIANT,
        nmd_predicted=True,
    )
    assert tc.nmd_predicted is True


# ---------------------------------------------------------------- PopulationEvidence

def test_population_evidence_observed_requires_overall_af():
    expect_schema_error(
        lambda: PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        )
    )


def test_population_evidence_observed_with_af_ok():
    pe = PopulationEvidence(
        source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        overall_af=0.001,
    )
    assert pe.overall_af == 0.001


def test_population_evidence_absent_requires_locus_coverage_adequate():
    expect_schema_error(
        lambda: PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
        )
    )


def test_population_evidence_unavailable_forbids_frequency_fields():
    expect_schema_error(
        lambda: PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.UNAVAILABLE,
            overall_af=0.001,
        )
    )


# ---------------------------------------------------------------- ComputationalEvidence

def test_computational_evidence_absent_status_rejected():
    expect_schema_error(
        lambda: ComputationalEvidence(
            tool="REVEL", tool_version="1.0", calibration_source="test",
            retrieval_status=PopulationRetrievalStatus.ABSENT,
        )
    )


def test_computational_evidence_observed_requires_calibrated_prediction():
    expect_schema_error(
        lambda: ComputationalEvidence(
            tool="REVEL", tool_version="1.0", calibration_source="test",
            retrieval_status=PopulationRetrievalStatus.OBSERVED,
        )
    )


def test_computational_evidence_observed_with_prediction_ok():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test",
        retrieval_status=PopulationRetrievalStatus.OBSERVED,
        score=0.9, calibrated_prediction=ComputationalPrediction.PATHOGENIC,
    )
    assert ce.calibrated_prediction == ComputationalPrediction.PATHOGENIC


# ---------------------------------------------------------------- CriterionResult

def test_criterion_result_rejects_unknown_code():
    expect_schema_error(
        lambda: CriterionResult(
            code="PZ99", status=CriterionStatus.NOT_MET, direction=EvidenceDirection.PATHOGENIC,
            rule_source="ACMG", rule_version="2015", rationale="n/a",
        )
    )


def test_criterion_result_met_requires_strength():
    expect_schema_error(
        lambda: CriterionResult(
            code="PM2", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC,
            rule_source="ACMG", rule_version="2015", rationale="absent from gnomAD",
        )
    )


def test_criterion_result_not_met_forbids_strength():
    expect_schema_error(
        lambda: CriterionResult(
            code="PM2", status=CriterionStatus.NOT_MET, direction=EvidenceDirection.PATHOGENIC,
            rule_source="ACMG", rule_version="2015", rationale="common variant",
            strength=CriterionStrength.MODERATE,
        )
    )


def test_criterion_result_met_with_strength_ok():
    cr = CriterionResult(
        code="PM2", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC,
        rule_source="ACMG", rule_version="2015", rationale="absent from gnomAD",
        strength=CriterionStrength.MODERATE,
    )
    assert cr.strength == CriterionStrength.MODERATE


# ---------------------------------------------------------------- ProvisionalClassification

def _met_pvs1():
    return CriterionResult(
        code="PVS1", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC,
        rule_source="ACMG", rule_version="2015", rationale="null variant, established LOF mechanism",
        strength=CriterionStrength.VERY_STRONG,
    )


def test_provisional_classification_rejects_empty_criteria():
    expect_schema_error(
        lambda: ProvisionalClassification(
            provisional_class=ProvisionalClass.VUS, status=ClassificationStatus.PROVISIONAL_AUTOMATED,
            criteria=[], combining_rule_source="ACMG", combining_rule_version="2015", rationale="n/a",
        )
    )


def test_provisional_classification_rejects_final_status_in_milestone_1():
    expect_schema_error(
        lambda: ProvisionalClassification(
            provisional_class=ProvisionalClass.VUS, status=ClassificationStatus.FINAL,
            criteria=[_met_pvs1()], combining_rule_source="ACMG", combining_rule_version="2015", rationale="n/a",
        )
    )


def test_provisional_classification_rejects_duplicate_codes():
    expect_schema_error(
        lambda: ProvisionalClassification(
            provisional_class=ProvisionalClass.LIKELY_PATHOGENIC, status=ClassificationStatus.PROVISIONAL_AUTOMATED,
            criteria=[_met_pvs1(), _met_pvs1()], combining_rule_source="ACMG", combining_rule_version="2015",
            rationale="n/a",
        )
    )


def test_provisional_classification_valid_construction():
    pc = ProvisionalClassification(
        provisional_class=ProvisionalClass.LIKELY_PATHOGENIC, status=ClassificationStatus.PROVISIONAL_AUTOMATED,
        criteria=[_met_pvs1()], combining_rule_source="ACMG", combining_rule_version="2015",
        rationale="1 Very Strong",
    )
    assert pc.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC


# ---------------------------------------------------------------- VariantEvidenceBundle

def _valid_transcript():
    from variant_classifier.models.enums import Consequence
    return TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT,
    )


def _valid_population():
    return PopulationEvidence(
        source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
    )


def test_evidence_bundle_rejects_gene_mismatch():
    expect_schema_error(
        lambda: VariantEvidenceBundle(
            variant=VariantIdentity(variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
            gene_disease_context=GeneDiseaseContext(
                gene="DMD", disease="x", inheritance=Inheritance.X_LINKED_RECESSIVE,
                mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_valid_spec(),
            ),
            transcript_consequences=[_valid_transcript()],
            population_evidence=[_valid_population()],
        )
    )


def test_evidence_bundle_requires_exactly_one_clinically_relevant_transcript():
    from variant_classifier.models.enums import Consequence
    not_relevant = TranscriptConsequence(
        transcript_id="NM_2", clinically_relevant=False, consequence=Consequence.INTRON_VARIANT,
    )
    expect_schema_error(
        lambda: VariantEvidenceBundle(
            variant=VariantIdentity(variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
            gene_disease_context=GeneDiseaseContext(
                gene="CAPN3", disease="x", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
                mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_valid_spec(),
            ),
            transcript_consequences=[not_relevant],
            population_evidence=[_valid_population()],
        )
    )


def test_evidence_bundle_valid_construction():
    bundle = VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="v1", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene="CAPN3", disease="x", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_valid_spec(),
        ),
        transcript_consequences=[_valid_transcript()],
        population_evidence=[_valid_population()],
    )
    assert bundle.variant.variant_id == "v1"


# ---------------------------------------------------------------- GoldenCase

def test_golden_case_rejects_empty_expected_criterion_status():
    expect_schema_error(
        lambda: GoldenCase(
            variant_id="v1", expected_provisional_class=ProvisionalClass.VUS,
            expected_criterion_status={}, source="hand-curated",
        )
    )


def test_golden_case_rejects_unknown_criterion_code():
    expect_schema_error(
        lambda: GoldenCase(
            variant_id="v1", expected_provisional_class=ProvisionalClass.VUS,
            expected_criterion_status={"PZ99": CriterionStatus.NOT_MET}, source="hand-curated",
        )
    )


def test_golden_case_valid_construction():
    gc = GoldenCase(
        variant_id="v1", expected_provisional_class=ProvisionalClass.VUS,
        expected_criterion_status={"PM2": CriterionStatus.MET}, source="hand-curated",
    )
    assert gc.variant_id == "v1"
