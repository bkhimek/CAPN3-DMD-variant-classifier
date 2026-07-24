"""Tests for the PP3 and BP4 evaluators — the single-calibrated-score
computational evidence criteria. Same two-part structure as the others.
"""

from variant_classifier import loader
from variant_classifier.evaluators import evaluate_bp4, evaluate_pp3
from variant_classifier.models import (
    ComputationalEvidence,
    GeneDiseaseContext,
    PopulationEvidence,
    Specification,
    TranscriptConsequence,
    VariantEvidenceBundle,
    VariantIdentity,
)
from variant_classifier.models.enums import (
    ComputationalPrediction,
    Consequence,
    CriterionStatus,
    CriterionStrength,
    DiseaseMechanism,
    GenomeBuild,
    Inheritance,
    PopulationRetrievalStatus,
    SpecificationType,
)


# ------------------------------------------------------- against real fixtures

def test_pp3_and_bp4_match_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        for code, evaluate in (("PP3", evaluate_pp3), ("BP4", evaluate_bp4)):
            if code not in golden.expected_criterion_status:
                continue
            result = evaluate(bundle)
            expected = golden.expected_criterion_status[code]
            assert result.status == expected, (
                f"{bundle.variant.variant_id}: {code} evaluated to {result.status}, "
                f"golden case expects {expected}. Rationale: {result.rationale}"
            )
            checked += 1
    assert checked == 10  # PP3 + BP4, across all five curated cases


# ------------------------------------------------------- hand-built edge cases

def _bundle(computational_evidence=None) -> VariantEvidenceBundle:
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
        population_evidence=[
            PopulationEvidence(
                source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
                locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
            )
        ],
        computational_evidence=computational_evidence,
    )


def test_pp3_no_computational_evidence_yields_not_evaluated():
    result = evaluate_pp3(_bundle(computational_evidence=None))
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_bp4_no_computational_evidence_yields_not_evaluated():
    result = evaluate_bp4(_bundle(computational_evidence=None))
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pp3_pathogenic_prediction_yields_met_supporting():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        score=0.9, calibrated_prediction=ComputationalPrediction.PATHOGENIC,
    )
    result = evaluate_pp3(_bundle(computational_evidence=ce))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_bp4_benign_prediction_yields_met_supporting():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        score=0.1, calibrated_prediction=ComputationalPrediction.BENIGN,
    )
    result = evaluate_bp4(_bundle(computational_evidence=ce))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_pp3_indeterminate_prediction_yields_not_met():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        score=0.5, calibrated_prediction=ComputationalPrediction.INDETERMINATE,
    )
    result = evaluate_pp3(_bundle(computational_evidence=ce))
    assert result.status == CriterionStatus.NOT_MET


def test_bp4_indeterminate_prediction_yields_not_met():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test", retrieval_status=PopulationRetrievalStatus.OBSERVED,
        score=0.5, calibrated_prediction=ComputationalPrediction.INDETERMINATE,
    )
    result = evaluate_bp4(_bundle(computational_evidence=ce))
    assert result.status == CriterionStatus.NOT_MET


def test_pp3_not_assessed_yields_not_evaluated():
    ce = ComputationalEvidence(
        tool="REVEL", tool_version="1.0", calibration_source="test", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED,
    )
    result = evaluate_pp3(_bundle(computational_evidence=ce))
    assert result.status == CriterionStatus.NOT_EVALUATED
