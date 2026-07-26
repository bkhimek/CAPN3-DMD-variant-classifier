"""Tests for the PM4 evaluator, added batch 14.

Same two-part structure as every other evaluator's tests: (1) run it
against every curated fixture that has a PM4 expectation and check the
result matches the golden case; (2) hand-built edge cases for branches
the curated fixtures don't happen to cover (repeat-region exclusion,
single-residue downgrade, unrecorded size, non-PM4 consequence types).
"""

from variant_classifier import loader
from variant_classifier.evaluators import evaluate_pm4
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


# ------------------------------------------------------- against real fixtures

def test_pm4_matches_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        if "PM4" not in golden.expected_criterion_status:
            continue
        result = evaluate_pm4(bundle)
        expected = golden.expected_criterion_status["PM4"]
        assert result.status == expected, (
            f"{bundle.variant.variant_id}: PM4 evaluated to {result.status}, "
            f"golden case expects {expected}. Rationale: {result.rationale}"
        )
        checked += 1
    assert checked == len(bundles) > 0


# ------------------------------------------------------- hand-built edge cases

def _spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="2015")


def _population():
    return PopulationEvidence(
        source="gnomAD", source_version="test", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True,
    )


def _bundle(transcript: TranscriptConsequence, gene="CAPN3", mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True):
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene=gene, genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene=gene, disease="test disease", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=mechanism, lof_established=lof_established, specification=_spec(),
        ),
        transcript_consequences=[transcript],
        population_evidence=[_population()],
    )


def test_pm4_missense_variant_is_not_applicable():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT,
    )
    result = evaluate_pm4(_bundle(transcript))
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_pm4_repeat_region_indel_is_not_met():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_DELETION,
        repeat_region=True,
    )
    result = evaluate_pm4(_bundle(transcript))
    assert result.status == CriterionStatus.NOT_MET


def test_pm4_single_residue_indel_is_supporting_not_moderate():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_DELETION,
        repeat_region=False, protein_length_change_aa=1,
    )
    result = evaluate_pm4(_bundle(transcript))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_pm4_multi_residue_indel_is_moderate():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_INSERTION,
        repeat_region=False, protein_length_change_aa=3,
    )
    result = evaluate_pm4(_bundle(transcript))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE


def test_pm4_stop_loss_with_unrecorded_size_defaults_to_moderate_not_a_guess():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_LOST,
        repeat_region=False,
    )
    result = evaluate_pm4(_bundle(transcript))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE
    assert "not curated" in result.rationale or "unrecorded" in result.rationale


def test_pm4_does_not_require_lof_mechanism_unlike_pvs1():
    # PM4 is not gated on disease mechanism the way PVS1 is (Richards et al.
    # 2015 doesn't condition it on mechanism) — confirm it still evaluates
    # normally even when lof_established=False / mechanism is GAIN_OF_FUNCTION.
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_DELETION,
        repeat_region=False, protein_length_change_aa=2,
    )
    result = evaluate_pm4(_bundle(transcript, mechanism=DiseaseMechanism.GAIN_OF_FUNCTION, lof_established=False))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE


def test_transcript_consequence_requires_repeat_region_for_pm4_relevant_consequence():
    from variant_classifier.errors import SchemaValidationError
    try:
        TranscriptConsequence(
            transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_LOST,
        )
        raise AssertionError("expected SchemaValidationError, none was raised")
    except SchemaValidationError:
        pass
