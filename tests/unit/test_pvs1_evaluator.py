"""Tests for the PVS1 evaluator.

Same two-part structure as the PM2 tests: (1) run it against the three
real curated CAPN3 fixtures and check the result matches the golden case
exactly — these golden-case PVS1 expectations were written in Milestone 1,
before this evaluator existed, so this is a genuine check; (2) hand-built
edge cases for branches the three curated fixtures don't happen to cover
(non-LOF mechanism, no-NMD truncation, splice/start-loss variants).
"""

from variant_classifier import loader
from variant_classifier.evaluators import evaluate_pvs1
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

def test_pvs1_matches_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        if "PVS1" not in golden.expected_criterion_status:
            continue
        result = evaluate_pvs1(bundle)
        expected = golden.expected_criterion_status["PVS1"]
        assert result.status == expected, (
            f"{bundle.variant.variant_id}: PVS1 evaluated to {result.status}, "
            f"golden case expects {expected}. Rationale: {result.rationale}"
        )
        checked += 1
    assert checked == 3  # all three curated cases have a PVS1 expectation


def test_pvs1_frameshift_case_is_very_strong():
    bundles, _ = loader.load_variant_evidence_bundles()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.VERY_STRONG


def test_pvs1_missense_case_is_not_applicable():
    bundles, _ = loader.load_variant_evidence_bundles()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_SYNTH_LIKELY_BENIGN_01")
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.NOT_APPLICABLE


# ------------------------------------------------------- hand-built edge cases

def _spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="2015")


def _population():
    return PopulationEvidence(
        source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
    )


def _bundle(transcript: TranscriptConsequence, mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True):
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene="CAPN3", disease="LGMDR1", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=mechanism, lof_established=lof_established, specification=_spec(),
        ),
        transcript_consequences=[transcript],
        population_evidence=[_population()],
    )


def test_pvs1_not_applicable_when_mechanism_not_loss_of_function():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_GAINED, nmd_predicted=True,
    )
    bundle = _bundle(transcript, mechanism=DiseaseMechanism.GAIN_OF_FUNCTION, lof_established=False)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_pvs1_not_applicable_when_lof_not_established_even_if_mechanism_is_lof():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_GAINED, nmd_predicted=True,
    )
    bundle = _bundle(transcript, mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=False)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_pvs1_no_nmd_yields_manual_review_not_a_guess():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.FRAMESHIFT_VARIANT,
        nmd_predicted=False,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MANUAL_REVIEW
    assert result.strength is None


def test_pvs1_stop_gained_with_nmd_true_is_met():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_GAINED, nmd_predicted=True,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.VERY_STRONG


def test_pvs1_splice_donor_yields_manual_review():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.SPLICE_DONOR_VARIANT,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MANUAL_REVIEW


def test_pvs1_splice_acceptor_yields_manual_review():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.SPLICE_ACCEPTOR_VARIANT,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MANUAL_REVIEW


def test_pvs1_start_lost_yields_manual_review():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.START_LOST,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.MANUAL_REVIEW


def test_pvs1_synonymous_variant_is_not_applicable():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.SYNONYMOUS_VARIANT,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_pvs1_inframe_deletion_is_not_applicable():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_DELETION,
    )
    bundle = _bundle(transcript)
    result = evaluate_pvs1(bundle)
    assert result.status == CriterionStatus.NOT_APPLICABLE
