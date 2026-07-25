"""Tests for the PM2 evaluator.

Two kinds of coverage: (1) run it against the three real curated CAPN3
fixtures and check the result matches the corresponding golden case
exactly — the fixtures and golden cases were written before this
evaluator existed, so this is a genuine check, not a tautology; (2)
hand-built edge cases for retrieval-status branches the three curated
fixtures don't happen to exercise.
"""

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators import evaluate_pm2
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

def test_pm2_matches_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()
    thresholds = loader.load_frequency_thresholds()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        if "PM2" not in golden.expected_criterion_status:
            continue
        result = evaluate_pm2(bundle, thresholds)
        expected = golden.expected_criterion_status["PM2"]
        assert result.status == expected, (
            f"{bundle.variant.variant_id}: PM2 evaluated to {result.status}, "
            f"golden case expects {expected}. Rationale: {result.rationale}"
        )
        checked += 1
    # Tied to len(bundles) rather than a hardcoded number: every curated
    # bundle is expected to have a PM2 expectation, and this way adding
    # fixtures later doesn't require remembering to update a count too (a
    # maintenance gap that bit this test suite during Milestone 4).
    assert checked == len(bundles) > 0


def test_pm2_founder_case_is_now_not_met_under_real_vcep_threshold():
    # This fixture (CAPN3_c.550del) is the case that originally motivated
    # the evaluator's ancestry-AF MANUAL_REVIEW branch, back when CAPN3
    # used a placeholder PM2 threshold (0.001). Since batch 4 adopted the
    # real ClinGen LGMD VCEP CAPN3 threshold (<0.0001), this variant's
    # overall AF (0.00023) now exceeds the threshold on its own, so PM2 is
    # decided as NOT_MET before the ancestry-specific check is even
    # reached -- see variant_golden_cases.yaml's curator_note for the full
    # explanation. The ancestry-AF branch itself is still real code with
    # real behavior; test_pm2_ancestry_enrichment_below_threshold_is_manual_review
    # below exercises it directly with hand-built data now that this real
    # fixture no longer does.
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = evaluate_pm2(bundle, thresholds)
    assert result.status == CriterionStatus.NOT_MET


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
            TranscriptConsequence(
                transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT,
            )
        ],
        population_evidence=[pe],
    )


def _thresholds():
    return {
        "ba1_stand_alone_af": 0.05,
        "genes": {"CAPN3": {"pm2_max_credible_af": 0.001, "bs1_min_af": 0.001, "threshold_source": "test"}},
    }


def test_pm2_not_assessed_yields_not_evaluated():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED)
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pm2_unavailable_yields_not_evaluated():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.UNAVAILABLE)
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pm2_not_applicable_yields_not_evaluated():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.NOT_APPLICABLE)
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pm2_absent_yields_met_moderate():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
            locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
        )
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.MET
    from variant_classifier.models.enums import CriterionStrength
    assert result.strength == CriterionStrength.MODERATE


def test_pm2_observed_above_threshold_yields_not_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
            overall_af=0.05,
        )
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.NOT_MET


def test_pm2_observed_below_threshold_no_ancestry_data_yields_met():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
            overall_af=0.00001,
        )
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.MET


def test_pm2_missing_gene_threshold_raises():
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
            overall_af=0.00001,
        )
    )
    expect_schema_error(lambda: evaluate_pm2(bundle, {}))


def test_pm2_ancestry_enrichment_below_threshold_is_manual_review():
    # Direct, hand-built coverage of the ancestry-AF MANUAL_REVIEW branch,
    # using the local _thresholds() (pm2_max_credible_af=0.001) rather than
    # the real CAPN3 config -- kept independent of config/population_thresholds.yaml
    # so this test does not silently stop covering the branch if the real
    # threshold changes again later. Overall AF is below 0.001, but the
    # ancestry-specific AF clears it -- the founder-enrichment pattern.
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.OBSERVED,
            overall_af=0.0005, ancestry_specific_max_af=0.002,
        )
    )
    result = evaluate_pm2(bundle, _thresholds())
    assert result.status == CriterionStatus.MANUAL_REVIEW
    assert "ancestry" in result.rationale.lower()


def test_pm2_gene_specific_strength_is_honored_when_configured():
    # Covers the batch-4 addition: pm2_strength is read per gene, defaulting
    # to MODERATE, but a gene can be configured (as CAPN3 now is, via the
    # real ClinGen LGMD VCEP spec) to SUPPORTING instead.
    from variant_classifier.models.enums import CriterionStrength

    thresholds = {
        "ba1_stand_alone_af": 0.05,
        "genes": {
            "CAPN3": {
                "pm2_max_credible_af": 0.0001,
                "pm2_strength": "SUPPORTING",
                "bs1_min_af": 0.001,
                "threshold_source": "test",
            }
        },
    }
    bundle = _bundle_with_population_evidence(
        PopulationEvidence(
            source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
            locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
        )
    )
    result = evaluate_pm2(bundle, thresholds)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_pm2_gene_without_strength_override_defaults_to_moderate():
    # A gene present in config but with no pm2_strength key (e.g. DMD, which
    # has no VCEP spec adopted here) should keep the generic-ACMG default.
    thresholds = {
        "ba1_stand_alone_af": 0.05,
        "genes": {"DMD": {"pm2_max_credible_af": 0.001, "bs1_min_af": 0.001, "threshold_source": "test"}},
    }
    bundle_pe = PopulationEvidence(
        source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True, allele_count=0, allele_number=1000000,
    )
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
    from variant_classifier.models.enums import CriterionStrength

    result = evaluate_pm2(bundle, thresholds)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE


def test_pm2_rejects_multiple_population_evidence_entries():
    bundle = VariantEvidenceBundle(
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
            PopulationEvidence(source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED),
            PopulationEvidence(source="TOPMed", source_version="freeze8", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED),
        ],
    )
    expect_schema_error(lambda: evaluate_pm2(bundle, _thresholds()))
