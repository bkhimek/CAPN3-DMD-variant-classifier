"""Tests for the PM3 evaluator, added batch 28 -- closing the case-level
gap this project has disclosed since Milestone 4 (see clinical.py's and
engine.py's module docstrings, and pm3_evidence.py's own extensive
docstring for the full design writeup).

Same two-part structure as every other evaluator's tests: (1) run against
every curated fixture that has a PM3 expectation and check it matches the
golden case; (2) hand-built edge cases for branches the curated fixtures
don't happen to cover -- VERY_STRONG, SUPPORTING, and the below-Supporting
NOT_MET band in particular, none of which any real or synthetic curated
fixture in this batch happens to reach (CAPN3_c.550del reaches STRONG,
the two new synthetic fixtures reach MODERATE and the cis-override
NOT_MET) -- plus Pm3ProbandObservation/Pm3Evidence model validation.
"""

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators import evaluate_pm3
from variant_classifier.models import (
    GeneDiseaseContext,
    Pm3Evidence,
    Pm3ProbandObservation,
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
    PhaseRelationship,
    Pm3Zygosity,
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


# ------------------------------------------------------- against real fixtures

def test_pm3_matches_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        if "PM3" not in golden.expected_criterion_status:
            continue
        result = evaluate_pm3(bundle)
        expected = golden.expected_criterion_status["PM3"]
        assert result.status == expected, (
            f"{bundle.variant.variant_id}: PM3 evaluated to {result.status}, golden case "
            f"expects {expected}. Rationale: {result.rationale}"
        )
        checked += 1
    assert checked == len(bundles) > 0


def test_pm3_met_strong_for_capn3_c_550del_real_homozygous_cohort_data():
    # Two independent real published cohorts (Czech, Polish) reporting
    # c.550delA homozygous LGMD2A patients -- see the fixture's notes in
    # data/curated/variant_evidence.json for the full citations. Each
    # capped at the real ACGS 2024 homozygous 1.0-point ceiling, summing
    # to 2.0 = Strong under the CAPN3 VCEP's own real threshold table.
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = evaluate_pm3(bundle)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STRONG


def test_pm3_met_moderate_for_compound_het_confirmed_trans_synthetic_fixture():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_SYNTH_PM3_MODERATE_01")
    result = evaluate_pm3(bundle)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE


def test_pm3_not_met_when_cis_cooccurrence_overrides_synthetic_fixture():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_SYNTH_PM3_CIS_OVERRIDE_01")
    result = evaluate_pm3(bundle)
    assert result.status == CriterionStatus.NOT_MET
    assert "co-occur" in result.rationale


# ------------------------------------------------------- hand-built edge cases

def _spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="2015")


def _population():
    return PopulationEvidence(
        source="gnomAD", source_version="v4.1.0", retrieval_status=PopulationRetrievalStatus.NOT_ASSESSED,
    )


def _bundle(pm3_evidence=None):
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene="CAPN3", genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene="CAPN3", disease="LGMDR1", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_spec(),
        ),
        transcript_consequences=[
            TranscriptConsequence(transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT),
        ],
        population_evidence=[_population()],
        pm3_evidence=pm3_evidence,
    )


def _proband(points, zygosity=Pm3Zygosity.COMPOUND_HETEROZYGOUS, phase=PhaseRelationship.TRANS, proband_id="p1"):
    return Pm3ProbandObservation(
        proband_id=proband_id,
        zygosity=zygosity,
        other_allele_classification=ProvisionalClass.PATHOGENIC,
        points=points,
        phase=phase,
    )


def test_pm3_not_evaluated_when_no_pm3_evidence():
    result = evaluate_pm3(_bundle(pm3_evidence=None))
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pm3_met_very_strong_at_or_above_four_points():
    ev = Pm3Evidence(probands=[_proband(2.0, proband_id="a"), _proband(2.0, proband_id="b")])
    result = evaluate_pm3(_bundle(pm3_evidence=ev))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.VERY_STRONG


def test_pm3_met_supporting_between_half_and_one_point():
    ev = Pm3Evidence(probands=[_proband(0.5)])
    result = evaluate_pm3(_bundle(pm3_evidence=ev))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_pm3_not_met_below_supporting_threshold():
    ev = Pm3Evidence(probands=[_proband(0.25)])
    result = evaluate_pm3(_bundle(pm3_evidence=ev))
    assert result.status == CriterionStatus.NOT_MET
    assert result.strength is None


def test_pm3_cis_cooccurrence_overrides_even_with_no_probands():
    ev = Pm3Evidence(probands=[], cis_cooccurrence_observed=True)
    result = evaluate_pm3(_bundle(pm3_evidence=ev))
    assert result.status == CriterionStatus.NOT_MET


# ------------------------------------------------------- model validation

def test_pm3_proband_rejects_points_above_homozygous_cap():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.HOMOZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=1.01,
    ))


def test_pm3_proband_allows_points_at_exactly_the_homozygous_cap():
    obs = Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.HOMOZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=1.0,
    )
    assert obs.points == 1.0


def test_pm3_proband_rejects_phase_set_on_homozygous():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.HOMOZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=1.0,
        phase=PhaseRelationship.TRANS,
    ))


def test_pm3_proband_rejects_missing_phase_on_compound_het():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.COMPOUND_HETEROZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=1.0,
    ))


def test_pm3_proband_rejects_cis_phase():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.COMPOUND_HETEROZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=1.0,
        phase=PhaseRelationship.CIS,
    ))


def test_pm3_proband_rejects_non_qualifying_other_allele_classification():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.COMPOUND_HETEROZYGOUS,
        other_allele_classification=ProvisionalClass.VUS, points=1.0,
        phase=PhaseRelationship.TRANS,
    ))


def test_pm3_proband_rejects_zero_or_negative_points():
    expect_schema_error(lambda: Pm3ProbandObservation(
        proband_id="p", zygosity=Pm3Zygosity.COMPOUND_HETEROZYGOUS,
        other_allele_classification=ProvisionalClass.PATHOGENIC, points=0.0,
        phase=PhaseRelationship.TRANS,
    ))


def test_pm3_evidence_rejects_empty_probands_with_no_cis_flag():
    expect_schema_error(lambda: Pm3Evidence(probands=[], cis_cooccurrence_observed=False))


def test_pm3_evidence_allows_empty_probands_when_cis_flag_set():
    ev = Pm3Evidence(probands=[], cis_cooccurrence_observed=True)
    assert ev.probands == []
