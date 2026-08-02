"""Tests for clinical.py — Milestone 4's case-level interpretation layer.

The headline test, test_interpret_case_matches_golden_case_for_all_curated_cases,
runs the full pipeline (evaluate every variant a case references via
engine.classify(), then clinical.interpret_case()) against all curated
ClinicalCase fixtures and checks the result against
case_interpretation_golden_cases.yaml — golden expectations written
independently of this code, same philosophy as every other golden-case
test in this project.
"""

from variant_classifier import loader
from variant_classifier.bayesian import classify_bayesian
from variant_classifier.clinical import interpret_case, interpret_recessive_case, interpret_x_linked_case
from variant_classifier.engine import classify
from variant_classifier.errors import SchemaValidationError
from variant_classifier.models import ClinicalCase, GeneDiseaseContext, Specification
from variant_classifier.models.enums import (
    CaseInterpretationStatus,
    DiseaseMechanism,
    Inheritance,
    KaryotypicSex,
    PhaseRelationship,
    SpecificationType,
)


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


def _classify_all_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    return {b.variant.variant_id: classify(b, thresholds) for b in bundles}


# ------------------------------------------------------- against real fixtures

def test_interpret_case_matches_golden_case_for_all_curated_cases():
    classifications = _classify_all_bundles()
    contexts = loader.load_gene_disease_contexts()
    cases = loader.load_clinical_cases()
    goldens = loader.load_case_interpretation_goldens()

    assert len(cases) == 11
    for case in cases:
        golden = goldens[case.case_id]
        result = interpret_case(case, classifications, contexts[case.gene])
        assert result.status == golden["expected_status"], (
            f"{case.case_id}: interpreted as {result.status}, golden case expects "
            f"{golden['expected_status']}. Rationale: {result.rationale}"
        )


def test_capn3_trans_vs_cis_differ_only_by_phase():
    # The specific pair this fixture set was built to prove: identical
    # variants, identical individual classifications, opposite case-level
    # outcome purely because of phase.
    classifications = _classify_all_bundles()
    contexts = loader.load_gene_disease_contexts()
    cases = {c.case_id: c for c in loader.load_clinical_cases()}

    trans_result = interpret_case(cases["CASE_CAPN3_BIALLELIC_TRANS"], classifications, contexts["CAPN3"])
    cis_result = interpret_case(cases["CASE_CAPN3_BIALLELIC_CIS"], classifications, contexts["CAPN3"])
    assert trans_result.status == CaseInterpretationStatus.EXPLAINED
    assert cis_result.status == CaseInterpretationStatus.INSUFFICIENT


def test_dmd_male_vs_female_differ_only_by_karyotypic_sex():
    classifications = _classify_all_bundles()
    contexts = loader.load_gene_disease_contexts()
    cases = {c.case_id: c for c in loader.load_clinical_cases()}

    male_result = interpret_case(cases["CASE_DMD_HEMIZYGOUS_MALE"], classifications, contexts["DMD"])
    female_result = interpret_case(cases["CASE_DMD_FEMALE_CARRIER"], classifications, contexts["DMD"])
    assert male_result.status == CaseInterpretationStatus.EXPLAINED
    assert female_result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_dmd_male_vs_female_differ_only_by_karyotypic_sex_real_variant():
    # Same proof as above, but end-to-end on a real ClinVar-sourced DMD
    # variant (DMD_c.2302C>T) instead of the synthetic one — added batch 12,
    # the project's first real-variant Milestone 4 clinical case.
    classifications = _classify_all_bundles()
    contexts = loader.load_gene_disease_contexts()
    cases = {c.case_id: c for c in loader.load_clinical_cases()}

    male_result = interpret_case(cases["CASE_DMD_HEMIZYGOUS_MALE_REAL"], classifications, contexts["DMD"])
    female_result = interpret_case(cases["CASE_DMD_FEMALE_CARRIER_REAL"], classifications, contexts["DMD"])
    assert male_result.status == CaseInterpretationStatus.EXPLAINED
    assert female_result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_clinical_case_interpretation_agnostic_to_combining_system():
    """Added batch 20 (Milestone 5), alongside bayesian.py.

    clinical.py's functions all take a pre-computed `classifications` dict
    as a parameter -- they never call engine.classify() themselves (see
    clinical.py's own module docstring: "It does not re-derive
    variant-level evidence; it consumes engine.classify()'s output").
    That's a deliberate separation of concerns, and this test is the
    concrete proof of it: feed every curated ClinicalCase the exact same
    real evidence, but classified via classify_bayesian() instead of
    classify(), and every case-level result should be identical to the
    Table-5-fed golden expectation.

    This matters concretely for DMD_c.2302C>T and DMD_c.8944C>T, both of
    which move from LIKELY_PATHOGENIC (Table 5) to PATHOGENIC (Bayesian) --
    see variant_golden_cases_bayesian.yaml. clinical.py's own
    _QUALIFYING set treats both tiers identically for case-level purposes,
    so CASE_DMD_HEMIZYGOUS_MALE_REAL still resolves EXPLAINED either way,
    and CASE_DMD_FEMALE_CARRIER_REAL still resolves MANUAL_REVIEW either
    way (karyotypic sex, not variant tier, drives that one). Nothing here
    required touching case_interpretation_golden_cases.yaml.
    """
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    contexts = loader.load_gene_disease_contexts()
    cases = loader.load_clinical_cases()
    goldens = loader.load_case_interpretation_goldens()

    bayesian_classifications = {b.variant.variant_id: classify_bayesian(b, thresholds) for b in bundles}

    for case in cases:
        golden = goldens[case.case_id]
        result = interpret_case(case, bayesian_classifications, contexts[case.gene])
        assert result.status == golden["expected_status"], (
            f"{case.case_id}: Bayesian-fed interpretation is {result.status}, golden case (written "
            f"against Table 5) expects {golden['expected_status']} -- if this genuinely differs by "
            "combining system, that's a real finding worth its own fixture, not a silent mismatch. "
            f"Rationale: {result.rationale}"
        )


# ------------------------------------------------------- hand-built edge cases

def _capn3_context():
    return GeneDiseaseContext(
        gene="CAPN3", disease="LGMDR1", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
        mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True,
        specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
    )


def _dmd_context():
    return GeneDiseaseContext(
        gene="DMD", disease="dystrophinopathy", inheritance=Inheritance.X_LINKED_RECESSIVE,
        mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True,
        specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
    )


def test_interpret_recessive_case_rejects_non_recessive_context():
    case = ClinicalCase(case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.XY, variant_ids=["v1"])
    expect_schema_error(lambda: interpret_recessive_case(case, {}, _dmd_context()))


def test_interpret_x_linked_case_rejects_non_x_linked_context():
    case = ClinicalCase(case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XY, variant_ids=["v1"])
    expect_schema_error(lambda: interpret_x_linked_case(case, {}, _capn3_context()))


def test_interpret_x_linked_case_rejects_two_variants():
    case = ClinicalCase(
        case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.XY,
        variant_ids=["v1", "v2"], phase=PhaseRelationship.TRANS,
    )
    expect_schema_error(lambda: interpret_x_linked_case(case, {}, _dmd_context()))


def test_dmd_xx_biallelic_trans_vs_cis_differ_only_by_phase():
    # Batch 29's parallel to test_capn3_trans_vs_cis_differ_only_by_phase --
    # except, unlike autosomal recessive, cis does NOT resolve to
    # INSUFFICIENT for X-linked XX (see clinical.py's
    # _interpret_xx_biallelic docstring for why: a wild-type X copy is not
    # always active the way a wild-type autosome copy is).
    classifications = _classify_all_bundles()
    contexts = loader.load_gene_disease_contexts()
    cases = {c.case_id: c for c in loader.load_clinical_cases()}

    trans_result = interpret_case(cases["CASE_DMD_XX_BIALLELIC_TRANS"], classifications, contexts["DMD"])
    cis_result = interpret_case(cases["CASE_DMD_XX_BIALLELIC_CIS"], classifications, contexts["DMD"])
    assert trans_result.status == CaseInterpretationStatus.EXPLAINED
    assert cis_result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_interpret_x_linked_case_xx_unknown_phase_is_manual_review():
    case = ClinicalCase(
        case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.XX,
        variant_ids=["DMD_c.2302C>T", "DMD_c.8944C>T"], phase=PhaseRelationship.UNKNOWN,
    )
    classifications = _classify_all_bundles()
    result = interpret_x_linked_case(case, classifications, _dmd_context())
    assert result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_interpret_x_linked_case_xx_trans_but_not_both_qualifying_is_manual_review():
    # DMD_c.10103A>G is real and VUS-classified (see CASE_DMD_HEMIZYGOUS_MALE_VUS_REAL) --
    # trans with a genuinely Likely Pathogenic partner still isn't enough.
    case = ClinicalCase(
        case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.XX,
        variant_ids=["DMD_c.2302C>T", "DMD_c.10103A>G"], phase=PhaseRelationship.TRANS,
    )
    classifications = _classify_all_bundles()
    result = interpret_x_linked_case(case, classifications, _dmd_context())
    assert result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_interpret_x_linked_case_rejects_two_variants_for_other_karyotype():
    case = ClinicalCase(
        case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.OTHER,
        variant_ids=["v1", "v2"], phase=PhaseRelationship.TRANS,
    )
    expect_schema_error(lambda: interpret_x_linked_case(case, {}, _dmd_context()))


def test_interpret_x_linked_case_xx_single_variant_still_manual_review():
    # Unchanged since Milestone 4 -- batch 29 only added the two-variant path.
    case = ClinicalCase(case_id="c1", gene="DMD", karyotypic_sex=KaryotypicSex.XX, variant_ids=["DMD_c.2302C>T"])
    classifications = _classify_all_bundles()
    result = interpret_x_linked_case(case, classifications, _dmd_context())
    assert result.status == CaseInterpretationStatus.MANUAL_REVIEW


def test_interpret_case_dispatches_not_applicable_for_unhandled_inheritance():
    context = GeneDiseaseContext(
        gene="SOMEGENE", disease="some dominant disease", inheritance=Inheritance.AUTOSOMAL_DOMINANT,
        mechanism=DiseaseMechanism.GAIN_OF_FUNCTION, lof_established=False,
        specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
    )
    case = ClinicalCase(case_id="c1", gene="SOMEGENE", karyotypic_sex=KaryotypicSex.XY, variant_ids=["v1"])
    result = interpret_case(case, {}, context)
    assert result.status == CaseInterpretationStatus.NOT_APPLICABLE
