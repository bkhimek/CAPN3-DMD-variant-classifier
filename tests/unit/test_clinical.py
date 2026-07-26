"""Tests for clinical.py — Milestone 4's case-level interpretation layer.

The headline test, test_interpret_case_matches_golden_case_for_all_curated_cases,
runs the full pipeline (evaluate every variant a case references via
engine.classify(), then clinical.interpret_case()) against all six curated
ClinicalCase fixtures and checks the result against
case_interpretation_golden_cases.yaml — golden expectations written
independently of this code, same philosophy as every other golden-case
test in this project.
"""

from variant_classifier import loader
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

    assert len(cases) == 8
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


def test_interpret_case_dispatches_not_applicable_for_unhandled_inheritance():
    context = GeneDiseaseContext(
        gene="SOMEGENE", disease="some dominant disease", inheritance=Inheritance.AUTOSOMAL_DOMINANT,
        mechanism=DiseaseMechanism.GAIN_OF_FUNCTION, lof_established=False,
        specification=Specification(type=SpecificationType.GENERIC_ACMG, version="2015"),
    )
    case = ClinicalCase(case_id="c1", gene="SOMEGENE", karyotypic_sex=KaryotypicSex.XY, variant_ids=["v1"])
    result = interpret_case(case, {}, context)
    assert result.status == CaseInterpretationStatus.NOT_APPLICABLE
