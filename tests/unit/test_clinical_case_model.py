"""Schema-validation tests for the two Milestone 4 models: ClinicalCase
and CaseInterpretation. Same style as test_models.py.
"""

from variant_classifier.errors import SchemaValidationError
from variant_classifier.models import CaseInterpretation, ClinicalCase
from variant_classifier.models.enums import CaseInterpretationStatus, KaryotypicSex, PhaseRelationship


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


def test_clinical_case_valid_single_variant():
    case = ClinicalCase(case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XY, variant_ids=["v1"])
    assert case.phase is None


def test_clinical_case_valid_two_variants_with_phase():
    case = ClinicalCase(
        case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XX,
        variant_ids=["v1", "v2"], phase=PhaseRelationship.TRANS,
    )
    assert case.phase == PhaseRelationship.TRANS


def test_clinical_case_rejects_empty_variant_ids():
    expect_schema_error(lambda: ClinicalCase(case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XY, variant_ids=[]))


def test_clinical_case_rejects_more_than_two_variant_ids():
    expect_schema_error(
        lambda: ClinicalCase(
            case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XY,
            variant_ids=["v1", "v2", "v3"], phase=PhaseRelationship.TRANS,
        )
    )


def test_clinical_case_rejects_duplicate_variant_ids():
    expect_schema_error(
        lambda: ClinicalCase(
            case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XX,
            variant_ids=["v1", "v1"], phase=PhaseRelationship.TRANS,
        )
    )


def test_clinical_case_rejects_two_variants_without_phase():
    expect_schema_error(
        lambda: ClinicalCase(case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XX, variant_ids=["v1", "v2"])
    )


def test_clinical_case_rejects_phase_with_single_variant():
    expect_schema_error(
        lambda: ClinicalCase(
            case_id="c1", gene="CAPN3", karyotypic_sex=KaryotypicSex.XY,
            variant_ids=["v1"], phase=PhaseRelationship.UNKNOWN,
        )
    )


def test_clinical_case_from_dict_round_trip():
    data = {
        "case_id": "c1", "gene": "DMD", "karyotypic_sex": "XY", "variant_ids": ["v1"],
    }
    case = ClinicalCase.from_dict(data)
    assert case.case_id == "c1"
    assert case.karyotypic_sex == KaryotypicSex.XY


def test_case_interpretation_valid_construction():
    ci = CaseInterpretation(case_id="c1", gene="CAPN3", status=CaseInterpretationStatus.EXPLAINED, rationale="test")
    assert ci.status == CaseInterpretationStatus.EXPLAINED


def test_case_interpretation_rejects_empty_rationale():
    expect_schema_error(
        lambda: CaseInterpretation(case_id="c1", gene="CAPN3", status=CaseInterpretationStatus.EXPLAINED, rationale="  ")
    )
