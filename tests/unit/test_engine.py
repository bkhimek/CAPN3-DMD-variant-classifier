"""Tests for engine.py — the Milestone 3 combining engine.

The headline test here is test_classify_matches_golden_case_provisional_class:
running the full pipeline (all six evaluators, then combine()) against the
three real curated CAPN3 fixtures should reproduce their golden-case
expected_provisional_class exactly. Those golden cases were written in
Milestone 1, long before any evaluator or engine code existed, so this is
the strongest end-to-end check in the whole test suite — it's the first
time anything in this project goes from raw evidence to a classification.
"""

from variant_classifier import loader
from variant_classifier.engine import classify, combine, evaluate_all
from variant_classifier.models import CriterionResult
from variant_classifier.models.enums import ClassificationStatus, CriterionStatus, CriterionStrength, EvidenceDirection, ProvisionalClass


# ------------------------------------------------------- against real fixtures

def test_classify_matches_golden_case_provisional_class_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()
    thresholds = loader.load_frequency_thresholds()

    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        result = classify(bundle, thresholds)
        assert result.provisional_class == golden.expected_provisional_class, (
            f"{bundle.variant.variant_id}: classified as {result.provisional_class}, "
            f"golden case expects {golden.expected_provisional_class}. Rationale: {result.rationale}"
        )


def test_evaluate_all_returns_exactly_six_criteria_in_fixed_order():
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = bundles[0]
    results = evaluate_all(bundle, thresholds)
    assert [r.code for r in results] == ["PVS1", "PM2", "PP3", "BP4", "BA1", "BS1"]


def test_pathogenic_case_has_no_manual_review_and_no_conflict():
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_SYNTH_PATHOGENIC_01")
    result = classify(bundle, thresholds)
    assert result.provisional_class == ProvisionalClass.PATHOGENIC
    assert result.manual_review_required is False
    assert result.conflicting_evidence_flag is False


def test_founder_case_flags_manual_review_but_not_conflict():
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = classify(bundle, thresholds)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.manual_review_required is True  # PM2 and BS1 are both MANUAL_REVIEW
    assert result.conflicting_evidence_flag is False  # no combining rule satisfied on either side


# ------------------------------------------------------- hand-built combine() cases

def _result(code, status, direction, strength=None):
    return CriterionResult(
        code=code, status=status, direction=direction, strength=strength,
        rule_source="test", rule_version="test", rationale="hand-built test case",
    )


def test_combine_no_evidence_at_all_yields_vus_no_conflict():
    criteria = [
        _result("PVS1", CriterionStatus.NOT_APPLICABLE, EvidenceDirection.PATHOGENIC),
        _result("PM2", CriterionStatus.NOT_MET, EvidenceDirection.PATHOGENIC),
        _result("PP3", CriterionStatus.NOT_EVALUATED, EvidenceDirection.PATHOGENIC),
        _result("BP4", CriterionStatus.NOT_EVALUATED, EvidenceDirection.BENIGN),
        _result("BA1", CriterionStatus.NOT_MET, EvidenceDirection.BENIGN),
        _result("BS1", CriterionStatus.NOT_MET, EvidenceDirection.BENIGN),
    ]
    result = combine(criteria)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.conflicting_evidence_flag is False
    assert result.status == ClassificationStatus.PROVISIONAL_AUTOMATED


def test_combine_two_strong_pathogenic_yields_pathogenic():
    # Table 5: >=2 Strong -> Pathogenic. Using hypothetical codes since no
    # evaluator in this project currently produces two independent Strong
    # pathogenic criteria on the same variant — combine() doesn't care.
    criteria = [
        _result("PS1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.STRONG),
        _result("PS3", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.STRONG),
    ]
    result = combine(criteria)
    assert result.provisional_class == ProvisionalClass.PATHOGENIC


def test_combine_conflicting_evidence_yields_vus_with_conflict_flag():
    # Contrived: enough pathogenic AND enough benign evidence at once.
    criteria = [
        _result("PVS1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.VERY_STRONG),
        _result("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _result("BS1", CriterionStatus.MET, EvidenceDirection.BENIGN, CriterionStrength.STRONG),
        _result("BP4", CriterionStatus.MET, EvidenceDirection.BENIGN, CriterionStrength.SUPPORTING),
    ]
    result = combine(criteria)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.conflicting_evidence_flag is True


def test_combine_three_moderate_yields_likely_pathogenic():
    criteria = [
        _result("PM1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _result("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _result("PM5", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
    ]
    result = combine(criteria)
    assert result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC


def test_combine_manual_review_flag_set_when_any_criterion_needs_review():
    criteria = [
        _result("PVS1", CriterionStatus.NOT_APPLICABLE, EvidenceDirection.PATHOGENIC),
        _result("PM2", CriterionStatus.MANUAL_REVIEW, EvidenceDirection.PATHOGENIC),
    ]
    result = combine(criteria)
    assert result.manual_review_required is True
