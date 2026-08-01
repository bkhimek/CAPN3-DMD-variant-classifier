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
from variant_classifier.bayesian import classify_bayesian
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


def test_no_real_capn3_variant_currently_reaches_pathogenic_tier():
    """Regression/documentation test, added batch 13.

    This is a structural consequence of two facts, not a coincidence of
    which fixtures happen to be curated: (1) the real ClinGen LGMD VCEP
    threshold adopted for CAPN3 (batch 4) fixes PM2 at SUPPORTING
    strength, never MODERATE; (2) PP3's strength is hardcoded SUPPORTING
    everywhere in this engine (pp3.py has no gene-specific override). For
    a loss-of-function CAPN3 variant, the only pathogenic-direction
    criteria that can be MET are PVS1 (Very Strong) and PM2 (Supporting) —
    Table 5 has no rule for "1 Very Strong + 1 Supporting" alone (see
    CAPN3_c.1939G>T). For a missense CAPN3 variant, only PM2 and PP3 can
    be MET, both Supporting — Table 5 has no rule for "2 Supporting"
    alone without at least one Moderate or Strong criterion either. So,
    as currently configured, no real CAPN3 variant — however strong its
    individual evidence — can reach PATHOGENIC or LIKELY_PATHOGENIC
    through this engine. This test locks that fact in as an explicit,
    intentional claim rather than an implicit one: if it ever starts
    failing (e.g. a gene-specific PM2 strength override, or a
    higher-strength PP3), that's a real, deliberate change to notice and
    document — not something that should happen silently.

    Extended in batch 14 with a third shape, in-frame indel/stop-loss
    CAPN3 variants: only PM2 and PM4 can be MET (PVS1 and PP3 are both
    out of scope for that consequence class). This evaluator's PM4
    defaults a single-residue change to Supporting, so the best case is
    "2 Supporting" — same insufficient shape as the missense case, and
    still no Table 5 rule for it. Even in the case where PM4 landed at
    its unmodified Moderate strength (as the real ClinGen LGMD VCEP
    appears to have applied it for CAPN3_c.1401_1403del — see that
    fixture's curator_note), "1 Moderate + 1 Supporting" alone still
    doesn't satisfy any Likely-Pathogenic rule in Table 5. So PM4's
    addition does not create a new way for a real CAPN3 variant to clear
    this bar; it was checked, not assumed.
    """
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    real_capn3_ids = [
        b.variant.variant_id for b in bundles
        if b.variant.gene == "CAPN3" and "SYNTH" not in b.variant.variant_id
    ]
    assert len(real_capn3_ids) >= 5, "expected several real CAPN3 fixtures to check this against"

    for variant_id in real_capn3_ids:
        expected = golden_cases[variant_id].expected_provisional_class
        assert expected not in (ProvisionalClass.PATHOGENIC, ProvisionalClass.LIKELY_PATHOGENIC), (
            f"{variant_id}: golden case now expects {expected}, contradicting this test's "
            "documented claim that no real CAPN3 variant can currently reach that tier — if "
            "this is an intentional config/evaluator change, update or remove this test "
            "deliberately rather than leaving it failing."
        )


def test_capn3_c1939_reaches_likely_pathogenic_under_bayesian_combining():
    """Closes the loop test_no_real_capn3_variant_currently_reaches_pathogenic_tier's
    own docstring anticipated, added batch 20 (Milestone 5).

    That test's claim is scoped to *this engine as currently configured* --
    i.e. Table 5 combining (engine.combine()). It explicitly names one of
    the two ways that scoping could someday change: "a gene-specific PM2
    strength override, or a higher-strength PP3." Milestone 5 adds a third
    way that test's docstring didn't anticipate by name but its own
    reasoning already implied: keep every evaluator and threshold
    identical, and change the *combining system* instead. CAPN3_c.1939G>T
    (PVS1 Very Strong + PM2 Supporting) has no Table 5 rule and stays VUS
    under classify() -- confirmed unchanged by
    test_no_real_capn3_variant_currently_reaches_pathogenic_tier still
    passing -- but reaches 9 points (Likely Pathogenic, 6-9) under
    classify_bayesian(), because Tavtigian et al. 2020's point system does
    define a combination for "1 Very Strong + 1 Supporting" that Table 5
    simply never enumerated. This is not a contradiction between the two
    tests: one says "Table 5, as configured, cannot reach this tier for
    real CAPN3 evidence"; this one says "a different, equally real
    combining system, given the identical evidence, can." Both are true
    at once, which is exactly bayesian.py's own point.
    """
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.1939G>T")

    table5_result = classify(bundle, thresholds)
    bayesian_result = classify_bayesian(bundle, thresholds)

    assert table5_result.provisional_class == ProvisionalClass.VUS
    assert bayesian_result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC
    assert bayesian_result.points == 9


def test_evaluate_all_returns_exactly_nine_criteria_in_fixed_order():
    # Was "exactly six" through Milestone 3; PM4 added batch 14; PS1/PM5 added this round.
    bundles, _ = loader.load_variant_evidence_bundles()
    thresholds = loader.load_frequency_thresholds()
    bundle = bundles[0]
    results = evaluate_all(bundle, thresholds)
    assert [r.code for r in results] == ["PVS1", "PM2", "PM4", "PS1", "PM5", "PP3", "BP4", "BA1", "BS1"]


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
    # As of batch 4 (real ClinGen LGMD VCEP CAPN3 thresholds), it's BA1 and
    # BS1 that are MANUAL_REVIEW here, not PM2 and BS1 -- PM2 is now decided
    # (NOT_MET) since this variant's overall AF alone exceeds the real,
    # stricter PM2 threshold. See variant_golden_cases.yaml's curator_note.
    assert result.manual_review_required is True  # BA1 and BS1 are both MANUAL_REVIEW
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
