"""Tests for bayesian.py -- Milestone 5's Bayesian point-based combining
(Tavtigian et al. 2020), added alongside (not replacing) engine.py's
classic Table 5 combine(). See bayesian.py's module docstring for the
full citation and point-value/threshold quotes this module implements
against.

The headline test is
test_bayesian_matches_hand_derivation_for_all_curated_bundles: running
evaluate_all() (unchanged from engine.py) through combine_bayesian()
against all 22 curated fixtures should reproduce the independently
hand-derived expected_provisional_class and expected_points in
variant_golden_cases_bayesian.yaml exactly -- same golden-case philosophy
as every other headline test in this project, applied to a second
combining system.
"""

from variant_classifier import loader
from variant_classifier.bayesian import BENIGN_POINTS, PATHOGENIC_POINTS, classify_bayesian, combine_bayesian
from variant_classifier.engine import classify, evaluate_all
from variant_classifier.models import CriterionResult
from variant_classifier.models.enums import CriterionStatus, CriterionStrength, EvidenceDirection, ProvisionalClass


def _criterion(code, status, direction=None, strength=None):
    return CriterionResult(
        code=code,
        status=status,
        direction=direction or EvidenceDirection.PATHOGENIC,
        strength=strength,
        rule_source="test fixture",
        rule_version="1",
        rationale="test",
    )


# ------------------------------------------------------- against real fixtures

def test_bayesian_matches_hand_derivation_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    goldens = loader.load_golden_cases_bayesian()

    assert len(goldens) == 27
    for bundle in bundles:
        golden = goldens[bundle.variant.variant_id]
        result = classify_bayesian(bundle, thresholds)
        assert result.provisional_class == golden["expected_provisional_class"], (
            f"{bundle.variant.variant_id}: Bayesian-classified as {result.provisional_class}, "
            f"golden case expects {golden['expected_provisional_class']}. Rationale: {result.rationale}"
        )
        if golden["expected_points"] is not None:
            assert result.points == golden["expected_points"], (
                f"{bundle.variant.variant_id}: {result.points} points, "
                f"golden case expects {golden['expected_points']}"
            )


def test_bayesian_diverges_from_table5_for_exactly_the_five_documented_fixtures():
    """Locks in the full set of divergences found while building this
    milestone, so a future evaluator/threshold change that silently
    creates (or removes) a divergence gets noticed rather than passing
    quietly. Five fixtures, two distinct shapes:

    - CAPN3_c.1939G>T, CAPN3_c.550del: PVS1 Very Strong + one Supporting
      criterion (PM2 for the former, PS3 for the latter, added batch 25)
      = 9 points (Likely Pathogenic), but no Table 5 rule for "1 Very
      Strong + 1 Supporting" alone -> VUS. Both real fixtures, real data.
    - DMD_SYNTH_PATHOGENIC_01, DMD_c.2302C>T, DMD_c.8944C>T: PVS1 Very
      Strong + PM2 Moderate = 10 points (Pathogenic), but Table 5's flat
      PATHOGENIC tier needs "1 Very Strong + >=2 Moderate", so "1 Very
      Strong + 1 Moderate" alone only reaches LIKELY_PATHOGENIC there.
      Two of these three are real fixtures (DMD_c.2302C>T,
      DMD_c.8944C>T).
    """
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()

    expected_divergent_ids = {
        "CAPN3_c.1939G>T",
        "CAPN3_c.550del",
        "DMD_SYNTH_PATHOGENIC_01",
        "DMD_c.2302C>T",
        "DMD_c.8944C>T",
    }

    actual_divergent_ids = set()
    for bundle in bundles:
        t5 = classify(bundle, thresholds)
        bayes = classify_bayesian(bundle, thresholds)
        if t5.provisional_class != bayes.provisional_class:
            actual_divergent_ids.add(bundle.variant.variant_id)

    assert actual_divergent_ids == expected_divergent_ids, (
        f"divergence set changed: now {actual_divergent_ids}, expected {expected_divergent_ids}. "
        "If this is an intentional evaluator/threshold change, update this test and the "
        "corresponding golden-case curator_notes deliberately rather than leaving it failing."
    )


# ------------------------------------------------------- point-value spot checks

def test_pathogenic_point_values_match_tavtigian_2020():
    assert PATHOGENIC_POINTS[CriterionStrength.VERY_STRONG] == 8
    assert PATHOGENIC_POINTS[CriterionStrength.STRONG] == 4
    assert PATHOGENIC_POINTS[CriterionStrength.MODERATE] == 2
    assert PATHOGENIC_POINTS[CriterionStrength.SUPPORTING] == 1


def test_benign_point_values_match_tavtigian_2020():
    assert BENIGN_POINTS[CriterionStrength.STRONG] == -4
    assert BENIGN_POINTS[CriterionStrength.MODERATE] == -2
    assert BENIGN_POINTS[CriterionStrength.SUPPORTING] == -1


# ------------------------------------------------------- hand-built edge cases

def test_ba1_stand_alone_bypasses_point_summing():
    criteria = [
        _criterion("BA1", CriterionStatus.MET, EvidenceDirection.BENIGN, CriterionStrength.STAND_ALONE),
        _criterion("PVS1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.VERY_STRONG),
    ]
    result = combine_bayesian(criteria)
    assert result.provisional_class == ProvisionalClass.BENIGN
    assert "Stand-Alone" in result.rationale


def test_single_very_strong_criterion_alone_is_vus_not_likely_pathogenic():
    # The exact ACGS 2024 example: PVS1_vstr (8 points) alone -> VUS, not
    # Likely Pathogenic, because only one criterion contributed.
    criteria = [_criterion("PVS1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.VERY_STRONG)]
    result = combine_bayesian(criteria)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.points == 8


def test_single_supporting_benign_criterion_alone_is_vus_not_likely_benign():
    # The other exact ACGS 2024 example: BP4_sup (-1 points) alone -> VUS.
    criteria = [_criterion("BP4", CriterionStatus.MET, EvidenceDirection.BENIGN, CriterionStrength.SUPPORTING)]
    result = combine_bayesian(criteria)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.points == -1


def test_two_supporting_criteria_reach_likely_pathogenic_at_exactly_six_points():
    # 1 Very Strong (8) + ... no; use a combination that lands exactly at
    # the 6-point Likely Pathogenic floor: 1 Very Strong is too big alone
    # (needs a second criterion per the 2-minimum rule) -- use Moderate x3
    # instead (2+2+2=6, 3 contributing criteria).
    criteria = [
        _criterion("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _criterion("PM4", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _criterion("PM5", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
    ]
    result = combine_bayesian(criteria)
    assert result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC
    assert result.points == 6


def test_ten_points_is_pathogenic_not_likely_pathogenic():
    criteria = [
        _criterion("PVS1", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.VERY_STRONG),
        _criterion("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
    ]
    result = combine_bayesian(criteria)
    assert result.points == 10
    assert result.provisional_class == ProvisionalClass.PATHOGENIC


def test_no_met_criteria_is_vus_with_zero_points():
    criteria = [
        _criterion("PVS1", CriterionStatus.NOT_APPLICABLE),
        _criterion("PM2", CriterionStatus.NOT_MET),
    ]
    result = combine_bayesian(criteria)
    assert result.provisional_class == ProvisionalClass.VUS
    assert result.points == 0


def test_conflicting_evidence_flag_is_always_false_for_bayesian():
    # Unlike Table 5's combine(), a net point sum can't land in two bands
    # at once, so there's no "conflict" state to flag -- even when both
    # directions had MET criteria (net sum just resolves it).
    criteria = [
        _criterion("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.MODERATE),
        _criterion("BS1", CriterionStatus.MET, EvidenceDirection.BENIGN, CriterionStrength.STRONG),
    ]
    result = combine_bayesian(criteria)
    assert result.conflicting_evidence_flag is False
    assert result.points == 2 + (-4)
    assert "both directions had MET criteria" in result.rationale


def test_manual_review_required_propagates_same_as_table5():
    criteria = [
        _criterion("PVS1", CriterionStatus.MANUAL_REVIEW),
        _criterion("PM2", CriterionStatus.MET, EvidenceDirection.PATHOGENIC, CriterionStrength.SUPPORTING),
    ]
    result = combine_bayesian(criteria)
    assert result.manual_review_required is True


def test_classify_bayesian_uses_same_evaluate_all_as_classify():
    # Confirms the two combining systems really do share the exact same
    # per-criterion evidence -- classify_bayesian() isn't quietly running
    # a different evaluation pipeline.
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    bundle = bundles[0]

    direct = evaluate_all(bundle, thresholds)
    via_table5 = classify(bundle, thresholds).criteria
    via_bayesian = classify_bayesian(bundle, thresholds).criteria

    assert [c.code for c in direct] == [c.code for c in via_table5] == [c.code for c in via_bayesian]
    assert [c.status for c in direct] == [c.status for c in via_table5] == [c.status for c in via_bayesian]
