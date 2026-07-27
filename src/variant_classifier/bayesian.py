"""Milestone 5 -- Bayesian point-based combining (Tavtigian et al. 2020),
"Fitting a naturally scaled point system to the ACMG/AMP variant
classification guidelines," Human Mutation 41(6):1023-1041.

Offered ALONGSIDE engine.py's classic Table 5 combine(), not as a
replacement -- both are real, and this project's own design notes (see
README, "Real ClinGen LGMD VCEP thresholds for CAPN3") already document a
concrete case where they disagree: CAPN3_c.1939G>T (PVS1 Very Strong + PM2
Supporting) lands VUS under Table 5, because Table 5 has no combining rule
for "1 Very Strong + 1 Supporting" alone, but reaches 9 points (Likely
Pathogenic, 6-9) under this module -- the exact discrepancy this milestone
was built to make concrete and testable rather than just described in
prose.

Point values and thresholds below are quoted directly from the ACGS 2024
UK Practice Guidelines for Variant Classification (v1, Aug 2024), which
cites Tavtigian et al. 2020 verbatim:

    "Criteria can be easily combined using Bayesian-derived evidence
    points for pathogenicity (Very Strong=8, Strong=4, Moderate=2,
    Supporting=1) or benignity (Strong=-4, Moderate=-2, Supporting=-1).
    Thresholds are: >=10 (Pathogenic), 6-9 (Likely Pathogenic), 0-5
    (VUS), -1 to -5 (Likely Benign), <=-6 (Benign) (Tavtigian et al.,
    2020). With the exception of BA1 stand-alone, a minimum of two
    criteria are required to classify a variant as (likely) benign or
    (likely) pathogenic; therefore, variants with only one piece of
    evidence e.g. PVS1_vstr (8 points) or BP4_sup (-1 points) are
    classified as a VUS pending a second corroborating piece of
    evidence."

Two things worth naming explicitly, both carried over deliberately rather
than re-derived from scratch:

- BA1 Stand-Alone bypasses point-summing entirely and resolves straight
  to BENIGN, exactly like Table 5's own treatment (engine.py's
  _benign_tier() has the identical "sa >= 1 -> BENIGN" short-circuit).
  This project's real BA1 fixture (DMD_c.5234G>A, batch 18) is expected
  to land on BENIGN under both systems for exactly this reason.
- The "minimum of two criteria" rule reproduces, as a side effect, a
  design decision this project already made independently for Table 5:
  the "Expanding the curated set" note for CAPN3_SYNTH_LIKELY_BENIGN_01
  and the DMD_c.5163G>C golden-case curator_note both explain that BS1
  or BP4 alone can never reach LIKELY_BENIGN. That's not a coincidence --
  it's the same real-world combining constraint showing up in both
  systems, and this module's own fixtures re-confirm it rather than
  silently assuming it carries over.

One further discrepancy, documented but not currently exercisable by this
project's evaluators: Tavtigian et al. 2020 itself reports that Table 5's
"1 Strong + >=1 Strong" (i.e. >=2 Strong) pathogenic rule is the weakest of
Table 5's eight Pathogenic combining paths, with a posterior probability
of 0.975 versus >0.99 for the others -- meaning 2 Strong criteria (4+4=8
points) fall short of this module's own >=10 Pathogenic threshold, landing
Likely Pathogenic instead, a real published inconsistency between Table 5
and the naturally-scaled point system rather than a bug introduced here.
This project has no Strong-strength pathogenic-direction evaluator (PVS1
is Very Strong; PM2/PM4 are Moderate or Supporting; PP3 is Supporting) so
no current fixture can actually reach this discrepancy -- noted for
completeness, in the same "disclose known gaps rather than let them stay
invisible" spirit as everything else in this project.

Conflicting-evidence handling also differs from Table 5, and this is a
genuine design difference, not an oversight: Table 5's combine() checks
pathogenic-direction and benign-direction combining rules independently
and can find both satisfied at once (conflicting_evidence_flag=True).
A single net point sum cannot land in two bands simultaneously -- when
evidence points both ways, the sum quantitatively resolves which
direction dominates rather than flagging an unresolved conflict. This
module's ProvisionalClassification.conflicting_evidence_flag is therefore
always False; the rationale text still states both sub-totals whenever
both directions contributed a MET criterion, so nothing is hidden, just
resolved differently -- one of the actual, published motivations for the
Bayesian framework's existence, not this project's own invention.
"""

from typing import List, Tuple

from .engine import evaluate_all
from .models import CriterionResult, ProvisionalClassification, VariantEvidenceBundle
from .models.enums import ClassificationStatus, CriterionStatus, CriterionStrength, EvidenceDirection, ProvisionalClass

RULE_SOURCE = (
    "Bayesian point-based combining (Tavtigian et al. 2020, Human Mutation 41(6):1023-1041), "
    "point values and thresholds as quoted in ACGS 2024 UK Practice Guidelines for Variant "
    "Classification v1"
)
RULE_VERSION = "2020"

PATHOGENIC_POINTS = {
    CriterionStrength.VERY_STRONG: 8,
    CriterionStrength.STRONG: 4,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.SUPPORTING: 1,
}

# No Very Strong tier on the benign side, matching both Tavtigian 2020 and
# Richards et al. 2015 -- no published ACMG/AMP criterion (this project's or
# otherwise) is defined at Very-Strong-Benign. MODERATE is defined here for
# completeness (Tavtigian's own table includes it) even though this
# project's BS1/BP4 evaluators never emit it -- BS1 is always Strong, BP4
# always Supporting.
BENIGN_POINTS = {
    CriterionStrength.STRONG: -4,
    CriterionStrength.MODERATE: -2,
    CriterionStrength.SUPPORTING: -1,
}


def _sum_points(met: List[CriterionResult]) -> Tuple[int, int, int, int]:
    """Returns (net_points, contributing_count, pathogenic_subtotal, benign_subtotal)."""
    pathogenic_subtotal = 0
    benign_subtotal = 0
    contributing = 0
    for c in met:
        if c.direction == EvidenceDirection.PATHOGENIC:
            pts = PATHOGENIC_POINTS.get(c.strength)
            if pts is not None:
                pathogenic_subtotal += pts
                contributing += 1
        else:
            pts = BENIGN_POINTS.get(c.strength)
            if pts is not None:
                benign_subtotal += pts
                contributing += 1
    return pathogenic_subtotal + benign_subtotal, contributing, pathogenic_subtotal, benign_subtotal


def _classify_points(points: int, contributing: int) -> Tuple[ProvisionalClass, str]:
    if contributing < 2:
        # Matches "variants with only one piece of evidence ... are classified as a
        # VUS pending a second corroborating piece of evidence" (ACGS 2024 / Tavtigian
        # et al. 2020) -- applies even if that one criterion's own points would
        # otherwise land in a Likely-Pathogenic/Pathogenic/Likely-Benign/Benign band.
        return ProvisionalClass.VUS, f"{points} points from only {contributing} contributing criterion/criteria (below the 2-criterion minimum)"
    if points >= 10:
        return ProvisionalClass.PATHOGENIC, f"{points} points (>=10)"
    if points >= 6:
        return ProvisionalClass.LIKELY_PATHOGENIC, f"{points} points (6-9)"
    if points >= 0:
        return ProvisionalClass.VUS, f"{points} points (0-5)"
    if points >= -5:
        return ProvisionalClass.LIKELY_BENIGN, f"{points} points (-1 to -5)"
    return ProvisionalClass.BENIGN, f"{points} points (<=-6)"


def combine_bayesian(criteria: List[CriterionResult]) -> ProvisionalClassification:
    met = [c for c in criteria if c.status == CriterionStatus.MET]

    manual_review_required = any(
        c.status in (CriterionStatus.MANUAL_REVIEW, CriterionStatus.CONFLICTING_EVIDENCE) for c in criteria
    )
    met_summary = ", ".join(f"{c.code} ({c.direction.value}/{c.strength.value})" for c in met) or "none"

    net_points, contributing, path_subtotal, benign_subtotal = _sum_points(met)

    ba1_stand_alone = any(
        c.code == "BA1" and c.status == CriterionStatus.MET and c.strength == CriterionStrength.STAND_ALONE
        for c in criteria
    )

    both_directions_contributed = path_subtotal > 0 and benign_subtotal < 0
    direction_note = (
        f" (pathogenic-direction subtotal +{path_subtotal}, benign-direction subtotal {benign_subtotal} -- "
        "both directions had MET criteria; the net sum resolves this quantitatively rather than flagging "
        "an unresolved conflict, unlike Table 5)"
        if both_directions_contributed
        else ""
    )

    if ba1_stand_alone:
        provisional_class = ProvisionalClass.BENIGN
        rationale = (
            f"BENIGN via BA1 Stand-Alone -- bypasses point-summing entirely, the one exception Tavtigian "
            f"et al. 2020 carries over unchanged from Table 5 (net point sum of all other MET criteria: "
            f"{net_points}, shown for information only, not used to decide this result). "
            f"MET criteria: {met_summary}."
        )
    else:
        provisional_class, band_note = _classify_points(net_points, contributing)
        rationale = (
            f"{provisional_class.value} via Bayesian point sum: {net_points} points from {contributing} "
            f"contributing criteria ({band_note}){direction_note}. MET criteria: {met_summary}."
        )

    return ProvisionalClassification(
        provisional_class=provisional_class,
        status=ClassificationStatus.PROVISIONAL_AUTOMATED,
        criteria=criteria,
        combining_rule_source=RULE_SOURCE,
        combining_rule_version=RULE_VERSION,
        rationale=rationale,
        conflicting_evidence_flag=False,
        manual_review_required=manual_review_required,
        points=net_points,
    )


def classify_bayesian(bundle: VariantEvidenceBundle, thresholds: dict) -> ProvisionalClassification:
    return combine_bayesian(evaluate_all(bundle, thresholds))
