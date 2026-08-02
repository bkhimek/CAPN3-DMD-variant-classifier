"""PM3 evaluator -- "for recessive disorders, detected in trans with a
pathogenic variant" (Richards et al. 2015, Table 3). Added batch 28,
closing the case-level gap this project disclosed since Milestone 4 (see
clinical.py's and engine.py's module docstrings).

See models/pm3_evidence.py for the full design writeup: why the partner
allele's classification is a curated fact rather than something this
evaluator re-derives (avoids the circularity a live, per-variant PM3
evaluator would otherwise have), which real ACGS 2024 rules are enforced
directly (the homozygous 1-point cap, the cis-cooccurrence override), and
why the exact per-scenario points table is deliberately NOT hardcoded --
that table's own primary source (ClinGen SVI "Recommendation for the in
trans Criterion (PM3)" Version 1.0) was unreachable via web_fetch on every
attempt this session, so `points` on each proband observation is itself a
curated fact rather than a value this project asserts it has verified.

Threshold bands (Very Strong >=4, Strong >=2 but <4, Moderate >=1 but <2,
Supporting >=0.5 but <1) are the CAPN3 LGMD VCEP's own real, gene-specific
PM3 points table (cspec.genome.network/cspec/ui/svi/doc/GN187) -- the same
real source already used for this project's CAPN3-specific PM2/BA1/BS1
thresholds (config/population_thresholds.yaml). Scope, disclosed rather
than silently assumed: this evaluator applies CAPN3's own confirmed
thresholds project-wide (including to DMD fixtures), since no
DMD-specific PM3 points table was located this session and DMD is
predominantly X-linked (PM3 does not apply to the X-linked fixtures in
this project's curated set at all, since PM3 is defined for recessive
disorders) -- same "one gene's real VCEP numbers, applied as the best
available real-world anchor rather than a project-invented placeholder"
pattern already used for CAPN3's PM2/BA1/BS1 thresholds.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import CriterionStatus, CriterionStrength, EvidenceDirection

RULE_SOURCE = (
    "ACMG/AMP (Richards et al. 2015); points-based combining per ACGS 2024 UK Practice Guidelines "
    "for Variant Classification v1 (citing ClinGen SVI 'Recommendation for the in trans Criterion "
    "(PM3)' Version 1.0); strength thresholds per ClinGen LGMD VCEP CAPN3 specification v2.0"
)
RULE_VERSION = "2015 / 2024"

_VERY_STRONG_THRESHOLD = 4.0
_STRONG_THRESHOLD = 2.0
_MODERATE_THRESHOLD = 1.0
_SUPPORTING_THRESHOLD = 0.5


def evaluate_pm3(bundle: VariantEvidenceBundle) -> CriterionResult:
    pm3 = bundle.pm3_evidence
    evidence_id = f"pm3_evidence:{bundle.variant.variant_id}"

    if pm3 is None:
        return CriterionResult(
            code="PM3",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="No pm3_evidence recorded for this variant -- PM3 was never assessed.",
            evidence_ids=[evidence_id],
        )

    if pm3.cis_cooccurrence_observed:
        return CriterionResult(
            code="PM3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                "Population-level cis co-occurrence has been observed for this variant "
                f"({pm3.notes or 'see curated notes'}) -- per ACGS 2024, 'PM3 should not be applied "
                "at any level in the context of two variants that predominantly co-occur.' This "
                "overrides any proband-level trans observations that may also be recorded."
            ),
            evidence_ids=[evidence_id],
        )

    total_points = sum(p.points for p in pm3.probands)
    proband_summary = ", ".join(
        f"{p.proband_id} ({p.zygosity.value}, {p.points:g} pt)" for p in pm3.probands
    ) or "none"

    if total_points >= _VERY_STRONG_THRESHOLD:
        strength = CriterionStrength.VERY_STRONG
    elif total_points >= _STRONG_THRESHOLD:
        strength = CriterionStrength.STRONG
    elif total_points >= _MODERATE_THRESHOLD:
        strength = CriterionStrength.MODERATE
    elif total_points >= _SUPPORTING_THRESHOLD:
        strength = CriterionStrength.SUPPORTING
    else:
        return CriterionResult(
            code="PM3",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Summed proband points ({total_points:g}) fall below the Supporting threshold "
                f"({_SUPPORTING_THRESHOLD:g}) -- insufficient in-trans evidence to meet PM3 at any "
                f"strength. Probands: {proband_summary}."
            ),
            evidence_ids=[evidence_id],
        )

    return CriterionResult(
        code="PM3",
        status=CriterionStatus.MET,
        strength=strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Summed proband points ({total_points:g}) meet the {strength.value} threshold for PM3 "
            f"(CAPN3 LGMD VCEP v2.0 points table). Probands: {proband_summary}."
        ),
        evidence_ids=[evidence_id],
    )
