"""The Milestone 3 combining engine.

Two responsibilities, kept separate on purpose:

- evaluate_all(bundle, thresholds): run all six SUPPORTED_CRITERIA_MILESTONE_1
  evaluators against one VariantEvidenceBundle and return their six
  CriterionResults — always six, one per code, even when most are
  NOT_MET/NOT_EVALUATED/NOT_APPLICABLE. Nothing is dropped, so a report
  built from this list can show a scientist exactly what was and wasn't
  checked (the "preserved uncertainty" principle from the Workflow
  Architecture Guide).

- combine(criteria): apply the ACMG/AMP combining rules (Richards et al.
  2015, Table 5) to a list of CriterionResults and produce one
  ProvisionalClassification. This function does not know or care where
  the CriterionResults came from — it would work identically on results
  from criteria this project doesn't evaluate yet (PS1, PM3, ...), which
  is deliberate: extending criterion coverage later shouldn't require
  touching this function.

classify(bundle, thresholds) chains the two for convenience.

Two things this engine does NOT do, both by design for Milestone 3:
- It does not attempt PM3/PS4/other recessive-specific evidence (a second
  variant in trans, segregation, case-control counts) — that's Milestone 4,
  "clinical interpretation," which needs case-level information (parental
  testing, other observed alleles) this project's variant-level evidence
  bundle doesn't carry.
- Genuinely conflicting evidence (both a Pathogenic/Likely-Pathogenic
  combination AND a Benign/Likely-Benign combination satisfied at once)
  is reported as VUS with conflicting_evidence_flag=True, not silently
  resolved by picking a "winner" — ACMG/AMP treats this as its own
  category, not a tiebreak.
"""

from typing import List

from .evaluators import evaluate_ba1, evaluate_bp4, evaluate_bs1, evaluate_pm2, evaluate_pp3, evaluate_pvs1
from .models import CriterionResult, ProvisionalClassification, VariantEvidenceBundle
from .models.enums import ClassificationStatus, CriterionStatus, CriterionStrength, EvidenceDirection, ProvisionalClass

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015, Table 5)"
RULE_VERSION = "2015"


def evaluate_all(bundle: VariantEvidenceBundle, thresholds: dict) -> List[CriterionResult]:
    """Run all six Milestone-1-scope evaluators. Always returns exactly six
    CriterionResults, one per code in a fixed order — order doesn't matter
    for combine() but a fixed order makes output diffs/reports readable.
    """
    return [
        evaluate_pvs1(bundle),
        evaluate_pm2(bundle, thresholds),
        evaluate_pp3(bundle),
        evaluate_bp4(bundle),
        evaluate_ba1(bundle, thresholds),
        evaluate_bs1(bundle, thresholds),
    ]


def _pathogenic_tier(counts: dict):
    vs, s, m, sup = counts["VERY_STRONG"], counts["STRONG"], counts["MODERATE"], counts["SUPPORTING"]
    # Pathogenic (Table 5)
    if vs >= 1 and s >= 1:
        return ProvisionalClass.PATHOGENIC, "1 Very Strong + >=1 Strong"
    if vs >= 1 and m >= 2:
        return ProvisionalClass.PATHOGENIC, "1 Very Strong + >=2 Moderate"
    if vs >= 1 and m >= 1 and sup >= 1:
        return ProvisionalClass.PATHOGENIC, "1 Very Strong + 1 Moderate + >=1 Supporting"
    if vs >= 1 and sup >= 2:
        return ProvisionalClass.PATHOGENIC, "1 Very Strong + >=2 Supporting"
    if s >= 2:
        return ProvisionalClass.PATHOGENIC, ">=2 Strong"
    if s >= 1 and m >= 3:
        return ProvisionalClass.PATHOGENIC, "1 Strong + >=3 Moderate"
    if s >= 1 and m >= 2 and sup >= 2:
        return ProvisionalClass.PATHOGENIC, "1 Strong + 2 Moderate + >=2 Supporting"
    if s >= 1 and sup >= 4:
        return ProvisionalClass.PATHOGENIC, "1 Strong + >=4 Supporting"
    # Likely Pathogenic (Table 5)
    if vs >= 1 and m >= 1:
        return ProvisionalClass.LIKELY_PATHOGENIC, "1 Very Strong + 1 Moderate"
    if s >= 1 and 1 <= m <= 2:
        return ProvisionalClass.LIKELY_PATHOGENIC, "1 Strong + 1-2 Moderate"
    if s >= 1 and sup >= 2:
        return ProvisionalClass.LIKELY_PATHOGENIC, "1 Strong + >=2 Supporting"
    if m >= 3:
        return ProvisionalClass.LIKELY_PATHOGENIC, ">=3 Moderate"
    if m >= 2 and sup >= 2:
        return ProvisionalClass.LIKELY_PATHOGENIC, "2 Moderate + >=2 Supporting"
    if m >= 1 and sup >= 4:
        return ProvisionalClass.LIKELY_PATHOGENIC, "1 Moderate + >=4 Supporting"
    return None, None


def _benign_tier(counts: dict):
    sa, s, sup = counts["STAND_ALONE"], counts["STRONG"], counts["SUPPORTING"]
    # Benign (Table 5)
    if sa >= 1:
        return ProvisionalClass.BENIGN, "1 Stand-Alone (BA1)"
    if s >= 2:
        return ProvisionalClass.BENIGN, ">=2 Strong"
    # Likely Benign (Table 5)
    if s >= 1 and sup >= 1:
        return ProvisionalClass.LIKELY_BENIGN, "1 Strong + 1 Supporting"
    if sup >= 2:
        return ProvisionalClass.LIKELY_BENIGN, ">=2 Supporting"
    return None, None


def combine(criteria: List[CriterionResult]) -> ProvisionalClassification:
    met = [c for c in criteria if c.status == CriterionStatus.MET]

    pathogenic_counts = {"VERY_STRONG": 0, "STRONG": 0, "MODERATE": 0, "SUPPORTING": 0}
    benign_counts = {"STAND_ALONE": 0, "STRONG": 0, "SUPPORTING": 0}

    for c in met:
        if c.direction == EvidenceDirection.PATHOGENIC:
            pathogenic_counts[c.strength.value] += 1
        else:
            benign_counts[c.strength.value] += 1

    path_class, path_rule = _pathogenic_tier(pathogenic_counts)
    benign_class, benign_rule = _benign_tier(benign_counts)

    manual_review_required = any(
        c.status in (CriterionStatus.MANUAL_REVIEW, CriterionStatus.CONFLICTING_EVIDENCE) for c in criteria
    )

    met_summary = ", ".join(f"{c.code} ({c.direction.value}/{c.strength.value})" for c in met) or "none"

    if path_class is not None and benign_class is not None:
        provisional_class = ProvisionalClass.VUS
        conflicting_evidence_flag = True
        rationale = (
            f"Conflicting evidence: pathogenic-side criteria satisfy {path_class.value} "
            f"({path_rule}) while benign-side criteria simultaneously satisfy {benign_class.value} "
            f"({benign_rule}). Per Table 5 this is reported as Uncertain Significance with a conflict "
            f"flag, not resolved by picking a side. MET criteria: {met_summary}."
        )
    elif path_class is not None:
        provisional_class = path_class
        conflicting_evidence_flag = False
        rationale = f"{path_class.value} via combining rule: {path_rule}. MET criteria: {met_summary}."
    elif benign_class is not None:
        provisional_class = benign_class
        conflicting_evidence_flag = False
        rationale = f"{benign_class.value} via combining rule: {benign_rule}. MET criteria: {met_summary}."
    else:
        provisional_class = ProvisionalClass.VUS
        conflicting_evidence_flag = False
        rationale = (
            "No combining rule satisfied on either side (insufficient evidence, not conflicting "
            f"evidence). MET criteria: {met_summary}."
        )

    return ProvisionalClassification(
        provisional_class=provisional_class,
        status=ClassificationStatus.PROVISIONAL_AUTOMATED,
        criteria=criteria,
        combining_rule_source=RULE_SOURCE,
        combining_rule_version=RULE_VERSION,
        rationale=rationale,
        conflicting_evidence_flag=conflicting_evidence_flag,
        manual_review_required=manual_review_required,
    )


def classify(bundle: VariantEvidenceBundle, thresholds: dict) -> ProvisionalClassification:
    return combine(evaluate_all(bundle, thresholds))
