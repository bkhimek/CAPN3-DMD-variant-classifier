"""Fixture loading for Milestone 1.

Deliberately dumb: this module's only job is to turn the curated files
under data/curated/ into validated model instances, rejecting anything
malformed with a SchemaValidationError rather than silently coercing or
skipping it. No evaluators, no combining logic — that is Milestone 2/3
scope (see README.md).
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple

import yaml

from .errors import SchemaValidationError
from .models import CnvDeletionEvidence, ClinicalCase, GeneDiseaseContext, GoldenCase, VariantEvidenceBundle
from .models.enums import CaseInterpretationStatus, ProvisionalClass

REPO_ROOT = Path(__file__).resolve().parents[2]
CURATED_DIR = REPO_ROOT / "data" / "curated"
GOLDEN_CASES_DIR = REPO_ROOT / "validation" / "golden_cases"
CONFIG_DIR = REPO_ROOT / "config"


def load_gene_disease_contexts(path: Path = None) -> Dict[str, GeneDiseaseContext]:
    """Load data/curated/gene_disease_context.yaml into {gene: GeneDiseaseContext}.

    Raises SchemaValidationError on the first malformed entry — this file is
    small and hand-curated, so fail loudly rather than silently drop a gene.
    """
    path = path or (CURATED_DIR / "gene_disease_context.yaml")
    if not path.exists():
        raise FileNotFoundError(f"gene/disease context file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "genes" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'genes' mapping")
    contexts: Dict[str, GeneDiseaseContext] = {}
    for gene, entry in raw["genes"].items():
        context = GeneDiseaseContext.from_dict(entry)
        if context.gene != gene:
            raise SchemaValidationError(
                f"{path}: top-level key {gene!r} does not match embedded gene field {context.gene!r}"
            )
        contexts[gene] = context
    return contexts


_VALID_PM2_STRENGTHS = ("STAND_ALONE", "VERY_STRONG", "STRONG", "MODERATE", "SUPPORTING")


def load_frequency_thresholds(path: Path = None) -> dict:
    """Load config/population_thresholds.yaml into:

        {
            "ba1_stand_alone_af": float,
            "genes": {
                gene: {
                    "pm2_max_credible_af": float,
                    "pm2_strength": str,   # one of CriterionStrength's names; defaults to "MODERATE"
                    "bs1_min_af": float,
                    "ba1_af": float,       # optional per-gene BA1 override; absent means "use the global default"
                    "threshold_source": str,
                },
            },
        }

    Used by evaluators.pm2, evaluators.ba1, and evaluators.bs1 — kept in
    loader.py alongside the other "read a curated file, validate its
    shape" functions rather than inside the evaluators package, so
    evaluators stay focused on decision logic rather than file I/O.

    Named load_frequency_thresholds (not load_pm2_thresholds, its original
    Milestone-2 name) since it now serves three evaluators, not one.

    pm2_strength and ba1_af were added during the batch-4 curated-set
    expansion, when CAPN3 moved from a placeholder threshold to the real
    ClinGen LGMD VCEP specification (which uses PM2_Supporting, not
    PM2_Moderate, and a gene-specific BA1 threshold of 0.003 rather than
    the generic 0.05). Both are optional per gene, defaulting to the
    generic-ACMG conventions (Moderate strength, the global BA1 default)
    for any gene without a VCEP-sourced override — this is what keeps DMD
    (no VCEP spec adopted here) behaving exactly as before.
    """
    path = path or (CONFIG_DIR / "population_thresholds.yaml")
    if not path.exists():
        raise FileNotFoundError(f"Frequency threshold config not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "genes" not in raw or "ba1_stand_alone_af" not in raw:
        raise SchemaValidationError(f"{path}: expected top-level 'ba1_stand_alone_af' and 'genes' keys")

    ba1_af = raw["ba1_stand_alone_af"]
    if not isinstance(ba1_af, (int, float)) or isinstance(ba1_af, bool) or not (0.0 <= ba1_af <= 1.0):
        raise SchemaValidationError(f"{path}: ba1_stand_alone_af must be a number in [0, 1]")

    def _validate_af(key, af):
        if not isinstance(af, (int, float)) or isinstance(af, bool) or not (0.0 <= af <= 1.0):
            raise SchemaValidationError(f"{path}: genes.{gene}.{key} must be a number in [0, 1]")

    genes: Dict[str, dict] = {}
    for gene, entry in raw["genes"].items():
        if not isinstance(entry, dict) or "pm2_max_credible_af" not in entry or "bs1_min_af" not in entry:
            raise SchemaValidationError(
                f"{path}: genes.{gene} must include both 'pm2_max_credible_af' and 'bs1_min_af'"
            )
        pm2_af = entry["pm2_max_credible_af"]
        bs1_af = entry["bs1_min_af"]
        _validate_af("pm2_max_credible_af", pm2_af)
        _validate_af("bs1_min_af", bs1_af)

        pm2_strength = entry.get("pm2_strength", "MODERATE")
        if pm2_strength not in _VALID_PM2_STRENGTHS:
            raise SchemaValidationError(
                f"{path}: genes.{gene}.pm2_strength {pm2_strength!r} is not one of {_VALID_PM2_STRENGTHS}"
            )

        gene_config = {
            "pm2_max_credible_af": float(pm2_af),
            "pm2_strength": pm2_strength,
            "bs1_min_af": float(bs1_af),
            "threshold_source": entry.get("threshold_source", ""),
        }
        if "ba1_af" in entry:
            _validate_af("ba1_af", entry["ba1_af"])
            gene_config["ba1_af"] = float(entry["ba1_af"])
        genes[gene] = gene_config

    return {"ba1_stand_alone_af": float(ba1_af), "genes": genes}


def load_dosage_sensitivity(path: Path = None) -> Dict[str, dict]:
    """Load config/dosage_sensitivity.yaml into:

        {gene: {"hi_score": int, "hi_established": bool, "source": str}}

    Used by cnv_scoring.py's category 2A check ("complete overlap of an
    established dosage-sensitive gene"). Added batch 23 alongside
    cnv_scoring.py and models/cnv_deletion_evidence.py -- kept in loader.py
    next to load_frequency_thresholds() since both are the same shape of
    thing: a small, hand-curated, per-gene YAML config an evaluator/scoring
    module needs before it can reason about a specific variant/CNV.

    hi_established must equal (hi_score == 3) -- ClinGen's own "sufficient
    evidence" bar for category 2A -- checked here rather than trusted
    blindly from the YAML, so a future curator typo (hi_score: 2,
    hi_established: true) is rejected loudly instead of silently
    misclassifying a gene as dosage-sensitive-established.
    """
    path = path or (CONFIG_DIR / "dosage_sensitivity.yaml")
    if not path.exists():
        raise FileNotFoundError(f"Dosage sensitivity config not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "genes" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'genes' mapping")

    genes: Dict[str, dict] = {}
    for gene, entry in raw["genes"].items():
        if not isinstance(entry, dict) or "hi_score" not in entry or "hi_established" not in entry:
            raise SchemaValidationError(f"{path}: genes.{gene} must include both 'hi_score' and 'hi_established'")
        hi_score = entry["hi_score"]
        hi_established = entry["hi_established"]
        if not isinstance(hi_score, int) or isinstance(hi_score, bool) or not (0 <= hi_score <= 3):
            raise SchemaValidationError(f"{path}: genes.{gene}.hi_score must be an integer in [0, 3]")
        if not isinstance(hi_established, bool):
            raise SchemaValidationError(f"{path}: genes.{gene}.hi_established must be true/false")
        if hi_established != (hi_score == 3):
            raise SchemaValidationError(
                f"{path}: genes.{gene}.hi_established={hi_established} is inconsistent with hi_score="
                f"{hi_score} -- hi_established must be true iff hi_score == 3"
            )
        genes[gene] = {
            "hi_score": hi_score,
            "hi_established": hi_established,
            "source": entry.get("source", ""),
        }
    return genes


def load_variant_evidence_bundles(path: Path = None) -> Tuple[List[VariantEvidenceBundle], List[Tuple[dict, str]]]:
    """Load data/curated/variant_evidence.json into a list of validated
    VariantEvidenceBundle instances.

    Returns (bundles, rejected) rather than raising on the first bad record:
    a batch fixture file can have one bad case among several good ones, and
    the point of Milestone 1's validation layer is to *reject malformed
    records safely*, not to abort the whole load. Each rejected entry is
    (raw_dict, error_message).
    """
    path = path or (CURATED_DIR / "variant_evidence.json")
    if not path.exists():
        raise FileNotFoundError(f"variant evidence file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict) or "cases" not in raw or not isinstance(raw["cases"], list):
        raise SchemaValidationError(f"{path}: expected a top-level 'cases' list")

    bundles: List[VariantEvidenceBundle] = []
    rejected: List[Tuple[dict, str]] = []
    for i, entry in enumerate(raw["cases"]):
        try:
            bundles.append(VariantEvidenceBundle.from_dict(entry, f"cases[{i}]"))
        except SchemaValidationError as exc:
            rejected.append((entry, str(exc)))
    return bundles, rejected


def load_golden_cases(path: Path = None) -> Dict[str, GoldenCase]:
    """Load validation/golden_cases/variant_golden_cases.yaml into
    {variant_id: GoldenCase}. Golden cases are curated separately from
    data/curated/ on purpose — see golden_case.py docstring.
    """
    path = path or (GOLDEN_CASES_DIR / "variant_golden_cases.yaml")
    if not path.exists():
        raise FileNotFoundError(f"golden case file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "golden_cases" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'golden_cases' list")
    cases_raw = raw["golden_cases"]
    if not isinstance(cases_raw, list):
        raise SchemaValidationError(f"{path}: 'golden_cases' must be a list")
    golden_cases: Dict[str, GoldenCase] = {}
    for i, entry in enumerate(cases_raw):
        case = GoldenCase.from_dict(entry, f"golden_cases[{i}]")
        if case.variant_id in golden_cases:
            raise SchemaValidationError(f"{path}: duplicate golden case for variant_id={case.variant_id!r}")
        golden_cases[case.variant_id] = case
    return golden_cases


def load_golden_cases_bayesian(path: Path = None) -> Dict[str, dict]:
    """Load validation/golden_cases/variant_golden_cases_bayesian.yaml into
    {variant_id: {"expected_provisional_class": ProvisionalClass, "expected_points":
    Optional[int], "source": str, "curator_note": str}}.

    Deliberately a plain dict, not a GoldenCase (added batch 20 / Milestone 5,
    alongside bayesian.py): GoldenCase requires a non-empty
    expected_criterion_status mapping, which doesn't apply here -- the
    per-criterion results are identical to the Table 5 file (same
    evaluate_all() output feeds both combining systems), so this file only
    records what's new: the point total and resulting class. Same
    lightweight-dict pattern as load_case_interpretation_goldens().
    """
    path = path or (GOLDEN_CASES_DIR / "variant_golden_cases_bayesian.yaml")
    if not path.exists():
        raise FileNotFoundError(f"Bayesian golden case file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "golden_cases" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'golden_cases' list")
    entries = raw["golden_cases"]
    if not isinstance(entries, list):
        raise SchemaValidationError(f"{path}: 'golden_cases' must be a list")
    goldens: Dict[str, dict] = {}
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or "variant_id" not in entry or "expected_provisional_class" not in entry:
            raise SchemaValidationError(
                f"{path}: golden_cases[{i}] must include 'variant_id' and 'expected_provisional_class'"
            )
        variant_id = entry["variant_id"]
        try:
            expected_provisional_class = ProvisionalClass(entry["expected_provisional_class"])
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in ProvisionalClass))
            raise SchemaValidationError(
                f"{path}: golden_cases[{i}].expected_provisional_class="
                f"{entry['expected_provisional_class']!r} invalid; expected one of {valid}"
            ) from exc
        if variant_id in goldens:
            raise SchemaValidationError(f"{path}: duplicate variant_id {variant_id!r}")
        goldens[variant_id] = {
            "expected_provisional_class": expected_provisional_class,
            "expected_points": entry.get("expected_points"),
            "source": entry.get("source", ""),
            "curator_note": entry.get("curator_note", ""),
        }
    return goldens


def load_cnv_deletion_evidence(path: Path = None) -> Tuple[List[CnvDeletionEvidence], List[Tuple[dict, str]]]:
    """Load data/curated/cnv_deletion_evidence.json into a list of validated
    CnvDeletionEvidence instances. Mirrors load_variant_evidence_bundles():
    returns (evidence, rejected) rather than raising on the first bad
    record, for the same reason -- a growing curated set can have one bad
    entry among several good ones.
    """
    path = path or (CURATED_DIR / "cnv_deletion_evidence.json")
    if not path.exists():
        raise FileNotFoundError(f"CNV deletion evidence file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict) or "cases" not in raw or not isinstance(raw["cases"], list):
        raise SchemaValidationError(f"{path}: expected a top-level 'cases' list")

    evidence: List[CnvDeletionEvidence] = []
    rejected: List[Tuple[dict, str]] = []
    for i, entry in enumerate(raw["cases"]):
        try:
            evidence.append(CnvDeletionEvidence.from_dict(entry, f"cases[{i}]"))
        except SchemaValidationError as exc:
            rejected.append((entry, str(exc)))
    return evidence, rejected


def load_cnv_deletion_golden_cases(path: Path = None) -> Dict[str, dict]:
    """Load validation/golden_cases/cnv_deletion_golden_cases.yaml into
    {cnv_id: {"expected_provisional_class": ProvisionalClass, "expected_points":
    Optional[float], "expected_category_code": Optional[str], "source": str,
    "curator_note": str}}.

    Same lightweight-dict pattern as load_golden_cases_bayesian() -- a full
    GoldenCase doesn't fit here (it requires a non-empty
    expected_criterion_status mapping keyed by ACMG codes, which a CNV
    result never has), so this records only what cnv_scoring.py actually
    produces: a provisional class, a net point total, and (for this
    project's single-category-per-CNV model) which category code drove it.
    """
    path = path or (GOLDEN_CASES_DIR / "cnv_deletion_golden_cases.yaml")
    if not path.exists():
        raise FileNotFoundError(f"CNV deletion golden case file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "golden_cases" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'golden_cases' list")
    entries = raw["golden_cases"]
    if not isinstance(entries, list):
        raise SchemaValidationError(f"{path}: 'golden_cases' must be a list")
    goldens: Dict[str, dict] = {}
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or "cnv_id" not in entry or "expected_provisional_class" not in entry:
            raise SchemaValidationError(
                f"{path}: golden_cases[{i}] must include 'cnv_id' and 'expected_provisional_class'"
            )
        cnv_id = entry["cnv_id"]
        try:
            expected_provisional_class = ProvisionalClass(entry["expected_provisional_class"])
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in ProvisionalClass))
            raise SchemaValidationError(
                f"{path}: golden_cases[{i}].expected_provisional_class="
                f"{entry['expected_provisional_class']!r} invalid; expected one of {valid}"
            ) from exc
        if cnv_id in goldens:
            raise SchemaValidationError(f"{path}: duplicate cnv_id {cnv_id!r}")
        expected_points = entry.get("expected_points")
        goldens[cnv_id] = {
            "expected_provisional_class": expected_provisional_class,
            "expected_points": float(expected_points) if expected_points is not None else None,
            "expected_category_code": entry.get("expected_category_code"),
            "source": entry.get("source", ""),
            "curator_note": entry.get("curator_note", ""),
        }
    return goldens


def load_clinical_cases(path: Path = None) -> List[ClinicalCase]:
    """Load data/curated/clinical_cases.json into a list of validated
    ClinicalCase instances. Unlike load_variant_evidence_bundles, this
    raises on the first malformed entry rather than collecting rejects —
    this file is small and hand-curated (Milestone 4 scope: at most a
    handful of cases), so failing loudly is more useful than partial
    tolerance here.
    """
    path = path or (CURATED_DIR / "clinical_cases.json")
    if not path.exists():
        raise FileNotFoundError(f"clinical cases file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict) or "cases" not in raw or not isinstance(raw["cases"], list):
        raise SchemaValidationError(f"{path}: expected a top-level 'cases' list")
    cases = [ClinicalCase.from_dict(entry, f"cases[{i}]") for i, entry in enumerate(raw["cases"])]
    seen_ids = set()
    for case in cases:
        if case.case_id in seen_ids:
            raise SchemaValidationError(f"{path}: duplicate case_id {case.case_id!r}")
        seen_ids.add(case.case_id)
    return cases


def load_case_interpretation_goldens(path: Path = None) -> Dict[str, dict]:
    """Load validation/golden_cases/case_interpretation_golden_cases.yaml
    into {case_id: {"expected_status": CaseInterpretationStatus, "source": str, "curator_note": str}}.
    """
    path = path or (GOLDEN_CASES_DIR / "case_interpretation_golden_cases.yaml")
    if not path.exists():
        raise FileNotFoundError(f"case interpretation golden case file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict) or "case_interpretations" not in raw:
        raise SchemaValidationError(f"{path}: expected a top-level 'case_interpretations' list")
    entries = raw["case_interpretations"]
    if not isinstance(entries, list):
        raise SchemaValidationError(f"{path}: 'case_interpretations' must be a list")
    goldens: Dict[str, dict] = {}
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or "case_id" not in entry or "expected_status" not in entry:
            raise SchemaValidationError(f"{path}: case_interpretations[{i}] must include 'case_id' and 'expected_status'")
        case_id = entry["case_id"]
        try:
            expected_status = CaseInterpretationStatus(entry["expected_status"])
        except ValueError as exc:
            valid = ", ".join(sorted(v.value for v in CaseInterpretationStatus))
            raise SchemaValidationError(
                f"{path}: case_interpretations[{i}].expected_status={entry['expected_status']!r} invalid; expected one of {valid}"
            ) from exc
        if case_id in goldens:
            raise SchemaValidationError(f"{path}: duplicate case_id {case_id!r}")
        goldens[case_id] = {
            "expected_status": expected_status,
            "source": entry.get("source", ""),
            "curator_note": entry.get("curator_note", ""),
        }
    return goldens


def load_all(base: Path = None) -> dict:
    """Convenience entry point: load contexts, evidence bundles, and golden
    cases together, and cross-check that every bundle has a gene/disease
    context and every golden case matches a loaded bundle. Returns a dict
    with keys: gene_disease_contexts, evidence_bundles, rejected_evidence,
    golden_cases, cross_check_warnings.
    """
    base = base or REPO_ROOT
    contexts = load_gene_disease_contexts((base / "data" / "curated" / "gene_disease_context.yaml") if base else None)
    bundles, rejected = load_variant_evidence_bundles((base / "data" / "curated" / "variant_evidence.json") if base else None)
    golden_cases = load_golden_cases((base / "validation" / "golden_cases" / "variant_golden_cases.yaml") if base else None)

    warnings: List[str] = []
    bundle_ids = {b.variant.variant_id for b in bundles}
    for variant_id in golden_cases:
        if variant_id not in bundle_ids:
            warnings.append(f"golden case for {variant_id!r} has no matching evidence bundle")
    for bundle in bundles:
        if bundle.variant.variant_id not in golden_cases:
            warnings.append(f"evidence bundle {bundle.variant.variant_id!r} has no matching golden case")
        if bundle.gene_disease_context.gene not in contexts:
            warnings.append(
                f"evidence bundle {bundle.variant.variant_id!r} references gene "
                f"{bundle.gene_disease_context.gene!r} not present in gene_disease_context.yaml"
            )

    return {
        "gene_disease_contexts": contexts,
        "evidence_bundles": bundles,
        "rejected_evidence": rejected,
        "golden_cases": golden_cases,
        "cross_check_warnings": warnings,
    }


if __name__ == "__main__":
    result = load_all()
    print(f"Loaded {len(result['gene_disease_contexts'])} gene/disease context(s): "
          f"{sorted(result['gene_disease_contexts'])}")
    print(f"Loaded {len(result['evidence_bundles'])} evidence bundle(s), "
          f"{len(result['rejected_evidence'])} rejected")
    for entry, err in result["rejected_evidence"]:
        print(f"  REJECTED: {err}")
    print(f"Loaded {len(result['golden_cases'])} golden case(s): {sorted(result['golden_cases'])}")
    for w in result["cross_check_warnings"]:
        print(f"  WARNING: {w}")
