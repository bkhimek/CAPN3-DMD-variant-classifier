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
from .models import GeneDiseaseContext, GoldenCase, VariantEvidenceBundle

REPO_ROOT = Path(__file__).resolve().parents[2]
CURATED_DIR = REPO_ROOT / "data" / "curated"
GOLDEN_CASES_DIR = REPO_ROOT / "validation" / "golden_cases"


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
    """Load validation/golden_cases/capn3_milestone1.yaml into
    {variant_id: GoldenCase}. Golden cases are curated separately from
    data/curated/ on purpose — see golden_case.py docstring.
    """
    path = path or (GOLDEN_CASES_DIR / "capn3_milestone1.yaml")
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
    golden_cases = load_golden_cases((base / "validation" / "golden_cases" / "capn3_milestone1.yaml") if base else None)

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
