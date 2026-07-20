"""Tests for loader.py against the real curated fixtures, plus a
malformed-record test proving loading rejects bad records safely instead
of crashing the whole batch or silently dropping the failure.
"""

import json
import tempfile
from pathlib import Path

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError


def test_load_gene_disease_contexts_from_real_fixture():
    contexts = loader.load_gene_disease_contexts()
    assert "CAPN3" in contexts
    assert contexts["CAPN3"].inheritance.value == "AUTOSOMAL_RECESSIVE"


def test_load_variant_evidence_bundles_from_real_fixture():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert len(bundles) == 3
    assert rejected == []
    ids = sorted(b.variant.variant_id for b in bundles)
    assert ids == ["CAPN3_SYNTH_LIKELY_BENIGN_01", "CAPN3_SYNTH_PATHOGENIC_01", "CAPN3_c.550del"]


def test_load_golden_cases_from_real_fixture():
    golden_cases = loader.load_golden_cases()
    assert len(golden_cases) == 3
    assert golden_cases["CAPN3_SYNTH_PATHOGENIC_01"].expected_provisional_class.value == "PATHOGENIC"


def test_load_all_real_fixtures_have_no_cross_check_warnings():
    result = loader.load_all()
    assert result["cross_check_warnings"] == [], result["cross_check_warnings"]
    assert len(result["evidence_bundles"]) == 3
    assert len(result["golden_cases"]) == 3
    assert result["rejected_evidence"] == []


def test_load_all_flags_golden_case_with_no_matching_bundle(tmp_path=None):
    # Build a minimal scratch repo layout so this test doesn't depend on
    # mutating the real curated fixtures.
    tmp_path = Path(tempfile.mkdtemp())
    (tmp_path / "data" / "curated").mkdir(parents=True)
    (tmp_path / "validation" / "golden_cases").mkdir(parents=True)

    real_gene_context = loader.CURATED_DIR / "gene_disease_context.yaml"
    (tmp_path / "data" / "curated" / "gene_disease_context.yaml").write_text(real_gene_context.read_text())

    real_evidence = loader.CURATED_DIR / "variant_evidence.json"
    evidence_data = json.loads(real_evidence.read_text())
    evidence_data["cases"] = evidence_data["cases"][:1]  # keep only CAPN3_c.550del
    (tmp_path / "data" / "curated" / "variant_evidence.json").write_text(json.dumps(evidence_data))

    golden_data = {
        "golden_cases": [
            {
                "variant_id": "NOT_A_REAL_VARIANT_ID",
                "expected_provisional_class": "VUS",
                "expected_criterion_status": {"PM2": "NOT_MET"},
                "source": "synthetic test fixture",
            }
        ]
    }
    (tmp_path / "validation" / "golden_cases" / "capn3_milestone1.yaml").write_text(
        json.dumps(golden_data)  # valid YAML is a superset concern here: JSON is valid YAML
    )

    result = loader.load_all(base=tmp_path)
    warnings = result["cross_check_warnings"]
    assert any("NOT_A_REAL_VARIANT_ID" in w and "no matching evidence bundle" in w for w in warnings)
    assert any("CAPN3_c.550del" in w and "no matching golden case" in w for w in warnings)


def test_load_variant_evidence_bundles_rejects_malformed_record_without_losing_good_ones():
    tmp_dir = Path(tempfile.mkdtemp())
    real_evidence = loader.CURATED_DIR / "variant_evidence.json"
    data = json.loads(real_evidence.read_text())
    good_case = data["cases"][0]
    broken_case = json.loads(json.dumps(good_case))  # deep copy
    del broken_case["population_evidence"]  # required field -> should be rejected, not silently coerced
    broken_case["variant"]["variant_id"] = "BROKEN_CASE"

    scratch_file = tmp_dir / "variant_evidence.json"
    scratch_file.write_text(json.dumps({"cases": [good_case, broken_case]}))

    bundles, rejected = loader.load_variant_evidence_bundles(scratch_file)
    assert len(bundles) == 1
    assert bundles[0].variant.variant_id == good_case["variant"]["variant_id"]
    assert len(rejected) == 1
    entry, message = rejected[0]
    assert entry["variant"]["variant_id"] == "BROKEN_CASE"
    assert "population_evidence" in message


def test_load_gene_disease_contexts_missing_file_raises_file_not_found():
    try:
        loader.load_gene_disease_contexts(Path("/nonexistent/path/gene_disease_context.yaml"))
    except FileNotFoundError:
        return
    raise AssertionError("expected FileNotFoundError")


def test_load_gene_disease_contexts_rejects_key_mismatch():
    tmp_dir = Path(tempfile.mkdtemp())
    scratch_file = tmp_dir / "gene_disease_context.yaml"
    scratch_file.write_text(
        "genes:\n"
        "  CAPN3:\n"
        "    gene: DMD\n"  # deliberate mismatch with the top-level key
        "    disease: x\n"
        "    inheritance: AUTOSOMAL_RECESSIVE\n"
        "    mechanism: LOSS_OF_FUNCTION\n"
        "    lof_established: true\n"
        "    specification:\n"
        "      type: GENERIC_ACMG\n"
        "      version: '2015'\n"
    )
    try:
        loader.load_gene_disease_contexts(scratch_file)
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError")
