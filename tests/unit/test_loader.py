"""Tests for loader.py against the real curated fixtures, plus a
malformed-record test proving loading rejects bad records safely instead
of crashing the whole batch or silently dropping the failure.
"""

import json
import tempfile
from pathlib import Path

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators.pm2 import PM2_INDEL_DELINS_CONSEQUENCES, evaluate_pm2
from variant_classifier.models.enums import CriterionStatus
from variant_classifier.models.transcript_consequence import TranscriptConsequence


def test_load_gene_disease_contexts_from_real_fixture():
    contexts = loader.load_gene_disease_contexts()
    assert "CAPN3" in contexts
    assert contexts["CAPN3"].inheritance.value == "AUTOSOMAL_RECESSIVE"
    assert "DMD" in contexts
    assert contexts["DMD"].inheritance.value == "X_LINKED_RECESSIVE"


def test_load_variant_evidence_bundles_from_real_fixture():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert len(bundles) == 35
    assert rejected == []
    ids = sorted(b.variant.variant_id for b in bundles)
    assert ids == [
        "BRCA1_SYNTH_BENIGN_01",
        "BRCA1_SYNTH_PATHOGENIC_01",
        "BRCA1_c.181T>G",
        "BRCA1_c.5266dup",
        "BRCA1_c.5559C>G",
        "BRCA1_c.68_69delAG",
        "CAPN3_SYNTH_LIKELY_BENIGN_01",
        "CAPN3_SYNTH_PATHOGENIC_01",
        "CAPN3_SYNTH_PATHOGENIC_02",
        "CAPN3_SYNTH_PM3_CIS_OVERRIDE_01",
        "CAPN3_SYNTH_PM3_MODERATE_01",
        "CAPN3_SYNTH_PS1_01",
        "CAPN3_SYNTH_PVS1_SPLICE_RNA_01",
        "CAPN3_SYNTH_PVS1_STARTLOSS_NO_ALT_01",
        "CAPN3_SYNTH_PVS1_STARTLOSS_SUPPORTING_01",
        "CAPN3_c.1132T>C",
        "CAPN3_c.1342C>T",
        "CAPN3_c.1343G>A",
        "CAPN3_c.1401_1403del",
        "CAPN3_c.1939G>T",
        "CAPN3_c.1A>G",
        "CAPN3_c.2050+1G>A",
        "CAPN3_c.2120A>G",
        "CAPN3_c.2257G>A",
        "CAPN3_c.550del",
        "CAPN3_c.598_612del",
        "CAPN3_c.946-1G>A",
        "DMD_SYNTH_PATHOGENIC_01",
        "DMD_c.10103A>G",
        "DMD_c.11041A>T",
        "DMD_c.2302C>T",
        "DMD_c.5163G>C",
        "DMD_c.5234G>A",
        "DMD_c.8944C>T",
        "DMD_c.93+1G>A",
    ]


def test_load_golden_cases_from_real_fixture():
    golden_cases = loader.load_golden_cases()
    assert len(golden_cases) == 35
    assert golden_cases["CAPN3_SYNTH_PATHOGENIC_01"].expected_provisional_class.value == "PATHOGENIC"


def test_load_all_real_fixtures_have_no_cross_check_warnings():
    result = loader.load_all()
    assert result["cross_check_warnings"] == [], result["cross_check_warnings"]
    assert len(result["evidence_bundles"]) == 35
    assert len(result["golden_cases"]) == 35
    assert result["rejected_evidence"] == []


def test_load_golden_cases_bayesian_from_real_fixture():
    # Added batch 20 (Milestone 5) alongside bayesian.py.
    goldens = loader.load_golden_cases_bayesian()
    assert len(goldens) == 35
    assert goldens["CAPN3_SYNTH_PATHOGENIC_01"]["expected_provisional_class"].value == "PATHOGENIC"
    assert goldens["CAPN3_SYNTH_PATHOGENIC_01"]["expected_points"] == 10
    # BA1 stand-alone fixtures bypass point-summing -- expected_points is null/None.
    assert goldens["DMD_c.5234G>A"]["expected_points"] is None
    # every variant with a Bayesian golden case must also have a real evidence bundle
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle_ids = {b.variant.variant_id for b in bundles}
    assert set(goldens) == bundle_ids


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
    (tmp_path / "validation" / "golden_cases" / "variant_golden_cases.yaml").write_text(
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


def test_brca1_fixtures_do_not_exercise_pm4_pm5_or_unenforced_pm2_indel_gaps():
    """Batch 31 guard test (design doc v2 section 9, refined during
    implementation -- see README.md's BRCA1 design note for the full
    writeup). BRCA1's real ENIGMA spec excludes PM4 (folded into PP3/BP4)
    and missense-based PM5 (its real PM5 is a structurally different,
    PTC/exon-only code) entirely, and neither evaluators/pm4.py nor
    evaluators/pm5.py has a gene gate -- this batch avoided curating a
    BRCA1 fixture shape that would exercise either gap rather than adding
    one to shared evaluator code CAPN3/DMD also run through (design doc
    sections 1.4/1.5). "Never silently guess" is this project's whole
    ethos; relying on that discipline alone for every future BRCA1
    fixture was flagged as insufficient, hence this test.

    PM2's indel/delins exclusion is different: it IS enforced in code
    (evaluators/pm2.py's pm2_excludes_indel_delins gate, added this batch
    after discovering fixture-shape avoidance cannot work for PM2 --
    population_evidence is a required field the founder fixtures also
    need, OBSERVED, for BA1/BS1's founder-frequency handling). This test
    still checks it explicitly rather than trusting the config alone: it
    confirms the gate is actually configured for BRCA1 and actually fires
    for every BRCA1 indel/delins fixture, so a future accidental removal
    of the config flag is caught here, not by chance.
    """
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    thresholds = loader.load_frequency_thresholds()
    brca1_bundles = [b for b in bundles if b.variant.gene == "BRCA1"]
    assert brca1_bundles, "expected at least one BRCA1 fixture -- this guard test would be vacuous otherwise"

    for bundle in brca1_bundles:
        variant_id = bundle.variant.variant_id
        transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)

        assert transcript.consequence not in TranscriptConsequence.PM4_RELEVANT_CONSEQUENCES, (
            f"{variant_id}: consequence {transcript.consequence.value} is PM4-relevant, but "
            "evaluators/pm4.py has no gene gate for BRCA1 (design doc section 1.4) -- this "
            "fixture would incorrectly exercise PM4 for a gene whose real spec folds it into "
            "PP3/BP4 instead. Curate a non-PM4-relevant consequence, or add a real PM4 gene gate "
            "before adding this fixture."
        )

        sre = bundle.same_residue_evidence
        assert sre is None or not sre.pm5_precedent_established, (
            f"{variant_id}: same_residue_evidence.pm5_precedent_established is set, but "
            "evaluators/pm5.py has no gene gate for BRCA1 (design doc section 1.5) -- ENIGMA's "
            "real BRCA1 PM5 is a structurally different, PTC/exon-only code, not the classic "
            "missense-residue PM5 this evaluator implements. Remove this precedent, or add a "
            "real PM5 gene gate before adding this fixture."
        )

        if transcript.consequence in PM2_INDEL_DELINS_CONSEQUENCES:
            gene_config = thresholds["genes"].get(bundle.variant.gene, {})
            assert gene_config.get("pm2_excludes_indel_delins", False), (
                f"{variant_id}: has an indel/delins consequence ({transcript.consequence.value}) "
                f"but {bundle.variant.gene}'s population_thresholds.yaml entry does not set "
                "pm2_excludes_indel_delins=True -- PM2 would incorrectly evaluate a real "
                "frequency comparison against a variant type the real BRCA1 spec excludes "
                "entirely (design doc section 6)."
            )
            result = evaluate_pm2(bundle, thresholds)
            assert result.status == CriterionStatus.NOT_APPLICABLE, (
                f"{variant_id}: pm2_excludes_indel_delins=True is configured, but evaluate_pm2 "
                f"returned {result.status.value}, not NOT_APPLICABLE -- the gate is not firing "
                "for this indel/delins consequence."
            )

