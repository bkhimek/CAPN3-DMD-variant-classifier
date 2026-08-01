"""Tests for batch 23's DMD CNV/deletion scoring AND batch 24's DMD CNV/
duplication scoring -- cnv_scoring.py, models/cnv_deletion_evidence.py,
models/cnv_duplication_evidence.py, models/cnv_category_result.py,
models/cnv_provisional_classification.py, and the loader functions for
both curated CNV sets.

Same two-part structure as every other evaluator's tests in this project:
(1) run against every curated CNV fixture and check it matches its golden
case; (2) hand-built edge cases for every decision branch and validation
rule, since the curated sets (small and real/literature-grounded on
purpose, see cnv_deletion_evidence.py / cnv_duplication_evidence.py) do
not happen to exercise every branch of either decision tree (2C, both 2D
point values, and 2F have no real/literature deletion fixture; the
not-tandem and whole-gene-duplication branches have no real/literature
duplication fixture -- see README.md).
"""

from variant_classifier import cnv_scoring, loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.models import CnvDeletionEvidence, CnvDuplicationEvidence
from variant_classifier.models.enums import (
    CnvDuplicationOrientation,
    CnvReadingFrameEffect,
    CriterionStatus,
    GenomeBuild,
    ProvisionalClass,
)


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


# ------------------------------------------------------- against real fixtures

def test_cnv_scoring_matches_golden_case_for_all_curated_evidence():
    evidence, rejected = loader.load_cnv_deletion_evidence()
    assert rejected == []
    dosage_config = loader.load_dosage_sensitivity()
    golden_cases = loader.load_cnv_deletion_golden_cases()

    assert len(evidence) == 3
    checked = 0
    for ev in evidence:
        golden = golden_cases[ev.cnv_id]
        result = cnv_scoring.score_cnv_deletion(ev, dosage_config)
        assert result.provisional_class == golden["expected_provisional_class"], (
            f"{ev.cnv_id}: scored {result.provisional_class}, golden case expects "
            f"{golden['expected_provisional_class']}. Rationale: {result.rationale}"
        )
        if golden["expected_points"] is not None:
            assert result.points == golden["expected_points"], (
                f"{ev.cnv_id}: scored {result.points} points, golden case expects {golden['expected_points']}"
            )
        if golden["expected_category_code"] is not None:
            assert result.categories[0].code == golden["expected_category_code"], (
                f"{ev.cnv_id}: category {result.categories[0].code}, golden case expects "
                f"{golden['expected_category_code']}"
            )
        checked += 1
    assert checked == 3


def test_cnv_deletion_golden_cases_cover_every_curated_cnv_id():
    evidence, _ = loader.load_cnv_deletion_evidence()
    golden_cases = loader.load_cnv_deletion_golden_cases()
    evidence_ids = {ev.cnv_id for ev in evidence}
    assert evidence_ids == set(golden_cases.keys())


# ---------------------------------------------------------- dosage sensitivity config

def test_dosage_sensitivity_config_has_established_dmd():
    dosage_config = loader.load_dosage_sensitivity()
    assert dosage_config["DMD"]["hi_score"] == 3
    assert dosage_config["DMD"]["hi_established"] is True
    assert "CAPN3" not in dosage_config  # autosomal recessive -- see config file's own docstring


# ---------------------------------------------------------- hand-built decision-tree edge cases

def _dmd_dosage_config():
    return loader.load_dosage_sensitivity()


def test_whole_gene_deletion_of_established_gene_is_2a_pathogenic():
    ev = CnvDeletionEvidence(cnv_id="EDGE_2A", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=True)
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2A"
    assert result.categories[0].status == CriterionStatus.MET
    assert result.points == 1.0
    assert result.provisional_class == ProvisionalClass.PATHOGENIC


def test_whole_gene_deletion_of_unestablished_gene_is_none_applicable():
    ev = CnvDeletionEvidence(cnv_id="EDGE_2A_NONE", gene="UNLISTED_GENE", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=True)
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.points == 0.0
    assert result.provisional_class == ProvisionalClass.VUS


def test_benign_region_overlap_is_2f_benign():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2F", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        overlaps_benign_region=True,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2F"
    assert result.points == -1.0
    assert result.provisional_class == ProvisionalClass.BENIGN


def test_intragenic_out_of_frame_with_nmd_is_2e_likely_pathogenic():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2E", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        reading_frame_effect=CnvReadingFrameEffect.OUT_OF_FRAME, nmd_predicted=True,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2E"
    assert result.points == 0.9
    assert result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC


def test_intragenic_out_of_frame_but_nmd_escaped_is_none_applicable():
    # An out-of-frame deletion confined to the last exon can escape NMD --
    # the same real caveat TranscriptConsequence.nmd_predicted encodes for
    # point-variant frameshifts. 2E requires BOTH out-of-frame AND NMD.
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2E_NO_NMD", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        reading_frame_effect=CnvReadingFrameEffect.OUT_OF_FRAME, nmd_predicted=False,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_intragenic_in_frame_is_none_applicable_vus():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_INFRAME", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        reading_frame_effect=CnvReadingFrameEffect.IN_FRAME,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_five_prime_end_plus_cds_is_2c_likely_pathogenic():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2C", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        five_prime_end_deleted=True, cds_involved=True,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2C"
    assert result.points == 0.9
    assert result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC


def test_five_prime_end_without_cds_is_none_applicable():
    # UTR-only 5' end deletion, no CDS involvement -- real Riggs 2020 2C
    # definitionally requires CDS involvement, so this is a legitimate
    # zero-points outcome, not an error.
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2C_UTR_ONLY", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        five_prime_end_deleted=True, cds_involved=False,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_three_prime_end_with_other_exons_is_2d_major_likely_pathogenic():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2D_MAJOR", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        three_prime_end_deleted=True, other_exons_involved=True,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2D"
    assert result.points == 0.9
    assert result.provisional_class == ProvisionalClass.LIKELY_PATHOGENIC


def test_three_prime_end_last_exon_cds_only_is_2d_minor_vus():
    ev = CnvDeletionEvidence(
        cnv_id="EDGE_2D_MINOR", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        three_prime_end_deleted=True, cds_involved=True, other_exons_involved=False,
    )
    result = cnv_scoring.score_cnv_deletion(ev, _dmd_dosage_config())
    assert result.categories[0].code == "2D"
    assert result.points == 0.3
    assert result.provisional_class == ProvisionalClass.VUS


# ---------------------------------------------------------- model validation

def test_whole_gene_and_benign_region_is_contradictory():
    expect_schema_error(lambda: CnvDeletionEvidence(
        cnv_id="BAD1", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_deleted=True, overlaps_benign_region=True,
    ))


def test_purely_intragenic_deletion_requires_reading_frame_effect():
    expect_schema_error(lambda: CnvDeletionEvidence(
        cnv_id="BAD2", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
    ))


def test_out_of_frame_requires_nmd_predicted_stated():
    expect_schema_error(lambda: CnvDeletionEvidence(
        cnv_id="BAD3", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=False,
        reading_frame_effect=CnvReadingFrameEffect.OUT_OF_FRAME,
    ))


def test_unverified_coordinates_cannot_carry_position_fields():
    expect_schema_error(lambda: CnvDeletionEvidence(
        cnv_id="BAD4", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=True,
        coordinate_verified=False, chromosome="X",
    ))


def test_end_must_be_after_start():
    expect_schema_error(lambda: CnvDeletionEvidence(
        cnv_id="BAD5", gene="DMD", genome_build=GenomeBuild.GRCH38, whole_gene_deleted=True,
        start=100, end=50,
    ))


def test_cnv_category_result_rejects_unknown_code():
    from variant_classifier.models import CnvCategoryResult
    from variant_classifier.models.enums import CriterionStatus, EvidenceDirection

    expect_schema_error(lambda: CnvCategoryResult(
        code="9Z", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC, points=1.0,
        rule_source="x", rule_version="1", rationale="test",
    ))


def test_cnv_category_result_met_requires_nonzero_points():
    from variant_classifier.models import CnvCategoryResult
    from variant_classifier.models.enums import CriterionStatus, EvidenceDirection

    expect_schema_error(lambda: CnvCategoryResult(
        code="2A", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC, points=0.0,
        rule_source="x", rule_version="1", rationale="test",
    ))


# ==================================================== batch 24: duplication (gain) scoring

# ------------------------------------------------------- against real fixtures

def test_cnv_duplication_scoring_matches_golden_case_for_all_curated_evidence():
    evidence, rejected = loader.load_cnv_duplication_evidence()
    assert rejected == []
    golden_cases = loader.load_cnv_duplication_golden_cases()

    assert len(evidence) == 2
    checked = 0
    for ev in evidence:
        golden = golden_cases[ev.cnv_id]
        result = cnv_scoring.score_cnv_duplication(ev)
        assert result.provisional_class == golden["expected_provisional_class"], (
            f"{ev.cnv_id}: scored {result.provisional_class}, golden case expects "
            f"{golden['expected_provisional_class']}. Rationale: {result.rationale}"
        )
        if golden["expected_points"] is not None:
            assert result.points == golden["expected_points"], (
                f"{ev.cnv_id}: scored {result.points} points, golden case expects {golden['expected_points']}"
            )
        if golden["expected_category_code"] is not None:
            assert result.categories[0].code == golden["expected_category_code"], (
                f"{ev.cnv_id}: category {result.categories[0].code}, golden case expects "
                f"{golden['expected_category_code']}"
            )
        checked += 1
    assert checked == 2


def test_cnv_duplication_golden_cases_cover_every_curated_cnv_id():
    evidence, _ = loader.load_cnv_duplication_evidence()
    golden_cases = loader.load_cnv_duplication_golden_cases()
    evidence_ids = {ev.cnv_id for ev in evidence}
    assert evidence_ids == set(golden_cases.keys())


# ---------------------------------------------------------- hand-built decision-tree edge cases

def test_benign_region_overlap_is_gain_benign():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_GAIN_BENIGN", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=False, overlaps_benign_region=True,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "GAIN_BENIGN"
    assert result.points == -1.0
    assert result.provisional_class == ProvisionalClass.BENIGN


def test_whole_gene_duplication_is_none_applicable_not_triplosensitivity():
    # This project deliberately does not implement TS scoring -- see
    # cnv_duplication_evidence.py's module docstring for why.
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_WHOLE_DUP", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=True, breakpoint_within_gene=False,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.points == 0.0
    assert result.provisional_class == ProvisionalClass.VUS


def test_no_breakpoint_within_gene_and_no_other_overlap_is_none_applicable():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_NO_OVERLAP", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=False,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_breakpoint_within_gene_not_tandem_is_none_applicable():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_NOT_TANDEM", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.NOT_TANDEM_OR_COMPLEX,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_breakpoint_within_gene_tandem_unknown_orientation_is_none_applicable():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_UNKNOWN_TANDEM", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.UNKNOWN,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "NONE_APPLICABLE"
    assert result.provisional_class == ProvisionalClass.VUS


def test_tandem_out_of_frame_is_gain_2k_equiv():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_2K", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.TANDEM, reading_frame_effect=CnvReadingFrameEffect.OUT_OF_FRAME,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "GAIN_2K_EQUIV"
    assert result.points == 0.45
    assert result.provisional_class == ProvisionalClass.VUS  # 0.45 < 0.90 Likely Pathogenic cutoff


def test_tandem_in_frame_is_gain_2j_equiv():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_2J", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.TANDEM, reading_frame_effect=CnvReadingFrameEffect.IN_FRAME,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "GAIN_2J_EQUIV"
    assert result.points == 0.0
    assert result.provisional_class == ProvisionalClass.VUS


def test_tandem_unknown_frame_effect_is_gain_2j_equiv():
    ev = CnvDuplicationEvidence(
        cnv_id="EDGE_2J_UNKNOWN_FRAME", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.TANDEM, reading_frame_effect=CnvReadingFrameEffect.UNKNOWN,
    )
    result = cnv_scoring.score_cnv_duplication(ev)
    assert result.categories[0].code == "GAIN_2J_EQUIV"
    assert result.points == 0.0


# ---------------------------------------------------------- model validation

def test_whole_gene_duplicated_and_benign_region_is_contradictory():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD1", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=True, breakpoint_within_gene=False, overlaps_benign_region=True,
    ))


def test_whole_gene_duplicated_and_breakpoint_within_gene_is_contradictory():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD2", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=True, breakpoint_within_gene=True,
    ))


def test_breakpoint_within_gene_requires_is_tandem_stated():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD3", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
    ))


def test_confirmed_tandem_requires_reading_frame_effect_stated():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD4", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=False, breakpoint_within_gene=True,
        is_tandem=CnvDuplicationOrientation.TANDEM,
    ))


def test_duplication_unverified_coordinates_cannot_carry_position_fields():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD5", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=True, breakpoint_within_gene=False,
        coordinate_verified=False, chromosome="X",
    ))


def test_duplication_end_must_be_after_start():
    expect_schema_error(lambda: CnvDuplicationEvidence(
        cnv_id="DBAD6", gene="DMD", genome_build=GenomeBuild.GRCH38,
        whole_gene_duplicated=True, breakpoint_within_gene=False,
        start=100, end=50,
    ))


def test_cnv_category_result_accepts_gain_codes_too():
    from variant_classifier.models import CnvCategoryResult
    from variant_classifier.models.enums import CriterionStatus, EvidenceDirection

    # Should not raise -- GAIN_2K_EQUIV is a valid code shared via the
    # combined CNV_LOSS_CATEGORY_CODES | CNV_GAIN_CATEGORY_CODES check.
    result = CnvCategoryResult(
        code="GAIN_2K_EQUIV", status=CriterionStatus.MET, direction=EvidenceDirection.PATHOGENIC,
        points=0.45, rule_source="x", rule_version="1", rationale="test",
    )
    assert result.code == "GAIN_2K_EQUIV"
