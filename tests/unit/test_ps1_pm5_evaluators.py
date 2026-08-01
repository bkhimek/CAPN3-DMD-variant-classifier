"""Tests for the PS1 and PM5 evaluators, added this round -- the first
new criteria since PM4 (batch 14). Same two-part structure as every other
evaluator's tests: (1) run against every curated fixture that has a
PS1/PM5 expectation and check it matches the golden case; (2) hand-built
edge cases for branches the curated fixtures don't happen to cover.
"""

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators import evaluate_pm5, evaluate_ps1
from variant_classifier.models import (
    GeneDiseaseContext,
    PopulationEvidence,
    SameResidueEvidence,
    Specification,
    TranscriptConsequence,
    VariantEvidenceBundle,
    VariantIdentity,
)
from variant_classifier.models.enums import (
    Consequence,
    CriterionStatus,
    CriterionStrength,
    DiseaseMechanism,
    GenomeBuild,
    Inheritance,
    PopulationRetrievalStatus,
    SpecificationType,
)


def expect_schema_error(callable_):
    try:
        callable_()
    except SchemaValidationError:
        return
    raise AssertionError("expected SchemaValidationError, none was raised")


# ------------------------------------------------------- against real fixtures

def test_ps1_and_pm5_match_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        for code, evaluate in (("PS1", evaluate_ps1), ("PM5", evaluate_pm5)):
            if code not in golden.expected_criterion_status:
                continue
            result = evaluate(bundle)
            expected = golden.expected_criterion_status[code]
            assert result.status == expected, (
                f"{bundle.variant.variant_id}: {code} evaluated to {result.status}, "
                f"golden case expects {expected}. Rationale: {result.rationale}"
            )
            checked += 1
    assert checked == 2 * len(bundles) > 0


# ------------------------------------------------------- hand-built edge cases

def _spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="2015")


def _population():
    return PopulationEvidence(
        source="gnomAD", source_version="test", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True,
    )


def _bundle(transcript: TranscriptConsequence, same_residue_evidence=None, gene="CAPN3"):
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene=gene, genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene=gene, disease="test disease", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_spec(),
        ),
        transcript_consequences=[transcript],
        population_evidence=[_population()],
        same_residue_evidence=same_residue_evidence,
    )


def _missense(hgvs_p="p.Arg1His"):
    return TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT, hgvs_p=hgvs_p,
    )


def test_ps1_non_missense_is_not_applicable():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.STOP_GAINED, nmd_predicted=True,
    )
    result = evaluate_ps1(_bundle(transcript))
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_pm5_non_missense_is_not_applicable():
    transcript = TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.INFRAME_DELETION, repeat_region=False,
    )
    result = evaluate_pm5(_bundle(transcript))
    assert result.status == CriterionStatus.NOT_APPLICABLE


def test_ps1_missense_with_no_same_residue_evidence_is_not_evaluated():
    result = evaluate_ps1(_bundle(_missense()))
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_pm5_missense_with_same_residue_evidence_present_but_pm5_field_unset_is_not_evaluated():
    # ps1 curated, pm5 deliberately left unset -- each field is independent.
    sre = SameResidueEvidence(
        ps1_precedent_established=False,
    )
    result = evaluate_pm5(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_ps1_precedent_checked_and_not_found_is_not_met():
    sre = SameResidueEvidence(ps1_precedent_established=False)
    result = evaluate_ps1(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.NOT_MET


def test_pm5_precedent_checked_and_not_found_is_not_met():
    sre = SameResidueEvidence(pm5_precedent_established=False)
    result = evaluate_pm5(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.NOT_MET


def test_ps1_precedent_found_but_splice_impact_not_excluded_is_manual_review():
    sre = SameResidueEvidence(
        ps1_precedent_established=True,
        ps1_precedent_classification="PATHOGENIC",
        ps1_precedent_variant="TESTGENE c.1A>G p.(Arg1His), synthetic precedent",
        splice_impact_excluded=False,
    )
    result = evaluate_ps1(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MANUAL_REVIEW


def test_pm5_precedent_found_but_splice_impact_not_excluded_is_manual_review():
    sre = SameResidueEvidence(
        pm5_precedent_established=True,
        pm5_precedent_classification="PATHOGENIC",
        pm5_precedent_variant="TESTGENE c.1A>G p.(Arg1Cys), synthetic precedent",
        splice_impact_excluded=False,
    )
    result = evaluate_pm5(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MANUAL_REVIEW


def test_ps1_pathogenic_precedent_is_met_strong():
    sre = SameResidueEvidence(
        ps1_precedent_established=True,
        ps1_precedent_classification="PATHOGENIC",
        ps1_precedent_variant="TESTGENE c.1A>G p.(Arg1His), synthetic precedent",
        splice_impact_excluded=True,
    )
    result = evaluate_ps1(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STRONG


def test_ps1_likely_pathogenic_precedent_is_met_moderate_downgrade():
    sre = SameResidueEvidence(
        ps1_precedent_established=True,
        ps1_precedent_classification="LIKELY_PATHOGENIC",
        ps1_precedent_variant="TESTGENE c.1A>G p.(Arg1His), synthetic precedent",
        splice_impact_excluded=True,
    )
    result = evaluate_ps1(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE
    assert "downgraded" in result.rationale


def test_pm5_pathogenic_precedent_is_met_moderate():
    sre = SameResidueEvidence(
        pm5_precedent_established=True,
        pm5_precedent_classification="PATHOGENIC",
        pm5_precedent_variant="TESTGENE c.1A>G p.(Arg1Cys), synthetic precedent",
        splice_impact_excluded=True,
    )
    result = evaluate_pm5(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.MODERATE


def test_pm5_likely_pathogenic_precedent_is_met_supporting_downgrade():
    sre = SameResidueEvidence(
        pm5_precedent_established=True,
        pm5_precedent_classification="LIKELY_PATHOGENIC",
        pm5_precedent_variant="TESTGENE c.1A>G p.(Arg1Cys), synthetic precedent",
        splice_impact_excluded=True,
    )
    result = evaluate_pm5(_bundle(_missense(), same_residue_evidence=sre))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING
    assert "downgraded" in result.rationale


# ------------------------------------------------------- SameResidueEvidence model validation

def test_same_residue_evidence_requires_classification_when_ps1_precedent_established():
    expect_schema_error(lambda: SameResidueEvidence(
        ps1_precedent_established=True,
        ps1_precedent_variant="TESTGENE c.1A>G p.(Arg1His)",
        splice_impact_excluded=True,
    ))


def test_same_residue_evidence_requires_variant_citation_when_pm5_precedent_established():
    expect_schema_error(lambda: SameResidueEvidence(
        pm5_precedent_established=True,
        pm5_precedent_classification="PATHOGENIC",
        splice_impact_excluded=True,
    ))


def test_same_residue_evidence_rejects_classification_without_established_true():
    expect_schema_error(lambda: SameResidueEvidence(
        ps1_precedent_established=False,
        ps1_precedent_classification="PATHOGENIC",
    ))


def test_same_residue_evidence_requires_splice_impact_excluded_when_any_precedent_established():
    expect_schema_error(lambda: SameResidueEvidence(
        pm5_precedent_established=True,
        pm5_precedent_classification="PATHOGENIC",
        pm5_precedent_variant="TESTGENE c.1A>G p.(Arg1Cys)",
        # splice_impact_excluded deliberately omitted
    ))


def test_same_residue_evidence_rejects_invalid_classification_value():
    expect_schema_error(lambda: SameResidueEvidence(
        ps1_precedent_established=True,
        ps1_precedent_classification="BENIGN",  # not a valid PS1/PM5 precedent classification
        ps1_precedent_variant="TESTGENE c.1A>G p.(Arg1His)",
        splice_impact_excluded=True,
    ))
