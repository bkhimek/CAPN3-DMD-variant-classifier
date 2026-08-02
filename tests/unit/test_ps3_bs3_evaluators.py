"""Tests for the PS3 and BS3 evaluators, added batch 25 -- the first new
criteria since PS1/PM5 (batch 22). Same two-part structure as every other
evaluator's tests: (1) run against every curated fixture that has a
PS3/BS3 expectation and check it matches the golden case; (2) hand-built
edge cases for branches the curated fixtures don't happen to cover --
notably BS3's MET branch, which no real curated fixture reaches (see
functional_evidence.py's docstring and the batch-25 README design note
for why no real "clean normal-WB-confirms-benign" CAPN3/DMD example was
used).
"""

from variant_classifier import loader
from variant_classifier.errors import SchemaValidationError
from variant_classifier.evaluators import evaluate_bs3, evaluate_ps3
from variant_classifier.models import (
    FunctionalEvidence,
    GeneDiseaseContext,
    PopulationEvidence,
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
    FunctionalAssayResult,
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

def test_ps3_and_bs3_match_golden_case_for_all_curated_bundles():
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    golden_cases = loader.load_golden_cases()

    checked = 0
    for bundle in bundles:
        golden = golden_cases[bundle.variant.variant_id]
        for code, evaluate in (("PS3", evaluate_ps3), ("BS3", evaluate_bs3)):
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


def test_ps3_met_supporting_for_capn3_c_550del_real_functional_data():
    # The real, ClinVar-Pathogenic founder allele: Czech LGMD2A cohort
    # Western blot data (Chrobakova et al. 2004 / Hermanova et al. 2006)
    # showing total absence of calpain-3 protein. See the fixture's
    # notes in data/curated/variant_evidence.json for the full citation.
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.550del")
    result = evaluate_ps3(bundle)
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.SUPPORTING


def test_ps3_and_bs3_both_not_met_for_capn3_c_2257g_a_indeterminate_assay():
    # The real, genuinely-contested p.Asp753Asn variant: Bruno et al.
    # 2025's own Western blot data for this exact variant was mixed
    # (Normal in 3/5, Reduced in 2/5) -- curated as INDETERMINATE, which
    # supports neither PS3 nor BS3.
    bundles, rejected = loader.load_variant_evidence_bundles()
    assert rejected == []
    bundle = next(b for b in bundles if b.variant.variant_id == "CAPN3_c.2257G>A")
    assert evaluate_ps3(bundle).status == CriterionStatus.NOT_MET
    assert evaluate_bs3(bundle).status == CriterionStatus.NOT_MET


# ------------------------------------------------------- hand-built edge cases

def _spec():
    return Specification(type=SpecificationType.GENERIC_ACMG, version="2015")


def _population():
    return PopulationEvidence(
        source="gnomAD", source_version="test", retrieval_status=PopulationRetrievalStatus.ABSENT,
        locus_coverage_adequate=True,
    )


def _transcript():
    return TranscriptConsequence(
        transcript_id="NM_1", clinically_relevant=True, consequence=Consequence.MISSENSE_VARIANT, hgvs_p="p.Arg1His",
    )


def _bundle(functional_evidence=None, gene="CAPN3"):
    return VariantEvidenceBundle(
        variant=VariantIdentity(variant_id="EDGE_CASE", gene=gene, genome_build=GenomeBuild.GRCH38),
        gene_disease_context=GeneDiseaseContext(
            gene=gene, disease="test disease", inheritance=Inheritance.AUTOSOMAL_RECESSIVE,
            mechanism=DiseaseMechanism.LOSS_OF_FUNCTION, lof_established=True, specification=_spec(),
        ),
        transcript_consequences=[_transcript()],
        population_evidence=[_population()],
        functional_evidence=functional_evidence,
    )


def test_ps3_no_functional_evidence_is_not_evaluated():
    result = evaluate_ps3(_bundle())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_bs3_no_functional_evidence_is_not_evaluated():
    result = evaluate_bs3(_bundle())
    assert result.status == CriterionStatus.NOT_EVALUATED


def test_ps3_and_bs3_both_not_met_when_assay_indeterminate():
    fe = FunctionalEvidence(assay_result=FunctionalAssayResult.INDETERMINATE)
    bundle = _bundle(functional_evidence=fe)
    assert evaluate_ps3(bundle).status == CriterionStatus.NOT_MET
    assert evaluate_bs3(bundle).status == CriterionStatus.NOT_MET


def test_ps3_not_met_when_assay_normal():
    fe = FunctionalEvidence(
        assay_result=FunctionalAssayResult.NORMAL, validation_strength=CriterionStrength.SUPPORTING,
    )
    result = evaluate_ps3(_bundle(functional_evidence=fe))
    assert result.status == CriterionStatus.NOT_MET


def test_bs3_not_met_when_assay_abnormal():
    fe = FunctionalEvidence(
        assay_result=FunctionalAssayResult.ABNORMAL, validation_strength=CriterionStrength.SUPPORTING,
    )
    result = evaluate_bs3(_bundle(functional_evidence=fe))
    assert result.status == CriterionStatus.NOT_MET


def test_ps3_met_strength_matches_curated_validation_strength():
    for strength in (CriterionStrength.SUPPORTING, CriterionStrength.MODERATE, CriterionStrength.STRONG):
        fe = FunctionalEvidence(assay_result=FunctionalAssayResult.ABNORMAL, validation_strength=strength)
        result = evaluate_ps3(_bundle(functional_evidence=fe))
        assert result.status == CriterionStatus.MET
        assert result.strength == strength


def test_bs3_met_at_strong_no_real_curated_fixture_reaches_this_branch():
    # No real CAPN3/DMD fixture in this project's curated set has a
    # clean, non-caveated normal-WB-confirms-benign result (see
    # functional_evidence.py's docstring) -- covered here as a hand-built
    # edge case instead, same convention used for CNV categories 2C/2D/2F
    # in batch 23.
    fe = FunctionalEvidence(assay_result=FunctionalAssayResult.NORMAL, validation_strength=CriterionStrength.STRONG)
    result = evaluate_bs3(_bundle(functional_evidence=fe))
    assert result.status == CriterionStatus.MET
    assert result.strength == CriterionStrength.STRONG


# ------------------------------------------------------- FunctionalEvidence model validation

def test_functional_evidence_requires_validation_strength_when_abnormal():
    expect_schema_error(lambda: FunctionalEvidence(assay_result=FunctionalAssayResult.ABNORMAL))


def test_functional_evidence_requires_validation_strength_when_normal():
    expect_schema_error(lambda: FunctionalEvidence(assay_result=FunctionalAssayResult.NORMAL))


def test_functional_evidence_rejects_validation_strength_when_indeterminate():
    expect_schema_error(lambda: FunctionalEvidence(
        assay_result=FunctionalAssayResult.INDETERMINATE, validation_strength=CriterionStrength.SUPPORTING,
    ))


def test_functional_evidence_rejects_very_strong():
    expect_schema_error(lambda: FunctionalEvidence(
        assay_result=FunctionalAssayResult.ABNORMAL, validation_strength=CriterionStrength.VERY_STRONG,
    ))


def test_functional_evidence_rejects_stand_alone():
    expect_schema_error(lambda: FunctionalEvidence(
        assay_result=FunctionalAssayResult.NORMAL, validation_strength=CriterionStrength.STAND_ALONE,
    ))
