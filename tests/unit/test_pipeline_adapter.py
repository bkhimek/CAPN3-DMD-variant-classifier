"""Tests for pipeline_adapter.py against a real (trimmed, not synthetic)
sample of CAPN3-DMD-variant-calling-pipeline's actual
ANNOTATE_CALLS output — see data/source/pipeline_annotate_calls/README.md for
exactly which records are in the fixture and why each one is there.
"""

from pathlib import Path

from variant_classifier import loader
from variant_classifier.pipeline_adapter import _split_per_allele, _vep_allele, build_bundles_from_pipeline_output
from variant_classifier.models.enums import Consequence, PopulationRetrievalStatus

REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_VCF = REPO_ROOT / "data" / "source" / "pipeline_annotate_calls" / "adapter_test_sample.vcf.gz"


def _load_sample_bundles():
    contexts = loader.load_gene_disease_contexts()
    return build_bundles_from_pipeline_output(SAMPLE_VCF, contexts)


def test_vep_allele_trimming_matches_real_vep_output():
    # Confirmed against CAPN3-DMD-variant-calling-pipeline's actual VEP
    # output for these exact REF/ALT pairs (chr15:42262438, not itself in this fixture) before
    # writing this rule into the adapter.
    assert _vep_allele("CTTT", "C") == "-"
    assert _vep_allele("CTTT", "CTT") == "TT"
    assert _vep_allele("G", "A") == "A"


def test_split_per_allele_handles_missing_token():
    assert _split_per_allele("95102,.", 2) == ["95102", None]
    assert _split_per_allele(None, 2) == [None, None]


def test_no_rejected_records_in_real_sample():
    bundles, rejected = _load_sample_bundles()
    assert rejected == []


def test_record_outside_mane_transcript_span_produces_no_bundle():
    # chr15:42259883 falls in CAPN3's padded flanking region but has no CSQ
    # hit on either gene's MANE transcript -- correctly produces zero
    # bundles rather than a fabricated one (see evidence_bundle.py's
    # requirement of exactly one clinically_relevant transcript).
    bundles, _ = _load_sample_bundles()
    assert not any(b.variant.position == 42259883 for b in bundles)


def test_real_capn3_missense_bundle():
    bundles, _ = _load_sample_bundles()
    matches = [b for b in bundles if b.variant.position == 42389001]
    assert len(matches) == 1
    bundle = matches[0]
    assert bundle.variant.gene == "CAPN3"
    assert bundle.variant.chromosome == "chr15"
    assert bundle.variant.reference == "G"
    assert bundle.variant.alternate == "A"
    assert bundle.gene_disease_context.gene == "CAPN3"
    tc = bundle.transcript_consequences[0]
    assert tc.transcript_id == "NM_000070.3"
    assert tc.clinically_relevant is True
    assert tc.mane_select is True
    assert tc.consequence == Consequence.MISSENSE_VARIANT
    assert tc.hgvs_c == "ENST00000397163.8:c.706G>A"
    assert tc.hgvs_p == "ENSP00000380349.1:p.Ala236Thr"
    pop = bundle.population_evidence[0]
    assert pop.retrieval_status == PopulationRetrievalStatus.OBSERVED
    assert pop.overall_af == 0.222151
    assert pop.allele_count == 33763
    assert pop.allele_number == 151982


def test_real_dmd_missense_bundle():
    bundles, _ = _load_sample_bundles()
    matches = [b for b in bundles if b.variant.position == 31478233]
    assert len(matches) == 1
    bundle = matches[0]
    assert bundle.variant.gene == "DMD"
    tc = bundle.transcript_consequences[0]
    assert tc.transcript_id == "NM_004006.3"
    assert tc.consequence == Consequence.MISSENSE_VARIANT
    assert tc.hgvs_c == "ENST00000357033.9:c.8810G>A"
    pop = bundle.population_evidence[0]
    assert pop.retrieval_status == PopulationRetrievalStatus.OBSERVED
    assert pop.overall_af == 0.938899


def test_multiallelic_capn3_site_splits_into_two_bundles_with_distinct_gnomad_status():
    bundles, _ = _load_sample_bundles()
    matches = [b for b in bundles if b.variant.position == 42398590]
    assert len(matches) == 2
    statuses = {b.variant.alternate: b.population_evidence[0].retrieval_status for b in matches}
    assert statuses["T"] == PopulationRetrievalStatus.OBSERVED
    assert statuses["TACACAC"] == PopulationRetrievalStatus.ABSENT
    for b in matches:
        assert b.variant.gene == "CAPN3"
        assert b.transcript_consequences[0].consequence == Consequence.INTRON_VARIANT


def test_multiallelic_dmd_site_splits_into_two_bundles():
    bundles, _ = _load_sample_bundles()
    matches = [b for b in bundles if b.variant.position == 32456578]
    assert len(matches) == 2
    assert all(b.variant.gene == "DMD" for b in matches)
    alts = {b.variant.alternate for b in matches}
    assert alts == {"T", "TTG"}
