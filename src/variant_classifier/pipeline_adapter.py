"""Adapter: turns CAPN3-DMD-variant-calling-pipeline's VCF output into this
project's VariantEvidenceBundle instances.

Lives here (CAPN3-DMD-variant-classifier), not in CAPN3-DMD-variant-calling-pipeline,
deliberately: the pipeline stays a generic, Python-classifier-agnostic
Nextflow/Docker project (no runtime dependency on this package, consistent
with its "Docker only, no pip installs" convention), while this project owns
its own dataclasses and already has the dependency-free test runner to
verify against. See ~/projects/HANDOFF.md's "Projects 4x5 integration" entry
for the fuller rationale.

Input is CAPN3-DMD-variant-calling-pipeline's `results/annotate_calls/annotated_calls.vcf.gz` — the
GATK/DeepVariant-concordant call set, annotated with VEP (transcript
consequence) and gnomAD v4.1 (population frequency). Reading it needs no
extra dependency: bgzip is gzip-compatible for sequential reads, so stdlib
`gzip` suffices (no pysam/htslib binding required, keeping this project's
"stdlib + PyYAML only" tech-stack rule intact).

What this adapter deliberately does NOT do:
- ComputationalEvidence: never populated. Project 4's own design requires a
  single *calibrated* predictor score (computational_evidence.py's docstring
  is explicit that raw per-tool scores don't qualify); that calibration is a
  research task in its own right, the same bucket as PS3/BS3 on the Roadmap.
  Bundles built here are still fully valid (the field is Optional) — PP3/BP4
  just won't fire on them.
- NMD prediction / repeat-region calls for consequences that need them
  (frameshift_variant, stop_gained, inframe_deletion, inframe_insertion,
  stop_lost): none of the real HG002 concordant calls on either gene's MANE
  transcript hit these consequences (checked directly against the real
  annotated output before writing this module), so this adapter has never
  needed to answer that question. TranscriptConsequence's own __post_init__
  already refuses to construct such a record without an explicit
  nmd_predicted/repeat_region value, so a future variant that does hit one of
  these consequences will be correctly rejected (surfaced in this module's
  `rejected` list) rather than silently accepted with a guessed value.

HG002 caveat: per ~/projects/HANDOFF.md, HG002's CAPN3/DMD variation is
confirmed ordinary benign background polymorphism (zero overlap with any
ClinVar pathogenic/likely-pathogenic position in either gene). Bundles built
from it are a real, mechanical schema/plumbing test of the pipeline-to-bundle
mapping — not a demonstration of the classifier flagging a genuine pathogenic
finding. `notes` on every bundle built here says so explicitly.
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import unquote

from .errors import SchemaValidationError
from .models import GeneDiseaseContext, PopulationEvidence, TranscriptConsequence, VariantEvidenceBundle, VariantIdentity
from .models.enums import Consequence, GenomeBuild, PopulationRetrievalStatus

# MANE Select transcripts (confirmed live against NCBI/EBI's MANE summary
# file — see CAPN3-DMD-variant-calling-pipeline's docs/data_sources.md — not carried over from
# memory). Ensembl IDs are version-stripped since VEP's --gff CSQ output
# reports bare "ENST..." Feature IDs without the trailing ".N".
MANE_TRANSCRIPTS = {
    "CAPN3": {"ensembl": "ENST00000397163", "refseq": "NM_000070.3"},
    "DMD": {"ensembl": "ENST00000357033", "refseq": "NM_004006.3"},
}

GNOMAD_SOURCE = "gnomAD"
GNOMAD_SOURCE_VERSION = "v4.1.0"

_MISSING_VCF_TOKENS = (".", "")


def _vep_allele(ref: str, alt: str) -> str:
    """Reproduce VEP's minimal-representation Allele string for an indel so
    it can be matched against CSQ's own Allele field. VEP (and the VCF
    convention it follows) trims exactly one shared leading anchor base for
    indels, representing a pure deletion as "-". Confirmed against real
    output: REF=CTTT/ALT=C -> "-", REF=CTTT/ALT=CTT -> "TT".
    """
    if len(ref) == 1 and len(alt) == 1:
        return alt
    if ref and alt and ref[0] == alt[0]:
        trimmed = alt[1:]
        return trimmed if trimmed else "-"
    return alt


def _parse_csq_format(header_lines: List[str]) -> List[str]:
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ,"):
            marker = 'Format: '
            idx = line.find(marker)
            if idx == -1:
                raise SchemaValidationError("VCF CSQ INFO header found but has no 'Format: ...' description")
            fields_part = line[idx + len(marker):].rstrip('">')
            return fields_part.split("|")
    raise SchemaValidationError("VCF has no ##INFO=<ID=CSQ,...> header — was this file run through VEP_ANNOTATE?")


def _parse_info(info_str: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            key, _, value = entry.partition("=")
            info[key] = value
        else:
            info[entry] = "true"
    return info


def _split_per_allele(value: Optional[str], n_alts: int) -> List[Optional[str]]:
    """gnomAD's per-ALT INFO fields (Number=A) are comma-separated, aligned
    1:1 with the VCF ALT list; "." marks an allele gnomAD's region query
    didn't find. gnomAD_AN (Number=1, site-wide) never goes through this."""
    if value is None:
        return [None] * n_alts
    parts = value.split(",")
    if len(parts) != n_alts:
        raise SchemaValidationError(f"expected {n_alts} comma-separated values, got {len(parts)}: {value!r}")
    return [None if p in _MISSING_VCF_TOKENS else p for p in parts]


def _consequence_from_vep_term(raw: str) -> Consequence:
    """VEP joins multiple SO terms with '&', ranked most-to-least severe;
    take the first. Falls back to OTHER for any term outside this project's
    deliberately small Consequence enum (see enums.py's own docstring — this
    is the intended behavior, not a gap)."""
    primary = raw.split("&")[0]
    try:
        return Consequence(primary)
    except ValueError:
        return Consequence.OTHER


def _build_transcript_consequence(csq_fields: List[str], csq_entry: Dict[str, str], mane_refseq: str) -> TranscriptConsequence:
    consequence = _consequence_from_vep_term(csq_entry.get("Consequence", ""))
    nmd_predicted = None
    repeat_region = None
    if consequence in TranscriptConsequence.NMD_RELEVANT_CONSEQUENCES or consequence in TranscriptConsequence.PM4_RELEVANT_CONSEQUENCES:
        # Deliberately left unset — see module docstring. TranscriptConsequence's
        # own __post_init__ will raise SchemaValidationError for this case,
        # which the caller (build_bundles_from_pipeline_output) catches into
        # `rejected` rather than fabricating a value here.
        pass
    hgvs_c = csq_entry.get("HGVSc") or None
    hgvs_p = csq_entry.get("HGVSp") or None
    if hgvs_p:
        hgvs_p = unquote(hgvs_p)  # VEP percent-encodes "=" as %3D for p.Xxx123= synonymous notation
    exon = csq_entry.get("EXON") or None
    return TranscriptConsequence(
        transcript_id=mane_refseq,
        clinically_relevant=True,
        consequence=consequence,
        mane_select=True,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        exon=exon,
        nmd_predicted=nmd_predicted,
        repeat_region=repeat_region,
    )


def _build_population_evidence(info: Dict[str, str], allele_idx: int, n_alts: int) -> PopulationEvidence:
    af_values = _split_per_allele(info.get("gnomAD_AF"), n_alts)
    af = af_values[allele_idx]
    if af is None:
        return PopulationEvidence(
            source=GNOMAD_SOURCE,
            source_version=GNOMAD_SOURCE_VERSION,
            retrieval_status=PopulationRetrievalStatus.ABSENT,
            locus_coverage_adequate=True,
        )
    ac_values = _split_per_allele(info.get("gnomAD_AC"), n_alts)
    grpmax_values = _split_per_allele(info.get("gnomAD_AF_grpmax"), n_alts)
    nhomalt_values = _split_per_allele(info.get("gnomAD_nhomalt"), n_alts)
    an = info.get("gnomAD_AN")  # site-wide (Number=1), not per-allele
    return PopulationEvidence(
        source=GNOMAD_SOURCE,
        source_version=GNOMAD_SOURCE_VERSION,
        retrieval_status=PopulationRetrievalStatus.OBSERVED,
        overall_af=float(af),
        ancestry_specific_max_af=float(grpmax_values[allele_idx]) if grpmax_values[allele_idx] is not None else None,
        allele_count=int(ac_values[allele_idx]) if ac_values[allele_idx] is not None else None,
        allele_number=int(an) if an not in (None, *_MISSING_VCF_TOKENS) else None,
        homozygote_count=int(nhomalt_values[allele_idx]) if nhomalt_values[allele_idx] is not None else None,
    )


def build_bundles_from_pipeline_output(
    vcf_path: Path,
    gene_disease_contexts: Dict[str, GeneDiseaseContext],
    sample_label: str = "HG002",
) -> Tuple[List[VariantEvidenceBundle], List[Tuple[dict, str]]]:
    """Parse CAPN3-DMD-variant-calling-pipeline's annotated_calls.vcf.gz into VariantEvidenceBundle
    instances. Returns (bundles, rejected) like loader.load_variant_evidence_bundles
    — one bad/unsupported record shouldn't abort the whole batch.

    A VCF row is split per-ALT allele (VCF's own unit-of-truth is one row per
    site, potentially multiple ALTs; VariantEvidenceBundle is one bundle per
    variant). A row/allele is silently skipped (not counted as rejected) only
    when NEITHER gene's MANE transcript has any CSQ entry for it — i.e. it
    falls outside both padded regions' actual gene bodies (pure flanking/
    intergenic sequence), so there is no clinically-relevant transcript
    consequence to build a bundle from at all, by design of this project's
    schema (evidence_bundle.py requires exactly one).
    """
    vcf_path = Path(vcf_path)
    bundles: List[VariantEvidenceBundle] = []
    rejected: List[Tuple[dict, str]] = []

    header_lines: List[str] = []
    csq_fields: Optional[List[str]] = None

    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                csq_fields = _parse_csq_format(header_lines)
                continue
            if csq_fields is None:
                raise SchemaValidationError(f"{vcf_path}: data line encountered before #CHROM header")

            chrom, pos, _vid, ref, alt_field, *_rest = line.split("\t")
            info_str = line.split("\t")[7]
            info = _parse_info(info_str)
            alts = alt_field.split(",")

            csq_entries = [dict(zip(csq_fields, entry.split("|"))) for entry in info.get("CSQ", "").split(",") if entry]

            for allele_idx, alt in enumerate(alts):
                vep_allele = _vep_allele(ref, alt)
                gene = None
                mane_refseq = None
                matched_csq = None
                for gene_name, ids in MANE_TRANSCRIPTS.items():
                    hit = next(
                        (c for c in csq_entries if c.get("Feature") == ids["ensembl"] and c.get("Allele") == vep_allele),
                        None,
                    )
                    if hit is not None:
                        gene, mane_refseq, matched_csq = gene_name, ids["refseq"], hit
                        break
                if matched_csq is None:
                    continue  # outside both gene bodies — no bundle to build, by design (see docstring)

                raw_record = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene,
                }
                variant_id = f"{sample_label}_{gene}_{chrom}:{pos}{ref}>{alt}"
                try:
                    variant = VariantIdentity(
                        variant_id=variant_id,
                        gene=gene,
                        genome_build=GenomeBuild.GRCH38,
                        chromosome=chrom,
                        position=int(pos),
                        reference=ref,
                        alternate=alt,
                        coordinate_verified=True,
                    )
                    gene_disease_context = gene_disease_contexts[gene]
                    transcript_consequence = _build_transcript_consequence(csq_fields, matched_csq, mane_refseq)
                    population_evidence = _build_population_evidence(info, allele_idx, len(alts))
                    bundle = VariantEvidenceBundle(
                        variant=variant,
                        gene_disease_context=gene_disease_context,
                        transcript_consequences=[transcript_consequence],
                        population_evidence=[population_evidence],
                        computational_evidence=None,
                        notes=(
                            f"Real pipeline-derived variant from CAPN3-DMD-variant-calling-pipeline's "
                            f"ANNOTATE_CALLS output (sample={sample_label}, GIAB reference individual). "
                            "Mechanical/schema-level verification of the VCF-to-VariantEvidenceBundle mapping "
                            "only -- HG002 is confirmed clinically empty for CAPN3/DMD pathogenic variation "
                            "(see ~/projects/HANDOFF.md), so this is NOT a demonstration of the classifier "
                            "flagging a genuine pathogenic finding."
                        ),
                    )
                    bundles.append(bundle)
                except (SchemaValidationError, KeyError) as exc:
                    rejected.append((raw_record, str(exc)))

    return bundles, rejected
