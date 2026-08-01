"""CnvDeletionEvidence -- curated evidence for a single-gene copy-number
deletion, scoped to DMD only this batch. See cnv_scoring.py for the
scoring logic this feeds, and README.md ("DMD CNV/structural-variant
scoring, batch 23") for the full design writeup of what is and is not
implemented.

This is a deliberately NEW, separate evidence shape, not an extension of
VariantEvidenceBundle -- batch 22's sizing note (see README, "DMD
CNV/structural-variant representation: sized, not implemented") concluded
that CNV interpretation is a structurally different framework (ACMG/ClinGen
Technical Standards for the Interpretation and Reporting of Constitutional
Copy-Number Variants, Riggs et al. 2020, Genetics in Medicine 22:245-257)
from the Richards et al. 2015 point-mutation framework every other model in
this package represents: different evidence codes (Section 1/2/3/4/5
categories, not PVS1/PM2/etc.), different combining math (summed
floating-point scores against fixed cutoffs, not Table 5 or Tavtigian et al.
2020's integer points), and a different identity shape (a genomic interval /
exon range, not a single HGVS change). Reusing VariantIdentity or
TranscriptConsequence for this would misrepresent a CNV as if it had a
single nucleotide-level HGVS description, which it does not.

Scope actually implemented (batch 23), and why:

This project could not obtain the primary Riggs et al. 2020 paper directly
(reCAPTCHA-blocked) or the official ClinGen CNV calculator (cnvcalc.
clinicalgenome.org -- a JS app that timed out on fetch). The exact category
point values used here are instead sourced from ClassifyCNV (Gurbich TA,
Ilinsky VV. "ClassifyCNV: a tool for clinical annotation of copy-number
variants." Sci Rep 10, 20375 (2020), DOI 10.1038/s41598-020-76425-3), an
open-source, published, actively-used reimplementation of the same rubric
-- code fetched directly from github.com/Genotek/ClassifyCNV
(resources.py, ClassifyCNV.py) during this batch's research. This is a
disclosed reliance on a secondary (but primary-adjacent, peer-reviewed,
executable) source for the exact numbers, not the Riggs paper's own text
-- the same "disclosed simplification, not an invented rule" treatment
this project has given other gaps (e.g. same_residue_evidence.py's PS1/PM5
precedent-strength downgrade).

Of Riggs 2020's five evidence sections, this batch implements only a slice
of Section 2 (dosage-sensitivity / haploinsufficiency-overlap categories
for LOSS/deletion only, as reimplemented by ClassifyCNV):

  2A (1.0 pts)  -- complete overlap of an established (ClinGen HI=3)
                   dosage-sensitive gene -- i.e. the whole gene is deleted.
  2C (0.9 pts)  -- the deletion includes the gene's 5' end (5'UTR or first
                   exon) AND coding sequence.
  2D (0.9/0.3)  -- the deletion includes the gene's 3' end (3'UTR or last
                   exon); 0.9 if other exons are also involved, 0.3 if the
                   deletion is confined to the last exon's CDS.
  2E (0.9 pts)  -- an intragenic (both ends of the gene intact) deletion
                   that disrupts the reading frame and is predicted to
                   trigger nonsense-mediated decay -- this project's direct
                   link to the Aartsma-Rus DMD reading-frame rule (Aartsma-
                   Rus et al. 2006, PMID 16770791; 2019 update, Human
                   Mutation): an out-of-frame internal deletion is exactly
                   a 2E-shaped deletion.
  2F (-1.0 pts) -- the deletion falls completely within an established
                   ClinGen benign copy-number region.

Explicitly NOT implemented, named rather than hidden:

- Section 1 (genomic content / "does the CNV overlap any gene at all") --
  every evidence record in this project's scope is, by construction, a
  deletion of the DMD gene, so this category can never fire here. Real
  for CNVs with no gene overlap; not modeled.
- Section 2H (predicted-but-not-established haploinsufficiency, via
  DECIPHER/pLI/LOEUF automated prediction) -- requires infrastructure
  (three external prediction datasets) this project does not have, and
  CAPN3 -- the project's other gene -- is autosomal recessive, for which
  ClinGen's own haploinsufficiency dosage-sensitivity framework does not
  straightforwardly apply the way it does to DMD's X-linked-hemizygous
  mechanism (a single-copy deletion of a recessive gene does not, by
  itself, cause the recessive disease). No real fixture would exercise
  this category honestly, so it is deferred rather than guessed at.
- Duplications (the ACMG/ClinGen gain-side rubric) entirely -- deferred,
  mirroring the deletion-only decision confirmed with the user before
  this batch began.
- Sections 3 (gene count -- moot for a single-gene CNV), 4 (case/
  case-control/population evidence), and 5 (inheritance/family history) --
  all deferred, exactly as sized in batch 22.

One real consequence of implementing Section 2 only, disclosed here and in
the golden cases: a real out-of-frame DMD deletion is very often clinically
classified Pathogenic once case/family evidence (Sections 4/5) is folded
in, but this milestone's Section-2-only score for the same deletion can
land Likely Pathogenic (2E = 0.9, just under the 0.99 Pathogenic cutoff) --
an intentional, disclosed gap, the same "documented partial result, not a
silent misrepresentation" treatment CAPN3's PVS1+PM2 Table-5-vs-Bayesian
discrepancy already has in bayesian.py.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_bool, optional_int, optional_str, require_bool, require_dict, require_str
from .enums import CnvReadingFrameEffect, GenomeBuild


@dataclass(frozen=True)
class CnvDeletionEvidence:
    cnv_id: str
    gene: str
    genome_build: GenomeBuild
    whole_gene_deleted: bool
    chromosome: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    coordinate_verified: bool = True
    five_prime_end_deleted: bool = False
    three_prime_end_deleted: bool = False
    cds_involved: bool = False
    other_exons_involved: bool = False
    overlaps_benign_region: bool = False
    reading_frame_effect: Optional[CnvReadingFrameEffect] = None
    nmd_predicted: Optional[bool] = None
    exon_description: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"CnvDeletionEvidence[{self.cnv_id}]"
        if self.start is not None and self.start <= 0:
            raise SchemaValidationError(f"{context}: start must be a positive integer if provided")
        if self.end is not None and self.end <= 0:
            raise SchemaValidationError(f"{context}: end must be a positive integer if provided")
        if self.start is not None and self.end is not None and self.end <= self.start:
            raise SchemaValidationError(f"{context}: end ({self.end}) must be greater than start ({self.start})")
        if not self.coordinate_verified and (self.chromosome or self.start or self.end):
            raise SchemaValidationError(
                f"{context}: coordinate_verified=False but chromosome/start/end were provided -- "
                "either verify them and set coordinate_verified=True, or omit them "
                "(same convention as VariantIdentity)"
            )
        if self.whole_gene_deleted and self.overlaps_benign_region:
            raise SchemaValidationError(
                f"{context}: whole_gene_deleted=True and overlaps_benign_region=True are contradictory -- "
                "a deletion cannot both remove the entire established dosage-sensitive gene (2A) and lie "
                "completely inside an established benign region (2F)"
            )
        is_purely_intragenic = (
            not self.whole_gene_deleted
            and not self.overlaps_benign_region
            and not self.five_prime_end_deleted
            and not self.three_prime_end_deleted
        )
        if is_purely_intragenic and self.reading_frame_effect is None:
            raise SchemaValidationError(
                f"{context}: reading_frame_effect must be explicitly stated for a deletion that does not "
                "touch either end of the gene and is not a whole-gene or established-benign-region overlap "
                "-- category 2E (intragenic frameshift + NMD) cannot be safely evaluated from an unstated "
                "reading-frame effect. Never left unstated, same 'never silently guess' convention as "
                "TranscriptConsequence.nmd_predicted and SameResidueEvidence.splice_impact_excluded."
            )
        if self.reading_frame_effect == CnvReadingFrameEffect.OUT_OF_FRAME and self.nmd_predicted is None:
            raise SchemaValidationError(
                f"{context}: nmd_predicted must be explicitly true/false when reading_frame_effect is "
                "OUT_OF_FRAME -- an out-of-frame deletion could still land in the gene's last exon and "
                "escape nonsense-mediated decay, exactly the same real caveat "
                "TranscriptConsequence.nmd_predicted already encodes for point-variant frameshifts."
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CnvDeletionEvidence":
        data = require_dict(data, context or "CnvDeletionEvidence")
        cnv_id = require_str(data, "cnv_id", context or "CnvDeletionEvidence")
        ctx = f"CnvDeletionEvidence[{cnv_id}]"
        gene = require_str(data, "gene", ctx)
        genome_build = coerce_enum(GenomeBuild, data.get("genome_build"), "genome_build", ctx)
        whole_gene_deleted = require_bool(data, "whole_gene_deleted", ctx)
        reading_frame_raw = data.get("reading_frame_effect")
        reading_frame_effect = (
            coerce_enum(CnvReadingFrameEffect, reading_frame_raw, "reading_frame_effect", ctx)
            if reading_frame_raw is not None
            else None
        )
        return cls(
            cnv_id=cnv_id,
            gene=gene,
            genome_build=genome_build,
            whole_gene_deleted=whole_gene_deleted,
            chromosome=optional_str(data, "chromosome"),
            start=optional_int(data, "start", ctx, minimum=1),
            end=optional_int(data, "end", ctx, minimum=1),
            coordinate_verified=bool(data.get("coordinate_verified", True)),
            five_prime_end_deleted=bool(data.get("five_prime_end_deleted", False)),
            three_prime_end_deleted=bool(data.get("three_prime_end_deleted", False)),
            cds_involved=bool(data.get("cds_involved", False)),
            other_exons_involved=bool(data.get("other_exons_involved", False)),
            overlaps_benign_region=bool(data.get("overlaps_benign_region", False)),
            reading_frame_effect=reading_frame_effect,
            nmd_predicted=optional_bool(data, "nmd_predicted", ctx),
            exon_description=optional_str(data, "exon_description"),
            notes=optional_str(data, "notes"),
        )
