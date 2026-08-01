"""CnvDuplicationEvidence -- curated evidence for a single-gene copy-number
duplication, scoped to DMD only, batch 24. Parallel to
CnvDeletionEvidence (batch 23), but NOT a variant of it -- the ACMG/ClinGen
gain-side Section 2 rubric asks fundamentally different questions than the
loss side, per this batch's research:

Whole-gene triplosensitivity (TS) scoring is deliberately NOT implemented
here. ClinGen's own DMD-specific dosage curation states "whole gene
duplications have not been reported in association with clinical
phenotypes" for DMD -- the real, clinically relevant DMD duplication
mechanism is not a "triple dose" triplosensitivity effect at all. Instead:
"intragenic DMD duplications and triplications have been reported in
patients with DMD and BMD, presumably by a loss-of-function-type
mechanism" (ClinGen Dosage Sensitivity Curation, DMD). A tandem duplication
with a breakpoint inside the gene can disrupt it exactly the way a
deletion can -- via the same Aartsma-Rus reading-frame rule already used
for `CnvDeletionEvidence`, not via a separate triplosensitivity pathway.

This is why `whole_gene_duplicated` exists here as a representable but
NOT_MET-by-design input (see cnv_scoring.py): the field lets a curator
state a real whole-gene duplication if one is ever found, but this
project has no verified triplosensitivity scoring to apply to it, and no
real DMD fixture would exercise it honestly anyway (per the ClinGen
finding above) -- an explicit, disclosed gap, not an oversight.

Tandem confirmation is a hard prerequisite for any functional call,
sourced from real CNV-interpretation literature found during this batch's
research: a breakpoint study of 119 gain CNVs found 83% were tandem and
direct, "with the majority of the remainder interpreted as VUS because
the effect could not be determined." `is_tandem` is therefore required
whenever a breakpoint falls inside the gene, and `reading_frame_effect`
is only asked for (and only matters) once tandem orientation is
CONFIRMED -- an unconfirmed or non-tandem/complex insertion is treated as
functionally unpredictable regardless of any frame prediction a curator
might otherwise guess at.

Point-value disclosure, carried over from this batch's user-facing
scoping discussion: a 2019/2020 inter-laboratory CNV-concordance study
(PMC8960312) discusses real disagreement over "the use of 2K (0.45
points) or 2J (0 point) when a copy number gain breakpoint was observed
for the established HI genes" -- confirming that a real Riggs et al. 2020
gain-side category exists for exactly this scenario, but NOT which
condition (out-of-frame vs in-frame/unknown) maps to which of the two
point values, nor the real letter code (almost certainly NOT "2A"/"2C" as
ClassifyCNV's own internal dict-key reuse for the gain side might
suggest -- see CNV_GAIN_CATEGORY_CODES in enums.py). This project infers
the higher value (0.45) for the out-of-frame/disruptive case and the
lower value (0) for the in-frame/uncertain case, by direct analogy to
every other pathogenic-vs-uncertain split already used on the loss side
(2E's frameshift+NMD vs NONE_APPLICABLE's in-frame). That inference is
NOT independently confirmed against the primary paper -- disclosed here,
in cnv_scoring.py, and in every golden case that exercises it, per the
user's explicit choice to proceed on this basis (batch 24 scoping
discussion) rather than block on further primary-source access.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_bool, require_dict, require_str
from .enums import CnvDuplicationOrientation, CnvReadingFrameEffect, GenomeBuild


@dataclass(frozen=True)
class CnvDuplicationEvidence:
    cnv_id: str
    gene: str
    genome_build: GenomeBuild
    whole_gene_duplicated: bool
    breakpoint_within_gene: bool
    chromosome: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    coordinate_verified: bool = True
    overlaps_benign_region: bool = False
    is_tandem: Optional[CnvDuplicationOrientation] = None
    reading_frame_effect: Optional[CnvReadingFrameEffect] = None
    exon_description: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"CnvDuplicationEvidence[{self.cnv_id}]"
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
                "(same convention as VariantIdentity/CnvDeletionEvidence)"
            )
        if self.whole_gene_duplicated and self.overlaps_benign_region:
            raise SchemaValidationError(
                f"{context}: whole_gene_duplicated=True and overlaps_benign_region=True are contradictory"
            )
        if self.whole_gene_duplicated and self.breakpoint_within_gene:
            raise SchemaValidationError(
                f"{context}: whole_gene_duplicated=True and breakpoint_within_gene=True are contradictory -- "
                "a whole-gene duplication's breakpoints are, by definition, outside the gene"
            )
        needs_tandem_call = (
            self.breakpoint_within_gene
            and not self.whole_gene_duplicated
            and not self.overlaps_benign_region
        )
        if needs_tandem_call and self.is_tandem is None:
            raise SchemaValidationError(
                f"{context}: is_tandem must be explicitly stated when a breakpoint falls inside the gene -- "
                "no functional call (frameshift or otherwise) can be safely made from an unconfirmed "
                "tandem/direct orientation. Never left unstated, same 'never silently guess' convention "
                "as CnvDeletionEvidence.reading_frame_effect."
            )
        if self.is_tandem == CnvDuplicationOrientation.TANDEM and self.reading_frame_effect is None:
            raise SchemaValidationError(
                f"{context}: reading_frame_effect must be explicitly stated once is_tandem=TANDEM is "
                "confirmed -- this is exactly the scenario category GAIN_2K_EQUIV/GAIN_2J_EQUIV distinguish."
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "CnvDuplicationEvidence":
        data = require_dict(data, context or "CnvDuplicationEvidence")
        cnv_id = require_str(data, "cnv_id", context or "CnvDuplicationEvidence")
        ctx = f"CnvDuplicationEvidence[{cnv_id}]"
        gene = require_str(data, "gene", ctx)
        genome_build = coerce_enum(GenomeBuild, data.get("genome_build"), "genome_build", ctx)
        whole_gene_duplicated = require_bool(data, "whole_gene_duplicated", ctx)
        breakpoint_within_gene = require_bool(data, "breakpoint_within_gene", ctx)
        is_tandem_raw = data.get("is_tandem")
        is_tandem = (
            coerce_enum(CnvDuplicationOrientation, is_tandem_raw, "is_tandem", ctx)
            if is_tandem_raw is not None
            else None
        )
        reading_frame_raw = data.get("reading_frame_effect")
        reading_frame_effect = (
            coerce_enum(CnvReadingFrameEffect, reading_frame_raw, "reading_frame_effect", ctx)
            if reading_frame_raw is not None
            else None
        )
        start = data.get("start")
        if start is not None and (not isinstance(start, int) or isinstance(start, bool)):
            raise SchemaValidationError(f"{ctx}: 'start' must be an integer if provided")
        end = data.get("end")
        if end is not None and (not isinstance(end, int) or isinstance(end, bool)):
            raise SchemaValidationError(f"{ctx}: 'end' must be an integer if provided")
        return cls(
            cnv_id=cnv_id,
            gene=gene,
            genome_build=genome_build,
            whole_gene_duplicated=whole_gene_duplicated,
            breakpoint_within_gene=breakpoint_within_gene,
            chromosome=optional_str(data, "chromosome"),
            start=start,
            end=end,
            coordinate_verified=bool(data.get("coordinate_verified", True)),
            overlaps_benign_region=bool(data.get("overlaps_benign_region", False)),
            is_tandem=is_tandem,
            reading_frame_effect=reading_frame_effect,
            exon_description=optional_str(data, "exon_description"),
            notes=optional_str(data, "notes"),
        )
