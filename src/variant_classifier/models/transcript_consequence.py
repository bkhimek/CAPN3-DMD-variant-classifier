"""TranscriptConsequence — what a variant does at a specific transcript.
Kept separate from VariantIdentity because the same genomic change can be
described against multiple transcripts, and PVS1 in particular depends on
which transcript is clinically relevant (see ACMG Engine Detailed Design
Guide, Section 8, "PVS1 design in depth").

Batch 27 adds three more optional, START_LOST-only fields implementing
the real alternative-start-codon rule for PVS1's initiation-codon branch
(ACGS 2024 UK Practice Guidelines, quoting Abou Tayoun et al. 2018's
decision tree): "If there is a potential in-frame initiation codon
downstream, the missing N-terminal region of the protein should be
assessed according to the principles described in the decision tree
(i.e. is the missing region critical to protein function / is it >10%
of the entire protein length / are there any reported pathogenic
variants upstream of the potential initiation codon) and apply PVS1 at
either reduced strength or n/a, as appropriate. If no alternative
in-frame start codon is identified, use PVS1 at maximum strength." See
evaluators/pvs1.py's module docstring for exactly how these three facts
combine into a decision, and README's "PVS1 start-loss: the alternative-
start-codon rule (batch 27)" design note for the full writeup, including
how CAPN3_c.1A>G's real values were independently derived from the
primary RefSeq/Ensembl CDS sequence (not the blocked Abou Tayoun/Walker
PDFs) plus this project's own curated fixture set.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import (
    coerce_enum,
    optional_bool,
    optional_float,
    optional_int,
    optional_str,
    require_bool,
    require_dict,
    require_str,
)
from .enums import Consequence, SplicingRnaEvidence


@dataclass(frozen=True)
class TranscriptConsequence:
    transcript_id: str
    clinically_relevant: bool
    consequence: Consequence
    mane_select: bool = False
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    exon: Optional[str] = None
    nmd_predicted: Optional[bool] = None
    repeat_region: Optional[bool] = None
    protein_length_change_aa: Optional[int] = None
    splicing_rna_evidence: Optional[SplicingRnaEvidence] = None
    alternative_start_codon_identified: Optional[bool] = None
    alternative_start_codon_percent_protein_lost: Optional[float] = None
    alternative_start_codon_preceded_by_pathogenic_variant: Optional[bool] = None

    NMD_RELEVANT_CONSEQUENCES = (Consequence.FRAMESHIFT_VARIANT, Consequence.STOP_GAINED)
    PM4_RELEVANT_CONSEQUENCES = (Consequence.INFRAME_DELETION, Consequence.INFRAME_INSERTION, Consequence.STOP_LOST)
    SPLICE_RELEVANT_CONSEQUENCES = (Consequence.SPLICE_DONOR_VARIANT, Consequence.SPLICE_ACCEPTOR_VARIANT)

    def __post_init__(self) -> None:
        context = f"TranscriptConsequence[{self.transcript_id}]"
        if self.consequence in TranscriptConsequence.NMD_RELEVANT_CONSEQUENCES and self.nmd_predicted is None:
            raise SchemaValidationError(
                f"{context}: nmd_predicted must be explicitly true or false for a "
                f"{self.consequence.value} — PVS1 cannot be evaluated safely from an unstated "
                "NMD prediction. (Originally only enforced for frameshift_variant; widened to "
                "stop_gained too during the PVS1 evaluator build, since nonsense variants are "
                "subject to the same NMD-vs-last-exon logic.)"
            )
        if self.consequence in TranscriptConsequence.PM4_RELEVANT_CONSEQUENCES and self.repeat_region is None:
            raise SchemaValidationError(
                f"{context}: repeat_region must be explicitly true or false for a "
                f"{self.consequence.value} — PM4 excludes repeat-region indels/stop-losses by "
                "definition (Richards et al. 2015), so this cannot be evaluated safely from an "
                "unstated repeat-region status. Same 'never silently guess' convention as "
                "nmd_predicted for frameshift/stop_gained."
            )
        if (
            self.consequence not in TranscriptConsequence.SPLICE_RELEVANT_CONSEQUENCES
            and self.splicing_rna_evidence is not None
        ):
            raise SchemaValidationError(
                f"{context}: splicing_rna_evidence is only meaningful for "
                f"{'/'.join(c.value for c in TranscriptConsequence.SPLICE_RELEVANT_CONSEQUENCES)}, "
                f"not {self.consequence.value} — batch 26's PVS1 splice-RNA-evidence branch does not "
                "apply to this consequence class."
            )
        start_loss_fields_set = (
            self.alternative_start_codon_identified is not None
            or self.alternative_start_codon_percent_protein_lost is not None
            or self.alternative_start_codon_preceded_by_pathogenic_variant is not None
        )
        if self.consequence != Consequence.START_LOST and start_loss_fields_set:
            raise SchemaValidationError(
                f"{context}: alternative_start_codon_* fields are only meaningful for "
                f"{Consequence.START_LOST.value}, not {self.consequence.value} — batch 27's PVS1 "
                "alternative-start-codon branch does not apply to this consequence class."
            )
        if self.alternative_start_codon_identified is True:
            if self.alternative_start_codon_percent_protein_lost is None:
                raise SchemaValidationError(
                    f"{context}: alternative_start_codon_percent_protein_lost must be stated whenever "
                    "alternative_start_codon_identified=True — the ACGS 2024 / Abou Tayoun et al. 2018 "
                    "rule requires knowing what fraction of the protein the missing N-terminal region "
                    "represents to decide between a reduced-strength PVS1 call and a full "
                    "protein-domain-criticality assessment. Never left unstated, same 'never silently "
                    "guess' convention as nmd_predicted/repeat_region/splicing_rna_evidence."
                )
            if self.alternative_start_codon_preceded_by_pathogenic_variant is None:
                raise SchemaValidationError(
                    f"{context}: alternative_start_codon_preceded_by_pathogenic_variant must be stated "
                    "whenever alternative_start_codon_identified=True — whether a known pathogenic "
                    "variant falls within the region that would be lost is one of the three explicit "
                    "factors the real decision tree assesses, and cannot be safely left unstated."
                )
        elif self.alternative_start_codon_identified in (False, None):
            if (
                self.alternative_start_codon_percent_protein_lost is not None
                or self.alternative_start_codon_preceded_by_pathogenic_variant is not None
            ):
                raise SchemaValidationError(
                    f"{context}: alternative_start_codon_percent_protein_lost and "
                    "alternative_start_codon_preceded_by_pathogenic_variant are only meaningful when "
                    "alternative_start_codon_identified=True — neither question applies if no "
                    "alternative start codon exists (or its existence hasn't been assessed)."
                )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "TranscriptConsequence":
        data = require_dict(data, context or "TranscriptConsequence")
        transcript_id = require_str(data, "transcript_id", context or "TranscriptConsequence")
        ctx = f"TranscriptConsequence[{transcript_id}]"
        clinically_relevant = require_bool(data, "clinically_relevant", ctx)
        consequence = coerce_enum(Consequence, data.get("consequence"), "consequence", ctx)
        nmd_value = data.get("nmd_predicted")
        if nmd_value is not None and not isinstance(nmd_value, bool):
            raise SchemaValidationError(f"{ctx}: nmd_predicted must be true/false if provided")
        return cls(
            transcript_id=transcript_id,
            clinically_relevant=clinically_relevant,
            consequence=consequence,
            mane_select=bool(data.get("mane_select", False)),
            hgvs_c=optional_str(data, "hgvs_c"),
            hgvs_p=optional_str(data, "hgvs_p"),
            exon=optional_str(data, "exon"),
            nmd_predicted=nmd_value,
            repeat_region=optional_bool(data, "repeat_region", ctx),
            protein_length_change_aa=optional_int(data, "protein_length_change_aa", ctx, minimum=1),
            splicing_rna_evidence=(
                coerce_enum(SplicingRnaEvidence, data.get("splicing_rna_evidence"), "splicing_rna_evidence", ctx)
                if data.get("splicing_rna_evidence") is not None
                else None
            ),
            alternative_start_codon_identified=optional_bool(data, "alternative_start_codon_identified", ctx),
            alternative_start_codon_percent_protein_lost=optional_float(
                data, "alternative_start_codon_percent_protein_lost", ctx, minimum=0.0, maximum=100.0
            ),
            alternative_start_codon_preceded_by_pathogenic_variant=optional_bool(
                data, "alternative_start_codon_preceded_by_pathogenic_variant", ctx
            ),
        )
