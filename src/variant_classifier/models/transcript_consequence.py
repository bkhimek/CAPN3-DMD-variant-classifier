"""TranscriptConsequence — what a variant does at a specific transcript.
Kept separate from VariantIdentity because the same genomic change can be
described against multiple transcripts, and PVS1 in particular depends on
which transcript is clinically relevant (see ACMG Engine Detailed Design
Guide, Section 8, "PVS1 design in depth").
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_bool, require_dict, require_str
from .enums import Consequence


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

    def __post_init__(self) -> None:
        context = f"TranscriptConsequence[{self.transcript_id}]"
        if self.consequence == Consequence.FRAMESHIFT_VARIANT and self.nmd_predicted is None:
            raise SchemaValidationError(
                f"{context}: nmd_predicted must be explicitly true or false for a frameshift_variant — "
                "PVS1 cannot be evaluated safely from an unstated NMD prediction"
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
        )
