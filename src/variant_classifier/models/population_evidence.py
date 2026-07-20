"""PopulationEvidence — frequency data for one variant from one source
(gnomAD for Milestone 1). The retrieval_status field is the whole point of
this model: "gnomAD result unavailable" must never collapse into "variant
absent from gnomAD" (see ACMG Engine Detailed Design Guide, Table 13,
"Missing-data safeguard", and the Reporting and Dashboard Design Guide,
Section 4.2).
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_float, optional_int, optional_str, require_dict, require_str
from .enums import PopulationRetrievalStatus


@dataclass(frozen=True)
class PopulationEvidence:
    source: str
    source_version: str
    retrieval_status: PopulationRetrievalStatus
    overall_af: Optional[float] = None
    ancestry_specific_max_af: Optional[float] = None
    allele_count: Optional[int] = None
    allele_number: Optional[int] = None
    homozygote_count: Optional[int] = None
    locus_coverage_adequate: Optional[bool] = None

    def __post_init__(self) -> None:
        context = f"PopulationEvidence[{self.source} {self.source_version}]"
        frequency_fields_present = any(
            v is not None for v in (self.overall_af, self.allele_count, self.allele_number, self.homozygote_count)
        )
        if self.retrieval_status == PopulationRetrievalStatus.OBSERVED and self.overall_af is None:
            raise SchemaValidationError(
                f"{context}: retrieval_status=OBSERVED requires overall_af to be populated — "
                "OBSERVED means the source was queried successfully and returned this variant"
            )
        if self.retrieval_status == PopulationRetrievalStatus.ABSENT and self.locus_coverage_adequate is not True:
            raise SchemaValidationError(
                f"{context}: retrieval_status=ABSENT requires locus_coverage_adequate=True — "
                "a variant cannot be recorded as absent unless the locus was actually assessed "
                "(see the Milestone-1 gnomAD-outage synthetic case for the alternative: UNAVAILABLE)"
            )
        if self.retrieval_status in (
            PopulationRetrievalStatus.UNAVAILABLE,
            PopulationRetrievalStatus.NOT_ASSESSED,
            PopulationRetrievalStatus.UNKNOWN,
        ) and frequency_fields_present:
            raise SchemaValidationError(
                f"{context}: retrieval_status={self.retrieval_status.value} but frequency fields are populated — "
                "a source that was not successfully queried cannot also report numbers from it"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "PopulationEvidence":
        data = require_dict(data, context or "PopulationEvidence")
        source = require_str(data, "source", context or "PopulationEvidence")
        source_version = require_str(data, "source_version", f"PopulationEvidence[{source}]")
        ctx = f"PopulationEvidence[{source} {source_version}]"
        retrieval_status = coerce_enum(
            PopulationRetrievalStatus, data.get("retrieval_status"), "retrieval_status", ctx
        )
        coverage = data.get("locus_coverage_adequate")
        if coverage is not None and not isinstance(coverage, bool):
            raise SchemaValidationError(f"{ctx}: locus_coverage_adequate must be true/false if provided")
        return cls(
            source=source,
            source_version=source_version,
            retrieval_status=retrieval_status,
            overall_af=optional_float(data, "overall_af", ctx, minimum=0.0, maximum=1.0),
            ancestry_specific_max_af=optional_float(data, "ancestry_specific_max_af", ctx, minimum=0.0, maximum=1.0),
            allele_count=optional_int(data, "allele_count", ctx, minimum=0),
            allele_number=optional_int(data, "allele_number", ctx, minimum=0),
            homozygote_count=optional_int(data, "homozygote_count", ctx, minimum=0),
            locus_coverage_adequate=coverage,
        )
