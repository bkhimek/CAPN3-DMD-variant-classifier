"""VariantIdentity — what the variant is, independent of any transcript,
disease, or evidence. See ACMG Engine Detailed Design Guide, Section 4,
"Variant identity" evidence domain.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_int, optional_str, require_dict, require_str
from .enums import GenomeBuild


@dataclass(frozen=True)
class VariantIdentity:
    variant_id: str
    gene: str
    genome_build: GenomeBuild
    chromosome: Optional[str] = None
    position: Optional[int] = None
    reference: Optional[str] = None
    alternate: Optional[str] = None
    clinvar_variation_id: Optional[int] = None
    coordinate_verified: bool = True

    def __post_init__(self) -> None:
        context = f"VariantIdentity[{self.variant_id}]"
        if not isinstance(self.genome_build, GenomeBuild):
            raise SchemaValidationError(f"{context}: genome_build must be a GenomeBuild value")
        if self.position is not None and self.position <= 0:
            raise SchemaValidationError(f"{context}: position must be a positive integer if provided")
        if self.clinvar_variation_id is not None and self.clinvar_variation_id <= 0:
            raise SchemaValidationError(f"{context}: clinvar_variation_id must be a positive integer if provided")
        if not self.coordinate_verified and (self.chromosome or self.position or self.reference or self.alternate):
            # A record can legitimately omit coordinates entirely (HGVS-only,
            # e.g. this project's real ClinVar case, see README). It should
            # not claim unverified coordinates as if they were verified.
            raise SchemaValidationError(
                f"{context}: coordinate_verified=False but chromosome/position/reference/alternate "
                "were provided — either verify them and set coordinate_verified=True, or omit them"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "VariantIdentity":
        data = require_dict(data, context or "VariantIdentity")
        variant_id = require_str(data, "variant_id", context or "VariantIdentity")
        ctx = f"VariantIdentity[{variant_id}]"
        gene = require_str(data, "gene", ctx)
        genome_build = coerce_enum(GenomeBuild, data.get("genome_build"), "genome_build", ctx)
        return cls(
            variant_id=variant_id,
            gene=gene,
            genome_build=genome_build,
            chromosome=optional_str(data, "chromosome"),
            position=optional_int(data, "position", ctx, minimum=1),
            reference=optional_str(data, "reference"),
            alternate=optional_str(data, "alternate"),
            clinvar_variation_id=optional_int(data, "clinvar_variation_id", ctx, minimum=1),
            coordinate_verified=bool(data.get("coordinate_verified", True)),
        )
