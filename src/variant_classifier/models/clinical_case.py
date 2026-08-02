"""ClinicalCase — Milestone 4's new layer: what was found in one patient,
as opposed to what's true about one variant in isolation.

Everything through Milestone 3 (VariantEvidenceBundle, CriterionResult,
ProvisionalClassification) describes a single variant on its own. Real
clinical interpretation for a recessive or X-linked disease can't be
decided from one variant alone — it depends on how many variants this
specific patient has, and how they relate to each other. ClinicalCase is
deliberately minimal: just enough to answer that question for the two
inheritance patterns Milestone 4 covers (autosomal recessive, X-linked),
not a general patient/case record.

Scope, stated plainly: at most two variant_ids (a third variant, or
compound scenarios beyond simple biallelic recessive, hemizygous X-linked,
or biallelic XX X-linked (batch 29 — see clinical.py's
_interpret_xx_biallelic), are out of scope). phase is required when there
are two variants and forbidden when there's one, so a case can never be
silently ambiguous about whether phase was even considered. This model
itself stays gene/inheritance-agnostic on purpose — it is clinical.py's
interpret_x_linked_case, not this dataclass, that decides which
variant_ids-count/karyotypic_sex combinations are actually meaningful for
a given inheritance pattern (e.g. rejecting two variant_ids for a
hemizygous XY case).
"""

from dataclasses import dataclass, field
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_dict, require_list, require_str
from .enums import KaryotypicSex, PhaseRelationship


@dataclass(frozen=True)
class ClinicalCase:
    case_id: str
    gene: str
    karyotypic_sex: KaryotypicSex
    variant_ids: List[str] = field(default_factory=list)
    phase: Optional[PhaseRelationship] = None
    phase_evidence_note: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"ClinicalCase[{self.case_id}]"
        if not self.variant_ids:
            raise SchemaValidationError(f"{context}: variant_ids must not be empty")
        if len(self.variant_ids) > 2:
            raise SchemaValidationError(
                f"{context}: {len(self.variant_ids)} variant_ids given; Milestone 4 only supports "
                "one or two variants per case (simple biallelic recessive, hemizygous X-linked, or "
                "biallelic XX X-linked)"
            )
        if len(set(self.variant_ids)) != len(self.variant_ids):
            raise SchemaValidationError(f"{context}: variant_ids contains duplicates")
        if len(self.variant_ids) == 2 and self.phase is None:
            raise SchemaValidationError(
                f"{context}: two variant_ids given but phase is not set — phase must be explicitly "
                "TRANS, CIS, or UNKNOWN, never silently omitted"
            )
        if len(self.variant_ids) == 1 and self.phase is not None:
            raise SchemaValidationError(f"{context}: phase is meaningless for a single-variant case")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "ClinicalCase":
        data = require_dict(data, context or "ClinicalCase")
        case_id = require_str(data, "case_id", context or "ClinicalCase")
        ctx = f"ClinicalCase[{case_id}]"
        gene = require_str(data, "gene", ctx)
        karyotypic_sex = coerce_enum(KaryotypicSex, data.get("karyotypic_sex"), "karyotypic_sex", ctx)
        raw_variant_ids = require_list(data, "variant_ids", ctx)
        for vid in raw_variant_ids:
            if not isinstance(vid, str) or not vid.strip():
                raise SchemaValidationError(f"{ctx}: every entry in variant_ids must be a non-empty string")
        phase_raw = data.get("phase")
        phase = coerce_enum(PhaseRelationship, phase_raw, "phase", ctx) if phase_raw is not None else None
        return cls(
            case_id=case_id,
            gene=gene,
            karyotypic_sex=karyotypic_sex,
            variant_ids=list(raw_variant_ids),
            phase=phase,
            phase_evidence_note=optional_str(data, "phase_evidence_note"),
        )
