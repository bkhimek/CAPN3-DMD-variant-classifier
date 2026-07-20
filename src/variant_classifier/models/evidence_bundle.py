"""VariantEvidenceBundle — the single-source-of-truth container that
gathers everything an evaluator needs for one variant in one gene-disease
context: identity, gene/disease context, transcript consequences, and
population evidence. Mirrors the CaseRecord single-source-of-truth pattern
from the Workflow Architecture Guide, scoped down to variant-level (not
case/patient-level) for Milestone 1.
"""

from dataclasses import dataclass
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import require_dict, require_list, require_str
from .computational_evidence import ComputationalEvidence
from .gene_disease_context import GeneDiseaseContext
from .population_evidence import PopulationEvidence
from .transcript_consequence import TranscriptConsequence
from .variant_identity import VariantIdentity


@dataclass(frozen=True)
class VariantEvidenceBundle:
    variant: VariantIdentity
    gene_disease_context: GeneDiseaseContext
    transcript_consequences: List[TranscriptConsequence]
    population_evidence: List[PopulationEvidence]
    computational_evidence: Optional[ComputationalEvidence] = None
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"VariantEvidenceBundle[{self.variant.variant_id}]"
        if self.variant.gene != self.gene_disease_context.gene:
            raise SchemaValidationError(
                f"{context}: variant.gene ({self.variant.gene!r}) does not match "
                f"gene_disease_context.gene ({self.gene_disease_context.gene!r})"
            )
        if not self.transcript_consequences:
            raise SchemaValidationError(f"{context}: at least one transcript_consequences entry is required")
        relevant = [tc for tc in self.transcript_consequences if tc.clinically_relevant]
        if len(relevant) != 1:
            raise SchemaValidationError(
                f"{context}: exactly one transcript_consequences entry must have "
                f"clinically_relevant=True (found {len(relevant)}) — PVS1 in particular "
                "requires an unambiguous clinically relevant transcript"
            )
        if not self.population_evidence:
            raise SchemaValidationError(f"{context}: at least one population_evidence entry is required")

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "VariantEvidenceBundle":
        data = require_dict(data, context or "VariantEvidenceBundle")
        variant = VariantIdentity.from_dict(data.get("variant", {}), context)
        ctx = f"VariantEvidenceBundle[{variant.variant_id}]"
        gene_disease_context = GeneDiseaseContext.from_dict(data.get("gene_disease_context", {}))
        raw_transcripts = require_list(data, "transcript_consequences", ctx)
        transcript_consequences = [
            TranscriptConsequence.from_dict(t, f"{ctx}.transcript_consequences[{i}]")
            for i, t in enumerate(raw_transcripts)
        ]
        raw_population = require_list(data, "population_evidence", ctx)
        population_evidence = [
            PopulationEvidence.from_dict(p, f"{ctx}.population_evidence[{i}]") for i, p in enumerate(raw_population)
        ]
        raw_computational = data.get("computational_evidence")
        computational_evidence = (
            ComputationalEvidence.from_dict(raw_computational, f"{ctx}.computational_evidence")
            if raw_computational is not None
            else None
        )
        notes = data.get("notes")
        if notes is not None and not isinstance(notes, str):
            raise SchemaValidationError(f"{ctx}: notes must be a string if provided")
        return cls(
            variant=variant,
            gene_disease_context=gene_disease_context,
            transcript_consequences=transcript_consequences,
            population_evidence=population_evidence,
            computational_evidence=computational_evidence,
            notes=notes,
        )
