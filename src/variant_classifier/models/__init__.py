"""Typed data models for the CAPN3/DMD variant classification prototype.

Seven models, matching the schemas defined across the ACMG Engine Detailed
Design Guide and the Workflow Architecture Guide:

- VariantIdentity            (variant_identity.py)
- GeneDiseaseContext          (gene_disease_context.py, + Specification)
- TranscriptConsequence       (transcript_consequence.py)
- PopulationEvidence          (population_evidence.py)
- CriterionResult             (criterion_result.py)
- ProvisionalClassification   (provisional_classification.py)

Plus one supporting evidence-domain model added while wiring up fixtures
(PP3/BP4 need something to evaluate from):

- ComputationalEvidence       (computational_evidence.py)

Plus two containers used only by this prototype's fixture/golden-case
loading, not part of the architecture guides themselves:

- VariantEvidenceBundle       (evidence_bundle.py)
- GoldenCase                  (golden_case.py)
"""

from .computational_evidence import ComputationalEvidence
from .criterion_result import CriterionResult
from .evidence_bundle import VariantEvidenceBundle
from .gene_disease_context import GeneDiseaseContext, Specification
from .golden_case import GoldenCase
from .population_evidence import PopulationEvidence
from .provisional_classification import ProvisionalClassification
from .transcript_consequence import TranscriptConsequence
from .variant_identity import VariantIdentity

__all__ = [
    "ComputationalEvidence",
    "CriterionResult",
    "GeneDiseaseContext",
    "GoldenCase",
    "PopulationEvidence",
    "ProvisionalClassification",
    "Specification",
    "TranscriptConsequence",
    "VariantEvidenceBundle",
    "VariantIdentity",
]
