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

Plus one more evidence-domain model added for PS1/PM5 (same-residue
precedent evidence):

- SameResidueEvidence          (same_residue_evidence.py)

Plus two containers used only by this prototype's fixture/golden-case
loading, not part of the architecture guides themselves:

- VariantEvidenceBundle       (evidence_bundle.py)
- GoldenCase                  (golden_case.py)

Plus two Milestone 4 models for case-level (not variant-level) reasoning:

- ClinicalCase                (clinical_case.py)
- CaseInterpretation           (case_interpretation.py)

Plus three batch-23 models for DMD CNV/deletion scoring -- a deliberately
separate evidence family, not an extension of VariantEvidenceBundle (see
cnv_deletion_evidence.py's docstring for why):

- CnvDeletionEvidence          (cnv_deletion_evidence.py)
- CnvCategoryResult            (cnv_category_result.py)
- CnvProvisionalClassification (cnv_provisional_classification.py)

Plus one batch-24 model for DMD CNV/duplication scoring -- parallel to
CnvDeletionEvidence, sharing CnvCategoryResult/CnvProvisionalClassification
rather than forking them (see cnv_duplication_evidence.py's docstring):

- CnvDuplicationEvidence       (cnv_duplication_evidence.py)

Plus one batch-25 model for PS3/BS3 functional-assay evidence, attached
to VariantEvidenceBundle alongside ComputationalEvidence and
SameResidueEvidence (see functional_evidence.py's docstring):

- FunctionalEvidence           (functional_evidence.py)

Plus two batch-28 models for PM3 (in-trans-with-a-pathogenic-variant)
evidence, attached to VariantEvidenceBundle alongside FunctionalEvidence
-- the case-level circularity PM3 has carried since Milestone 4 is
resolved by curating the partner allele's classification as a known fact
rather than re-deriving it from this engine (see pm3_evidence.py's
docstring for the full writeup):

- Pm3Evidence                  (pm3_evidence.py)
- Pm3ProbandObservation        (pm3_evidence.py)
"""

from .case_interpretation import CaseInterpretation
from .clinical_case import ClinicalCase
from .cnv_category_result import CnvCategoryResult
from .cnv_deletion_evidence import CnvDeletionEvidence
from .cnv_duplication_evidence import CnvDuplicationEvidence
from .cnv_provisional_classification import CnvProvisionalClassification
from .computational_evidence import ComputationalEvidence
from .criterion_result import CriterionResult
from .evidence_bundle import VariantEvidenceBundle
from .functional_evidence import FunctionalEvidence
from .gene_disease_context import GeneDiseaseContext, Specification
from .golden_case import GoldenCase
from .pm3_evidence import Pm3Evidence, Pm3ProbandObservation
from .population_evidence import PopulationEvidence
from .provisional_classification import ProvisionalClassification
from .same_residue_evidence import SameResidueEvidence
from .transcript_consequence import TranscriptConsequence
from .variant_identity import VariantIdentity

__all__ = [
    "CaseInterpretation",
    "ClinicalCase",
    "CnvCategoryResult",
    "CnvDeletionEvidence",
    "CnvDuplicationEvidence",
    "CnvProvisionalClassification",
    "ComputationalEvidence",
    "CriterionResult",
    "FunctionalEvidence",
    "GeneDiseaseContext",
    "GoldenCase",
    "Pm3Evidence",
    "Pm3ProbandObservation",
    "PopulationEvidence",
    "ProvisionalClassification",
    "SameResidueEvidence",
    "Specification",
    "TranscriptConsequence",
    "VariantEvidenceBundle",
    "VariantIdentity",
]
