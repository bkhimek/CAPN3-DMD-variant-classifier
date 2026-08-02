"""Milestone 4: case-level clinical interpretation.

Deliberately its own module, not folded into engine.py. engine.py answers
"how strong is the evidence for this one variant" from raw evidence
alone. clinical.py answers a different question — "does what was found in
this patient explain their disease" — and needs case-level information
(ClinicalCase) plus each involved variant's already-computed
ProvisionalClassification as input. It does not re-derive variant-level
evidence; it consumes engine.classify()'s output.

PM3 ("detected in trans with a pathogenic variant") was not bolted onto
the per-variant evaluator pattern from Milestone 2/3 either -- it has a
structural circularity a naive per-variant evaluator can't cleanly
express (variant A's PM3 would depend on variant B's classification,
which could itself depend on evidence entirely outside variant A's own
VariantEvidenceBundle). Batch 28's evaluators/pm3.py resolves this a
different way than this module does: rather than reasoning about it here
at the case level after both variants already have their own
classification, it curates the partner allele's classification as a
known fact directly on the variant's own bundle (models/pm3_evidence.py),
the same way every other "can't be safely re-derived by this engine"
value in this project is handled. This module (ClinicalCase /
CaseInterpretation) remains the place for phase/hemizygosity reasoning
that genuinely needs case-level information PM3 itself doesn't (does a
qualifying genotype in this specific patient explain their disease) --
the two are complementary, not overlapping.

Scope, stated plainly:
- Autosomal recessive: handles exactly one or two variants, with phase
  required whenever there are two (see ClinicalCase). Compound scenarios
  beyond that are out of scope.
- X-linked: only handles a single variant. A hemizygous male
  (karyotypic_sex=XY) with a convincing variant is the clean case this
  module resolves confidently. Any non-XY case (XX, OTHER, UNKNOWN) is
  deferred to MANUAL_REVIEW rather than reasoned about — female X-linked
  carrier interpretation involves X-inactivation biology this project
  does not model, and guessing would be worse than admitting the gap.
"""

from typing import Dict

from .errors import SchemaValidationError
from .models import CaseInterpretation, ClinicalCase, GeneDiseaseContext, ProvisionalClassification
from .models.enums import CaseInterpretationStatus, Inheritance, KaryotypicSex, PhaseRelationship, ProvisionalClass

_QUALIFYING = (ProvisionalClass.PATHOGENIC, ProvisionalClass.LIKELY_PATHOGENIC)
_BENIGN_SIDE = (ProvisionalClass.BENIGN, ProvisionalClass.LIKELY_BENIGN)


def interpret_recessive_case(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
    gene_disease_context: GeneDiseaseContext,
) -> CaseInterpretation:
    if gene_disease_context.inheritance != Inheritance.AUTOSOMAL_RECESSIVE:
        raise SchemaValidationError(
            f"interpret_recessive_case[{case.case_id}]: gene_disease_context.inheritance="
            f"{gene_disease_context.inheritance.value}, not AUTOSOMAL_RECESSIVE"
        )

    if len(case.variant_ids) == 1:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.INSUFFICIENT,
            rationale=(
                "Only one variant identified in this gene. Autosomal recessive disease requires two "
                "affected alleles (one from each parent); a single variant, however strong its own "
                "evidence, does not on its own explain a recessive disease presentation."
            ),
        )

    v1_id, v2_id = case.variant_ids
    c1, c2 = classifications[v1_id], classifications[v2_id]

    if case.phase == PhaseRelationship.CIS:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.INSUFFICIENT,
            rationale=(
                f"{v1_id} and {v2_id} are confirmed in cis (same chromosome copy) — one copy of the "
                "gene carries both variants, but the other copy is untouched and remains functional. "
                "This does not explain autosomal recessive disease regardless of how pathogenic either "
                "variant looks individually."
            ),
        )

    if case.phase == PhaseRelationship.UNKNOWN:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.MANUAL_REVIEW,
            rationale=(
                f"{v1_id} and {v2_id} were both identified, but their phase relationship (trans vs cis) "
                "has not been established — via parental testing or otherwise. Whether this patient "
                "has one broken copy of the gene or two cannot be determined without it."
            ),
        )

    # phase == TRANS
    if c1.provisional_class in _BENIGN_SIDE or c2.provisional_class in _BENIGN_SIDE:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.INSUFFICIENT,
            rationale=(
                f"{v1_id} and {v2_id} are confirmed in trans, but at least one is classified "
                f"({c1.provisional_class.value} / {c2.provisional_class.value}) as (Likely) Benign — "
                "it does not count as a disease-causing allele, so biallelic involvement is not "
                "established even though phase is known."
            ),
        )

    if c1.provisional_class in _QUALIFYING and c2.provisional_class in _QUALIFYING:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.EXPLAINED,
            rationale=(
                f"{v1_id} ({c1.provisional_class.value}) and {v2_id} ({c2.provisional_class.value}) are "
                "confirmed in trans, and both are individually classified Pathogenic or Likely "
                "Pathogenic — biallelic loss of function is consistent with this patient's autosomal "
                "recessive disease."
            ),
        )

    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"{v1_id} ({c1.provisional_class.value}) and {v2_id} ({c2.provisional_class.value}) are "
            "confirmed in trans, but at least one is not yet classified Pathogenic or Likely "
            "Pathogenic — biallelic involvement is plausible but not conclusively established."
        ),
    )


def interpret_x_linked_case(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
    gene_disease_context: GeneDiseaseContext,
) -> CaseInterpretation:
    if gene_disease_context.inheritance not in (Inheritance.X_LINKED_RECESSIVE, Inheritance.X_LINKED_DOMINANT):
        raise SchemaValidationError(
            f"interpret_x_linked_case[{case.case_id}]: gene_disease_context.inheritance="
            f"{gene_disease_context.inheritance.value}, not X-linked"
        )
    if len(case.variant_ids) != 1:
        raise SchemaValidationError(
            f"interpret_x_linked_case[{case.case_id}]: Milestone 4 X-linked interpretation only "
            f"supports a single variant per case, found {len(case.variant_ids)}"
        )

    variant_id = case.variant_ids[0]
    classification = classifications[variant_id]

    if case.karyotypic_sex != KaryotypicSex.XY:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.MANUAL_REVIEW,
            rationale=(
                f"karyotypic_sex={case.karyotypic_sex.value}, not XY. Non-hemizygous X-linked "
                "interpretation (carrier status, X-inactivation effects on expressivity) is out of "
                "scope for this evaluator — deferred to manual review rather than reasoned about "
                "incorrectly."
            ),
        )

    if classification.provisional_class in _BENIGN_SIDE:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.INSUFFICIENT,
            rationale=(
                f"{variant_id} is classified {classification.provisional_class.value} — does not "
                "explain this patient's disease regardless of hemizygosity."
            ),
        )

    if classification.provisional_class in _QUALIFYING:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.EXPLAINED,
            rationale=(
                f"karyotypic_sex=XY (hemizygous — a single copy of the X chromosome), and {variant_id} "
                f"is classified {classification.provisional_class.value}. A single variant is "
                "sufficient to explain X-linked disease in a hemizygous individual; no second allele "
                "is needed the way it would be for an autosomal recessive gene."
            ),
        )

    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"karyotypic_sex=XY, but {variant_id} is classified "
            f"{classification.provisional_class.value} — not yet conclusive enough to explain disease "
            "on its own."
        ),
    )


def interpret_case(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
    gene_disease_context: GeneDiseaseContext,
) -> CaseInterpretation:
    """Dispatches to interpret_recessive_case or interpret_x_linked_case based on
    gene_disease_context.inheritance. Any other inheritance pattern (autosomal
    dominant, mitochondrial, unknown) is NOT_APPLICABLE — Milestone 4 only
    covers the two patterns CAPN3 and DMD illustrate.
    """
    if gene_disease_context.inheritance == Inheritance.AUTOSOMAL_RECESSIVE:
        return interpret_recessive_case(case, classifications, gene_disease_context)
    if gene_disease_context.inheritance in (Inheritance.X_LINKED_RECESSIVE, Inheritance.X_LINKED_DOMINANT):
        return interpret_x_linked_case(case, classifications, gene_disease_context)
    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.NOT_APPLICABLE,
        rationale=(
            f"inheritance={gene_disease_context.inheritance.value} is not one of the patterns "
            "Milestone 4 covers (autosomal recessive, X-linked)."
        ),
    )
