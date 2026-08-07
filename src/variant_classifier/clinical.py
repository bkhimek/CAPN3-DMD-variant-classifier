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
- X-linked, karyotypic_sex=XY: handles exactly one variant. A hemizygous
  male with a convincing variant is the clean case this module resolves
  confidently; two variant_ids for XY is rejected outright (only one X
  chromosome exists, so it isn't a real genotype), not silently ignored.
- X-linked, karyotypic_sex=XX: handles one OR two variants (extended in
  batch 29 — see "X-linked female/other-karyotype case interpretation
  (batch 29)" in the README for the full design). A single variant is
  still MANUAL_REVIEW, unchanged since Milestone 4: heterozygous carrier
  phenotype depends on X-inactivation biology the field itself cannot
  reliably predict (Brioschi et al. 2012). Two variants confirmed in
  trans, both independently Pathogenic/Likely Pathogenic, resolve to
  EXPLAINED — a real, documented, X-inactivation-independent mechanism
  (both X copies affected, so no genuinely functional copy exists for
  inactivation to preferentially spare). Every other two-variant XX
  combination (cis, unknown phase, or trans-but-not-both-qualifying)
  stays MANUAL_REVIEW for the same real reason as the single-variant
  case — see _interpret_xx_biallelic's docstring for exactly why cis
  isn't treated the way autosomal recessive cis is.
- X-linked, karyotypic_sex=OTHER or UNKNOWN: only a single variant is
  modeled, always MANUAL_REVIEW — deferred rather than reasoned about,
  since OTHER lumps together karyotypes with genuinely different X-linked
  dosage biology (X0/Turner is functionally hemizygous like XY; XXY is
  diploid-X like XX; mosaicism is neither uniformly) that this project
  does not attempt to disambiguate.
- Autosomal dominant (batch 31 — see "BRCA1 extension (Batch 31)" in the
  README for the full design): handles exactly one variant. BRCA1 is this
  module's first risk-conferring (as opposed to deterministic-diagnosis)
  gene — a qualifying monoallelic variant does not explain a diagnosis the
  way AR/X-linked qualifying genotypes do, it confers elevated,
  penetrance-dependent risk. See interpret_dominant_case and
  CaseInterpretationStatus.RISK_CONFERRING. Two-variant AD case reasoning
  (e.g. a suspected biallelic BRCA1/BRCA2 Fanconi-anemia-phenotype
  presentation) is out of scope this batch, rejected outright rather than
  reasoned about incorrectly, the same treatment X-linked XY gives a
  two-variant case that isn't a real genotype.
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

    if case.karyotypic_sex == KaryotypicSex.XY:
        if len(case.variant_ids) != 1:
            raise SchemaValidationError(
                f"interpret_x_linked_case[{case.case_id}]: karyotypic_sex=XY is hemizygous — there is "
                f"only one X chromosome, so {len(case.variant_ids)} variant_ids for an X-linked gene "
                "does not correspond to a real genotype"
            )
        return _interpret_hemizygous_male(case, classifications)

    if case.karyotypic_sex == KaryotypicSex.XX:
        if len(case.variant_ids) == 1:
            return _interpret_xx_single_variant(case)
        return _interpret_xx_biallelic(case, classifications)

    # OTHER (XXY, X0/Turner, mosaicism, ...) and UNKNOWN: genuinely heterogeneous karyotypes with
    # different X-linked dosage biology each (X0 is hemizygous like XY; XXY is diploid-X like XX;
    # mosaicism is neither cleanly) that this project does not attempt to disambiguate -- batch 29
    # research (see "X-linked female/other-karyotype case interpretation (batch 29)" in the README)
    # found real literature describing X0/Turner DMD carriers as functionally hemizygous, but
    # OTHER does not distinguish X0 from XXY/mosaic, so guessing which sub-case applies would be
    # worse than admitting the gap. Batch 29 relaxed XX from one variant to one-or-two; OTHER/UNKNOWN
    # remain exactly as narrow as before that batch -- single variant only, always MANUAL_REVIEW.
    if len(case.variant_ids) != 1:
        raise SchemaValidationError(
            f"interpret_x_linked_case[{case.case_id}]: karyotypic_sex={case.karyotypic_sex.value} with "
            f"{len(case.variant_ids)} variant_ids is not modeled -- multi-variant X-linked reasoning is "
            "only implemented for XX (see _interpret_xx_biallelic), not for OTHER/UNKNOWN karyotypes, "
            "whose X-linked dosage biology varies too much by specific karyotype to reason about "
            "generically"
        )
    variant_id = case.variant_ids[0]
    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"karyotypic_sex={case.karyotypic_sex.value}, not XY or XX. This karyotype's X-linked "
            "dosage biology (e.g. X0/Turner is functionally hemizygous, XXY is diploid-X, mosaicism "
            "is neither uniformly) is out of scope for this evaluator — deferred to manual review "
            "rather than reasoned about incorrectly."
        ),
    )


def _interpret_hemizygous_male(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
) -> CaseInterpretation:
    variant_id = case.variant_ids[0]
    classification = classifications[variant_id]

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


def _interpret_xx_single_variant(case: ClinicalCase) -> CaseInterpretation:
    variant_id = case.variant_ids[0]
    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"karyotypic_sex=XX with a single variant ({variant_id}) — a heterozygous X-linked "
            "carrier. Unchanged since Milestone 4: real literature (Brioschi et al. 2012, BMC Med "
            "Genet 13:73, cited in this project's CASE_DMD_FEMALE_CARRIER_REAL fixture) found skewed "
            "X-inactivation in only 2 of 6 symptomatic DMD carriers and in 5 of 11 asymptomatic ones "
            "-- the field itself cannot reliably predict a heterozygous carrier's phenotype from "
            "genotype plus X-inactivation pattern alone, so this evaluator does not either."
        ),
    )


def _interpret_xx_biallelic(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
) -> CaseInterpretation:
    """Batch 29. karyotypic_sex=XX gives two copies of an X-linked gene, so a second, genuinely
    different scenario from the single-heterozygous-carrier case above becomes representable:
    both copies affected (homozygous or compound heterozygous), documented in real literature as a
    real, if less common, cause of a fully or near-fully manifesting DMD phenotype in females --
    see "X-linked female/other-karyotype case interpretation (batch 29)" in the README for the
    real citations (Ulm et al. 2022; Fujii et al. 2009; Takeshita et al. 2017) and why this is
    mechanistically distinct from (and does not contradict) the single-carrier XCI-unpredictability
    finding above: when BOTH X copies carry a qualifying variant, there is no genuinely functional
    copy anywhere for X-inactivation to preferentially silence or spare -- every cell's active X is
    a qualifying one, regardless of which X gets inactivated in that cell. X-inactivation mosaicism
    only creates unpredictability when one copy is truly wild-type (or benign), which is exactly why
    every other combination below (CIS, unknown phase, or trans-but-not-both-qualifying) still
    resolves to MANUAL_REVIEW rather than being upgraded the way autosomal recessive's equivalent
    branches are: a benign or wild-type second X copy is still subject to the same real,
    unpredictable inactivation pattern as the single-carrier case, even though this project's own
    ClinicalCase model can now represent two variant_ids for XX cases.
    """
    v1_id, v2_id = case.variant_ids
    c1, c2 = classifications[v1_id], classifications[v2_id]

    if case.phase == PhaseRelationship.TRANS and c1.provisional_class in _QUALIFYING and c2.provisional_class in _QUALIFYING:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.EXPLAINED,
            rationale=(
                f"karyotypic_sex=XX, {v1_id} ({c1.provisional_class.value}) and {v2_id} "
                f"({c2.provisional_class.value}) confirmed in trans -- both copies of this X-linked "
                "gene carry a qualifying variant, so no genuinely functional copy exists in any cell "
                "for X-inactivation to preferentially spare. This biallelic mechanism is real and "
                "documented in the literature (see this function's docstring), and does not depend "
                "on X-inactivation pattern the way a single-heterozygous-variant carrier does."
            ),
        )

    if case.phase == PhaseRelationship.CIS:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.MANUAL_REVIEW,
            rationale=(
                f"karyotypic_sex=XX, {v1_id} and {v2_id} confirmed in cis (same X chromosome copy) -- "
                "unlike the autosomal recessive case, this is NOT resolved to INSUFFICIENT: the other "
                "X chromosome copy is genuinely wild-type, but (unlike an autosome) it is not always "
                "active -- X-inactivation randomly silences one X per cell, so a fully wild-type "
                "second copy can still be silenced in some cells, carrying the same real, "
                "unpredictable-phenotype uncertainty as a single-heterozygous carrier."
            ),
        )

    if case.phase == PhaseRelationship.UNKNOWN:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.MANUAL_REVIEW,
            rationale=(
                f"karyotypic_sex=XX, {v1_id} and {v2_id} were both identified, but their phase "
                "relationship (trans vs cis) has not been established. Whether this patient has a "
                "genuinely functional X copy anywhere (which would reintroduce the same real "
                "X-inactivation unpredictability as the single-carrier case) cannot be determined "
                "without it."
            ),
        )

    # phase == TRANS but not both qualifying
    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"karyotypic_sex=XX, {v1_id} ({c1.provisional_class.value}) and {v2_id} "
            f"({c2.provisional_class.value}) are confirmed in trans, but at least one is not yet "
            "classified Pathogenic or Likely Pathogenic -- a non-qualifying variant on the other X "
            "copy is functionally similar to a wild-type copy for this purpose, so the same real "
            "X-inactivation unpredictability as the single-carrier case still applies rather than "
            "being resolved either way."
        ),
    )


def interpret_dominant_case(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
    gene_disease_context: GeneDiseaseContext,
) -> CaseInterpretation:
    """Batch 31 (BRCA1 extension). Autosomal dominant, risk-conferring case-level
    interpretation -- this module's first use of Inheritance.AUTOSOMAL_DOMINANT.

    Distinct in kind, not just in name, from interpret_recessive_case and
    interpret_x_linked_case: those resolve to EXPLAINED because a qualifying
    genotype deterministically accounts for a Mendelian diagnosis. A single
    monoallelic BRCA1 pathogenic variant does not deterministically cause
    disease -- it confers elevated, penetrance-dependent cancer risk.
    Conflating the two under EXPLAINED would misrepresent what this classifier
    is actually claiming, so this branch resolves to the distinct
    CaseInterpretationStatus.RISK_CONFERRING instead. See README.md, "BRCA1
    extension (Batch 31)" for the full design writeup.

    A Benign/Likely Benign classification here resolves to INSUFFICIENT, the
    same status interpret_recessive_case's own trans-phase branch already uses
    when a co-variant is (Likely) Benign ("it does not count as a
    disease-causing allele"). That's a deliberate, considered reuse, not an
    accidental one: INSUFFICIENT already spans both "not enough evidence
    quantity yet" (e.g. only one of two AR variants found) and "this specific
    variant is confirmed not to matter" (a refuted hypothesis, not a quantity
    problem) elsewhere in this module, so extending that same generalization
    to AD's single-variant case is consistent with existing precedent rather
    than a new kind of overload. A distinct status (e.g. REFUTED) would be
    equally defensible and was considered -- not introduced this batch to
    avoid growing CaseInterpretationStatus by two new values in one pass when
    only RISK_CONFERRING has a clear, load-bearing semantic justification.

    Scope, stated plainly, matching this module's existing convention:
    - Handles exactly one variant_id. A patient with two BRCA1 variant_ids
      (e.g. a suspected biallelic/Fanconi-anemia-phenotype presentation) is
      out of scope this batch -- deferred alongside PM3/BS2 (see README) --
      and raises SchemaValidationError rather than silently reasoning about
      it incorrectly, the same "reject outright, don't guess" treatment
      interpret_x_linked_case gives two variant_ids for XY.
    """
    if gene_disease_context.inheritance != Inheritance.AUTOSOMAL_DOMINANT:
        raise SchemaValidationError(
            f"interpret_dominant_case[{case.case_id}]: gene_disease_context.inheritance="
            f"{gene_disease_context.inheritance.value}, not AUTOSOMAL_DOMINANT"
        )

    if len(case.variant_ids) != 1:
        raise SchemaValidationError(
            f"interpret_dominant_case[{case.case_id}]: {len(case.variant_ids)} variant_ids given, but "
            "this batch only handles a single monoallelic variant. Two-variant biallelic/"
            "Fanconi-anemia-phenotype case reasoning for BRCA1 is deferred to a later batch (see "
            "README.md, 'BRCA1 extension (Batch 31)') -- rejected outright rather than reasoned "
            "about incorrectly."
        )

    variant_id = case.variant_ids[0]
    classification = classifications[variant_id]

    if classification.provisional_class in _BENIGN_SIDE:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.INSUFFICIENT,
            rationale=(
                f"{variant_id} is classified {classification.provisional_class.value} — does not "
                "confer elevated disease risk regardless of zygosity."
            ),
        )

    if classification.provisional_class in _QUALIFYING:
        return CaseInterpretation(
            case_id=case.case_id,
            gene=case.gene,
            status=CaseInterpretationStatus.RISK_CONFERRING,
            rationale=(
                f"inheritance=AUTOSOMAL_DOMINANT, and {variant_id} is classified "
                f"{classification.provisional_class.value}. A single monoallelic variant is "
                "sufficient to confer elevated, penetrance-dependent disease risk in an autosomal "
                "dominant, risk-conferring gene -- distinct from EXPLAINED, which this project "
                "reserves for genotypes that deterministically account for a Mendelian diagnosis "
                "(as in CAPN3's autosomal recessive or DMD's X-linked hemizygous mechanisms)."
            ),
        )

    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.MANUAL_REVIEW,
        rationale=(
            f"inheritance=AUTOSOMAL_DOMINANT, but {variant_id} is classified "
            f"{classification.provisional_class.value} — not yet conclusive enough to establish "
            "elevated risk on its own."
        ),
    )


def interpret_case(
    case: ClinicalCase,
    classifications: Dict[str, ProvisionalClassification],
    gene_disease_context: GeneDiseaseContext,
) -> CaseInterpretation:
    """Dispatches to interpret_recessive_case, interpret_x_linked_case, or
    interpret_dominant_case based on gene_disease_context.inheritance. Any
    other inheritance pattern (mitochondrial, unknown) remains NOT_APPLICABLE.
    """
    if gene_disease_context.inheritance == Inheritance.AUTOSOMAL_RECESSIVE:
        return interpret_recessive_case(case, classifications, gene_disease_context)
    if gene_disease_context.inheritance in (Inheritance.X_LINKED_RECESSIVE, Inheritance.X_LINKED_DOMINANT):
        return interpret_x_linked_case(case, classifications, gene_disease_context)
    if gene_disease_context.inheritance == Inheritance.AUTOSOMAL_DOMINANT:
        return interpret_dominant_case(case, classifications, gene_disease_context)
    return CaseInterpretation(
        case_id=case.case_id,
        gene=case.gene,
        status=CaseInterpretationStatus.NOT_APPLICABLE,
        rationale=(
            f"inheritance={gene_disease_context.inheritance.value} is not one of the patterns "
            "this project covers (autosomal recessive, X-linked, autosomal dominant)."
        ),
    )
