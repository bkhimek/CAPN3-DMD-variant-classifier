"""PS1 evaluator -- "same amino acid change as a previously established
pathogenic variant, regardless of nucleotide change" (Richards et al.
2015, Table 3). Strong pathogenic evidence in the base ACMG/AMP framework.

Added this round alongside PM5 (evaluators/pm5.py) -- see
models/same_residue_evidence.py for the full design writeup (why this
needed a small new evidence field rather than a new evidence-domain
model, why it isn't circular the way PM3 would be, the splice caveat, and
the disclosed precedent-strength downgrade convention).

Scope, disclosed rather than silently assumed:

1. Only MISSENSE_VARIANT is in scope here -- everything else is
   NOT_APPLICABLE. Richards et al. 2015's own PS1 example ("Val->Leu
   caused by either G>C or G>T in the same codon") is a missense case,
   and every real ClinGen VCEP specification found while researching this
   evaluator (CAPN3's LGMD VCEP, RYR1's) restricts PS1 to missense
   variants for which "the amino acid change is the expected mechanism of
   disease." The ClinGen SVI Splicing Subgroup (Walker et al. 2023, PMID
   37352859) also defines a PS1 usage for splice-altering variants, which
   is a genuinely different evaluation (comparing predicted splicing
   impact, not amino acid identity) -- out of scope here, same as this
   project's PVS1 evaluator not resolving splice donor/acceptor variants.
2. Depends on bundle.same_residue_evidence, a curated fact about a
   DIFFERENT variant's own established classification -- never computed
   by scanning this project's other curated fixtures, and never this
   project's own engine output (no circularity).
3. The real ClinGen LGMD VCEP specification for CAPN3 (v2.0, 2025-07-09)
   requires, in addition to what's implemented here: a minimum REVEL
   score for the variant under curation, an excluded SpliceAI score for
   both variants, no benign missense variation at the residue, exclusion
   of missense changes encoded by the first/last 3 nucleotides of an
   exon (a splice-region proxy), and counting multiple pathogenic/likely
   pathogenic precedents toward Strength (2P or 3LP = Strong). None of
   that VCEP-specific machinery is implemented here -- this evaluator
   applies the generic Richards et al. 2015 definition plus the base
   splice-vs-protein-level caveat only, the same "deliberately partial"
   treatment PVS1 has had since Milestone 2.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import Consequence, CriterionStatus, CriterionStrength, EvidenceDirection

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_ps1(bundle: VariantEvidenceBundle) -> CriterionResult:
    transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)
    evidence_id = f"transcript:{transcript.transcript_id}"

    if transcript.consequence != Consequence.MISSENSE_VARIANT:
        return CriterionResult(
            code="PS1",
            status=CriterionStatus.NOT_APPLICABLE,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{transcript.consequence.value} is not a missense consequence -- PS1 (\"same amino "
                "acid change as an established pathogenic variant\") only applies where the amino acid "
                "change itself is the expected disease mechanism."
            ),
            evidence_ids=[evidence_id],
        )

    sre = bundle.same_residue_evidence
    if sre is None or sre.ps1_precedent_established is None:
        return CriterionResult(
            code="PS1",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="No same_residue_evidence.ps1_precedent_established recorded -- PS1 was never assessed for this variant.",
            evidence_ids=[evidence_id],
        )

    if sre.ps1_precedent_established is False:
        return CriterionResult(
            code="PS1",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                "Curated precedent search found no previously established pathogenic or likely "
                "pathogenic variant with this exact same amino acid change."
            ),
            evidence_ids=[evidence_id],
        )

    # ps1_precedent_established is True.
    if sre.splice_impact_excluded is not True:
        return CriterionResult(
            code="PS1",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"A same-amino-acid-change precedent is recorded ({sre.ps1_precedent_variant}), but "
                "a splice-driven mechanism has not been excluded for this variant and/or the precedent "
                "variant. Per Richards et al. 2015's own caveat and the ClinGen SVI Splicing Subgroup "
                "(Walker et al. 2023, PMID 37352859), PS1's amino-acid-level comparison is not valid "
                "when the nucleotide change may instead be acting through splicing -- flagged for "
                "manual review rather than guessed either way."
            ),
            evidence_ids=[evidence_id],
        )

    if sre.ps1_precedent_classification == "PATHOGENIC":
        strength = CriterionStrength.STRONG
    else:  # LIKELY_PATHOGENIC, the only other value SameResidueEvidence permits here
        strength = CriterionStrength.MODERATE

    return CriterionResult(
        code="PS1",
        status=CriterionStatus.MET,
        strength=strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Same amino acid change ({transcript.hgvs_p}) as {sre.ps1_precedent_variant}, established "
            f"{sre.ps1_precedent_classification}, with a splice-driven mechanism excluded. Strength "
            f"{strength.value}" + (" (downgraded one level: precedent is Likely Pathogenic, not Pathogenic)" if strength == CriterionStrength.MODERATE else "") + "."
        ),
        evidence_ids=[evidence_id],
    )
