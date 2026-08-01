"""PM5 evaluator -- "novel missense change at an amino acid residue where
a different missense change determined to be pathogenic has been seen
before" (Richards et al. 2015, Table 3). Moderate pathogenic evidence in
the base ACMG/AMP framework. PS1's direct counterpart: same evidence
model (models/same_residue_evidence.py), same splice caveat, same
disclosed precedent-strength downgrade, different comparison (a DIFFERENT
resulting amino acid change at the same residue, not the same one -- see
ps1.py for that side).

Scope notes below intentionally mirror ps1.py's -- read that module's
docstring for the fuller design rationale (why a new evidence field
rather than a new evidence-domain model, why this isn't circular the way
PM3 would be).

1. Only MISSENSE_VARIANT is in scope -- everything else is NOT_APPLICABLE.
2. Depends on bundle.same_residue_evidence.pm5_precedent_established, a
   curated fact about a different variant's own established
   classification, not computed from this project's own fixtures/engine.
3. The real ClinGen LGMD VCEP specification for CAPN3 (v2.0, 2025-07-09)
   requires considerably more than this evaluator checks: REVEL >0.7 for
   the variant under curation, an excluded SpliceAI score for both
   variants, no benign missense variation at the residue, exclusion of
   missense changes encoded by the first/last 3 nucleotides of an exon,
   and counting multiple precedents toward Strength (2P or 3LP at
   Strong). None of that gene-specific machinery is implemented here --
   this evaluator applies the generic Richards et al. 2015 definition and
   the base splice-vs-protein-level caveat only, disclosed rather than
   silently assumed, the same "deliberately partial" treatment PVS1 and
   PS1 both have.
"""

from ..models import CriterionResult, VariantEvidenceBundle
from ..models.enums import Consequence, CriterionStatus, CriterionStrength, EvidenceDirection

RULE_SOURCE = "ACMG/AMP (Richards et al. 2015)"
RULE_VERSION = "2015"


def evaluate_pm5(bundle: VariantEvidenceBundle) -> CriterionResult:
    transcript = next(tc for tc in bundle.transcript_consequences if tc.clinically_relevant)
    evidence_id = f"transcript:{transcript.transcript_id}"

    if transcript.consequence != Consequence.MISSENSE_VARIANT:
        return CriterionResult(
            code="PM5",
            status=CriterionStatus.NOT_APPLICABLE,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"{transcript.consequence.value} is not a missense consequence -- PM5 (\"different "
                "amino acid change at a residue with an established pathogenic change\") only applies "
                "where the amino acid change itself is the expected disease mechanism."
            ),
            evidence_ids=[evidence_id],
        )

    sre = bundle.same_residue_evidence
    if sre is None or sre.pm5_precedent_established is None:
        return CriterionResult(
            code="PM5",
            status=CriterionStatus.NOT_EVALUATED,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale="No same_residue_evidence.pm5_precedent_established recorded -- PM5 was never assessed for this variant.",
            evidence_ids=[evidence_id],
        )

    if sre.pm5_precedent_established is False:
        return CriterionResult(
            code="PM5",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                "Curated precedent search found no previously established pathogenic or likely "
                "pathogenic variant with a different amino acid change at this residue."
            ),
            evidence_ids=[evidence_id],
        )

    # pm5_precedent_established is True.
    if sre.splice_impact_excluded is not True:
        return CriterionResult(
            code="PM5",
            status=CriterionStatus.MANUAL_REVIEW,
            direction=EvidenceDirection.PATHOGENIC,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"A different-amino-acid-change precedent is recorded at this residue "
                f"({sre.pm5_precedent_variant}), but a splice-driven mechanism has not been excluded "
                "for this variant and/or the precedent variant. Per Richards et al. 2015's own caveat "
                "and the ClinGen SVI Splicing Subgroup (Walker et al. 2023, PMID 37352859), PM5's "
                "amino-acid-level comparison is not valid when the nucleotide change may instead be "
                "acting through splicing -- flagged for manual review rather than guessed either way."
            ),
            evidence_ids=[evidence_id],
        )

    if sre.pm5_precedent_classification == "PATHOGENIC":
        strength = CriterionStrength.MODERATE
        downgrade_note = ""
    else:  # LIKELY_PATHOGENIC
        strength = CriterionStrength.SUPPORTING
        downgrade_note = " (downgraded one level: precedent is Likely Pathogenic, not Pathogenic)"

    return CriterionResult(
        code="PM5",
        status=CriterionStatus.MET,
        strength=strength,
        direction=EvidenceDirection.PATHOGENIC,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Different amino acid change ({transcript.hgvs_p}) at the same residue as "
            f"{sre.pm5_precedent_variant}, established {sre.pm5_precedent_classification}, with a "
            f"splice-driven mechanism excluded. Strength {strength.value}{downgrade_note}."
        ),
        evidence_ids=[evidence_id],
    )
