"""Batch 23 -- DMD CNV/structural-variant scoring: a deliberately partial
implementation of the ACMG/ClinGen Technical Standards for the
Interpretation and Reporting of Constitutional Copy-Number Variants (Riggs
et al. 2020, Genetics in Medicine 22:245-257), scoped to DMD deletions
only. See models/cnv_deletion_evidence.py for the full scope writeup
(what's implemented, what's deferred, and why) and README.md ("DMD
CNV/structural-variant scoring, batch 23") for the batch-level narrative.

Offered as an entirely separate scoring module from engine.py/bayesian.py,
never mixed with either -- a CNV category result (CnvCategoryResult) is
not an ACMG/AMP CriterionResult, and this module's combining math (summed
category points against fixed cutoffs) has its own real source, distinct
from both Table 5 (Richards et al. 2015) and the Tavtigian et al. 2020
Bayesian point system engine.py/bayesian.py already implement for
point-mutation evidence.

Category point values and pathogenicity cutoffs below are quoted directly
from ClassifyCNV (Gurbich TA, Ilinsky VV. "ClassifyCNV: a tool for
clinical annotation of copy-number variants." Sci Rep 10, 20375 (2020),
DOI 10.1038/s41598-020-76425-3) -- an open-source, peer-reviewed
reimplementation of the Riggs et al. 2020 rubric, fetched directly from
github.com/Genotek/ClassifyCNV (resources.py, ClassifyCNV.py) during this
batch's research, since the primary paper and the official ClinGen CNV
calculator (cnvcalc.clinicalgenome.org) were both unreachable. See
cnv_deletion_evidence.py for the disclosure of this secondary-source
reliance.

Decision order, mirroring ClassifyCNV's own assign_del_points_s2() /
analyze_intragenic_deletions() control flow (a whole-gene deletion always
wins over every other check; an intragenic frameshift+NMD deletion is
checked before the 5'/3' end categories; benign-region overlap is checked
before any pathogenic-direction category):

    1. whole_gene_deleted (and gene is an established HI=3 gene)  -> 2A
    2. overlaps_benign_region                                     -> 2F
    3. intragenic + reading_frame_effect == OUT_OF_FRAME + NMD    -> 2E
    4. five_prime_end_deleted and cds_involved                    -> 2C
    5. three_prime_end_deleted and other_exons_involved           -> 2D (0.9)
    6. three_prime_end_deleted and cds_involved (last exon only)  -> 2D (0.3)
    7. otherwise (established HI gene, but none of the above)     -> NONE_APPLICABLE

Because this project's evidence model represents exactly one CNV
overlapping exactly one gene (no multi-gene CNVs, per the Section-3
gene-count category being out of scope), at most one category can ever
fire per CnvDeletionEvidence record -- unlike ClassifyCNV itself, which
sums points across every gene a real multi-gene CNV overlaps.
"""

from typing import Dict

from .models import CnvCategoryResult, CnvDeletionEvidence, CnvDuplicationEvidence, CnvProvisionalClassification
from .models.enums import (
    ClassificationStatus,
    CnvDuplicationOrientation,
    CnvReadingFrameEffect,
    CriterionStatus,
    EvidenceDirection,
    ProvisionalClass,
)

RULE_SOURCE = (
    "ACMG/ClinGen Technical Standards for Constitutional Copy-Number Variants (Riggs et al. 2020, "
    "Genetics in Medicine 22:245-257), Section 2 (loss/deletion) category point values as reimplemented "
    "by ClassifyCNV (Gurbich & Ilinsky 2020, Sci Rep 10:20375)"
)
RULE_VERSION = "2020"

# Pathogenicity cutoffs, quoted directly from ClassifyCNV's resources.py.
PATHOGENIC_CUTOFF = 0.99
LIKELY_PATHOGENIC_CUTOFF = 0.90
LIKELY_BENIGN_CUTOFF = -0.90
BENIGN_CUTOFF = -0.99


def _evidence_id(evidence: CnvDeletionEvidence) -> str:
    return f"cnv:{evidence.cnv_id}"


def _score_category(evidence: CnvDeletionEvidence, dosage_config: Dict[str, dict]) -> CnvCategoryResult:
    eid = _evidence_id(evidence)
    gene_config = dosage_config.get(evidence.gene, {})
    hi_established = bool(gene_config.get("hi_established", False))

    if evidence.whole_gene_deleted:
        if hi_established:
            return CnvCategoryResult(
                code="2A",
                status=CriterionStatus.MET,
                direction=EvidenceDirection.PATHOGENIC,
                points=1.0,
                rule_source=RULE_SOURCE,
                rule_version=RULE_VERSION,
                rationale=(
                    f"Complete deletion of {evidence.gene}, an established ClinGen dosage-sensitive "
                    f"(HI=3) gene ({gene_config.get('source', 'no source recorded')}) -- category 2A."
                ),
                evidence_ids=[eid],
            )
        return CnvCategoryResult(
            code="NONE_APPLICABLE",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.0,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Whole-gene deletion of {evidence.gene}, but this gene is not recorded as an "
                "established (HI=3) dosage-sensitive gene in config/dosage_sensitivity.yaml -- category "
                "2A requires an established gene. The predicted-but-not-established fallback (category "
                "2H) is an explicit, disclosed gap in this milestone (see cnv_deletion_evidence.py), so "
                "no points are assigned rather than guessed at."
            ),
            evidence_ids=[eid],
        )

    if evidence.overlaps_benign_region:
        return CnvCategoryResult(
            code="2F",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.BENIGN,
            points=-1.0,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=f"Deletion falls completely within an established ClinGen benign region for {evidence.gene} -- category 2F.",
            evidence_ids=[eid],
        )

    if (
        evidence.reading_frame_effect == CnvReadingFrameEffect.OUT_OF_FRAME
        and evidence.nmd_predicted is True
    ):
        return CnvCategoryResult(
            code="2E",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.9,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Intragenic deletion of {evidence.gene}"
                + (f" ({evidence.exon_description})" if evidence.exon_description else "")
                + " disrupts the reading frame (Aartsma-Rus reading-frame rule) and is predicted to "
                "trigger nonsense-mediated decay -- category 2E."
            ),
            evidence_ids=[eid],
        )

    if evidence.five_prime_end_deleted and evidence.cds_involved:
        return CnvCategoryResult(
            code="2C",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.9,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Deletion removes the 5' end (5'UTR/first exon) and coding sequence of {evidence.gene} -- "
                "category 2C."
            ),
            evidence_ids=[eid],
        )

    if evidence.three_prime_end_deleted and evidence.other_exons_involved:
        return CnvCategoryResult(
            code="2D",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.9,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Deletion removes the 3' end (3'UTR/last exon) of {evidence.gene} and involves other "
                "exons besides the last one -- category 2D (0.9 pts)."
            ),
            evidence_ids=[eid],
        )

    if evidence.three_prime_end_deleted and evidence.cds_involved:
        return CnvCategoryResult(
            code="2D",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.3,
            rule_source=RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Deletion removes the 3' end of {evidence.gene}, confined to the last exon's coding "
                "sequence with no other exons involved -- category 2D (0.3 pts)."
            ),
            evidence_ids=[eid],
        )

    return CnvCategoryResult(
        code="NONE_APPLICABLE",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.PATHOGENIC,
        points=0.0,
        rule_source=RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Intragenic deletion of {evidence.gene}"
            + (f" ({evidence.exon_description})" if evidence.exon_description else "")
            + f" with reading_frame_effect={evidence.reading_frame_effect.value if evidence.reading_frame_effect else 'None'} "
            "-- does not match any Section 2 category this project implements (2A/2C/2D/2E/2F). Real "
            "in-frame internal deletions of this shape (e.g. the classic DMD exon 45-47 Becker deletion) "
            "likely correspond to a real Riggs 2020 category this project has not independently verified "
            "(candidates: 2B, 2G) -- see cnv_deletion_evidence.py. Reported as zero Section-2 points "
            "rather than guessed at."
        ),
        evidence_ids=[eid],
    )


def _classify_points(points: float):
    if points >= PATHOGENIC_CUTOFF:
        return ProvisionalClass.PATHOGENIC, f"{points:.2f} points (>= {PATHOGENIC_CUTOFF})"
    if points >= LIKELY_PATHOGENIC_CUTOFF:
        return ProvisionalClass.LIKELY_PATHOGENIC, f"{points:.2f} points ({LIKELY_PATHOGENIC_CUTOFF} to {PATHOGENIC_CUTOFF})"
    if points > LIKELY_BENIGN_CUTOFF:
        return ProvisionalClass.VUS, f"{points:.2f} points ({LIKELY_BENIGN_CUTOFF} to {LIKELY_PATHOGENIC_CUTOFF})"
    if points > BENIGN_CUTOFF:
        return ProvisionalClass.LIKELY_BENIGN, f"{points:.2f} points ({BENIGN_CUTOFF} to {LIKELY_BENIGN_CUTOFF})"
    return ProvisionalClass.BENIGN, f"{points:.2f} points (<= {BENIGN_CUTOFF})"


def score_cnv_deletion(evidence: CnvDeletionEvidence, dosage_config: Dict[str, dict]) -> CnvProvisionalClassification:
    """Score a single DMD deletion against this project's implemented
    slice of Riggs et al. 2020 Section 2, and classify the result via
    ClassifyCNV's published cutoffs. dosage_config is the dict returned by
    loader.load_dosage_sensitivity().
    """
    category = _score_category(evidence, dosage_config)
    provisional_class, band_note = _classify_points(category.points)

    rationale = (
        f"{provisional_class.value} via CNV Section-2 scoring: {category.points:.2f} points from category "
        f"{category.code} ({band_note}). {category.rationale}"
    )

    return CnvProvisionalClassification(
        provisional_class=provisional_class,
        status=ClassificationStatus.PROVISIONAL_AUTOMATED,
        categories=[category],
        combining_rule_source=RULE_SOURCE,
        combining_rule_version=RULE_VERSION,
        rationale=rationale,
        points=category.points,
        manual_review_required=False,
    )


# --------------------------------------------------------------------------
# Batch 24: DMD duplication (gain) scoring.
#
# A deliberately narrower slice than the deletion side, per this batch's
# research and the scope confirmed with the user before writing any of
# this: no whole-gene triplosensitivity scoring (ClinGen's own DMD dosage
# curation says whole-gene DMD duplications aren't clinically reported;
# this project has no verified TS category or config to apply anyway).
# What IS implemented is the real, clinically-relevant DMD mechanism --
# an intragenic tandem duplication disrupting the gene via the same
# Aartsma-Rus reading-frame rule already used for deletions -- with the
# exact point-value-to-condition mapping explicitly disclosed as this
# project's own inference, not an independently verified primary-source
# fact. See models/cnv_duplication_evidence.py's module docstring and
# models/enums.py's CNV_GAIN_CATEGORY_CODES comment for the full
# disclosure this rationale text below summarizes.

GAIN_RULE_SOURCE = (
    "ACMG/ClinGen Technical Standards for Constitutional Copy-Number Variants (Riggs et al. 2020, "
    "Genetics in Medicine 22:245-257), gain-side Section 2 category existence and point-value "
    "fragments (2K=0.45, 2J=0, for a gain CNV breakpoint observed within an established HI gene) "
    "per PMC8960312 (inter-laboratory CNV concordance study); benign-region-overlap point value "
    "(-1.0) per ClassifyCNV (Gurbich & Ilinsky 2020, Sci Rep 10:20375). The condition-to-value "
    "mapping (out-of-frame -> 0.45, in-frame/unknown -> 0) is this project's own disclosed "
    "inference by analogy to the loss side, NOT independently verified against the primary paper."
)


def _duplication_evidence_id(evidence: CnvDuplicationEvidence) -> str:
    return f"cnv:{evidence.cnv_id}"


def _score_duplication_category(evidence: CnvDuplicationEvidence) -> CnvCategoryResult:
    eid = _duplication_evidence_id(evidence)

    if evidence.overlaps_benign_region:
        return CnvCategoryResult(
            code="GAIN_BENIGN",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.BENIGN,
            points=-1.0,
            rule_source=GAIN_RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Duplication falls completely within an established ClinGen benign region for "
                f"{evidence.gene} -- point value confirmed directly from ClassifyCNV's "
                "assign_dup_points_s2()."
            ),
            evidence_ids=[eid],
        )

    if evidence.whole_gene_duplicated:
        return CnvCategoryResult(
            code="NONE_APPLICABLE",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.0,
            rule_source=GAIN_RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Whole-gene duplication of {evidence.gene}. This project does not implement "
                "triplosensitivity (TS) scoring -- ClinGen's own DMD dosage curation states whole-"
                "gene DMD duplications are not clinically reported, and no real DMD fixture would "
                "exercise a TS category honestly. Explicit, disclosed gap, not a guess."
            ),
            evidence_ids=[eid],
        )

    if not evidence.breakpoint_within_gene:
        return CnvCategoryResult(
            code="NONE_APPLICABLE",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.0,
            rule_source=GAIN_RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Duplication of {evidence.gene} has no breakpoint inside the gene and does not "
                "overlap a benign region or the whole gene -- no implemented gain category applies."
            ),
            evidence_ids=[eid],
        )

    if evidence.is_tandem != CnvDuplicationOrientation.TANDEM:
        return CnvCategoryResult(
            code="NONE_APPLICABLE",
            status=CriterionStatus.NOT_MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.0,
            rule_source=GAIN_RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Breakpoint inside {evidence.gene}, but tandem/direct orientation is "
                f"{evidence.is_tandem.value if evidence.is_tandem else 'not stated'} -- no functional "
                "call can be made without confirmed tandem orientation (a breakpoint study of 119 "
                "gain CNVs found 83% tandem/direct, with the majority of the remainder interpreted "
                "as VUS because the effect could not be determined)."
            ),
            evidence_ids=[eid],
        )

    if evidence.reading_frame_effect == CnvReadingFrameEffect.OUT_OF_FRAME:
        return CnvCategoryResult(
            code="GAIN_2K_EQUIV",
            status=CriterionStatus.MET,
            direction=EvidenceDirection.PATHOGENIC,
            points=0.45,
            rule_source=GAIN_RULE_SOURCE,
            rule_version=RULE_VERSION,
            rationale=(
                f"Confirmed-tandem duplication with a breakpoint inside {evidence.gene}"
                + (f" ({evidence.exon_description})" if evidence.exon_description else "")
                + ", predicted out-of-frame (Aartsma-Rus reading-frame rule) -- this project's "
                "disclosed GAIN_2K_EQUIV label (0.45 points, the higher of two secondary-sourced "
                "gain-side breakpoint-in-HI-gene point values; see module docstring for why this "
                "specific condition-to-value mapping is an inference, not a verified fact)."
            ),
            evidence_ids=[eid],
        )

    return CnvCategoryResult(
        code="GAIN_2J_EQUIV",
        status=CriterionStatus.NOT_MET,
        direction=EvidenceDirection.PATHOGENIC,
        points=0.0,
        rule_source=GAIN_RULE_SOURCE,
        rule_version=RULE_VERSION,
        rationale=(
            f"Confirmed-tandem duplication with a breakpoint inside {evidence.gene}"
            + (f" ({evidence.exon_description})" if evidence.exon_description else "")
            + f", reading_frame_effect={evidence.reading_frame_effect.value if evidence.reading_frame_effect else 'None'} "
            "(in-frame or unknown) -- this project's disclosed GAIN_2J_EQUIV label (0 points, the "
            "lower of two secondary-sourced gain-side breakpoint-in-HI-gene point values)."
        ),
        evidence_ids=[eid],
    )


def score_cnv_duplication(evidence: CnvDuplicationEvidence) -> CnvProvisionalClassification:
    """Score a single DMD duplication against this project's implemented
    slice of gain-side CNV evidence, and classify via the same
    ClassifyCNV-published cutoffs used for deletions. Unlike
    score_cnv_deletion(), no dosage_config is needed -- the implemented
    categories here don't depend on an established-gene lookup (whole-gene
    TS scoring, which would need one, is explicitly not implemented).
    """
    category = _score_duplication_category(evidence)
    provisional_class, band_note = _classify_points(category.points)

    rationale = (
        f"{provisional_class.value} via CNV Section-2 gain scoring: {category.points:.2f} points from "
        f"category {category.code} ({band_note}). {category.rationale}"
    )

    return CnvProvisionalClassification(
        provisional_class=provisional_class,
        status=ClassificationStatus.PROVISIONAL_AUTOMATED,
        categories=[category],
        combining_rule_source=GAIN_RULE_SOURCE,
        combining_rule_version=RULE_VERSION,
        rationale=rationale,
        points=category.points,
        manual_review_required=False,
    )
