"""Controlled vocabularies shared by every model in this package.

These enums are the internal representation the Evidence Normalisation
module (see the Workflow Architecture Guide) is responsible for producing.
Curated fixtures for this milestone are hand-normalised into these values
directly; nothing here should ever be a source-specific raw string
(ClinVar's "Likely pathogenic" vs "likely_pathogenic" vs "LP" all become one
value here: LIKELY_PATHOGENIC).
"""

from enum import Enum


class GenomeBuild(str, Enum):
    GRCH38 = "GRCh38"


class Inheritance(str, Enum):
    AUTOSOMAL_DOMINANT = "AUTOSOMAL_DOMINANT"
    AUTOSOMAL_RECESSIVE = "AUTOSOMAL_RECESSIVE"
    X_LINKED_RECESSIVE = "X_LINKED_RECESSIVE"
    X_LINKED_DOMINANT = "X_LINKED_DOMINANT"
    MITOCHONDRIAL = "MITOCHONDRIAL"
    UNKNOWN = "UNKNOWN"


class DiseaseMechanism(str, Enum):
    LOSS_OF_FUNCTION = "LOSS_OF_FUNCTION"
    GAIN_OF_FUNCTION = "GAIN_OF_FUNCTION"
    DOMINANT_NEGATIVE = "DOMINANT_NEGATIVE"
    UNKNOWN = "UNKNOWN"


class SpecificationType(str, Enum):
    GENERIC_ACMG = "GENERIC_ACMG"
    VCEP = "VCEP"


class Consequence(str, Enum):
    """A deliberately small subset of Sequence Ontology consequence terms —
    just enough to cover Milestone 1's CAPN3/DMD cases. Extend as needed;
    do not silently accept an arbitrary string here (see loader.py)."""

    FRAMESHIFT_VARIANT = "frameshift_variant"
    STOP_GAINED = "stop_gained"
    MISSENSE_VARIANT = "missense_variant"
    SYNONYMOUS_VARIANT = "synonymous_variant"
    SPLICE_DONOR_VARIANT = "splice_donor_variant"
    SPLICE_ACCEPTOR_VARIANT = "splice_acceptor_variant"
    INFRAME_DELETION = "inframe_deletion"
    INFRAME_INSERTION = "inframe_insertion"
    START_LOST = "start_lost"
    STOP_LOST = "stop_lost"
    INTRON_VARIANT = "intron_variant"
    FIVE_PRIME_UTR_VARIANT = "5_prime_UTR_variant"
    THREE_PRIME_UTR_VARIANT = "3_prime_UTR_variant"
    OTHER = "other"


class PopulationRetrievalStatus(str, Enum):
    """Mirrors the missing-data states defined in the Reporting and
    Dashboard Design Guide (Section 4.2), plus OBSERVED for the ordinary
    successful-retrieval case: a variant that was never found is not the
    same as a variant whose locus could not be assessed, and neither is the
    same as a variant that WAS found and has real frequency data attached.

    OBSERVED was added during Milestone-1 fixture-writing after finding
    that the original five states were all "missing/negative" outcomes,
    with no status for a normal successful lookup that returns a
    nonzero-AF variant (e.g. the real CAPN3 c.550del founder-mutation
    case) — a genuine gap, not a stylistic choice."""

    OBSERVED = "OBSERVED"
    ABSENT = "ABSENT"
    NOT_ASSESSED = "NOT_ASSESSED"
    UNAVAILABLE = "UNAVAILABLE"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    UNKNOWN = "UNKNOWN"


class CriterionStatus(str, Enum):
    """The six-state criterion model settled across the whole guide set
    (Workflow Architecture Guide Rev. 6+, ACMG Engine Detailed Design Guide,
    Validation and Verification Guide)."""

    MET = "MET"
    NOT_MET = "NOT_MET"
    NOT_EVALUATED = "NOT_EVALUATED"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    MANUAL_REVIEW = "MANUAL_REVIEW"
    CONFLICTING_EVIDENCE = "CONFLICTING_EVIDENCE"


class CriterionStrength(str, Enum):
    STAND_ALONE = "STAND_ALONE"
    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"


class EvidenceDirection(str, Enum):
    PATHOGENIC = "PATHOGENIC"
    BENIGN = "BENIGN"


class AutomationConfidence(str, Enum):
    AUTOMATED_HIGH = "AUTOMATED_HIGH"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"


class ComputationalPrediction(str, Enum):
    """The calibrated output of a single computational-evidence source for
    PP3/BP4 — a prediction bucket, not a raw tool score. Per the ACMG
    Engine Detailed Design Guide, "Computational and splicing evidence":
    one calibrated call per variant, not several correlated tool votes."""

    PATHOGENIC = "PATHOGENIC"
    BENIGN = "BENIGN"
    INDETERMINATE = "INDETERMINATE"


class ProvisionalClass(str, Enum):
    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    VUS = "VUS"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    BENIGN = "BENIGN"


class ClassificationStatus(str, Enum):
    """A classification is PROVISIONAL_AUTOMATED until Scientist Review &
    Sign-off marks it FINAL. Milestone 1 has no review step, so every
    ProvisionalClassification produced here stays PROVISIONAL_AUTOMATED —
    Reporting (out of scope for this project) would refuse to render
    anything else, and this prototype makes the same rule true by construction."""

    PROVISIONAL_AUTOMATED = "PROVISIONAL_AUTOMATED"
    FINAL = "FINAL"


class KaryotypicSex(str, Enum):
    """Chromosomal sex as relevant to X-linked inheritance reasoning
    (hemizygous vs heterozygous X), not a stand-in for gender. OTHER
    covers real karyotypic variation (e.g. XXY, X0, mosaicism) that the
    Milestone 4 X-linked interpretation logic does not attempt to reason
    about — it is deferred to manual review, the same as XX."""

    XY = "XY"
    XX = "XX"
    OTHER = "OTHER"
    UNKNOWN = "UNKNOWN"


class PhaseRelationship(str, Enum):
    """Whether two variants in the same gene, found in the same patient,
    are confirmed on different copies of the chromosome (TRANS — each
    parent contributed one broken copy) or the same copy (CIS — one copy
    is doubly hit, the other is untouched). Established via parental
    testing or phasing sequencing; UNKNOWN when neither is available."""

    TRANS = "TRANS"
    CIS = "CIS"
    UNKNOWN = "UNKNOWN"


class CaseInterpretationStatus(str, Enum):
    """The outcome of Milestone 4's case-level (not variant-level)
    reasoning: does what was found in this patient explain their disease,
    given how the gene's disease is inherited?"""

    EXPLAINED = "EXPLAINED"
    INSUFFICIENT = "INSUFFICIENT"
    MANUAL_REVIEW = "MANUAL_REVIEW"
    NOT_APPLICABLE = "NOT_APPLICABLE"


# The full 28-code ACMG/AMP controlled vocabulary (Richards et al. 2015).
# Milestone 1 only *evaluates* a subset (see SUPPORTED_CRITERIA_MILESTONE_1)
# but fixtures and models accept any code from this set, so adding an
# evaluator later never requires a schema change.
ACMG_CRITERION_CODES = frozenset(
    {
        "PVS1",
        "PS1", "PS2", "PS3", "PS4",
        "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
        "PP1", "PP2", "PP3", "PP4", "PP5",
        "BA1",
        "BS1", "BS2", "BS3", "BS4",
        "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
    }
)

# Criteria this Milestone 1 prototype actually evaluates against curated
# fixtures. PP3/BP4 are included as a pair deliberately: they are the same
# calibrated computational-evidence family evaluated in opposite directions
# (see the ACMG Engine Detailed Design Guide, "Computational and splicing
# evidence"), and without BP4 the combination engine cannot reach
# LIKELY_BENIGN (BS1 alone is only one Strong-Benign criterion; Likely
# Benign needs 1 Strong + 1 Supporting, or 2 Supporting).
SUPPORTED_CRITERIA_MILESTONE_1 = frozenset({"PVS1", "PM2", "PP3", "BP4", "BA1", "BS1"})


class CnvReadingFrameEffect(str, Enum):
    """Whether an intragenic (both gene ends intact) deletion is predicted
    to shift the reading frame, per the Aartsma-Rus DMD reading-frame rule
    (Aartsma-Rus et al. 2006, PMID 16770791; 2019 update, Human Mutation).
    Curated explicitly per CnvDeletionEvidence record -- never computed by
    this project from raw deleted-base counts, the same "state the fact,
    don't derive it" convention TranscriptConsequence.nmd_predicted and
    SameResidueEvidence.splice_impact_excluded already use. UNKNOWN is a
    real, legitimate state (e.g. breakpoints not precisely resolved to the
    base pair) distinct from simply omitting the field."""

    OUT_OF_FRAME = "OUT_OF_FRAME"
    IN_FRAME = "IN_FRAME"
    UNKNOWN = "UNKNOWN"


# The Section 2 (loss/deletion, dosage-sensitivity) category codes this
# project's cnv_scoring.py actually evaluates, per the ACMG/ClinGen
# Technical Standards for Copy-Number Variants (Riggs et al. 2020), point
# values as reimplemented by ClassifyCNV (Gurbich & Ilinsky 2020, Sci Rep
# 10:20375). NONE_APPLICABLE is this project's OWN bookkeeping label (not
# a Riggs/ClassifyCNV code) for a deletion that is intragenic, established
# to be in an HI gene, but whose reading-frame effect is IN_FRAME or
# UNKNOWN -- i.e. none of 2A/2C/2D/2E/2F match. The real Riggs rubric
# likely has a specific code for this shape (candidates seen in secondary
# sources include 2B and 2G) but this project has not independently
# verified either one's exact definition or point value, so rather than
# guess, it reports zero Section-2 points under this disclosed
# project-internal label. See cnv_deletion_evidence.py and cnv_scoring.py
# for the full scope writeup.
CNV_LOSS_CATEGORY_CODES = frozenset({"2A", "2C", "2D", "2E", "2F", "NONE_APPLICABLE"})


class CnvDuplicationOrientation(str, Enum):
    """Whether a duplication with a breakpoint inside a gene has been
    confirmed to be a tandem, direct-orientation insertion adjacent to the
    original copy -- the configuration needed before ANY functional
    interpretation (reading-frame effect included) can be applied at all.
    Riggs et al. 2020 / the CNV-interpretation literature this project
    found while researching duplications (a breakpoint study of 119 gain
    CNVs: 83% tandem and direct, with "the majority of the remainder ...
    interpreted as VUS because the effect could not be determined") treats
    an unconfirmed or non-tandem/complex insertion as functionally
    unpredictable -- NOT a synonym for "no evidence," but a real, distinct
    reason no functional call can be made. Curated explicitly per
    CnvDuplicationEvidence record, never assumed tandem by default."""

    TANDEM = "TANDEM"
    NOT_TANDEM_OR_COMPLEX = "NOT_TANDEM_OR_COMPLEX"
    UNKNOWN = "UNKNOWN"


# The Section 2 (gain/duplication) category labels this project's
# cnv_scoring.py evaluates. Unlike CNV_LOSS_CATEGORY_CODES, NONE of these
# are asserted to be the real Riggs et al. 2020 letter codes -- research
# for batch 24 found secondary-source evidence (an inter-laboratory
# concordance study, PMC8960312, discussing "the use of 2K (0.45 points)
# or 2J (0 point) when a copy number gain breakpoint was observed for the
# established HI genes") that the real rubric numbers gain-side
# haploinsufficiency-breakpoint categories distinctly from the loss side
# (which occupies 2A-2H) -- almost certainly continuing the letter
# sequence (2I onward), NOT reusing 2A/2C/2D/2F the way ClassifyCNV's own
# code does internally (a dict-key-reuse simplification in that tool,
# confirmed by reading its source directly, not proof of the true Riggs
# lettering for gains). Rather than assert an unverified real code, this
# project uses its own disclosed "_EQUIV" labels naming which secondary-
# sourced point value each maps to:
#   GAIN_2K_EQUIV (0.45 pts) -- a confirmed-tandem duplication with a
#     breakpoint inside an established HI gene, predicted out-of-frame.
#   GAIN_2J_EQUIV (0 pts) -- same, but in-frame or frame effect unknown.
#   GAIN_BENIGN (-1.0 pts) -- falls completely within an established
#     benign copy-number region (point value IS confirmed directly from
#     ClassifyCNV's assign_dup_points_s2()).
#   NONE_APPLICABLE (0 pts) -- whole-gene duplication (this project does
#     not implement triplosensitivity scoring, and ClinGen's own DMD
#     dosage curation states whole-gene DMD duplications are not
#     clinically reported anyway), or a duplication whose tandem/direct
#     orientation is not confirmed (per the 83%-tandem literature above).
# See models/cnv_duplication_evidence.py and cnv_scoring.py for the full
# writeup and the exact condition-to-value mapping's disclosed uncertainty.
CNV_GAIN_CATEGORY_CODES = frozenset({"GAIN_2K_EQUIV", "GAIN_2J_EQUIV", "GAIN_BENIGN", "NONE_APPLICABLE"})


class FunctionalAssayResult(str, Enum):
    """The three real outcomes a curated functional assay (PS3/BS3, see
    models/functional_evidence.py) can report. ABNORMAL and NORMAL are
    self-explanatory; INDETERMINATE is a distinct, real third state --
    an assay that was performed but did not clearly discriminate
    pathogenic from benign for this specific variant (e.g. a Western
    blot showing "variably reduced" or ambiguous protein expression).
    INDETERMINATE is NOT the same as having no functional evidence at
    all (functional_evidence absent from the bundle) -- it is curated
    explicitly, same "never silently guess" convention used throughout
    this project, and results in NOT_MET for both PS3 and BS3 rather
    than NOT_EVALUATED."""

    ABNORMAL = "ABNORMAL"
    NORMAL = "NORMAL"
    INDETERMINATE = "INDETERMINATE"
