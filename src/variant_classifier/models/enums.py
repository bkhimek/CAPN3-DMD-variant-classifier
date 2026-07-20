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
