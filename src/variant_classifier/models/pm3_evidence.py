"""Pm3Evidence / Pm3ProbandObservation -- curated evidence for PM3
("for recessive disorders, detected in trans with a pathogenic variant"),
Richards et al. 2015, Table 3. Added batch 28, closing a gap this project
has disclosed since Milestone 4 (see clinical.py's module docstring and
engine.py's "what this engine does NOT do" section): PM3 could not be
bolted onto the per-variant evaluator pattern the way PS1/PM5/PS3/BS3
were, because a per-variant evaluator taking only its own
VariantEvidenceBundle cannot see a *second* variant's classification --
and re-deriving that second classification live, from this same engine,
inside the first variant's own evaluation would be circular (variant A's
PM3 depends on variant B's ProvisionalClassification, which in a fully
general case could itself depend on variant B's own PM3, depending on
variant A).

This module resolves that circularity the same way the project has
resolved every other "the exact computation isn't independently
verifiable" gap since Milestone 1: by curating the fact a human already
knows, rather than having the engine re-derive it. Real-world PM3
curation already works this way in practice -- a curator citing a
published compound-heterozygous case report is trusting that report's
own stated classification of the partner allele (often itself a ClinVar
entry, a VCEP call, or the original paper's own ACMG interpretation), not
re-running this project's own 12-evaluator engine on a variant that may
not even be in this project's own curated set. other_allele_classification
below is exactly that already-known, curator-supplied fact.

Real, quoted rules this module implements directly, sourced from the
ACGS 2024 UK Practice Guidelines for Variant Classification (v1, Aug
2024) -- the same document already used for PVS1's start-loss rules
(batch 27):

  - Points-based system: "each proband is awarded a point value based
    upon phasing of the two variants in question (confirmed in trans
    versus unknown) and classification of the variant on the other
    allele. The combined point value of all proband occurrences is then
    summed."
  - Homozygous cap: homozygous observations of the variant being
    evaluated are capped at a maximum of 1 point each, regardless of how
    the homozygosity was confirmed -- enforced in
    Pm3ProbandObservation.__post_init__ below.
  - Phase can be established via direct parental testing OR proxy
    methods (informative SNP linkage on NGS reads); "unknown" (neither
    available) is a real, distinct, lower-weighted third state, not an
    error -- this is why phase is Optional[PhaseRelationship] restricted
    to TRANS/UNKNOWN rather than a bool.
  - Cis-cooccurrence override: "PM3 should not be applied at any level
    in the context of two variants that predominantly co-occur" (i.e.
    are found together on the same chromosome copy across the population,
    typically established via the gnomAD variant co-occurrence tool,
    https://gnomad.broadinstitute.org/variant-cooccurrence) -- modeled
    here as Pm3Evidence.cis_cooccurrence_observed, a hard override the
    evaluator checks before summing any proband points at all.

What this module deliberately does NOT do: assert an exact points-per-
scenario table (e.g. "confirmed trans + Pathogenic partner = 1.0 point,
unknown phase + Likely Pathogenic partner = 0.25 points"). That table is
published in the ClinGen SVI "Recommendation for the in trans Criterion
(PM3)" Version 1.0, whose primary PDF
(clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf)
returned empty/unreadable on every fetch attempt this session (as did its
clinicalgenome.org doc page) -- consistent with this project's repeated
experience of ClinGen SVI primary documents being unreachable via
web_fetch (Walker et al. 2023's PVS1-splicing paper and CAPN3 PVS1
flowchart, batch 26, had the identical problem). Rather than hardcode a
recalled-but-unverified numeric table, `points` on each
Pm3ProbandObservation is itself a curated fact -- the curator states the
real point value directly (from their own knowledge of the true SVI
table), exactly the same "never silently guess, state the decision-
relevant fact explicitly" pattern already used for
FunctionalEvidence.validation_strength (batch 25). The only numeric
constraint this project enforces itself is the one rule that WAS directly
quoted and confirmed (the homozygous 1-point cap); everything else about
how many points a given proband is worth is left to the curator, same as
Brnich et al. 2019's OddsPath tier is left to the curator for PS3/BS3
rather than this project computing it from raw assay data.

Thresholds for combining summed points into a PM3 strength ARE a
confirmed, real, gene-specific fact -- the ClinGen LGMD VCEP CAPN3
specification v2.0's own PM3 points table (cspec.genome.network/cspec/ui/svi/doc/GN187):
Very Strong >=4 points, Strong >=2 but <4, Moderate >=1 but <2,
Supporting >=0.5 but <1. See evaluators/pm3.py for how these are applied.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_dict, require_list
from .enums import PhaseRelationship, Pm3Zygosity, ProvisionalClass

_QUALIFYING_OTHER_ALLELE = (ProvisionalClass.PATHOGENIC, ProvisionalClass.LIKELY_PATHOGENIC)
_HOMOZYGOUS_POINT_CAP = 1.0


@dataclass(frozen=True)
class Pm3ProbandObservation:
    proband_id: str
    zygosity: Pm3Zygosity
    other_allele_classification: ProvisionalClass
    points: float
    phase: Optional[PhaseRelationship] = None
    source: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        context = f"Pm3ProbandObservation[{self.proband_id}]"
        if not self.proband_id.strip():
            raise SchemaValidationError("Pm3ProbandObservation: proband_id must not be empty")
        if not isinstance(self.zygosity, Pm3Zygosity):
            raise SchemaValidationError(f"{context}: zygosity must be a Pm3Zygosity value")
        if self.other_allele_classification not in _QUALIFYING_OTHER_ALLELE:
            raise SchemaValidationError(
                f"{context}: other_allele_classification="
                f"{self.other_allele_classification.value} is not valid for PM3 -- the partner "
                "allele must itself be PATHOGENIC or LIKELY_PATHOGENIC (Richards et al. 2015: "
                "'detected in trans with a pathogenic variant'); a partner allele that isn't at "
                "least Likely Pathogenic does not support PM3 for this proband, so it should not "
                "be recorded as a contributing observation at all"
            )
        if self.points <= 0:
            raise SchemaValidationError(
                f"{context}: points={self.points} must be > 0 -- a non-contributing observation "
                "should simply be omitted from Pm3Evidence.probands rather than recorded at zero"
            )
        if self.zygosity == Pm3Zygosity.HOMOZYGOUS:
            if self.phase is not None:
                raise SchemaValidationError(
                    f"{context}: phase must not be set when zygosity=HOMOZYGOUS -- phase (trans vs. "
                    "cis vs. unknown) describes the relationship between two *different* variants; "
                    "a homozygous proband carries the identical variant on both copies, so phase is "
                    "not a meaningful fact here"
                )
            if self.points > _HOMOZYGOUS_POINT_CAP:
                raise SchemaValidationError(
                    f"{context}: points={self.points} exceeds the real ACGS 2024 homozygous cap of "
                    f"{_HOMOZYGOUS_POINT_CAP} -- 'homozygous variants should be capped to a maximum "
                    "of 1 point ... regardless of whether the variant is confirmed via parental "
                    "testing or not'"
                )
        else:  # COMPOUND_HETEROZYGOUS
            if self.phase is None:
                raise SchemaValidationError(
                    f"{context}: phase must be explicitly TRANS or UNKNOWN when "
                    "zygosity=COMPOUND_HETEROZYGOUS -- never silently omitted, same convention as "
                    "ClinicalCase.phase"
                )
            if self.phase == PhaseRelationship.CIS:
                raise SchemaValidationError(
                    f"{context}: phase=CIS is not valid for a Pm3ProbandObservation -- a proband "
                    "confirmed in cis contributes zero PM3 points and should simply be omitted from "
                    "Pm3Evidence.probands, not recorded with phase=CIS. (Population-level cis "
                    "co-occurrence, e.g. via gnomAD, is a separate fact -- "
                    "Pm3Evidence.cis_cooccurrence_observed -- that overrides PM3 for the variant as "
                    "a whole, not per-proband.)"
                )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "Pm3ProbandObservation":
        data = require_dict(data, context or "Pm3ProbandObservation")
        proband_id = data.get("proband_id")
        if not isinstance(proband_id, str) or not proband_id.strip():
            raise SchemaValidationError(f"{context or 'Pm3ProbandObservation'}: 'proband_id' must be a non-empty string")
        ctx = f"Pm3ProbandObservation[{proband_id}]"
        zygosity = coerce_enum(Pm3Zygosity, data.get("zygosity"), "zygosity", ctx)
        other_allele_classification = coerce_enum(
            ProvisionalClass, data.get("other_allele_classification"), "other_allele_classification", ctx
        )
        points_raw = data.get("points")
        if not isinstance(points_raw, (int, float)) or isinstance(points_raw, bool):
            raise SchemaValidationError(f"{ctx}: 'points' must be a number, got {points_raw!r}")
        points = float(points_raw)
        phase_raw = data.get("phase")
        phase = coerce_enum(PhaseRelationship, phase_raw, "phase", ctx) if phase_raw is not None else None
        return cls(
            proband_id=proband_id,
            zygosity=zygosity,
            other_allele_classification=other_allele_classification,
            points=points,
            phase=phase,
            source=optional_str(data, "source"),
            notes=optional_str(data, "notes"),
        )


@dataclass(frozen=True)
class Pm3Evidence:
    probands: List[Pm3ProbandObservation] = field(default_factory=list)
    cis_cooccurrence_observed: bool = False
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.probands and not self.cis_cooccurrence_observed:
            raise SchemaValidationError(
                "Pm3Evidence: probands is empty and cis_cooccurrence_observed is False -- an empty, "
                "no-op Pm3Evidence record should simply not be attached to the bundle "
                "(pm3_evidence=None) rather than curated as an empty object"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "Pm3Evidence":
        data = require_dict(data, context or "Pm3Evidence")
        ctx = context or "Pm3Evidence"
        raw_probands = require_list(data, "probands", ctx)
        probands = [
            Pm3ProbandObservation.from_dict(p, f"{ctx}.probands[{i}]") for i, p in enumerate(raw_probands)
        ]
        cis_raw = data.get("cis_cooccurrence_observed", False)
        if not isinstance(cis_raw, bool):
            raise SchemaValidationError(f"{ctx}: 'cis_cooccurrence_observed' must be true/false if provided")
        return cls(
            probands=probands,
            cis_cooccurrence_observed=cis_raw,
            notes=optional_str(data, "notes"),
        )
