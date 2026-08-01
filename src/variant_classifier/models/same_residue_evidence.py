"""SameResidueEvidence -- curated precedent evidence for PS1 ("same amino
acid change as a previously established pathogenic variant, regardless of
nucleotide change") and PM5 ("novel missense change at an amino acid
residue where a different missense change determined to be pathogenic has
been seen before"), Richards et al. 2015, Table 3.

Added alongside evaluators/ps1.py and evaluators/pm5.py -- the first new
evidence added since PM4 (batch 14), and, per the same "no new evidence
type needed" scoping that made PS1/PM5 this round's first priority, this
is a small addition to the existing evidence-bundle shape rather than a
new bundle-level evidence domain the way ComputationalEvidence was for
PP3/BP4.

Both PS1 and PM5 depend on a fact about a DIFFERENT, previously classified
variant (sourced from ClinVar or a VCEP curation) -- never on this
project's own engine output. That is a real, deliberate difference from
PM3 ("detected in trans with a pathogenic variant"), which clinical.py's
own docstring explains cannot be a per-variant evaluator because variant
A's PM3 would depend on variant B's classification from this same engine,
a structural circularity. PS1/PM5 reference an external authority's
already-established classification instead, so no such circularity
exists here -- the precedent is curated data, exactly like PM2's
population-frequency thresholds or PP3/BP4's calibrated computational
score, not something this evaluator computes by cross-referencing other
curated bundles in this project's own fixture set.

Splice caveat: both criteria are explicitly not safe to apply when the
nucleotide change under evaluation (or the precedent variant's nucleotide
change) might itself be acting through altered splicing rather than
through the amino acid substitution -- Richards et al. 2015's own caveat
for both criteria, refined by the ClinGen SVI Splicing Subgroup (Walker
et al. 2023, PMID 37352859), which specifically lists PS1 and PM5 among
the codes its splicing-evidence framework covers. `splice_impact_excluded`
must be stated explicitly (True or False) whenever a precedent is
recorded, never left unstated -- same "never silently guess" convention
as `nmd_predicted` (PVS1) and `repeat_region` (PM4).

Precedent strength: the base Richards et al. 2015 framework does not
formally vary PS1/PM5 strength by how strongly the PRECEDENT variant
itself was classified (PS1 is always Strong, PM5 is always Moderate).
This project downgrades one level (PS1 -> Moderate, PM5 -> Supporting)
when the precedent is only Likely Pathogenic rather than Pathogenic -- a
disclosed simplification, not an invented rule: at least one real ClinGen
VCEP specification (RYR1, Malignant Hyperthermia Susceptibility) already
documents exactly this downgrade convention, though this project has not
verified it is a universal ClinGen SVI mandate. The real ClinGen LGMD
VCEP specification for CAPN3 (v2.0, 2025-07-09) goes considerably further
than either of these -- requiring a minimum REVEL score, an excluded
SpliceAI score, no benign missense variation at the residue, exclusion of
first/last-3-nucleotide-of-exon codons, and counting MULTIPLE precedent
variants toward a Strong-level PS1/PM5 -- none of which is implemented
here. That gap is disclosed in ps1.py/pm5.py and README.md, the same
"deliberately partial, gap named rather than hidden" treatment PVS1 has
had since Milestone 2.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import optional_bool, optional_str, require_dict

_VALID_PRECEDENT_CLASSIFICATIONS = ("PATHOGENIC", "LIKELY_PATHOGENIC")


@dataclass(frozen=True)
class SameResidueEvidence:
    ps1_precedent_established: Optional[bool] = None
    ps1_precedent_classification: Optional[str] = None
    ps1_precedent_variant: Optional[str] = None
    pm5_precedent_established: Optional[bool] = None
    pm5_precedent_classification: Optional[str] = None
    pm5_precedent_variant: Optional[str] = None
    splice_impact_excluded: Optional[bool] = None

    def __post_init__(self) -> None:
        context = "SameResidueEvidence"
        self._validate_precedent(context, "ps1", self.ps1_precedent_established, self.ps1_precedent_classification, self.ps1_precedent_variant)
        self._validate_precedent(context, "pm5", self.pm5_precedent_established, self.pm5_precedent_classification, self.pm5_precedent_variant)
        if (self.ps1_precedent_established or self.pm5_precedent_established) and self.splice_impact_excluded is None:
            raise SchemaValidationError(
                f"{context}: splice_impact_excluded must be explicitly true/false whenever a PS1 or "
                "PM5 precedent is established=True -- PS1/PM5 cannot be safely evaluated without "
                "excluding a splice-driven mechanism (Richards et al. 2015 caveat; Walker et al. 2023, "
                "PMID 37352859). Never left unstated, same convention as nmd_predicted/repeat_region."
            )

    @staticmethod
    def _validate_precedent(context: str, code: str, established: Optional[bool], classification: Optional[str], variant: Optional[str]) -> None:
        if established is True:
            if classification not in _VALID_PRECEDENT_CLASSIFICATIONS:
                raise SchemaValidationError(
                    f"{context}: {code}_precedent_established=True requires {code}_precedent_classification "
                    f"to be one of {_VALID_PRECEDENT_CLASSIFICATIONS}, got {classification!r}"
                )
            if not variant or not variant.strip():
                raise SchemaValidationError(
                    f"{context}: {code}_precedent_established=True requires a non-empty "
                    f"{code}_precedent_variant citing the precedent (never assert a precedent without citing it)"
                )
        elif classification is not None or variant is not None:
            raise SchemaValidationError(
                f"{context}: {code}_precedent_classification/{code}_precedent_variant must not be set "
                f"unless {code}_precedent_established=True"
            )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "SameResidueEvidence":
        data = require_dict(data, context or "SameResidueEvidence")
        ctx = context or "SameResidueEvidence"
        return cls(
            ps1_precedent_established=optional_bool(data, "ps1_precedent_established", ctx),
            ps1_precedent_classification=optional_str(data, "ps1_precedent_classification"),
            ps1_precedent_variant=optional_str(data, "ps1_precedent_variant"),
            pm5_precedent_established=optional_bool(data, "pm5_precedent_established", ctx),
            pm5_precedent_classification=optional_str(data, "pm5_precedent_classification"),
            pm5_precedent_variant=optional_str(data, "pm5_precedent_variant"),
            splice_impact_excluded=optional_bool(data, "splice_impact_excluded", ctx),
        )
