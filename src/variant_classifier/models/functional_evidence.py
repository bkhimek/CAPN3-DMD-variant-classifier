"""FunctionalEvidence -- curated evidence for PS3 ("well-established
in vitro or in vivo functional studies supportive of a damaging effect on
the gene or gene product") and BS3 ("well-established in vitro or in vivo
functional studies show no damaging effect on protein function or
splicing"), Richards et al. 2015, Table 3. Added batch 25 alongside
evaluators/ps3.py and evaluators/bs3.py.

Base Richards et al. 2015 gives PS3/BS3 no further structure: "well-
established" is left to the curator's judgment, both criteria default to
Strong, and nothing distinguishes a rigorously controlled, statistically
validated assay from a single uncontrolled observation. The ClinGen SVI
Working Group's own refinement -- Brnich et al. 2019, "Recommendations
for application of the functional evidence PS3/BS3 criterion using the
ACMG/AMP sequence variant interpretation framework," Genome Medicine
11:98 -- replaces that judgment call with an explicit validation-tier
ladder this model follows directly:

  - Supporting: an assay with experimental controls and replicates but
    <=10 validation controls (a mix of known pathogenic/benign variants
    used to show the assay discriminates), OR a historically/broadly
    accepted assay class or commercial kit where controls aren't
    documented for this specific instance.
  - Moderate: >=11 total validation controls (a mix of pathogenic and
    benign), but no formal statistical analysis of discriminating power.
  - Strong (the base framework's ceiling for PS3/BS3): a rigorous
    statistical analysis produces a formal OddsPath value, with strength
    scaling to how extreme that OddsPath is (Brnich et al. 2019, Table 3).

This project does not compute OddsPath itself -- validation_strength is a
curated fact (the calibrated tier a curator/VCEP already assigned),
exactly like ComputationalEvidence's "one calibrated call per variant,
not several correlated tool votes" and SameResidueEvidence's precedent
classification. Only SUPPORTING, MODERATE, and STRONG are valid here --
Richards et al. 2015 never defines a Very-Strong or Stand-Alone tier for
PS3/BS3, and no VCEP specification this project has adopted (CAPN3's
LGMD VCEP) extends it further, so accepting those values would silently
overclaim a strength the base framework doesn't offer.

Assay direction and mechanism relevance are also curated facts, not
derived: `assay_result` states whether the assay showed a damaging
(ABNORMAL), non-damaging (NORMAL), or non-discriminating (INDETERMINATE)
effect -- INDETERMINATE is a real, distinct state from "no functional
data at all" (functional_evidence absent from the bundle entirely),
covering the common real case of an assay that was performed but did not
clearly distinguish pathogenic from benign for this specific variant
(see CAPN3_c.2257G>A's real fixture, batch 25 -- a directly-relevant
real-world illustration of exactly this Brnich caveat: Western blot
calpain-3 expression is well documented to not reliably track
pathogenicity for this gene, and this project's own curated notes quote
the primary literature saying so directly).

Per Brnich et al. 2019 Recommendation 6 (conflicting assay results):
"If the results are conflicting, the assay that most closely reflects
the disease mechanism and is more well-validated can override... If the
assays are essentially at the same level of validation, conflicting
functional evidence should not be used." This project represents ONE
already-reconciled FunctionalEvidence record per variant -- exactly like
ComputationalEvidence -- on the assumption that a curator has already
applied this reconciliation before curating a fixture, not something the
evaluator re-derives from multiple raw assay observations.
"""

from dataclasses import dataclass
from typing import Optional

from ..errors import SchemaValidationError
from ._coerce import coerce_enum, optional_str, require_dict
from .enums import CriterionStrength
from .enums import FunctionalAssayResult

_VALID_FUNCTIONAL_STRENGTHS = (CriterionStrength.SUPPORTING, CriterionStrength.MODERATE, CriterionStrength.STRONG)


@dataclass(frozen=True)
class FunctionalEvidence:
    assay_result: FunctionalAssayResult
    validation_strength: Optional[CriterionStrength] = None
    assay_description: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        context = "FunctionalEvidence"
        if not isinstance(self.assay_result, FunctionalAssayResult):
            raise SchemaValidationError(f"{context}: assay_result must be a FunctionalAssayResult value")
        if self.assay_result == FunctionalAssayResult.INDETERMINATE:
            if self.validation_strength is not None:
                raise SchemaValidationError(
                    f"{context}: validation_strength must not be set when assay_result=INDETERMINATE -- "
                    "an indeterminate assay supports neither PS3 nor BS3 at any strength"
                )
        else:
            if self.validation_strength is None:
                raise SchemaValidationError(
                    f"{context}: validation_strength must be explicitly stated whenever assay_result is "
                    "ABNORMAL or NORMAL -- PS3/BS3 cannot be safely evaluated from an unstated validation "
                    "tier (Brnich et al. 2019). Never left unstated, same 'never silently guess' convention "
                    "as nmd_predicted/repeat_region/splice_impact_excluded."
                )
            if self.validation_strength not in _VALID_FUNCTIONAL_STRENGTHS:
                raise SchemaValidationError(
                    f"{context}: validation_strength={self.validation_strength.value} is not valid for "
                    f"PS3/BS3 -- expected one of {[s.value for s in _VALID_FUNCTIONAL_STRENGTHS]} "
                    "(Richards et al. 2015 defines no Very-Strong or Stand-Alone tier for PS3/BS3)"
                )

    @classmethod
    def from_dict(cls, data: dict, context: Optional[str] = None) -> "FunctionalEvidence":
        data = require_dict(data, context or "FunctionalEvidence")
        ctx = context or "FunctionalEvidence"
        assay_result = coerce_enum(FunctionalAssayResult, data.get("assay_result"), "assay_result", ctx)
        strength_raw = data.get("validation_strength")
        validation_strength = (
            coerce_enum(CriterionStrength, strength_raw, "validation_strength", ctx)
            if strength_raw is not None
            else None
        )
        return cls(
            assay_result=assay_result,
            validation_strength=validation_strength,
            assay_description=optional_str(data, "assay_description"),
            notes=optional_str(data, "notes"),
        )
