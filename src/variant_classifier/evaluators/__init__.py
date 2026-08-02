"""Milestone 2/3 (plus PM4, added batch 14; plus PS1/PM5, added batch 22;
plus PS3/BS3, added batch 25; plus PM3, added batch 28): the criterion
evaluators.

An evaluator takes a VariantEvidenceBundle plus whatever configuration it
needs, and returns a CriterionResult for one ACMG/AMP code. These twelve
cover SUPPORTED_CRITERIA_MILESTONE_1 plus PM4, PS1, PM5, PS3, BS3, and PM3
(see models/enums.py) -- combining their results into a classification is
engine.py, not this package.
"""

from .ba1 import evaluate_ba1
from .bp4 import evaluate_bp4
from .bs1 import evaluate_bs1
from .bs3 import evaluate_bs3
from .pm2 import evaluate_pm2
from .pm3 import evaluate_pm3
from .pm4 import evaluate_pm4
from .pm5 import evaluate_pm5
from .pp3 import evaluate_pp3
from .ps1 import evaluate_ps1
from .ps3 import evaluate_ps3
from .pvs1 import evaluate_pvs1

__all__ = [
    "evaluate_ba1",
    "evaluate_bp4",
    "evaluate_bs1",
    "evaluate_bs3",
    "evaluate_pm2",
    "evaluate_pm3",
    "evaluate_pm4",
    "evaluate_pm5",
    "evaluate_pp3",
    "evaluate_ps1",
    "evaluate_ps3",
    "evaluate_pvs1",
]
