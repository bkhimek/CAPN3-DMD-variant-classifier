"""Milestone 2/3: the criterion evaluators.

An evaluator takes a VariantEvidenceBundle plus whatever configuration it
needs, and returns a CriterionResult for one ACMG/AMP code. These six
cover SUPPORTED_CRITERIA_MILESTONE_1 (see models/enums.py) — combining
their results into a classification is engine.py, not this package.
"""

from .ba1 import evaluate_ba1
from .bp4 import evaluate_bp4
from .bs1 import evaluate_bs1
from .pm2 import evaluate_pm2
from .pp3 import evaluate_pp3
from .pvs1 import evaluate_pvs1

__all__ = ["evaluate_ba1", "evaluate_bp4", "evaluate_bs1", "evaluate_pm2", "evaluate_pp3", "evaluate_pvs1"]
