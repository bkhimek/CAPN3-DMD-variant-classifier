"""Milestone 2: the first criterion evaluators.

An evaluator takes a VariantEvidenceBundle plus whatever configuration it
needs, and returns a CriterionResult for one ACMG/AMP code. Nothing here
combines multiple criteria into a classification yet — that's Milestone 3.
"""

from .pm2 import evaluate_pm2
from .pvs1 import evaluate_pvs1

__all__ = ["evaluate_pm2", "evaluate_pvs1"]
