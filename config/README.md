# config/

Placeholder for runtime configuration once there is any: population
frequency thresholds per gene/disease, tool version pins for computational
evidence sources, combining-rule specification versions to use by default,
etc. — the kind of thing the ACMG Engine Detailed Design Guide says should
be data, not hard-coded constants in the evaluators.

Empty in Milestone 1: there are no evaluators yet to configure. Milestone 2
(PM2 evaluator) is expected to be the first thing that needs a file here,
most likely a per-gene population frequency threshold table.
