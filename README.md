# CAPN3/DMD Variant Classifier

A small, scoped prototype implementing the ACMG Engine architecture
described in the companion design-guide set, starting with two genes:
**CAPN3** (autosomal recessive, LGMDR1/calpainopathy) and **DMD**
(X-linked, out of schema scope until Milestone 4 — see Roadmap).

## Status: Milestone 2 complete — PM2 and PVS1 evaluators

Milestone 1 built the schema and fixtures; there was no evaluator and no
combining engine. Milestone 2 adds evaluators — code that reads a
VariantEvidenceBundle and decides whether one ACMG/AMP criterion is MET.
**PM2 and PVS1 are done; a combining engine that turns MET criteria into a
classification is not (Milestone 3).** What exists:

- Seven typed data models (`src/variant_classifier/models/`) matching the
  schemas in the *Building an ACMG Engine* and *Clinical Variant Pipeline
  Workflow Architecture* design guides, each validating its own invariants
  and rejecting malformed input with a single `SchemaValidationError`.
- Three curated CAPN3 evidence bundles (`data/curated/variant_evidence.json`):
  one real ClinVar-grounded case and two synthetic cases constructed to
  exercise combining-rule paths the real case doesn't reach. A CAPN3
  gene/disease context (`data/curated/gene_disease_context.yaml`).
- Three golden cases (`validation/golden_cases/capn3_milestone1.yaml`) —
  expected results, hand-derived from the ACMG/AMP combining rules and
  curated *separately* from the evidence bundles, per the golden-case
  philosophy in the Validation and Verification design guide. These are
  what a future evaluator will be checked against; nothing in this repo
  checks them yet.
- 37 schema-validation tests (`tests/unit/`) covering both valid and
  invalid records for each model.
- A **PM2 evaluator** (`src/variant_classifier/evaluators/pm2.py`), driven
  by a per-gene frequency threshold in `config/population_thresholds.yaml`.
- A **PVS1 evaluator** (`src/variant_classifier/evaluators/pvs1.py`),
  deliberately partial — see "PVS1 scope" below.
- Both evaluators verified against all three curated CAPN3 cases, matching
  their golden cases exactly, plus edge-case tests for branches the three
  curated cases don't happen to cover
  (`tests/unit/test_pm2_evaluator.py`, `tests/unit/test_pvs1_evaluator.py`).
  59 tests pass in total.

## Design notes

**Criterion coverage.** Per Richards et al. 2015 (Table 5), reaching
LIKELY_BENIGN requires either 1 Strong-Benign + 1 Supporting-Benign
criterion, or 2 Supporting-Benign criteria. `SUPPORTED_CRITERIA_MILESTONE_1`
therefore includes both **BS1** (Strong-Benign) and **BP4**
(Supporting-Benign) — BS1 alone can never combine to LIKELY_BENIGN. The
`CAPN3_SYNTH_LIKELY_BENIGN_01` golden case exercises this path directly.

**Population retrieval status.** `PopulationRetrievalStatus` distinguishes
several "missing/negative" outcomes (`ABSENT`, `NOT_ASSESSED`,
`UNAVAILABLE`, `NOT_APPLICABLE`, `UNKNOWN`) from `OBSERVED` — a variant
successfully retrieved from a source with real frequency data attached.
Collapsing "not found" and "not assessed" into one state would silently
treat missing evidence as negative evidence, which the Reporting and
Dashboard design guide specifically warns against. `ComputationalEvidence`
(the model backing PP3/BP4) reuses the same enum for the same reason.

**PVS1 scope.** The full PVS1 decision tree (Abou Tayoun et al. 2018)
branches on protein-domain criticality and constitutive-exon-splicing
information this project doesn't model. This evaluator only ever returns
MET for the one case it can defend end-to-end: an early frameshift or
nonsense variant predicted to trigger nonsense-mediated decay, in a gene
with an established loss-of-function mechanism. Everything harder —
truncations that escape NMD (typically last-exon), splice donor/acceptor
variants, and start-loss variants — returns MANUAL_REVIEW with a rationale
explaining why, rather than a guessed MET or NOT_MET. Non-null-variant
consequence types (missense, synonymous, etc.) return NOT_APPLICABLE.
`TranscriptConsequence` requires an explicit `nmd_predicted` value for
both frameshift and stop-gained variants for exactly this reason — this
requirement was originally frameshift-only in Milestone 1 and widened here
once the evaluator needed it for stop-gained variants too.

**PM2 and founder mutations.** PM2 asks whether a variant is absent or at
extremely low frequency in the general population. A single global allele
frequency threshold isn't enough to answer that safely: `CAPN3_c.550del` is
rare overall (0.023%) but a known founder mutation enriched to 0.75% in
specific ancestries. The evaluator does not silently pass PM2 using the
lower, reassuring global number — when an ancestry-specific frequency
clears the threshold while the overall frequency doesn't, it returns
MANUAL_REVIEW rather than guessing, because whether "extremely low" holds
depends on the tested individual's ancestry, which isn't available here.

**Dataclasses instead of pydantic.** All seven models use the Python
standard library's `dataclasses` module with hand-written `from_dict()`
validation rather than pydantic. This keeps the dependency footprint to
just PyYAML for fixture loading. Converting to pydantic later, if its
validation machinery becomes useful, is a contained, mechanical change
scoped to these seven files.

## Repository layout

```
src/variant_classifier/
  errors.py                  SchemaValidationError — the one exception type
  models/
    enums.py                 controlled vocabularies + ACMG_CRITERION_CODES
    _coerce.py                shared from_dict() validation helpers
    variant_identity.py        VariantIdentity
    gene_disease_context.py    GeneDiseaseContext, Specification
    transcript_consequence.py  TranscriptConsequence
    population_evidence.py     PopulationEvidence
    computational_evidence.py  ComputationalEvidence
    criterion_result.py        CriterionResult
    provisional_classification.py  ProvisionalClassification
    evidence_bundle.py         VariantEvidenceBundle (container, this repo only)
    golden_case.py             GoldenCase (container, this repo only)
  loader.py                  loads/validates the curated fixtures below
  evaluators/
    pm2.py                    evaluate_pm2()
    pvs1.py                   evaluate_pvs1() — see "PVS1 scope" above

config/
  population_thresholds.yaml per-gene PM2 frequency thresholds (see Design notes)

data/
  curated/                   hand-curated fixtures (the 3 CAPN3 cases live here)
  source/                    placeholder — raw pulls from ClinVar/gnomAD/VEP (empty)
  synthetic/                 placeholder — larger generated datasets (empty)

validation/golden_cases/     expected results, curated separately from data/curated/

tests/unit/                  pytest tests
tests/run_tests.py           dependency-free runner (see below)
```

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running the tests

```bash
pytest
```

`pytest.ini` sets `pythonpath = src`, so this works out of the box with no
extra environment variables. All 37 tests currently pass.

A dependency-free alternative is also included, useful in environments
without PyPI access:

```bash
PYTHONPATH=src python3 tests/run_tests.py
```

It discovers `tests/unit/test_*.py`, runs every `test_*` function, and
prints a PASS/FAIL summary — the same tests, no `pytest` package required.

## Sanity-checking the fixtures

```bash
PYTHONPATH=src python3 -m variant_classifier.loader
```

Prints how many gene/disease contexts, evidence bundles, and golden cases
loaded, any rejected records, and any cross-check warnings (e.g. a golden
case with no matching evidence bundle).

## Roadmap

- **Milestone 2** — done. PM2 and PVS1 evaluators (see above; PVS1 is
  intentionally partial).
- **Milestone 3** — combination engine (limited scope: the six
  Milestone-1 criteria only). Also the natural point to revisit whether
  PP3/BP4 evaluators are needed, since the LIKELY_BENIGN and full
  PATHOGENIC combining-rule paths depend on them.
- **Milestone 4** — clinical interpretation layer: CAPN3 recessive
  allele-count handling, DMD X-linked hemizygous handling.
- Later: expand curated fixtures to the full 20–30 ClinVar variant set;
  add PM3/PS1/PM5/PS3/BS3.
