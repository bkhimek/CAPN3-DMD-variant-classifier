# CAPN3/DMD Variant Classifier — Project 1

A small, scoped prototype implementing the ACMG Engine architecture
described in the companion design-guide set, starting with two genes:
**CAPN3** (autosomal recessive, LGMDR1/calpainopathy) and **DMD**
(X-linked, out of schema scope until Milestone 4 — see below).

This is Krzysztof's first human-genomics project, following several
bacterial-genomics projects. It's meant to be pushed to
[github.com/bkhimek](https://github.com/bkhimek) as a new repo — see
"Moving this into your WSL projects folder" at the bottom.

## What Milestone 1 actually is

Milestone 1 is schema and fixtures only. There is **no evaluator and no
combining engine yet** — nothing in this repo computes a classification
from raw evidence. What exists:

- Seven typed data models (`src/variant_classifier/models/`) matching the
  schemas in the *Building an ACMG Engine* and *Clinical Variant Pipeline
  Workflow Architecture* guides, each validating its own invariants and
  rejecting malformed input with a single `SchemaValidationError`.
- Three curated CAPN3 evidence bundles (`data/curated/variant_evidence.json`):
  one real ClinVar-grounded case and two synthetic cases built specifically
  to exercise combining-rule paths the real case doesn't reach.
  A CAPN3 gene/disease context (`data/curated/gene_disease_context.yaml`).
- Three golden cases (`validation/golden_cases/capn3_milestone1.yaml`) —
  expected results, hand-derived from the ACMG/AMP combining rules and
  curated *separately* from the evidence bundles, per the golden-case
  philosophy in the Validation and Verification guide. These are what a
  future evaluator will be checked against; nothing in this repo checks
  them today.
- 37 schema-validation tests (`tests/unit/`) proving the models accept
  valid records and reject the invalid ones they're supposed to.

Milestone 1 is deliberately "something running," not "something that
classifies variants." The next milestones build the parts that do that.

## Why this milestone's plan changed slightly while building it

Two real gaps turned up during construction, both fixed rather than
deferred, since they were cheap to fix now and expensive to fix once
fixtures depended on the old shape:

1. **The originally-proposed 5-criterion set (PVS1/PM2/PP3/BA1/BS1) could
   never reach LIKELY_BENIGN.** Per Richards et al. 2015 Table 5, Likely
   Benign requires 1 Strong-Benign + 1 Supporting-Benign, or 2
   Supporting-Benign — and that 5-criterion set has only BS1 (a single
   Strong-Benign criterion) with no Supporting-Benign criterion at all.
   Fixed by adding **BP4** to `SUPPORTED_CRITERIA_MILESTONE_1`. The
   `CAPN3_SYNTH_LIKELY_BENIGN_01` golden case exists specifically to prove
   this path now works.

2. **`PopulationRetrievalStatus` had five states, all "missing/negative"
   outcomes (`ABSENT`, `NOT_ASSESSED`, `UNAVAILABLE`, `NOT_APPLICABLE`,
   `UNKNOWN`) — none of them meant "successfully found, here's the AF."**
   That gap only became visible when writing the real CAPN3 fixture, which
   has a genuine nonzero gnomAD frequency. Fixed by adding an `OBSERVED`
   state. `ComputationalEvidence` (a model not in the original six-model
   plan, needed because PP3/BP4 need *something* to evaluate from) reuses
   the same enum for the same reason.

Both are noted inline in the relevant source files, not just here.

## Why dataclasses instead of pydantic

The build environment had no PyPI access (`pip install pydantic` failed
with a proxy 403), so all seven models use the Python standard library's
`dataclasses` module with hand-written `from_dict()` validation instead.
This is a real constraint of *this build environment*, not a recommendation
for yours — `requirements.txt` lists `pydantic` nowhere because the models
don't need it to function, but if you'd rather have pydantic's validation
machinery going forward, converting these seven small files is a
contained, mechanical piece of work. Worth deciding before Milestone 2
adds more models, not after.

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

data/
  curated/                   hand-curated fixtures (the 3 CAPN3 cases live here)
  source/                    placeholder — raw pulls from ClinVar/gnomAD/VEP (empty)
  synthetic/                 placeholder — larger generated datasets (empty)

validation/golden_cases/     expected results, curated separately from data/curated/

config/                      placeholder — runtime thresholds etc. (empty)

tests/unit/                  pytest-convention tests
tests/run_tests.py           standalone runner (see below)
```

## Running the tests

This environment couldn't install `pytest` (no PyPI access), so the test
files are plain pytest-convention (`test_*.py`, `test_` functions, bare
`assert`) but were verified with a small dependency-free runner instead:

```bash
cd CAPN3-DMD-variant-classifier
PYTHONPATH=src python3 tests/run_tests.py
```

This prints PASS/FAIL per test and a summary; all 37 currently pass.

Once you have real `pytest` available (you will, on your own machine —
this was purely a sandbox limitation):

```bash
pip install -r requirements.txt
pytest
```

No changes to the test files are needed either way.

You can also load and sanity-check the fixtures directly:

```bash
PYTHONPATH=src python3 -m variant_classifier.loader
```

This prints how many gene/disease contexts, evidence bundles, and golden
cases loaded, any rejected records, and any cross-check warnings (e.g. a
golden case with no matching evidence bundle).

## Roadmap (not started)

- **Milestone 2** — first evaluators: PM2, then PVS1.
- **Milestone 3** — combination engine (limited scope: the six
  Milestone-1 criteria only).
- **Milestone 4** — clinical interpretation layer: CAPN3 recessive
  allele-count handling, DMD X-linked hemizygous handling.
- Later: expand curated fixtures to the full 20–30 ClinVar variant set;
  add PM3/PS1/PM5/PS3/BS3.

## Moving this into your WSL projects folder

From WSL, with this folder's contents copied somewhere accessible (e.g.
your Windows Downloads, then into WSL):

```bash
cp -r /mnt/c/Users/<you>/Downloads/CAPN3-DMD-variant-classifier ~/projects/
cd ~/projects/CAPN3-DMD-variant-classifier
git init
git add .
git commit -m "Milestone 1: repo scaffold, 7 data models, curated CAPN3 fixtures, golden cases"
git branch -M main
git remote add origin https://github.com/bkhimek/CAPN3-DMD-variant-classifier.git
git push -u origin main
```

(Create the empty `CAPN3-DMD-variant-classifier` repo on
[github.com/bkhimek](https://github.com/bkhimek) first if it doesn't
already exist — GitHub won't auto-create it on push.)
