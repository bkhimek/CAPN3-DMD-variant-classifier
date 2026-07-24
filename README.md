# CAPN3/DMD Variant Classifier

A small, scoped prototype implementing the ACMG Engine architecture
described in the companion design-guide set, starting with two genes:
**CAPN3** (autosomal recessive, LGMDR1/calpainopathy) and **DMD**
(X-linked, out of schema scope until Milestone 4 — see Roadmap).

## Status: Milestone 4 complete, curated variant set expanding toward 20-30 (10 done)

Milestone 1 built the schema and fixtures. Milestone 2 added the first two
evaluators (PM2, PVS1). Milestone 3 added the remaining four (BA1, BS1,
PP3, BP4) and a combining engine (`classify()`) that turns a full set of
evaluated criteria into a classification — Pathogenic / Likely Pathogenic
/ VUS / Likely Benign / Benign, for one variant at a time.

Milestone 4 adds something qualitatively different: **case-level**
reasoning. Every prior milestone answers "how strong is the evidence for
this one variant." Milestone 4 answers "does what was found in this
patient explain their disease" — which depends on how many variants they
have and how the gene's disease is inherited, not on any single variant in
isolation. `clinical.interpret_case(case, classifications, gene_disease_context)`
now goes from a `ClinicalCase` (what was found in a specific patient) plus
each variant's own classification to EXPLAINED / INSUFFICIENT /
MANUAL_REVIEW / NOT_APPLICABLE. What exists:

- Seven typed data models (`src/variant_classifier/models/`) matching the
  schemas in the *Building an ACMG Engine* and *Clinical Variant Pipeline
  Workflow Architecture* design guides, each validating its own invariants
  and rejecting malformed input with a single `SchemaValidationError`.
- **Ten curated evidence bundles** (`data/curated/variant_evidence.json`):
  six real ClinVar-grounded variants (five CAPN3, one DMD) and four
  synthetic cases (three CAPN3, one DMD) constructed to exercise specific
  combining-rule and case-level paths. This is two increments toward the
  ~20-30 ClinVar variant set from the original project plan — see
  "Expanding the curated set" below for what each real variant adds and
  where this is headed. Gene/disease context for both CAPN3 and DMD
  (`data/curated/gene_disease_context.yaml`).
- Golden cases for every curated fixture, curated *separately* from the
  evidence they judge, per the golden-case philosophy in the Validation
  and Verification design guide —
  `validation/golden_cases/variant_golden_cases.yaml` for per-variant
  results and `case_interpretation_golden_cases.yaml` for case-level
  results (see below).
- Schema-validation tests (`tests/unit/`) covering both valid and
  invalid records for every model.
- Six evaluators, one per Milestone-1-scope criterion
  (`src/variant_classifier/evaluators/`): `pvs1.py`, `pm2.py`, `pp3.py`,
  `bp4.py`, `ba1.py`, `bs1.py`. PVS1 is deliberately partial — see "PVS1
  scope" below.
- A **combining engine** (`src/variant_classifier/engine.py`) implementing
  the ACMG/AMP combining rules (Richards et al. 2015, Table 5), including
  a genuine conflict path (`conflicting_evidence_flag`) rather than
  silently picking a side when evidence satisfies both a pathogenic and a
  benign combining rule at once.
- Two new Milestone 4 models — `ClinicalCase` (what was found in one
  patient: gene, karyotypic sex, one or two variant_ids, phase between
  them if two) and `CaseInterpretation` (the case-level verdict) — plus
  **`src/variant_classifier/clinical.py`**, which reasons about autosomal
  recessive cases (trans/cis/unknown phase, single-variant insufficiency)
  and X-linked cases (hemizygous male vs everything else) separately. See
  "Case-level scope" below for exactly what is and isn't covered.
- Ten curated variants (CAPN3 and DMD) and six curated `ClinicalCase`
  fixtures covering every branch above. Every evaluator, the combining
  engine, and the case interpretation layer are each verified against
  golden cases written independently of the code — including, for the
  case-level tests, that trans vs cis and male vs female change the
  outcome despite everything else being identical
  (`tests/unit/test_*.py`). 103 tests pass in total (the "matches golden
  case" tests iterate over however many fixtures exist rather than a
  hardcoded count, so this number grows automatically as the curated set
  does).

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

**BA1/BS1 mirror PM2's founder-frequency handling, deliberately.**
While building BS1, the same founder-enrichment ambiguity that makes PM2
return MANUAL_REVIEW for `CAPN3_c.550del` turned out to apply to BS1 too:
a variant locally common in one ancestry (0.75%) but rare overall (0.023%)
is exactly as ambiguous for "is this too common to be pathogenic" as it is
for "is this rare enough to be pathogenic." The Milestone 1 golden case
had BS1 as a confident NOT_MET, which would have been an inconsistency —
resolving the ambiguity automatically on the benign side while refusing to
on the pathogenic side, with no principled reason for the difference. The
golden case was corrected to BS1: MANUAL_REVIEW (see its curator_note for
the full explanation); the final classification for that variant doesn't
change (still VUS), but two criteria are now honestly flagged for human
judgment instead of one. BA1 got the same three-way branch for
consistency, even though its 5% threshold makes the branch unlikely to
matter in practice.

**Expanding the curated set.** Three real, ClinVar-documented CAPN3
variants were added alongside the original `CAPN3_c.550del`, chosen
deliberately for diversity rather than for easy wins:

- `CAPN3_c.1939G>T` (p.Glu647Ter) — a clean nonsense variant, no
  ancestry-enrichment complication. Reaches LIKELY_PATHOGENIC, not
  PATHOGENIC, because no computational evidence was gathered for it — an
  honest limit of the current criterion set, not a bug.
- `CAPN3_c.2257G>A` (p.Asp753Asn) — a real, still-debated variant.
  Population frequency clears this project's BS1 threshold, but BS1 alone
  isn't sufficient for LIKELY_BENIGN per Table 5, so the engine lands on
  VUS. Recent literature suggests this variant might reflect a
  low-penetrance *dominant* mechanism distinct from CAPN3's usual
  recessive model — something this project's one-inheritance-pattern-per-
  gene design can't represent, noted explicitly rather than silently
  mismodeled.
- `CAPN3_c.946-1G>A` — a real, experimentally-confirmed pathogenic splice
  acceptor variant. Because PVS1 doesn't auto-resolve splice variants (see
  "PVS1 scope" below), this variant lands on VUS despite being confidently
  Pathogenic in ClinVar with wet-lab support — the clearest real-world
  case yet of the PVS1 scope gap actually costing something.

A second round added two more, chosen to round out the picture rather
than repeat the same story:

- `DMD_c.2302C>T` (p.Arg768Ter) — the project's first real (non-synthetic)
  DMD variant, replacing the placeholder role `DMD_SYNTH_PATHOGENIC_01`
  had been playing alone. Same LIKELY_PATHOGENIC-not-PATHOGENIC story as
  `CAPN3_c.1939G>T`, for the same reason (no computational evidence
  gathered). Not yet wired into a `ClinicalCase` fixture — the Milestone 4
  case-level tests still use the synthetic DMD variant.
- `CAPN3_c.1343G>A` (p.Arg448His) — a real, extremely rare missense
  variant that ClinVar itself calls Uncertain Significance, and so does
  this engine (PM2 MET alone, no other criterion applies). Added
  deliberately because the three prior real variants were all "engine
  under-calls a real Pathogenic/Benign-leaning variant due to scope
  gaps" stories — worth having at least one case in the set where the
  limited criterion set's answer simply agrees with the real-world
  consensus, so the pattern of gaps doesn't look like the whole picture.

A true common (>5% gnomAD) CAPN3 variant to exercise BA1 MET against real
data was searched for but not found in this round — CAPN3's recessive
disease model means genuinely common variants in this gene are scarce in
the literature searched so far. Still a gap; a hand-built synthetic edge
case covers the BA1 MET logic branch in the meantime
(`tests/unit/test_ba1_bs1_evaluators.py`).

Reaching the full ~20-30 variant set is expected to take several more
rounds of this same process (research a real variant, ground its
evidence, hand-derive the expected result, verify against the engine).

**Case-level scope.** clinical.py deliberately does not extend the
per-variant evaluator pattern from Milestone 2/3 — PM3 ("detected in trans
with a pathogenic variant") has a structural circularity a per-variant
evaluator can't cleanly express (variant A's PM3 depends on variant B's
classification, which lives entirely outside variant A's own
VariantEvidenceBundle). Reasoning about it at the case level, after both
variants already have their own classification, sidesteps that rather than
working around it. Scope, stated plainly: autosomal recessive handles
exactly one or two variants, with phase required whenever there are two;
X-linked only handles a single variant, and only resolves the hemizygous
male case (karyotypic_sex=XY) confidently — any other karyotypic sex is
deferred to MANUAL_REVIEW rather than reasoned about, since female X-linked
carrier interpretation involves X-inactivation biology this project does
not model.

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
  engine.py                  evaluate_all() + combine() + classify() — the combining engine
  clinical.py                 interpret_case() — case-level reasoning, see "Case-level scope" above
  evaluators/
    pvs1.py                   evaluate_pvs1() — see "PVS1 scope" above
    pm2.py                    evaluate_pm2()
    pp3.py                    evaluate_pp3()
    bp4.py                    evaluate_bp4()
    ba1.py                    evaluate_ba1()
    bs1.py                    evaluate_bs1()

config/
  population_thresholds.yaml per-gene PM2 frequency thresholds (see Design notes)

data/
  curated/
    gene_disease_context.yaml   CAPN3 and DMD
    variant_evidence.json       10 curated variants (CAPN3 + DMD) -- growing toward ~20-30
    clinical_cases.json         6 curated ClinicalCase fixtures (Milestone 4)
  source/                    placeholder — raw pulls from ClinVar/gnomAD/VEP (empty)
  synthetic/                 placeholder — larger generated datasets (empty)

validation/golden_cases/
  variant_golden_cases.yaml            expected per-variant results (renamed from
                                        capn3_milestone1.yaml once DMD variants existed)
  case_interpretation_golden_cases.yaml expected case-level results (Milestone 4)

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
extra environment variables. All 103 tests currently pass.

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

- **Milestone 2** — done. PM2 and PVS1 evaluators (PVS1 intentionally
  partial).
- **Milestone 3** — done. BA1/BS1/PP3/BP4 evaluators and the combining
  engine.
- **Milestone 4** — done. `ClinicalCase`/`CaseInterpretation` models and
  `clinical.py`'s case-level reasoning (see "Case-level scope" above).
- Later: continue expanding curated fixtures toward the full 20-30
  ClinVar variant set (10 of ~20-30 done, see "Expanding the curated set"
  above; a true common/BA1-level CAPN3 variant is a known remaining gap);
  add PM3/PS1/PM5/PS3/BS3 as real per-variant criteria (distinct from how
  clinical.py currently handles trans/cis reasoning at the case level);
  revisit PVS1's partial scope (protein-domain criticality,
  constitutive-exon data); revisit DMD's CNV representation gap (see
  gene_disease_context.yaml) if real DMD variants are ever added; extend
  X-linked case interpretation beyond the hemizygous-male case once
  X-inactivation is worth modeling properly.
