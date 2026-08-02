# CAPN3/DMD Variant Classifier

A small, scoped prototype implementing the ACMG Engine architecture
described in the companion design-guide set, starting with two genes:
**CAPN3** (autosomal recessive, LGMDR1/calpainopathy) and **DMD**
(X-linked, out of schema scope until Milestone 4 — see Roadmap).

## Status: COMPLETE (as of Batch 30). Milestone 5 complete, PS1/PM5 evaluators added (Batch 22), DMD CNV/structural-variant deletion scoring (Batch 23) and duplication scoring (Batch 24) both implemented as deliberately partial slices, PS3/BS3 functional-evidence evaluators added (Batch 25), PVS1 extended with a real splice-RNA-evidence branch (Batch 26) and a real start-loss alternative-start-codon branch (Batch 27), PM3 (in-trans-with-a-pathogenic-variant) implemented (Batch 28) via curated per-proband points -- this project's twelfth evaluator and first Strong-strength pathogenic-direction evaluator, and the change that finally moves the real `CAPN3_c.550del` founder-allele fixture to PATHOGENIC, matching its real ClinVar call -- curated variant set at 19 real of ~20-30 (29 point-mutation fixtures total, plus 3 CNV deletion + 2 CNV duplication fixtures in separate curated sets), `clinical.py` extended to handle biallelic XX X-linked case interpretation (Batch 29) -- a real, X-inactivation-independent mechanism, closing the last structural gap in case-level reasoning -- CAPN3-DMD-variant-calling-pipeline integration adapter done (Batch 21). Batch 30 closed the project out with a documentation-only pass (no code changes) -- see "Final status and scope boundary" at the end of the Roadmap section for the definitive summary of what's built and what's deliberately left out of scope.

Milestone 5 (batch 20) added a second combining system -- Bayesian
point-based combining (Tavtigian et al. 2020), offered alongside, not
replacing, Milestone 3's classic Table 5 rules -- and was expected to be
this project's last milestone before development moved to a new,
complementary project starting one stage earlier in the pipeline. That
move hadn't happened yet by batch 22, which returned to this project for
one more round: PS1 and PM5, the two criteria the Roadmap had already
flagged as "most tractable next," are now implemented (see "PS1 and PM5:
same-residue precedent evidence" below), and DMD's long-disclosed CNV/
structural-variant representation gap was researched and sized against
the real ClinGen technical standard for CNV interpretation -- deliberately
NOT implemented in batch 22, per its own sizing conclusion that it needs a
new identity representation and a new, parallel scoring system, not an
incremental extension of the existing evaluator pattern. Batch 23 picked
that sizing conclusion back up and built the scoped slice it identified:
a new `CnvDeletionEvidence` evidence family, a new `cnv_scoring.py`
combining module (parallel to `bayesian.py`, not part of it), and a
deletion-only implementation of Section 2 (dosage-sensitivity) of the
Riggs et al. 2020 CNV rubric, grounded in real ClinVar and literature DMD
deletion examples -- see "DMD CNV/structural-variant scoring (batch 23)"
below for the full writeup, including exactly what is and isn't
implemented. Batch 24 extended this to duplications (`CnvDuplicationEvidence`,
`score_cnv_duplication()`), but research surfaced a real wrinkle first:
ClinGen's own DMD dosage curation says whole-gene DMD duplications aren't
clinically reported, and the exact gain-side point values needed for the
real, clinically-relevant intragenic-disruption mechanism could only be
partially confirmed from secondary sources -- see "DMD CNV/structural-
variant duplication scoring (batch 24)" below for the full writeup,
including the specific, disclosed inference this project made to proceed
anyway (confirmed with the user before writing any code). Batch 25 added
two more point-mutation criteria, PS3 and BS3 (well-established functional
studies supportive of, or showing no, damaging effect) -- a new
`FunctionalEvidence` model, `evaluate_ps3()`/`evaluate_bs3()`, and two
real curated fixture enrichments, one showing PS3's MET path (a founder
CAPN3 allele with real Western-blot protein-absence data) and one showing
the framework's INDETERMINATE assay-result path (a real, genuinely mixed
Western-blot result for a separately-contested CAPN3 missense variant) --
see "PS3 and BS3: functional-evidence criteria (batch 25)" below. Batch
25 also researched PM3 (in trans with a second pathogenic variant) and
explicitly sized it as a real, larger architectural change requiring
multi-proband evidence aggregation, deliberately not implemented this
round -- see "PM3: sized, not implemented (batch 25)" below. Batch 26
extended PVS1 (see "PVS1 scope" and "Splice-RNA evidence feeds PVS1
directly (batch 26)" below): a new `splicing_rna_evidence` field on
`TranscriptConsequence` lets a real RNA/splicing assay result resolve
splice donor/acceptor variants directly, mirroring the real ClinGen LGMD
VCEP CAPN3 specification's own instruction that experimental splicing
evidence should be scored under PVS1, not PS3. The exact percentage-based
thresholds from the full Abou Tayoun et al. 2018 / Walker et al. 2023
decision trees remain unimplemented -- both the primary paper and the
CAPN3-specific PVS1 flowchart PDF were unreachable during this batch's
research -- so only the two threshold-free branches (a confirmed null-
equivalent transcript, or confirmed normal splicing directly contradicting
the predicted disruption) are implemented; none of the three real curated
splice fixtures had precise enough primary data to exercise the new
MET path, which is instead demonstrated on one disclosed synthetic
fixture. Batch 27 continued PVS1's start-loss branch: this time the
governing rule (ACGS 2024 UK Practice Guidelines, quoting Abou Tayoun et
al. 2018's decision tree for initiation-codon variants) WAS fully
reachable, and the real alternative-start-codon question for CAPN3's own
start-loss fixture (`CAPN3_c.1A>G`) was independently answered from
primary RefSeq/Ensembl CDS sequence data plus this project's own curated
fixture set -- see "PVS1 start-loss: the alternative-start-codon rule
(batch 27)" below for the full writeup. Batch 28 came back to PM3, sized
but deliberately left unimplemented in batch 25, and built it: a new
`Pm3Evidence`/`Pm3ProbandObservation` model family (curated per-proband
points, resolving the same-engine circularity `clinical.py` had flagged
by treating the partner allele's classification as a known fact rather
than something re-derived live) and a new `evaluate_pm3()` -- this
project's twelfth evaluator, and its first Strong-strength pathogenic-
direction one. The real founder-allele fixture, `CAPN3_c.550del`, gained
two independent real published-cohort homozygous PM3 observations and now
classifies PATHOGENIC (PVS1 Very Strong + PM3 Strong), finally matching
its real ClinVar call and dropping out of the Table-5-vs-Bayesian
divergence set -- see "PM3: implemented (batch 28)" below for the full
writeup, including exactly which primary source was unreachable (again)
and how that shaped the design. Batch 29 returned to `clinical.py`,
untouched since Milestone 4: `interpret_x_linked_case` handled only the
hemizygous-male case, deferring every other karyotypic sex to
MANUAL_REVIEW since female X-linked carrier phenotype depends on
X-inactivation biology this project found (via real, cited literature)
the field itself cannot reliably predict. Batch 29's research found a
second, genuinely different, real mechanism this project's existing
model could newly represent without touching that finding: biallelic XX
X-linked involvement (both copies affected), which is X-inactivation-
*independent* and resolves confidently to EXPLAINED, while every other
XX combination -- including cis, which is deliberately NOT treated the
way autosomal recessive cis is -- correctly stays MANUAL_REVIEW for the
same real reason as before. See "X-linked female/other-karyotype case
interpretation (batch 29)" below for the full design and citations. See
"Milestone 5: Bayesian point-based combining" below for what that
milestone found and why batch 20 felt like the right place to pause
before batch 22/23/24/25/26/27/28/29 picked back up. Batch 30 is where
this project actually stops -- a deliberate, user-decided close rather
than another pause, documented (not just implied) in "Final status and
scope boundary" at the end of the Roadmap section.

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

- Eight typed data models (`src/variant_classifier/models/`) matching the
  schemas in the *Building an ACMG Engine* and *Clinical Variant Pipeline
  Workflow Architecture* design guides, each validating its own invariants
  and rejecting malformed input with a single `SchemaValidationError`.
  `SameResidueEvidence` (batch 22) was added for PS1/PM5 -- see "PS1 and
  PM5: same-residue precedent evidence" below. Batch 23 adds three more,
  in a deliberately separate family for CNV evidence rather than an
  extension of the point-mutation models above: `CnvDeletionEvidence`,
  `CnvCategoryResult`, `CnvProvisionalClassification` -- see "DMD CNV/
  structural-variant scoring (batch 23)" below. Batch 25 adds one more to
  the point-mutation family alongside `SameResidueEvidence`:
  `FunctionalEvidence`, attached to `VariantEvidenceBundle` for PS3/BS3 --
  see "PS3 and BS3: functional-evidence criteria (batch 25)" below. Batch
  28 adds two more, also attached to `VariantEvidenceBundle`:
  `Pm3Evidence` and `Pm3ProbandObservation`, for PM3 -- see "PM3:
  implemented (batch 28)" below.
- **Twenty-nine curated evidence bundles** (`data/curated/variant_evidence.json`):
  nineteen real ClinVar/VCEP-grounded variants (twelve CAPN3, seven DMD)
  and ten synthetic cases (nine CAPN3, one DMD -- `CAPN3_SYNTH_PS1_01`
  added batch 22 to exercise PS1's MET path, `CAPN3_SYNTH_PVS1_SPLICE_RNA_01`
  added batch 26 for the PVS1 splice-RNA-evidence MET path,
  `CAPN3_SYNTH_PVS1_STARTLOSS_NO_ALT_01`/`CAPN3_SYNTH_PVS1_STARTLOSS_SUPPORTING_01`
  added batch 27 for PVS1's two new start-loss MET branches,
  `CAPN3_SYNTH_PM3_MODERATE_01`/`CAPN3_SYNTH_PM3_CIS_OVERRIDE_01` added
  batch 28 for PM3's compound-heterozygous branch and cis-cooccurrence
  override, see "Expanding the curated set" below) constructed to
  exercise specific combining-rule and case-level paths. This is nineteen
  increments toward the ~20-30 ClinVar variant set from the original
  project plan — see "Expanding the curated set" below for what each real
  variant adds and where this is headed. Gene/disease context for both
  CAPN3 and DMD (`data/curated/gene_disease_context.yaml`).
- Golden cases for every curated fixture, curated *separately* from the
  evidence they judge, per the golden-case philosophy in the Validation
  and Verification design guide —
  `validation/golden_cases/variant_golden_cases.yaml` for per-variant
  results and `case_interpretation_golden_cases.yaml` for case-level
  results (see below).
- Schema-validation tests (`tests/unit/`) covering both valid and
  invalid records for every model.
- Twelve evaluators, one per implemented ACMG/AMP point-mutation criterion
  (`src/variant_classifier/evaluators/`): `pvs1.py`, `pm2.py`, `pm4.py`,
  `ps1.py`, `pm5.py`, `pm3.py`, `ps3.py`, `pp3.py`, `bp4.py`, `ba1.py`,
  `bs1.py`, `bs3.py`. PVS1 is deliberately partial — see "PVS1 scope"
  below. PM4 was added in batch 14 — see "PM4: a second new criterion"
  below. PS1 and PM5 were added in batch 22 — see "PS1 and PM5:
  same-residue precedent evidence" below. PS3 and BS3 were added in batch
  25 — see "PS3 and BS3: functional-evidence criteria (batch 25)" below.
  PM3 was added in batch 28, this project's first Strong-strength
  pathogenic-direction evaluator — see "PM3: implemented (batch 28)"
  below. Batch 23's CNV deletion scoring and batch 24's CNV duplication
  scoring (both in `cnv_scoring.py`) are a separate, parallel scoring
  module over the ACMG/ClinGen CNV rubric (Riggs et al. 2020), not one
  more evaluator feeding this list -- see "DMD CNV/structural-variant
  scoring (batch 23)" and "DMD CNV/structural-variant duplication scoring
  (batch 24)" below.
- A **combining engine** (`src/variant_classifier/engine.py`) implementing
  the ACMG/AMP combining rules (Richards et al. 2015, Table 5), including
  a genuine conflict path (`conflicting_evidence_flag`) rather than
  silently picking a side when evidence satisfies both a pathogenic and a
  benign combining rule at once.
- A **second combining system** (`src/variant_classifier/bayesian.py`),
  added Milestone 5 (batch 20): Bayesian point-based combining (Tavtigian
  et al. 2020), run over the exact same `evaluate_all()` output as Table
  5, offered alongside it rather than replacing it. See "Milestone 5:
  Bayesian point-based combining" below.
- Two new Milestone 4 models — `ClinicalCase` (what was found in one
  patient: gene, karyotypic sex, one or two variant_ids, phase between
  them if two) and `CaseInterpretation` (the case-level verdict) — plus
  **`src/variant_classifier/clinical.py`**, which reasons about autosomal
  recessive cases (trans/cis/unknown phase, single-variant insufficiency)
  and X-linked cases (hemizygous male, biallelic XX as of batch 29, and
  everything else) separately. See "Case-level scope" below for exactly
  what is and isn't covered.
- Twenty-nine curated variants (CAPN3 and DMD) and eleven curated
  `ClinicalCase` fixtures covering every branch above — including, as of
  batch 12, a pair built on a real ClinVar-sourced DMD variant rather
  than only the original synthetic ones, as of batch 17, the case-
  level layer's previously-untested MANUAL_REVIEW catch-all branch (a
  hemizygous male whose only variant is VUS, not yet qualifying or
  benign), and, as of batch 29, a trans/cis pair on a biallelic XX
  X-linked case (`CASE_DMD_XX_BIALLELIC_TRANS`/`_CIS`) — see "X-linked
  female/other-karyotype case interpretation (batch 29)" below. Every
  evaluator, the combining
  engine, and the case interpretation layer are each verified against
  golden cases written independently of the code — including, for the
  case-level tests, that trans vs cis and male vs female change the
  outcome despite everything else being identical
  (`tests/unit/test_*.py`). 250 tests pass in total (the "matches golden
  case" tests iterate over however many fixtures exist rather than a
  hardcoded count, so this number grows automatically as the curated set
  does; also includes two hand-built tests pairing the real DMD
  male/female cases (batch 12), proving no real CAPN3 variant *without
  PM3 evidence* can currently reach Pathogenic/Likely Pathogenic *under
  Table 5* (batch 13, extended in batch 14 to cover PM4 too, given an
  explicit Bayesian counterpart in batch 20 rather than left to look
  contradicted by it, re-verified in batch 22 to still hold with PM5
  implemented and MET on a real fixture, and narrowed in batch 28 to
  explicitly exclude and separately verify `CAPN3_c.550del`, now a
  documented, intentional exception via real PM3 evidence), a dedicated
  `test_pm4_evaluator.py` suite added in batch 14, a dedicated
  `test_bayesian.py` suite (13 tests) added in batch 20, a dedicated
  `test_ps1_pm5_evaluators.py` suite (18 tests) added in batch 22, a
  dedicated `test_cnv_scoring.py` suite (20 tests in batch 23, growing to
  37 in batch 24) for the separate CNV deletion and duplication
  evidence/scoring families described below, a dedicated
  `test_ps3_bs3_evaluators.py` suite (15 tests) added in batch 25,
  `test_pvs1_evaluator.py` grew from 12 to 19 tests in batch 26 (new
  splice-RNA-evidence branch) and again to 27 tests in batch 27 (new
  start-loss alternative-start-codon branch), a dedicated
  `test_pm3_evaluator.py` suite (18 tests) added in batch 28, and 5 new
  `test_clinical.py` tests added in batch 29 for the biallelic-XX branch.

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

**Real ClinGen LGMD VCEP thresholds for CAPN3 (batch 4), and what adopting
them without a matching combining-system upgrade costs.** Through batch 3,
CAPN3's PM2/BA1/BS1 thresholds were hand-picked placeholders, explicitly
labeled as such. Batch 4 found the real, published ClinGen Limb-Girdle
Muscular Dystrophy (LGMD) Variant Curation Expert Panel specification for
CAPN3 (v2.0, released 2025-07-09, DOI 10.5281/zenodo.21434844) and adopted
its PM2/BA1/BS1 thresholds and PM2's strength into
`config/population_thresholds.yaml` and `gene_disease_context.yaml`
(CAPN3's specification is now `VCEP`, not `GENERIC_ACMG`; DMD, with no
adopted VCEP spec, is unaffected and stays `GENERIC_ACMG`). Three changes
worth naming explicitly:

- **PM2's threshold got much stricter** (0.01% instead of the placeholder
  0.1%), and its **strength changed from Moderate to Supporting** — the
  real VCEP only defines a Supporting-strength PM2 for CAPN3 at all.
- **BA1's threshold got much stricter for CAPN3** (0.3% instead of the
  generic ACMG/AMP 5% default) — a real, gene-specific override, now
  configurable per gene (`ba1_af` in `population_thresholds.yaml`) rather
  than hardcoded as one universal value.
- **BS1's threshold (0.1%) was unchanged** — it turned out the earlier
  placeholder had already guessed the right number, confirmed rather than
  corrected once the real spec was found.

This is not free: this project's combining engine (`engine.py`) still uses
the classic Richards et al. 2015 Table 5 categorical rules, but the real
VCEP spec pairs its thresholds with a **Bayesian point-based combining
system** instead (Tavtigian et al. 2020 adaptation — Very Strong=8,
Strong=4, Moderate=2, Supporting=1 point; Pathogenic >=10, Likely
Pathogenic 6-9, VUS 0-5, Likely Benign -6 to -1, Benign <=-7). Table 5 has
no rule at all for "1 Very Strong + 1 Supporting" — so a CAPN3 variant
whose only pathogenic evidence is PVS1 (Very Strong) plus PM2 (now
Supporting, not Moderate) no longer reaches any tier here, even though the
VCEP's own point system would score that exact evidence 8+1=9 points,
squarely in their Likely Pathogenic range. `CAPN3_c.1939G>T` — a real,
ClinVar-Pathogenic nonsense variant — shows this honestly: it dropped from
LIKELY_PATHOGENIC to VUS in this batch, left unmodified deliberately (see
its golden-case curator_note). A synthetic case-level fixture,
`CAPN3_SYNTH_PATHOGENIC_02`, hit the same gap and was deliberately
strengthened (a PP3 computational-evidence entry was added) to keep
demonstrating a clean case-level EXPLAINED result — a disclosed fixture
design choice, not a silent one, and a different treatment than the real
variant got on purpose. This is the strongest concrete argument this
project has produced so far for eventually adopting Bayesian point-based
combining rather than Table 5 (see Roadmap) — batch 20 (Milestone 5,
"Bayesian point-based combining" below) followed through on this: run
through `bayesian.py` instead of `engine.py`, `CAPN3_c.1939G>T`'s exact
same evidence (PVS1 Very Strong + PM2 Supporting) sums to 9 points,
Likely Pathogenic, closing this specific gap without touching Table 5's
own behavior at all.

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

A third round added two more:

- `CAPN3_c.2050+1G>A` — a real, Pathogenic splice DONOR variant,
  complementing `CAPN3_c.946-1G>A` (a splice ACCEPTOR variant). Same
  PVS1-scope-gap story, and having two independent real splice variants
  land the same way shows it's a systematic limitation, not a one-off.
- `CAPN3_c.1132T>C` — a real missense VUS where several individual
  in-silico tools (SIFT, PolyPhen-2, Align-GVGD) lean benign. Deliberately
  *not* wired into a `ComputationalEvidence` record: this project's
  PP3/BP4 model expects one calibrated ensemble score, not several raw
  per-tool votes stacked together, and retrofitting one from these would
  contradict that design principle. Lands on VUS, agreeing with ClinVar.

A fourth round (batch 5) added one more:

- `CAPN3_c.1A>G` (p.Met1Val) — a real, ClinVar-Pathogenic start-loss
  (initiator codon) variant, the project's first fixture exercising PVS1's
  third documented scope gap. PVS1 treats start-loss as within scope in
  principle but doesn't auto-resolve it (correctly doing so requires
  checking for a plausible downstream alternative start codon, which this
  project doesn't model), so it returns MANUAL_REVIEW here just as it does
  for the two splice-site fixtures. With PVS1 unresolved, PM2 alone (MET,
  Supporting under the real VCEP threshold) never reaches a Table 5 tier —
  VUS, the same systematic-limitation story as the splice-site fixtures,
  now shown a third independent way.

A fifth round (batch 6) added one more, sourced directly from an official
VCEP expert curation rather than an aggregate ClinVar record:

- `CAPN3_c.2120A>G` (p.Asp707Gly) — a real variant with an unusually rich
  paper trail: the ClinGen LGMD VCEP's own Evidence Repository curation
  (approved 2025-04-22) classifies it PATHOGENIC using PM3_Strong,
  PP1_Strong, and PP4 — all case-level/segregation criteria this project's
  variant-only engine doesn't implement. The one criterion the two share,
  PP3, uses this fixture's REVEL score of 0.966 — **the project's first
  real, precisely-sourced calibrated computational score**, closing a gap
  open since Milestone 1 (see below). More interestingly: the real VCEP
  found this variant's East Asian-specific frequency technically crosses
  their BS1 threshold, but explicitly declined to apply BS1 anyway (an
  documented "BS1 exception") because it's a well-established recurrent
  founder-pathogenic variant in that population, not a benign one. This
  project's BS1 evaluator can't know that context, but its
  founder-enrichment MANUAL_REVIEW branch (built in Milestone 3, see
  "BA1/BS1 mirror PM2's founder-frequency handling" above) flags the exact
  same case for human judgment anyway — a real-world confirmation that
  design decision catches what it was built to catch. Final result: VUS
  with manual review required, an honest reflection of how much of this
  variant's real classification lives outside a 6-criterion engine's
  reach.

Batch 8 didn't add a new fixture — instead it enriched an existing one
once new real evidence turned up, and made a sixth unsuccessful attempt
at the BA1 gap (ClinVar search for benign-classified CAPN3 variants, an
erepo.clinicalgenome.org affiliation-level browse, still nothing above
0.3%). While searching for a real BP4-benign example, `CAPN3_c.2257G>A`
(already in the set since batch 2) turned up two new pieces of real
information: its actual REVEL score (0.649 — added as computational
evidence, correctly landing as INDETERMINATE since it's below the real
VCEP's PP3 threshold but above its BP4 one, so both criteria now
correctly return NOT_MET instead of one of them being left NOT_EVALUATED)
and the real ClinGen LGMD VCEP's own current classification for this
exact variant: **Likely Benign, via BS1 + BP2** — a case-level criterion
this project doesn't implement at all. That's the benign-side mirror of
the `CAPN3_c.1939G>T` gap: a real VCEP reaching a confident classification
through evidence entirely outside this project's variant-only, Table-5
engine. (This particular finding came from a web search summary rather
than the primary erepo record, unlike `CAPN3_c.2120A>G` — flagged in the
fixture's notes as secondhand rather than primary-source-confirmed.)

A ninth round (batch 9) added one more, closing the last of PVS1's three
documented scope gaps with real data:

- `DMD_c.11041A>T` (p.Arg3681Ter) — a real ClinVar entry (VCV000641807.15,
  Uncertain significance, 4 concordant submitters), a nonsense variant in
  DMD's coding exon 78 of 79, the project's first REAL fixture exercising
  PVS1's NMD-escape branch (`nmd_predicted=false`) — previously only
  covered by a hand-built edge-case test. Coordinates weren't found by
  direct search; they were derived independently first (mapping DMD's
  MANE Select transcript's protein residue 3681 to genomic coordinates via
  the Ensembl REST API, retrieving the codon on both strands, and
  enumerating all nine possible single-nucleotide substitutions to find
  the one that produces a stop codon), then confirmed against the primary
  ClinVar record afterward, which returned the identical HGVS —
  independent derivation matching the authoritative source exactly, not
  built from it. `nmd_predicted=false` is itself independently
  corroborated: three of the four contributing ClinVar submitters state in
  their own comments that NMD is not expected for this variant, applying
  the same last-exon/near-final-junction rule this evaluator implements.
  This fixture also surfaced a second, latent gap: DMD has no published
  ClinGen VCEP variant-curation specification (checked this batch), so
  PM2/BS1's OBSERVED-branch requirement for a configured per-gene
  threshold — previously invisible, since every prior real DMD fixture
  happened to be genuinely absent from population databases — finally got
  exercised. Rather than inventing a fabricated "real" number or
  misrepresenting the variant's real, low-but-nonzero gnomAD frequency as
  ABSENT, `config/population_thresholds.yaml` now carries an
  explicitly-labeled DMD placeholder reusing CAPN3's own pre-batch-4
  placeholder values (0.0001/0.001) — the same disclosed-placeholder role,
  ready to be replaced the same way batch 4 replaced CAPN3's once a real
  DMD VCEP spec exists. With that placeholder, PM2 lands on the same
  founder-enrichment MANUAL_REVIEW pattern as `CAPN3_c.550del` and
  `CAPN3_c.2120A>G`: all 6 observed gnomAD alleles are African-ancestry,
  a real founder-enrichment signal. Final result: VUS, manual review
  required — two independent, disclosed modeling gaps (not conflicting
  evidence) converging on the same "needs a human" conclusion the real
  multi-lab ClinVar consensus also reached.

Batch 10 didn't add a new fixture either — a seventh unsuccessful BA1
search (see below) led back to `CAPN3_c.2257G>A` (p.Asp753Asn), already
in the set since batch 2 and enriched once already in batch 8, and this
round upgraded it substantially on two fronts:

- **Coordinates and the VCEP classification are now primary-source
  confirmed.** Batch 8 had the real LGMD VCEP's Likely-Benign-via-BS1+BP2
  call from a web search summary only, flagged as unverified. Batch 10
  fetched the primary ClinVar VCV000281081.42 page directly: it confirms
  that call word-for-word (3-star "reviewed by expert panel" status,
  evaluated 2025-10-28) and supplies real coordinates (GRCh38
  chr15:42410660 G>A), so `coordinate_verified` is now `true`.
- **A genuinely opposing real source, not just an unimplemented
  criterion.** A peer-reviewed paper published one month after that VCEP
  evaluation (Bruno et al. 2025, *Int J Mol Sci* 26(23):11384, PMID
  41373542) reports this exact variant as the single most frequent
  finding — 8 of 59 patients (13.5%) — in a multicenter Italian cohort of
  CAPN3-heterozygous patients, with real phenotype data and new
  structural modeling, concluding "the aggregated evidence supports a
  pathogenic role." This doesn't contradict the VCEP's numbers (both cite
  essentially the same gnomAD frequency and REVEL score) — it answers a
  different question. The VCEP classified this variant Likely Benign FOR
  AUTOSOMAL RECESSIVE calpainopathy; the paper argues for pathogenicity in
  a distinct, less common AUTOSOMAL DOMINANT calpainopathy (LGMDD4) model.
  `gene_disease_context.yaml` fixes CAPN3 as one inheritance pattern for
  every variant, so this project structurally cannot represent "benign
  for recessive disease, possibly pathogenic for dominant disease" as one
  variant's answer — a schema-level gap, not a missing criterion,
  previously only mentioned here in passing (batch 8) and now fully
  primary-source-documented on both sides. The engine's own VUS call ends
  up being the most epistemically honest of the three real answers on the
  table (VCEP: Likely Benign; paper: pathogenic-leaning; this project:
  uncertain) specifically because it commits to neither disease model — a
  fortunate byproduct of missing criteria, not a designed strength.

Batch 11 added one more, again while searching (unsuccessfully) for
something else:

- `DMD_c.8944C>T` (p.Arg2982Ter) — a real, cleanly-sourced DMD nonsense
  variant (ClinVar VCV000011211.24, Pathogenic since a 1992 literature
  report and reconfirmed by seven independent submitters since),
  unanimously absent from gnomAD/1000 Genomes/TOPMed. PVS1 Very Strong +
  PM2 Moderate lands on LIKELY_PATHOGENIC — the same shape as
  `DMD_c.2302C>T`, just at a different exon (60 of 79). Found while
  looking for a real BP4-MET fixture (a calibrated computational score
  landing BENIGN from real data), which this project still doesn't have
  — the only real benign-leaning score on record is
  `CAPN3_c.2257G>A`'s INDETERMINATE REVEL (batch 8). That gap remains
  open, genuinely hard to source via search rather than merely
  unexamined, in the same way the BA1 gap is. This fixture doesn't close
  it, but does grow the set (17 of ~20-30) with another fully real,
  unsurprising, easy-to-verify addition.

A sixth round (batch 7) added one more:

- `DMD_c.93+1G>A` — a real, ClinVar-Pathogenic splice donor variant in DMD
  intron 2, the project's second real DMD fixture and the first real DMD
  variant to exercise PVS1's splice-site scope gap (previously only shown
  on CAPN3). Same MANUAL_REVIEW-then-VUS story as the CAPN3 splice
  fixtures, confirming the limitation is a property of the shared
  evaluator, not something specific to one gene's data.

Before landing on that, batch 7 made a further, unsuccessful attempt at
the long-standing BA1-common-variant gap (see below) — tried ClinVar
molecular-consequence searches for synonymous variants, general dbSNP/
gnomAD polymorphism searches, and pulling the full CAPN3 genomic-region
variant list from Ensembl's REST API to check further candidates. None
of these turned up a real CAPN3 variant crossing the 0.3% BA1 threshold
either. This has now been searched for across five separate rounds
(batches 3, 4, 6, and 7) without success.

The real, precisely-sourced calibrated computational score gap mentioned
in every prior round is now closed (`CAPN3_c.2120A>G`'s REVEL=0.966). One
gap remains open, searched for across seven rounds and still not found
rather than silently skipped: a true common CAPN3 variant to exercise
BA1 MET against real data. Searched for repeatedly (batches 3, 4, 6, 7,
8, and 10) via ClinVar molecular-consequence queries, general
population-frequency web searches, an Ensembl REST API pull of the full
CAPN3 genomic-region variant list, and a ClinVar/erepo search
specifically for benign-classified CAPN3 variants. Batch 10 tried the
"fundamentally different approach" this note used to suggest: pulling
every cataloged variant directly from the CAPN3 gene region via
Ensembl's `overlap/region` endpoint (the same technique that found
`DMD_c.11041A>T` in batch 9) instead of ClinVar-searching for one. It
didn't work here — a 53kb gene region returns thousands of dbSNP
entries, the vast majority ultra-rare/singleton, and Ensembl's overlap
response doesn't include population frequency at all (a separate
per-variant lookup is needed for that), so brute-force enumeration
isn't actually more tractable than the ClinVar-anchored searches it was
meant to replace; gnomAD's own GraphQL API — the tool actually built for
"which variants in this gene cross X% frequency" — would likely work,
but is POST-only and out of reach of this project's available fetch
tooling. Batch 10 instead found real value elsewhere: a web search for
"CAPN3 common variant" surfaced a brand-new (Nov 2025) peer-reviewed
paper on `CAPN3_c.2257G>A` (p.Asp753Asn) itself — still, as before, the
closest real candidate to the BA1 threshold at 0.2457% overall AF, just
under the 0.3% cutoff — and that paper turned into a substantial
enrichment of the existing fixture (see below) rather than a new one.
At this point the likeliest explanation for the BA1 gap isn't that the
search strategy is missing something, but that CAPN3 genuinely doesn't
have many easily-indexed common coding variants — consistent with it
being a disease-relevant gene under negative selection even for
variants short of clinical significance. Remains covered only by
hand-built synthetic tests (`tests/unit/test_ba1_bs1_evaluators.py`) as
of this point in the project -- batch 18 (see below) eventually closed
the real-data gap with a DMD variant instead, once it became clear the
gap's own framing ("a true common variant to exercise BA1 MET") was
never actually CAPN3-specific.

Batch 12 didn't add a new variant fixture, but closed a different,
longstanding gap: every Milestone 4 `ClinicalCase` fixture had used
synthetic variants (`DMD_SYNTH_PATHOGENIC_01`, `CAPN3_SYNTH_PATHOGENIC_01`/
`02`), even after the curated set grew to include several real, cleanly
`LIKELY_PATHOGENIC` DMD variants across batches 7-11 — none of them had
ever actually been wired into a case-level fixture. Added
`CASE_DMD_HEMIZYGOUS_MALE_REAL` and `CASE_DMD_FEMALE_CARRIER_REAL`, both
built on the real `DMD_c.2302C>T`, reproducing the existing
male-EXPLAINED / female-MANUAL_REVIEW pair end-to-end on real data for
the first time. The female-carrier case is grounded in real biology, not
just an abstract "what if XX": DMD manifesting carriers are a documented
phenomenon (2.5-7.8% report muscle weakness, up to ~8% present with
dilated cardiomyopathy), and Brioschi et al. (BMC Med Genet 2012;13:73,
DOI 10.1186/1471-2350-13-73) found X-inactivation skewing — the usual
proposed mechanism — did not reliably track with phenotype (skewed in 2
of 6 symptomatic carriers vs. 5 of 11 asymptomatic ones). That real
unpredictability is exactly why this project's MANUAL_REVIEW here is
honest rather than a shortfall: the published literature itself can't
cleanly resolve genotype-to-phenotype for a carrier from X-inactivation
data alone. `CAPN3`'s biallelic-case fixtures remain synthetic for a
structural reason, not an oversight: no real CAPN3 variant in this
project's curated set reaches PATHOGENIC/LIKELY_PATHOGENIC on its own
(every real CAPN3 fixture lands VUS, per the scope gaps documented
throughout this section) — case_interpretation's EXPLAINED branch
requires two independently-QUALIFYING classifications, so a real
biallelic CAPN3 EXPLAINED case isn't buildable from this project's
current real fixtures without also closing one of those variant-level
gaps first.

Batch 13 turned that loose end into a proven claim rather than an
open question. Given (1) the real ClinGen LGMD VCEP threshold this
project adopted for CAPN3 (batch 4) fixes PM2 at SUPPORTING strength,
never MODERATE, and (2) PP3's strength is hardcoded SUPPORTING
everywhere in `pp3.py` with no gene override — the only pathogenic-
direction criteria any single real CAPN3 variant can ever have MET are,
at most, PVS1 (Very Strong, loss-of-function variants only) + PM2
(Supporting), or PM2 (Supporting) + PP3 (Supporting, missense variants
only). Table 5 has no combining rule for "1 Very Strong + 1 Supporting"
alone, nor for "2 Supporting" alone without at least one Moderate or
Strong criterion. So, as this project is currently configured, **no
real CAPN3 variant can reach PATHOGENIC or LIKELY_PATHOGENIC through
this engine, regardless of how strong its individual real-world
evidence is** — not a fixture-search problem, a structural one. This is
now locked in by
`tests/unit/test_engine.py::test_no_real_capn3_variant_currently_reaches_pathogenic_tier`,
which checks every real CAPN3 fixture's golden case against this claim
and is designed to fail loudly (not silently pass) if a future change
(e.g. a gene-specific PM2 Moderate override, or a higher-strength PP3)
ever makes it false. Concretely, this also means `CASE_CAPN3_BIALLELIC_TRANS`-style
real biallelic clinical cases (see "Case-level scope" and batch 12
above) cannot be built from real CAPN3 fixtures until this changes —
closing it would need either implementing more of the real VCEP's
criteria (PM3/PP1/PP4/PS1/PS3), or adopting the VCEP's Bayesian point
system instead of Table 5 (see the "Real ClinGen LGMD VCEP thresholds"
note above) — both already-identified, larger pieces of future work,
not new ones.

Batch 14 added a new criterion, not just a new fixture — see "PM4: a
second new criterion" below for how the gap was found. Its real fixture:

- `CAPN3_c.1401_1403del` (p.Glu467del) — a real, expert-panel-reviewed
  in-frame single-glutamate deletion (ClinVar VCV000553852.19, ClinGen
  LGMD VCEP, Pathogenic since 2025-06-24). Chosen specifically to source
  PM4's first real fixture, and interesting for a reason beyond that: its
  official HGVS (`c.1395GGA[2]`) uses short-tandem-repeat notation — the
  reference is a 3-copy GGA run, contracting to 2 copies — which reads,
  naively, like exactly the "repeat region" PM4 is defined to exclude
  (Richards et al. 2015: "...in a nonrepeat region"). The real VCEP's own
  Pathogenic comment applies PM4 to this variant anyway, with no
  repeat-region caveat. This project's fixture follows that real-world
  precedent (`repeat_region: false`, with the reasoning laid out in full
  in the fixture's own `notes` field) rather than mechanically reading
  the HGVS notation as disqualifying — the distinction being drawn is
  between an incidental 2-3-copy in-gene run (what this variant is) and
  a genuinely repeat-prone locus with a documented high population rate
  of benign in-frame indels (what the PM4 exclusion is actually aimed
  at). PM2 and PM4 both MET (Supporting strength each — this evaluator
  downgrades single-residue in-frame changes by default), landing on
  VUS: a second, independent proof point for batch 13's structural
  finding that no real CAPN3 variant currently reaches Pathogenic/Likely
  Pathogenic through this engine (now covering all three consequence
  shapes: loss-of-function, missense, and in-frame indel/stop-loss).

Batch 15 added one more, in the same round as a wider search for real
BA1 (8th attempt) and BP4-MET (3rd attempt) fixtures hit real tooling
walls (client-rendered ClinVar/ClinVar-Miner search pages, a blocked
gnomAD GraphQL endpoint, a temporary fetch rate limit) rather than
turning up nothing — both gaps remain open, now documented as
tooling-blocked rather than merely unsearched this round:

- `CAPN3_c.598_612del` (p.Phe200_Leu204del) — a second real PM4 fixture,
  chosen specifically to complement `CAPN3_c.1401_1403del`. That one is a
  single-residue in-frame deletion (PM4 downgraded to Supporting by this
  evaluator's single-residue rule); this one deletes five residues, so
  it's the first real fixture to exercise PM4's default Moderate-strength
  branch. Real, expert-panel-reviewed (ClinGen LGMD VCEP, evaluated
  2025-03-18), Pathogenic via PM4 + PM3_Very Strong + PP4_Strong +
  PM2_Supporting. Its `repeat_region: false` curation is more
  straightforward than the first PM4 fixture's: the VCEP's own comment
  states outright that this is "an in-frame deletion of five amino acids
  in a non-repeat region," with no HGVS tandem-repeat notation to reason
  through the way `c.1395GGA[2]` required. PM2 and PM4 both MET (Moderate
  this time, not Supporting), still landing on VUS — a third independent
  confirmation of batch 13's structural finding, now also showing that
  PM4 at full Moderate strength doesn't change the outcome either.

Batch 16 added one more, this project's first DMD missense fixture:

- `DMD_c.10103A>G` (p.Asp3368Gly) — a real, de novo variant found by exome
  sequencing in a woman with X-linked dilated cardiomyopathy (DMD's rare
  cardiac-predominant phenotype), hemizygous in her son, who has elevated
  CK and calf hypertrophy but no cardiac involvement yet — a cardiac-axis
  echo of the same phenotypic-spectrum variability the batch 12 DMD
  female-carrier case already documents on the skeletal-muscle axis.
  Every prior real DMD fixture was truncating or splice; this is the
  project's first DMD missense variant and first DMD fixture with
  computational evidence. Its REVEL score (0.975) is cited directly from
  the source publication (d'Apolito et al. 2024, Int J Mol Sci
  25(5):2787, PMID 38474032) rather than any ClinVar submitter comment.
  DMD has no VCEP-specific REVEL threshold the way CAPN3 does, so this
  fixture uses the real, published, gene-agnostic Pejaver et al. 2022
  calibration instead -- the first real fixture in this project to do so.
  PM2 MET Moderate (ABSENT from gnomAD/ExAC/1000G, DMD's generic-default
  PM2 strength since no VCEP override exists) + PP3 MET Supporting
  (REVEL clears even the Moderate calibration tier) still lands on VUS --
  the same "1 Moderate + 1 Supporting is insufficient" shape already
  established for CAPN3, now shown for DMD too. Unlike most of this
  project's real fixtures, though, this isn't a story about the engine
  under-calling a confidently-classified variant: ClinVar's own aggregate
  is "conflicting classifications" (Invitae: Pathogenic, via de novo/
  segregation/same-residue reasoning this project doesn't implement;
  GeneDx: Likely Pathogenic; Revvity: Uncertain Significance), and the
  discovery paper's own authors, applying ACMG/AMP criteria themselves
  with exactly PM2 and PP3, independently classified it a "rare clinical
  VUS" too -- this engine's VUS call matches the primary literature's own
  classification directly, not just a defensible fallback against a scope
  gap.

Batch 17 didn't add a new variant fixture -- it closed a case-level gap
found the same way PM4's gap was found (batch 14): checking what the
curated set actually exercises, not assuming coverage exists because the
code does. `interpret_x_linked_case`'s final branch -- karyotypic_sex=XY,
but the variant's own classification is neither qualifying (Likely
Pathogenic/Pathogenic) nor benign-side -- returns MANUAL_REVIEW, and had
existed in `clinical.py` since Milestone 4, but no curated `ClinicalCase`
fixture, real or synthetic, had ever actually reached it; every prior
X-linked case used a cleanly (Likely) Pathogenic or (Likely) Benign
variant. `CASE_DMD_HEMIZYGOUS_MALE_VUS_REAL` closes it using
`DMD_c.10103A>G`'s real hemizygous case directly: the source paper's
8-year-old hemizygous son has elevated CK and calf hypertrophy but was
diagnosed only with "unspecified myopathy due to DMD defect" -- his own
clinicians didn't resolve his case to a confident diagnosis either. This
engine's MANUAL_REVIEW answer lands in the same place real clinical
practice did for this exact patient.

Batch 18 closed the longest-standing open gap in this project: a real
common variant to exercise BA1 MET, searched for across nine rounds since
Milestone 1 without success (see above) -- the search had stayed
CAPN3-scoped for that entire time even though the gap's own stated
framing never required that gene.

- `DMD_c.5234G>A` (p.Arg1745His) — a real, extremely common DMD missense
  polymorphism (ClinVar VCV000094657.45, Benign, 16 of 19 contributing
  submissions). gnomAD exomes overall AF 48.59%, corroborated by TOPMed
  (38.1%), 1000 Genomes (46.5%), and ExAC (51.4%) — a near-50/50 common
  variant in most populations, lowest in the African/African American
  gnomAD subgroup at ~6.7%, still comfortably above the 5% generic BA1
  threshold on its own. Found via a ClinVar Miner per-submitter
  benign-variant listing page, which rendered as static HTML where the
  gene-level search/listing pages tried in earlier rounds (and again in
  batches 14/15) needed JavaScript — the same tooling obstacle, worked
  around by finding one page shape that doesn't require it. BA1 MET
  (Stand-Alone) resolves this directly to BENIGN per Table 5 — this
  project's first real fixture to land there; every real fixture before
  it landed somewhere in the VUS/Likely-Pathogenic/Likely-Benign middle
  band, a side effect of which scope gaps the real fixtures happened to
  probe rather than anything about BENIGN being harder to reach.

Batch 19 closed that last open real-data gap:

- `DMD_c.5163G>C` (p.Lys1721Asn) — a real, ClinVar-documented DMD
  missense variant (Variation ID 455905, rs72468630), classified Benign
  or Likely Benign across six independent RCV records spanning multiple
  conditions and submitters with no conflicts. Found via
  `myvariant.info`'s public REST API (`dbnsfp.revel` field), a new
  discovery this round — the first working source in this project for
  precomputed REVEL scores after both NCBI eutils and gnomAD's GraphQL
  API were confirmed blocked/non-functional in this environment. REVEL
  = 0.167, clearing even DMD's real Pejaver et al. 2022 Moderate-tier
  BP4 threshold (<=0.183), not just the Supporting one (<=0.290) —
  unambiguously benign-zone. Population frequency (0.021%) was
  deliberately chosen to sit above DMD's PM2 threshold but below its
  BS1 placeholder, so BP4 is the sole benign-direction criterion MET,
  not overshadowed by frequency evidence the way BA1 dominated
  `DMD_c.5234G>A` (batch 18). BP4 MET (Supporting) alone doesn't
  satisfy Table 5's Likely-Benign rule (needs 1 Strong + 1 Supporting,
  or 2 Supporting), so this lands on VUS — the direct real-data
  complement to `CAPN3_SYNTH_LIKELY_BENIGN_01`, which was built with
  both BS1 and BP4 MET specifically because neither alone suffices.

Every real-data gap the project set out to find — BA1 MET, BS1 MET
(pre-existing), PM4 at both strengths, and now BP4 MET — has a real
fixture behind it. The curated variant set has reached the low end of
its original ~20-30 target with room to keep growing, but no longer has
a specific "still missing" criterion-level story driving what gets added
next.

Reaching the full ~20-30 variant set is expected to take several more
rounds of this same process (research a real variant, ground its
evidence, hand-derive the expected result, verify against the engine).

Batch 22 added one more real fixture, found organically while
researching real PS1/PM5 precedent pairs rather than through a dedicated
search:

- `CAPN3_c.1342C>T` (p.Arg448Cys) — a real CAPN3 missense variant
  (ClinVar VCV000280038.68, RCV006646441.1), classified Pathogenic by the
  ClinGen LGMD VCEP itself (reviewed by expert panel, evaluated
  2026-01-20), via PM3_Strong, PP4_Strong, PM2_Supporting, PP3, and PM5 —
  the PM5 citation points directly at this project's own pre-existing
  `CAPN3_c.1343G>A` (p.Arg448His) fixture as its precedent. This project's
  engine implements PM2, PP3, and (as of this batch) PM5 of those five —
  PM3_Strong and PP4_Strong remain case-level/segregation criteria this
  project's variant-only design doesn't reach, and they are doing most of
  the real work in the VCEP's own Pathogenic call. PM2 Supporting + PP3
  Supporting + PM5 Moderate doesn't clear any Table 5 pathogenic-tier
  rule — VUS here, a sixth independent confirmation of batch 13's
  structural finding (no real CAPN3 variant reaches Pathogenic/Likely
  Pathogenic through this engine), now shown to persist even with PM5
  newly implemented and MET on real data.

  Finding this pair also surfaced that `CAPN3_c.1343G>A` itself needed
  updating: when it was added (batch 5), ClinVar's aggregate call was
  Uncertain Significance and no computational evidence had been gathered.
  The same 2026-01-20 LGMD VCEP evaluation that classified
  `CAPN3_c.1342C>T` also reclassified `CAPN3_c.1343G>A` Pathogenic (via
  PM3_Strong, PP4_Strong, PM2_Supporting, PP3 — notably *not* citing PM5
  in this direction, even though `CAPN3_c.1342C>T` is a real, established
  different-amino-acid-change precedent at the same residue). This
  project's fixture for `CAPN3_c.1343G>A` was updated with the real REVEL
  score (0.904, now a real `computational_evidence` entry closing the gap
  the original notes flagged as missing) and verified coordinates, but
  deliberately does *not* also add reciprocal PM5 evidence citing
  `CAPN3_c.1342C>T` — mirroring the real VCEP's own choice not to invoke
  it, rather than second-guessing an unexplained-but-plausibly-principled
  omission (both variants were curated and reached Pathogenic from the
  same joint evaluation event without needing PM5, which the VCEP's own
  PM5 specification text says is fine: "PM5 can potentially be applied to
  multiple amino acid changes at the same residue as long as the variant
  classification that determines the strength level does not depend on
  PM5 application" — applying it mutually from a single joint curation
  risks the same kind of circularity `clinical.py`'s own docstring already
  flags for PM3). An honest "not curated one way or the other," not a
  guess, and disclosed as an open discussion point rather than resolved
  by fiat. See "PS1 and PM5: same-residue precedent evidence" above for
  the full criterion-level design story, and `CAPN3_SYNTH_PS1_01`'s notes
  in `data/curated/variant_evidence.json` for the disclosed synthetic
  fixture built to cover PS1's MET path, since no real "same amino acid
  change via a different nucleotide" precedent pair was found this round.

**Milestone 5: Bayesian point-based combining (batch 20).** Every prior
milestone changed what evidence gets gathered or what criteria get
evaluated; Milestone 5 changes neither -- it adds a second way to combine
the exact same evaluated criteria into a classification.
`src/variant_classifier/bayesian.py` implements Tavtigian et al. 2020's
"naturally scaled point system" (Human Mutation 41(6):1023-1041), offered
*alongside* `engine.py`'s classic Table 5 rules rather than replacing
them -- `combine_bayesian()`/`classify_bayesian()` sit next to
`combine()`/`classify()`, both callable, both tested, both real. Point
values and thresholds are quoted directly from the ACGS 2024 UK Practice
Guidelines for Variant Classification (which cites Tavtigian et al. 2020
verbatim; see `bayesian.py`'s module docstring for the exact quote):
pathogenic-direction points (Very Strong=8, Strong=4, Moderate=2,
Supporting=1) and benign-direction points (Strong=-4, Moderate=-2,
Supporting=-1) sum to a single net total, classified >=10 Pathogenic, 6-9
Likely Pathogenic, 0-5 VUS, -1 to -5 Likely Benign, <=-6 Benign -- with
two carried-over exceptions, not new inventions: BA1 Stand-Alone still
bypasses point-summing entirely and resolves straight to BENIGN, exactly
as it does in Table 5; and, per Tavtigian et al. 2020 (via ACGS 2024), a
minimum of two contributing criteria is required to reach any
Likely/definite Pathogenic or Likely/definite Benign result, so a single
criterion's points alone -- however large -- caps out at VUS.

Why this was worth building now rather than left as a roadmap bullet:
batch 4 already produced a concrete, real discrepancy (`CAPN3_c.1939G>T`
dropping from LIKELY_PATHOGENIC to VUS once CAPN3's real VCEP threshold
set PM2 to Supporting strength) and explicitly named the Bayesian point
system as the real-world fix the VCEP itself uses. Rather than continue
describing that gap in prose, Milestone 5 made it a running, tested
comparison. Hand-deriving Bayesian point totals for all 22 curated
fixtures (independently, before running any new code, same convention as
every other batch) turned up four real divergences, not just the one
already anticipated:

- `CAPN3_c.1939G>T` (PVS1 Very Strong + PM2 Supporting = 9 points) — VUS
  under Table 5 (no rule for "1 Very Strong + 1 Supporting" alone),
  LIKELY_PATHOGENIC under Bayesian. The anticipated case, now concrete.
- `DMD_SYNTH_PATHOGENIC_01`, `DMD_c.2302C>T`, `DMD_c.8944C>T` (PVS1 Very
  Strong + PM2 Moderate = 10 points each) — LIKELY_PATHOGENIC under Table
  5 (its flat PATHOGENIC tier needs "1 Very Strong + >=2 Moderate", not
  just one), PATHOGENIC under Bayesian. A second, independent divergence
  shape found only by checking every fixture rather than just the one
  already suspected — two of these three are real fixtures.

UPDATED in batch 25: adding real PS3 evidence to `CAPN3_c.550del` (see
"PS3 and BS3: functional-evidence criteria (batch 25)" above) turned it
into a fifth real divergence, the same shape as `CAPN3_c.1939G>T`'s —
PVS1 Very Strong + PS3 Supporting = 9 points, VUS under Table 5, LIKELY_PATHOGENIC
under Bayesian. Not a new divergence *shape*, but a second independent
real-data example of the first one, found as a direct side effect of
adding unrelated evidence rather than a targeted search for more
divergences.

UPDATED in batch 28: adding real PM3 evidence to that same
`CAPN3_c.550del` (see "PM3: implemented (batch 28)" above) removed this
divergence rather than adding one. PM3 Strong (a Strong-strength
criterion, this project's first) plus the existing PVS1 Very Strong now
satisfies Table 5's own "1 Very Strong + >=1 Strong" pathogenic rule
directly, so both systems reach the same result -- PATHOGENIC, 13
Bayesian points. Back down to four documented divergences; see
`test_bayesian_diverges_from_table5_for_exactly_the_four_documented_fixtures`.

Both shapes share the same root cause: Table 5's combining rules were
built as a discrete, hand-enumerated list (Table 5 in Richards et al.
2015), not derived from the point values later fit to approximate it, so
a few real evidence combinations that clear the Bayesian point thresholds
were simply never enumerated as their own Table 5 rule. This project's
own further literature check (Tavtigian et al. 2020 itself) surfaces a
third, opposite-direction inconsistency neither this project's fixtures
happen to exercise: Table 5's "2 Strong" Pathogenic rule is reported as
the weakest of its eight Pathogenic combining paths (posterior probability
0.975 vs >0.99 for the rest), and 2 Strong criteria (4+4=8 points) fall
*short* of the 10-point Bayesian Pathogenic threshold — a case where Table
5 is arguably more permissive than the point system, not less. Documented
in `bayesian.py`'s module docstring for completeness; this project has no
Strong-strength pathogenic-direction evaluator, so no current fixture can
actually reach it.

One more real design difference, not a bug: Table 5's `combine()` can
report `conflicting_evidence_flag=True` when pathogenic-direction and
benign-direction criteria satisfy separate rules simultaneously.
`combine_bayesian()` cannot produce that state by construction — a single
net point sum always resolves to exactly one band. This is one of the
Bayesian framework's own claimed advantages, not something this project
invented; `combine_bayesian()`'s rationale still states both directions'
subtotals whenever both contributed, so nothing is hidden, just resolved
differently.

Case-level ripple effects were checked, not assumed. `clinical.py` never
calls `engine.classify()` directly — every function takes a pre-computed
`classifications` dict as a parameter, a separation of concerns baked in
since Milestone 4 specifically so case-level reasoning wouldn't need to
know or care how a variant's classification was produced. Feeding every
curated `ClinicalCase` fixture Bayesian-derived classifications instead of
Table 5 ones (`test_clinical_case_interpretation_agnostic_to_combining_system`)
reproduces every existing case-level golden expectation exactly, including
for `CASE_DMD_HEMIZYGOUS_MALE_REAL` and `CASE_DMD_FEMALE_CARRIER_REAL`,
built on `DMD_c.2302C>T` — even though that variant's own tier moves from
LIKELY_PATHOGENIC to PATHOGENIC between the two systems, `clinical.py`'s
`_QUALIFYING` set treats both tiers identically for case-level purposes,
so nothing downstream changes. No case-level golden case needed touching.

Golden cases for this milestone live in a new, parallel file,
`validation/golden_cases/variant_golden_cases_bayesian.yaml`, rather than
added as extra fields inside the existing Table 5 file — the per-criterion
evidence and evaluator results are identical between the two systems (only
the combining step differs), so this file only records what's new: each
fixture's Bayesian point total and resulting class, with a curator_note
explaining any divergence from the Table 5 result. `test_bayesian.py` is
the dedicated test suite (13 tests): the headline
`test_bayesian_matches_hand_derivation_for_all_curated_bundles` check, a
`test_bayesian_diverges_from_table5_for_exactly_the_five_documented_fixtures`
regression lock (so a future evaluator or threshold change that silently
creates or removes a divergence gets noticed), point-value spot checks
against the Tavtigian 2020 values quoted above, and hand-built edge cases
for the 2-criterion minimum rule (including the two exact examples ACGS's
2024 guidelines themselves use: `PVS1_vstr` alone and `BP4_sup` alone —
`DMD_c.5163G>C` is a real fixture for the `BP4_sup`-alone case, and
`PVS1_vstr` alone was illustrated by the real fixture `CAPN3_c.550del`
through batch 24; as of batch 25 that fixture has real PS3 evidence too
(see above), so `PVS1_vstr` alone is now covered only by the hand-built
`test_single_very_strong_criterion_alone_is_vus_not_likely_pathogenic`,
not a real fixture — disclosed here rather than left silently inaccurate).

What Milestone 5 deliberately does not do: it does not change
`engine.classify()`'s behavior, does not change which combining system
`clinical.py` is fed by default (that remains a caller's choice, still
Table 5 in every existing test), and does not re-derive or replace any
existing Table 5 golden case. Adding a second, real, working combining
system was the goal — deciding to make it the project's default, if that
ever happens, is a separate decision for a separate day, not bundled into
this milestone.

**Case-level scope.** clinical.py deliberately does not extend the
per-variant evaluator pattern from Milestone 2/3 for case-level questions
— reasoning about "does what was found in this patient explain their
disease" after each involved variant already has its own classification
sidesteps circularity rather than working around it (PM3 itself is
handled differently — see "PM3: implemented (batch 28)" above — but the
same instinct against re-deriving classification-dependent facts live
applies here too). Scope, stated plainly: autosomal recessive handles
exactly one or two variants, with phase required whenever there are two.
X-linked, karyotypic_sex=XY, handles exactly one variant (hemizygous male
— the clean, confidently-resolved case since Milestone 4). X-linked,
karyotypic_sex=XX, handles one OR two variants as of batch 29 — see "X-
linked female/other-karyotype case interpretation (batch 29)" below for
the full design, including the real, documented, X-inactivation-
independent biallelic mechanism this newly resolves to EXPLAINED, and why
a single heterozygous variant (or a two-variant case that doesn't clear
that specific bar) still correctly stays MANUAL_REVIEW. Any other
karyotypic sex (OTHER, UNKNOWN) is still deferred to MANUAL_REVIEW with a
single variant only — OTHER lumps together karyotypes with genuinely
different X-linked dosage biology (X0/Turner functionally hemizygous;
XXY diploid-X; mosaicism neither uniformly) this project does not
attempt to disambiguate.

**PVS1 scope.** The full PVS1 decision tree (Abou Tayoun et al. 2018)
branches on protein-domain criticality and constitutive-exon-splicing
information this project doesn't fully model. This evaluator returns MET
for the cases it can defend end-to-end: an early frameshift or nonsense
variant predicted to trigger nonsense-mediated decay in a gene with an
established loss-of-function mechanism; (added batch 26 — see "Splice-RNA
evidence feeds PVS1 directly (batch 26)" below) a splice donor/acceptor
variant with a real RNA/splicing assay confirming a null-equivalent
transcript; or (added batch 27 — see "PVS1 start-loss: the
alternative-start-codon rule (batch 27)" below) a start-loss variant with
either no downstream in-frame alternative start codon at all, or one that
clears the real automatic-downgrade bar. Everything harder — truncations
that escape NMD (typically last-exon), splice donor/acceptor variants
without threshold-free RNA evidence, and start-loss variants whose
alternative start codon doesn't clear that bar (or hasn't been assessed)
— returns MANUAL_REVIEW with a rationale explaining why, rather than a
guessed MET or NOT_MET. Non-null-variant consequence types (missense,
synonymous, etc.) return NOT_APPLICABLE. `TranscriptConsequence` requires
an explicit `nmd_predicted` value for both frameshift and stop-gained
variants for exactly this reason — this requirement was originally
frameshift-only in Milestone 1 and widened here once the evaluator needed
it for stop-gained variants too. All three documented scope gaps (splice
donor/acceptor, start-loss, NMD-escaping truncations) now have at least
one real fixture exercising them, as of batch 9's `DMD_c.11041A>T` — see
"Expanding the curated set" below, and two of the three (splice-site,
start-loss) have been substantially narrowed by batches 26/27
respectively. The exact protein-domain-criticality thresholds within the
NMD-escape branch, and the remaining percentage-based thresholds Walker
et al. 2023's full splicing decision tree defines, remain open — see
"Splice-RNA evidence feeds PVS1 directly (batch 26)" and "PVS1 start-loss:
the alternative-start-codon rule (batch 27)" below for exactly what is and
isn't covered by each extension.

**PM4: a second new criterion (batch 14).** This project's first six
evaluators were all built in Milestones 1-3; PM4 is the first one added
since. It was found by accident of investigation, not by plan: batch 14
went looking for a real stop-loss variant to fill in what looked like a
PVS1 scope gap (`PVS1` returns `NOT_APPLICABLE` for `stop_lost`). Before
writing a fix, that assumption was checked against the primary criterion
definitions (Abou Tayoun et al. 2018's PVS1 decision tree; Richards et
al. 2015 Table 3) — and turned out to be wrong. PVS1 is specifically
about null/loss-of-function variants; a stop-loss variant produces an
elongated protein, not a null one, so `NOT_APPLICABLE` was already
correct. Richards et al. 2015 Table 3 assigns stop-loss variants (and
in-frame indels) to **PM4** — "protein length changes as a result of
in-frame deletions/insertions in a nonrepeat region or stop-loss
variants" (Moderate) — a criterion this project had never implemented,
even though `Consequence` has always included `INFRAME_DELETION`,
`INFRAME_INSERTION`, and `STOP_LOST`. No curated fixture had ever used
any of them, so the gap stayed invisible until batch 14 went looking.
Confirmed via a full-suite check before writing any evaluator code: zero
existing fixtures use a PM4-relevant consequence, so implementing it
carries zero risk of silently changing any existing classification (all
17 pre-existing golden cases gained only an additive `PM4: NOT_APPLICABLE`
line). `src/variant_classifier/evaluators/pm4.py`: in scope only for
`INFRAME_DELETION`/`INFRAME_INSERTION`/`STOP_LOST` (everything else is
`NOT_APPLICABLE`); not gated on an established loss-of-function
mechanism the way PVS1 is (Richards et al. 2015 doesn't condition PM4 on
mechanism); `NOT_MET` if `TranscriptConsequence.repeat_region` is `true`
(required to be stated explicitly for PM4-relevant consequences, same
never-guess convention as `nmd_predicted`); otherwise `MET`, at
Supporting strength if `protein_length_change_aa == 1` (the ClinGen SVI's
general caution against full Moderate strength for single-residue
indels absent gene-specific evidence) or Moderate otherwise, including
when the size isn't recorded at all (a disclosed simplification, not a
silent guess — the rationale always states which case applied). See
"Expanding the curated set" above for `CAPN3_c.1401_1403del`, its real
fixture.

**PS1 and PM5: same-residue precedent evidence (batch 22).** The
Roadmap's "Later" bullet had named PS1 and PM5 as the most tractable next
criteria since Milestone 3 ("same-residue pathogenic/benign lookups, no
new evidence type needed") -- batch 22 followed through on that. Real
definitions confirmed against Richards et al. 2015 Table 3 before writing
any code: **PS1** ("same amino acid change as a previously established
pathogenic variant, regardless of nucleotide change," Strong) and **PM5**
("novel missense change at an amino acid residue where a different
missense change determined to be pathogenic has been seen before,"
Moderate). Both carry a caveat this project takes seriously rather than
glossing over: the ClinGen SVI Splicing Subgroup (Walker et al. 2023,
*Am J Hum Genet*, PMID 37352859) explicitly lists PS1 and PM5 among the
codes needing care when the nucleotide change under evaluation might
itself be acting through altered splicing rather than through the amino
acid substitution -- in which case the whole "same/different amino acid
change" comparison isn't valid.

Design: a new model, `SameResidueEvidence` (`src/variant_classifier/models/same_residue_evidence.py`),
added to `VariantEvidenceBundle` as one more optional field alongside
`computational_evidence` -- not a new evidence-domain the way
`ComputationalEvidence` was for PP3/BP4, confirming the Roadmap's
original "no new evidence type needed" framing. Each precedent
(`ps1_precedent_established`/`pm5_precedent_established`, plus a required
classification and citation whenever `True`) is a curated fact about a
**different**, previously classified variant, sourced externally
(ClinVar/a VCEP curation) -- never this project's own engine output, and
never computed by scanning this project's other curated fixtures. That
distinction matters: `clinical.py`'s own docstring explains PM3 ("detected
in trans with a pathogenic variant") can't be a per-variant evaluator
because variant A's PM3 would depend on variant B's classification from
this same engine, a structural circularity. PS1/PM5 reference an external
authority's already-established classification instead, so no such
circularity exists here. `splice_impact_excluded` must be stated
explicitly (`True`/`False`) whenever a precedent is recorded, never left
unstated -- the same "never silently guess" convention as `nmd_predicted`
(PVS1) and `repeat_region` (PM4); left unstated or `False`, both
evaluators return `MANUAL_REVIEW` rather than guessing which way the
caveat cuts.

Scope, disclosed rather than assumed: only `MISSENSE_VARIANT` is in scope
(everything else is `NOT_APPLICABLE`), and precedent strength downgrades
one level when the precedent itself is only Likely Pathogenic rather than
Pathogenic (PS1 Strong->Moderate, PM5 Moderate->Supporting) -- a
convention at least one real VCEP (RYR1, Malignant Hyperthermia
Susceptibility) documents explicitly, though this project has not
verified it as a universal ClinGen SVI mandate. While researching a real
fixture, the actual ClinGen LGMD VCEP specification for CAPN3 (v2.0,
cspec.genome.network/cspec/ui/svi/doc/GN187) turned out to already define
PS1/PM5 for CAPN3 in far more depth than implemented here: a minimum
REVEL score and excluded SpliceAI score for the variant under curation,
no benign missense variation permitted at the residue, exclusion of
missense changes encoded by the first/last 3 nucleotides of an exon (a
splice-region proxy), and counting *multiple* precedent variants toward
Strength (2 Pathogenic or 3 Likely Pathogenic = Strong PS1; similarly for
PM5). None of that gene-specific machinery is implemented here -- this
evaluator applies the generic Richards et al. 2015 definition and the
base splice-vs-protein-level caveat only, the same "deliberately partial"
treatment PVS1 has had since Milestone 2, disclosed in `ps1.py`/`pm5.py`'s
own module docstrings rather than silently narrower than it looks.

A real fixture pair fell directly out of this research, not a separate
search: `CAPN3_c.1342C>T` (p.Arg448Cys) and the pre-existing
`CAPN3_c.1343G>A` (p.Arg448His) turned out to be the exact same-residue
pair the real ClinGen LGMD VCEP itself classified Pathogenic on the same
day (2026-01-20), citing each other reciprocally as PM5 evidence in one
direction. Re-checking `CAPN3_c.1343G>A` against the primary ClinVar
record (a routine step, same discipline as batches 8/10) found its own
classification had changed since this fixture was first added -- it was
Uncertain Significance at the time, now real-world Pathogenic -- so this
batch updated that fixture's `computational_evidence` (a REVEL score of
0.904 the earlier notes had flagged as missing) and notes accordingly,
rather than leaving stale information next to the new PM5 fixture. See
"Expanding the curated set" below for the full story, including why this
project deliberately did *not* also curate reciprocal PM5 evidence for
`CAPN3_c.1343G>A` itself (the real VCEP didn't either, for a documented
reason worth preserving rather than second-guessing). No real "same amino
acid change via a different nucleotide" precedent pair (a PS1 example) was
found for CAPN3 or DMD this round -- disclosed as an open real-data gap,
the same "searched, not found, still open" treatment other gaps have had
in this project (e.g. BA1 across nine rounds before batch 18 closed it).
PS1's MET/Strong path is instead exercised by a disclosed synthetic
fixture, `CAPN3_SYNTH_PS1_01`, modeled directly on Richards et al. 2015's
own PS1 example ("Val->Leu caused by either G>C or G>T in the same
codon").

Tests: `tests/unit/test_ps1_pm5_evaluators.py` (18 tests) -- golden-case
cross-check plus hand-built edge cases for every branch (`NOT_APPLICABLE`,
`NOT_EVALUATED`, `NOT_MET`, `MANUAL_REVIEW`, `MET` at both strength tiers
for both criteria) and `SameResidueEvidence`'s own schema validation.
`bayesian.py` needed no changes at all to pick up PS1/PM5 -- it calls
`engine.evaluate_all()` directly, so the two new evaluators are already
included in its point totals; `test_bayesian_diverges_from_table5_for_exactly_the_four_documented_fixtures`
(batch 20's regression lock) still passes unchanged, confirming PS1/PM5
didn't quietly create or remove a Table-5-vs-Bayesian divergence.

**DMD CNV/structural-variant scoring: a deliberately partial
implementation (batch 23).** Batch 22 sized this gap (`gene_disease_context.yaml`
has disclosed since Milestone 4 that "DMD pathogenic variants are very
often multi-exon deletions/duplications ... which this project's
`TranscriptConsequence` model does not represent") and concluded it needed
a new identity representation and a new, parallel scoring system -- not an
incremental evaluator. Batch 23 built exactly that, scoped to DMD
deletions only (confirmed with the user before writing any code):

- **`models/cnv_deletion_evidence.py` (`CnvDeletionEvidence`)** -- a new,
  separate evidence shape (not an extension of `VariantEvidenceBundle`):
  a genomic interval/exon-range identity plus which parts of the gene are
  affected (whole gene; 5' end + CDS; 3' end +/- other exons; purely
  intragenic) and, for intragenic deletions, an explicit
  `reading_frame_effect` (`CnvReadingFrameEffect`: `OUT_OF_FRAME` /
  `IN_FRAME` / `UNKNOWN`) and `nmd_predicted` -- this project's direct
  link to the Aartsma-Rus DMD reading-frame rule (Aartsma-Rus et al. 2006,
  PMID 16770791). Both fields follow the same "never silently guess"
  convention as `TranscriptConsequence.nmd_predicted` and
  `SameResidueEvidence.splice_impact_excluded`: required exactly when the
  scoring decision actually depends on them, never left unstated.
- **`models/cnv_category_result.py` / `cnv_provisional_classification.py`**
  -- parallel to `CriterionResult`/`ProvisionalClassification`, deliberately
  NOT a reuse: a CNV category code (`2A`, `2C`, `2D`, `2E`, `2F`) is not a
  member of the ACMG/AMP 28-code vocabulary `CriterionResult.code`
  enforces, and carries a raw point value rather than a
  `CriterionStrength` tier. The five-tier output vocabulary
  (`ProvisionalClass`) and `ClassificationStatus` ARE reused as-is -- both
  frameworks genuinely produce the same five-tier scale.
- **`config/dosage_sensitivity.yaml` + `loader.load_dosage_sensitivity()`**
  -- a small per-gene config (mirroring `population_thresholds.yaml`'s
  pattern) recording ClinGen's own Haploinsufficiency curation. Only DMD
  is populated (HI score 3, "sufficient evidence") -- CAPN3 is deliberately
  absent: it's autosomal recessive, and ClinGen's haploinsufficiency
  dosage-sensitivity framework is scoped to single-copy-loss-causes-disease
  mechanisms (dominant, or X-linked hemizygous as in DMD), not to a
  recessive gene where one copy's loss alone doesn't cause disease -- a
  real biological reason, not an oversight.
- **`cnv_scoring.py` (`score_cnv_deletion()`)** -- a new combining module,
  parallel to `bayesian.py`/`engine.py` but never mixed with either. Point
  values and cutoffs are quoted from ClassifyCNV (Gurbich & Ilinsky 2020,
  *Sci Rep* 10:20375), an open-source, peer-reviewed reimplementation of
  the Riggs et al. 2020 rubric -- fetched directly from
  github.com/Genotek/ClassifyCNV during this batch's research, since the
  primary paper (reCAPTCHA-blocked) and the official ClinGen CNV
  calculator (cnvcalc.clinicalgenome.org, a JS app that timed out on
  fetch) were both unreachable. A disclosed reliance on a secondary but
  primary-adjacent, executable source, not an invented rubric.

What is actually implemented -- only Section 2 (dosage-sensitivity /
haploinsufficiency-overlap), loss/deletion side only:

  - **2A** (1.0 pts): the whole gene is deleted and it's an established
    (HI=3) dosage-sensitive gene.
  - **2C** (0.9 pts): the deletion removes the gene's 5' end (5'UTR/first
    exon) and coding sequence.
  - **2D** (0.9 pts, or 0.3 if confined to the last exon's CDS with no
    other exons involved): the deletion removes the gene's 3' end.
  - **2E** (0.9 pts): an intragenic (both ends intact) deletion that is
    out-of-frame and predicted to trigger nonsense-mediated decay -- the
    direct Aartsma-Rus reading-frame-rule link.
  - **2F** (-1.0 pts): the deletion falls completely inside an established
    ClinGen benign copy-number region.
  - **`NONE_APPLICABLE`** (0 pts): this project's OWN bookkeeping label
    (not a Riggs/ClassifyCNV code) for an intragenic deletion that is
    in-frame or has unknown frame effect -- none of 2A/2C/2D/2E/2F apply.
    The real Riggs rubric likely has a specific code for this shape
    (candidates seen in secondary sources: 2B, 2G) but this project has
    not independently verified either one's definition or point value, so
    it reports zero points under a disclosed label rather than guessing.

Explicitly deferred, named rather than hidden: Section 1 (genomic content
-- moot, every fixture here is by construction a DMD-overlapping
deletion), Section 2H (predicted-but-not-established HI via
DECIPHER/pLI/LOEUF -- needs three external prediction datasets this
project doesn't have), duplications entirely (the gain-side rubric),
and Sections 3 (gene count -- moot for a single-gene CNV), 4 (case/
case-control/population evidence), and 5 (inheritance/family history).

**Curated CNV fixtures** (`data/curated/cnv_deletion_evidence.json`, 3
cases, validated separately from the point-mutation curated set) and
**golden cases** (`validation/golden_cases/cnv_deletion_golden_cases.yaml`,
hand-derived before running `cnv_scoring.py`, same discipline as every
other golden-case file in this project):

  - `DMD_CNV_del_ex47_50` -- REAL: ClinVar RCV000813350.9 / VCV000656845
    (`NC_000023.11:g.(?_31819965)_(31929755_?)del`), Labcorp Genetics
    (formerly Invitae), germline classification Pathogenic (1 star,
    evaluated 2018-11-08). Submitter-confirmed out-of-frame deletion of
    DMD exons 47-50. Scores category 2E (0.9 pts) -> **LIKELY_PATHOGENIC**
    under this project's Section-2-only rubric -- intentionally less
    severe than ClinVar's real Pathogenic call, since that call also
    draws on case/family evidence (Riggs 2020 Sections 4/5) this milestone
    doesn't implement. A disclosed gap, not an error, the same treatment
    CAPN3's PVS1+PM2 Table-5-vs-Bayesian discrepancy already has.
  - `DMD_CNV_del_whole_gene_Xp21` -- literature-grounded (not a single
    pinned ClinVar accession): the real Xp21 contiguous gene deletion
    syndrome / complex glycerol kinase deficiency, in which a deletion
    spanning the entire DMD gene causes dystrophinopathy (over 100 male
    patients reported; e.g. PMID 23739620, PMC8543963). Scores category 2A
    (1.0 pts) -> **PATHOGENIC**.
  - `DMD_CNV_del_ex45_47_inframe` -- literature-grounded: the classic
    in-frame DMD exon 45-47 deletion, the single most common in-frame
    deletion causing Becker muscular dystrophy (~25-30% of BMD in-frame
    deletions across multiple published cohorts). Scores
    `NONE_APPLICABLE` (0 pts) -> **VUS** under this project's Section-2-only
    rubric -- a deliberate, disclosed gap: real Becker in-frame deletions
    of this kind are typically classified Pathogenic/Likely Pathogenic in
    clinical practice via segregation/case-series evidence this milestone
    doesn't implement. This golden case exists specifically to make that
    gap visible and testable.

**Tests:** `tests/unit/test_cnv_scoring.py` (20 tests) -- golden-case
cross-check against all 3 curated CNV fixtures, plus hand-built edge
cases for every decision branch the curated set doesn't happen to
exercise (2C, both 2D point values, 2F, the out-of-frame-but-NMD-escaped
case, the unestablished-gene whole-gene-deletion case) and every
`CnvDeletionEvidence`/`CnvCategoryResult` schema-validation rule.

This remains a genuinely partial implementation by design, following the
same "deliberately partial first version, gaps named rather than hidden"
convention as PVS1 (Milestone 2), PM4 (batch 14), and PS1/PM5 (batch 22).
A full implementation would still need: duplications; Sections 1, 2H,
3, 4, and 5 of the Riggs rubric; and either independent verification of
the primary paper's exact 2B/2G definitions or a documented decision to
keep relying on ClassifyCNV's reimplementation.

**DMD CNV/structural-variant duplication scoring: a deliberately
narrower slice than deletions, with disclosed point-value uncertainty
(batch 24).** Extending batch 23's deletion-only CNV scoring to
duplications looked, at first, like a same-shape extension: mirror
`CnvDeletionEvidence` with a `CnvDuplicationEvidence`, mirror the category
decision tree. Research before writing any code (per this project's
standing discipline) found two real complications that changed the scope:

1. **Whole-gene triplosensitivity (TS) scoring doesn't apply to DMD.**
   ClinGen's own DMD-specific Dosage Sensitivity Curation states plainly:
   "whole gene duplications have not been reported in association with
   clinical phenotypes" for DMD. The real, clinically relevant DMD
   duplication mechanism is not a "triple dose" TS effect at all --
   ClinGen's curation continues: "intragenic DMD duplications and
   triplications have been reported in patients with DMD and BMD,
   presumably by a loss-of-function-type mechanism." So this batch does
   NOT implement whole-gene TS scoring (no config, no category) --
   `whole_gene_duplicated` is representable in the model for completeness,
   but always scores `NONE_APPLICABLE`, an explicit, disclosed gap rather
   than an invented TS category with no real DMD data to ground it.
2. **The real, relevant mechanism -- a tandem duplication with a
   breakpoint inside the gene, disrupting it via the same Aartsma-Rus
   reading-frame rule already used for deletions -- has a confirmed real
   category in the Riggs et al. 2020 rubric, but this project could only
   partially verify its numbers.** An inter-laboratory CNV-classification
   concordance study (PMC8960312) discusses real disagreement over "the
   use of 2K (0.45 points) or 2J (0 point) when a copy number gain
   breakpoint was observed for the established HI genes" -- confirming
   the category exists, but not which condition (out-of-frame vs
   in-frame/unknown) maps to which value, nor the real letter code (very
   likely NOT "2A" the way ClassifyCNV's own gain-side code might
   suggest -- see below). A breakpoint study of 119 gain CNVs also
   confirmed a hard prerequisite this project adopted directly: 83% were
   tandem and direct, "with the majority of the remainder interpreted as
   VUS because the effect could not be determined" -- so `is_tandem` must
   be confirmed before any frame-effect call is made at all.

Presented to the user as an explicit choice before writing any scoring
code (see the batch 24 scoping conversation): implement the real DMD
mechanism anyway, inferring the higher point value (0.45) for the
out-of-frame/disruptive case and the lower value (0) for the in-frame/
uncertain case by direct analogy to the loss side's own pathogenic-vs-
uncertain split (2E vs `NONE_APPLICABLE`) -- disclosed as an inference,
not a verified fact, in the model docstring, the scoring rationale text,
and every golden case that exercises it. The alternative (implementing
only the fully-confirmed benign-region-overlap category) was available
and declined in favor of covering DMD's actually-relevant mechanism.

What is implemented, in `models/cnv_duplication_evidence.py` +
`cnv_scoring.py`'s `score_cnv_duplication()`:

  - **`GAIN_2K_EQUIV`** (0.45 pts, this project's own disclosed label, not
    a verified Riggs code): a confirmed-tandem duplication with a
    breakpoint inside the gene, predicted out-of-frame.
  - **`GAIN_2J_EQUIV`** (0 pts): same, but in-frame or frame effect
    unknown.
  - **`GAIN_BENIGN`** (-1.0 pts, point value confirmed directly from
    ClassifyCNV's `assign_dup_points_s2()`): the duplication falls
    completely within an established benign copy-number region.
  - **`NONE_APPLICABLE`** (0 pts): whole-gene duplication (TS scoring not
    implemented, see above), no breakpoint inside the gene, or a
    breakpoint inside the gene whose tandem/direct orientation is not
    confirmed.

One real, striking consequence, visible in the golden cases: because 0.45
points falls well short of the 0.90 Likely Pathogenic cutoff, NO
DMD duplication can reach Likely Pathogenic or Pathogenic through this
project's Section-2-only gain scoring alone, no matter how clearly
out-of-frame and disruptive it is -- a real ClinVar-Pathogenic exon 2
duplication (`DMD_CNV_dup_ex2`, see below) scores VUS here. This is a
substantially larger, more conservative gap than the deletion side's
equivalent (2E's 0.9 points at least reaches Likely Pathogenic) --
disclosed rather than smoothed over, because the gain-side point value
itself is this project's own inference, not an independently verified
number, and a bigger disclosed gap is more honest than a confident-looking
score built on an unconfirmed value.

**Curated duplication fixtures**
(`data/curated/cnv_duplication_evidence.json`, 2 cases) and **golden
cases** (`validation/golden_cases/cnv_duplication_golden_cases.yaml`):

  - `DMD_CNV_dup_ex2` -- REAL: ClinVar RCV000240212.2 / VCV000254069
    (`NM_004006.2(DMD):c.32-?_93+?dup62`), Labcorp Genetics (formerly
    Invitae), germline classification Pathogenic (1 star, evaluated
    2016-12-03). Submitter-confirmed likely-tandem duplication of DMD
    exon 2, "results in an absent or disrupted protein product."
    Scores `GAIN_2K_EQUIV` (0.45 pts) -> **VUS** here.
  - `DMD_CNV_dup_ex42_inframe` -- literature-grounded (not a single pinned
    ClinVar accession): the real, published in-frame DMD exon 42
    duplication, PCR-confirmed tandem (not inverted). Scores
    `GAIN_2J_EQUIV` (0 pts) -> **VUS**.

**Tests:** `test_cnv_scoring.py` grew by 17 tests (batch 24 additions) --
golden-case cross-check against both curated duplication fixtures, plus
hand-built edge cases for every decision branch (benign overlap, whole-
gene duplication, no breakpoint, not-tandem, tandem-unknown-orientation,
both frame-effect outcomes) and every `CnvDuplicationEvidence` schema-
validation rule.

Explicitly deferred, same spirit as batch 23: whole-gene triplosensitivity
scoring (no real DMD data would exercise it honestly); independent
verification of the primary Riggs paper's exact gain-side letter codes
and point-value-to-condition mapping; and everything already deferred on
the deletion side (Sections 1, 2H, 3, 4, 5).

**PS3 and BS3: functional-evidence criteria (batch 25).** Richards et al.
2015 (Table 3) defines **PS3** ("well-established in vitro or in vivo
functional studies supportive of a damaging effect on the gene or gene
product," Strong) and **BS3** ("well-established in vitro or in vivo
functional studies show no damaging effect on protein function or
splicing," Strong) but leaves "well-established" entirely to curator
judgment -- no further structure. The ClinGen SVI Working Group's own
refinement, Brnich et al. 2019 (*Genome Medicine* 11:98, "Recommendations
for application of the functional evidence PS3/BS3 criterion"), replaces
that judgment call with an explicit validation-tier ladder based on how
many pathogenic/benign control variants an assay was validated against,
and whether a formal OddsPath statistic was computed: Supporting (few or
undocumented controls), Moderate (>=11 mixed controls, no formal
statistics), up to Strong (a calculated OddsPath). This project does not
compute OddsPath itself -- `validation_strength` on the new
`FunctionalEvidence` model (`src/variant_classifier/models/functional_evidence.py`)
is a curated fact, the calibrated tier a curator already assigned, exactly
like `ComputationalEvidence`'s "one calibrated call per variant" and
`SameResidueEvidence`'s externally-sourced precedent classification.
Only Supporting/Moderate/Strong are accepted -- Richards et al. 2015
defines no Very-Strong or Stand-Alone tier for PS3/BS3, so `__post_init__`
rejects those values rather than silently overclaiming a strength the
base framework doesn't offer.

`assay_result` is a three-way `FunctionalAssayResult` enum
(`ABNORMAL`/`NORMAL`/`INDETERMINATE`), not a boolean. `INDETERMINATE` is a
real, distinct, and — per this batch's own fixture research — common
state: an assay that was performed but did not clearly discriminate
pathogenic from benign for this specific variant. It is curated
explicitly (never left absent, which instead means "no functional data at
all," yielding `NOT_EVALUATED`), and results in `NOT_MET` for *both* PS3
and BS3 rather than guessing a direction. `evaluate_ps3()`/`evaluate_bs3()`
(`src/variant_classifier/evaluators/ps3.py`/`bs3.py`) are structural
mirror images of each other: `NOT_EVALUATED` if no `functional_evidence`
is recorded, `NOT_MET` on `INDETERMINATE` or the opposite-direction
result, `MET` at the curated `validation_strength` on a matching result.
Wiring both into `engine.evaluate_all()` required no changes to
`combine()`'s Table 5 logic or to `bayesian.py`'s point-summing at all --
both already operate generically over any `CriterionResult` list by
counting strength tiers per direction, exactly as the module docstrings
for both files already claimed before this batch put it to the test.

Two real fixture enrichments in `data/curated/variant_evidence.json`
(hand-derived against the new evaluators before being written into golden
cases, same discipline as every other batch): `CAPN3_c.550del` (the real,
ClinVar-Pathogenic founder LGMDR1 allele) gets `ABNORMAL`/`SUPPORTING`,
citing a Czech LGMD2A cohort Western blot study (Chrobakova et al. 2004,
*Neuromuscul Disord*; Hermanova et al. 2006, *Muscle Nerve*) reporting
total absence of calpain-3 protein in patients carrying this allele.
Supporting (rather than a higher tier) is a disclosed conservative choice
-- the primary papers were reCAPTCHA-blocked during curation, so the
exact validation-control counts Brnich et al. 2019 would use to justify
Moderate or Strong could not be independently verified. `CAPN3_c.2257G>A`
(the already-extensively-documented, genuinely contested p.Asp753Asn
variant -- see "Expanding the curated set" below) gets `INDETERMINATE`,
citing the *same* 2025 Bruno et al. paper already cited for that fixture:
its own Table 1 reports mixed Western blot results for this exact variant
(Normal in 3 of 5 patients tested, Reduced in 2), and the paper explicitly
warns that even a Normal result cannot be trusted as benign for this gene
("approximately 20% of individuals with pathogenic CAPN3 variants may have
normal protein levels despite functional impairment"). This is a genuine
real-world illustration of the Brnich framework's INDETERMINATE-outcome
caveat, not a curation shortfall -- and it changes nothing about this
fixture's already-VUS classification, which continues to rest on the same
schema-level single-inheritance-pattern gap documented there. No real
"clean, non-caveated normal-WB-confirms-benign" CAPN3/DMD fixture was
found this round, so BS3's MET branch is covered only by a hand-built
unit test, the same treatment CNV categories 2C/2D/2F got in batch 23.

Adding real PS3 evidence to `CAPN3_c.550del` also created this project's
second real Table-5-vs-Bayesian divergence case (see "Milestone 5"
below): PVS1 Very Strong + PS3 Supporting is "1 Very Strong + 1
Supporting," a shape Table 5 has no combining rule for (still VUS), while
Tavtigian et al. 2020's point system reaches 8+1=9 points, squarely
Likely Pathogenic. `CAPN3_c.1939G>T` showed this same shape already (with
PM2 as the Supporting criterion); this fixture now joins it as a second,
independent real-data example. The `test_bayesian_diverges_from_table5_for...`
regression lock in `test_bayesian.py` was updated (renamed to "five", not
"four", documented fixtures) to catch any future evaluator change that
silently adds or removes a divergence.

Tests: `tests/unit/test_ps3_bs3_evaluators.py` (15 tests) -- golden-case
cross-check plus hand-built edge cases for every branch (`NOT_EVALUATED`,
both `NOT_MET` shapes, `MET` at all three strength tiers for PS3, and a
hand-built BS3 `MET` case since no real fixture reaches it) and
`FunctionalEvidence`'s own schema validation (required/forbidden
`validation_strength` per `assay_result`, rejected Very-Strong/Stand-Alone
tiers).

**PM3: implemented (batch 28).** PM3 ("detected in trans with a
pathogenic variant, for recessive disorders," Table 3) was researched and
explicitly sized-but-not-implemented in batch 25 -- `clinical.py`'s own
docstring had flagged, since an earlier milestone, that a per-variant PM3
evaluator is structurally circular (variant A's PM3 would depend on
variant B's classification from this same engine), and batch 25's
research found the real mechanism avoiding that circularity is not a
simple two-variant workaround: PM3 is a **points-based system that
aggregates evidence across multiple probands**, with points per proband
set by phasing confidence and the other allele's own independently-
established classification, summed to a total, homozygous observations
capped at a maximum of 1 point, and phase required to be either directly
confirmed or ruled out via a gnomAD variant co-occurrence check.

Batch 28 came back to this and built it, but not by aggregating live
across `ClinicalCase`/`CaseInterpretation` the way the batch-25 sizing
implied it might have to. The primary source for the *exact* points-per-
scenario table -- the ClinGen SVI "Recommendation for the in trans
Criterion (PM3)" Version 1.0 PDF
(`clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf`)
and its doc page -- returned empty/unreadable on every fetch attempt this
batch, the same failure mode this project has hit repeatedly for
ClinGen SVI primary documents (Walker et al. 2023's splicing paper and
CAPN3 flowchart, batch 26). Rather than hardcode a recalled-but-unverified
numeric table, or block on an unreachable source a third time, this batch
resolved the circularity a different, more direct way: `Pm3Evidence` /
`Pm3ProbandObservation` (`models/pm3_evidence.py`) attach directly to
`VariantEvidenceBundle`, the same as `FunctionalEvidence` (batch 25), and
the partner allele's classification (`other_allele_classification`, must
be Pathogenic or Likely Pathogenic) plus each proband's points are
**curated facts** rather than values this engine re-derives live from
another variant's own classification -- exactly the same "never silently
guess, state the decision-relevant fact explicitly" pattern already used
for `FunctionalEvidence.validation_strength`. Real-world PM3 curation
already works this way in practice: a curator citing a published
compound-heterozygous case report is trusting that report's own stated
classification of the partner allele, not re-running this project's own
12-evaluator engine on a variant that may not even be in its curated set.

What IS enforced directly, because it WAS confirmed and quoted (ACGS 2024
UK Practice Guidelines for Variant Classification): the homozygous
1-point cap (`Pm3ProbandObservation.__post_init__` rejects any
`HOMOZYGOUS` observation with `points > 1.0`), and the cis-cooccurrence
override ("PM3 should not be applied at any level in the context of two
variants that predominantly co-occur" -- `Pm3Evidence.cis_cooccurrence_observed`
forces `evaluate_pm3()` straight to `NOT_MET`, checked before any
proband-point summing at all). Threshold bands for summed points (Very
Strong >=4, Strong >=2 but <4, Moderate >=1 but <2, Supporting >=0.5 but
<1) are the CAPN3 LGMD VCEP's own real, confirmed, gene-specific table
(`cspec.genome.network/cspec/ui/svi/doc/GN187`) -- the same real source
already used for CAPN3's PM2/BA1/BS1 thresholds. `evaluators/pm3.py` is
this project's twelfth evaluator, and PM3 is this project's **first
Strong-strength pathogenic-direction evaluator** (PVS1 is Very Strong;
PM2/PM4 are Moderate or Supporting; PP3 is Supporting) -- see
`bayesian.py`'s module docstring for the previously-undemonstrated
Table-5-vs-Bayesian discrepancy this newly makes reachable in principle
(2 Strong criteria = 8 Bayesian points, below the 10-point Pathogenic
threshold, a real published Tavtigian et al. 2020 inconsistency this
project still has no fixture combination that actually reaches).

`clinical.py` (`ClinicalCase`/`CaseInterpretation`) is unchanged and
remains complementary, not overlapping: it still answers "does this
patient's genotype explain their disease," which needs case-level facts
(this specific patient's phase, this specific patient's karyotypic sex)
PM3 itself doesn't -- PM3 answers a narrower, variant-level question
("has this variant been seen in trans with a pathogenic variant across
the literature") using facts a curator states once per variant, not once
per patient encounter.

**Real fixture: `CAPN3_c.550del` reaches PATHOGENIC.** This batch
enriched the real founder-allele fixture with two independent, real,
published-cohort HOMOZYGOUS PM3 observations (Czech: Chrobakova et al.
2004 / Hermanova et al. 2006, already cited for this variant's PS3
evidence; Polish: "The Frequency of c.550delA Mutation of the CANP3 Gene
in the Polish LGMD2A Population"), each capped at the real 1.0-point
ceiling, summing to 2.0 = PM3 Strong. Combined with the existing PVS1
Very Strong, Table 5's "1 Very Strong + >=1 Strong" pathogenic rule is
satisfied directly -- this fixture moves from VUS (batches 3-27) to
**PATHOGENIC**, finally matching its real ClinVar call, and closing a gap
this project's own notes have flagged since Milestone 3 ("The real-world
Pathogenic classification for this variant also draws on segregation/
case-count evidence (PM3, PS4-adjacent) that is entirely out of scope for
this project's 6-criterion engine" -- no longer true). It also drops out
of the Table-5-vs-Bayesian divergence set (both systems now agree,
13 Bayesian points, Pathogenic) -- see
`test_bayesian_diverges_from_table5_for_exactly_the_four_documented_fixtures`,
shrunk back from five to four fixtures this batch.

Two new synthetic fixtures (`CAPN3_SYNTH_PM3_MODERATE_01`,
`CAPN3_SYNTH_PM3_CIS_OVERRIDE_01`) cover the `COMPOUND_HETEROZYGOUS`
zygosity path and the cis-cooccurrence override respectively, since no
real curated fixture in this project's set has documented compound-
heterozygous phase confirmation or a gnomAD co-occurrence analysis on
record. VERY_STRONG and SUPPORTING/NOT_MET-below-threshold are covered by
hand-built unit tests in `tests/unit/test_pm3_evaluator.py` (18 tests
total), since no fixture happens to land there either -- the same
"golden-case coverage plus hand-built edge cases for the rest" pattern
every other evaluator's tests already follow.

**X-linked female/other-karyotype case interpretation (batch 29).**
`clinical.py` had handled exactly one X-linked shape since Milestone 4:
a hemizygous male (karyotypic_sex=XY) with a single variant. Every other
karyotypic sex was deferred to MANUAL_REVIEW with an explicit, disclosed
rationale — female X-linked carrier interpretation depends on
X-inactivation biology this project didn't model, and batch 12's own
research (Brioschi et al. 2012, BMC Med Genet 13:73, cited in
`CASE_DMD_FEMALE_CARRIER_REAL`) had already found skewed X-inactivation
in only 2 of 6 symptomatic DMD carriers and 5 of 11 asymptomatic ones —
real, published evidence that the field itself cannot reliably predict a
heterozygous carrier's phenotype from genotype plus X-inactivation
pattern alone. Batch 29's research confirmed that finding is still real
and still the right reason to leave the single-heterozygous-variant case
at MANUAL_REVIEW — but also found a second, genuinely different, real
mechanism this project's existing evidence model could newly represent:
**biallelic** X-linked involvement in an XX individual (both DMD copies
affected, via homozygosity or compound heterozygosity), documented in
real, independently published cases (Ulm et al. 2022, *Molecular
Genetics & Genomic Medicine* 11:e2088, DOI 10.1002/mgg3.2088 — "the
fifth published case ... of a female with multiple DMD variants
confirmed in trans"; Fujii et al. 2009, *Am J Med Genet A* 149A:1052,
DOI 10.1002/ajmg.a.32808, homozygous in-frame deletion, Becker
phenotype; Takeshita et al. 2017, biallelic exon 48–50/51–53 deletions,
DMD phenotype). Crucially, this mechanism does **not** depend on
X-inactivation the unpredictable way the single-carrier case does: when
both X copies carry a qualifying variant, every cell's active X is a
qualifying one, regardless of which X gets randomly inactivated in that
cell — there is no genuinely functional copy anywhere for inactivation
to preferentially spare. That is the real, biological reason this branch
can be resolved confidently while the single-carrier branch still can't,
not an inconsistency between the two.

`ClinicalCase` itself needed no change — it already generically allowed
one or two `variant_ids` with phase required whenever there are two (the
Milestone-4 design was gene/inheritance-agnostic on purpose, see
`models/clinical_case.py`'s docstring). Only `interpret_x_linked_case`'s
own hardcoded "exactly one variant" check needed relaxing, and only for
`karyotypic_sex=XX`: `karyotypic_sex=XY` now explicitly *rejects* two
variant_ids (`SchemaValidationError` — a hemizygous individual has only
one X chromosome, so two variant_ids for an X-linked gene isn't a real
genotype, not something to silently misinterpret), and `OTHER`/`UNKNOWN`
remain exactly as narrow as before (single variant only, always
MANUAL_REVIEW) since those karyotypes are too heterogeneous to reason
about generically — batch 29's research also surfaced that X0/Turner
carriers are functionally hemizygous (like XY) while XXY carriers are
diploid-X (like XX), and `OTHER` doesn't distinguish which, so extending
it would mean guessing which sub-case applies.

The new `_interpret_xx_biallelic` branch resolves EXPLAINED only for the
narrow, XCI-independent case: phase confirmed TRANS *and* both variants
independently classified Pathogenic or Likely Pathogenic. Every other
two-variant XX combination — CIS, UNKNOWN phase, or TRANS-but-not-both-
qualifying — stays MANUAL_REVIEW, and deliberately does **not** get the
same treatment autosomal recessive's equivalent branches get (CIS and
non-qualifying resolve to INSUFFICIENT there). The reasoning, worked out
and confirmed before writing any code: autosomal recessive CIS leaves a
genuinely wild-type *autosome* copy, which is always transcriptionally
active — a real protective fact. X-linked CIS (or a non-qualifying
partner variant) also leaves a genuinely wild-type *X* copy, but X
chromosomes are not always active — X-inactivation randomly silences one
X per cell, so that wild-type copy can still be silenced in some cells,
reintroducing exactly the same real, published phenotype
unpredictability (Brioschi et al. 2012 again) as the single-heterozygous-
variant case. Getting this asymmetry right — not just mechanically
mirroring `interpret_recessive_case`'s CIS/benign-side handling — was the
main design decision confirmed with the user before coding.

No real published case pairs two specific point-mutation DMD variants in
one biallelic-XX patient — every real case found this batch is CNV-based
or a CNV+point-mutation mix (the Ulm/Fujii/Takeshita citations above),
which doesn't fit this project's point-mutation-only `ClinicalCase`
model (CNV scoring, batches 23/24, is a deliberately separate, parallel
system with no case-level wiring). Two new fixtures
(`CASE_DMD_XX_BIALLELIC_TRANS`, `CASE_DMD_XX_BIALLELIC_CIS`) are
therefore synthetic at the case level, but built entirely from two REAL,
already-curated DMD variants (`DMD_c.2302C>T`, `DMD_c.8944C>T`, both
real ClinVar-grounded, both LIKELY_PATHOGENIC) rather than inventing new
synthetic point mutations — the same "synthetic case, real variant
evidence" pattern `CASE_DMD_FEMALE_CARRIER_REAL` already established in
batch 12. `UNKNOWN` phase and TRANS-but-not-both-qualifying are covered
by hand-built unit tests only (`test_clinical.py`), since no fixture
happens to land there either. 5 new tests plus the golden-case checks
already covered by the shared `test_interpret_case_matches_golden_case_for_all_curated_cases`
loop, 250 total.

**Splice-RNA evidence feeds PVS1 directly (batch 26).** Batch 26's
research into extending PVS1 (Abou Tayoun et al. 2018's full decision
tree branches on protein-domain criticality and constitutive-exon-
splicing information this project has never modeled — see "PVS1 scope"
above) turned up a real, useful, previously-unknown-to-this-project fact:
the ClinGen LGMD VCEP's own CAPN3 specification (v2.0,
cspec.genome.network/cspec/ui/svi/doc/GN187, fetched directly) states
that PS3 is "not applicable at this time" for CAPN3 in vitro functional
assays, and that "for any variant type, experimental evidence for altered
splicing should be scored under PVS1 in accordance with the decision tree
for RNA splicing assay results outlined in Walker et al. 2023 (PMID:
37352859)" — i.e. real RNA/splicing assay evidence sets PVS1's own
strength directly for this gene, it does not feed a separate PS3
criterion the way batch 25 might have suggested. This reframes
"experimental splicing evidence" as PVS1-domain evidence, not
functional-evidence-domain evidence, and this batch implements it that
way rather than folding it into `FunctionalEvidence`.

Both the primary Walker et al. 2023 paper (PMC, ScienceDirect, and
Cell.com all reCAPTCHA-blocked or JS-rendered-empty) and the CAPN3-
specific PVS1 flowchart PDF the VCEP spec attaches (a binary file this
project's fetch tooling cannot parse as text, and which this project's
rules do not permit downloading via a bash-level workaround) were
unreachable during this batch's research. So the exact percentage-of-
transcript and protein-region-criticality thresholds that govern most of
Walker et al. 2023's decision tree remain unimplemented — this is
NOT a full implementation of that tree. What is implemented is narrower
and threshold-free by construction: a new `SplicingRnaEvidence` enum
(`src/variant_classifier/models/enums.py`) and a new optional
`TranscriptConsequence.splicing_rna_evidence` field (restricted to
`SPLICE_DONOR_VARIANT`/`SPLICE_ACCEPTOR_VARIANT`, unlike `nmd_predicted`/
`repeat_region` this field is NOT required whenever its consequence class
applies — most splice fixtures simply have no RNA assay at all, the
common, unset default case, which falls through to the original,
unchanged predicted-only `MANUAL_REVIEW` path):

- `CONFIRMED_NULL_EQUIVALENT`: a real assay confirms the aberrant
  transcript is functionally equivalent to a null allele (out-of-frame
  exon skip, intron retention producing a frameshift, or a confirmed PTC
  with no normal transcript detected) — treated identically to a
  confirmed-NMD frameshift/nonsense variant — `evaluate_pvs1()` returns
  MET, Very Strong. No percentage threshold needed: the assay directly
  establishes the null-equivalent outcome.
- `CONFIRMED_NORMAL_SPLICING`: a real assay directly contradicts the
  predicted splice disruption — NOT_MET, not a guess in either remaining
  direction, since the null-variant mechanism this consequence class
  assumes is refuted by direct evidence. Also threshold-free.
- `CONFIRMED_IN_FRAME_OR_PARTIAL_FUNCTION` and `INCONCLUSIVE` both
  return MANUAL_REVIEW — the same open protein-domain-criticality
  question (for the in-frame case) or genuine ambiguity (for the
  inconclusive case) that the rest of this evaluator already declines to
  guess through.

Real fixture research surfaced a specific, honest limitation: none of the
three real curated splice fixtures (`CAPN3_c.946-1G>A`,
`CAPN3_c.2050+1G>A`, `DMD_c.93+1G>A`) has experimental data precise
enough to populate `CONFIRMED_NULL_EQUIVALENT`. `CAPN3_c.946-1G>A` came
closest — its real ClinVar submitter comment cites "Experimental studies
have shown that this variant disrupts mRNA splicing (PMID: 7720071)," a
genuine RNA-level study, not just a computational canonical-splice-site
prediction — but the available secondary source doesn't state the
resulting transcript's frame or NMD status, so it's now curated as
`INCONCLUSIVE` (a real, honest distinction from "no experimental data at
all," even though both currently route to the same MANUAL_REVIEW status).
The other two fixtures' records mention only the predicted consequence,
so `splicing_rna_evidence` is left unset for them, unchanged. The new
MET/Very-Strong path is instead demonstrated end-to-end on one disclosed
synthetic fixture, `CAPN3_SYNTH_PVS1_SPLICE_RNA_01` — modeled on a
hypothetical minigene/RT-PCR assay — the same "searched, not found, still
open, demonstrate on a labeled synthetic instead" treatment PS1's MET
path got in batch 22. Its `population_evidence` is deliberately
`NOT_ASSESSED` so PVS1 MET stands as its only contributing criterion,
keeping the fixture focused on the new branch alone rather than
incidentally becoming a sixth Table-5-vs-Bayesian divergence case (that
pattern is already well covered by real fixtures — see "Milestone 5"
above).

Tests: `tests/unit/test_pvs1_evaluator.py` grew from 12 to 19 tests --
one new golden-case-covered fixture plus six hand-built edge cases (all
four `SplicingRnaEvidence` branches, confirmation that `START_LOST` is
unaffected, and `TranscriptConsequence`'s new validation rule rejecting
`splicing_rna_evidence` on non-splice consequence types). `engine.py`
and `bayesian.py` needed no changes at all — same story as every prior
criterion addition, since both operate generically over whatever
`evaluate_pvs1()` returns.

**PVS1 start-loss: the alternative-start-codon rule (batch 27).** Where
batch 26's splice-RNA-evidence research hit blocked primary sources,
batch 27's start-loss research did not: the exact rule for PVS1's
initiation-codon branch is quoted directly in the ACGS 2024 UK Practice
Guidelines for Variant Classification (already a trusted, successfully-
fetched source in this project), citing Abou Tayoun et al. 2018's
decision tree verbatim: "If there is a potential in-frame initiation
codon downstream, the missing N-terminal region of the protein should be
assessed according to the principles described in the decision tree (i.e.
is the missing region critical to protein function / is it >10% of the
entire protein length / are there any reported pathogenic variants
upstream of the potential initiation codon) and apply PVS1 at either
reduced strength or n/a, as appropriate. If no alternative in-frame start
codon is identified, use PVS1 at maximum strength." A second source (a
web search summary of Abou Tayoun et al. 2018 itself) supplied the exact
downgrade rule: "the only null variant that should be applied the
downgraded criteria PVS1_Supporting is a variant of initiation codon when
a new methionine is not preceded by a pathogenic variant."

Three new optional `TranscriptConsequence` fields, restricted to
`START_LOST` and validated together the same way `nmd_predicted`/
`repeat_region`/`splicing_rna_evidence` already are:
`alternative_start_codon_identified` (bool), and, required together
whenever that's `True`, `alternative_start_codon_percent_protein_lost`
(float, 0-100) and `alternative_start_codon_preceded_by_pathogenic_variant`
(bool). `evaluate_pvs1()`'s new logic follows the quoted rule directly:
`identified=False` -> MET, Very Strong (no rescue possible); `identified`
unset -> MANUAL_REVIEW (unchanged default); `identified=True` with
`percent_protein_lost>10` or `preceded_by_pathogenic_variant=True` ->
MANUAL_REVIEW (falls outside the automatic downgrade, needs a
protein-domain-criticality judgment this evaluator still doesn't make);
`identified=True` with neither disqualifying factor -> MET, Supporting.

Unlike batch 26, this batch's real fixture (`CAPN3_c.1A>G`) could be
fully, independently resolved -- not to MET, but to a far more specific
MANUAL_REVIEW than before. NM_000070.3's real CDS sequence was fetched
directly from Ensembl's REST API (matched to the correct canonical
transcript, ENST00000397163, by confirming its translation-start genomic
coordinate, 42359806, is exactly this variant's own position) and scanned
programmatically for the first downstream in-frame ATG: codon 228,
meaning re-initiation there would lose the protein's first 227 of 821
residues -- 27.6%, independently computed from primary sequence data, not
estimated or guessed. That alone exceeds the >10% threshold. It is also
independently corroborated by this project's own curated fixture set,
without needing any new external lookup: `CAPN3_c.550del` (p.Thr184fs,
real, ClinVar Pathogenic) sits at residue 184, within the very region
(1-227) that would be lost -- direct, real evidence that region is not
dispensable padding. Both factors agree, so `CAPN3_c.1A>G` now
demonstrates the "outside the automatic downgrade" branch specifically,
with a rationale that names exactly which factors triggered it, rather
than a generic "not implemented" message.

Neither of the other two new branches had a real fixture available:
`CAPN3_c.1A>G`'s own alternative start codon is disqualified, and no
other real CAPN3/DMD start-loss variant is currently curated. Both are
instead demonstrated on new disclosed synthetic fixtures --
`CAPN3_SYNTH_PVS1_STARTLOSS_NO_ALT_01` (no alternative start codon at
all -> MET Very Strong) and `CAPN3_SYNTH_PVS1_STARTLOSS_SUPPORTING_01`
(a clean, small, unencumbered alternative start codon -> MET Supporting)
-- the same "searched, not found, still open, demonstrate on a labeled
synthetic instead" treatment used throughout this project. The first is
deliberately built with `population_evidence` `NOT_ASSESSED` so its PVS1
MET stands alone, avoiding another incidental Table-5-vs-Bayesian
divergence (same discipline as `CAPN3_SYNTH_PVS1_SPLICE_RNA_01`, batch
26); the second uses ordinary ABSENT population evidence and lands on
"PVS1 Supporting + PM2 Supporting = 2 Supporting alone," a combination
neither Table 5 nor the Bayesian point system has a pathogenic rule for,
so it stays VUS under both without needing special construction.

What remains open, disclosed rather than silently narrower than it looks:
the real protein-domain-criticality judgment inside the MANUAL_REVIEW
branch (is the lost region *actually* functionally critical, not just
">10%" or "a pathogenic variant happens to be there") is still not
modeled -- the same boundary the NMD-escape and in-frame-splice branches
already have. The >10% cutoff's interpretation (`>10.0` triggers review,
`<=10.0` is eligible for the automatic path) is this project's own
disclosed reading of the quoted threshold, not verbatim from a source
that specifies which side of exactly 10% falls where.

Tests: `tests/unit/test_pvs1_evaluator.py` grew again, from 19 to 27
tests -- three new golden-case-covered fixtures plus eight hand-built
edge cases (all four start-loss outcomes, including that a pathogenic
variant upstream disqualifies the downgrade even when the percentage
alone would not, and three new `TranscriptConsequence` validation rules).
`engine.py`/`bayesian.py` again needed zero changes.

**PM2 and founder mutations.** PM2 asks whether a variant is absent or at
extremely low frequency in the general population. A single global allele
frequency threshold isn't enough to answer that safely: `CAPN3_c.550del` is
rare overall (0.023%) but a known founder mutation enriched to 0.75% in
specific ancestries. The evaluator does not silently pass PM2 using the
lower, reassuring global number — when an ancestry-specific frequency
clears the threshold while the overall frequency doesn't, it returns
MANUAL_REVIEW rather than guessing, because whether "extremely low" holds
depends on the tested individual's ancestry, which isn't available here.

**Dataclasses instead of pydantic.** All eight models use the Python
standard library's `dataclasses` module with hand-written `from_dict()`
validation rather than pydantic. This keeps the dependency footprint to
just PyYAML for fixture loading. Converting to pydantic later, if its
validation machinery becomes useful, is a contained, mechanical change
scoped to these eight files.

## Repository layout

```
src/variant_classifier/
  errors.py                  SchemaValidationError — the one exception type
  models/
    enums.py                 controlled vocabularies + ACMG_CRITERION_CODES
    _coerce.py                shared from_dict() validation helpers
    variant_identity.py        VariantIdentity
    gene_disease_context.py    GeneDiseaseContext, Specification
    transcript_consequence.py  TranscriptConsequence — splicing_rna_evidence field added batch 26,
                                alternative_start_codon_* fields added batch 27, see "PVS1 scope" above
    population_evidence.py     PopulationEvidence
    computational_evidence.py  ComputationalEvidence
    same_residue_evidence.py   SameResidueEvidence — batch 22, see "PS1 and PM5" above
    functional_evidence.py     FunctionalEvidence — batch 25, see "PS3 and BS3" above
    pm3_evidence.py             Pm3Evidence, Pm3ProbandObservation — batch 28, see "PM3: implemented (batch 28)" above
    criterion_result.py        CriterionResult
    provisional_classification.py  ProvisionalClassification
    cnv_deletion_evidence.py    CnvDeletionEvidence — batch 23, see "DMD CNV/structural-variant scoring" above
    cnv_duplication_evidence.py CnvDuplicationEvidence — batch 24, see "DMD CNV/structural-variant duplication scoring" above
    cnv_category_result.py      CnvCategoryResult — batch 23/24, shared by both deletion and duplication scoring
    cnv_provisional_classification.py  CnvProvisionalClassification — batch 23/24, shared by both
    evidence_bundle.py         VariantEvidenceBundle (container, this repo only)
    golden_case.py             GoldenCase (container, this repo only)
  loader.py                  loads/validates the curated fixtures below
  engine.py                  evaluate_all() + combine() + classify() — Table 5 combining engine
  bayesian.py                 combine_bayesian() + classify_bayesian() — Milestone 5, see above
  cnv_scoring.py               score_cnv_deletion() (batch 23) + score_cnv_duplication() (batch 24), see above
  clinical.py                 interpret_case() — case-level reasoning, see "Case-level scope" above
  evaluators/
    pvs1.py                   evaluate_pvs1() — see "PVS1 scope" above; splice-RNA-evidence branch added batch 26
    pm2.py                    evaluate_pm2()
    pm4.py                    evaluate_pm4() — see "PM4: a second new criterion" above
    ps1.py                    evaluate_ps1() — batch 22, see "PS1 and PM5" above
    pm5.py                    evaluate_pm5() — batch 22, see "PS1 and PM5" above
    pm3.py                    evaluate_pm3() — batch 28, see "PM3: implemented (batch 28)" above
    ps3.py                    evaluate_ps3() — batch 25, see "PS3 and BS3" above
    pp3.py                    evaluate_pp3()
    bp4.py                    evaluate_bp4()
    ba1.py                    evaluate_ba1()
    bs1.py                    evaluate_bs1()
    bs3.py                    evaluate_bs3() — batch 25, see "PS3 and BS3" above

config/
  population_thresholds.yaml per-gene PM2/BA1/BS1 frequency thresholds -- CAPN3's
                              are the real ClinGen LGMD VCEP values as of batch 4
                              (see Design notes); DMD's remain generic ACMG defaults
  dosage_sensitivity.yaml    per-gene ClinGen Haploinsufficiency curation -- batch 23,
                              DMD only (HI=3); CAPN3 deliberately absent (autosomal
                              recessive, see "DMD CNV/structural-variant scoring" above)

data/
  curated/
    gene_disease_context.yaml   CAPN3 and DMD
    variant_evidence.json       29 curated variants (CAPN3 + DMD) -- growing toward ~20-30;
                                 two enriched with functional_evidence in batch 25 (see "PS3 and BS3" above);
                                 one new synthetic PVS1 splice-RNA-evidence fixture added batch 26;
                                 CAPN3_c.1A>G enriched, plus two new synthetic PVS1 start-loss fixtures,
                                 added batch 27; CAPN3_c.550del enriched with real PM3 evidence, plus two
                                 new synthetic PM3 fixtures, added batch 28 (see "PM3: implemented
                                 (batch 28)" above)
    clinical_cases.json         11 curated ClinicalCase fixtures (Milestone 4) -- 2 added batch 29 for biallelic XX X-linked
    cnv_deletion_evidence.json  3 curated CNV deletion fixtures (DMD only) -- batch 23,
                                 a separate curated set from variant_evidence.json above
    cnv_duplication_evidence.json  2 curated CNV duplication fixtures (DMD only) -- batch 24,
                                 a separate curated set again (see "DMD CNV/structural-variant
                                 duplication scoring" above for why it isn't folded into the
                                 deletion set)
  source/
    pipeline_annotate_calls/  real (trimmed) sample of CAPN3-DMD-variant-calling-pipeline's
                              ANNOTATE_CALLS VCF output, used by test_pipeline_adapter.py
  synthetic/                 placeholder — larger generated datasets (empty)

src/variant_classifier/
  pipeline_adapter.py         build_bundles_from_pipeline_output() — CAPN3-DMD-variant-calling-pipeline's
                              VCF -> VariantEvidenceBundle (Batch 21, see Status above)

validation/golden_cases/
  variant_golden_cases.yaml            expected per-variant results, Table 5 (renamed from
                                        capn3_milestone1.yaml once DMD variants existed)
  variant_golden_cases_bayesian.yaml   expected per-variant results, Bayesian (Milestone 5)
  case_interpretation_golden_cases.yaml expected case-level results (Milestone 4)
  cnv_deletion_golden_cases.yaml        expected CNV deletion results (batch 23), curated
                                        separately from cnv_deletion_evidence.json above
  cnv_duplication_golden_cases.yaml     expected CNV duplication results (batch 24), curated
                                        separately from cnv_duplication_evidence.json above

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
extra environment variables. All 160 tests currently pass.

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

Complete historical log, in order, of every milestone and batch that
built this project — see "Final status and scope boundary" at the end
of this section for where things stand now.

- **Milestone 2** — done. PM2 and PVS1 evaluators (PVS1 intentionally
  partial).
- **Milestone 3** — done. BA1/BS1/PP3/BP4 evaluators and the combining
  engine.
- **Milestone 4** — done. `ClinicalCase`/`CaseInterpretation` models and
  `clinical.py`'s case-level reasoning (see "Case-level scope" above).
- **Milestone 5** — done (batch 20). Bayesian point-based combining
  (`bayesian.py`, Tavtigian et al. 2020), offered alongside Table 5's
  `engine.py` rather than replacing it. See "Milestone 5: Bayesian
  point-based combining" above for the four real divergences found and
  the case-level agnosticism proof.
- **Batch 4** — done. Adopted the real ClinGen LGMD VCEP CAPN3
  specification's PM2/BA1/BS1 thresholds and PM2 strength (see "Real
  ClinGen LGMD VCEP thresholds for CAPN3" above); re-derived every golden
  case this touched.
- **Batch 5** — done. Added `CAPN3_c.1A>G`, a real start-loss variant —
  PVS1's third documented scope gap now has a real fixture, alongside the
  two splice-site variants.
- **Batch 6** — done. Added `CAPN3_c.2120A>G`, sourced from an official
  ClinGen LGMD VCEP Evidence Repository curation — the project's first
  real calibrated computational (REVEL) score, and a real-world
  confirmation of the founder-enrichment MANUAL_REVIEW design (see
  "Expanding the curated set" above).
- **Batch 7** — done. Added `DMD_c.93+1G>A`, the project's second real DMD
  variant (first real DMD splice-site fixture). Also made a further
  unsuccessful attempt at the BA1-common-variant gap (see "Expanding the
  curated set" above) — still open, now searched five separate times.
- **Batch 8** — done. No new fixture, but enriched `CAPN3_c.2257G>A` with
  a real REVEL score (first real INDETERMINATE-prediction fixture) and
  its real LGMD VCEP classification (Likely Benign via BS1+BP2 — a
  case-level gap mirroring `CAPN3_c.1939G>T`'s pathogenic-side one). Sixth
  unsuccessful BA1-gap search round.
- **Batch 9** — done. Added `DMD_c.11041A>T`, a real DMD nonsense variant
  in the final exon-exon junction region — the project's first real
  fixture for PVS1's NMD-escape scope gap (all three documented PVS1
  scope gaps now have real fixtures). Coordinates independently derived
  via the Ensembl REST API, then confirmed against the primary ClinVar
  record. Also added an explicitly-labeled DMD PM2/BS1 placeholder
  threshold to `config/population_thresholds.yaml` (no real DMD VCEP spec
  exists yet), surfaced by this being the first real, non-absent DMD
  population-frequency observation in the set.
- **Batch 10** — done. No new fixture; a seventh unsuccessful BA1 search
  (this time via an Ensembl gene-region variant pull instead of a
  ClinVar-anchored search — didn't pan out either, see "Expanding the
  curated set" above) led back to enriching `CAPN3_c.2257G>A` instead:
  primary-source-confirmed its VCEP Likely-Benign classification and real
  coordinates (previously secondhand), and added a brand-new (Nov 2025)
  peer-reviewed paper arguing the same variant may be pathogenic for a
  distinct autosomal-dominant calpainopathy this project's
  one-inheritance-pattern-per-gene schema cannot represent.
- **Batch 11** — done. Added `DMD_c.8944C>T`, a real, cleanly-sourced DMD
  nonsense variant (LIKELY_PATHOGENIC, same shape as `DMD_c.2302C>T` at a
  different exon). Found while searching unsuccessfully for a real
  BP4-MET fixture, which remains an open gap alongside BA1.
- **Batch 12** — done. No new variant fixture; added the project's first
  real-variant Milestone 4 clinical cases (`CASE_DMD_HEMIZYGOUS_MALE_REAL`
  / `CASE_DMD_FEMALE_CARRIER_REAL`, both on `DMD_c.2302C>T`), replacing
  reliance on synthetic variants for that layer's real-data coverage and
  grounding the female-carrier MANUAL_REVIEW case in real DMD
  manifesting-carrier literature.
- **Batch 13** — done. No new fixture; proved (and locked in with a
  regression test) that no real CAPN3 variant can currently reach
  PATHOGENIC/LIKELY_PATHOGENIC through this engine, a structural
  consequence of CAPN3's real VCEP PM2-Supporting-only threshold plus
  PP3's hardcoded Supporting strength — turning batch 12's open question
  about real biallelic CAPN3 cases into a documented, well-understood
  limitation rather than an unexplained gap.
- **Batch 14** — done. Implemented PM4 (`evaluate_pm4()`), this project's
  first new criterion since Milestone 3 — found while investigating what
  looked like a PVS1 stop-loss gap and turned out to be a missing PM4
  instead (see "PM4: a second new criterion" above). Added its first real
  fixture, `CAPN3_c.1401_1403del`, and extended batch 13's structural
  proof to cover the in-frame indel/stop-loss shape too.
- **Batch 15** — done. Added `CAPN3_c.598_612del`, a second real PM4
  fixture exercising the evaluator's default Moderate-strength branch
  (complementing batch 14's single-residue Supporting-strength one). An
  8th BA1 attempt and a 3rd BP4-MET attempt both hit real tooling walls
  this round (client-rendered search pages, a blocked gnomAD API, a
  temporary fetch rate limit) rather than turning up nothing — both gaps
  remain open and are now documented as tooling-blocked.
- **Batch 16** — done. Added `DMD_c.10103A>G`, this project's first DMD
  missense fixture and first DMD fixture with computational evidence
  (REVEL, cited from the source publication, calibrated against the
  real Pejaver et al. 2022 generic thresholds rather than CAPN3's
  VCEP-specific ones). Notable for matching a real published ACMG
  classification (the discovery paper's own VUS call) directly, not just
  a defensible fallback against a scope gap — see "Expanding the curated
  set" above. Curated variant set reaches the low end of the original
  ~20-30 target (20).
- **Batch 17** — done. No new variant fixture; closed a case-level gap
  instead. `CASE_DMD_HEMIZYGOUS_MALE_VUS_REAL` is the first curated
  `ClinicalCase` to exercise `interpret_x_linked_case`'s MANUAL_REVIEW
  catch-all branch (hemizygous male, non-qualifying/non-benign variant),
  which existed in code since Milestone 4 but had never actually been
  reached by any curated fixture. Uses `DMD_c.10103A>G`'s real hemizygous
  son, whose own real diagnosis was similarly unresolved.
- **Batch 18** — done. Added `DMD_c.5234G>A`, closing the project's
  longest-standing open gap (a real common variant to exercise BA1 MET),
  searched for across nine rounds since Milestone 1. The search had
  stayed CAPN3-scoped the whole time even though the gap was never
  actually CAPN3-specific; found by trying a different ClinVar Miner page
  shape (per-submitter listing, static HTML) after gene-level
  search/listing pages kept failing on client-side rendering. First real
  fixture to land on flat BENIGN. One real-data gap remains open: a real
  BP4-MET example — DMD's real Pejaver calibration (batch 16) opens an
  easier avenue for this than CAPN3's strict VCEP threshold, not yet
  searched.
- **Batch 19** — done. Added `DMD_c.5163G>C`, closing the project's last
  open real-data gap (a real BP4-MET example). REVEL score (0.167)
  sourced from `myvariant.info`'s dbNSFP-backed public API — a new
  working data source discovered this round after NCBI eutils and
  gnomAD's GraphQL API were both confirmed non-functional in this
  environment. Real ClinVar Benign/Likely-Benign consensus across six
  RCVs; BP4 MET alone lands on VUS (Table 5 needs a second
  benign-direction criterion for Likely Benign), the real-data
  complement to `CAPN3_SYNTH_LIKELY_BENIGN_01`'s synthetic BS1+BP4 pair.
- **Batch 20 (Milestone 5)** — done. Added Bayesian point-based combining
  (`bayesian.py`, Tavtigian et al. 2020) as a second, fully tested
  combining system alongside Table 5's `engine.py`. Found four real
  divergences by hand-deriving Bayesian results for all 22 fixtures
  before writing any code: `CAPN3_c.1939G>T` (VUS under Table 5,
  LIKELY_PATHOGENIC under Bayesian — the case batch 4 originally
  flagged as the motivating example) plus a second, independently
  discovered shape in `DMD_SYNTH_PATHOGENIC_01`, `DMD_c.2302C>T`, and
  `DMD_c.8944C>T` (LIKELY_PATHOGENIC under Table 5, PATHOGENIC under
  Bayesian). Proved `clinical.py`'s case-level layer is fully agnostic
  to which combining system produced its input classifications — no
  case-level golden case needed touching. This was chosen deliberately
  as this project's likely last milestone before development moves to a
  complementary project starting one stage earlier in the pipeline (raw
  sequencing reads rather than an already-identified variant) — see
  "Milestone 5: Bayesian point-based combining" above.
- **Batch 21 (Projects 4x5 integration)** — done. Added
  `src/variant_classifier/pipeline_adapter.py`: turns
  CAPN3-DMD-variant-calling-pipeline's `ANNOTATE_CALLS` VCF output (VEP
  transcript consequence + gnomAD v4.1 population frequency, GATK/DeepVariant-
  concordant HG002 calls) into real `VariantEvidenceBundle` instances, reusing
  `loader.load_gene_disease_contexts()` rather than duplicating it. Lives here,
  not in the pipeline repo, so that repo stays free of a runtime dependency on
  this package (see the adapter's own module docstring). Verified against a real,
  hand-trimmed 5-record sample of CAPN3-DMD-variant-calling-pipeline's actual output
  (`data/source/pipeline_annotate_calls/`, see its README for exactly which
  records and why): every VEP/gnomAD field hand-checked against an independent
  `bcftools`/`vep` run matches exactly, multi-allelic sites correctly split
  into one bundle per allele (including reproducing VEP's own indel-allele-
  trimming convention), and records with no MANE-transcript CSQ hit (i.e.
  outside both genes' actual transcript span) correctly produce no bundle
  rather than a fabricated one. `ComputationalEvidence` and NMD/repeat-region
  determination for consequences that need them are explicitly left absent,
  not guessed (see the adapter docstring for why) — neither path is exercised
  by real HG002 data, which is confirmed clinically empty for CAPN3/DMD
  pathogenic variation. 8 new tests (`test_pipeline_adapter.py`), 142 total.
- **Batch 22 (PS1/PM5, and DMD CNV sizing)** — done. Implemented
  `evaluate_ps1()`/`evaluate_pm5()` (`src/variant_classifier/evaluators/ps1.py`,
  `pm5.py`) and the new `SameResidueEvidence` model — the two criteria the
  Roadmap had named "most tractable next." See "PS1 and PM5: same-residue
  precedent evidence" above for the full design (why no new evidence-domain
  model was needed, the splice caveat, the disclosed precedent-strength
  downgrade, and the real ClinGen LGMD VCEP machinery deliberately left
  unimplemented). Added a real fixture pair found directly through this
  research (`CAPN3_c.1342C>T`, and an update to the pre-existing
  `CAPN3_c.1343G>A` once its real classification turned out to have
  changed) and a disclosed synthetic fixture (`CAPN3_SYNTH_PS1_01`) for
  PS1's MET path, since no real PS1 precedent pair was found this round —
  see "Expanding the curated set" above. Also researched and sized (but
  deliberately did not implement) DMD's long-disclosed CNV/structural-
  variant representation gap against the real ClinGen technical standard
  for CNV interpretation (Riggs et al. 2020) and DMD's own ClinGen dosage-
  sensitivity curation and reading-frame-rule literature — see "DMD CNV/
  structural-variant representation" above for the full sizing writeup and
  why it needs a new identity representation and a new parallel scoring
  system, not an incremental evaluator. 18 new tests
  (`test_ps1_pm5_evaluators.py`), 160 total.
- **Batch 23 (DMD CNV deletion scoring)** — done. Implemented the slice
  batch 22 sized: `models/cnv_deletion_evidence.py`, `cnv_category_result.py`,
  `cnv_provisional_classification.py`, `config/dosage_sensitivity.yaml`,
  and `src/variant_classifier/cnv_scoring.py` (`score_cnv_deletion()`),
  covering Section 2 (dosage-sensitivity) of the Riggs et al. 2020 CNV
  rubric for DMD deletions only (whole-gene, 5'/3' end, intragenic
  frameshift+NMD via the Aartsma-Rus reading-frame rule, and benign-region
  overlap) — see "DMD CNV/structural-variant scoring (batch 23)" above for
  the full design, exact point values/sources, and everything explicitly
  deferred. Added 3 curated CNV fixtures (`data/curated/cnv_deletion_evidence.json`):
  one REAL ClinVar-anchored out-of-frame deletion (`DMD_CNV_del_ex47_50`,
  RCV000813350/VCV000656845), and two literature-grounded composite
  fixtures (`DMD_CNV_del_whole_gene_Xp21`, on the real Xp21 contiguous
  gene deletion syndrome; `DMD_CNV_del_ex45_47_inframe`, on the classic
  in-frame Becker deletion) — disclosed as composites rather than single
  pinned ClinVar accessions, since none was found during this batch's
  research. Duplications remain entirely out of scope, confirmed with the
  user before coding. 20 new tests (`test_cnv_scoring.py`), 180 total.
- **Batch 24 (DMD CNV duplication scoring)** — done. Extended batch 23's
  CNV scoring to duplications: `models/cnv_duplication_evidence.py` and
  `cnv_scoring.py`'s `score_cnv_duplication()`. Research surfaced a real
  wrinkle before any code was written -- ClinGen's own DMD dosage curation
  says whole-gene DMD duplications aren't clinically reported (so no
  triplosensitivity scoring is implemented), and the real mechanism (an
  intragenic tandem duplication disrupting the gene, LOF-type, via the
  Aartsma-Rus reading-frame rule again) has a confirmed real category in
  Riggs et al. 2020 whose exact point values could only be partially
  verified from secondary sources (2K=0.45, 2J=0, condition-to-value
  mapping inferred by this project, not confirmed) -- presented to the
  user as an explicit choice before coding; see "DMD CNV/structural-
  variant duplication scoring (batch 24)" above for the full writeup. One
  real, disclosed consequence: since 0.45 points falls short of the 0.90
  Likely Pathogenic cutoff, no DMD duplication can reach Likely
  Pathogenic/Pathogenic through this scoring alone, however clearly
  out-of-frame -- visible directly in `DMD_CNV_dup_ex2`'s golden case (a
  real ClinVar-Pathogenic exon 2 duplication that scores VUS here). Added
  2 curated duplication fixtures (`data/curated/cnv_duplication_evidence.json`):
  one REAL ClinVar-anchored out-of-frame duplication (`DMD_CNV_dup_ex2`,
  RCV000240212/VCV000254069) and one literature-grounded in-frame
  composite (`DMD_CNV_dup_ex42_inframe`). Whole-gene TS scoring and
  independent verification of the primary paper's exact gain-side letter
  codes remain explicit, disclosed gaps. 17 new tests
  (`test_cnv_scoring.py` grew to 37 CNV tests total), 197 total.
- **Batch 25 (PS3/BS3 functional-evidence criteria; PM3 sized)** — done.
  Implemented PS3 and BS3 (`evaluators/ps3.py`/`bs3.py`, new
  `FunctionalEvidence` model) per the base Richards et al. 2015 definition
  and the Brnich et al. 2019 validation-strength ladder — see "PS3 and
  BS3: functional-evidence criteria (batch 25)" above for the full design,
  including the three-way `assay_result` (`ABNORMAL`/`NORMAL`/`INDETERMINATE`)
  and why `combine()`/`bayesian.py` needed zero changes to pick up the two
  new evaluators. Enriched two real curated fixtures with
  `functional_evidence`: `CAPN3_c.550del` (real Czech LGMD2A cohort
  Western blot data, `ABNORMAL`/`SUPPORTING`, disclosed as a conservative
  strength choice) and `CAPN3_c.2257G>A` (real, genuinely mixed Western
  blot data for this exact variant from the same Bruno et al. 2025 paper
  already cited for it, `INDETERMINATE`) — no fabricated BS3-MET fixture
  was added since no real "clean normal-WB-confirms-benign" CAPN3/DMD
  example was found; that branch is covered by a hand-built test only.
  Adding real PS3 evidence to `CAPN3_c.550del` created this project's
  second real Table-5-vs-Bayesian divergence (see "Milestone 5" above),
  locked at the time into `test_bayesian_diverges_from_table5_for_exactly_the_five_documented_fixtures`
  (renamed/shrunk back to "four" in batch 28, once real PM3 evidence gave
  this same fixture a Strong-strength criterion and it stopped
  diverging — see below). Also researched PM3 (in trans with a
  pathogenic variant) and explicitly sized it as a real multi-proband
  points-aggregation system (ClinGen SVI Recommendation for the in trans
  Criterion), a genuinely larger architectural change than this batch's
  other criteria needed — sized and disclosed, not implemented this
  round (see "Batch 28" below for when it was). 15 new tests
  (`test_ps3_bs3_evaluators.py`), 212 total.
- **Batch 26 (PVS1 splice-RNA-evidence branch)** — done. Extended PVS1
  (`evaluators/pvs1.py`) with a new `SplicingRnaEvidence`-driven branch
  for splice donor/acceptor variants, grounded in a real fact found this
  batch: the ClinGen LGMD VCEP's own CAPN3 specification states
  experimental splicing evidence should be scored under PVS1, not PS3,
  per the ClinGen SVI Splicing Subgroup (Walker et al. 2023) — see
  "Splice-RNA evidence feeds PVS1 directly (batch 26)" above for the full
  design, including exactly which two branches are threshold-free and
  implemented (`CONFIRMED_NULL_EQUIVALENT` -> MET Very Strong,
  `CONFIRMED_NORMAL_SPLICING` -> NOT_MET) versus which remain open (the
  full decision tree's percentage/protein-region-criticality thresholds
  — both the primary paper and the CAPN3-specific PVS1 flowchart PDF were
  unreachable this batch). New field `TranscriptConsequence.splicing_rna_evidence`,
  optional and unset by default. Enriched `CAPN3_c.946-1G>A`'s curated
  notes with `INCONCLUSIVE` (a real cited experimental study exists, but
  doesn't state frame/NMD outcome precisely enough) rather than leaving
  it silently unset; the other two real splice fixtures remain unset
  (genuinely no experimental data in their records) — none of the three
  had data precise enough to reach `CONFIRMED_NULL_EQUIVALENT`, so that
  MET path is demonstrated on one new disclosed synthetic fixture,
  `CAPN3_SYNTH_PVS1_SPLICE_RNA_01`, deliberately constructed to avoid
  also becoming a sixth Table-5-vs-Bayesian divergence case. `engine.py`/
  `bayesian.py` needed zero changes. 7 new tests (`test_pvs1_evaluator.py`
  grew from 12 to 19), 219 total.
- **Batch 27 (PVS1 start-loss alternative-start-codon branch)** — done.
  Continued PVS1's scope completion: unlike batch 26, this rule's exact
  text WAS reachable (ACGS 2024 UK Practice Guidelines, quoting Abou
  Tayoun et al. 2018's initiation-codon decision tree verbatim) — see
  "PVS1 start-loss: the alternative-start-codon rule (batch 27)" above
  for the full design and exact quotes. Three new `TranscriptConsequence`
  fields (`alternative_start_codon_identified`,
  `alternative_start_codon_percent_protein_lost`,
  `alternative_start_codon_preceded_by_pathogenic_variant`), and a new
  four-way `evaluate_pvs1()` start-loss branch: no alternative found ->
  MET Very Strong; a clean alternative (<=10% of protein, no pathogenic
  variant in the lost region) -> MET Supporting; a disqualified
  alternative (>10%, or a pathogenic variant present) -> MANUAL_REVIEW
  with a rationale naming exactly why; unassessed -> unchanged
  MANUAL_REVIEW default. `CAPN3_c.1A>G`, previously an unresolved
  generic-MANUAL_REVIEW real fixture, is now independently and fully
  resolved: its downstream alternative start codon (found by fetching
  NM_000070.3's real CDS sequence directly from Ensembl's REST API and
  scanning it programmatically) sits at codon 228, losing 27.6% of the
  protein -- doubly disqualified, since this project's own
  `CAPN3_c.550del` fixture (a real, ClinVar-Pathogenic variant at residue
  184) independently confirms the lost region isn't dispensable. Two new
  disclosed synthetic fixtures (`CAPN3_SYNTH_PVS1_STARTLOSS_NO_ALT_01`,
  `CAPN3_SYNTH_PVS1_STARTLOSS_SUPPORTING_01`) demonstrate the two MET
  branches no real fixture could reach. `engine.py`/`bayesian.py` again
  needed zero changes. 8 new tests (`test_pvs1_evaluator.py` grew from
  19 to 27), 227 total.
- **Batch 28 (PM3 implemented)** — done. Implemented PM3
  (`evaluators/pm3.py`, new `Pm3Evidence`/`Pm3ProbandObservation` models)
  — see "PM3: implemented (batch 28)" above for the full design. The
  primary ClinGen SVI points-table PDF was unreachable again this batch
  (same failure mode as batch 26's Walker et al. 2023 sources), so rather
  than block a second time or hardcode a recalled-but-unverified numeric
  table, per-proband points became a curated fact (like
  `FunctionalEvidence.validation_strength`), while the two rules that
  WERE confirmed and quoted (the homozygous 1-point cap, the
  cis-cooccurrence override) are enforced directly. This project's
  twelfth evaluator, and its first Strong-strength pathogenic-direction
  one. Enriched `CAPN3_c.550del` with two independent real published-
  cohort homozygous PM3 observations (Czech, Polish) — the fixture
  reaches PM3 Strong, and combined with its existing PVS1 Very Strong,
  now classifies PATHOGENIC under Table 5 (previously VUS since
  Milestone 3), finally matching its real ClinVar call and dropping out
  of the Table-5-vs-Bayesian divergence set. Two new synthetic fixtures
  (`CAPN3_SYNTH_PM3_MODERATE_01`, `CAPN3_SYNTH_PM3_CIS_OVERRIDE_01`)
  cover the compound-heterozygous phase branch and the cis-cooccurrence
  override, since no real curated fixture has either on record.
  `clinical.py` unchanged — remains complementary, not overlapping. 18
  new tests (`test_pm3_evaluator.py`), 245 total.
- **Batch 29 (X-linked female/biallelic-XX case interpretation)** — done.
  Extended `clinical.py`'s `interpret_x_linked_case` beyond the
  hemizygous-male-only case — see "X-linked female/other-karyotype case
  interpretation (batch 29)" above for the full design. Batch 12's
  finding that X-inactivation cannot reliably predict a single
  heterozygous carrier's phenotype (Brioschi et al. 2012) still holds and
  is unchanged; batch 29 instead found and implemented a second, real,
  genuinely X-inactivation-*independent* mechanism this project's
  existing evidence model could newly represent: biallelic XX
  involvement (Ulm et al. 2022; Fujii et al. 2009; Takeshita et al.
  2017). `karyotypic_sex=XX` now accepts one OR two variant_ids;
  confirmed-trans plus both variants independently Pathogenic/Likely
  Pathogenic resolves to EXPLAINED, while every other two-variant XX
  combination (notably including cis, deliberately NOT handled the way
  autosomal recessive cis is — see the design note for the real
  X-inactivation-mosaicism reasoning) stays MANUAL_REVIEW alongside the
  unchanged single-variant case. `karyotypic_sex=XY` now explicitly
  rejects two variant_ids (a hemizygous individual only has one X) rather
  than silently mishandling them, and `OTHER`/`UNKNOWN` remain exactly as
  narrow as before — genuinely heterogeneous karyotypes (X0/Turner
  functionally hemizygous, XXY diploid-X, mosaicism neither) this project
  still does not attempt to disambiguate. Two new fixtures
  (`CASE_DMD_XX_BIALLELIC_TRANS`/`_CIS`), synthetic at the case level but
  built from two real, already-curated DMD variants, since no real
  published case pairs two point mutations in one biallelic-XX patient
  (real cases found are all CNV-based, which doesn't fit this project's
  point-mutation-only `ClinicalCase`). 5 new tests (`test_clinical.py`),
  250 total.
- **Batch 30 (project wrap-up)** — done. Documentation-only, no code
  changes: rewrote this section's closing bullet into the "Final status
  and scope boundary" section below, updated the top-of-file Status line
  to mark the project COMPLETE, and added this entry so the batch-by-
  batch history stays complete through the actual last batch. The
  decision to close here rather than keep going was the user's, made
  explicitly at the start of this batch rather than assumed.

## Final status and scope boundary

This project is complete, not paused. Batch 29 (X-linked female/biallelic
case interpretation) was the last batch of new functionality; batch 30
closed out the documentation to reflect that rather than leaving the
Roadmap reading like an open todo list. Nothing below is a promise of
future work — it's a disclosed inventory of exactly where the scope
boundary sits, in the same spirit as every other "researched, sized, and
deliberately not implemented" decision this project made along the way,
so a future reader (in this repo or a similar one) knows precisely what
was checked versus what was left honestly unresolved.

What's here, end to end: two genes (CAPN3 autosomal recessive, DMD
X-linked), twelve ACMG/AMP point-mutation criterion evaluators (PVS1,
PM2, PM4, PS1, PM5, PM3, PS3, PP3, BP4, BA1, BS1, BS3), two independent,
fully tested combining systems (Table 5 and Bayesian point-based),
DMD CNV deletion and duplication scoring as a separate parallel system,
case-level interpretation for autosomal recessive and X-linked
(hemizygous male and biallelic-XX) inheritance, a real pipeline-adapter
integration, 29 point-mutation fixtures (19 real ClinVar/VCEP-grounded)
plus 5 CNV fixtures plus 11 case-level fixtures, all checked against
golden cases curated independently of the code that's judged against
them, and 250 tests passing under both a pytest suite and a
dependency-free runner.

What's deliberately out of scope, and why, grouped by area:

- **Curated fixture breadth.** 19 real variants toward the original
  ~20-30 target. Every criterion-level real-data gap that was searched
  for is now closed except five specific branch examples (a real PS1
  precedent pair, a real BS3-MET example, a real PVS1
  CONFIRMED_NULL_EQUIVALENT example, real PVS1 start-loss "no
  alternative"/"clean automatic-downgrade" examples, and a real
  compound-heterozygous-confirmed-trans PM3 example) — each searched for
  and not found, not simply unconsidered. Further additions beyond these
  five would mostly add volume, not new logic.
- **PVS1's full decision tree.** The real protein-domain-criticality
  judgment inside the NMD-escape and start-loss MANUAL_REVIEW branches,
  and Walker et al. 2023's exact percentage-of-transcript splicing
  thresholds, remain open. Both require primary sources (the Walker
  paper itself, the CAPN3-specific PVS1 flowchart PDF) that were
  unreachable via this project's tooling on every attempt across three
  separate batches (26, 27, 28's PM3 research) — a consistent, disclosed
  tooling limitation, not a research gap.
- **CNV scoring beyond deletions/duplications.** Whole-gene
  triplosensitivity scoring for genes where it's real (confirmed not
  DMD, per batch 24), and Riggs et al. 2020 Sections 1 (genomic
  content), 2H (predicted-but-not-established HI), 3 (gene count), 4
  (case/case-control/population evidence), and 5 (inheritance/family
  history) are all unimplemented. The primary paper's exact loss-side
  2B/2G definitions and gain-side letter codes/point values also rely on
  ClassifyCNV's reimplementation and a disclosed secondary-source
  inference (batch 24) rather than independent primary verification.
- **X-linked karyotypes beyond XY and biallelic XX.** `karyotypic_sex=OTHER`
  remains a single deferred bucket (X0/Turner, XXY, and mosaicism all
  return MANUAL_REVIEW without being distinguished from each other),
  even though batch 29's research found that X0/Turner is functionally
  hemizygous (the XY-shaped logic would apply) and XXY is diploid-X (the
  XX-shaped logic would apply) — splitting `OTHER` into specific real
  karyotypes was identified but not built, since `KaryotypicSex.OTHER`
  itself would need to become several distinct values first, a schema
  change this project stopped short of.
- **A real point-mutation-only biallelic-XX case.** No real published
  case pairs two point mutations (rather than CNVs) in one biallelic-XX
  patient, so `CASE_DMD_XX_BIALLELIC_TRANS`/`_CIS` are synthetic cases
  built from two real, already-curated variants rather than an
  end-to-end real fixture.
- **Table 5 vs. Bayesian default.** `clinical.py` (and any other caller)
  takes a pre-computed `classifications` dict and is agnostic to which
  combining system produced it — Milestone 5 built both as equally real,
  equally tested options and deliberately left the choice of which one a
  caller should default to unmade, rather than picking a "winner" this
  project has no authority to declare.

None of the above is silently missing — each is named here, in the
design note it originated from, and in the code (a docstring, a
rationale string, or both) that a reader would actually encounter while
using the affected feature.
