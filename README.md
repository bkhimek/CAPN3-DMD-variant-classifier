# CAPN3/DMD Variant Classifier

A small, scoped prototype implementing the ACMG Engine architecture
described in the companion design-guide set, starting with two genes:
**CAPN3** (autosomal recessive, LGMDR1/calpainopathy) and **DMD**
(X-linked, out of schema scope until Milestone 4 — see Roadmap).

## Status: Milestone 4 complete, curated variant set at 21 of ~20-30

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
- **Twenty-one curated evidence bundles** (`data/curated/variant_evidence.json`):
  seventeen real ClinVar-grounded variants (eleven CAPN3, six DMD) and four
  synthetic cases (three CAPN3, one DMD) constructed to exercise specific
  combining-rule and case-level paths. This is seventeen increments toward
  the ~20-30 ClinVar variant set from the original project plan — see
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
- Seven evaluators, one per implemented criterion
  (`src/variant_classifier/evaluators/`): `pvs1.py`, `pm2.py`, `pm4.py`,
  `pp3.py`, `bp4.py`, `ba1.py`, `bs1.py`. PVS1 is deliberately partial —
  see "PVS1 scope" below. PM4 was added in batch 14 — see "PM4: a second
  new criterion" below.
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
- Twenty-one curated variants (CAPN3 and DMD) and nine curated
  `ClinicalCase` fixtures covering every branch above — including, as of
  batch 12, a pair built on a real ClinVar-sourced DMD variant rather
  than only the original synthetic ones, and, as of batch 17, the case-
  level layer's previously-untested MANUAL_REVIEW catch-all branch (a
  hemizygous male whose only variant is VUS, not yet qualifying or
  benign). Every evaluator, the combining
  engine, and the case interpretation layer are each verified against
  golden cases written independently of the code — including, for the
  case-level tests, that trans vs cis and male vs female change the
  outcome despite everything else being identical
  (`tests/unit/test_*.py`). 118 tests pass in total (the "matches golden
  case" tests iterate over however many fixtures exist rather than a
  hardcoded count, so this number grows automatically as the curated set
  does; also includes two hand-built tests pairing the real DMD
  male/female cases (batch 12), proving no real CAPN3 variant can
  currently reach Pathogenic/Likely Pathogenic (batch 13, extended in
  batch 14 to cover PM4 too), and a dedicated `test_pm4_evaluator.py`
  suite added in batch 14).

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
combining rather than Table 5 (see Roadmap).

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

One real-data gap remains open: a real, precisely-sourced BP4-MET example
(a calibrated computational score landing BENIGN from real data). The
project's one real benign-leaning computational score remains
`CAPN3_c.2257G>A`'s INDETERMINATE REVEL (batch 8); DMD's real Pejaver
calibration (batch 16) opens a second, easier-to-clear avenue for this
(BP4 Supporting at REVEL <=0.290, vs. CAPN3's real VCEP threshold of
<=0.10) that hasn't been searched yet.

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
once the evaluator needed it for stop-gained variants too. All three
documented scope gaps (splice donor/acceptor, start-loss, NMD-escaping
truncations) now have at least one real fixture exercising them, as of
batch 9's `DMD_c.11041A>T` — see "Expanding the curated set" below.

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
    pm4.py                    evaluate_pm4() — see "PM4: a second new criterion" above
    pp3.py                    evaluate_pp3()
    bp4.py                    evaluate_bp4()
    ba1.py                    evaluate_ba1()
    bs1.py                    evaluate_bs1()

config/
  population_thresholds.yaml per-gene PM2/BA1/BS1 frequency thresholds -- CAPN3's
                              are the real ClinGen LGMD VCEP values as of batch 4
                              (see Design notes); DMD's remain generic ACMG defaults

data/
  curated/
    gene_disease_context.yaml   CAPN3 and DMD
    variant_evidence.json       21 curated variants (CAPN3 + DMD) -- growing toward ~20-30
    clinical_cases.json         9 curated ClinicalCase fixtures (Milestone 4)
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
extra environment variables. All 118 tests currently pass.

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
- Later: continue expanding curated fixtures toward the full 20-30
  ClinVar variant set (21 of ~20-30 done, see "Expanding the curated set"
  above; a real BP4-MET example is the one known remaining gap); add
  PM3/PS1/PM5/PS3/BS3
  as real per-variant criteria (distinct from how clinical.py currently
  handles trans/cis reasoning at the case level);
  revisit PVS1's partial scope (protein-domain criticality,
  constitutive-exon data — the real CAPN3 VCEP spec includes a gene-specific
  PVS1 flowchart that isn't implemented here either); revisit DMD's CNV
  representation gap (see gene_disease_context.yaml) if real DMD variants
  are ever added; extend X-linked case interpretation beyond the
  hemizygous-male case once X-inactivation is worth modeling properly;
  **consider switching `engine.py`'s combining rules from classic Table 5
  to the VCEP's Bayesian point-based system** — batch 4 produced a
  concrete, real example (`CAPN3_c.1939G>T`) of Table 5 under-calling a
  variant the point system would resolve, which is a much stronger
  motivating case than existed before batch 4.
