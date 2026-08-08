# How to Analyze Any Other Variant — A Practical Guide

*Written after Batch 31 (adding BRCA1), which was deliberately treated as a stress test of this question: does the classifier's architecture actually generalize beyond the two gene/inheritance patterns it was built against? The answer was yes, but not for free — PM2 needed a real, gene-gated code fix along the way, and separately, BA1/BS1's founder-frequency handling required curating the right fixture data against evaluator logic that was already correct and never changed — two different categories of "not for free," worth telling apart. This guide turns that experience into a repeatable procedure.*

---

## First, decide which question you're actually asking

There are two different jobs that both sound like "analyze a variant in a new gene," and they take very different amounts of work:

**A. "Is this one variant, in a gene we don't currently support, pathogenic or benign?"**
You need steps 1–4 below, at minimum, before the existing 12 evaluators can be trusted to produce a meaningful answer. Skipping straight to "just run it through the code" will produce a number, but not a defensible one — the evaluators are gene-agnostic in their *logic*, not in the *thresholds and applicability rules* they depend on, and those come from the new gene's real specification.

**B. "Add durable support for this gene, so future variants in it can be classified routinely."**
You need the full seven-step procedure below, ending in curated fixtures, hand-derived expected results, and guard tests — the same discipline every one of the project's 31 batches has followed.

If you're not sure which one you're doing, do (A) first. It's a strict subset of (B) and tells you quickly whether the gene is a clean fit for the existing architecture or something that will need new case-level logic (as BRCA1 did).

---

## Steps 1–4: the real prerequisite (needed even for a single variant)

### 1. Confirm the gene's disease mechanism and inheritance pattern against a real, citable source

Not a generic reading of the base ACMG/AMP guidelines — those are gene-agnostic by design and don't tell you the numbers that matter for *this* gene. Look for a published expert-panel specification first:

- A **ClinGen Variant Curation Expert Panel (VCEP)** specification for the gene, if one exists. This is what BRCA1 used (ENIGMA BRCA1/BRCA2 VCEP, v1.2). VCEPs override the generic ACMG/AMP thresholds with gene-specific ones — different `PM2` allele-frequency cutoffs, different `PVS1` applicability rules, sometimes entirely new criteria (PM3/BS2's points-based Fanconi-anemia logic for BRCA1, currently out of scope here).
- If no VCEP exists, the next-best sources are the original disease-gene literature and OMIM/GeneReviews for inheritance pattern, with any quantitative thresholds sourced from population databases (gnomAD) or functional-assay literature directly, and disclosed as a secondary source if the primary one is unreachable — the same pattern this project used three times for CAPN3/DMD (CNV rubric, PVS1 splice thresholds, PM3 points table).

Record the answer to one specific question before moving on: **does a single damaging allele cause disease (dominant), does it require two (recessive), or is the risk elevated but not deterministic (incomplete penetrance, as with BRCA1)?** This determines almost everything else.

### 2. Audit every existing evaluator against that real specification — one criterion at a time

Don't assume the shared evaluator code is safe for a new gene just because it was gene-agnostic for the genes tried so far. Go through all twelve, criterion by criterion, and ask two questions for each: does this criterion even apply to this gene's variant types and mechanism, and if so, does the existing generic logic already implement it correctly, or does the new gene's spec define a different rule?

| Criterion | What to check for a new gene |
|---|---|
| PVS1 | Does loss-of-function actually cause this gene's disease? (Not all genes are LOF-mechanism genes.) |
| PM2 / BA1 / BS1 | What allele-frequency thresholds does the gene's own spec define? (This is exactly where BRCA1 needed a fix — see the worked example below.) |
| PM4 | Does the gene have a known structurally important region where in-frame indels matter? |
| PS1 / PM5 | Same-residue precedent logic is generally mechanism-agnostic, but confirm the gene's spec doesn't add extra conditions. |
| PM3 / BS2 | Only applies to recessive/biallelic-relevant genes with a real points table or trans/cis logic — BRCA1's version of this is a structurally different, points-based Fanconi-anemia code and was explicitly deferred rather than force-fit. |
| PS3 / BS3 | Does a validated functional assay exist for this gene? (BRCA1 used the Findlay et al. 2018 saturation genome editing dataset.) |
| PP3 / BP4 | Are the computational predictors already in use validated for this gene, or does the spec recommend different ones/different thresholds? |

The BRCA1 audit found eleven of twelve needed zero changes. Expect a similarly short list of real changes for most genes — but confirm it by checking, not by assuming the last gene's result generalizes.

### 3. Add the gene's real, cited thresholds to configuration, not code, wherever possible

Follow the pattern `PM2`/`BA1`/`BS1` already use in `population_thresholds.yaml`: gene-specific numeric thresholds and boolean flags live in config, keyed by gene, so a shared evaluator function stays gene-agnostic in its *logic* while varying its *behavior* per gene through data it's handed. This is what made the actual BRCA1 fix minimal — one new opt-in config key (`pm2_excludes_indel_delins: true`), not a rewritten evaluator.

### 4. Decide whether the gene's inheritance pattern needs a new case-level outcome

If the new gene fits an inheritance pattern the classifier already models (autosomal recessive like CAPN3, X-linked like DMD), the existing case-level statuses (`EXPLAINED` / does-not-explain / flag-for-review) likely already fit. If it doesn't — as BRCA1's incomplete penetrance didn't — don't overload an existing status to mean something subtly different. Add a new, clearly distinct one. BRCA1 got `RISK_CONFERRING`: "elevated, but not deterministic, disease risk," kept separate from the deterministic `EXPLAINED` outcome the other two genes use, because collapsing the two would have been a real misrepresentation, not a simplification.

**Steps 1–4 are where the actual thinking happens.** Once they're done, running a single variant through the existing evaluators is genuinely mechanical. Skipping them and running the variant anyway will produce output that looks like a real classification but isn't grounded in anything gene-specific — the single most likely way to get a confidently wrong answer out of this system.

---

## Steps 5–7: turning that into durable, tested support for the gene

Only needed if you're adding the gene properly (path B above), not just answering one variant's question.

### 5. Curate a small set of real, cited example variants

Prioritize real, publicly sourced examples — ClinVar accessions, gnomAD frequency data, published functional-assay datasets — over invented ones, and cover the gene's most clinically important evidence types (a clear pathogenic LOF case, a clear benign/common case, at least one case that exercises whatever the gene's distinctive logic is). Add synthetic fixtures only where no real example demonstrates a specific logic path you need to test, and always label them as synthetic in the data itself — never let a synthetic fixture pass as if it were real. BRCA1 landed on 4 real + a handful of synthetic point-mutation fixtures, plus 3 real-mechanism case-level genotype scenarios (monoallelic pathogenic / monoallelic benign / VUS).

### 6. Write the expected result for every example independently, by hand, before running any code against it

This is the discipline that actually catches errors — a test suite that's curated from what the code already outputs will certify bugs along with everything else. Hand-derive each fixture's expected classification from the gene's real spec and the evidence you cited, write it down, *then* run the evaluator and compare. A mismatch is either a bug in the code or an error in your hand-derivation — both are worth finding before this ships.

### 7. Write a guard test for any known limitation the design deliberately routes around rather than fixes

If you decide a criterion doesn't need a code change for this gene because none of your fixtures happen to exercise it (the way BRCA1 initially assumed for PM4, PM5, and PM2-on-indels), don't leave that as an implicit assumption — write an automated test that checks it explicitly, so a future contributor adding one more fixture to this gene can't silently reintroduce the gap. This is exactly the test that, for BRCA1, surfaced that the PM4/PM5 avoidance-by-fixture-shape trick worked cleanly but the PM2 one didn't — turning a would-be silent bug into a caught, fixed, and disclosed one before release.

---

## Worked example 1: what this looked like for BRCA1's PM2

Worth walking through once end-to-end, because it's the one place the "eleven of twelve need nothing" story wasn't quite true, and shows what step 7 catching a real problem actually looks like in practice:

1. **Step 1** established that BRCA1's real spec (ENIGMA VCEP) says PM2 simply does not apply to insertion/deletion variants — a rule the base ACMG/AMP guidelines don't state, because it's gene-specific.
2. **Step 2**'s audit noted this, but the original plan assumed it could be handled by *fixture choice* alone — just don't curate a BRCA1 indel fixture that needs PM2 to fire — rather than a code change, since the shared PM2 evaluator has no concept of variant type at all.
3. That assumption held for two other criteria in the same batch (PM4, PM5) but not for PM2: two of the batch's own real, cited example variants were Ashkenazi-founder frameshift indels that *needed* PM2's underlying population-frequency data anyway, for their own BA1/BS1 founder-frequency handling. So PM2 always ran a real frequency comparison for those fixtures, contradicting the spec's "does not apply" rule.
4. **Step 7**'s guard test — written to lock in the fixture-shape-avoidance assumption for PM4/PM5/PM2 — is what surfaced this. It failed for PM2, which is exactly what a guard test is for.
5. The fix followed **step 3**'s pattern exactly: one small, additive, opt-in config flag (`pm2_excludes_indel_delins`) on the shared PM2 evaluator, checked before any other logic, `true` only for BRCA1, leaving CAPN3/DMD provably unaffected.
6. Both affected fixtures' golden-case results and the project's own documentation were corrected to match — not left as a known-wrong intermediate result just because it didn't change the final classification tier.

---

## Worked example 2: BA1/BS1 and the founder-population trap

A second, differently-shaped lesson from the same batch — worth telling apart from PM2's, because it looks similar on the surface but the fix was a completely different kind:

1. Two of BRCA1's real fixtures are Ashkenazi-founder pathogenic frameshifts (185delAG, 5382insC) — variants that are locally common in one specific population but real, ENIGMA-reviewed Pathogenic calls, not benign polymorphisms.
2. `ba1.py`/`bs1.py` — the "too common to be pathogenic" evaluators — compare a variant's `overall_af` against the gene's threshold *before* ever consulting `ancestry_specific_max_af`. That ordering is correct, gene-agnostic logic, and it was never changed for BRCA1.
3. What was wrong, initially, was the *data*: curating `overall_af` as gnomAD's raw, undifferentiated total allele frequency let a founder-population spike in one subgroup inflate the overall number past BRCA1's BS1 threshold — which would have made BS1 wrongly return MET for a known-Pathogenic variant.
4. The fix was to curate `overall_af` itself as gnomAD's own reported **Grpmax Filtering AF** — a statistic gnomAD computes specifically by finding the highest-frequency *adequately-powered* ancestry group, which by gnomAD's own methodology already excludes small, bottlenecked founder groups like Ashkenazi Jewish. No evaluator code and no config value changed — only the fixture data was corrected to use the right upstream statistic.
5. The practical distinction to carry forward: PM2's fix was a real, disclosed **code change** (new gate, new config key). BA1/BS1's fix was a real, disclosed **fixture-curation correction** against evaluator logic that was already right. Both are legitimate "not for free" findings from generalizing to a new gene — but conflating them (as an early draft of this project's own documentation briefly did) overstates what changed in the code and understates the importance of getting fixture data curation right.

---

## Quick checklist

- [ ] Real, citable disease-mechanism and inheritance-pattern source identified (VCEP preferred)
- [ ] All twelve evaluators individually audited against that source
- [ ] Any new thresholds/flags added to config, not hardcoded into evaluator logic
- [ ] New case-level status added if the inheritance pattern genuinely needs one — not reused from an existing one
- [ ] Real, cited example variants curated; synthetic ones clearly labeled
- [ ] Expected results hand-derived and written down before running any code
- [ ] Guard tests written for every "we're avoiding this by fixture choice, not fixing it" assumption
- [ ] Full test suite run and passing before merging

If you're only analyzing one variant (path A), the first four boxes are the ones that matter — but check them for real, on this gene, rather than assuming the last gene's answer still applies.
