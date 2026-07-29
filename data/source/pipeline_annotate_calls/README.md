# pipeline_annotate_calls/

`adapter_test_sample.vcf.gz` is a real (not synthetic), hand-trimmed slice of
CAPN3-DMD-variant-calling-pipeline's actual
`results/annotate_calls/annotated_calls.vcf.gz` output — the GATK/DeepVariant-
concordant HG002 call set, annotated with VEP (transcript consequence) and
gnomAD v4.1 (population frequency). Extracted 2026-07-29 via
`bcftools view -r <5 positions>` from a real, fully-verified pipeline run
(see CAPN3-DMD-variant-calling-pipeline's README.md "Status" section and this workspace's
`~/projects/HANDOFF.md`) — not fabricated.

Five records, chosen to exercise every branch of `pipeline_adapter.py`:

| Position | REF/ALT | Why it's here |
|----------|---------|----------------|
| chr15:42259883 | C/T | Falls in CAPN3's padded flanking region but outside the MANE transcript's own span (no CSQ hit on either gene's MANE transcript) — tests the "no bundle, by design" skip path. |
| chr15:42389001 | G/A | Real CAPN3 missense (`p.Ala236Thr`) on the MANE transcript, gnomAD-observed common variant (AF=0.222) — the main OBSERVED-in-gnomAD path, hand-verified independently during this adapter's build (see CAPN3-DMD-variant-calling-pipeline's README/HANDOFF). |
| chr15:42398590 | TACACACACACACAC/T,TACACAC | Multi-allelic, both alleles land in CAPN3's intron — tests per-allele CSQ/gnomAD splitting (one allele OBSERVED, one ABSENT in gnomAD) and the VEP indel-allele-trimming logic (`_vep_allele`). |
| chrX:31478233 | C/T | Real DMD missense (`p.Arg2937Gln`) on the MANE transcript, gnomAD-observed near-fixed common variant (AF=0.939) — same OBSERVED path as the CAPN3 case, other gene/chromosome. |
| chrX:32456578 | TTGTG/T,TTG | Multi-allelic DMD intron site — second multi-allelic exercise, confirms the splitting logic isn't CAPN3-specific. |

HG002 is confirmed clinically empty for CAPN3/DMD pathogenic variation (see
`~/projects/HANDOFF.md`) — every variant here is ordinary benign background
polymorphism. This fixture is for adapter/schema testing only, never cited as
evidence of the classifier detecting a real pathogenic finding.
