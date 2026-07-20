# data/synthetic/

Placeholder for larger generated/synthetic test datasets (e.g. randomly
constructed variants for stress-testing the loader or combining engine)
that are too voluminous or too clearly non-biological to belong in
`data/curated/` alongside the small, individually-reviewed fixture set.

Empty in Milestone 1: the two synthetic cases built so far
(`CAPN3_SYNTH_PATHOGENIC_01`, `CAPN3_SYNTH_LIKELY_BENIGN_01`) are small and
hand-reviewed enough to live directly in `data/curated/variant_evidence.json`
with the real case, each clearly labelled as synthetic in its `notes` field.
This directory is here for when that stops being true.
