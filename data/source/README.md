# data/source/

Placeholder for raw, unmodified pulls from external sources (ClinVar XML/VCF
exports, gnomAD API responses, VEP annotation output, ClinGen gene-disease
validity exports, etc.) — kept separate from `data/curated/` so it is always
possible to tell what a human curator added, changed, or normalised versus
what a source system originally returned.

Empty in Milestone 1: the one real case in this milestone (`CAPN3_c.550del`)
was curated by hand from published ClinVar/gnomAD information rather than
from a saved raw API/export file. Milestone 2+ (batch ClinVar/gnomAD/VEP
retrieval) is expected to populate this directory.
