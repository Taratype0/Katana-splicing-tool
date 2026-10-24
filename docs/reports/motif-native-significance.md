# Motif significance in the report

The report shows at most ten distinct motif models with a recorded native test
FDR at or below 0.05. It selects one actual window or full-sequence test per
model, ordered by its native FDR, its paired native P value, and stable model
identity. Both values and the effect come from the same recorded test row.
The displayed FDR belongs to the complete comparison test family; it is not a
new whole-model FDR or a test of the illustrated RNA sequence.

When no model has an eligible significant test, the comparison displays a
short **no significant motif** notice. Sequence panels and diagnostic plots
do not fill the report body. If tests cannot be estimated, the notice says
**not estimable** instead of treating missing values as a negative finding.
All tested models, source P/FDR values, sequences and complete output tables
remain available in the analysis outputs.

Identical adjusted values can be legitimate. Native paired-swap P values are
discrete, and Benjamini–Hochberg correction across the complete set of finite
window and full-sequence tests can assign identical FDR values to distinct
tests. Report selection does not recalculate these values or correct only the
displayed Top 10.

MAPP is a separate result family. Its original per-fit/window correction is
labelled by its native method; Bonferroni-adjusted values are not relabelled as
MotifTrace FDR. RBP identity, boundary motif enrichment, CLIP evidence and local
RNA sequence consequences retain their separate statistical scopes.
