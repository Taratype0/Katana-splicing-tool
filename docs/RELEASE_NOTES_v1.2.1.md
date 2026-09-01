# Katana CoSplice v1.2.1

## Confidence optimization

- Kept the frozen posterior-mean WGCNA/TOM partition and rejected further module-boundary tuning.
- Added fixed-module bootstrap kME as a complementary membership-stability signal.
- Added split-half confirmation: `core` and `peripheral` confidence must reproduce in two disjoint bootstrap halves.
- Added event-level evidence provenance, both split-half roles, kME lower-tail stability and module-level kME summaries.
- Added a membership-confidence plot and evidence summary table to the real-count ablation.
- Fixed direct execution of `scripts/run_real_count_four_arm_ablation.py`.

## Frozen Setd2 result

Using the same 3,000 Stable-AS events, 36 samples, seed 37 and 10 condition-stratified posterior-draw bootstraps:

- the production partition remains six modules, 1,624 assigned events and a largest-module fraction of 0.327;
- posterior-mean sample-bootstrap ARI remains 0.559 versus 0.546 for raw-PSI WGCNA;
- split-half confirmed confidence contains 844 core, 422 peripheral and 358 uncertain assigned events;
- core membership therefore rises from 138/1,624 (8.5%) under co-assignment-only annotation to 844/1,624 (52.0%) under the dual-evidence, split-half-confirmed model;
- 716 core events are supported by bootstrap kME, 109 by both kME and co-assignment, and 19 by co-assignment.

The improvement is in stable event-membership identification, not in the WGCNA partition ARI. CoSplice does not claim universal clustering superiority over WGCNA.
