# Development milestones

| Stage | Deliverable | Current status |
|---|---|---|
| M01 | T01–T08 scope, source inventory, method boundaries and protected baseline | Complete as planning/baseline work |
| M02 | T09–T16 and T57–T60: data contracts, eligibility, adapters and headless core | T09 contracts implemented; remainder pending |
| M03 | T17–T48: core science, MotifTrace, networks and qualified single-cell/nucleus adapters | Existing local source; generic package migration and acceptance pending |
| M04 | T49–T56 and T61: English reports and interactive product integration | Existing local report/UI source; new integration pending |
| M05 | T62–T64: builds, dependencies, containers and selected Nextflow execution DAG | Pending |
| M06 | T65–T70: clean-install, scientific, performance and route acceptance | Pending for the new distribution |
| M07 | T71–T72: user/developer documentation and release review | Pending |

All 72 tasks are retained in the development plan. Existing scientific source,
importing previous results, executing new calculations and acceptance are tracked
separately. A GitHub commit does not mark any pending task complete.

Implementation proceeds in task-sized commits with explicit scope and checks.
Algorithms/defaults/statistical scope require versioned changes and affected
acceptance; reports and UI have their own presentation versions. Reuse compatible
completed results and write new runs to new locations. Recovery creates a clean
checkout at a known baseline without deleting data or cleaning a dirty worktree.

Community scope/name approval is a separate gate from internal engineering.
This repository is an independent development checkpoint, not an nf-core release.
