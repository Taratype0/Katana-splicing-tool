# CoSplice and MotifTrace integration review

Review date: 2026-11-02. Status: **not ready for main promotion**.

This review compares the selected local sealed application, `katana-cosplice-event-options-fix-packaged-20261102`, with the pipeline development checkout and GitHub main `4f046758cefe244082c2acf4436c3215d7ebb6b3`. The pipeline entry files listed below are identical between that main commit and the reviewed development baseline. No production implementation was replaced in this review.

## What is connected

| Entry | Selected local application | Published pipeline |
| --- | --- | --- |
| MotifTrace | EventGraph 1.1.0; workflow `prioritize.motiftrace_consequences` 2.0.0; standalone CLI and workflow dispatcher exist | Workflow still imports the legacy execution adapter |
| CoSplice | `programs.cosplice` and `programs.psi_network` resolve to the whole-study workflow 4.3.3 | Native network task adapter still invokes legacy `CoSplicePipeline` |
| Dedicated network command | `katana-cosplice` itself still invokes the legacy implementation, even in the selected local application | Not a verified standalone entry for the new engine |
| Report | Latest renderer produces the seven-chapter synthetic report, motif sequence cards and whole-study network browser | Updated computation, portable installation and report collection have not been jointly verified |

The native Nextflow DAG already contains network, motif and report processes. Their existence does not establish that they execute the newly selected algorithms. An ordinary `nextflow run` of the current published version is **not** an accepted execution route for the new CoSplice/MotifTrace pair.

## Integration work still required

1. Route both standalone and Nextflow entries to the same versioned implementation. Preserve explicit legacy compatibility rather than silently routing names to different algorithms.
2. Package the new `cosplice_core`, its annotation resources and native dependencies. It currently sits outside the Python application package discovery.
3. Replace machine-specific runtime discovery: the network core fixes a Windows R executable path; the formal local report configuration fixes an external Python environment and sibling release.
4. Resolve the new network input contract before claiming support for a six-sample study. The actual `prepare` function rejects six generated independent units with `Fewer than 15 independent biological units`. This is a whole-study minimum, not six samples per comparison. Memory/time safeguards also remain in this implementation. This review does not change these scientific/resource rules.
5. Verify the published wheel and matching Linux execution environment through the new native network and motif routes, collection and report commands. A local report render is insufficient evidence.

The selected new network version builds **one whole-study network** and displays its module associations across conditions. Changing a comparison in the report does not reconstruct the network. This is the current sealed implementation, distinct from the historical per-comparison pipeline route.

## Checks performed

- Verified all 926 files in the selected application's source hash manifest; no mismatch.
- Loaded both workflow registries in isolated subprocesses and recorded their real implementation modules, versions and defaults.
- Inspected the Nextflow DAG, planner, native task adapter, standalone CLI and packaging declarations.
- Exercised the current network preflight with generated data; confirmed the six-unit rejection without running a network fit.
- Ran 12 existing EventGraph/frozen-method checks: 11 passed; one renderer test failed because it asserts an older exact sentence. Its preceding assertions preserving regional significance and zero significant windows passed. The failed test is retained, not reported as a pass.
- Verified the current synthetic report's regional-FWER and window-BH distinction, 10/5/0 qualifying motif-card cases, ordinary-DAS/RS-DAS browser switching, embedded network and complete tables.
- Checked all 77 displayed figures have data links, 1,365 local links resolve inside the package, and all 413 ZIP members equal their on-disk files.

These are source, generated-data and presentation checks. There was no real-data rerun, Slurm job, matching container acceptance or full new native scientific execution.

## Synthetic report reproduction

The review scripts are `scripts/audit_updated_module_entries.py`, `scripts/build_integration_review_dummy.py` and `scripts/check_integration_review_dummy.py`.

For the report builder, use Python with the selected release's `app` directory on `PYTHONPATH`, and supply `--outdir NEW_DIRECTORY --network-demo SYNTHETIC_NETWORK_PREVIEW_DIRECTORY`. The script intentionally refuses an existing output directory.

All panels and statistical values are synthetic display fixtures; no biological interpretation is implied. MAPP and matched CLIP are optional calculations not illustrated in this dummy. The unmodified renderer page is retained outside the delivery directory as `renderer-original.html`. In the review copy only, the header is explicitly synthetic and missing-real-input notices are grouped in the appendix; their original records remain in the report manifest and `DUMMY_REVIEW.json`. No production renderer was edited.

Approval of this layout does not replace the remaining integration and installation checks. Main promotion remains subject to user approval and passing those checks.
