# Katana v1.0.2 GUI and Service Hardening Acceptance

Date: 2026-07-28

## Scope

This acceptance pass covers desktop behavior across displays, plot interaction,
responsive navigation, detail/data dialogs, page-open execution behavior, and a
non-blocking incremental split of `ProjectService`.

Representative project:

- `E:\Project\RNASEQ\Result\Setd2_results`
- 24 configured comparisons

## Manual GUI Acceptance

Tested window configurations:

- Laptop display, approximately `1363 x 891`
- External display, approximately `1360 x 890`
- Maximized external display, `1920 x 1032`
- Window movement between displays with different origins and scaling

Verified:

- Project loading completed without starting analysis modules.
- The sidebar and page content remained scrollable on the smaller display.
- The layout remained usable after moving and maximizing the application across
  displays.
- 5.2 opened without triggering 5.3 or another downstream analysis.
- Event-type navigation exposed SE, RI, A3SS, A5SS, and MXE.
- 5.4 and 5.4.1 opened as lightweight state/cache views and did not run analysis.
- Plot wheel zoom and mouse pan worked without first changing width or height.
- Plot canvases relaid out after tab/window size changes.
- Data preview dialogs could be maximized.
- Detail dialogs were scrollable and now expose minimize/maximize controls.
- Long plot titles reserve additional top margin and wrap before exceeding the
  preview width.

The last visual spot-check after reducing the shared title wrap width to 46
characters was interrupted by the user with Escape. Static verification confirms
the representative titles wrap to lines of 40-42 characters, but this final
single visual confirmation remains outstanding.

## GUI Fixes

- Shared plot layout now reserves more title space (`top=0.78`).
- Shared plot titles wrap at 46 characters.
- Comparison detail dialogs use normal top-level window flags with minimize and
  maximize buttons.

These changes update display layout only. They do not trigger analysis,
generation, cache rebuilding, or downstream execution.

## Service Hardening

The following responsibilities were extracted from `ProjectService` without
changing their public service contract:

- Replicate-aware cross-comparison interaction logic:
  `src/services/cross_interaction_mixin.py`
- Native Jutils plot generation:
  `src/services/native_jutils_plot_mixin.py`
- Cross-comparison directional pattern table construction:
  `src/services/cross_comparison_service.py`

`ProjectService` is approximately 1,800 lines smaller. It remains a large
orchestration service and can be split further incrementally, but this does not
block the current release.

## Automated Verification

Commands:

```text
python -m compileall -q app src katana_cosplice tests
git diff --check
pytest -q
```

Result:

- Compilation passed.
- Diff validation passed.
- `85 passed, 5 warnings`.
- Warnings are limited to existing statsmodels mixed-model covariance/convergence
  diagnostics and one seaborn pending deprecation warning.

## Release Assessment

The analyzed code paths and representative desktop interactions are release
ready. The only remaining manual item is the final visual spot-check of the
46-character wrapped title on both displays; no analysis or data correctness
risk is associated with that check.
