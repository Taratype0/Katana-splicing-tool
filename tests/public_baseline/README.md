# Usable application acceptance

Run from the repository root after installing the application and test dependencies:

```bash
python -m pytest -q tests/public_baseline
```

The original workflow inventory remains preserved in `workflow_inventory.json`.
`workflow_inventory_read_support_v1.json` records the approved every-sample
ten-count bulk policy. Tests permit method-version changes only to the approved
input-policy workflows: Quality, bulk Landscape, DAS, exact cross-comparison and
count-bearing CoSplice. The approved module controls add two nullable runtime
overrides: blank inherits the YAML configuration; standard effective defaults
remain 20 events and merge height 0.25. The module kernel is unchanged. Both
report entrypoints have explicit version changes to record actual parameters,
without inferring historical values or rerunning scientific stages. The approved
network annotation policy adds `apply_threshold_labels=false` to CoSplice, the
imported-PSI workflow and explicit hub analysis. Continuous metrics and complete
ranks remain; opt-in labels cannot delete events or change fitted assignments.
Other methods, including MotifTrace and DEG runtime parameters, remain unchanged.

The approved DAS selector defaults to RS-DAS and preserves the ordinary-DAS
flag, native effects and P/FDR. Tests explicitly allow the downstream selection
metadata and the detected/testable-background controls without weakening the
read-support contract. Synthetic MotifTrace selections declare the separate
native-significance flag and the every-sample ten-count policy; historical
stable flags alone remain invalid. Both earlier workflow inventory files remain
unchanged, with the permitted version/default differences listed in the test.
Interface hashes in the MotifTrace integration record are revisioned separately
from its 23 frozen scientific sources, retaining prior interface digests.

The tests use generated synthetic inputs. The API test opens a project, executes
a small DEG result-reader workflow, verifies reuse of its completed result,
assembles and downloads a report containing HTML, PDF, provenance and result
files, and checks that input files remain unchanged. The resource checks also
verify the retained network R sources, original vendor license and report
templates.

Passing this suite does not establish execution of optional R tools, reference
databases, new network fits or an experimental study. Those need their declared
software and inputs. Synthetic narrative-report, MotifTrace and network
mathematics checks are retained separately in the main test suite.
