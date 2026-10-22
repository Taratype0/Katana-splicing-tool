# Usable application acceptance

Run from the repository root after installing the application and test dependencies:

```bash
python -m pytest -q tests/public_baseline
```

The workflow inventory is copied from the preserved English-report application
baseline. It fixes the 54 workflow identifiers, method versions and public
parameter defaults, including MotifTrace, network analysis and both report
entrypoints.

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
