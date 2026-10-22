# Katana T09 contract development checkpoint

This allow-listed source snapshot publishes the new read-only exchange-contract
boundary plus its unchanged import dependencies. It does **not** contain the
complete new analysis package or a Nextflow runner. The root-level application
is legacy code; commands below use only `development/source`.

## Run the validator from source

Run commands from `development/`.

Use Python 3.11 or newer and an environment of your choice. Install these small
contract/test dependencies rather than installing the legacy root distribution:

```sh
python -m pip install -r requirements-contracts.txt
```

PowerShell:

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) 'source'
$env:PYTHONDONTWRITEBYTECODE = '1'
python -B -m katana_next.contracts validate source/examples/contracts/v1/bulk-native-import.project.json --verify-hashes
python -B -m pytest -q -p no:cacheprovider source/tests/contracts
```

Linux/macOS:

```sh
PYTHONPATH=source python -B -m katana_next.contracts validate source/examples/contracts/v1/bulk-native-import.project.json --verify-hashes
PYTHONPATH=source PYTHONDONTWRITEBYTECODE=1 python -B -m pytest -q -p no:cacheprovider source/tests/contracts
```

Exit code 0 means the requested metadata/file checks passed; 2 means invalid
input or failed checks. A consistent blocked/missing-result manifest can validate
successfully. Validation is not biological completion or execution permission.

## Included

Exchange schema `1.0.0`, project/asset/run/result identity and provenance checks,
sample/contrast/method/unit/path/missing-state rules, three synthetic projects,
and the original 88 contract/API/CLI tests. `domain/` and the parent initializer
are unchanged import dependencies; their legacy platform schema `2.0` remains
separate from the exchange schema. The source package version label identifies
the development source; this subset is not an installable scientific release.

The validator performs no project creation, catalog write, GUI launch, remote
fetch, Nextflow submission, model fitting or scientific testing. Optional local
file/size/hash verification is explicitly requested. Examples contain synthetic
identities and tiny fixtures, never experimental data.

This subset is MIT-licensed under the root LICENSE. It does not vendor scientific
tools, annotation databases, front-end fonts or other external resources.
