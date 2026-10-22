# Development input contracts

These modules add read-only checks for project manifests, sample identity and
user-selected comparisons to the existing Katana software. Established nf-core
pipelines retain responsibility for reads, lane merging, alignment, quantification
and native input checks. Katana records experimental conditions, biological units,
donors and requested contrasts alongside those outputs.

The code lives in `katana_next/contracts/` in the main package. It is development
preflight support; the existing GUI and scientific workflows keep their current
interfaces. Passing these checks establishes metadata consistency. Native model,
reference, measurement and execution eligibility require their own checks.

## Run from a source checkout

Run these commands from the repository root. For contract-only source commands,
the dependencies are:

```bash
python -m pip install 'jsonschema>=4.18,<5' 'pydantic>=2.9,<3' 'numpy>=1.26,<3' 'pytest>=8.3,<10'
```

Full application installation follows the repository software instructions.
The `contracts` install extra adds validation dependencies to the current package;
it does not remove the package's existing application dependencies.

```bash
python -B -m katana_next.contracts validate examples/contracts/v1/bulk-native-import.project.json --verify-hashes
python -B -m katana_next.contracts samplesheets --input examples/contracts/sample-identity-v1/fastq-multilane/upstream.csv --metadata examples/contracts/sample-identity-v1/fastq-multilane/metadata.csv --format fastq
python -B -m katana_next.contracts contrasts --metadata examples/contracts/selected-contrasts-v1/selected-subset.metadata.csv --contrasts examples/contracts/selected-contrasts-v1/selected-subset.contrasts.csv --route bulk_complete
python -B -m pytest -q -p no:cacheprovider tests/contracts
```

The examples are synthetic. A condition list does not generate every pairwise
comparison: only the user-supplied contrast rows are checked. Multiple lanes do
not add biological replicates, and donor pairing requires explicit identity.

- [Samples, lanes and biological units](sample-identity-v1.md)
- [Selected contrasts and method boundaries](selected-contrasts-v1.md)

The complete bulk route currently declares independent `~ condition` designs.
More elaborate designs retain separate adapter and integration states in the
versioned method-support table. These checks run no alignment or statistical fit,
and validation never authorizes execution or certifies a completed scientific run.
