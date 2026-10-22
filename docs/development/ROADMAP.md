# Development roadmap

This roadmap separates the usable application baseline from work needed for the package/pipeline product.

| Stage | Status and intended outcome |
| --- | --- |
| Software baseline | Complete application modules, MotifTrace, networks and existing report routes preserved in source control; private datasets and workstation configuration excluded. |
| Scope and architecture | Function inventory, compatibility boundaries, defaults and preservation rules documented. These are design decisions, not completed package migration. |
| Shared input contracts | Project, sample/lane/donor and selected-contrast validation developed separately. These validators do not replace mature upstream algorithms. |
| Analysis package extraction | Separate reusable scientific services from application/GUI dependencies; retain native method identities and existing behavior. Pending module-by-module migration and acceptance. |
| Scientific adapters | Make declared input/reference requirements explicit and qualify new computation versus imported evidence, including MotifTrace, networks and single-cell routes. |
| English reports and interface | Preserve the existing report structure and connect consumers to the same result contracts; qualify complete portable deliveries. |
| Packaging and Nextflow integration | Build reproducible distributions and connect mature bulk upstream processing to the reusable downstream package. Not completed merely by publishing the source baseline. |
| Acceptance and release | Validate representative workflows, both FASTQ/BAM entry routes where supported, portable outputs, installation and reproducibility before claiming a new release. |

Main presents the usable software baseline. Package restructuring continues on the development branch. Changes to algorithms, defaults, result identities and presentation are reviewed separately; completed scientific outputs are reused only when their compatibility is established.

Internal task numbers are bookkeeping identifiers for this work. They are not software features or evidence that an entire product stage has shipped.
