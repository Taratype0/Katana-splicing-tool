# Current development container — 23 October 2026

This records installation of the current packaged workflow, including the
per-comparison network/MotifTrace report changes. It uses software fixtures and
an already accepted result collection. No study, public sequencing dataset or
completed scientific calculation was rerun.

## Software identity

| Item | Identity |
| --- | --- |
| Software source | `b0fbd7223b392867b4098d59a2a75dc992a80170` on `codex/package-development` |
| Wheel | `katana_splicing_tool-2.0.0b5.dev20261023-py3-none-any.whl` |
| Wheel SHA256 | `2d979b97613e01b7cc6f79ddbc62058ecc55c9dbdcd80cc7abc1f7e16a229663` |
| Local image tag | `katana-downstream:dev-b0fbd72` |
| Image ID | `sha256:8313b9da1c43ecd8d2df6f7576ccdb8d17327cda7b42069874e0d2a4afe6d199` |
| Export archive | `katana-downstream-dev-b0fbd72-linux-amd64.tar.gz` (1,459,686,331 bytes) |
| Archive SHA256 | `cc56d2654285180d398a126e6772b6695dd7d84c6e5136721ab799e4872acabe` |
| Platform | `linux/amd64` |
| Native lock SHA256 | `cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99` |
| Pip lock SHA256 | `4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1` |

The fixed base and dependency locks are unchanged from the
[earlier native-runtime acceptance](downstream-runtime-20261023.md). Verified
dependency build layers were reused; the current wheel was installed into a
new image. Acceptance ran as UID 1000 with networking disabled and without a
host Python/R environment mount. This is not a fresh dependency solver run.

## Observed checks

| Check | Result |
| --- | --- |
| Installed current wheel outside the source checkout | 41 portable checks passed; CLI, pipeline export and browser assets checked |
| Current container native runtime | 51 passed, zero failures |
| Installed package content | 522 package/resource files matched the current wheel byte-for-byte |
| Offline browser service | Health check and six HTTP routes/assets passed; served JavaScript/CSS matched the wheel |
| Actual Nextflow 25.10.4 / Docker execution | Three independent runtime probes completed with exit zero |
| Nextflow resume | All three probes cached; output hashes preserved |
| Container report rendering from existing results | Positive report acceptance; all 98 task records verified, no calculations rerun |
| Exported image round trip | Archive hashed, loaded with Docker and image ID verified unchanged |

The report check retained the approved seven-section specification. It used the
existing two-comparison fixture collection, including independently fitted
networks, multiple MotifTrace models, Top 10, native significance values,
compact enrichment summaries and complete hit-gene exports. The verified task
records include 42 completed calculations/adapters and 56 explicitly
not-applicable MAPP preparation/fit tasks. They are not 98 new successful fits.
The [generic workflow record](generic-pipeline-20261023.md) describes the
underlying result acceptance and its limits.

## GitHub CI repair

The native scheduling check initially failed because its simulated report
omitted `flow.acceptance.json`, which the production wrapper now requires.
The fixture now supplies a clearly marked transport-only acceptance record.
The production report gate was preserved. The repaired
[Native Nextflow integration run](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40025766609)
passed, including task dependency scheduling, failure propagation, cache resume
and the reference-exclusion fixture.

The same source commit also passed the
[Katana tests run](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40025766589):
Linux and Windows application/report/browser regression, interface tests/build,
and installed-wheel checks outside the checkout. The following documentation
update changes no software or image inputs and retains these tested identities.

## Use and limits

Follow the [installation guide](../INSTALL.md) or the
[container build/load instructions](../../containers/downstream/README.md).
The image is a local development artifact, not a published registry release.
No PyPI release is implied. Container installation requires no host Conda;
standalone downstream users can install the matching Conda/pip locks.

The tested container platform is Linux x86-64 using Docker. This run does not
certify native ARM/macOS execution, Apptainer/Singularity image conversion,
every HPC site's scheduler configuration, new FASTQ/BAM scientific runs or a
real-study report. References and databases remain separately configured
inputs. The paused Primary culture analysis and the original application were
not modified or executed.
