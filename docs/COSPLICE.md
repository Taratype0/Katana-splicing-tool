# CoSplice

CoSplice builds one network from all comparable samples in a Katana project.
Each sample remains an individual PSI column. WGCNA discovers modules; kME
describes event membership. Module associations use the comparisons selected
by the user. Changing the report comparison does not rebuild the network.

## Run on an existing project

Use the matching Katana container or the standalone Linux Conda environment.
The Python wheel includes the engine and report code; native computation also
requires R and WGCNA from that environment.

```bash
katana-cosplice --project project.json plan
katana-cosplice --project project.json run
katana-cosplice --project project.json artifacts
```

The project can supply a PSI matrix and event annotation, native rMATS results
with their sample identities, or an existing compatible CoSplice input config.
Gene-set resources are supplied separately; no example-study database is used.

## Pairing and sample counts

Pairing is optional. A contrast declared `paired=true` retains the matching
`biological_sample_id` values. For a standalone project whose comparisons have
not declared pairing, explicitly enable it with `run --paired`.

Repeated observations count once when reporting biological units. Six samples
from three paired individuals are three biological units. There is no minimum
unit-count gate for network construction. Twelve independent units are a project
recommendation, not a significance threshold or a literature-derived guarantee.
A statistic that cannot be estimated remains NA with its reason.

## Parameters

```json
{
  "event_limit_mode": "5000",
  "min_module_size": 20,
  "merge_cut_height": 0.25,
  "bootstrap_iterations": 0,
  "association_contrast_ids": ["Treated_vs_Control"]
}
```

```bash
katana-cosplice --project project.json run --parameters network.json
```

Event presets are 2,000, 5,000 or 10,000; `Custom` uses `max_events`.
Selection follows input QC and variance, not DAS significance. The association
list limits downstream comparisons only: discovery still uses all project
samples. Omitting the list uses all declared comparisons. Bootstrap is off by
default. The report records actual parameters and retains continuous membership
values, complete gene/event lists and association statistics.

## Python and Nextflow

```python
from katana_cosplice.study_api import CoSpliceStudy

with CoSpliceStudy("project.json") as network:
    print(network.plan())
    result = network.run()
```

The standard [Nextflow command](NEXTFLOW_SIMPLE.md) uses this same workflow and
installed engine. `network.engine=whole_study_4_3` is the default. The older
per-comparison implementation remains explicitly available as `engine=legacy`
and `katana-cosplice-legacy`; it is not silently substituted for the new engine.

Results include the native WGCNA model, module assignments, kME, module scores,
comparison associations, enrichment, network browser and complete tables. The
approved seven-chapter report ends with this shared study network.
