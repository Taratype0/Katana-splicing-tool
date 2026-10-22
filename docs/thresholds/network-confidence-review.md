# Network confidence, hubs and statistical labels

**Status: stage 3 approved on 2026-10-22; labels are optional and off by default.**
The approved output policy is `network-descriptive-labels-opt-in-v1`:
retain continuous measurements and complete rankings, with no default rejection
by the heuristic membership or hub cutoffs below. `apply_threshold_labels=true`
adds descriptive labels only; it does not remove rows or change the fitted
partition. A disabled label is `not_applied`, not a failed candidate or an absent
hub. The standard statistical reporting level remains FDR 0.05, with the actual
test family and native values retained.

The tables below document the **legacy optional labeling presets**, not newly
validated scientific boundaries. The stage-2 defaults, **20 splice events** and
**merge height 0.25**, remain approved and configurable. Previously completed
results retain their original policy. See the [active policy register](network-labels.v1.json).

## Member confidence annotates the fixed modules

The production default is `backbone_confidence_only`. Assigned events remain in
their full-data WGCNA modules even when their confidence label is `uncertain`.
Confidence is not a P-value, FDR, posterior probability or proof of function.
See [`consensus.py:509–525`](../../katana_cosplice/consensus.py#L509).

| Existing setting | Default | Meaning |
| --- | ---: | --- |
| Core co-assignment | >=0.70 | Mean bootstrap co-assignment with other members of the fixed module; self-pairs are excluded. |
| Peripheral co-assignment | >=0.45, margin >=0.05 | Own-module support must exceed the strongest alternative by at least 0.05. |
| Alternative core kME evidence | mean absolute kME >=0.70; 10th percentile >=0.50; sign consistency >=0.80 | All three conditions and a positive co-assignment margin are required. |
| Alternative peripheral kME evidence | mean absolute kME >=0.50; 10th percentile >=0.30; sign consistency >=0.70 | All three conditions and a positive co-assignment margin are required; this does not downgrade a co-assignment core. |
| Production confidence bootstrap count | 50 | The application overrides the standalone `ConsensusConfig` default of 20. The imported-PSI workflow also uses 50. |

The co-assignment and kME routes are **alternatives (OR)**, not two mandatory
tests. An initially supported label must then pass confirmation on two disjoint
halves of the bootstrap runs. A core label must be core in both halves;
otherwise a supported event can be peripheral only if both halves support at
least peripheral membership. Failure gives `uncertain`, without removing the
event or refitting its module. Fewer than four runs gives `not_evaluated` for
the half-run check; this is not independent-cohort validation.

Sources: [`config.py:211–239`](../../katana_cosplice/config.py#L211),
[`consensus.py:113–190`](../../katana_cosplice/consensus.py#L113),
[`consensus.py:208–288`](../../katana_cosplice/consensus.py#L208),
[`consensus.py:930–1003`](../../katana_cosplice/consensus.py#L930),
[`programs.py:237,324–325`](../../katana_next/workflows/programs.py#L237),
[`psi_network.py:115`](../../katana_next/science/psi_network.py#L115).

A separate continuous confidence score combines available co-assignment, kME,
split-half and optional read-reliability components geometrically. Its evidence
floor is 0.05; split-half component scores are 0.05/0.35/0.75/1.0. These are
Katana scoring choices, not probabilities or additional inclusion cutoffs. Missing
domains are recorded, not presented as measured evidence. See
[`membership_validation.py:15–116`](../../katana_cosplice/membership_validation.py#L15).

## Hub rankings and independently assessed hub candidates differ

The built-in module-interpretation path ranks **all** assigned events using the
mean of within-module percentiles of absolute kME, co-assignment confidence and
connectivity. The gene score averages percentiles of maximum absolute kME,
maximum event score and the number of module events. These are Katana descriptive
scores, not published composite statistics. A plot's **top 10 events** is a
display limit, not a hub significance threshold. No hard kME cutoff is applied
by this ranking path. See
[`module_interpretation.py:146–209,485–490`](../../katana_cosplice/module_interpretation.py#L146).

The explicit `programs.hub_events` workflow assesses candidates against fixed
modules using observed PSI, official R WGCNA correlation/PC1 and complete
within-module connectivity. Existing defaults are:

| Candidate criterion | Current default |
| --- | --- |
| Connectivity rank | Top 10% within each module and sample scope; average tied ranks; cutoff `max(1, ceil(module_events * 0.10))`. |
| Module membership | Signed kME >=0.70; absolute kME >=0.70 for an unsigned source network. |
| Biological-block bootstrap | 100 iterations; requested invalid iterations remain in the selection-frequency denominator. |
| Bootstrap selection frequency | >=0.70. |
| Valid bootstrap fraction | >=0.80. |
| Estimable module-pair fraction | >=0.80. |
| Supported connections | Positive connectivity and at least one estimable connection at the source adjacency threshold, default 0.10. |
| Same-gene sensitivity | Both all-edge and cross-gene criteria must pass; cross-gene candidates require complete host-gene annotation and possible cross-gene connections. |
| Biological units per scope | At least `max(15, source network minimum)`. Insufficient scopes are not assessed as supported hubs. |
| Requested scopes | Overall and each condition, assessed separately; these are not independent replications. |

Sources: [`hub_events.py:52–82`](../../katana_next/workflows/hub_events.py#L52),
[`science/hub_events.py:69–133`](../../katana_next/science/hub_events.py#L69),
[`hub_metrics.R:20–73`](../../katana_next/science/hub_metrics.R#L20).

**Adjacency 0.10 is not an edge P-value.** Full kWithin sums every estimable
within-module adjacency without this truncation. The threshold nevertheless
affects the supported-degree requirement above, as well as drawing exports.
The drawing exports at most eight connections per event; that limit does not
alter the full connectivity calculation. The general module dashboard does not
infer validated hubs from its thinned drawing. See
[`module_dashboard.py:73–84`](../../katana_next/services/module_dashboard.py#L73).

The separate discovery-validation table also reports how often an event is in
the top kME decile. That diagnostic is not the explicit hub-candidate rule above
([`validation.py:127–140,190–194`](../../katana_cosplice/validation.py#L127)).

## Statistical labels retain their test families

The existing conventional significance label is adjusted P/FDR **<=0.05** where
the corresponding method computes a valid test. This does not make membership,
stability, connectivity or candidate ranks statistically significant. All tested
rows remain available; display ordering is not a significance filter.

| Analysis | Actual correction family |
| --- | --- |
| Module-trait association | Within one observed/residual network run, BH across all returned module x factor x omnibus/level-vs-rest tests together. |
| Factorial module association | BH across all returned module x model-term tests; coefficient tests form a separate family. |
| Count-model event condition effect | BH across successfully fitted event tests; separate from module-score models. |
| Module-score condition effect | BH across successfully fitted module omnibus tests; `condition_driven` uses `statistics.fdr_cutoff`, default 0.05. |
| Module-expression association | Overall and condition-adjusted results form separate BH families. Each includes all module x pathway/RBP-feature tests; pathway and RBP features are not corrected separately. |
| Explicit module GMT enrichment | BH across all module x term tests within each library, including zero-overlap terms. Background includes all mapped, measured network genes, including unassigned events. |
| Built-in module interpretation | Membership correlations, motif enrichments, RBP-expression correlations and GMT enrichments each have their own table-wide BH family. A configured external enrichment provider retains its own correction semantics. |

Sources: [`association.py:167,258–262`](../../katana_cosplice/association.py#L167),
[`module_models.py:262–264,330–332`](../../katana_cosplice/module_models.py#L262),
[`module_expression.py:32–35`](../../katana_next/science/module_expression.py#L32),
[`module_functions.py:21–40`](../../katana_next/workflows/module_functions.py#L21),
[`module_interpretation.py:174,297–299,346–348,376–389`](../../katana_cosplice/module_interpretation.py#L174).

Functional labels require enrichment FDR <=0.05. A report may also show the ten
lowest-FDR terms even if none is significant; it states this explicitly
([`module_explanation.py:73–77,100–104`](../../katana_next/workflows/module_explanation.py#L73)).
The dashboard currently colors terms with **FDR <0.05**, whereas the functional
label uses **<=0.05**; this boundary mismatch is recorded here, not silently
changed by this review
([`ModuleDashboard.tsx:99–100`](../../apps/web/src/components/ModuleDashboard.tsx#L99)).

The existing module-validation score is the minimum of mean matched-module
Jaccard similarity and mean event-membership stability. Labels are `robust`
at >=0.75, `moderate` at >=0.50 and `unstable` below 0.50; they are resampling
descriptions, not P/FDR. This optional validation stage has a standalone default
of 20 runs, distinct from the production confidence count of 50 and hub count
of 100. See [`config.py:201–205`](../../katana_cosplice/config.py#L201),
[`validation.py:89–94,195–237`](../../katana_cosplice/validation.py#L89).

## Basis and decision boundary

- [WGCNA signedKME documentation](https://search.r-project.org/CRAN/refmans/WGCNA/html/signedKME.html)
  defines kME as correlation with a module eigengene.
- [WGCNA intramodular connectivity documentation](https://search.r-project.org/CRAN/refmans/WGCNA/html/intramodularConnectivity.html)
  defines connectivity as the sum of within-module adjacency, excluding the diagonal.
- [R `p.adjust` documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
  describes BH adjustment; it does not require the chosen reporting alpha of 0.05.
- [fpc `clusterboot` documentation](https://search.r-project.org/CRAN/refmans/fpc/html/clusterboot.html)
  gives bootstrap Jaccard interpretation guidance, including 0.75 and 0.50.
  Katana applies these numbers to its own combined score; its `moderate` label
  is not identical to the author's recommended stable-cluster interpretation.
- [WGCNA methods paper](https://doi.org/10.1186/1471-2105-9-559)
  supports module, hub and trait-association concepts. It does not establish
  Katana's exact OR rule, split-half rule, composite scores, hub cutoff values
  or co-splicing-specific validity.

The exact confidence, hub and composite-score values above are **Katana
operational heuristics**, not universal literature standards. Existing external
membership checks and hub simulation evaluations do not establish that these
specific numeric cutoffs are calibrated or optimal. Under the approved policy
they are optional labels, not default inclusion rules. Association is not
causation; a hub event is not automatically a regulatory protein or a causal
driver. Same-cohort resampling is not independent replication.

The legacy composite confidence score includes components derived from discrete
split-half labels. It must therefore remain identified as threshold-dependent
legacy evidence, rather than being presented as a threshold-free probability or
the required default ranking. Individual raw measurements and their full ranks
remain available. This policy does not change the WGCNA fit, resampling formulas,
native statistical tests or completed studies.
