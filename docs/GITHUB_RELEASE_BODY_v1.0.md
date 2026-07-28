# Katana v1.0

Katana is a comparison-first, direction-aware desktop workflow for RNA-seq and alternative splicing candidate screening.

This release establishes a modular GUI workflow for:

- project-specific input discovery
- pairing-aware direction control
- per-comparison transcript + splicing integration
- candidate gene ranking and selection
- cross-comparison candidate follow-up
- event-driven sashimi review
- Jutils output browsing
- RBP support and cross-comparison RBP support plots
- isoform candidate support and cross-comparison isoform support plots

## Highlights

- Comparison-first workflow with per-comparison pages, outputs, and state
- Canonical splicing-side direction inherited by rMATS, DEXSeq, and DTU
- Read-only canonical AS and DEG direction audit in Pairing, with downstream standardized effects
- Backward-compatible loading of legacy saved direction transforms
- Independent module execution with no automatic chained downstream runs
- Gene-centered candidate screening and follow-up
- Event selection before sashimi generation
- Project-specific annotation loading instead of hard-coded dataset mapping
- Cached TSV/PNG preview tabs for downstream RBP and isoform support modules

## Included in this release

- Complete README with installation, workflow overview, expected inputs, outputs, testing, and known limitations
- Release notes for v1.0
- GitHub repository description suggestions

## Current known limitations

- Some PDF outputs are opened externally instead of embedded inline
- Sashimi plotting depends on the local `rmats2sashimiplot` runtime and Python dependencies such as `pysam`
- Optional modules still depend on project-specific resource availability
- RBP and isoform support plots are prioritization views, not causal regulatory models

## Author

Prepared for release by **Tara**.
