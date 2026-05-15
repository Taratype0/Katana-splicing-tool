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

## Highlights

- Comparison-first workflow with per-comparison pages, outputs, and state
- Canonical splicing-side direction inherited by rMATS, DEXSeq, and DTU
- Manually selectable DEG final direction when DESeq2 naming is unreliable
- Independent module execution with no automatic chained downstream runs
- Gene-centered candidate screening and follow-up
- Event selection before sashimi generation
- Project-specific annotation loading instead of hard-coded dataset mapping

## Included in this release

- Complete README with installation, workflow overview, expected inputs, outputs, testing, and known limitations
- Release notes for v1.0
- GitHub repository description suggestions

## Current known limitations

- Some PDF outputs are opened externally instead of embedded inline
- Sashimi plotting depends on the local `rmats2sashimiplot` runtime and Python dependencies such as `pysam`
- Optional modules still depend on project-specific resource availability

## Author

Prepared for release by **Tara**.
