# Katana Agreed Requirements

This file captures the currently agreed workflow and UI expectations so the
implementation can be checked against the user's requirements without losing
earlier requests.

## Workflow

1. Project
   - Load project inputs.
   - Confirm project before moving forward.
   - Input-file checking must be available and easy to read.
2. Pairing
   - Review all detected comparisons.
   - Confirm per-comparison splicing/DEG direction.
   - Confirm pairing before moving forward.
   - Changes here must propagate downstream.
3. Comparison Sets
   - Define comparison-vs-comparison analyses.
   - Each set needs custom name support.
   - Meanings and labels must inherit from Pairing.
   - Confirm comparison sets before analysis.
4. Analysis
   - Not one monolithic page.
   - Each major analysis branch must be separate.
   - Each branch can have sub-pages / sub-options.
   - Each sub-analysis should be runnable separately.
5. Cards / Sashimi / Candidate Hooks
   - Card generation must depend on user-selected genes and comparison groups.
   - Sashimi should be grouped similarly.
   - A dedicated candidate-gene list hook is required for downstream AI /
     enrichment / biology interpretation.

## Analysis branches

### 1. Landscape

- Parent branch only, no mixed content.
- Separate runnable sub-pages:
  - Bar plot
  - Pie chart
  - Jutils heatmap
  - Jutils PCA

### 2. Group-vs-group splicing comparison

- Separate branch.
- Should support the user's comparison-set design.
- Needs a spreadsheet / Excel-like preview like the provided scripts.
- Used for candidate-gene selection.

### 3. Transcript-vs-splicing integration

- Separate branch.
- User can choose how many experiments/comparisons participate.
- Each selected experiment should appear as its own sub-item.
- This is where thresholds matter most.
- Must be interactive:
  - show gene identity
  - show direction
  - show DEG significance / fold change
  - show splicing PSI / significance
- Must inherit direction from Pairing.

### 4. Cards / Sashimi grouped visualisation

- Separate branch.
- User defines experiment groups, e.g.:
  - 1/2/3/4 as one group
  - 2/3 as another group
- Each group should then expose sub-actions:
  - card visualisation
  - sashimi visualisation
- Cards depend on chosen genes.
- One gene = one card output.

## Input validation

- Input checking must be readable and scrollable.
- It should show loaded files per comparison / experiment.
- It must help verify that detected files are complete and correct.

## Thresholds

- Thresholds must refresh immediately after project load.
- Thresholds should use literature-style presets.
- Thresholds are not equally important for every branch; branch-specific use is preferred.

## Non-negotiable UX points

- Do not silently omit requested visualisations.
- Do not merge unrelated analyses into a single mixed page.
- Do not run every analysis automatically when the user only wants one.
- Preserve inheritance from:
  Project -> Pairing -> Comparison Sets -> downstream analysis branches.
