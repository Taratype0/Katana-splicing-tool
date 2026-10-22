# Historical bulk policy — not active

These files preserve the exact documentation bytes from commit `8c368b74cba5a7a7b25ca9edea8c40231b9d598f`. They describe the superseded replicate-consistency presets and optional approximate-SE filter. They are an archive of what that revision documented, not current defaults, user choices or validation claims.

The active policy is **`bulk-all-samples-counts-v1`**: every expected sample in both arms of a selected contrast must have finite informative counts `IJC + SJC >= 10`, with unambiguous sample mapping. See [the current threshold reference](../../../THRESHOLDS.md) and [current register](../../bulk-defaults.v1.json).

The archive manifest records the source commit and SHA256 of each unchanged historical file. Reading these archived files does not change a current run or migrate historical results.

Archived originals: [threshold reference as unmodified text](THRESHOLDS.original.txt), [machine-readable register](bulk-defaults.v1.json), and [archive manifest](archive-manifest.json). Historical relative source links in the original text refer to its original commit, not this archive folder.
