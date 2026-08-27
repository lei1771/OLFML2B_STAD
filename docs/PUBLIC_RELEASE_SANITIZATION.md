# Public-release sanitation note

This GitHub release is designed to be portable and safe to publish. The scientific analysis contract is unchanged relative to the manuscript release.

Environment-only cleanup in the public package:

1. **Part9 v8.10.1** explicitly namespaces all ten `bind_rows()` calls as `dplyr::bind_rows()` to prevent masking by a Part2 project-local helper during a full Part0-Part9 session. No statistical/data/figure rule changed.
2. **Part8** no longer discovers data from an unrelated historical cross-project cache fallback. Canonical project-local discovery and the explicit `OLFML2B_STAD_SHARED_DATA_ROOT` option remain. No Part8 analytical method changed.
3. **From-zero launcher v1.0.2** resolves the project root from the launcher's own directory by default, making the release portable across Windows/macOS/Linux directory layouts. It also removes a duplicated error-message string from the prior launcher. No analysis step changed.
4. Raw data, restricted/manual inputs, large caches, generated R objects, logs and local outputs are excluded from the public repository by design.

The SHA256 manifest in the repository records the exact public-release files.
