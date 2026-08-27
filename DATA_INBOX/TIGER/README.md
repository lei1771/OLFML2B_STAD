# PRJEB25780 / TIGER manual input

Do **not** commit the source RDS payloads in this directory.

The frozen Part8 source semantics require the two TIGER files:

- `STAD-PRJEB25780.Response.Rds`
- `STAD-PRJEB25780.Response (1).Rds`

Place them here before running the from-zero launcher. The launcher copies them into `data/cache/ICI_PRJEB25780/TIGER/` during local execution.

Do not silently substitute a separately reprocessed ENA expression matrix, because that would change the frozen source semantics used by Part8.
