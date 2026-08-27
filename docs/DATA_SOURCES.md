# Data sources and acquisition contract

| Layer | Source | Role | Acquisition / redistribution note |
|---|---|---|---|
| Bulk discovery | TCGA-STAD | Tumour-normal expression, ecology, OS | Acquired programmatically through the GDC/TCGAbiolinks workflow |
| Bulk validation | GSE62254 | Ecology, stage, OS, DFS | Public GEO source |
| Bulk validation | GSE15459 | Ecology, stage, OS | Public GEO source |
| Bulk validation | GSE26253 | Ecology, stage, RFS | Public GEO source |
| Bulk sensitivity | GSE84437 | Ecology, OS sensitivity | Public GEO source; frozen subseries mapping is audited |
| scRNA-seq | GSE134520 | Source localisation | Public GEO supplementary archive |
| scRNA-seq | GSE150290 | Source localisation and paired patient evidence | Public GEO supplementary archive |
| scRNA-seq | GSE167297 | Source localisation | Public GEO supplementary archive |
| scRNA-seq | GSE183904 | Source localisation and paired patient evidence | Public GEO supplementary archive |
| Spatial | GSE251950 | Spatial contextualisation | Public GEO supplementary archive; 10 sections / 9 patients |
| Proteomics | PDC000614 | Orthogonal protein direction | Obtain from PDC; source files are not redistributed in this repository |
| Anti-PD-1 | PRJEB25780/TIGER | Treatment-response claim boundary | Obtain the frozen TIGER RDS files from the original source; not redistributed here |

## Heavy GEO archives recovered by the launcher

The current from-zero recovery layer is configured for the supplementary archives of GSE150290, GSE183904, GSE167297, GSE134520 and GSE251950, with resumable downloads.

## Manual inputs

See:

- `DATA_INBOX/PDC000614/README.md`
- `DATA_INBOX/TIGER/README.md`

The `.gitignore` rules intentionally prevent actual PDC/TIGER payloads from being committed.
