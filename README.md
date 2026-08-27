# OLFML2B-STAD

Reproducible analysis workflow supporting the manuscript:

**Cross-cohort single-cell and spatial triangulation positions OLFML2B as a fibroblast-dominant stromal-ecology marker in gastric adenocarcinoma**

## Authors

Lei Huai¹*, Xin-Zhen Cai²*, Zhang-Hui Long¹*, Yan-Xi Fu¹*, Yuan-Sheng Zang¹

¹ Department of Medical Oncology, Changzheng Hospital, Naval Medical University, Shanghai, China  
² Department of Hematopathology, People's Hospital of GaoYou County, Yangzhou, Jiangsu, China  
*Equal contribution. Correspondence: Yuan-Sheng Zang (doctorzangys@163.com).

## Scope

This repository contains the complete Part0-Part9 R workflow used for the manuscript's bulk-transcriptomic, survival, tumour-microenvironment, proteomic, single-cell, spatial-transcriptomic, and exploratory immunotherapy analyses.

The public release is intentionally code-centric. Raw public datasets, restricted/redistribution-sensitive source files, large caches, intermediate R objects, logs, and generated output are **not committed**. Public datasets are downloaded by the workflow where feasible; inputs whose redistribution is inappropriate must be obtained from their original repositories and placed in `DATA_INBOX/` as described below.

## Repository structure

```text
OLFML2B_STAD/
├── 00_RECOVER_AND_RUN_OLFML2B_FROM_ZERO_v1_0_2.R
├── 00_INSTALL_REQUIRED_PACKAGES_ONCE.R
├── R/
│   ├── 00_OLFML2B_PART0_CONFIG_CORE.R
│   ├── 01_OLFML2B_PART1_TCGA.R
│   ├── 02_OLFML2B_PART2_GEO.R
│   ├── 03_OLFML2B_PART3_BULK_SURVIVAL.R
│   ├── 04_OLFML2B_PART4_IMMUNE_TME_PRODUCTION.R
│   ├── 05_OLFML2B_PART5_PDC_PRODUCTION.R
│   ├── 06_OLFML2B_PART6_SCRNA_PRODUCTION.R
│   ├── 07_OLFML2B_PART7_SPATIAL_TRANSCRIPTOMICS.R
│   ├── 08_OLFML2B_PART8_ICI_PRJEB25780_TIGER_ONLY_MOLECULAR_CONTEXT.R
│   └── 09_OLFML2B_PART9_PUBLICATION_OUTPUTS_v8_10_1.R
├── DATA_INBOX/
│   ├── PDC000614/README.md
│   └── TIGER/README.md
├── docs/
├── CITATION.cff
├── .zenodo.json
├── LICENSE
└── .gitignore
```

## Quick start

The manuscript analyses were completed with R 4.5.2; R 4.5.x is recommended.

1. Clone or download this repository.
2. If the launcher reports missing R packages, run once:

```r
source("00_INSTALL_REQUIRED_PACKAGES_ONCE.R", encoding = "UTF-8", local = FALSE)
```

3. Place the manual PDC000614 and TIGER inputs in the corresponding `DATA_INBOX/` directories (see their README files).
4. Start the full recovery and Part0-Part9 workflow:

```r
source(
  "00_RECOVER_AND_RUN_OLFML2B_FROM_ZERO_v1_0_2.R",
  encoding = "UTF-8",
  local = FALSE
)
```

The launcher resolves the project root from its own location by default. A custom root may be supplied with the `OLFML2B_STAD_ROOT` environment variable.

## Main data layers

- TCGA-STAD: bulk RNA and clinical data
- GSE62254, GSE15459, GSE26253, GSE84437: external bulk cohorts
- GSE134520, GSE150290, GSE167297, GSE183904: single-cell RNA-seq
- GSE251950: spatial transcriptomics
- PDC000614: proteomics
- PRJEB25780/TIGER: exploratory anti-PD-1 response boundary

See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) for acquisition and redistribution notes.

## Reproducibility principles

- Patient/cohort-level inference is used where the data structure permits it.
- Cells and spatial spots are not treated as independent patient replicates.
- OLFML2B is removed from any programme score containing it before target-programme association analyses.
- Continuous within-cohort standardisation is used for survival models; median splits are display-only.
- Proportional-hazards diagnostics, random-effects synthesis, prediction intervals, multiplicity control, leave-one-out analyses, and claim ceilings are retained.
- Negative and precision-limited results are preserved rather than rescued by post hoc cut-points.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) and [`docs/PUBLIC_RELEASE_SANITIZATION.md`](docs/PUBLIC_RELEASE_SANITIZATION.md).

## Outputs

Generated data, figures, tables, audit files and logs are written under project-local `output/`, `recovery/` and `logs/` directories and are excluded from Git by default. Part9 produces the publication figures, supplementary tables, source-data exports and audit files from the upstream frozen result tables.

## Citation

Please use the metadata in [`CITATION.cff`](CITATION.cff). After a GitHub release is archived in Zenodo, update the repository metadata with the Zenodo DOI and cite the frozen release corresponding to the manuscript.

## License

Code in this repository is released under the MIT License unless a third-party dependency or source dataset specifies different terms. Source datasets remain subject to the terms of their original repositories.
