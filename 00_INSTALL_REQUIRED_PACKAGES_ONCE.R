# ==============================================================================
# OLFML2B-STAD environment recovery helper
# Run ONLY if the main launcher reports missing R packages.
# This script is deliberately separate from the scientific analysis runner.
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "data.table","ggplot2","matrixStats","readxl","R.utils","jsonlite","curl","httr",
  "yaml","survival","metafor","sandwich","lmtest","mice","digest","logistf","Matrix",
  "dplyr","tidyr","readr","forcats","stringr","tibble","purrr","scales","patchwork",
  "ggrepel"
)

bioc_pkgs <- c(
  "TCGAbiolinks","SummarizedExperiment","S4Vectors","Biobase","GEOquery",
  "AnnotationDbi","org.Hs.eg.db","hgu133plus2.db","edgeR","GSVA","limma"
)

missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran)) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc)) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

all_pkgs <- c(cran_pkgs, bioc_pkgs)
status <- data.frame(
  package = all_pkgs,
  available = vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    all_pkgs,
    function(p) if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else NA_character_,
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(status)
if (!all(status$available)) {
  stop(
    "Some packages are still missing: ",
    paste(status$package[!status$available], collapse = ", "),
    call. = FALSE
  )
}
message("Environment recovery complete. Now rerun 00_RECOVER_AND_RUN_OLFML2B_FROM_ZERO_v1_0_2.R")
