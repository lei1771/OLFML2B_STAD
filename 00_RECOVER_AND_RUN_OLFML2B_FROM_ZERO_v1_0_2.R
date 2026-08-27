# ==============================================================================
# OLFML2B-STAD TRUE FROM-ZERO RECOVERY + PART0-PART9 RUNNER
# Version: v1.0.2_20260827_GITHUB_PUBLIC_RELEASE_PORTABILITY_REPAIR
#
# One command:
#   source("D:/OLFML2B_STAD/00_RECOVER_AND_RUN_OLFML2B_FROM_ZERO_v1_0_2.R",
#          encoding = "UTF-8", local = FALSE)
#
# Design:
#   - project-local and reproducible
#   - no reuse of deleted output/cache is assumed
#   - downloads public heavy GEO scRNA/spatial archives with resume
#   - TCGA is acquired by Part1; bulk GEO is acquired by Part2
#   - PDC000614 and exact TIGER RDS retain their frozen source semantics
#   - if PDC/TIGER are missing, recovery stops BEFORE formal Part0-Part9 analysis
#   - rerunning the same launcher resumes/skips completed downloads
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
options(timeout = max(getOption("timeout", 60), 7200))

ENTRY <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
ENTRY_DIR <- if (!is.null(ENTRY) && file.exists(ENTRY)) {
  dirname(normalizePath(ENTRY, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

PROJECT_ROOT <- Sys.getenv(
  "OLFML2B_STAD_ROOT",
  unset = ENTRY_DIR
)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = FALSE)
CODE_ROOT <- file.path(PROJECT_ROOT, "R")

DOWNLOAD_HEAVY_GEO <- !identical(toupper(Sys.getenv("OLFML2B_DOWNLOAD_HEAVY", "TRUE")), "FALSE")
FORCE_CLEAN_OUTPUT <- identical(toupper(Sys.getenv("OLFML2B_FORCE_CLEAN_OUTPUT", "FALSE")), "TRUE")
MAKE_FIGURES <- !identical(toupper(Sys.getenv("OLFML2B_MAKE_FIGURES", "TRUE")), "FALSE")

# ------------------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------------------

dir.create(PROJECT_ROOT, recursive = TRUE, showWarnings = FALSE)
RECOVERY_DIR <- file.path(PROJECT_ROOT, "recovery")
AUDIT_DIR <- file.path(RECOVERY_DIR, "audit")
LOG_DIR <- file.path(PROJECT_ROOT, "logs", "RUN_ALL")
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_LOG <- file.path(
  LOG_DIR,
  paste0("RUN_FROM_ZERO_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logm <- function(level = "INFO", ...) {
  line <- sprintf(
    "[%s] [OLFML2B-FROM-ZERO] [%s] %s",
    ts(), toupper(level), paste0(..., collapse = "")
  )
  message(line)
  cat(line, "\n", file = RUN_LOG, append = TRUE, sep = "")
  invisible(line)
}
abort <- function(...) {
  msg <- paste0(..., collapse = "")
  logm("ERROR", msg)
  stop(msg, call. = FALSE)
}
assert <- function(ok, ...) {
  if (!isTRUE(ok)) abort(...)
  invisible(TRUE)
}

logm("INFO", "Launcher version: v1.0.2_20260827_GITHUB_PUBLIC_RELEASE_PORTABILITY_REPAIR")
logm("INFO", "Project root: ", PROJECT_ROOT)
logm("INFO", "Code root: ", CODE_ROOT)

# ------------------------------------------------------------------------------
# 1. Canonical code preflight
# ------------------------------------------------------------------------------

SCRIPTS <- c(
  Part0 = "00_OLFML2B_PART0_CONFIG_CORE.R",
  Part1 = "01_OLFML2B_PART1_TCGA.R",
  Part2 = "02_OLFML2B_PART2_GEO.R",
  Part3 = "03_OLFML2B_PART3_BULK_SURVIVAL.R",
  Part4 = "04_OLFML2B_PART4_IMMUNE_TME_PRODUCTION.R",
  Part5 = "05_OLFML2B_PART5_PDC_PRODUCTION.R",
  Part6 = "06_OLFML2B_PART6_SCRNA_PRODUCTION.R",
  Part7 = "07_OLFML2B_PART7_SPATIAL_TRANSCRIPTOMICS.R",
  Part8 = "08_OLFML2B_PART8_ICI_PRJEB25780_TIGER_ONLY_MOLECULAR_CONTEXT.R",
  Part9 = "09_OLFML2B_PART9_PUBLICATION_OUTPUTS_v8_10_1.R"
)

script_paths <- file.path(CODE_ROOT, unname(SCRIPTS))
missing_scripts <- script_paths[!file.exists(script_paths)]
if (length(missing_scripts)) {
  abort(
    "Canonical R source preflight failed. Missing:\n  - ",
    paste(missing_scripts, collapse = "\n  - ")
  )
}

parse_errors <- character()
for (fp in script_paths) {
  er <- tryCatch({
    parse(file = fp, encoding = "UTF-8")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(er)) parse_errors <- c(parse_errors, paste0(basename(fp), ": ", er))
}
if (length(parse_errors)) {
  abort("R syntax preflight failed:\n  - ", paste(parse_errors, collapse = "\n  - "))
}

code_manifest <- data.frame(
  part = names(SCRIPTS),
  file = unname(SCRIPTS),
  path = normalizePath(script_paths, winslash = "/", mustWork = FALSE),
  size_bytes = file.info(script_paths)$size,
  md5 = unname(tools::md5sum(script_paths)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  code_manifest,
  file.path(AUDIT_DIR, "00_ACTIVE_CODE_MANIFEST.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 2. Frozen-environment preflight
# ------------------------------------------------------------------------------

required_packages <- unique(c(
  "data.table","ggplot2","matrixStats","readxl","R.utils","jsonlite","curl","httr",
  "yaml","survival","metafor","sandwich","lmtest","mice","digest","logistf","Matrix",
  "dplyr","tidyr","readr","forcats","stringr","tibble","purrr","scales","patchwork",
  "ggrepel",
  "TCGAbiolinks","SummarizedExperiment","S4Vectors","Biobase","GEOquery",
  "AnnotationDbi","org.Hs.eg.db","hgu133plus2.db","edgeR","GSVA","limma"
))

pkg_audit <- data.frame(
  package = required_packages,
  available = vapply(required_packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    required_packages,
    function(p) if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else NA_character_,
    character(1)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  pkg_audit,
  file.path(AUDIT_DIR, "01_R_PACKAGE_PREFLIGHT.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (any(!pkg_audit$available)) {
  missing <- pkg_audit$package[!pkg_audit$available]
  writeLines(
    c(
      "The scientific pipeline does not install packages at runtime.",
      "Run this environment-recovery helper once, then rerun the main launcher:",
      file.path(PROJECT_ROOT, "00_INSTALL_REQUIRED_PACKAGES_ONCE.R"),
      "",
      "Missing packages:",
      paste0("- ", missing)
    ),
    file.path(AUDIT_DIR, "01_MISSING_R_PACKAGES.txt")
  )
  abort(
    "Frozen-environment preflight failed. Missing packages: ",
    paste(missing, collapse = ", "),
    ". Run 00_INSTALL_REQUIRED_PACKAGES_ONCE.R once, then rerun this same launcher."
  )
}

# Part9-specific version gates.
if (utils::compareVersion(as.character(utils::packageVersion("ggplot2")), "3.5.0") < 0L) {
  abort("Part9 requires ggplot2 >= 3.5.0.")
}
if (!"free" %in% getNamespaceExports("patchwork")) {
  abort("Part9 requires a patchwork version exporting free().")
}

# ------------------------------------------------------------------------------
# 3. Directory layout
# ------------------------------------------------------------------------------

dirs <- c(
  file.path(PROJECT_ROOT, "data", "raw", "single_cell", "GSE150290"),
  file.path(PROJECT_ROOT, "data", "raw", "single_cell", "GSE183904"),
  file.path(PROJECT_ROOT, "data", "raw", "single_cell", "GSE167297"),
  file.path(PROJECT_ROOT, "data", "raw", "single_cell", "GSE134520"),
  file.path(PROJECT_ROOT, "data", "raw", "spatial", "GSE251950"),
  file.path(PROJECT_ROOT, "data", "PDC_STAD", "PDC000614_standardized", "protein_assembly"),
  file.path(PROJECT_ROOT, "data", "PDC_STAD", "PDC000614_standardized", "metadata"),
  file.path(PROJECT_ROOT, "data", "cache", "ICI_PRJEB25780", "TIGER"),
  file.path(PROJECT_ROOT, "DATA_INBOX", "PDC000614"),
  file.path(PROJECT_ROOT, "DATA_INBOX", "TIGER")
)
invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

if (isTRUE(FORCE_CLEAN_OUTPUT)) {
  for (d in c(file.path(PROJECT_ROOT, "output"), file.path(PROJECT_ROOT, "data", "derived"))) {
    if (dir.exists(d)) unlink(d, recursive = TRUE, force = TRUE)
  }
  logm("WARN", "FORCE_CLEAN_OUTPUT=TRUE: prior output and derived data were removed.")
}

# ------------------------------------------------------------------------------
# 4. Robust resumable downloader for large public GEO supplementary archives
# ------------------------------------------------------------------------------

valid_download <- function(path, min_bytes = 1e6) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size >= min_bytes
}

download_resume <- function(url, dest, min_bytes = 1e6) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (valid_download(dest, min_bytes)) {
    logm("INFO", "Download already present; skipping: ", dest)
    return(invisible(dest))
  }

  part <- paste0(dest, ".part")
  curl_bin <- Sys.which("curl")
  ok <- FALSE

  if (nzchar(curl_bin)) {
    logm("INFO", "Downloading/resuming with curl: ", basename(dest))
    args <- c(
      "-L", "--fail",
      "--retry", "8",
      "--retry-delay", "5",
      "--connect-timeout", "60",
      "--continue-at", "-",
      "--output", shQuote(part),
      shQuote(url)
    )
    if (.Platform$OS.type == "windows") {
      # Prevent the Windows Schannel offline-revocation failure previously
      # observed on this project while retaining HTTPS certificate validation.
      args <- c("--ssl-no-revoke", args)
    }
    status <- suppressWarnings(system2(curl_bin, args = args))
    ok <- identical(as.integer(status), 0L) && valid_download(part, min_bytes)

    # Fallback for curl builds that do not recognize --ssl-no-revoke.
    if (!ok && .Platform$OS.type == "windows") {
      logm("WARN", "First curl attempt failed; retrying without --ssl-no-revoke.")
      args2 <- args[args != "--ssl-no-revoke"]
      status <- suppressWarnings(system2(curl_bin, args = args2))
      ok <- identical(as.integer(status), 0L) && valid_download(part, min_bytes)
    }
  }

  if (!ok) {
    logm("WARN", "curl executable path failed/unavailable; trying R download.file(libcurl).")
    er <- tryCatch({
      utils::download.file(url, part, mode = "wb", method = "libcurl", quiet = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
    ok <- is.null(er) && valid_download(part, min_bytes)
    if (!ok && !is.null(er)) logm("ERROR", "download.file error: ", er)
  }

  if (!ok) {
    abort(
      "Download failed or file is unexpectedly small: ", url,
      "\nPartial file is retained for resume: ", part
    )
  }

  if (file.exists(dest)) unlink(dest, force = TRUE)
  renamed <- file.rename(part, dest)
  if (!renamed) {
    copied <- file.copy(part, dest, overwrite = TRUE)
    if (copied) unlink(part, force = TRUE)
    renamed <- copied
  }
  assert(renamed && valid_download(dest, min_bytes), "Failed to finalize download: ", dest)
  logm("INFO", "Download complete: ", dest, " | bytes=", file.info(dest)$size)
  invisible(dest)
}

geo_raw_url <- function(accession) {
  family <- sub("[0-9]{3}$", "nnn", accession)
  sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s_RAW.tar",
    family, accession, accession
  )
}

heavy_geo <- data.frame(
  accession = c("GSE150290","GSE183904","GSE167297","GSE134520","GSE251950"),
  type = c("single_cell","single_cell","single_cell","single_cell","spatial"),
  min_bytes = c(100e6, 50e6, 10e6, 10e6, 100e6),
  stringsAsFactors = FALSE
)

if (isTRUE(DOWNLOAD_HEAVY_GEO)) {
  for (i in seq_len(nrow(heavy_geo))) {
    acc <- heavy_geo$accession[i]
    parent <- if (heavy_geo$type[i] == "spatial") {
      file.path(PROJECT_ROOT, "data", "raw", "spatial", acc)
    } else {
      file.path(PROJECT_ROOT, "data", "raw", "single_cell", acc)
    }
    dest <- file.path(parent, paste0(acc, "_RAW.tar"))
    download_resume(geo_raw_url(acc), dest, min_bytes = heavy_geo$min_bytes[i])
  }
}

download_audit <- do.call(rbind, lapply(seq_len(nrow(heavy_geo)), function(i) {
  acc <- heavy_geo$accession[i]
  parent <- if (heavy_geo$type[i] == "spatial") {
    file.path(PROJECT_ROOT, "data", "raw", "spatial", acc)
  } else {
    file.path(PROJECT_ROOT, "data", "raw", "single_cell", acc)
  }
  fp <- file.path(parent, paste0(acc, "_RAW.tar"))
  data.frame(
    accession = acc,
    type = heavy_geo$type[i],
    file = fp,
    exists = file.exists(fp),
    size_bytes = if (file.exists(fp)) file.info(fp)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  download_audit,
  file.path(AUDIT_DIR, "02_HEAVY_GEO_RECOVERY_AUDIT.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 5. Unpack spatial GSE251950 now because the final Part7 entry only discovers
#    already unpacked sample directories.
# ------------------------------------------------------------------------------

unpack_gse251950 <- function() {
  raw_dir <- file.path(PROJECT_ROOT, "data", "raw", "spatial", "GSE251950")
  outer <- file.path(raw_dir, "GSE251950_RAW.tar")
  if (!file.exists(outer)) return(data.frame())

  sample_archives <- list.files(
    raw_dir,
    pattern = "^GSM[0-9]+.*\\.tar(\\.gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  rows <- list()

  if (!length(sample_archives)) {
    logm("INFO", "Unpacking GSE251950 outer archive.")
    err <- tryCatch({
      utils::untar(outer, exdir = raw_dir)
      NA_character_
    }, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <- data.frame(
      archive = outer,
      destination = raw_dir,
      status = ifelse(is.na(err), "UNPACKED_OUTER", "FAIL_OUTER"),
      error = err,
      stringsAsFactors = FALSE
    )
  }

  sample_archives <- list.files(
    raw_dir,
    pattern = "^GSM[0-9]+.*\\.tar(\\.gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  for (a in sample_archives) {
    dest <- file.path(raw_dir, sub("\\.tar(\\.gz)?$", "", basename(a), ignore.case = TRUE))
    existing <- dir.exists(dest) && length(list.files(dest, recursive = TRUE, all.files = FALSE)) > 0L
    if (existing) {
      rows[[length(rows) + 1L]] <- data.frame(
        archive = a, destination = dest, status = "SKIPPED_EXISTING", error = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    logm("INFO", "Unpacking spatial section archive: ", basename(a))
    err <- tryCatch({
      utils::untar(a, exdir = dest)
      NA_character_
    }, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <- data.frame(
      archive = a,
      destination = dest,
      status = ifelse(is.na(err), "UNPACKED_SAMPLE", "FAIL_SAMPLE"),
      error = err,
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

spatial_unpack_audit <- unpack_gse251950()
utils::write.csv(
  spatial_unpack_audit,
  file.path(AUDIT_DIR, "03_GSE251950_UNPACK_AUDIT.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 6. Import source-semantic manual inputs from DATA_INBOX / optional env inboxes
# ------------------------------------------------------------------------------

choose_largest <- function(paths) {
  paths <- unique(paths[file.exists(paths) & !dir.exists(paths)])
  if (!length(paths)) return(NA_character_)
  paths[order(file.info(paths)$size, decreasing = TRUE)][1L]
}

find_recursive <- function(root, pattern) {
  if (is.null(root) || !length(root) || is.na(root) || !nzchar(root) || !dir.exists(root)) {
    return(character())
  }
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

pdc_inboxes <- unique(c(
  file.path(PROJECT_ROOT, "DATA_INBOX", "PDC000614"),
  Sys.getenv("OLFML2B_PDC_INBOX", unset = "")
))
pdc_inboxes <- pdc_inboxes[nzchar(pdc_inboxes) & dir.exists(pdc_inboxes)]

pdc_protein_dir <- file.path(
  PROJECT_ROOT, "data", "PDC_STAD", "PDC000614_standardized", "protein_assembly"
)
pdc_metadata_dir <- file.path(
  PROJECT_ROOT, "data", "PDC_STAD", "PDC000614_standardized", "metadata"
)
dir.create(pdc_protein_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdc_metadata_dir, recursive = TRUE, showWarnings = FALSE)

pdc_specs <- list(
  tmt18 = list(pattern = "tmt18\\.tsv$", dest = file.path(pdc_protein_dir, "CPTAC4_Gastric_Cancer_JHU_Proteome.tmt18.tsv"), required = TRUE),
  sample = list(pattern = "sample\\.txt$", dest = file.path(pdc_protein_dir, "CPTAC4_Gastric_Cancer_JHU_Proteome.sample.txt"), required = TRUE),
  label = list(pattern = "label\\.txt$", dest = file.path(pdc_protein_dir, "CPTAC4_Gastric_Cancer_JHU_Proteome.label.txt"), required = FALSE),
  summary = list(pattern = "summary\\.tsv$", dest = file.path(pdc_protein_dir, "CPTAC4_Gastric_Cancer_JHU_Proteome.summary.tsv"), required = FALSE),
  biospecimen = list(pattern = "biospecimen.*\\.tsv$", dest = file.path(pdc_metadata_dir, "PDC000614_biospecimen_latest.tsv"), required = TRUE)
)

pdc_import_rows <- list()
for (nm in names(pdc_specs)) {
  spec <- pdc_specs[[nm]]
  if (!file.exists(spec$dest)) {
    candidates <- unlist(lapply(pdc_inboxes, find_recursive, pattern = spec$pattern), use.names = FALSE)
    src <- choose_largest(candidates)
    if (!is.na(src)) {
      copied <- file.copy(src, spec$dest, overwrite = TRUE)
      if (!copied) abort("Failed to import PDC input: ", src)
      logm("INFO", "PDC import: ", nm, " -> ", spec$dest)
    }
  }
  pdc_import_rows[[nm]] <- data.frame(
    component = nm,
    required = spec$required,
    canonical_file = spec$dest,
    exists = file.exists(spec$dest),
    size_bytes = if (file.exists(spec$dest)) file.info(spec$dest)$size else NA_real_,
    md5 = if (file.exists(spec$dest)) unname(tools::md5sum(spec$dest)) else NA_character_,
    stringsAsFactors = FALSE
  )
}
pdc_audit <- do.call(rbind, pdc_import_rows)
utils::write.csv(
  pdc_audit,
  file.path(AUDIT_DIR, "04_PDC000614_INPUT_AUDIT.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

tiger_dir <- file.path(PROJECT_ROOT, "data", "cache", "ICI_PRJEB25780", "TIGER")
ici_cache_dir <- dirname(tiger_dir)
dir.create(tiger_dir, recursive = TRUE, showWarnings = FALSE)

tiger_inboxes <- unique(c(
  file.path(PROJECT_ROOT, "DATA_INBOX", "TIGER"),
  Sys.getenv("OLFML2B_TIGER_INBOX", unset = "")
))
tiger_inboxes <- tiger_inboxes[nzchar(tiger_inboxes) & dir.exists(tiger_inboxes)]

tiger_exact <- c(
  "STAD-PRJEB25780.Response.Rds",
  "STAD-PRJEB25780.Response (1).Rds"
)

tiger_rows <- list()
for (fn in tiger_exact) {
  dest <- file.path(tiger_dir, fn)
  if (!file.exists(dest)) {
    hits <- unlist(lapply(tiger_inboxes, function(d) {
      list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    }), use.names = FALSE)
    hits <- hits[basename(hits) == fn]
    src <- choose_largest(hits)
    if (!is.na(src)) {
      copied <- file.copy(src, dest, overwrite = TRUE)
      if (!copied) abort("Failed to import TIGER input: ", src)
      logm("INFO", "TIGER import: ", fn, " -> ", dest)
    }
  }
  tiger_rows[[fn]] <- data.frame(
    file = fn,
    canonical_file = dest,
    exists = file.exists(dest),
    size_bytes = if (file.exists(dest)) file.info(dest)$size else NA_real_,
    md5 = if (file.exists(dest)) unname(tools::md5sum(dest)) else NA_character_,
    stringsAsFactors = FALSE
  )
}
tiger_audit <- do.call(rbind, tiger_rows)
utils::write.csv(
  tiger_audit,
  file.path(AUDIT_DIR, "05_TIGER_EXACT_INPUT_AUDIT.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Optional Kim/Bagaev molecular context import. Part8 can run without it.
molecular_names <- c(
  "PRJEB25780_Kim2018_Bagaev2021_clinical_curated.csv",
  "Bagaev2021_Kim_PRJEB25780_clinical_rows_raw.csv",
  "Bagaev2021_mmc5_hits_Gastric.csv"
)
for (fn in molecular_names) {
  dest <- file.path(ici_cache_dir, fn)
  if (file.exists(dest)) next
  hits <- unlist(lapply(tiger_inboxes, function(d) {
    allf <- list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    allf[basename(allf) == fn]
  }), use.names = FALSE)
  src <- choose_largest(hits)
  if (!is.na(src)) {
    file.copy(src, dest, overwrite = TRUE)
    logm("INFO", "Optional molecular-context import: ", fn)
  }
}

# ------------------------------------------------------------------------------
# 7. Recovery readiness gate
# ------------------------------------------------------------------------------

heavy_ready <- all(download_audit$exists & download_audit$size_bytes >= heavy_geo$min_bytes)
pdc_ready <- all(pdc_audit$exists[pdc_audit$required])
tiger_ready <- all(tiger_audit$exists)

readiness <- data.frame(
  component = c(
    "active_code_00_to_09",
    "R_environment",
    "heavy_GEO_scRNA_spatial",
    "PDC000614_required_inputs",
    "TIGER_exact_two_RDS"
  ),
  ready = c(
    TRUE,
    all(pkg_audit$available),
    heavy_ready,
    pdc_ready,
    tiger_ready
  ),
  action_if_false = c(
    "Restore the canonical R bundle.",
    "Run 00_INSTALL_REQUIRED_PACKAGES_ONCE.R, then rerun.",
    "Rerun this same launcher; *.part files are retained for resume.",
    "Download PDC000614 TMT18 + sample map + biospecimen manifest and place them in DATA_INBOX/PDC000614.",
    "Obtain the two exact TIGER RDS files and place them in DATA_INBOX/TIGER."
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  readiness,
  file.path(AUDIT_DIR, "06_FROM_ZERO_READINESS.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

FULL_READY <- all(readiness$ready)

writeLines(
  c(
    paste0("TIMESTAMP=", ts()),
    paste0("FULL_PART0_TO_PART9_LOCAL_PREFLIGHT_READY=", toupper(as.character(FULL_READY))),
    paste0("PDC000614_READY=", toupper(as.character(pdc_ready))),
    paste0("TIGER_READY=", toupper(as.character(tiger_ready))),
    paste0("HEAVY_GEO_READY=", toupper(as.character(heavy_ready)))
  ),
  file.path(AUDIT_DIR, "99_RECOVERY_STATUS.txt")
)

if (!FULL_READY) {
  missing_manual <- character()
  if (!pdc_ready) {
    missing_manual <- c(
      missing_manual,
      paste0(
        "PDC000614: place required files in ",
        file.path(PROJECT_ROOT, "DATA_INBOX", "PDC000614"),
        "\n  required: *tmt18.tsv, *sample.txt, *biospecimen*.tsv",
        "\n  official study: https://pdc.cancer.gov/pdc/study/PDC000614"
      )
    )
  }
  if (!tiger_ready) {
    missing_manual <- c(
      missing_manual,
      paste0(
        "PRJEB25780/TIGER: place BOTH exact files in ",
        file.path(PROJECT_ROOT, "DATA_INBOX", "TIGER"),
        "\n  STAD-PRJEB25780.Response.Rds",
        "\n  STAD-PRJEB25780.Response (1).Rds",
        "\n  TIGER portal: http://tiger.canceromics.org/"
      )
    )
  }

  writeLines(
    c(
      "RECOVERY STAGE COMPLETE, BUT FORMAL PART0-PART9 HAS NOT STARTED.",
      "",
      missing_manual,
      "",
      "After adding the missing files, run the SAME launcher again.",
      "Existing GEO downloads will be skipped/resumed automatically."
    ),
    file.path(AUDIT_DIR, "MANUAL_INPUTS_REQUIRED.txt")
  )

  abort(
    "Recovery is not yet fully ready; formal Part0-Part9 was deliberately NOT started.\n",
    paste(missing_manual, collapse = "\n\n"),
    "\n\nAfter adding the missing files, rerun this same launcher."
  )
}

logm("INFO", "FULL_PART0_TO_PART9_LOCAL_PREFLIGHT_READY=TRUE")

# ------------------------------------------------------------------------------
# 8. Formal Part0-Part9 execution
# ------------------------------------------------------------------------------

Sys.setenv(
  OLFML2B_STAD_ROOT = PROJECT_ROOT,
  OLFML2B_ROOT = PROJECT_ROOT,
  OLFML2B_STAD_CODE_ROOT = CODE_ROOT,
  OLFML2B_TABLES = file.path(PROJECT_ROOT, "output", "tables"),
  OLFML2B_AUTO_INSTALL = "FALSE",
  OLFML2B_PART6_REFRESH_CACHE = "FALSE"
)
setwd(PROJECT_ROOT)

source_part <- function(part) {
  fp <- file.path(CODE_ROOT, SCRIPTS[[part]])
  logm("INFO", "Sourcing ", part, ": ", fp)
  source(fp, encoding = "UTF-8", local = FALSE, chdir = FALSE)
  invisible(fp)
}

run_stage <- function(stage, expr) {
  started <- Sys.time()
  logm("INFO", "========== ", stage, " START ==========")
  out <- tryCatch(
    force(expr),
    error = function(e) {
      logm("ERROR", stage, " FAILED: ", conditionMessage(e))
      stop(e)
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  logm("INFO", "========== ", stage, " COMPLETE | ", sprintf("%.2f", elapsed), " min ==========")
  invisible(out)
}

# Part0
source_part("Part0")
assert(exists("run_part0", mode = "function"), "Part0 entry run_part0() not found.")
ctx <- run_stage(
  "PART0",
  run_part0(
    root = PROJECT_ROOT,
    allow_internet = TRUE,
    overwrite_downloads = FALSE,
    overwrite_results = TRUE,
    run_scrna = FALSE,
    run_cnv = FALSE,
    run_methylation = FALSE,
    auto_install = FALSE,
    save_large_objects = FALSE
  )
)

# Part1: TCGA-STAD acquisition + processing
source_part("Part1")
assert(exists("run_olfml2b_part1_tcga", mode = "function"), "Part1 entry not found.")
res1 <- run_stage("PART1", run_olfml2b_part1_tcga(ctx))

# Part2: GEO bulk acquisition + processing
source_part("Part2")
assert(exists("run_olfml2b_part2_geo", mode = "function"), "Part2 entry not found.")
res2 <- run_stage("PART2", run_olfml2b_part2_geo(ctx = ctx, include_gse84437 = TRUE))

# Part3: bulk survival
source_part("Part3")
assert(exists("run_olfml2b_part3_bulk_survival", mode = "function"), "Part3 entry not found.")
res3 <- run_stage("PART3", run_olfml2b_part3_bulk_survival(ctx = ctx, make_figures = MAKE_FIGURES))

# Part4: immune/TME
source_part("Part4")
assert(exists("run_olfml2b_part4_immune_tme_production", mode = "function"), "Part4 entry not found.")
res4 <- run_stage(
  "PART4",
  run_olfml2b_part4_immune_tme_production(
    root = PROJECT_ROOT,
    geo_validation_cohorts = c("GSE26253","GSE84437","GSE62254","GSE15459"),
    output_subdir = "Part4",
    make_figures = MAKE_FIGURES
  )
)

# Part5: PDC000614
source_part("Part5")
assert(exists("run_olfml2b_part5_pdc_production", mode = "function"), "Part5 entry not found.")
res5 <- run_stage(
  "PART5",
  run_olfml2b_part5_pdc_production(
    root = PROJECT_ROOT,
    pdc614_dir = file.path(PROJECT_ROOT, "data", "PDC_STAD", "PDC000614_standardized"),
    output_subdir = "Part5",
    make_figures = MAKE_FIGURES,
    reset_part5_outputs = TRUE
  )
)

# Part6: four scRNA GEO cohorts; force_unpack is required from RAW.tar clean-room inputs.
source_part("Part6")
assert(exists("run_olfml2b_part6_singlecell_production", mode = "function"), "Part6 entry not found.")
res6 <- run_stage(
  "PART6",
  run_olfml2b_part6_singlecell_production(
    root = PROJECT_ROOT,
    sc_dir = file.path(PROJECT_ROOT, "data", "raw", "single_cell"),
    force_unpack = TRUE,
    make_figures = MAKE_FIGURES,
    refresh_cache = FALSE
  )
)

# Part7: spatial; recovery has already unpacked the outer/sample archives.
source_part("Part7")
assert(exists("run_olfml2b_part7_spatial_transcriptomics", mode = "function"), "Part7 entry not found.")
spatial_dir <- file.path(PROJECT_ROOT, "data", "raw", "spatial", "GSE251950")
if (exists("o2b_p7_unpack_local_archives", mode = "function")) {
  invisible(o2b_p7_unpack_local_archives(spatial_dir, force_unpack = TRUE))
}
res7 <- run_stage(
  "PART7",
  run_olfml2b_part7_spatial_transcriptomics(
    root = PROJECT_ROOT,
    spatial_dir = spatial_dir,
    standardized_dir = file.path(PROJECT_ROOT, "data", "processed", "spatial", "GSE251950_standardized"),
    output_subdir = "Part7",
    make_figures = MAKE_FIGURES,
    bootstrap_B = 2000L,
    moran_permutation_B = 999L,
    distance_permutation_B = 999L,
    reuse_existing_validated_results = FALSE
  )
)

# Part8: restored production TIGER molecular-context implementation.
source_part("Part8")
assert(exists("run_part8_ici_prjeb25780_tiger_only", mode = "function"), "Production Part8 entry not found.")
res8 <- run_stage(
  "PART8",
  run_part8_ici_prjeb25780_tiger_only(
    root = PROJECT_ROOT,
    target_gene = "OLFML2B",
    tiger_dir = tiger_dir,
    molecular_file = NA_character_,
    expression_transform = "auto",
    make_figures = MAKE_FIGURES
  )
)

# Verify exact Part8 -> Part9 handoff before publication rendering.
part8_required <- c(
  "06_TIGER_tumor_only_response_analysis_input.csv",
  "08_ORR_response_wilcoxon_all_features.csv",
  "09_ORR_response_logistic_per_1SD_all_features.csv",
  "23_OLFML2B_module_spearman_correlations.csv",
  "27_final_molecular_context_claim_ceiling.csv",
  "28_part8_primary_dataset_publication_summary.csv"
)
part8_out <- file.path(PROJECT_ROOT, "output", "tables", "Part8_Immunotherapy")
missing_p8 <- part8_required[!file.exists(file.path(part8_out, part8_required))]
if (length(missing_p8)) {
  abort(
    "Part8 completed without the full Part9 handoff contract. Missing: ",
    paste(missing_p8, collapse = ", ")
  )
}

# Part9: source-executed publication renderer.
logm("INFO", "========== PART9 START ==========")
source_part("Part9")
logm("INFO", "========== PART9 COMPLETE ==========")

# ------------------------------------------------------------------------------
# 9. Completion marker
# ------------------------------------------------------------------------------

completion <- c(
  "STATUS=COMPLETE",
  "PIPELINE=OLFML2B_STAD_PART0_TO_PART9",
  "LAUNCHER=v1.0.2_20260827_GITHUB_PUBLIC_RELEASE_PORTABILITY_REPAIR",
  paste0("COMPLETED_AT=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("PROJECT_ROOT=", PROJECT_ROOT),
  "FROM_ZERO_PUBLIC_DATA_RECOVERY=PASS",
  "PDC000614_SOURCE_SEMANTIC_GATE=PASS",
  "TIGER_EXACT_RDS_GATE=PASS",
  "PART8_PRODUCTION_TO_PART9_HANDOFF=PASS"
)
completion_path <- file.path(AUDIT_DIR, "100_FULL_PART0_TO_PART9_COMPLETE.txt")
writeLines(completion, completion_path)

logm("INFO", "============================================================")
logm("INFO", "FULL PART0-PART9 COMPLETE")
logm("INFO", "Completion marker: ", completion_path)
logm("INFO", "============================================================")

invisible(list(
  context = ctx,
  Part1 = res1,
  Part2 = res2,
  Part3 = res3,
  Part4 = res4,
  Part5 = res5,
  Part6 = res6,
  Part7 = res7,
  Part8 = res8,
  completion_marker = completion_path,
  log = RUN_LOG
))
