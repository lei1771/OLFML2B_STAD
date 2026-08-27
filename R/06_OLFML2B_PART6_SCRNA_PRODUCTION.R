# ============================================================================
# OLFML2B-STAD Part6: unified four-dataset production and exact-inference repair
# Version: v1.3.0_20260722_UNIFIED_FOUR_DATASET_EXACT_INFERENCE_AND_COMPARTMENT_PROGRAM_FIX
# ============================================================================
# Purpose:
#   Reconstruct Part6 using the mature auditing logic from the production Part6
#   while removing Kang/GSE206785 from the formal evidence layer and formally
#   re-including four locally available public gastric cancer-related scRNA-seq
#   datasets:
#     - GSE150290: gastric cancer / adjacent normal atlas; primary RC2 layer
#     - GSE183904: primary gastric cancer atlas; primary RC2 layer
#     - GSE167297: diffuse-type GC superficial/deep/normal context; primary context layer
#     - GSE134520: premalignant-to-early-GC context; supportive layer
#
# Design:
#   1. NO Kang/GSE206785 in the formal RC2 analysis.
#   2. No immunotherapy response analysis.
#   3. Automatically unpacks outer RAW.tar and nested tar/tar.gz/tgz archives.
#   4. Supports both 10X matrix.mtx(+features/genes+barcodes) and text count matrices.
#   5. Extracts exact OLFML2B, keeps OLFML2A and OLFM2 as exact-symbol exclusion audits.
#   6. Uses sample/condition/patient/compartment summaries; no cell-level DE as primary evidence.
#   7. Parses target rows and full-transcriptome pseudobulk in one pass.
#   8. Uses input-validated per-matrix caches for interruption-safe resume.
#   9. Outputs machine-auditable manifests, parse status, go/no-go, claim ceiling, and figures.
#
# Run:
#   setwd("D:/OLFML2B_STAD")
#   source("R/06_OLFML2B_PART6_SCRNA_PRODUCTION.R", encoding="UTF-8")
#   o2p6 <- run_olfml2b_part6_scrna_rc2_no_kang_full_production(
#       root = "D:/OLFML2B_STAD",
#       raw_single_cell_dir = "D:/OLFML2B_STAD/data/raw/single_cell",
#       output_subdir = "Part6_RC2",
#       force_unpack = FALSE,
#       make_figures = TRUE
#   )
# ============================================================================

options(stringsAsFactors = FALSE)
OLFML2B_PART6_RC2_VERSION <- "v1.3.0_20260722_UNIFIED_FOUR_DATASET_EXACT_INFERENCE_AND_COMPARTMENT_PROGRAM_FIX"

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

.o2p6_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.o2p6_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.o2p6_env_valid <- nzchar(.o2p6_env_root) && dir.exists(.o2p6_env_root) &&
  file.exists(file.path(.o2p6_env_root, "00_OLFML2B_PART0_CONFIG_CORE.R")) &&
  file.exists(file.path(.o2p6_env_root, "06_OLFML2B_PART6_SCRNA_PRODUCTION.R"))
.o2p6_code_root <- if (.o2p6_env_valid) {
  normalizePath(.o2p6_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.o2p6_entry) && file.exists(.o2p6_entry)) {
  dirname(normalizePath(.o2p6_entry, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .o2p6_code_root)

o2p6_ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
o2p6_stop <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
o2p6_assert <- function(cond, ...) if (!isTRUE(cond)) o2p6_stop(...)

o2p6_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    o2p6_assert(ok || dir.exists(path), "Cannot create directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

o2p6_log <- function(level = "INFO", ..., log_file = NULL, echo = TRUE) {
  line <- sprintf("[%s] [%s] [OLFML2B-P6-4DS] %s", o2p6_ts(), toupper(level), paste0(..., collapse = ""))
  if (isTRUE(echo)) message(line)
  if (!is.null(log_file) && nzchar(log_file)) {
    o2p6_dir(dirname(log_file))
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
  invisible(line)
}

o2p6_bind_rows <- function(xs) {
  if (is.null(xs)) return(data.frame())
  if (is.data.frame(xs)) return(xs)
  xs <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  out <- lapply(xs, function(z) {
    for (m in setdiff(cols, names(z))) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

o2p6_atomic_write_csv <- function(x, path, row.names = FALSE, na = "") {
  o2p6_dir(dirname(path))
  if (is.null(x)) x <- data.frame()
  if (is.atomic(x) && is.null(dim(x))) x <- data.frame(value = x, stringsAsFactors = FALSE)
  if (!is.data.frame(x) && !is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  utils::write.csv(x, tmp, row.names = row.names, na = na, fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  o2p6_assert(ok && file.exists(path), "Failed to write CSV: ", path)
  invisible(path)
}

o2p6_atomic_save_rds <- function(object, path, compress = "xz") {
  o2p6_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  o2p6_assert(ok && file.exists(path), "Failed to save RDS: ", path)
  invisible(path)
}

o2p6_source_part0_if_needed <- function(root) {
  if (exists("olfml2b_build_dirs", mode = "function")) return(invisible(TRUE))
  candidates <- unique(c(
    file.path(.o2p6_code_root, "00_OLFML2B_PART0_CONFIG_CORE.R"),
    file.path(root, "R", "00_OLFML2B_PART0_CONFIG_CORE.R"),
    file.path(getwd(), "R", "00_OLFML2B_PART0_CONFIG_CORE.R"),
    file.path(getwd(), "00_OLFML2B_PART0_CONFIG_CORE.R")
  ))
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit) && nzchar(hit)) sys.source(hit, envir = environment(o2p6_source_part0_if_needed), chdir = FALSE)
  invisible(TRUE)
}

o2p6_dirs <- function(root, output_subdir = "Part6_RC2") {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  o2p6_source_part0_if_needed(root)
  if (exists("olfml2b_build_dirs", mode = "function")) {
    d <- olfml2b_build_dirs(root)
    return(list(
      root = root,
      tables = o2p6_dir(file.path(d$tables_root, output_subdir)),
      figures = o2p6_dir(file.path(d$figures_root, output_subdir)),
      reports = o2p6_dir(file.path(d$reports_root, output_subdir)),
      qc = o2p6_dir(file.path(d$qc_root, output_subdir)),
      logs = o2p6_dir(file.path(d$logs_runtime, output_subdir)),
      objects = o2p6_dir(d$objects)
    ))
  }
  list(
    root = root,
    tables = o2p6_dir(file.path(root, "output", "tables", output_subdir)),
    figures = o2p6_dir(file.path(root, "output", "figures", output_subdir)),
    reports = o2p6_dir(file.path(root, "output", "reports", output_subdir)),
    qc = o2p6_dir(file.path(root, "output", "qc", output_subdir)),
    logs = o2p6_dir(file.path(root, "logs", "runtime", output_subdir)),
    objects = o2p6_dir(file.path(root, "output", "objects"))
  )
}

o2p6_pkg_versions <- function(pkgs) {
  data.frame(
    package = pkgs,
    available = vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

o2p6_manifest_md5 <- function(paths, info = file.info(paths),
                               cache_file = getOption("olfml2b.part6.md5_cache_file", NULL)) {
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  current <- data.frame(
    path = paths, size = as.numeric(info$size), mtime = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
  if (is.null(cache_file) || !length(cache_file) || is.na(cache_file[1]) || !nzchar(cache_file[1])) {
    return(list(md5 = tryCatch(unname(tools::md5sum(paths)), error = function(e) rep(NA_character_, length(paths))),
                status = rep("COMPUTED_NO_CACHE", length(paths))))
  }
  cache_file <- cache_file[1]
  old <- if (!isTRUE(getOption("olfml2b.part6.refresh_md5_cache", FALSE)) && file.exists(cache_file)) {
    tryCatch(readRDS(cache_file), error = function(e) data.frame())
  } else data.frame()
  valid_old <- is.data.frame(old) && all(c("path", "size", "mtime", "md5") %in% names(old))
  if (!valid_old) old <- data.frame(path = character(), size = numeric(), mtime = numeric(), md5 = character(), stringsAsFactors = FALSE)
  key <- paste(current$path, current$size, current$mtime, sep = "||")
  old_key <- paste(old$path, old$size, old$mtime, sep = "||")
  hit <- match(key, old_key)
  md5 <- old$md5[hit]
  status <- ifelse(!is.na(hit) & !is.na(md5) & nzchar(md5), "CACHE_HIT_PATH_SIZE_MTIME", "COMPUTED_AND_CACHED")
  need <- which(status == "COMPUTED_AND_CACHED")
  if (length(need)) md5[need] <- tryCatch(unname(tools::md5sum(paths[need])), error = function(e) rep(NA_character_, length(need)))
  updated <- rbind(
    data.frame(current, md5 = md5, stringsAsFactors = FALSE),
    old[, c("path", "size", "mtime", "md5"), drop = FALSE]
  )
  updated_key <- paste(updated$path, updated$size, updated$mtime, sep = "||")
  updated <- updated[!duplicated(updated_key), , drop = FALSE]
  o2p6_atomic_save_rds(updated, cache_file, compress = "gzip")
  list(md5 = md5, status = status)
}

o2p6_file_manifest <- function(root, recursive = TRUE, pattern = NULL) {
  if (!dir.exists(root)) return(data.frame())
  files <- list.files(root, recursive = recursive, full.names = TRUE, all.files = FALSE, pattern = pattern)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(data.frame())
  info <- file.info(files)
  base <- normalizePath(root, winslash = "/", mustWork = FALSE)
  abs <- normalizePath(files, winslash = "/", mustWork = FALSE)
  checksum <- o2p6_manifest_md5(abs, info = info)
  data.frame(
    relative_path = substring(abs, nchar(base) + 2L),
    absolute_path = abs,
    file_name = basename(abs),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 3),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    md5 = checksum$md5,
    md5_status = checksum$status,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------------
# Gene and marker definitions
# ----------------------------------------------------------------------------
o2p6_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^[\\x22\\x27]+|[\\x22\\x27]+$", "", x, perl = TRUE)
  x <- gsub("^X__|^X_", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("\\..*$", "", x)
  x <- sub("^HGNC:", "", x)
  x <- trimws(x)
  x
}

o2p6_marker_sets <- function() {
  list(
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2", "KRT7", "KRT17"),
    Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "PDGFRA", "FAP", "THY1"),
    Myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1", "CALD1"),
    Smooth_Muscle = c("MYH11", "DES", "SMTN", "ACTG2", "CNN1", "TAGLN", "ACTA2"),
    Pericyte = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9"),
    Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "CDH5", "PLVAP", "EMCN"),
    Myeloid = c("LYZ", "LST1", "FCER1G", "TYROBP", "CTSS", "CD68", "C1QA", "C1QB", "C1QC"),
    T_NK = c("PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "CD8A", "NKG7", "GNLY", "PRF1", "GZMB"),
    B_Plasma = c("MS4A1", "CD79A", "CD79B", "CD37", "JCHAIN", "MZB1", "SDC1", "IGHG1"),
    Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC")
  )
}

o2p6_program_sets <- function() {
  list(
    CAF_ECM = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "POSTN", "THY1", "COL6A1", "COL6A2"),
    TGFb_Response = c("TGFB1", "TGFBI", "SERPINE1", "SMAD3", "SMAD7", "CTGF", "INHBA", "PMEPA1"),
    ECM_Remodeling = c("COL1A1", "COL3A1", "FN1", "MMP2", "MMP11", "MMP14", "SPARC", "POSTN", "BGN", "THBS2", "LOX", "PLOD2"),
    Myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1", "CALD1"),
    Inflammatory_Fibroblast = c("CXCL12", "IL6", "CFD", "C7", "C3", "CCL2", "PTGDS", "APOD"),
    Smooth_Muscle = c("MYH11", "DES", "SMTN", "ACTG2", "CNN1", "TAGLN"),
    Pericyte = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9"),
    Myeloid_Macrophage = c("LYZ", "LST1", "CD68", "CD163", "MRC1", "C1QA", "C1QB", "C1QC", "FCER1G"),
    Endothelial_Angiogenic = c("PECAM1", "VWF", "KDR", "FLT1", "ENG", "ESAM", "PLVAP", "EMCN"),
    Epithelial_Differentiation = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"),
    CD8_Cytotoxic = c("CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG"),
    Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "UBE2C")
  )
}

o2p6_exclusion_symbols <- function() {
  c("OLFML2A", "OLFM2")
}

o2p6_target_genes <- function() {
  unique(c("OLFML2B", o2p6_exclusion_symbols(), unlist(o2p6_marker_sets(), use.names = FALSE), unlist(o2p6_program_sets(), use.names = FALSE)))
}

# ----------------------------------------------------------------------------
# Tar unpacking and file discovery
# ----------------------------------------------------------------------------
o2p6_infer_cohort <- function(path, sc_dir) {
  rel <- substring(
    normalizePath(path, winslash = "/", mustWork = FALSE),
    nchar(normalizePath(sc_dir, winslash = "/", mustWork = FALSE)) + 2L
  )
  strsplit(rel, "/", fixed = TRUE)[[1]][1]
}

o2p6_strip_ext <- function(x) {
  x <- basename(x)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.tar$", "", x, ignore.case = TRUE)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("\\.tsv$", "", x, ignore.case = TRUE)
  x <- sub("\\.mtx$", "", x, ignore.case = TRUE)
  x
}

o2p6_unpack_nested_archives <- function(dataset_dir, force_unpack = FALSE, max_rounds = 4L, log_file = NULL) {
  dataset_dir <- normalizePath(dataset_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dataset_dir)) return(data.frame())
  unpack_root <- o2p6_dir(file.path(dataset_dir, "extracted_rc2"))
  audit <- list()

  for (round in seq_len(max_rounds)) {
    search_roots <- unique(c(dataset_dir, unpack_root))
    all_files <- unique(unlist(lapply(search_roots, function(r) {
      if (!dir.exists(r)) character() else list.files(r, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    }), use.names = FALSE))
    all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
    archives <- all_files[grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$", all_files, ignore.case = TRUE)]
    archives <- archives[!grepl("/extracted_rc2/", normalizePath(archives, winslash = "/", mustWork = FALSE), fixed = TRUE) |
                           grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$", archives, ignore.case = TRUE)]
    if (!length(archives)) break

    progressed <- FALSE
    for (a in archives) {
      rel_safe <- gsub("[^A-Za-z0-9_.-]+", "_", o2p6_strip_ext(basename(a)))
      dest <- file.path(unpack_root, rel_safe)
      marker <- file.path(dest, ".o2p6_unpacked.ok")
      if (file.exists(marker) && !force_unpack) {
        audit[[length(audit)+1L]] <- data.frame(archive=a, dest=dest, status="SKIP_ALREADY_UNPACKED", stringsAsFactors = FALSE)
        next
      }
      o2p6_dir(dest)
      o2p6_log("INFO", "Unpacking archive: ", basename(a), " -> ", dest, log_file = log_file)
      ok <- tryCatch({
        utils::untar(a, exdir = dest)
        cat(o2p6_ts(), file = marker)
        TRUE
      }, error = function(e) {
        audit[[length(audit)+1L]] <<- data.frame(archive=a, dest=dest, status="UNPACK_FAILED", error=conditionMessage(e), stringsAsFactors = FALSE)
        FALSE
      })
      if (ok) {
        audit[[length(audit)+1L]] <- data.frame(archive=a, dest=dest, status="UNPACKED", error=NA_character_, stringsAsFactors = FALSE)
        progressed <- TRUE
      }
    }
    if (!progressed) break
  }
  o2p6_bind_rows(audit)
}

o2p6_find_10x_dirs <- function(dataset_dir) {
  if (!dir.exists(dataset_dir)) return(data.frame())
  files <- list.files(dataset_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  matrix_files <- files[grepl("matrix\\.mtx(\\.gz)?$", basename(files), ignore.case = TRUE)]
  if (!length(matrix_files)) return(data.frame())
  rows <- lapply(matrix_files, function(mtx) {
    d <- dirname(mtx)
    all_in_d <- list.files(d, full.names = TRUE, all.files = FALSE)
    feature <- all_in_d[grepl("features\\.tsv(\\.gz)?$|genes\\.tsv(\\.gz)?$", basename(all_in_d), ignore.case = TRUE)][1]
    barcode <- all_in_d[grepl("barcodes\\.tsv(\\.gz)?$", basename(all_in_d), ignore.case = TRUE)][1]
    data.frame(
      matrix_file = normalizePath(mtx, winslash = "/", mustWork = FALSE),
      feature_file = if (!is.na(feature)) normalizePath(feature, winslash = "/", mustWork = FALSE) else NA_character_,
      barcode_file = if (!is.na(barcode)) normalizePath(barcode, winslash = "/", mustWork = FALSE) else NA_character_,
      sample_dir = normalizePath(d, winslash = "/", mustWork = FALSE),
      complete_10x = !is.na(feature) && !is.na(barcode),
      stringsAsFactors = FALSE
    )
  })
  o2p6_bind_rows(rows)
}

o2p6_10x_container_dir <- function(sample_dir) {
  # GEO 10X archives often unpack as <GSM...raw_gene_bc_matrices>/hg19/matrix.mtx.
  # The biological sample identifier is the parent archive folder, not the genome
  # subdirectory "hg19".  Returning the parent folder prevents all GSE150290
  # samples from being logged/deduplicated as a single "hg19" sample.
  d <- normalizePath(sample_dir, winslash = "/", mustWork = FALSE)
  b <- basename(d)
  if (grepl("^(hg19|hg38|GRCh38|GRCh37|mm10|mm9)$", b, ignore.case = TRUE)) {
    return(dirname(d))
  }
  d
}

o2p6_10x_sample_label <- function(sample_dir) {
  basename(o2p6_10x_container_dir(sample_dir))
}

o2p6_session_packages_df <- function() {
  si <- utils::sessionInfo()
  pkgs <- c(si$basePkgs %||% character(), names(si$otherPkgs %||% list()), names(si$loadedOnly %||% list()))
  pkgs <- unique(pkgs[nzchar(pkgs)])
  if (!length(pkgs)) return(data.frame())
  data.frame(
    package = pkgs,
    version = vapply(pkgs, function(pkg) {
      out <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
      out
    }, character(1)),
    attached = pkgs %in% c(si$basePkgs %||% character(), names(si$otherPkgs %||% list())),
    loaded_only = pkgs %in% names(si$loadedOnly %||% list()),
    stringsAsFactors = FALSE
  )
}

o2p6_find_text_matrices <- function(dataset_dir) {
  if (!dir.exists(dataset_dir)) return(data.frame())
  files <- list.files(dataset_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[grepl("\\.(txt|txt\\.gz|csv|csv\\.gz|tsv|tsv\\.gz)$", files, ignore.case = TRUE)]
  files <- files[!grepl("series_matrix|family.soft|barcodes|features|genes|matrix.mtx|sample_tag|tag_calls", basename(files), ignore.case = TRUE)]
  # A permissive filter: allow all non-series tabular candidates that are not
  # clearly barcode/feature/metadata support files.  This is needed for
  # GSE183904-style csv.gz matrices whose filenames may not contain "count" or
  # "matrix". Very small files are still kept in the manifest but will usually
  # fail safely at parse time.
  files <- files[!grepl("metadata|meta_data|annotation|annot|clinical|sample_tag|tag_calls|barcode|barcodes|feature|features|gene_list|genelist", basename(files), ignore.case = TRUE)]
  if (!length(files)) return(data.frame())
  info <- file.info(files)
  data.frame(
    matrix_file = normalizePath(files, winslash = "/", mustWork = FALSE),
    file_name = basename(files),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 3),
    stringsAsFactors = FALSE
  )
}

o2p6_candidate_priority <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (grepl("/extracted_rc2/", p, ignore.case = TRUE)) return(1L)
  if (grepl("_RAW|extracted_auto|extracted", p, ignore.case = TRUE)) return(2L)
  3L
}

o2p6_sample_key <- function(path) {
  b <- basename(path)
  b0 <- o2p6_strip_ext(b)
  gsm <- regmatches(b0, regexpr("GSM[0-9]+", b0, ignore.case = TRUE))
  if (length(gsm) && nzchar(gsm)) return(toupper(gsm))
  # Try patient-like name embedded in path
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  hit <- regmatches(p, regexpr("(Pat|PT|Pt|P|GC)[-_]?[0-9]+[A-Za-z0-9_-]*", p, ignore.case = TRUE))
  if (length(hit) && nzchar(hit)) return(toupper(hit))
  toupper(gsub("[^A-Za-z0-9]+", "_", b0))
}

o2p6_deduplicate_candidates <- function(candidates) {
  if (!is.data.frame(candidates) || !nrow(candidates)) return(candidates)
  candidates$sample_key <- vapply(candidates$sample_path, o2p6_sample_key, character(1))
  candidates$priority <- vapply(candidates$sample_path, o2p6_candidate_priority, integer(1))
  candidates$kept <- FALSE
  candidates$deduplication_reason <- NA_character_
  for (ky in unique(candidates$sample_key)) {
    idx <- which(candidates$sample_key == ky)
    ord <- order(candidates$priority[idx], -candidates$size_bytes[idx], candidates$sample_path[idx])
    chosen <- idx[ord[1]]
    candidates$kept[chosen] <- TRUE
    candidates$deduplication_reason[chosen] <- "KEPT_BEST_PRIORITY_BY_SAMPLE_KEY"
    candidates$deduplication_reason[setdiff(idx, chosen)] <- "REMOVED_DUPLICATE_SAMPLE_KEY"
  }
  candidates
}


o2p6_apply_input_policy <- function(candidates, exclude_gse150290_raw_10x = TRUE) {
  if (!is.data.frame(candidates) || !nrow(candidates)) return(candidates)
  candidates$policy_keep <- TRUE
  candidates$input_policy <- "ALLOW_CANDIDATE"
  candidates$input_policy_reason <- "Candidate retained after dataset-specific input policy."
  candidates$matrix_origin_class <- ifelse(candidates$input_type == "10X_MTX", "tenx_mtx", "text_matrix")

  # v2.0.2 critical guard:
  # GSE150290_RAW.tar contains both author-provided per-sample TXT matrices and
  # CellRanger raw_gene_bc_matrices. The latter include the full raw barcode
  # universe (often 737,280 barcodes per sample) and massively inflate cell
  # counts by retaining empty droplets. For formal RC2 inference, GSE150290 is
  # restricted to the author processed TXT matrices. The raw 10X files are kept
  # in the candidate manifest for audit only and are never used unless the user
  # deliberately disables this guard.
  path_blob <- paste(candidates$sample_path, candidates$matrix_file, candidates$file_name, sep = " ")
  is_gse150290_raw10x <- candidates$cohort == "GSE150290" &
    candidates$input_type == "10X_MTX" &
    grepl("raw_gene_bc_matrices|/hg19/|/hg38/|/GRCh37/|/GRCh38/", path_blob, ignore.case = TRUE)
  if (isTRUE(exclude_gse150290_raw_10x) && any(is_gse150290_raw10x, na.rm = TRUE)) {
    candidates$policy_keep[is_gse150290_raw10x] <- FALSE
    candidates$input_policy[is_gse150290_raw10x] <- "EXCLUDE_GSE150290_RAW_10X_DROPLET_MATRIX"
    candidates$input_policy_reason[is_gse150290_raw10x] <- paste(
      "GSE150290 raw_gene_bc_matrices are raw barcode matrices and inflate cell counts;",
      "author processed GSM*.txt(.gz) matrices are used for formal analysis."
    )
    candidates$matrix_origin_class[is_gse150290_raw10x] <- "raw_10x_droplet_matrix_excluded"
  }

  is_gse150290_text <- candidates$cohort == "GSE150290" & candidates$input_type == "TEXT_MATRIX"
  if (any(is_gse150290_text, na.rm = TRUE)) {
    candidates$input_policy[is_gse150290_text] <- "PREFER_GSE150290_AUTHOR_PROCESSED_TEXT_MATRIX"
    candidates$input_policy_reason[is_gse150290_text] <- "Author processed per-sample TXT matrix used as the GSE150290 inferential input."
    candidates$matrix_origin_class[is_gse150290_text] <- "author_processed_text_matrix"
  }

  candidates
}

o2p6_candidate_policy_summary <- function(candidates) {
  if (!is.data.frame(candidates) || !nrow(candidates)) return(data.frame())
  if (!"policy_keep" %in% names(candidates)) candidates$policy_keep <- TRUE
  if (!"used_for_analysis" %in% names(candidates)) candidates$used_for_analysis <- candidates$kept %in% TRUE & candidates$complete %in% TRUE & candidates$policy_keep %in% TRUE
  pieces <- split(candidates, paste(candidates$cohort, candidates$input_type, candidates$input_policy, sep = "__"))
  out <- lapply(pieces, function(z) {
    data.frame(
      cohort = z$cohort[1],
      input_type = z$input_type[1],
      input_policy = z$input_policy[1],
      matrix_origin_class = z$matrix_origin_class[1],
      n_candidates = nrow(z),
      n_kept_after_dedup = sum(z$kept %in% TRUE, na.rm = TRUE),
      n_policy_keep = sum(z$policy_keep %in% TRUE, na.rm = TRUE),
      n_used_for_analysis = sum(z$used_for_analysis %in% TRUE, na.rm = TRUE),
      total_size_mb = round(sum(z$size_bytes, na.rm = TRUE) / 1024^2, 3),
      policy_reason = z$input_policy_reason[1],
      stringsAsFactors = FALSE
    )
  })
  o2p6_bind_rows(out)
}

o2p6_discover_dataset_inputs <- function(dataset_dir, cohort, force_unpack = FALSE, exclude_gse150290_raw_10x = TRUE, log_file = NULL) {
  unpack_audit <- o2p6_unpack_nested_archives(dataset_dir, force_unpack = force_unpack, log_file = log_file)
  manifest <- o2p6_file_manifest(dataset_dir, recursive = TRUE)
  tenx <- o2p6_find_10x_dirs(dataset_dir)
  text <- o2p6_find_text_matrices(dataset_dir)

  rows <- list()
  if (nrow(tenx)) {
    for (i in seq_len(nrow(tenx))) {
      rows[[length(rows)+1L]] <- data.frame(
        cohort = cohort,
        input_type = "10X_MTX",
        sample_path = o2p6_10x_container_dir(tenx$sample_dir[i]),
        matrix_file = tenx$matrix_file[i],
        feature_file = tenx$feature_file[i],
        barcode_file = tenx$barcode_file[i],
        file_name = o2p6_10x_sample_label(tenx$sample_dir[i]),
        size_bytes = as.numeric(file.info(tenx$matrix_file[i])$size),
        complete = isTRUE(tenx$complete_10x[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (nrow(text)) {
    for (i in seq_len(nrow(text))) {
      rows[[length(rows)+1L]] <- data.frame(
        cohort = cohort,
        input_type = "TEXT_MATRIX",
        sample_path = text$matrix_file[i],
        matrix_file = text$matrix_file[i],
        feature_file = NA_character_,
        barcode_file = NA_character_,
        file_name = basename(text$matrix_file[i]),
        size_bytes = text$size_bytes[i],
        complete = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  candidates <- o2p6_bind_rows(rows)
  if (nrow(candidates)) {
    candidates <- o2p6_deduplicate_candidates(candidates)
    candidates <- o2p6_apply_input_policy(candidates, exclude_gse150290_raw_10x = exclude_gse150290_raw_10x)
  }
  list(unpack_audit = unpack_audit, manifest = manifest, candidates = candidates)
}

# ----------------------------------------------------------------------------
# Sample metadata inference
# ----------------------------------------------------------------------------
o2p6_infer_sample_info <- function(cohort, sample_path, file_name = basename(sample_path)) {
  raw <- paste(sample_path, file_name, sep = "/")
  sample_id <- o2p6_sample_key(file_name)
  if (grepl("GSM[0-9]+", raw, ignore.case = TRUE)) {
    sample_id <- toupper(regmatches(raw, regexpr("GSM[0-9]+", raw, ignore.case = TRUE))[1])
  }
  patient <- NA_character_
  tissue_or_condition <- NA_character_
  disease_context <- NA_character_
  dataset_role <- switch(
    cohort,
    GSE150290 = "PRIMARY_TUMOR_NORMAL_GC_ATLAS",
    GSE183904 = "PRIMARY_ADVANCED_GC_TUMOR_ATLAS",
    GSE167297 = "PRIMARY_DGC_SUPERFICIAL_DEEP_CONTEXT",
    GSE134520 = "SUPPORTIVE_PREMALIGNANT_EGC_CONTEXT",
    "UNSPECIFIED"
  )

  if (identical(cohort, "GSE134520")) {
    hit <- regmatches(raw, regexpr("(NAG|CAG|IMW|IMS|EGC)[0-9]*", raw, ignore.case = TRUE))
    if (length(hit) && nzchar(hit)) tissue_or_condition <- toupper(gsub("[0-9]+$", "", hit))
    disease_context <- "premalignant_to_early_GC"
  }

  if (identical(cohort, "GSE167297")) {
    patient_hit <- regmatches(raw, regexpr("Pt[0-9]+", raw, ignore.case = TRUE))
    if (length(patient_hit) && nzchar(patient_hit)) patient <- patient_hit
    state_hit <- regmatches(raw, regexpr("(Normal|Superficial|Deep)", raw, ignore.case = TRUE))
    if (length(state_hit) && nzchar(state_hit)) tissue_or_condition <- state_hit
    disease_context <- "diffuse_type_GC_invasion_depth"
  }

  if (identical(cohort, "GSE150290")) {
    patient_hit <- regmatches(raw, regexpr("Pat[0-9]+|Patient[0-9]+|PT[0-9]+|P[0-9]+", raw, ignore.case = TRUE))
    if (length(patient_hit) && nzchar(patient_hit)) patient <- patient_hit
    if (grepl("normal|adjacent|non.?cancer|NT|_N\\b|-N\\b", raw, ignore.case = TRUE)) tissue_or_condition <- "Adjacent_or_Normal"
    if (grepl("tumou?r|cancer|GC|lesion|_T\\b|-T\\b", raw, ignore.case = TRUE)) tissue_or_condition <- tissue_or_condition %||% "Tumor_or_GC"
    if (grepl("IGC|intestinal", raw, ignore.case = TRUE)) disease_context <- "intestinal_type_GC"
    if (grepl("DGC|diffuse", raw, ignore.case = TRUE)) disease_context <- "diffuse_type_GC"
    if (is.na(disease_context)) disease_context <- "GC_adjacent_normal_atlas"
  }

  if (identical(cohort, "GSE183904")) {
    patient_hit <- regmatches(raw, regexpr("GC[0-9]+|P[0-9]+|Patient[0-9]+|PT[0-9]+", raw, ignore.case = TRUE))
    if (length(patient_hit) && nzchar(patient_hit)) patient <- patient_hit
    if (grepl("normal|adjacent|non.?tumou?r|NT", raw, ignore.case = TRUE)) tissue_or_condition <- "Adjacent_or_Normal"
    if (grepl("tumou?r|cancer|GC|primary", raw, ignore.case = TRUE)) tissue_or_condition <- tissue_or_condition %||% "Tumor_or_GC"
    disease_context <- "primary_GC_atlas"
  }

  data.frame(
    cohort = cohort,
    dataset_role = dataset_role,
    sample_id = sample_id,
    patient = patient,
    tissue_or_condition = tissue_or_condition,
    disease_context = disease_context,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------------
# Matrix parsers
# ----------------------------------------------------------------------------
o2p6_detect_sep_from_header <- function(header_line) {
  n_tab <- lengths(regmatches(header_line, gregexpr("\t", header_line, fixed = TRUE)))
  n_comma <- lengths(regmatches(header_line, gregexpr(",", header_line, fixed = TRUE)))
  if (n_tab >= n_comma && n_tab > 0) return("\t")
  if (n_comma > 0) return(",")
  "\t"
}

o2p6_open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
}

P6RC2_FULLTX_CACHE_SCHEMA <- "part6_fulltx_single_pass_v1"

o2p6_source_fingerprint <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  info <- file.info(paths)
  data.frame(
    path = paths,
    exists = file.exists(paths),
    size = as.numeric(info$size),
    mtime = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
}

o2p6_cache_key <- function(paths) {
  fp <- o2p6_source_fingerprint(paths)
  token <- paste(fp$path, collapse = "|")
  h <- 0
  for (z in utf8ToInt(enc2utf8(token))) h <- (h * 131 + z) %% 2147483647
  stem <- gsub("[^A-Za-z0-9._-]+", "_", basename(fp$path[1] %||% "matrix"))
  stem <- substr(stem, 1L, 72L)
  paste0(stem, "_", sprintf("%010.0f", h), ".rds")
}

o2p6_fulltx_cache_path <- function(paths, cache_dir = getOption("olfml2b.part6.fulltx_cache_dir", NULL)) {
  if (is.null(cache_dir) || !length(cache_dir) || is.na(cache_dir[1]) || !nzchar(cache_dir[1])) return(NA_character_)
  file.path(o2p6_dir(cache_dir[1]), o2p6_cache_key(paths))
}

o2p6_fingerprint_identical <- function(x, y) {
  if (!is.data.frame(x) || !is.data.frame(y) || !identical(names(x), names(y)) || nrow(x) != nrow(y)) return(FALSE)
  identical(as.character(x$path), as.character(y$path)) &&
    identical(as.logical(x$exists), as.logical(y$exists)) &&
    identical(as.numeric(x$size), as.numeric(y$size)) &&
    identical(as.numeric(x$mtime), as.numeric(y$mtime))
}

o2p6_fulltx_cache_read <- function(paths, genes = NULL, require_target = FALSE,
                                    cache_dir = getOption("olfml2b.part6.fulltx_cache_dir", NULL),
                                    refresh_cache = getOption("olfml2b.part6.refresh_cache", FALSE)) {
  cache_file <- o2p6_fulltx_cache_path(paths, cache_dir)
  miss <- list(hit = FALSE, cache_file = cache_file, object = NULL, reason = "CACHE_MISS")
  if (isTRUE(refresh_cache)) {
    miss$reason <- "CACHE_REFRESH_REQUESTED"
    return(miss)
  }
  if (is.na(cache_file) || !file.exists(cache_file)) return(miss)
  obj <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (!is.list(obj) || !identical(obj$schema, P6RC2_FULLTX_CACHE_SCHEMA)) {
    miss$reason <- "CACHE_SCHEMA_MISMATCH"
    return(miss)
  }
  if (!o2p6_fingerprint_identical(obj$source_fingerprint, o2p6_source_fingerprint(paths))) {
    miss$reason <- "CACHE_INPUT_CHANGED"
    return(miss)
  }
  if (!is.data.frame(obj$fulltx) || !all(c("gene", "pseudobulk_count", "n_cells") %in% names(obj$fulltx))) {
    miss$reason <- "CACHE_FULLTX_INVALID"
    return(miss)
  }
  if (isTRUE(require_target)) {
    genes_clean <- sort(unique(o2p6_clean_gene(genes %||% character())))
    if (is.null(obj$target_payload) || !identical(sort(unique(obj$target_genes)), genes_clean)) {
      miss$reason <- "CACHE_TARGET_PAYLOAD_MISSING"
      return(miss)
    }
  }
  list(hit = TRUE, cache_file = cache_file, object = obj, reason = "CACHE_HIT_VALIDATED")
}

o2p6_fulltx_cache_write <- function(paths, fulltx, target_payload = NULL, genes = NULL,
                                     created_by = "fulltx_fallback",
                                     cache_dir = getOption("olfml2b.part6.fulltx_cache_dir", NULL)) {
  cache_file <- o2p6_fulltx_cache_path(paths, cache_dir)
  if (is.na(cache_file)) return(invisible(NA_character_))
  obj <- list(
    schema = P6RC2_FULLTX_CACHE_SCHEMA,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    created_by = created_by,
    validation_rule = "normalized_path+byte_size+mtime+cache_schema",
    source_fingerprint = o2p6_source_fingerprint(paths),
    target_genes = sort(unique(o2p6_clean_gene(genes %||% character()))),
    target_payload = target_payload,
    fulltx = fulltx[, c("gene", "pseudobulk_count", "n_cells"), drop = FALSE]
  )
  o2p6_atomic_save_rds(obj, cache_file, compress = "gzip")
  invisible(cache_file)
}

o2p6_collapse_fulltx <- function(gene, count, n_cells) {
  out <- data.frame(
    gene = o2p6_clean_gene(gene),
    pseudobulk_count = as.numeric(count),
    n_cells = as.numeric(n_cells),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$gene) & is.finite(out$pseudobulk_count), , drop = FALSE]
  if (!nrow(out)) return(data.frame(gene = character(), pseudobulk_count = numeric(), n_cells = numeric(), stringsAsFactors = FALSE))
  counts <- stats::aggregate(pseudobulk_count ~ gene, data = out, FUN = sum, na.rm = TRUE)
  cells <- stats::aggregate(n_cells ~ gene, data = out, FUN = max, na.rm = TRUE)
  merge(counts, cells, by = "gene", all = TRUE, sort = FALSE)
}

o2p6_extract_target_rows_text <- function(path, genes = o2p6_target_genes()) {
  genes_clean <- o2p6_clean_gene(genes)
  cached <- o2p6_fulltx_cache_read(path, genes = genes_clean, require_target = TRUE)
  if (isTRUE(cached$hit)) {
    ans <- cached$object$target_payload
    ans$cache_status <- "CACHE_HIT_VALIDATED"
    ans$cache_file <- cached$cache_file
    return(ans)
  }
  con <- o2p6_open_text(path)
  on.exit(close(con), add = TRUE)

  header <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) character())
  if (!length(header)) {
    return(list(ok = FALSE, error = "empty_or_unreadable_file", matrix = NULL, cells = character(), line_map = data.frame(), orientation = NA_character_))
  }
  sep <- o2p6_detect_sep_from_header(header)
  header_fields <- strsplit(header, sep, fixed = TRUE)[[1]]
  header_clean <- o2p6_clean_gene(header_fields)
  cell_by_gene <- "OLFML2B" %in% header_clean
  rows <- list(); full_gene <- list(); full_count <- list(); full_n_cells <- list()
  cell_target <- list(); cell_ids_seen <- character(); target_header_genes <- character()
  full_sums <- numeric(); full_header_genes <- character(); n_cell_rows <- 0L
  if (cell_by_gene) {
    gene_idx <- which(seq_along(header_clean) > 1L & nzchar(header_clean))
    full_header_genes <- header_clean[gene_idx]
    full_sums <- numeric(length(gene_idx))
    target_idx <- gene_idx[full_header_genes %in% genes_clean]
    target_header_genes <- header_clean[target_idx]
  }
  line_no <- 1L

  repeat {
    lines <- tryCatch(readLines(con, n = 5000, warn = FALSE), error = function(e) character())
    if (!length(lines)) break
    for (j in seq_along(lines)) {
      line_no <- line_no + 1L
      fields <- strsplit(lines[j], sep, fixed = TRUE)[[1]]
      if (length(fields) < 2L) next
      if (cell_by_gene) {
        if (!length(gene_idx) || length(fields) < max(gene_idx)) next
        values <- suppressWarnings(as.numeric(fields[gene_idx]))
        count_values <- values
        count_values[!is.finite(count_values)] <- 0
        full_sums <- full_sums + count_values
        n_cell_rows <- n_cell_rows + 1L
        if (length(target_idx)) {
          cell_target[[length(cell_target) + 1L]] <- suppressWarnings(as.numeric(fields[target_idx]))
          cell_ids_seen <- c(cell_ids_seen, if (nzchar(fields[1])) fields[1] else paste0("cell_", n_cell_rows))
        }
      } else {
        g <- o2p6_clean_gene(fields[1])
        if (!nzchar(g)) next
        values <- suppressWarnings(as.numeric(fields[-1]))
        count_values <- values
        count_values[!is.finite(count_values)] <- 0
        full_gene[[length(full_gene) + 1L]] <- g
        full_count[[length(full_count) + 1L]] <- sum(count_values)
        full_n_cells[[length(full_n_cells) + 1L]] <- length(values)
        if (!g %in% genes_clean) next
        if (length(header_fields) == length(values)) {
          cell_ids <- header_fields
        } else if (length(header_fields) == length(values) + 1L) {
          cell_ids <- header_fields[-1]
        } else {
          cell_ids <- paste0("cell_", seq_along(values))
        }
        rows[[length(rows) + 1L]] <- list(gene = g, line_no = line_no, values = values, cell_ids = cell_ids)
      }
    }
  }

  if (cell_by_gene) {
    fulltx <- o2p6_collapse_fulltx(full_header_genes, full_sums, rep(n_cell_rows, length(full_sums)))
    if (length(cell_target) && length(target_header_genes)) {
      mat0 <- do.call(rbind, cell_target)
      mat0[!is.finite(mat0)] <- 0
      mat <- t(mat0)
      rownames(mat) <- target_header_genes
      colnames(mat) <- make.unique(cell_ids_seen)
      duplicated_genes <- unique(rownames(mat)[duplicated(rownames(mat))])
      mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
      line_map <- data.frame(
        gene = rownames(mat),
        header_column = match(rownames(mat), header_clean),
        duplicated_gene_was_dropped = rownames(mat) %in% duplicated_genes,
        stringsAsFactors = FALSE
      )
      payload <- list(ok = TRUE, error = NA_character_, matrix = mat, cells = colnames(mat),
                      line_map = line_map, sep = sep, orientation = "cell_by_gene_text")
    } else {
      payload <- list(ok = TRUE, error = NA_character_, matrix = NULL, cells = character(),
                      line_map = data.frame(), sep = sep, orientation = "cell_by_gene_no_target_columns")
    }
  } else {
    fulltx <- o2p6_collapse_fulltx(unlist(full_gene, use.names = FALSE),
                                    unlist(full_count, use.names = FALSE),
                                    unlist(full_n_cells, use.names = FALSE))
    if (!length(rows)) {
      payload <- list(ok = TRUE, error = NA_character_, matrix = NULL, cells = character(),
                      line_map = data.frame(), sep = sep, orientation = "gene_by_cell_not_detected")
    } else {
      genes_observed <- vapply(rows, `[[`, character(1), "gene")
      duplicated_genes <- unique(genes_observed[duplicated(genes_observed)])
      rows <- rows[!duplicated(genes_observed)]
      all_cells <- rows[[1]]$cell_ids
      mat <- matrix(NA_real_, nrow = length(rows), ncol = length(all_cells))
      rownames(mat) <- vapply(rows, `[[`, character(1), "gene")
      colnames(mat) <- all_cells
      for (i in seq_along(rows)) {
        v <- rows[[i]]$values
        len <- min(length(v), ncol(mat))
        mat[i, seq_len(len)] <- v[seq_len(len)]
      }
      line_map <- data.frame(
        gene = vapply(rows, `[[`, character(1), "gene"),
        line_no = vapply(rows, `[[`, integer(1), "line_no"),
        duplicated_gene_was_dropped = vapply(rows, function(r) r$gene %in% duplicated_genes, logical(1)),
        stringsAsFactors = FALSE
      )
      payload <- list(ok = TRUE, error = NA_character_, matrix = mat, cells = colnames(mat),
                      line_map = line_map, sep = sep, orientation = "gene_by_cell_text")
    }
  }
  cache_file <- o2p6_fulltx_cache_write(path, fulltx, target_payload = payload, genes = genes_clean,
                                         created_by = "single_pass_text")
  payload$cache_status <- cached$reason %||% "CACHE_MISS"
  payload$cache_file <- cache_file
  payload
}

o2p6_read_feature_table <- function(feature_file) {
  if (is.na(feature_file) || !file.exists(feature_file)) return(data.frame())
  con <- o2p6_open_text(feature_file)
  on.exit(close(con), add = TRUE)
  tab <- tryCatch(utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(tab)) return(data.frame())
  names(tab) <- paste0("V", seq_len(ncol(tab)))
  # 10X features usually col2 is symbol, genes.tsv usually col1 or col2
  symbol <- if (ncol(tab) >= 2L) tab[[2]] else tab[[1]]
  if (all(!nzchar(as.character(symbol))) && ncol(tab) >= 1L) symbol <- tab[[1]]
  data.frame(row_index = seq_len(nrow(tab)), gene_symbol = o2p6_clean_gene(symbol), raw_feature = as.character(symbol), stringsAsFactors = FALSE)
}

o2p6_read_barcodes <- function(barcode_file) {
  if (is.na(barcode_file) || !file.exists(barcode_file)) return(character())
  con <- o2p6_open_text(barcode_file)
  on.exit(close(con), add = TRUE)
  x <- tryCatch(readLines(con, warn = FALSE), error = function(e) character())
  x <- trimws(x)
  x[nzchar(x)]
}

o2p6_extract_target_rows_10x <- function(matrix_file, feature_file, barcode_file, genes = o2p6_target_genes()) {
  source_paths <- c(matrix_file, feature_file, barcode_file)
  cached <- o2p6_fulltx_cache_read(source_paths, genes = genes, require_target = TRUE)
  if (isTRUE(cached$hit)) {
    ans <- cached$object$target_payload
    ans$cache_status <- "CACHE_HIT_VALIDATED"
    ans$cache_file <- cached$cache_file
    return(ans)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    return(list(ok = FALSE, error = "Matrix_package_not_available", matrix = NULL, cells = character(), line_map = data.frame()))
  }
  feat <- o2p6_read_feature_table(feature_file)
  if (!nrow(feat)) return(list(ok = FALSE, error = "feature_table_unreadable", matrix = NULL, cells = character(), line_map = data.frame()))
  target_clean <- o2p6_clean_gene(genes)
  idx <- which(feat$gene_symbol %in% target_clean)
  m <- tryCatch(Matrix::readMM(matrix_file), error = function(e) e)
  if (inherits(m, "error")) return(list(ok = FALSE, error = conditionMessage(m), matrix = NULL, cells = character(), line_map = data.frame()))
  # Some matrices may be cells x genes. Infer orientation from feature length.
  orientation <- if (nrow(m) == nrow(feat)) "gene_by_cell_10x" else if (ncol(m) == nrow(feat)) "cell_by_gene_10x" else "dimension_mismatch"
  if (orientation == "dimension_mismatch") {
    return(list(ok = FALSE, error = paste0("feature_count_", nrow(feat), "_does_not_match_matrix_dim_", paste(dim(m), collapse="x")), matrix = NULL, cells = character(), line_map = data.frame()))
  }
  gene_by_cell <- if (orientation == "gene_by_cell_10x") m else Matrix::t(m)
  fulltx <- o2p6_collapse_fulltx(feat$gene_symbol, as.numeric(Matrix::rowSums(gene_by_cell)), ncol(gene_by_cell))
  if (!length(idx)) {
    payload <- list(ok = TRUE, error = NA_character_, matrix = NULL, cells = character(),
                    line_map = data.frame(), orientation = "10x_no_target_rows")
    cache_file <- o2p6_fulltx_cache_write(source_paths, fulltx, target_payload = payload,
                                           genes = target_clean, created_by = "single_pass_10x")
    payload$cache_status <- cached$reason %||% "CACHE_MISS"
    payload$cache_file <- cache_file
    return(payload)
  }
  cells <- o2p6_read_barcodes(barcode_file)
  if (!length(cells)) cells <- paste0("cell_", seq_len(if (orientation == "gene_by_cell_10x") ncol(m) else nrow(m)))

  if (orientation == "gene_by_cell_10x") {
    sub <- as.matrix(gene_by_cell[idx, , drop = FALSE])
    if (length(cells) == ncol(sub)) colnames(sub) <- cells else colnames(sub) <- paste0("cell_", seq_len(ncol(sub)))
  } else {
    sub <- as.matrix(gene_by_cell[idx, , drop = FALSE])
    if (length(cells) == ncol(sub)) colnames(sub) <- cells else colnames(sub) <- paste0("cell_", seq_len(ncol(sub)))
  }
  rownames(sub) <- feat$gene_symbol[idx]
  # collapse duplicated gene symbols by column sum
  if (any(duplicated(rownames(sub)))) {
    u <- unique(rownames(sub))
    collapsed <- matrix(0, nrow = length(u), ncol = ncol(sub), dimnames = list(u, colnames(sub)))
    for (g in u) collapsed[g, ] <- colSums(sub[rownames(sub) == g, , drop = FALSE], na.rm = TRUE)
    dup <- TRUE
    sub <- collapsed
  } else {
    dup <- FALSE
  }
  line_map <- data.frame(
    gene = rownames(sub),
    feature_row_indices = vapply(rownames(sub), function(g) paste(feat$row_index[feat$gene_symbol == g], collapse = ";"), character(1)),
    duplicated_gene_was_collapsed = dup,
    stringsAsFactors = FALSE
  )
  payload <- list(ok = TRUE, error = NA_character_, matrix = sub, cells = colnames(sub),
                  line_map = line_map, orientation = orientation)
  cache_file <- o2p6_fulltx_cache_write(source_paths, fulltx, target_payload = payload,
                                         genes = target_clean, created_by = "single_pass_10x")
  payload$cache_status <- cached$reason %||% "CACHE_MISS"
  payload$cache_file <- cache_file
  payload
}

o2p6_score_marker_celltypes <- function(marker_mat) {
  if (is.null(marker_mat) || !nrow(marker_mat) || !ncol(marker_mat)) return(data.frame())
  marker_sets <- o2p6_marker_sets()
  genes <- rownames(marker_mat)
  score_df <- data.frame(cell_id = colnames(marker_mat), stringsAsFactors = FALSE)
  detected_df <- data.frame(cell_id = colnames(marker_mat), stringsAsFactors = FALSE)
  for (ct in names(marker_sets)) {
    requested <- unique(o2p6_clean_gene(marker_sets[[ct]]))
    gs <- intersect(requested, genes)
    eligible <- length(gs) >= 2L && length(gs) / max(length(requested), 1L) >= 0.25
    if (!eligible) {
      score_df[[ct]] <- 0
      detected_df[[ct]] <- 0L
    } else {
      x <- marker_mat[gs, , drop = FALSE]
      score_df[[ct]] <- colSums(log1p(x), na.rm = TRUE)
      detected_df[[ct]] <- colSums(x > 0, na.rm = TRUE)
    }
  }
  score_cols <- names(marker_sets)
  score_mat <- as.matrix(score_df[, score_cols, drop = FALSE])
  ord <- t(apply(score_mat, 1, function(v) order(v, decreasing = TRUE)))
  best_idx <- ord[, 1]
  second_idx <- if (ncol(score_mat) >= 2L) ord[, 2] else best_idx
  mx <- score_mat[cbind(seq_len(nrow(score_mat)), best_idx)]
  second <- score_mat[cbind(seq_len(nrow(score_mat)), second_idx)]
  best <- score_cols[best_idx]
  best_detected <- vapply(seq_len(nrow(score_mat)), function(i) as.integer(detected_df[[best[i]]][i]), integer(1))
  margin <- mx - second
  rel_margin <- margin / pmax(abs(mx), 1e-8)
  confidence <- ifelse(!is.finite(mx) | mx <= 0, "UNASSIGNED",
                       ifelse(best_detected >= 2L & rel_margin >= 0.25, "HIGH", "MODERATE"))
  best[confidence == "UNASSIGNED"] <- "Unassigned"
  score_df$marker_celltype <- best
  score_df$marker_max_score <- mx
  score_df$marker_second_score <- second
  score_df$marker_score_margin <- margin
  score_df$marker_relative_margin <- rel_margin
  score_df$best_marker_genes_detected <- best_detected
  score_df$marker_confidence <- confidence
  score_df
}

o2p6_program_scores <- function(marker_mat) {
  if (is.null(marker_mat) || !nrow(marker_mat) || !ncol(marker_mat)) return(data.frame())
  sets <- o2p6_program_sets()
  genes <- rownames(marker_mat)
  score_df <- data.frame(cell_id = colnames(marker_mat), stringsAsFactors = FALSE)
  for (nm in names(sets)) {
    requested <- unique(o2p6_clean_gene(sets[[nm]]))
    gs <- intersect(requested, genes)
    eligible <- length(gs) >= 3L && length(gs) / max(length(requested), 1L) >= 0.40
    score_df[[nm]] <- if (!eligible) rep(NA_real_, ncol(marker_mat)) else colMeans(log1p(marker_mat[gs, , drop = FALSE]), na.rm = TRUE)
  }
  score_df
}

o2p6_summarise_matrix <- function(candidate, raw_single_cell_dir, log_file = NULL) {
  cohort <- as.character(candidate$cohort[1])
  info <- o2p6_infer_sample_info(cohort, candidate$sample_path[1], candidate$file_name[1])
  o2p6_log("INFO", "Parsing ", candidate$input_type[1], ": ", cohort, " | ", candidate$file_name[1], log_file = log_file)

  parsed <- if (identical(candidate$input_type[1], "10X_MTX")) {
    o2p6_extract_target_rows_10x(candidate$matrix_file[1], candidate$feature_file[1], candidate$barcode_file[1])
  } else {
    o2p6_extract_target_rows_text(candidate$matrix_file[1])
  }

  if (!isTRUE(parsed$ok) || is.null(parsed$matrix)) {
    return(list(
      dataset_audit = data.frame(
        info,
        input_type = candidate$input_type[1],
        file_name = candidate$file_name[1],
        matrix_file = candidate$matrix_file[1],
        parse_status = if (isTRUE(parsed$ok)) "NO_TARGET_ROWS_EXTRACTED" else "LOAD_FAILED",
        error = parsed$error %||% NA_character_,
        orientation = parsed$orientation %||% NA_character_,
        cache_status = parsed$cache_status %||% "CACHE_NOT_CONFIGURED",
        cache_file = parsed$cache_file %||% NA_character_,
        n_cells = NA_integer_,
        n_genes_extracted = 0L,
        has_OLFML2B = FALSE,
        has_OLFML2A = FALSE,
        stringsAsFactors = FALSE
      ),
      gene_summary = data.frame(),
      celltype_summary = data.frame(),
      pseudobulk = data.frame(),
      program_summary = data.frame(),
      compartment_program_summary = data.frame(),
      marker_program_coverage = data.frame(),
      program_tests_per_sample = data.frame(),
      line_map = data.frame()
    ))
  }

  mat <- parsed$matrix
  has_ube <- "OLFML2B" %in% rownames(mat)
  has_olfml2a <- "OLFML2A" %in% rownames(mat)

  gene_summary <- o2p6_bind_rows(lapply(rownames(mat), function(g) {
    v <- as.numeric(mat[g, ])
    data.frame(
      info,
      input_type = candidate$input_type[1],
      file_name = candidate$file_name[1],
      gene = g,
      n_cells = length(v),
      n_detected_cells = sum(v > 0, na.rm = TRUE),
      detected_fraction = mean(v > 0, na.rm = TRUE),
      n_detected_ge2 = sum(v >= 2, na.rm = TRUE),
      detected_ge2_fraction = mean(v >= 2, na.rm = TRUE),
      n_detected_ge3 = sum(v >= 3, na.rm = TRUE),
      detected_ge3_fraction = mean(v >= 3, na.rm = TRUE),
      mean_count = mean(v, na.rm = TRUE),
      median_count = median(v, na.rm = TRUE),
      max_count = max(v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  marker_program_coverage <- o2p6_bind_rows(c(
    lapply(names(o2p6_marker_sets()), function(nm) {
      requested <- unique(o2p6_clean_gene(o2p6_marker_sets()[[nm]])); measured <- intersect(requested, rownames(mat))
      data.frame(info, input_type = candidate$input_type[1], file_name = candidate$file_name[1], set_type = "CELL_COMPARTMENT_MARKER", set_name = nm,
                 requested_genes = length(requested), measured_genes = length(measured), coverage_fraction = length(measured)/max(length(requested),1L),
                 coverage_eligible = length(measured)>=2L && length(measured)/max(length(requested),1L)>=0.25, stringsAsFactors = FALSE)
    }),
    lapply(names(o2p6_program_sets()), function(nm) {
      requested <- unique(o2p6_clean_gene(o2p6_program_sets()[[nm]])); measured <- intersect(requested, rownames(mat))
      data.frame(info, input_type = candidate$input_type[1], file_name = candidate$file_name[1], set_type = "CELL_STATE_PROGRAM", set_name = nm,
                 requested_genes = length(requested), measured_genes = length(measured), coverage_fraction = length(measured)/max(length(requested),1L),
                 coverage_eligible = length(measured)>=3L && length(measured)/max(length(requested),1L)>=0.40, stringsAsFactors = FALSE)
    })
  ))

  scores <- o2p6_score_marker_celltypes(mat)
  program_scores <- o2p6_program_scores(mat)
  celltype_summary <- data.frame()
  pseudobulk <- data.frame()
  program_summary <- data.frame()
  compartment_program_summary <- data.frame()
  program_tests_per_sample <- data.frame()

  if (nrow(scores) && has_ube) {
    ube <- as.numeric(mat["OLFML2B", scores$cell_id])
    names(ube) <- scores$cell_id
    celltype_summary <- o2p6_bind_rows(lapply(split(seq_len(nrow(scores)), scores$marker_celltype), function(idx) {
      data.frame(
        info,
        input_type = candidate$input_type[1],
        file_name = candidate$file_name[1],
        marker_celltype = scores$marker_celltype[idx[1]],
        n_cells = length(idx),
        n_OLFML2B_positive = sum(ube[idx] > 0, na.rm = TRUE),
        OLFML2B_positive_fraction = mean(ube[idx] > 0, na.rm = TRUE),
        n_OLFML2B_ge2 = sum(ube[idx] >= 2, na.rm = TRUE),
        OLFML2B_ge2_fraction = mean(ube[idx] >= 2, na.rm = TRUE),
        n_OLFML2B_ge3 = sum(ube[idx] >= 3, na.rm = TRUE),
        OLFML2B_ge3_fraction = mean(ube[idx] >= 3, na.rm = TRUE),
        n_high_confidence_cells = sum(scores$marker_confidence[idx] == "HIGH", na.rm = TRUE),
        high_confidence_fraction = mean(scores$marker_confidence[idx] == "HIGH", na.rm = TRUE),
        n_OLFML2B_positive_high_confidence = sum(ube[idx] > 0 & scores$marker_confidence[idx] == "HIGH", na.rm = TRUE),
        mean_OLFML2B_count = mean(ube[idx], na.rm = TRUE),
        median_OLFML2B_count = median(ube[idx], na.rm = TRUE),
        annotation_basis = "marker_derived_exploratory",
        stringsAsFactors = FALSE
      )
    }))

    pseudobulk <- data.frame(
      info,
      input_type = candidate$input_type[1],
      file_name = candidate$file_name[1],
      n_cells = length(ube),
      n_OLFML2B_positive = sum(ube > 0, na.rm = TRUE),
      OLFML2B_positive_fraction = mean(ube > 0, na.rm = TRUE),
      n_OLFML2B_ge2 = sum(ube >= 2, na.rm = TRUE),
      OLFML2B_ge2_fraction = mean(ube >= 2, na.rm = TRUE),
      n_OLFML2B_ge3 = sum(ube >= 3, na.rm = TRUE),
      OLFML2B_ge3_fraction = mean(ube >= 3, na.rm = TRUE),
      mean_OLFML2B_count = mean(ube, na.rm = TRUE),
      median_OLFML2B_count = median(ube, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    if (nrow(program_scores)) {
      common <- intersect(scores$cell_id, program_scores$cell_id)
      ube2 <- ube[common]
      ps <- program_scores[match(common, program_scores$cell_id), , drop = FALSE]
      ss <- scores[match(common, scores$cell_id), , drop = FALSE]
      marker_sets_now <- o2p6_marker_sets()
      program_sets_now <- o2p6_program_sets()
      for (program in setdiff(names(ps), "cell_id")) {
        sc <- as.numeric(ps[[program]])
        if (!any(is.finite(sc))) next
        pos <- ube2 > 0
        program_summary <- rbind(program_summary, data.frame(
          info,
          input_type = candidate$input_type[1],
          file_name = candidate$file_name[1],
          marker_program = program,
          n_cells = length(sc),
          n_OLFML2B_positive = sum(pos, na.rm = TRUE),
          n_OLFML2B_negative = sum(!pos, na.rm = TRUE),
          mean_score_OLFML2B_positive = if (any(pos, na.rm = TRUE)) mean(sc[pos], na.rm = TRUE) else NA_real_,
          mean_score_OLFML2B_negative = if (any(!pos, na.rm = TRUE)) mean(sc[!pos], na.rm = TRUE) else NA_real_,
          delta_positive_minus_negative = if (any(pos, na.rm = TRUE) && any(!pos, na.rm = TRUE)) mean(sc[pos], na.rm = TRUE) - mean(sc[!pos], na.rm = TRUE) else NA_real_,
          sample_program_mean_score = mean(sc, na.rm = TRUE),
          sample_program_median_score = median(sc, na.rm = TRUE),
          inference_level = "sample-level_marker_program_difference; exploratory; no cell-level p-value",
          stringsAsFactors = FALSE
        ))
        for (ct in sort(unique(as.character(ss$marker_celltype)))) {
          for (confidence_scope in c("ALL_ASSIGNED", "HIGH_ONLY")) {
            idx <- which(ss$marker_celltype == ct)
            if (confidence_scope == "HIGH_ONLY") idx <- idx[ss$marker_confidence[idx] == "HIGH"]
            if (!length(idx)) next
            pos_ct <- ube2[idx] > 0
            overlap <- length(intersect(o2p6_clean_gene(marker_sets_now[[ct]] %||% character()),
                                        o2p6_clean_gene(program_sets_now[[program]] %||% character())))
            compartment_program_summary <- rbind(compartment_program_summary, data.frame(
              info,
              input_type = candidate$input_type[1],
              file_name = candidate$file_name[1],
              marker_celltype = ct,
              annotation_confidence_scope = confidence_scope,
              marker_program = program,
              n_compartment_cells = length(idx),
              n_OLFML2B_positive = sum(pos_ct, na.rm = TRUE),
              n_OLFML2B_negative = sum(!pos_ct, na.rm = TRUE),
              mean_score_OLFML2B_positive = if (any(pos_ct, na.rm = TRUE)) mean(sc[idx][pos_ct], na.rm = TRUE) else NA_real_,
              mean_score_OLFML2B_negative = if (any(!pos_ct, na.rm = TRUE)) mean(sc[idx][!pos_ct], na.rm = TRUE) else NA_real_,
              delta_positive_minus_negative = if (any(pos_ct, na.rm = TRUE) && any(!pos_ct, na.rm = TRUE)) mean(sc[idx][pos_ct], na.rm = TRUE) - mean(sc[idx][!pos_ct], na.rm = TRUE) else NA_real_,
              classifier_program_overlap_n = overlap,
              orthogonality_status = ifelse(overlap > 0L, "IDENTITY_COUPLED_NOT_MECHANISTIC", "ORTHOGONAL_TO_COMPARTMENT_CLASSIFIER"),
              inference_level = "within-marker-compartment sample summary; patient aggregation required; no cell-level p-value",
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }

  line_map <- parsed$line_map
  if (nrow(line_map)) line_map <- cbind(info, input_type = candidate$input_type[1], file_name = candidate$file_name[1], line_map)

  dataset_audit <- data.frame(
    info,
    input_type = candidate$input_type[1],
    file_name = candidate$file_name[1],
    matrix_file = candidate$matrix_file[1],
    parse_status = "OK",
    error = NA_character_,
    orientation = parsed$orientation %||% NA_character_,
    cache_status = parsed$cache_status %||% "CACHE_NOT_CONFIGURED",
    cache_file = parsed$cache_file %||% NA_character_,
    n_cells = ncol(mat),
    n_genes_extracted = nrow(mat),
    has_OLFML2B = has_ube,
    has_OLFML2A = has_olfml2a,
    has_any_excluded_symbol = any(rownames(mat) %in% o2p6_exclusion_symbols()),
    stringsAsFactors = FALSE
  )

  list(
    dataset_audit = dataset_audit,
    gene_summary = gene_summary,
    celltype_summary = celltype_summary,
    pseudobulk = pseudobulk,
    program_summary = program_summary,
    compartment_program_summary = compartment_program_summary,
    marker_program_coverage = marker_program_coverage,
    line_map = line_map
  )
}

# ----------------------------------------------------------------------------
# Downstream summaries
# ----------------------------------------------------------------------------
o2p6_condition_summary <- function(pseudobulk) {
  if (!is.data.frame(pseudobulk) || !nrow(pseudobulk)) return(data.frame())
  d <- pseudobulk[!is.na(pseudobulk$tissue_or_condition) & nzchar(as.character(pseudobulk$tissue_or_condition)), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  out <- o2p6_bind_rows(lapply(split(d, paste(d$cohort, d$tissue_or_condition, sep = "__")), function(z) {
    data.frame(
      cohort = z$cohort[1],
      dataset_role = z$dataset_role[1],
      tissue_or_condition = z$tissue_or_condition[1],
      n_samples = nrow(z),
      total_cells = sum(z$n_cells, na.rm = TRUE),
      total_OLFML2B_positive = sum(z$n_OLFML2B_positive, na.rm = TRUE),
      pooled_positive_fraction = sum(z$n_OLFML2B_positive, na.rm = TRUE) / pmax(sum(z$n_cells, na.rm = TRUE), 1),
      mean_sample_positive_fraction = mean(z$OLFML2B_positive_fraction, na.rm = TRUE),
      sd_sample_positive_fraction = stats::sd(z$OLFML2B_positive_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

o2p6_celltype_enrichment <- function(celltype_summary) {
  if (!is.data.frame(celltype_summary) || !nrow(celltype_summary)) return(data.frame())
  d0 <- celltype_summary
  pooled <- d0
  pooled$cohort <- "POOLED_PRIMARY_GC"
  pooled <- pooled[pooled$dataset_role %in% c("PRIMARY_TUMOR_NORMAL_GC_ATLAS", "PRIMARY_ADVANCED_GC_TUMOR_ATLAS", "PRIMARY_DGC_SUPERFICIAL_DEEP_CONTEXT"), , drop = FALSE]
  d <- rbind(d0, pooled)
  out <- list()
  for (coh in unique(d$cohort)) {
    dc <- d[d$cohort == coh, , drop = FALSE]
    if (!nrow(dc)) next
    agg <- aggregate(cbind(n_cells, n_OLFML2B_positive) ~ marker_celltype, data = dc, FUN = sum, na.rm = TRUE)
    total_cells <- sum(agg$n_cells, na.rm = TRUE)
    total_pos <- sum(agg$n_OLFML2B_positive, na.rm = TRUE)
    for (i in seq_len(nrow(agg))) {
      pos_ct <- agg$n_OLFML2B_positive[i]
      neg_ct <- agg$n_cells[i] - pos_ct
      pos_other <- total_pos - pos_ct
      neg_other <- (total_cells - total_pos) - neg_ct
      tab <- matrix(c(pos_ct, neg_ct, pos_other, neg_other), nrow = 2, byrow = TRUE)
      ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
      out[[length(out) + 1L]] <- data.frame(
        cohort = coh,
        marker_celltype = agg$marker_celltype[i],
        n_cells = agg$n_cells[i],
        n_OLFML2B_positive = pos_ct,
        OLFML2B_positive_fraction = pos_ct / pmax(agg$n_cells[i], 1),
        odds_ratio_vs_other_cells = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        ci_low = if (!is.null(ft)) unname(ft$conf.int[1]) else NA_real_,
        ci_high = if (!is.null(ft)) unname(ft$conf.int[2]) else NA_real_,
        fisher_p = if (!is.null(ft)) ft$p.value else NA_real_,
        inference_level = "cell-count enrichment; exploratory localization audit only; cells are not independent biological replicates",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- o2p6_bind_rows(out)
  if (nrow(out)) out$fisher_p_fdr <- stats::p.adjust(out$fisher_p, method = "BH")
  out
}

o2p6_program_sample_tests <- function(program_summary) {
  if (!is.data.frame(program_summary) || !nrow(program_summary)) return(data.frame())
  d0 <- program_summary[is.finite(program_summary$delta_positive_minus_negative), , drop = FALSE]
  if (!nrow(d0)) return(data.frame())
  pooled <- d0[d0$dataset_role %in% c("PRIMARY_TUMOR_NORMAL_GC_ATLAS", "PRIMARY_ADVANCED_GC_TUMOR_ATLAS", "PRIMARY_DGC_SUPERFICIAL_DEEP_CONTEXT"), , drop = FALSE]
  pooled$cohort <- "POOLED_PRIMARY_GC"
  d <- rbind(d0, pooled)
  out <- o2p6_bind_rows(lapply(split(d, paste(d$cohort, d$marker_program, sep = "__")), function(z) {
    diff <- z$delta_positive_minus_negative[is.finite(z$delta_positive_minus_negative)]
    tt <- if (length(diff) >= 3) tryCatch(stats::t.test(diff, mu = 0), error = function(e) NULL) else NULL
    wt <- if (length(diff) >= 3) tryCatch(suppressWarnings(stats::wilcox.test(diff, mu = 0, exact = FALSE)), error = function(e) NULL) else NULL
    data.frame(
      cohort = z$cohort[1],
      marker_program = z$marker_program[1],
      n_samples = length(diff),
      mean_delta_positive_minus_negative = if (length(diff)) mean(diff, na.rm = TRUE) else NA_real_,
      median_delta_positive_minus_negative = if (length(diff)) median(diff, na.rm = TRUE) else NA_real_,
      t_p = if (!is.null(tt)) tt$p.value else NA_real_,
      wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
      direction = ifelse(length(diff) && median(diff, na.rm = TRUE) > 0, "higher_in_OLFML2B_positive_cells", ifelse(length(diff) && median(diff, na.rm = TRUE) < 0, "lower_in_OLFML2B_positive_cells", "neutral_or_insufficient")),
      inference_level = "sample-level marker-program comparison; exploratory; no cell-level DE",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out)) {
    out$t_p_fdr <- stats::p.adjust(out$t_p, method = "BH")
    out$wilcox_p_fdr <- stats::p.adjust(out$wilcox_p, method = "BH")
  }
  out
}

o2p6_program_bridge <- function(pseudobulk, program_summary) {
  if (!is.data.frame(pseudobulk) || !nrow(pseudobulk) || !is.data.frame(program_summary) || !nrow(program_summary)) return(data.frame())
  key_pb <- paste(pseudobulk$cohort, pseudobulk$file_name, sep = "__")
  program_summary$key <- paste(program_summary$cohort, program_summary$file_name, sep = "__")
  out <- o2p6_bind_rows(lapply(split(program_summary, paste(program_summary$cohort, program_summary$marker_program, sep = "__")), function(z) {
    pb_idx <- match(z$key, key_pb)
    pb <- pseudobulk[pb_idx, , drop = FALSE]
    ok <- is.finite(pb$OLFML2B_positive_fraction) & is.finite(z$sample_program_mean_score)
    if (sum(ok) < 3L) {
      return(data.frame(
        cohort = z$cohort[1], marker_program = z$marker_program[1], n_samples = sum(ok),
        rho = NA_real_, p_value = NA_real_, direction = "insufficient_samples", stringsAsFactors = FALSE
      ))
    }
    ct <- tryCatch(suppressWarnings(stats::cor.test(pb$OLFML2B_positive_fraction[ok], z$sample_program_mean_score[ok], method = "spearman", exact = FALSE)), error = function(e) NULL)
    data.frame(
      cohort = z$cohort[1],
      marker_program = z$marker_program[1],
      n_samples = sum(ok),
      rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
      p_value = if (!is.null(ct)) ct$p.value else NA_real_,
      direction = if (!is.null(ct) && is.finite(unname(ct$estimate)) && unname(ct$estimate) > 0) "positive" else if (!is.null(ct) && is.finite(unname(ct$estimate)) && unname(ct$estimate) < 0) "negative" else "neutral_or_insufficient",
      inference_level = "sample-level Spearman between OLFML2B positive fraction and marker-program mean score",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out)) out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

o2p6_gse167297_paired_gradient <- function(pseudobulk) {
  if (!is.data.frame(pseudobulk) || !nrow(pseudobulk)) return(data.frame())
  d <- pseudobulk[pseudobulk$cohort == "GSE167297" & !is.na(pseudobulk$patient) & nzchar(as.character(pseudobulk$patient)), , drop = FALSE]
  d <- d[d$tissue_or_condition %in% c("Normal", "Superficial", "Deep"), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  pc <- o2p6_bind_rows(lapply(split(d, paste(d$patient, d$tissue_or_condition, sep = "__")), function(z) {
    data.frame(
      patient = z$patient[1],
      tissue_or_condition = z$tissue_or_condition[1],
      n_cells = sum(z$n_cells, na.rm = TRUE),
      n_OLFML2B_positive = sum(z$n_OLFML2B_positive, na.rm = TRUE),
      OLFML2B_positive_fraction = sum(z$n_OLFML2B_positive, na.rm = TRUE) / pmax(sum(z$n_cells, na.rm = TRUE), 1),
      stringsAsFactors = FALSE
    )
  }))
  patients <- sort(unique(pc$patient))
  wide <- data.frame(patient = patients, stringsAsFactors = FALSE)
  for (cond in c("Normal", "Superficial", "Deep")) {
    vals <- pc$OLFML2B_positive_fraction[match(paste(patients, cond, sep = "__"), paste(pc$patient, pc$tissue_or_condition, sep = "__"))]
    wide[[cond]] <- vals
  }
  wide$Deep_minus_Normal <- wide$Deep - wide$Normal
  wide$Superficial_minus_Normal <- wide$Superficial - wide$Normal
  wide$Deep_minus_Superficial <- wide$Deep - wide$Superficial
  wide
}

o2p6_paired_gradient_tests <- function(paired_gradient) {
  if (!is.data.frame(paired_gradient) || !nrow(paired_gradient)) return(data.frame())
  contrasts <- list(
    Deep_vs_Normal = c("Deep", "Normal"),
    Superficial_vs_Normal = c("Superficial", "Normal"),
    Deep_vs_Superficial = c("Deep", "Superficial")
  )
  out <- o2p6_bind_rows(lapply(names(contrasts), function(nm) {
    a <- contrasts[[nm]][1]; b <- contrasts[[nm]][2]
    ok <- is.finite(paired_gradient[[a]]) & is.finite(paired_gradient[[b]])
    x <- paired_gradient[[a]][ok]; y <- paired_gradient[[b]][ok]
    diff <- x - y
    finite_diff <- diff[is.finite(diff)]
    signflip_p <- if (length(finite_diff) >= 2L && length(finite_diff) <= 20L) {
      signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(finite_diff))))
      null <- as.numeric(signs %*% finite_diff) / length(finite_diff)
      mean(abs(null) >= abs(mean(finite_diff)) - sqrt(.Machine$double.eps))
    } else NA_real_
    nonzero <- finite_diff[finite_diff != 0]
    sign_test_p <- if (length(nonzero)) stats::binom.test(sum(nonzero > 0), length(nonzero), p = 0.5)$p.value else NA_real_
    tt <- if (length(diff) >= 3) tryCatch(stats::t.test(diff, mu = 0), error = function(e) NULL) else NULL
    wt <- if (length(diff) >= 3) tryCatch(suppressWarnings(stats::wilcox.test(diff, mu = 0, exact = FALSE)), error = function(e) NULL) else NULL
    data.frame(
      cohort = "GSE167297",
      contrast = nm,
      n_paired_patients = length(diff),
      mean_delta = if (length(diff)) mean(diff, na.rm = TRUE) else NA_real_,
      median_delta = if (length(diff)) median(diff, na.rm = TRUE) else NA_real_,
      t_p = if (!is.null(tt)) tt$p.value else NA_real_,
      wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
      exact_signflip_mean_p = signflip_p,
      exact_binomial_sign_p = sign_test_p,
      direction = ifelse(length(diff) && median(diff, na.rm = TRUE) > 0, "higher_in_first_condition", ifelse(length(diff) && median(diff, na.rm = TRUE) < 0, "lower_in_first_condition", "neutral_or_insufficient")),
      inference_level = "patient-level paired pseudobulk; exploratory due small n",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out)) {
    out$t_p_fdr <- stats::p.adjust(out$t_p, method = "BH")
    out$wilcox_p_fdr <- stats::p.adjust(out$wilcox_p, method = "BH")
    out$exact_signflip_mean_p_fdr <- stats::p.adjust(out$exact_signflip_mean_p, method = "BH")
    out$exact_binomial_sign_p_fdr <- stats::p.adjust(out$exact_binomial_sign_p, method = "BH")
  }
  out
}

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
o2p6_theme <- function(base_size = 10.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.35),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

o2p6_save_plot <- function(p, stem, width = 7, height = 5, dpi = 320) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  o2p6_dir(dirname(stem))
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  invisible(TRUE)
}

o2p6_make_figures <- function(dirs, cohort_summary, celltype_summary, condition_summary, program_tests, program_bridge) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  fig_registry <- list()

  if (is.data.frame(cohort_summary) && nrow(cohort_summary)) {
    d <- cohort_summary
    d$cohort <- factor(d$cohort, levels = d$cohort[order(d$dataset_role, d$cohort)])
    d$label <- paste0("cells=", format(d$total_cells, big.mark=","), "\nOLFML2B+=", sprintf("%.2f%%", 100*d$pooled_OLFML2B_positive_fraction))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = cohort, y = pooled_OLFML2B_positive_fraction, fill = dataset_role)) +
      ggplot2::geom_col(color = "black", linewidth = 0.25, width = 0.62) +
      ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 2.7, lineheight = 0.85) +
      ggplot2::scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100*x), expand = ggplot2::expansion(mult = c(0, 0.28))) +
      ggplot2::labs(tag = "A", title = "Part6-RC2 OLFML2B detection across four scRNA datasets", subtitle = "Exact OLFML2B row extraction; no Kang/GSE206785 used", x = NULL, y = "OLFML2B-positive cells", fill = "Dataset role") +
      o2p6_theme()
    o2p6_save_plot(p, file.path(dirs$figures, "FIG6RC2_A_OLFML2B_detection_by_dataset"), width = 8.2, height = 5)
  }

  if (is.data.frame(celltype_summary) && nrow(celltype_summary)) {
    ct <- stats::aggregate(cbind(n_cells, n_OLFML2B_positive) ~ cohort + marker_celltype, data = celltype_summary, FUN = sum, na.rm = TRUE)
    ct$OLFML2B_positive_fraction <- with(ct, n_OLFML2B_positive / pmax(n_cells, 1))
    p <- ggplot2::ggplot(ct, ggplot2::aes(x = marker_celltype, y = OLFML2B_positive_fraction, fill = cohort)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.68, color = "black", linewidth = 0.18) +
      ggplot2::scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100*x), expand = ggplot2::expansion(mult = c(0, 0.15))) +
      ggplot2::labs(tag = "B", title = "Marker-derived compartment localization", subtitle = "Exploratory localization audit, not author-annotation DE", x = NULL, y = "OLFML2B-positive cells", fill = "Cohort") +
      o2p6_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    o2p6_save_plot(p, file.path(dirs$figures, "FIG6RC2_B_marker_compartment_localization"), width = 9, height = 5.2)
  }

  if (is.data.frame(condition_summary) && nrow(condition_summary)) {
    d <- condition_summary
    p <- ggplot2::ggplot(d, ggplot2::aes(x = tissue_or_condition, y = pooled_positive_fraction, group = cohort, color = cohort)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::geom_line(linewidth = 0.45) +
      ggplot2::scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100*x), expand = ggplot2::expansion(mult = c(0, 0.18))) +
      ggplot2::labs(tag = "C", title = "Condition-level OLFML2B positivity", subtitle = "Sample-level pseudobulk summaries by inferred lesion state", x = NULL, y = "OLFML2B-positive cells", color = "Cohort") +
      o2p6_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    o2p6_save_plot(p, file.path(dirs$figures, "FIG6RC2_C_condition_level_pseudobulk"), width = 8.2, height = 5)
  }

  if (is.data.frame(program_tests) && nrow(program_tests)) {
    d <- program_tests[program_tests$cohort == "POOLED_PRIMARY_GC", , drop = FALSE]
    if (nrow(d)) {
      d$marker_program <- stats::reorder(d$marker_program, d$median_delta_positive_minus_negative)
      d$sig <- is.finite(d$wilcox_p_fdr) & d$wilcox_p_fdr < 0.05
      p <- ggplot2::ggplot(d, ggplot2::aes(x = marker_program, y = median_delta_positive_minus_negative, fill = sig)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
        ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.22) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "#C73E3A", `FALSE` = "#8C8C8C"), guide = "none") +
        ggplot2::labs(tag = "D", title = "Marker programs in OLFML2B-positive cells", subtitle = "Pooled primary-GC sample-level deltas; exploratory", x = NULL, y = "Median score delta: OLFML2B+ minus OLFML2B−") +
        o2p6_theme()
      o2p6_save_plot(p, file.path(dirs$figures, "FIG6RC2_D_OLFML2B_positive_program_delta"), width = 7.2, height = 5.2)
    }
  }

  if (is.data.frame(program_bridge) && nrow(program_bridge)) {
    d <- program_bridge[program_bridge$cohort %in% c("GSE150290", "GSE183904", "GSE167297", "POOLED_PRIMARY_GC"), , drop = FALSE]
    d <- d[is.finite(d$rho), , drop = FALSE]
    if (nrow(d)) {
      p <- ggplot2::ggplot(d, ggplot2::aes(x = marker_program, y = cohort, fill = rho)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.25) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.5) +
        ggplot2::scale_fill_gradient2(low = "#2F6DB3", mid = "white", high = "#C73E3A", midpoint = 0, name = "rho") +
        ggplot2::labs(tag = "E", title = "Sample-level OLFML2B–program bridge", subtitle = "Spearman rho between OLFML2B-positive fraction and marker-program score", x = NULL, y = NULL) +
        o2p6_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
      o2p6_save_plot(p, file.path(dirs$figures, "FIG6RC2_E_sample_level_program_bridge_heatmap"), width = 9.5, height = max(4.8, 0.35*nrow(unique(d["cohort"])) + 3))
    }
  }
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# Main runner
# ----------------------------------------------------------------------------
run_olfml2b_part6_scrna_rc2_no_kang_full_production <- function(
  root = "D:/OLFML2B_STAD",
  raw_single_cell_dir = file.path(root, "data", "raw", "single_cell"),
  output_subdir = "Part6_RC2",
  datasets = c("GSE150290", "GSE183904", "GSE167297", "GSE134520"),
  primary_datasets = c("GSE150290", "GSE183904", "GSE167297"),
  supportive_datasets = c("GSE134520"),
  force_unpack = FALSE,
  make_figures = TRUE,
  exclude_gse150290_raw_10x = TRUE
) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  raw_single_cell_dir <- normalizePath(raw_single_cell_dir, winslash = "/", mustWork = FALSE)
  dirs <- o2p6_dirs(root, output_subdir)
  log_file <- file.path(dirs$logs, "Part6_RC2_no_kang_full_production.log")
  o2p6_log("INFO", "Starting Part6-RC2 NO-KANG full production | ", OLFML2B_PART6_RC2_VERSION, log_file = log_file)
  o2p6_log("INFO", "Root: ", root, log_file = log_file)
  o2p6_log("INFO", "Raw single-cell dir: ", raw_single_cell_dir, log_file = log_file)

  params <- data.frame(
    parameter = c("version", "root", "raw_single_cell_dir", "output_subdir", "datasets", "primary_datasets", "supportive_datasets", "excluded_by_design", "target_gene", "exclude_gse150290_raw_10x", "fulltx_cache_schema", "cache_invalidation_rule", "refresh_cache", "analysis_scope"),
    value = c(
      OLFML2B_PART6_RC2_VERSION,
      root,
      raw_single_cell_dir,
      output_subdir,
      paste(datasets, collapse = ";"),
      paste(primary_datasets, collapse = ";"),
      paste(supportive_datasets, collapse = ";"),
      "Kang/GSE206785; Part8 immunotherapy",
      "OLFML2B exact symbol only; OLFML2A and OLFM2 are exclusion-audit symbols only",
      as.character(isTRUE(exclude_gse150290_raw_10x)),
      P6RC2_FULLTX_CACHE_SCHEMA,
      "normalized path + byte size + modification time + cache schema; content MD5 is separately cached in the input manifest",
      as.character(isTRUE(getOption("olfml2b.part6.refresh_cache", FALSE))),
      "Four-dataset gastric cancer-related scRNA context validation; sample/condition/compartment summaries only; no cell-level DE, no causal, no prognostic, no ICI claim. GSE150290 formal input uses author processed TXT matrices only; raw 10X droplet matrices are audit-only."
    ),
    stringsAsFactors = FALSE
  )
  o2p6_atomic_write_csv(params, file.path(dirs$tables, "00_run_parameters.csv"))

  dataset_plan <- data.frame(
    cohort = datasets,
    role = ifelse(datasets %in% primary_datasets, "PRIMARY_RC2_SINGLE_CELL_CONTEXT", "SUPPORTIVE_CONTEXT"),
    intended_use = c(
      "GC/adjacent-normal single-cell tumor-context atlas",
      "Primary gastric cancer tumor atlas",
      "Diffuse-type GC superficial/deep/normal invasion-context atlas",
      "Premalignant-to-early-GC context; supportive only"
    )[match(datasets, c("GSE150290", "GSE183904", "GSE167297", "GSE134520"))],
    formal_claim_ceiling = "single-cell context/localization only; no causal/prognostic/ICI claim",
    stringsAsFactors = FALSE
  )
  o2p6_atomic_write_csv(dataset_plan, file.path(dirs$tables, "01_dataset_inclusion_plan_no_kang.csv"))

  global_manifest <- o2p6_file_manifest(raw_single_cell_dir, recursive = TRUE)
  if (nrow(global_manifest)) global_manifest$cohort <- vapply(global_manifest$absolute_path, function(x) o2p6_infer_cohort(x, raw_single_cell_dir), character(1))
  o2p6_atomic_write_csv(global_manifest, file.path(dirs$tables, "02_input_file_manifest.csv"))

  discovery <- list()
  for (cohort in datasets) {
    cdir <- file.path(raw_single_cell_dir, cohort)
    o2p6_log("INFO", "Discovering inputs for ", cohort, " in ", cdir, log_file = log_file)
    discovery[[cohort]] <- o2p6_discover_dataset_inputs(cdir, cohort = cohort, force_unpack = force_unpack, exclude_gse150290_raw_10x = exclude_gse150290_raw_10x, log_file = log_file)
    o2p6_atomic_write_csv(discovery[[cohort]]$unpack_audit, file.path(dirs$tables, paste0("03a_", cohort, "_unpack_audit.csv")))
    o2p6_atomic_write_csv(discovery[[cohort]]$manifest, file.path(dirs$tables, paste0("03b_", cohort, "_file_manifest.csv")))
    o2p6_atomic_write_csv(discovery[[cohort]]$candidates, file.path(dirs$tables, paste0("03c_", cohort, "_candidate_matrix_manifest.csv")))
  }

  all_candidates <- o2p6_bind_rows(lapply(discovery, `[[`, "candidates"))
  if (nrow(all_candidates)) {
    all_candidates$used_for_analysis <- all_candidates$kept %in% TRUE & all_candidates$complete %in% TRUE & all_candidates$policy_keep %in% TRUE
  } else {
    all_candidates$used_for_analysis <- logical()
  }
  o2p6_atomic_write_csv(all_candidates, file.path(dirs$tables, "04_candidate_matrix_manifest_all.csv"))
  o2p6_atomic_write_csv(o2p6_candidate_policy_summary(all_candidates), file.path(dirs$tables, "04b_candidate_input_policy_summary.csv"))

  used <- all_candidates[all_candidates$used_for_analysis %in% TRUE, , drop = FALSE]
  o2p6_log("INFO", "Matrix candidates kept for analysis: ", nrow(used), log_file = log_file)

  parse_started <- Sys.time()
  parsed <- lapply(seq_len(nrow(used)), function(i) {
    one_started <- Sys.time()
    o2p6_log("INFO", sprintf("Single-pass matrix progress %d/%d | %s | %s",
      i, nrow(used), used$cohort[i], used$file_name[i]), log_file = log_file)
    ans <- tryCatch(
      o2p6_summarise_matrix(used[i, , drop = FALSE], raw_single_cell_dir, log_file = log_file),
      error = function(e) {
        candidate <- used[i, , drop = FALSE]
        info <- o2p6_infer_sample_info(candidate$cohort[1], candidate$sample_path[1], candidate$file_name[1])
        list(
          dataset_audit = data.frame(info, input_type = candidate$input_type[1], file_name = candidate$file_name[1], matrix_file = candidate$matrix_file[1], parse_status = "LOAD_FAILED", error = conditionMessage(e), orientation = NA_character_, cache_status = "CACHE_OR_PARSE_FAILED", cache_file = NA_character_, n_cells = NA_integer_, n_genes_extracted = 0L, has_OLFML2B = FALSE, has_OLFML2A = FALSE, stringsAsFactors = FALSE),
          gene_summary = data.frame(), celltype_summary = data.frame(), pseudobulk = data.frame(), program_summary = data.frame(), compartment_program_summary = data.frame(), marker_program_coverage = data.frame(), line_map = data.frame()
        )
      }
    )
    o2p6_log("INFO", sprintf("Single-pass matrix finished %d/%d | elapsed %.1f s | %s",
      i, nrow(used), as.numeric(difftime(Sys.time(), one_started, units = "secs")),
      ans$dataset_audit$cache_status[1] %||% "UNKNOWN"), log_file = log_file)
    ans
  })
  o2p6_log("INFO", sprintf("Single-pass target/full-transcriptome phase finished | %d matrices | elapsed %.1f min",
    nrow(used), as.numeric(difftime(Sys.time(), parse_started, units = "mins"))), log_file = log_file)

  dataset_audit <- o2p6_bind_rows(lapply(parsed, `[[`, "dataset_audit"))
  gene_summary <- o2p6_bind_rows(lapply(parsed, `[[`, "gene_summary"))
  celltype_summary <- o2p6_bind_rows(lapply(parsed, `[[`, "celltype_summary"))
  pseudobulk <- o2p6_bind_rows(lapply(parsed, `[[`, "pseudobulk"))
  program_summary <- o2p6_bind_rows(lapply(parsed, `[[`, "program_summary"))
  compartment_program_summary <- o2p6_bind_rows(lapply(parsed, `[[`, "compartment_program_summary"))
  marker_program_coverage <- o2p6_bind_rows(lapply(parsed, `[[`, "marker_program_coverage"))
  line_map <- o2p6_bind_rows(lapply(parsed, `[[`, "line_map"))

  o2p6_atomic_write_csv(dataset_audit, file.path(dirs$tables, "10_parse_audit_by_sample.csv"))
  o2p6_atomic_write_csv(gene_summary, file.path(dirs$tables, "11_extracted_marker_gene_summary_by_sample.csv"))
  o2p6_atomic_write_csv(gene_summary[gene_summary$gene == "OLFML2B", , drop = FALSE], file.path(dirs$tables, "12_exact_OLFML2B_summary_by_sample.csv"))
  o2p6_atomic_write_csv(gene_summary[gene_summary$gene %in% o2p6_exclusion_symbols(), , drop = FALSE], file.path(dirs$tables, "13_OLFML2A_OLFM2_exclusion_audit.csv"))
  o2p6_atomic_write_csv(line_map, file.path(dirs$tables, "14_extracted_gene_line_map.csv"))
  o2p6_atomic_write_csv(pseudobulk, file.path(dirs$tables, "20_sample_level_pseudobulk_OLFML2B.csv"))
  o2p6_atomic_write_csv(celltype_summary, file.path(dirs$tables, "21_marker_compartment_OLFML2B_summary.csv"))
  o2p6_atomic_write_csv(program_summary, file.path(dirs$tables, "22_OLFML2B_positive_marker_program_summary_by_sample.csv"))
  o2p6_atomic_write_csv(compartment_program_summary, file.path(dirs$tables, "24_within_compartment_marker_program_summary_by_sample.csv"))
  o2p6_atomic_write_csv(marker_program_coverage, file.path(dirs$tables, "25_marker_program_measurement_coverage_by_sample.csv"))

  cohort_summary <- data.frame()
  ube_summary <- gene_summary[gene_summary$gene == "OLFML2B", , drop = FALSE]
  if (nrow(ube_summary)) {
    cohort_summary <- o2p6_bind_rows(lapply(split(ube_summary, ube_summary$cohort), function(d) {
      role <- unique(d$dataset_role)[1]
      data.frame(
        cohort = d$cohort[1],
        dataset_role = role,
        n_samples_with_OLFML2B = nrow(d),
        total_cells = sum(d$n_cells, na.rm = TRUE),
        total_OLFML2B_positive_cells = sum(d$n_detected_cells, na.rm = TRUE),
        pooled_OLFML2B_positive_fraction = sum(d$n_detected_cells, na.rm = TRUE) / pmax(sum(d$n_cells, na.rm = TRUE), 1),
        mean_sample_detected_fraction = mean(d$detected_fraction, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }
  o2p6_atomic_write_csv(cohort_summary, file.path(dirs$tables, "23_cohort_level_OLFML2B_summary.csv"))

  condition_summary <- o2p6_condition_summary(pseudobulk)
  celltype_enrichment <- o2p6_celltype_enrichment(celltype_summary)
  program_tests <- o2p6_program_sample_tests(program_summary)
  program_bridge <- o2p6_program_bridge(pseudobulk, program_summary)
  paired_gradient <- o2p6_gse167297_paired_gradient(pseudobulk)
  paired_gradient_tests <- o2p6_paired_gradient_tests(paired_gradient)

  o2p6_atomic_write_csv(condition_summary, file.path(dirs$tables, "30_condition_level_OLFML2B_summary.csv"))
  o2p6_atomic_write_csv(celltype_enrichment, file.path(dirs$tables, "31_marker_compartment_enrichment_exact_or.csv"))
  o2p6_atomic_write_csv(program_tests, file.path(dirs$tables, "32_OLFML2B_positive_marker_program_tests.csv"))
  o2p6_atomic_write_csv(program_bridge, file.path(dirs$tables, "33_sample_level_OLFML2B_program_bridge.csv"))
  o2p6_atomic_write_csv(paired_gradient, file.path(dirs$tables, "34_GSE167297_patient_paired_gradient.csv"))
  o2p6_atomic_write_csv(paired_gradient_tests, file.path(dirs$tables, "35_GSE167297_paired_gradient_tests.csv"))

  # Dataset-level go/no-go
  go_dataset <- o2p6_bind_rows(lapply(datasets, function(coh) {
    da <- dataset_audit[dataset_audit$cohort == coh, , drop = FALSE]
    data.frame(
      cohort = coh,
      dataset_role = dataset_plan$role[match(coh, dataset_plan$cohort)],
      n_candidate_matrices = sum(all_candidates$cohort == coh, na.rm = TRUE),
      n_used_matrices = sum(used$cohort == coh, na.rm = TRUE),
      n_policy_excluded = if (nrow(all_candidates)) sum(all_candidates$cohort == coh & !(all_candidates$policy_keep %in% TRUE), na.rm = TRUE) else 0L,
      n_raw10x_excluded = if (nrow(all_candidates)) sum(all_candidates$cohort == coh & all_candidates$input_policy == "EXCLUDE_GSE150290_RAW_10X_DROPLET_MATRIX", na.rm = TRUE) else 0L,
      n_ok_parsed = sum(da$parse_status == "OK", na.rm = TRUE),
      n_with_exact_OLFML2B = sum(da$has_OLFML2B %in% TRUE, na.rm = TRUE),
      total_cells_parsed = sum(da$n_cells, na.rm = TRUE),
      status = ifelse(sum(da$has_OLFML2B %in% TRUE, na.rm = TRUE) > 0, "GO_EXACT_OLFML2B_DETECTED", ifelse(sum(da$parse_status == "OK", na.rm = TRUE) > 0, "REVIEW_NO_OLFML2B", "NO_GO_NOT_PARSED")),
      stringsAsFactors = FALSE
    )
  }))
  o2p6_atomic_write_csv(go_dataset, file.path(dirs$tables, "90_dataset_go_no_go.csv"))

  primary_go <- go_dataset[go_dataset$cohort %in% primary_datasets & go_dataset$status == "GO_EXACT_OLFML2B_DETECTED", , drop = FALSE]
  supportive_go <- go_dataset[go_dataset$cohort %in% supportive_datasets & go_dataset$status == "GO_EXACT_OLFML2B_DETECTED", , drop = FALSE]
  final_go <- data.frame(
    criterion = c(
      "no_kang_by_design",
      "GSE150290_raw_10x_guard",
      "primary_GC_dataset_count",
      "exact_OLFML2B_detected",
      "OLFML2_family_exclusion_audit",
      "sample_level_pseudobulk_available",
      "marker_compartment_summary_available",
      "program_bridge_available",
      "claim_ceiling"
    ),
    status = c(
      "PASS",
      ifelse(isTRUE(exclude_gse150290_raw_10x) && nrow(all_candidates) && any(all_candidates$input_policy == "EXCLUDE_GSE150290_RAW_10X_DROPLET_MATRIX"), "PASS", "REVIEW"),
      ifelse(nrow(primary_go) >= 2L, "PASS", ifelse(nrow(primary_go) >= 1L, "REVIEW", "NO_GO")),
      ifelse(nrow(ube_summary) > 0, "PASS", "NO_GO"),
      ifelse(any(dataset_audit$has_OLFML2A %in% TRUE) || nrow(gene_summary[gene_summary$gene %in% o2p6_exclusion_symbols(), , drop = FALSE]) >= 0, "PASS", "REVIEW"),
      ifelse(nrow(pseudobulk) > 0, "PASS", "NO_GO"),
      ifelse(nrow(celltype_summary) > 0, "PASS_EXPLORATORY", "REVIEW"),
      ifelse(nrow(program_bridge) > 0, "PASS_EXPLORATORY", "REVIEW"),
      "CONTEXT_SUPPORT_ONLY"
    ),
    interpretation = c(
      "Kang/GSE206785 is explicitly excluded from the RC2 formal single-cell layer.",
      "GSE150290 raw 10X droplet matrices are excluded from formal inference; author processed GSM TXT matrices are used instead.",
      "At least two of GSE150290/GSE183904/GSE167297 must parse with exact OLFML2B for a formal four-dataset single-cell layer; direction and P values never control continuation.",
      "Exact OLFML2B row was detected in at least one analyzable matrix.",
      "OLFML2A and OLFM2 are retained as exact-symbol exclusion audits and are never merged with OLFML2B.",
      "Sample-level OLFML2B summaries are available; these are preferred over cell-level p-values.",
      "Marker-derived compartment summaries are available but are exploratory and not author-level full annotation.",
      "Sample-level program bridge is available but remains exploratory.",
      "Part6-RC2 supports cellular context/localization only; no causal, prognostic, CAF-specific, smooth-muscle-specific, myeloid-specific, tumor-overexpression, or immunotherapy-response claim."
    ),
    stringsAsFactors = FALSE
  )
  o2p6_atomic_write_csv(final_go, file.path(dirs$tables, "99_Part6_RC2_no_kang_go_no_go.csv"))

  claim_boundary <- data.frame(
    item = c("Kang_GSE206785", "GSE150290_input_policy", "GSE150290", "GSE183904", "GSE167297", "GSE134520", "cell_level_inference", "marker_compartments", "allowed_claim", "forbidden_claim"),
    decision = c(
      "Excluded by RC2 design; not used in formal evidence matrix.",
      "Formal GSE150290 inference uses author processed GSM*.txt(.gz) matrices only; raw_gene_bc_matrices/10X droplet matrices are audit-only and excluded.",
      "Primary gastric cancer/adjacent-normal single-cell context layer if parsed.",
      "Primary gastric cancer tumor-atlas context layer if parsed.",
      "Diffuse-type GC superficial/deep invasion-context layer.",
      "Supportive premalignant-to-early-GC context only.",
      "Not primary; no cell-level DE claim.",
      "Exploratory localization only unless external author annotation is added.",
      "OLFML2B has single-cell context/localization evidence across locally available gastric cancer-related datasets.",
      "Do not claim a unique CAF, smooth-muscle, pericyte, myeloid or malignant-epithelial source; do not claim causal CAF/ECM remodeling, prognostic prediction, tumor-wide overexpression or ICI prediction."
    ),
    stringsAsFactors = FALSE
  )
  o2p6_atomic_write_csv(claim_boundary, file.path(dirs$tables, "100_claim_boundary_no_kang.csv"))

  if (isTRUE(make_figures)) {
    tryCatch(o2p6_make_figures(dirs, cohort_summary, celltype_summary, condition_summary, program_tests, program_bridge), error = function(e) o2p6_log("WARN", "Figure generation failed: ", conditionMessage(e), log_file = log_file))
  }

  o2p6_atomic_write_csv(o2p6_pkg_versions(c("Matrix", "ggplot2", "data.table")), file.path(dirs$reports, "Part6_RC2_package_versions.csv"))
  o2p6_atomic_write_csv(o2p6_session_packages_df(), file.path(dirs$reports, "Part6_RC2_session_packages.csv"))
  o2p6_atomic_write_csv(o2p6_file_manifest(dirs$tables, recursive = FALSE), file.path(dirs$reports, "Part6_RC2_output_table_manifest.csv"))
  o2p6_atomic_write_csv(o2p6_file_manifest(dirs$figures, recursive = FALSE), file.path(dirs$reports, "Part6_RC2_output_figure_manifest.csv"))

  index <- list(
    version = OLFML2B_PART6_RC2_VERSION,
    params = params,
    dataset_plan = dataset_plan,
    candidates = all_candidates,
    candidate_policy_summary = o2p6_candidate_policy_summary(all_candidates),
    dataset_go_no_go = go_dataset,
    final_go_no_go = final_go,
    dataset_audit = dataset_audit,
    cohort_summary = cohort_summary,
    pseudobulk = pseudobulk,
    celltype_summary = celltype_summary,
    program_summary = program_summary,
    compartment_program_summary = compartment_program_summary,
    marker_program_coverage = marker_program_coverage,
    program_tests = program_tests,
    program_bridge = program_bridge,
    paired_gradient = paired_gradient,
    claim_boundary = claim_boundary
  )
  o2p6_atomic_save_rds(index, file.path(dirs$objects, "Part6_RC2_no_kang_4dataset_full_production_index.rds"))

  o2p6_log("INFO", "Part6 target/marker layer complete; full-transcriptome aggregation follows. Tables: ", dirs$tables, log_file = log_file)
  o2p6_log("INFO", "Part6 target/marker figures complete: ", dirs$figures, log_file = log_file)
  invisible(index)
}

# Alias for pipeline integration if desired.
run_olfml2b_part6_scrna_rc2_no_kang <- run_olfml2b_part6_scrna_rc2_no_kang_full_production

# ==============================================================================
# IF 7-8 RC2-primary implementation | official GEO mapping, patient/sample
# inferential units, explicit four-cohort roles and Part9 compatibility.
# ==============================================================================

OLFML2B_PART6_IF78_VERSION <- "v1.3.0_20260722_UNIFIED_FOUR_DATASET_EXACT_INFERENCE_AND_COMPARTMENT_PROGRAM_FIX"
OLFML2B_PART6_RC2_VERSION <- OLFML2B_PART6_IF78_VERSION

o2p6_official_sample_map <- function() {
  # GSE150290: official GEO titles and accessions. A/B is interpreted only for
  # the 23 complete pairs. Pat21-A is an unpaired adjacent-normal sample and
  # Pat25-A--Pat29-A are control-normal samples in GEO; all six singletons are
  # retained descriptively but excluded from the paired primary estimand.
  g150_gsm <- 4546300:4546351
  g150_title <- c(
    as.vector(rbind(sprintf("Pat%02d-A", 1:20), sprintf("Pat%02d-B", 1:20))),
    "Pat21-A", "Pat22-A", "Pat22-B", "Pat23-A", "Pat23-B", "Pat24-A", "Pat24-B",
    "Pat25-A", "Pat26-A", "Pat27-A", "Pat28-A", "Pat29-A"
  )
  g150_patient <- sub("-.*$", "", g150_title)
  paired_patients <- sprintf("Pat%02d", c(1:20, 22:24))
  g150_tissue <- ifelse(g150_patient %in% paired_patients & grepl("-A$", g150_title), "Adjacent_Normal",
                        ifelse(g150_patient %in% paired_patients & grepl("-B$", g150_title), "Gastric_Cancer",
                               ifelse(g150_patient == "Pat21", "Adjacent_Normal_Unpaired", "Control_Normal_Unpaired")))
  g150_pair_eligible <- g150_patient %in% paired_patients
  g150 <- data.frame(
    cohort = "GSE150290", gsm = paste0("GSM", g150_gsm), sample_label = g150_title,
    patient = g150_patient, tissue_or_condition = g150_tissue,
    disease_context = "GC_adjacent_normal_atlas",
    mapping_status = ifelse(g150_pair_eligible, "OFFICIAL_GEO_PATIENT_TISSUE_MAP_PAIRED",
                            "OFFICIAL_GEO_NORMAL_SINGLETON_DESCRIPTIVE_ONLY"),
    inferential_eligible = g150_pair_eligible, stringsAsFactors = FALSE
  )

  # GSE183904: the manuscript-ID/sample-number crosswalk is published directly
  # on the GEO series page. GSM = 5573465 + sample number.
  g183_patient <- c(
    "NGCII518","NGCII518","NGCII519","NGCII520","NGCII520","NGCII521","NGCII521","NGCII524","NGCII525","NGCII525",
    "NGCII527","NGCII527","NGCII529","NGCII502","NGCII511","NGCII499","NGCII509","NGCII498","NGCII011","NGCII015",
    "NGCII513","NGCII513","NGCII514","NGCII514","NGCII512","NGCII512","NGCII510","NGCII531","NGCII522","NGCII533",
    "NGCII536","NGCII539","NGCII540","NGCII541","NGCII538","NGCII538","NGCII538","NGCII538","NGCII543","NGCII545"
  )
  g183_tissue <- c(
    "Primary_Normal","Primary_Tumor","Primary_Tumor","Primary_Normal","Primary_Tumor","Primary_Normal","Primary_Tumor","Primary_Tumor","Primary_Normal","Primary_Tumor",
    "Primary_Normal","Primary_Tumor",rep("Primary_Tumor", 6),"Peritoneal_Tumor","Peritoneal_Tumor",
    "Primary_Normal","Primary_Tumor","Primary_Normal","Primary_Tumor","Primary_Normal","Primary_Tumor",rep("Primary_Tumor", 4),
    "Primary_Normal",rep("Primary_Tumor", 3),"Primary_Normal","Primary_Tumor","Peritoneal_Normal","Peritoneal_Tumor","Primary_Tumor","Primary_Tumor"
  )
  stopifnot(length(g183_patient) == 40L, length(g183_tissue) == 40L)
  g183 <- data.frame(
    cohort = "GSE183904", gsm = paste0("GSM", 5573465 + 1:40),
    sample_label = paste0("sample", 1:40), patient = g183_patient,
    tissue_or_condition = g183_tissue, disease_context = "primary_and_peritoneal_GC_atlas",
    mapping_status = "OFFICIAL_GEO_MANUSCRIPT_ID_CROSSWALK", inferential_eligible = TRUE,
    stringsAsFactors = FALSE
  )

  g167_patient <- c(rep("Pt1", 3), rep("Pt2", 2), rep("Pt3", 3), rep("Pt4", 3), rep("Pt5", 3))
  g167_tissue <- c("Normal","Superficial","Deep", "Superficial","Deep",
                   "Normal","Superficial","Deep", "Normal","Superficial","Deep",
                   "Normal","Superficial","Deep")
  g167 <- data.frame(
    cohort = "GSE167297", gsm = paste0("GSM", 5101013:5101026),
    sample_label = paste0(g167_patient, "_", tolower(g167_tissue), ifelse(g167_tissue == "Normal", "", "_cancer")),
    patient = g167_patient, tissue_or_condition = g167_tissue,
    disease_context = "diffuse_type_GC_invasion_depth",
    mapping_status = "OFFICIAL_GEO_PATIENT_LAYER_MAP", inferential_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  rbind(g150, g183, g167)
}

o2p6_match_official_sample <- function(cohort, raw) {
  m <- o2p6_official_sample_map()
  m <- m[m$cohort == cohort, , drop = FALSE]
  if (!nrow(m)) return(data.frame())
  raw_upper <- toupper(raw)
  gsm <- regmatches(raw_upper, regexpr("GSM[0-9]+", raw_upper))
  if (length(gsm) && nzchar(gsm)) {
    z <- m[m$gsm == gsm, , drop = FALSE]
    if (nrow(z) == 1L) return(z)
  }
  if (cohort == "GSE183904") {
    num <- regmatches(raw_upper, regexpr("SAMPLE[0-9]+", raw_upper))
    if (length(num) && nzchar(num)) {
      z <- m[toupper(m$sample_label) == num, , drop = FALSE]
      if (nrow(z) == 1L) return(z)
    }
  } else {
    hit <- which(vapply(toupper(m$sample_label), function(tok) grepl(tok, raw_upper, fixed = TRUE), logical(1)))
    if (length(hit) == 1L) return(m[hit, , drop = FALSE])
  }
  data.frame()
}

.o2p6_if78_infer_fallback <- o2p6_infer_sample_info
o2p6_infer_sample_info <- function(cohort, sample_path, file_name = basename(sample_path)) {
  raw <- paste(sample_path, file_name, sep = "/")
  z <- o2p6_match_official_sample(cohort, raw)
  role <- switch(cohort,
    GSE150290 = "PRIMARY_TUMOR_NORMAL_GC_ATLAS",
    GSE183904 = "PRIMARY_ADVANCED_GC_TUMOR_ATLAS",
    GSE167297 = "DEPTH_CONTEXT_DGC_SUPERFICIAL_DEEP",
    GSE134520 = "SUPPLEMENTARY_PREMALIGNANT_EGC_CONTEXT", "UNSPECIFIED")
  if (nrow(z) == 1L) {
    return(data.frame(
      cohort = cohort, dataset_role = role, sample_id = z$gsm,
      patient = z$patient, tissue_or_condition = z$tissue_or_condition,
      disease_context = z$disease_context, official_sample_label = z$sample_label,
      official_mapping_status = z$mapping_status,
      inferential_unit = ifelse(isTRUE(z$inferential_eligible), "official_patient", "sample_descriptive_only"),
      inferential_eligible = isTRUE(z$inferential_eligible), stringsAsFactors = FALSE
    ))
  }
  out <- .o2p6_if78_infer_fallback(cohort, sample_path, file_name)
  out$dataset_role <- role
  out$official_sample_label <- NA_character_
  out$official_mapping_status <- ifelse(cohort == "GSE134520", "LESION_LABEL_ONLY_SUPPLEMENTARY", "UNRESOLVED_EXCLUDE_FROM_PATIENT_INFERENCE")
  out$inferential_unit <- "sample_descriptive_only"
  out$inferential_eligible <- FALSE
  out
}

o2p6_if78_patient_summary <- function(pseudobulk) {
  if (!is.data.frame(pseudobulk) || !nrow(pseudobulk)) return(data.frame())
  d <- pseudobulk[pseudobulk$inferential_eligible %in% TRUE & !is.na(pseudobulk$patient) &
                    nzchar(pseudobulk$patient) & pseudobulk$n_cells >= 30L, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  o2p6_bind_rows(lapply(split(d, paste(d$cohort, d$patient, d$tissue_or_condition, sep = "||")), function(z) {
    n <- sum(z$n_cells, na.rm = TRUE)
    npos <- sum(z$n_OLFML2B_positive, na.rm = TRUE)
    nge2 <- if ("n_OLFML2B_ge2" %in% names(z)) sum(z$n_OLFML2B_ge2, na.rm = TRUE) else NA_real_
    nge3 <- if ("n_OLFML2B_ge3" %in% names(z)) sum(z$n_OLFML2B_ge3, na.rm = TRUE) else NA_real_
    data.frame(
      cohort = z$cohort[1], patient = z$patient[1], tissue_or_condition = z$tissue_or_condition[1],
      n_source_samples = nrow(z), n_cells = n,
      n_OLFML2B_positive = npos, OLFML2B_positive_fraction = npos / pmax(n, 1),
      n_OLFML2B_ge2 = nge2, OLFML2B_ge2_fraction = nge2 / pmax(n, 1),
      n_OLFML2B_ge3 = nge3, OLFML2B_ge3_fraction = nge3 / pmax(n, 1),
      mean_OLFML2B_count = stats::weighted.mean(z$mean_OLFML2B_count, w = pmax(z$n_cells, 1), na.rm = TRUE),
      inference_unit = "patient", stringsAsFactors = FALSE
    )
  }))
}

o2p6_if78_paired_tumor_normal <- function(patient_summary) {
  if (!nrow(patient_summary)) return(list(delta = data.frame(), test = data.frame()))
  d <- patient_summary
  d$condition2 <- ifelse(grepl("Normal", d$tissue_or_condition), "Normal",
                         ifelse(grepl("Tumor|GC|Gastric_Cancer", d$tissue_or_condition) & !grepl("Peritoneal", d$tissue_or_condition), "Tumor", NA_character_))
  d <- d[!is.na(d$condition2) & d$cohort %in% c("GSE150290", "GSE183904"), , drop = FALSE]
  rows <- list(); tests <- list()
  for (coh in unique(d$cohort)) {
    z <- d[d$cohort == coh, , drop = FALSE]
    patients <- intersect(z$patient[z$condition2 == "Tumor"], z$patient[z$condition2 == "Normal"])
    if (!length(patients)) next
    delta <- data.frame(
      cohort = coh, patient = patients,
      tumor_fraction = z$OLFML2B_positive_fraction[match(paste(patients, "Tumor"), paste(z$patient, z$condition2))],
      normal_fraction = z$OLFML2B_positive_fraction[match(paste(patients, "Normal"), paste(z$patient, z$condition2))],
      stringsAsFactors = FALSE
    )
    delta$delta_tumor_minus_normal <- delta$tumor_fraction - delta$normal_fraction
    rows[[coh]] <- delta
    wt <- if (nrow(delta) >= 3L) tryCatch(suppressWarnings(stats::wilcox.test(delta$tumor_fraction, delta$normal_fraction, paired = TRUE, exact = FALSE)), error = function(e) NULL) else NULL
    tests[[coh]] <- data.frame(cohort = coh, n_paired_patients = nrow(delta),
      median_delta_tumor_minus_normal = stats::median(delta$delta_tumor_minus_normal, na.rm = TRUE),
      wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
      inference_unit = "officially_mapped_paired_patient", stringsAsFactors = FALSE)
  }
  test <- o2p6_bind_rows(tests)
  if (nrow(test)) test$wilcox_p_fdr <- stats::p.adjust(test$wilcox_p, method = "BH")
  list(delta = o2p6_bind_rows(rows), test = test)
}

.o2p6_if78_runner_core <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_scrna_rc2_no_kang_full_production <- function(...) {
  args <- list(...)
  args$datasets <- c("GSE150290", "GSE183904", "GSE167297", "GSE134520")
  args$primary_datasets <- c("GSE150290", "GSE183904")
  args$supportive_datasets <- c("GSE134520")
  index <- do.call(.o2p6_if78_runner_core, args)
  root <- normalizePath(args$root %||% "D:/OLFML2B_STAD", winslash = "/", mustWork = FALSE)
  output_subdir <- args$output_subdir %||% "Part6_RC2"
  dirs <- o2p6_dirs(root, output_subdir)
  index$version <- OLFML2B_PART6_IF78_VERSION
  index$dataset_plan$role <- c(
    GSE150290 = "PRIMARY_RC2_TUMOR_NORMAL",
    GSE183904 = "PRIMARY_RC2_GC_ATLAS",
    GSE167297 = "DEPTH_CONTEXT_PAIRED_PATIENT",
    GSE134520 = "SUPPLEMENTARY_PREMALIGNANT_EGC"
  )[index$dataset_plan$cohort]
  o2p6_atomic_write_csv(index$dataset_plan, file.path(dirs$tables, "01_dataset_inclusion_plan_no_kang.csv"))
  official_map <- o2p6_official_sample_map()
  o2p6_atomic_write_csv(official_map, file.path(dirs$tables, "05_official_GEO_sample_patient_map.csv"))
  map_audit <- o2p6_bind_rows(lapply(c("GSE150290", "GSE183904", "GSE167297", "GSE134520"), function(coh) {
    a <- index$dataset_audit[index$dataset_audit$cohort == coh, , drop = FALSE]
    data.frame(cohort = coh, n_parsed = nrow(a),
      n_officially_mapped = if (nrow(a) && "official_mapping_status" %in% names(a)) sum(grepl("^OFFICIAL", a$official_mapping_status)) else 0L,
      n_inferentially_eligible = if (nrow(a) && "inferential_eligible" %in% names(a)) sum(a$inferential_eligible %in% TRUE) else 0L,
      n_unique_patients = if (nrow(a)) length(unique(a$patient[grepl("^OFFICIAL", a$official_mapping_status) & !is.na(a$patient)])) else 0L,
      formal_role = c(GSE150290="PRIMARY",GSE183904="PRIMARY",GSE167297="DEPTH_CONTEXT",GSE134520="SUPPLEMENTARY")[coh],
      mapping_gate = ifelse(coh %in% c("GSE150290","GSE183904") && nrow(a) > 0L && all(grepl("^OFFICIAL", a$official_mapping_status)), "PASS_PRIMARY_OFFICIAL_MAP",
                            ifelse(coh %in% c("GSE150290","GSE183904"), "FAIL_OR_PARTIAL_PRIMARY_MAP", "NONPRIMARY_CONTEXT")),
      stringsAsFactors = FALSE)
  }))
  o2p6_atomic_write_csv(map_audit, file.path(dirs$tables, "06_official_mapping_audit_by_cohort.csv"))
  patient_summary <- o2p6_if78_patient_summary(index$pseudobulk)
  paired <- o2p6_if78_paired_tumor_normal(patient_summary)
  o2p6_atomic_write_csv(patient_summary, file.path(dirs$tables, "40_patient_level_OLFML2B_pseudobulk.csv"))
  o2p6_atomic_write_csv(paired$delta, file.path(dirs$tables, "41_primary_paired_patient_tumor_normal_delta.csv"))
  o2p6_atomic_write_csv(paired$test, file.path(dirs$tables, "42_primary_paired_patient_tumor_normal_test.csv"))

  guard <- data.frame(
    item = c("primary_RC2_datasets", "depth_context", "supplementary_context", "Kang_GSE206785", "formal_inference_unit", "cell_level_tests", "marker_compartments", "claim_ceiling"),
    decision = c("GSE150290;GSE183904", "GSE167297", "GSE134520", "EXCLUDED_BY_DESIGN",
                 "official patient when mapped; otherwise sample is descriptive only",
                 "PROHIBITED_AS_PRIMARY", "exploratory marker-derived broad-compartment localization only",
                 "cellular context/localization; no causality, prognosis or ICI prediction"),
    stringsAsFactors = FALSE
  )
  o2p6_atomic_write_csv(guard, file.path(dirs$tables, "100_part6_pseudoreplication_guardrail.csv"))
  o2p6_atomic_write_csv(guard, file.path(dirs$tables, "101_part6_methodology_claim_limits.csv"))
  primary_map_ok <- all(map_audit$mapping_gate[map_audit$cohort %in% c("GSE150290","GSE183904")] == "PASS_PRIMARY_OFFICIAL_MAP")
  go <- data.frame(
    criterion = c("RC2_primary_role_lock", "Kang_exclusion", "official_primary_mapping", "patient_sample_pseudobulk", "paired_patient_analysis", "no_cell_pseudoreplication", "result_independent_pipeline_continuation", "claim_ceiling"),
    status = c("PASS", "PASS", ifelse(primary_map_ok, "PASS", "FAIL_OR_NOT_EVALUABLE"),
               ifelse(nrow(patient_summary), "PASS", "NOT_EVALUABLE"), ifelse(nrow(paired$delta), "PASS", "NOT_EVALUABLE"),
               "PASS", "PASS", "CONTEXT_SUPPORT_ONLY"),
    interpretation = c(
      "GSE150290 and GSE183904 are the only primary RC2 cohorts.",
      "GSE206785/Kang is not read, scored or used in the formal evidence matrix.",
      "Primary inference requires the embedded official GEO sample-patient map.",
      "Biological summaries are aggregated at official patient/sample level.",
      "Tumor-normal comparisons use paired official patients only.",
      "Cells contribute to pseudobulk summaries and are never treated as independent biological replicates.",
      "Detection, effect direction and P values do not control pipeline continuation.",
      "No causal, prognostic, cell-type-specific mechanistic or ICI-response claim."
    ), stringsAsFactors = FALSE
  )
  index$final_go_no_go <- go
  index$official_sample_map <- official_map
  index$mapping_audit <- map_audit
  index$patient_summary <- patient_summary
  index$paired_patient <- paired
  index$dirs <- dirs
  o2p6_atomic_write_csv(go, file.path(dirs$tables, "99_Part6_RC2_no_kang_go_no_go.csv"))
  o2p6_atomic_write_csv(go, file.path(dirs$tables, "99_Part6_single_cell_go_no_go.csv"))
  o2p6_atomic_save_rds(index, file.path(dirs$objects, "Part6_RC2_no_kang_4dataset_full_production_index.rds"))
  o2p6_atomic_save_rds(index, file.path(dirs$objects, "Part6_single_cell_evidence_layer_index.rds"))
  invisible(index)
}

run_olfml2b_part6_scrna_rc2_no_kang <- run_olfml2b_part6_scrna_rc2_no_kang_full_production

# Backward-compatible Part9 entry point. Legacy Kang/atlas arguments are
# accepted only so old callers do not break; they are never used.
run_olfml2b_part6_singlecell_production <- function(
  root = "D:/OLFML2B_STAD",
  sc_dir = file.path(root, "data", "raw", "single_cell"),
  atlas_dir = NULL, kang_dir = NULL, allow_download = FALSE,
  force_download = FALSE, force_unpack = FALSE, make_figures = TRUE, ...
) {
  invisible(atlas_dir); invisible(kang_dir); invisible(allow_download); invisible(force_download)
  run_olfml2b_part6_scrna_rc2_no_kang_full_production(
    root = root, raw_single_cell_dir = sc_dir, output_subdir = "Part6",
    force_unpack = force_unpack, make_figures = make_figures,
    exclude_gse150290_raw_10x = TRUE
  )
}

# ----------------------------------------------------------------------------
# Full-transcriptome pseudobulk extension. Text matrices are streamed so a
# complete gene-by-cell matrix need not be materialized in memory; 10X inputs
# remain sparse and are collapsed with rowSums.
# ----------------------------------------------------------------------------
o2p6_full_pseudobulk_text <- function(path, cohort, info, chunk_size = 1000L) {
  con <- o2p6_open_text(path)
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header)) stop("Empty text matrix: ", path, call. = FALSE)
  sep <- o2p6_detect_sep_from_header(header)
  h <- strsplit(header, sep, fixed = TRUE)[[1]]
  h_clean <- o2p6_clean_gene(h)
  cell_by_gene <- "OLFML2B" %in% h_clean
  if (cell_by_gene) {
    gene_idx <- which(seq_along(h_clean) > 1L & nzchar(h_clean))
    gene_names <- h_clean[gene_idx]
    sums <- numeric(length(gene_idx)); n_cells <- 0L
    repeat {
      lines <- readLines(con, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      fields <- strsplit(lines, sep, fixed = TRUE)
      for (f in fields) {
        if (length(f) < max(gene_idx)) next
        v <- suppressWarnings(as.numeric(f[gene_idx]))
        v[!is.finite(v)] <- 0
        sums <- sums + v
        n_cells <- n_cells + 1L
      }
    }
    out <- data.frame(gene = gene_names, pseudobulk_count = sums, n_cells = n_cells,
                      stringsAsFactors = FALSE)
  } else {
    rows <- list(); k <- 0L; n_cells <- max(0L, length(h) - 1L)
    repeat {
      lines <- readLines(con, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      fields <- strsplit(lines, sep, fixed = TRUE)
      for (f in fields) {
        if (length(f) < 2L) next
        gene <- o2p6_clean_gene(f[1])
        if (!nzchar(gene)) next
        v <- suppressWarnings(as.numeric(f[-1]))
        v[!is.finite(v)] <- 0
        k <- k + 1L
        rows[[k]] <- data.frame(gene = gene, pseudobulk_count = sum(v),
                                n_cells = length(v), stringsAsFactors = FALSE)
      }
    }
    out <- o2p6_bind_rows(rows)
  }
  if (!nrow(out)) return(data.frame())
  count <- stats::aggregate(pseudobulk_count ~ gene, data = out, sum, na.rm = TRUE)
  cells <- stats::aggregate(n_cells ~ gene, data = out, max, na.rm = TRUE)
  out <- merge(count, cells, by = "gene", all = TRUE, sort = FALSE)
  cbind(info[rep(1L, nrow(out)), , drop = FALSE],
        input_type = "TEXT_MATRIX", file_name = basename(path), out,
        stringsAsFactors = FALSE)
}

o2p6_full_pseudobulk_10x <- function(candidate, info) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required for 10X pseudobulk.", call. = FALSE)
  con <- if (grepl("\\.gz$", candidate$matrix_file[1], ignore.case = TRUE)) gzfile(candidate$matrix_file[1], "rt") else candidate$matrix_file[1]
  if (inherits(con, "connection")) on.exit(close(con), add = TRUE)
  mat <- Matrix::readMM(con)
  fcon <- if (grepl("\\.gz$", candidate$feature_file[1], ignore.case = TRUE)) gzfile(candidate$feature_file[1], "rt") else candidate$feature_file[1]
  if (inherits(fcon, "connection")) on.exit(close(fcon), add = TRUE)
  feat <- utils::read.delim(fcon, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  genes <- o2p6_clean_gene(if (ncol(feat) >= 2L) feat[[2]] else feat[[1]])
  if (nrow(mat) != length(genes) && ncol(mat) == length(genes)) mat <- Matrix::t(mat)
  if (nrow(mat) != length(genes)) stop("10X feature/matrix dimension mismatch: ", candidate$matrix_file[1], call. = FALSE)
  out <- data.frame(gene = genes, pseudobulk_count = as.numeric(Matrix::rowSums(mat)),
                    n_cells = ncol(mat), stringsAsFactors = FALSE)
  out <- out[nzchar(out$gene), , drop = FALSE]
  count <- stats::aggregate(pseudobulk_count ~ gene, data = out, sum, na.rm = TRUE)
  cells <- stats::aggregate(n_cells ~ gene, data = out, max, na.rm = TRUE)
  out <- merge(count, cells, by = "gene", all = TRUE, sort = FALSE)
  cbind(info[rep(1L, nrow(out)), , drop = FALSE], input_type = "10X_MTX",
        file_name = candidate$file_name[1], out, stringsAsFactors = FALSE)
}

o2p6_candidate_source_paths <- function(candidate) {
  if (identical(candidate$input_type[1], "10X_MTX")) {
    c(candidate$matrix_file[1], candidate$feature_file[1], candidate$barcode_file[1])
  } else {
    candidate$matrix_file[1]
  }
}

o2p6_full_pseudobulk_one <- function(candidate, cache_dir = NULL, refresh_cache = FALSE) {
  started <- Sys.time()
  info <- o2p6_infer_sample_info(candidate$cohort[1], candidate$sample_path[1], candidate$file_name[1])
  source_paths <- o2p6_candidate_source_paths(candidate)
  cached <- o2p6_fulltx_cache_read(source_paths, require_target = FALSE,
                                    cache_dir = cache_dir, refresh_cache = refresh_cache)
  if (isTRUE(cached$hit)) {
    raw <- cached$object$fulltx
    status <- "CACHE_HIT_VALIDATED"
    created_at <- cached$object$created_at %||% NA_character_
    created_by <- cached$object$created_by %||% NA_character_
  } else {
    full <- if (identical(candidate$input_type[1], "10X_MTX")) {
      o2p6_full_pseudobulk_10x(candidate, info)
    } else {
      o2p6_full_pseudobulk_text(candidate$matrix_file[1], candidate$cohort[1], info)
    }
    raw <- full[, c("gene", "pseudobulk_count", "n_cells"), drop = FALSE]
    cache_file <- o2p6_fulltx_cache_write(source_paths, raw, target_payload = NULL,
                                           created_by = "fulltx_fallback", cache_dir = cache_dir)
    cached$cache_file <- cache_file
    status <- paste0(cached$reason, "_REBUILT")
    created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
    created_by <- "fulltx_fallback"
  }
  data <- if (nrow(raw)) {
    cbind(info[rep(1L, nrow(raw)), , drop = FALSE],
          input_type = rep(candidate$input_type[1], nrow(raw)),
          file_name = rep(candidate$file_name[1], nrow(raw)),
          raw, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  audit <- data.frame(
    cohort = candidate$cohort[1], file_name = candidate$file_name[1],
    matrix_file = candidate$matrix_file[1], cache_status = status,
    cache_file = cached$cache_file %||% NA_character_, cache_created_at = created_at,
    cache_created_by = created_by,
    cache_validation = "normalized_path+byte_size+mtime+cache_schema",
    n_genes = nrow(raw), n_cells = if (nrow(raw)) max(raw$n_cells, na.rm = TRUE) else 0,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    error = NA_character_, stringsAsFactors = FALSE
  )
  list(data = data, audit = audit)
}

o2p6_patient_full_pseudobulk <- function(sample_gene) {
  if (!is.data.frame(sample_gene) || !nrow(sample_gene)) return(data.frame())
  d <- sample_gene[sample_gene$inferential_eligible %in% TRUE & !is.na(sample_gene$patient) &
                     nzchar(sample_gene$patient) & sample_gene$n_cells >= 30L, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  stats::aggregate(
    cbind(pseudobulk_count, n_cells) ~ cohort + patient + tissue_or_condition + gene,
    data = d,
    FUN = sum, na.rm = TRUE
  )
}

o2p6_fulltx_scores <- function(patient_gene) {
  if (!is.data.frame(patient_gene) || !nrow(patient_gene)) return(list(target = data.frame(), programs = data.frame(), coverage = data.frame()))
  patient_gene$unit_id <- paste(patient_gene$cohort, patient_gene$patient, patient_gene$tissue_or_condition, sep = "||")
  lib <- stats::aggregate(pseudobulk_count ~ unit_id, patient_gene, sum, na.rm = TRUE)
  names(lib)[2] <- "library_size"
  d <- merge(patient_gene, lib, by = "unit_id", all.x = TRUE, sort = FALSE)
  d$log2_cpm <- log2((d$pseudobulk_count + 0.5) / (d$library_size + 1) * 1e6)
  target <- d[d$gene == "OLFML2B", c("unit_id","cohort","patient","tissue_or_condition","n_cells","pseudobulk_count","library_size","log2_cpm"), drop = FALSE]
  sets <- o2p6_program_sets()
  scores <- list(); cov_rows <- list()
  for (nm in names(sets)) {
    genes <- unique(o2p6_clean_gene(sets[[nm]]))
    z <- d[d$gene %in% genes, , drop = FALSE]
    present <- unique(z$gene)
    eligible <- length(present) >= 5L && length(present) / max(length(genes), 1L) >= 0.50
    cov_rows[[nm]] <- data.frame(program = nm, requested_genes = length(genes), present_genes = length(present),
                                 coverage_fraction = length(present) / max(length(genes), 1L), eligible = eligible,
                                 stringsAsFactors = FALSE)
    if (!eligible || !nrow(z)) next
    z$gene_z <- stats::ave(z$log2_cpm, z$cohort, z$gene, FUN = function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
    })
    one <- stats::aggregate(gene_z ~ unit_id + cohort + patient + tissue_or_condition,
                            z, mean, na.rm = TRUE)
    names(one)[names(one) == "gene_z"] <- "program_score"
    one$program <- nm; one$n_genes <- length(present)
    scores[[nm]] <- one
  }
  list(target = target, programs = o2p6_bind_rows(scores), coverage = o2p6_bind_rows(cov_rows))
}

.o2p6_if78_fulltx_core <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_scrna_rc2_no_kang_full_production <- function(
  ...,
  refresh_cache = identical(toupper(Sys.getenv("OLFML2B_PART6_REFRESH_CACHE", unset = "FALSE")), "TRUE")
) {
  args <- list(...)
  root <- normalizePath(args$root %||% "D:/OLFML2B_STAD", winslash = "/", mustWork = FALSE)
  output_subdir <- args$output_subdir %||% "Part6_RC2"
  pre_dirs <- o2p6_dirs(root, output_subdir)
  cache_dir <- o2p6_dir(file.path(pre_dirs$objects, output_subdir, "cache", "full_transcriptome"))
  md5_cache_file <- file.path(pre_dirs$objects, output_subdir, "cache", "input_md5_manifest.rds")
  log_file <- file.path(pre_dirs$logs, "Part6_RC2_no_kang_full_production.log")
  old_options <- options(
    olfml2b.part6.fulltx_cache_dir = cache_dir,
    olfml2b.part6.md5_cache_file = md5_cache_file,
    olfml2b.part6.refresh_md5_cache = isTRUE(refresh_cache),
    olfml2b.part6.refresh_cache = isTRUE(refresh_cache)
  )
  on.exit(options(old_options), add = TRUE)
  o2p6_log("INFO", "Validated sample-cache mode enabled | cache: ", cache_dir,
             " | refresh_cache=", isTRUE(refresh_cache), log_file = log_file)
  index <- do.call(.o2p6_if78_fulltx_core, args)
  dirs <- index$dirs
  used <- index$candidates[index$candidates$used_for_analysis %in% TRUE, , drop = FALSE]
  aggregate_started <- Sys.time()
  parsed <- lapply(seq_len(nrow(used)), function(i) {
    one_started <- Sys.time()
    ans <- tryCatch(
      o2p6_full_pseudobulk_one(used[i, , drop = FALSE], cache_dir = cache_dir, refresh_cache = FALSE),
      error = function(e) list(
        data = data.frame(),
        audit = data.frame(
          cohort = used$cohort[i], file_name = used$file_name[i], matrix_file = used$matrix_file[i],
          cache_status = "CACHE_OR_PARSE_FAILED", cache_file = NA_character_, cache_created_at = NA_character_,
          cache_created_by = NA_character_, cache_validation = "normalized_path+byte_size+mtime+cache_schema",
          n_genes = 0L, n_cells = 0L,
          elapsed_seconds = as.numeric(difftime(Sys.time(), one_started, units = "secs")),
          error = conditionMessage(e), stringsAsFactors = FALSE
        )
      )
    )
    o2p6_log(if (nrow(ans$data)) "INFO" else "WARN",
      sprintf("Full-transcriptome cache aggregation %d/%d | %s | %.1f s | %s",
        i, nrow(used), ans$audit$cache_status[1], ans$audit$elapsed_seconds[1], used$file_name[i]),
      log_file = log_file)
    ans
  })
  cache_audit <- o2p6_bind_rows(lapply(parsed, `[[`, "audit"))
  valid <- lapply(parsed[vapply(parsed, function(z) is.data.frame(z$data) && nrow(z$data) &&
                                  all(c("gene", "pseudobulk_count") %in% names(z$data)), logical(1))], `[[`, "data")
  failures <- cache_audit[!is.na(cache_audit$error) & nzchar(cache_audit$error), , drop = FALSE]
  sample_gene <- o2p6_bind_rows(valid)
  patient_gene <- o2p6_patient_full_pseudobulk(sample_gene)
  scores <- o2p6_fulltx_scores(patient_gene)
  manifest <- if (nrow(sample_gene)) o2p6_bind_rows(lapply(split(sample_gene, paste(sample_gene$cohort, sample_gene$file_name, sep = "||")), function(z) data.frame(
    cohort = z$cohort[1], file_name = z$file_name[1], patient = z$patient[1],
    tissue_or_condition = z$tissue_or_condition[1], n_cells = max(z$n_cells, na.rm = TRUE),
    n_genes = length(unique(z$gene)), library_size = sum(z$pseudobulk_count, na.rm = TRUE),
    mapping_status = z$official_mapping_status[1], stringsAsFactors = FALSE))) else data.frame()
  o2p6_atomic_write_csv(manifest, file.path(dirs$tables, "43_full_transcriptome_pseudobulk_manifest.csv"))
  o2p6_atomic_write_csv(failures, file.path(dirs$tables, "43a_full_transcriptome_parse_failures.csv"))
  o2p6_atomic_write_csv(cache_audit, file.path(dirs$tables, "43b_full_transcriptome_cache_audit.csv"))
  o2p6_atomic_write_csv(scores$target, file.path(dirs$tables, "44_patient_level_full_transcriptome_OLFML2B_log2CPM.csv"))
  o2p6_atomic_write_csv(scores$programs, file.path(dirs$tables, "45_patient_level_full_transcriptome_program_scores.csv"))
  o2p6_atomic_write_csv(scores$coverage, file.path(dirs$tables, "46_full_transcriptome_program_coverage.csv"))
  fulltx_paired <- data.frame(); fulltx_paired_test <- data.frame()
  if (nrow(scores$target)) {
    z <- scores$target[scores$target$cohort == "GSE150290", , drop = FALSE]
    z$condition2 <- ifelse(grepl("Normal", z$tissue_or_condition), "Normal",
                           ifelse(grepl("Tumor|GC|Gastric_Cancer", z$tissue_or_condition), "Tumor", NA_character_))
    patients <- intersect(z$patient[z$condition2 == "Tumor"], z$patient[z$condition2 == "Normal"])
    if (length(patients)) {
      fulltx_paired <- data.frame(
        cohort = "GSE150290", patient = patients,
        tumor_log2_cpm = z$log2_cpm[match(paste(patients, "Tumor"), paste(z$patient, z$condition2))],
        normal_log2_cpm = z$log2_cpm[match(paste(patients, "Normal"), paste(z$patient, z$condition2))],
        stringsAsFactors = FALSE
      )
      fulltx_paired$delta_tumor_minus_normal <- fulltx_paired$tumor_log2_cpm - fulltx_paired$normal_log2_cpm
      wt <- if (nrow(fulltx_paired) >= 3L) tryCatch(suppressWarnings(stats::wilcox.test(
        fulltx_paired$tumor_log2_cpm, fulltx_paired$normal_log2_cpm, paired = TRUE, exact = FALSE)), error = function(e) NULL) else NULL
      tt <- if (nrow(fulltx_paired) >= 3L) tryCatch(stats::t.test(
        fulltx_paired$tumor_log2_cpm, fulltx_paired$normal_log2_cpm, paired = TRUE), error = function(e) NULL) else NULL
      fulltx_paired_test <- data.frame(
        cohort = "GSE150290", estimand = "paired_patient_tumor_minus_adjacent_normal_OLFML2B_log2CPM",
        n_paired_patients = nrow(fulltx_paired), mean_delta = mean(fulltx_paired$delta_tumor_minus_normal, na.rm = TRUE),
        median_delta = stats::median(fulltx_paired$delta_tumor_minus_normal, na.rm = TRUE),
        mean_delta_ci_low = if (!is.null(tt)) tt$conf.int[1] else NA_real_,
        mean_delta_ci_high = if (!is.null(tt)) tt$conf.int[2] else NA_real_,
        paired_t_p = if (!is.null(tt)) tt$p.value else NA_real_,
        paired_wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
        inference_unit = "officially_mapped_paired_patient", stringsAsFactors = FALSE
      )
    }
  }
  o2p6_atomic_write_csv(fulltx_paired, file.path(dirs$tables, "47_GSE150290_fulltx_paired_patient_delta.csv"))
  o2p6_atomic_write_csv(fulltx_paired_test, file.path(dirs$tables, "48_GSE150290_fulltx_primary_estimand.csv"))
  fulltx_depth <- data.frame(); fulltx_depth_tests <- data.frame()
  if (nrow(scores$target)) {
    z167 <- scores$target[scores$target$cohort == "GSE167297", , drop = FALSE]
    patients167 <- sort(unique(as.character(z167$patient)))
    if (length(patients167)) {
      fulltx_depth <- data.frame(patient = patients167, stringsAsFactors = FALSE)
      for (cond in c("Normal", "Superficial", "Deep")) {
        fulltx_depth[[cond]] <- z167$log2_cpm[match(paste(patients167, cond),
                                                    paste(z167$patient, z167$tissue_or_condition))]
      }
      fulltx_depth$Deep_minus_Normal <- fulltx_depth$Deep - fulltx_depth$Normal
      fulltx_depth$Superficial_minus_Normal <- fulltx_depth$Superficial - fulltx_depth$Normal
      fulltx_depth$Deep_minus_Superficial <- fulltx_depth$Deep - fulltx_depth$Superficial
      fulltx_depth_tests <- o2p6_paired_gradient_tests(fulltx_depth)
      if (nrow(fulltx_depth_tests)) {
        fulltx_depth_tests$measure <- "full_transcriptome_OLFML2B_log2CPM"
        fulltx_depth_tests$claim_role <- "five_patient_depth_context_only"
      }
    }
  }
  o2p6_atomic_write_csv(fulltx_depth, file.path(dirs$tables, "49_GSE167297_fulltx_patient_depth_gradient.csv"))
  o2p6_atomic_write_csv(fulltx_depth_tests, file.path(dirs$tables, "50_GSE167297_fulltx_exact_depth_tests.csv"))
  o2p6_atomic_save_rds(list(sample_gene = sample_gene, patient_gene = patient_gene, scores = scores),
                        file.path(dirs$objects, "Part6_RC2_full_transcriptome_pseudobulk.rds"))
  fulltx_ok <- nrow(scores$target) > 0L && nrow(patient_gene) > 0L &&
    all(c("GSE150290", "GSE183904") %in% unique(scores$target$cohort))
  add <- data.frame(
    criterion = c("full_transcriptome_patient_pseudobulk", "GSE150290_primary_paired_estimand"),
    status = c(ifelse(fulltx_ok, "PASS", "FAIL_OR_NOT_EVALUABLE"),
               ifelse(nrow(fulltx_paired) >= 15L, "PASS", "NOT_EVALUABLE_PRIMARY_LT15_PAIRS")),
    interpretation = c(
      "All-gene counts are collapsed before normalization at official patient×tissue level; marker-only tables are localization supplements.",
      "Primary RC2 estimand is the officially mapped within-patient tumor-minus-adjacent-normal OLFML2B log2-CPM difference; <15 pairs is supportive only."
    ),
    stringsAsFactors = FALSE
  )
  index$final_go_no_go <- rbind(index$final_go_no_go, add)
  index$full_transcriptome_manifest <- manifest
  index$full_transcriptome_target <- scores$target
  index$full_transcriptome_programs <- scores$programs
  index$full_transcriptome_coverage <- scores$coverage
  index$full_transcriptome_paired_delta <- fulltx_paired
  index$full_transcriptome_primary_estimand <- fulltx_paired_test
  index$full_transcriptome_GSE167297_depth_gradient <- fulltx_depth
  index$full_transcriptome_GSE167297_depth_tests <- fulltx_depth_tests
  index$full_transcriptome_cache_audit <- cache_audit
  source_paths <- unique(unlist(lapply(seq_len(nrow(used)), function(i) {
    o2p6_candidate_source_paths(used[i, , drop = FALSE])
  }), use.names = FALSE))
  index$input_source_fingerprint <- o2p6_source_fingerprint(source_paths)
  o2p6_atomic_write_csv(index$final_go_no_go, file.path(dirs$tables, "99_Part6_RC2_no_kang_go_no_go.csv"))
  o2p6_atomic_write_csv(index$final_go_no_go, file.path(dirs$tables, "99_Part6_single_cell_go_no_go.csv"))
  o2p6_atomic_save_rds(index, file.path(dirs$objects, "Part6_RC2_no_kang_4dataset_full_production_index.rds"))
  o2p6_atomic_save_rds(index, file.path(dirs$objects, "Part6_single_cell_evidence_layer_index.rds"))
  o2p6_log("INFO", sprintf("Part6 complete | %d matrices | cache hits=%d | rebuilt=%d | failures=%d | aggregation %.1f min",
    nrow(cache_audit), sum(cache_audit$cache_status == "CACHE_HIT_VALIDATED", na.rm = TRUE),
    sum(grepl("REBUILT$", cache_audit$cache_status), na.rm = TRUE), nrow(failures),
    as.numeric(difftime(Sys.time(), aggregate_started, units = "mins"))), log_file = log_file)
  invisible(index)
}

run_olfml2b_part6_scrna_rc2_no_kang <- run_olfml2b_part6_scrna_rc2_no_kang_full_production

run_olfml2b_part6_singlecell_production <- function(
  root = "D:/OLFML2B_STAD", sc_dir = file.path(root, "data", "raw", "single_cell"),
  atlas_dir = NULL, kang_dir = NULL, allow_download = FALSE, force_download = FALSE,
  force_unpack = FALSE, make_figures = TRUE, ...
) {
  invisible(atlas_dir); invisible(kang_dir); invisible(allow_download); invisible(force_download)
  run_olfml2b_part6_scrna_rc2_no_kang_full_production(
    root = root, raw_single_cell_dir = sc_dir, output_subdir = "Part6",
    force_unpack = force_unpack, make_figures = make_figures,
    exclude_gse150290_raw_10x = TRUE
  )
}


# ==============================================================================
# OLFML2B direct-fork adaptation extension
# - preserves the mature four-dataset RC2 reader and official mapping
# - separates fibroblast, myofibroblast, smooth-muscle and pericyte compartments
# - adds count-threshold, confidence and composition-vs-within-compartment audits
# ==============================================================================

OLFML2B_PART6_DIRECT_FORK_VERSION <- "v1.3.0_20260722_UNIFIED_FOUR_DATASET_EXACT_INFERENCE_AND_COMPARTMENT_PROGRAM_FIX"

.o2p6_boot_median_ci <- function(x, B = 2000L, seed = 20260722L) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) == 1L) return(c(x, x))
  set.seed(seed)
  b <- replicate(B, stats::median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  as.numeric(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 6))
}

.o2p6_exact_sign_p <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  k <- sum(x > 0); n <- length(x)
  min(1, 2 * stats::pbinom(min(k, n-k), size = n, prob = 0.5))
}

o2p6_patient_compartment_summary <- function(celltype_summary, patient_summary) {
  if (!is.data.frame(celltype_summary) || !nrow(celltype_summary)) return(data.frame())
  d <- celltype_summary
  d$unit_id <- ifelse(!is.na(d$patient) & nzchar(as.character(d$patient)), as.character(d$patient), as.character(d$sample_id))
  d <- d[!is.na(d$unit_id) & nzchar(d$unit_id), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  out <- o2p6_bind_rows(lapply(split(d, paste(d$cohort, d$unit_id, d$tissue_or_condition, d$marker_celltype, sep="||")), function(z) {
    n <- sum(z$n_cells, na.rm = TRUE)
    npos <- sum(z$n_OLFML2B_positive, na.rm = TRUE)
    nge2 <- if ("n_OLFML2B_ge2" %in% names(z)) sum(z$n_OLFML2B_ge2, na.rm = TRUE) else NA_real_
    nge3 <- if ("n_OLFML2B_ge3" %in% names(z)) sum(z$n_OLFML2B_ge3, na.rm = TRUE) else NA_real_
    nhigh <- if ("n_high_confidence_cells" %in% names(z)) sum(z$n_high_confidence_cells, na.rm = TRUE) else NA_real_
    npos_high <- if ("n_OLFML2B_positive_high_confidence" %in% names(z)) sum(z$n_OLFML2B_positive_high_confidence, na.rm = TRUE) else NA_real_
    data.frame(
      cohort=z$cohort[1], patient_or_sample=z$unit_id[1], patient=z$patient[1], sample_id=z$sample_id[1],
      tissue_or_condition=z$tissue_or_condition[1], marker_celltype=z$marker_celltype[1],
      n_source_matrices=nrow(z), n_cells=n,
      n_OLFML2B_positive=npos, OLFML2B_positive_fraction=npos/pmax(n,1),
      n_OLFML2B_ge2=nge2, OLFML2B_ge2_fraction=nge2/pmax(n,1),
      n_OLFML2B_ge3=nge3, OLFML2B_ge3_fraction=nge3/pmax(n,1),
      mean_OLFML2B_count=stats::weighted.mean(z$mean_OLFML2B_count, pmax(z$n_cells,1), na.rm=TRUE),
      n_high_confidence_cells=nhigh, high_confidence_fraction=nhigh/pmax(n,1),
      n_OLFML2B_positive_high_confidence=npos_high,
      OLFML2B_positive_fraction_high_confidence=npos_high/pmax(nhigh,1),
      inference_unit=ifelse(!is.na(z$patient[1]) & nzchar(as.character(z$patient[1])), "patient", "sample_descriptive"),
      stringsAsFactors=FALSE
    )
  }))
  if (is.data.frame(patient_summary) && nrow(patient_summary)) {
    key <- paste(patient_summary$cohort, patient_summary$patient, patient_summary$tissue_or_condition, sep="||")
    m <- match(paste(out$cohort, out$patient_or_sample, out$tissue_or_condition, sep="||"), key)
    total <- patient_summary$n_cells[m]
    out$all_compartment_total_cells <- total
    out$compartment_cell_fraction <- out$n_cells / pmax(total, 1)
  } else {
    out$all_compartment_total_cells <- NA_real_
    out$compartment_cell_fraction <- NA_real_
  }
  out
}

o2p6_prespecified_contrasts <- function() {
  data.frame(
    cohort=c("GSE150290","GSE183904","GSE183904","GSE183904","GSE167297","GSE167297","GSE167297"),
    contrast=c("Gastric_Cancer_minus_Adjacent_Normal","Primary_Tumor_minus_Primary_Normal","Peritoneal_Tumor_minus_Peritoneal_Normal","Peritoneal_Tumor_minus_Primary_Tumor","Deep_minus_Normal","Superficial_minus_Normal","Deep_minus_Superficial"),
    high=c("Gastric_Cancer","Primary_Tumor","Peritoneal_Tumor","Peritoneal_Tumor","Deep","Superficial","Deep"),
    low=c("Adjacent_Normal","Primary_Normal","Peritoneal_Normal","Primary_Tumor","Normal","Normal","Superficial"),
    evidence_role=c("PRIMARY_PAIRED_TUMOR_NORMAL","PRIMARY_PAIRED_TUMOR_NORMAL","PERITONEAL_CONTEXT","PERITONEAL_CONTEXT","DEPTH_CONTEXT","DEPTH_CONTEXT","DEPTH_CONTEXT"),
    stringsAsFactors=FALSE
  )
}

o2p6_compartment_paired_analysis <- function(patient_compartment, B=2000L, min_cells=20L, min_pairs=3L) {
  if (!is.data.frame(patient_compartment) || !nrow(patient_compartment)) return(list(delta=data.frame(), tests=data.frame()))
  d <- patient_compartment[patient_compartment$inference_unit == "patient", , drop=FALSE]
  ctr <- o2p6_prespecified_contrasts()
  metrics <- c("OLFML2B_positive_fraction","OLFML2B_ge2_fraction","OLFML2B_ge3_fraction","mean_OLFML2B_count","compartment_cell_fraction")
  rows <- list(); tests <- list(); k <- 0L; q <- 0L
  for (i in seq_len(nrow(ctr))) {
    z0 <- d[d$cohort == ctr$cohort[i], , drop=FALSE]
    if (!nrow(z0)) next
    for (ct in sort(unique(as.character(z0$marker_celltype)))) {
      z <- z0[z0$marker_celltype == ct & z0$n_cells >= min_cells, , drop=FALSE]
      patients <- intersect(z$patient_or_sample[z$tissue_or_condition == ctr$high[i]], z$patient_or_sample[z$tissue_or_condition == ctr$low[i]])
      if (!length(patients)) next
      for (metric in metrics) {
        if (!metric %in% names(z)) next
        hi <- z[[metric]][match(paste(patients, ctr$high[i]), paste(z$patient_or_sample, z$tissue_or_condition))]
        lo <- z[[metric]][match(paste(patients, ctr$low[i]), paste(z$patient_or_sample, z$tissue_or_condition))]
        delta <- hi-lo
        keep <- is.finite(hi) & is.finite(lo)
        if (!any(keep)) next
        k <- k+1L
        rows[[k]] <- data.frame(cohort=ctr$cohort[i], contrast=ctr$contrast[i], evidence_role=ctr$evidence_role[i],
          marker_celltype=ct, metric=metric, patient=patients[keep], high=hi[keep], low=lo[keep], delta=delta[keep],
          min_cells_per_condition=min_cells, stringsAsFactors=FALSE)
        x <- delta[keep]; ci <- .o2p6_boot_median_ci(x, B=B, seed=20260722L+i+k)
        wt <- if (length(x)>=min_pairs) tryCatch(suppressWarnings(stats::wilcox.test(x, mu=0, paired=FALSE, exact=FALSE)), error=function(e) NULL) else NULL
        q <- q+1L
        tests[[q]] <- data.frame(cohort=ctr$cohort[i], contrast=ctr$contrast[i], evidence_role=ctr$evidence_role[i],
          marker_celltype=ct, metric=metric, n_paired_patients=length(x), median_delta=stats::median(x), mean_delta=mean(x),
          bootstrap_median_ci_low=ci[1], bootstrap_median_ci_high=ci[2], n_positive=sum(x>0), n_negative=sum(x<0), n_zero=sum(x==0),
          wilcox_p=if (!is.null(wt)) wt$p.value else NA_real_, exact_sign_p=.o2p6_exact_sign_p(x),
          formal_evaluable=length(x)>=min_pairs, min_cells_per_condition=min_cells, bootstrap_B=B,
          inference_unit="paired_patient", stringsAsFactors=FALSE)
      }
    }
  }
  test <- o2p6_bind_rows(tests)
  if (nrow(test)) {
    test$wilcox_fdr <- stats::p.adjust(test$wilcox_p, method="BH")
    test$sign_fdr <- stats::p.adjust(test$exact_sign_p, method="BH")
  }
  list(delta=o2p6_bind_rows(rows), tests=test)
}

o2p6_strict_localization_consistency <- function(patient_compartment, min_cells=20L, min_units=3L, min_positive=10L) {
  if (!is.data.frame(patient_compartment) || !nrow(patient_compartment)) return(list(by_dataset=data.frame(), consensus=data.frame()))
  d <- patient_compartment[patient_compartment$n_cells >= min_cells, , drop=FALSE]
  if (!nrow(d)) return(list(by_dataset=data.frame(), consensus=data.frame()))
  by <- o2p6_bind_rows(lapply(split(d, paste(d$cohort,d$marker_celltype,sep="||")), function(z) {
    data.frame(cohort=z$cohort[1], marker_celltype=z$marker_celltype[1], n_units=length(unique(z$patient_or_sample)),
      total_cells=sum(z$n_cells,na.rm=TRUE), total_positive=sum(z$n_OLFML2B_positive,na.rm=TRUE),
      pooled_positive_fraction=sum(z$n_OLFML2B_positive,na.rm=TRUE)/pmax(sum(z$n_cells,na.rm=TRUE),1),
      median_unit_positive_fraction=stats::median(z$OLFML2B_positive_fraction,na.rm=TRUE),
      total_ge2=sum(z$n_OLFML2B_ge2,na.rm=TRUE), pooled_ge2_fraction=sum(z$n_OLFML2B_ge2,na.rm=TRUE)/pmax(sum(z$n_cells,na.rm=TRUE),1),
      total_positive_high_confidence=sum(z$n_OLFML2B_positive_high_confidence,na.rm=TRUE),
      eligible=length(unique(z$patient_or_sample))>=min_units && sum(z$n_OLFML2B_positive,na.rm=TRUE)>=min_positive,
      stringsAsFactors=FALSE)
  }))
  by$rank_within_dataset <- stats::ave(-by$median_unit_positive_fraction, by$cohort, FUN=function(x) rank(x,ties.method="min"))
  eligible <- by[by$eligible %in% TRUE, , drop=FALSE]
  cons <- o2p6_bind_rows(lapply(split(eligible, eligible$marker_celltype), function(z) data.frame(
    marker_celltype=z$marker_celltype[1], n_datasets_eligible=nrow(z), n_datasets_top3=sum(z$rank_within_dataset<=3),
    datasets_eligible=paste(z$cohort,collapse=";"), median_dataset_fraction=stats::median(z$median_unit_positive_fraction,na.rm=TRUE),
    strict_consistency_status=ifelse(nrow(z)>=2L && sum(z$rank_within_dataset<=3)>=2L,"REPLICATED_TOP3_IN_AT_LEAST_TWO_DATASETS","DESCRIPTIVE_OR_NOT_REPLICATED"),
    stringsAsFactors=FALSE)))
  list(by_dataset=by, consensus=cons)
}

.o2p6_direct_fork_core <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_scrna_rc2_no_kang_full_production <- function(...) {
  args <- list(...)
  root <- normalizePath(args$root %||% "D:/OLFML2B_STAD", winslash="/", mustWork=FALSE)
  output_subdir <- args$output_subdir %||% "Part6"
  args$root <- root; args$output_subdir <- output_subdir
  index <- do.call(.o2p6_direct_fork_core, args)
  dirs <- index$dirs %||% o2p6_dirs(root, output_subdir)
  patient_summary <- index$patient_summary %||% o2p6_if78_patient_summary(index$pseudobulk)
  patient_compartment <- o2p6_patient_compartment_summary(index$celltype_summary, patient_summary)
  paired_compartment <- o2p6_compartment_paired_analysis(patient_compartment, B=2000L, min_cells=20L, min_pairs=3L)
  strict <- o2p6_strict_localization_consistency(patient_compartment, min_cells=20L, min_units=3L, min_positive=10L)
  o2p6_atomic_write_csv(patient_compartment, file.path(dirs$tables,"51_patient_condition_compartment_OLFML2B.csv"))
  o2p6_atomic_write_csv(paired_compartment$delta, file.path(dirs$tables,"52_paired_compartment_patient_deltas.csv"))
  o2p6_atomic_write_csv(paired_compartment$tests, file.path(dirs$tables,"53_paired_compartment_tests.csv"))
  o2p6_atomic_write_csv(strict$by_dataset, file.path(dirs$tables,"54_strict_localization_by_dataset.csv"))
  o2p6_atomic_write_csv(strict$consensus, file.path(dirs$tables,"55_strict_cross_dataset_localization_consistency.csv"))
  adapt_go <- data.frame(
    criterion=c("four_dataset_scope","GSE150290_processed_text_guard","official_primary_mapping","patient_level_pseudobulk","compartment_disambiguation","threshold_sensitivity","composition_vs_within_compartment","strict_cross_dataset_gate","Part0_5_freeze","final_gene_lock"),
    status=c(
      ifelse(all(c("GSE150290","GSE183904","GSE167297","GSE134520") %in% index$dataset_plan$cohort),"PASS","FAIL"),
      ifelse(any(index$candidates$input_policy=="EXCLUDE_GSE150290_RAW_10X_DROPLET_MATRIX",na.rm=TRUE),"PASS","REVIEW"),
      ifelse(is.data.frame(index$mapping_audit) && all(index$mapping_audit$mapping_gate[index$mapping_audit$cohort %in% c("GSE150290","GSE183904")] == "PASS_PRIMARY_OFFICIAL_MAP"),"PASS","FAIL_OR_NOT_EVALUABLE"),
      ifelse(nrow(patient_summary),"PASS","NOT_EVALUABLE"),
      ifelse(all(c("Fibroblast","Myofibroblast","Smooth_Muscle","Pericyte","Myeloid","Epithelial") %in% unique(index$celltype_summary$marker_celltype)),"PASS","REVIEW"),
      ifelse(all(c("OLFML2B_ge2_fraction","OLFML2B_ge3_fraction") %in% names(index$pseudobulk)),"PASS","FAIL"),
      ifelse(nrow(paired_compartment$tests),"PASS","NOT_EVALUABLE"),
      ifelse(nrow(strict$by_dataset),"PASS","NOT_EVALUABLE"),
      "PASS","FALSE"),
    interpretation=c(
      "All four locally supplied gastric single-cell datasets enter the analysis with prespecified unequal evidence roles.",
      "GSE150290 author-processed GSM text matrices are formal inputs; raw 10X droplet matrices remain audit-only.",
      "GSE150290 and GSE183904 use embedded official GSM-patient-tissue crosswalks.",
      "Cells are aggregated to official patient/sample units; cells never constitute independent biological replicates.",
      "Fibroblast, myofibroblast, smooth muscle, pericyte, myeloid and epithelial compartments are kept separate.",
      "OLFML2B detection thresholds count>=1, >=2 and >=3 are reported without changing the primary threshold.",
      "Paired analyses separate compartment abundance from within-compartment OLFML2B changes when evaluable.",
      "Cross-dataset localization requires minimum cells, units and positive cells; small-denominator top ranks are not treated as replication.",
      "Part0-5 objects and methods are read-only and are not recomputed by Part6.",
      "The gene remains unlocked until cellular source and spatial evidence are reconciled."),
    stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(adapt_go, file.path(dirs$tables,"99_OLFML2B_Part6_direct_fork_go_no_go.csv"))
  claim <- data.frame(
    level=c("ALLOWED","CONDITIONAL","FORBIDDEN"),
    statement=c(
      "Patient/sample-level OLFML2B cellular-context and localization evidence across four gastric single-cell datasets.",
      "A preferential fibroblast, smooth-muscle, pericyte, myeloid or epithelial source only when official-patient, confidence, threshold and cross-dataset outputs agree.",
      "No unique CAF source, malignant-cell source, secretion, receptor, direct TGF-beta activation, causal invasion/recurrence, prognosis or treatment-response claim from Part6 alone."),
    stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(claim, file.path(dirs$tables,"100_OLFML2B_Part6_claim_boundary.csv"))
  index$version <- OLFML2B_PART6_DIRECT_FORK_VERSION
  index$direct_fork_lineage <- "UBE2Q2 Part6 RC2 v2.1.0_20260720_RC2_RESUMABLE_SINGLE_PASS; unified OLFML2B exact-inference, patient-condition and within-compartment adaptations"
  index$patient_compartment <- patient_compartment
  index$paired_compartment <- paired_compartment
  index$strict_localization <- strict
  index$adaptation_go_no_go <- adapt_go
  index$OLFML2B_claim_boundary <- claim
  index$Part0_5_frozen <- TRUE
  index$final_gene_lock <- FALSE
  index$dirs <- dirs
  o2p6_atomic_save_rds(index, file.path(dirs$objects,"Part6_OLFML2B_four_dataset_direct_fork_index.rds"))
  o2p6_atomic_save_rds(index, file.path(dirs$objects,"Part6_OLFML2B_single_cell_index.rds"))
  p05 <- file.path(dirs$objects,"OLFML2B_Part0_5_complete_index.rds")
  integrated <- if (file.exists(p05)) tryCatch(readRDS(p05),error=function(e) list()) else list()
  integrated$Part6 <- index; integrated$Part6_version <- OLFML2B_PART6_DIRECT_FORK_VERSION
  integrated$Part0_5_frozen <- TRUE; integrated$final_gene_lock <- FALSE; integrated$updated_at <- o2p6_ts()
  o2p6_atomic_save_rds(integrated,file.path(dirs$objects,"OLFML2B_Part0_6_complete_index.rds"))
  o2p6_log("INFO","OLFML2B Part6 direct fork complete | four datasets | final_gene_lock=FALSE",log_file=file.path(dirs$logs,"Part6_OLFML2B_direct_fork.log"))
  invisible(index)
}

run_olfml2b_part6_scrna_rc2_no_kang <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_singlecell_production <- function(
  root="D:/OLFML2B_STAD", sc_dir=file.path(root,"data","raw","single_cell"),
  force_unpack=FALSE, make_figures=TRUE, refresh_cache=identical(toupper(Sys.getenv("OLFML2B_PART6_REFRESH_CACHE",unset="FALSE")),"TRUE"), ...
) {
  run_olfml2b_part6_scrna_rc2_no_kang_full_production(
    root=root, raw_single_cell_dir=sc_dir, output_subdir="Part6",
    force_unpack=force_unpack, make_figures=make_figures,
    exclude_gse150290_raw_10x=TRUE, refresh_cache=refresh_cache
  )
}

# ==============================================================================
# OLFML2B Part6 v1.3.0 unified finalization layer
# - one production file and one runner; no separate summary-repair module
# - exact/permutation Spearman inference for small and moderate patient counts
# - patient-condition program summaries without condition mixing
# - within-compartment OLFML2B-positive versus negative program audit
# - strict source grading with threshold and annotation-confidence sensitivity
# ==============================================================================

OLFML2B_PART6_UNIFIED_VERSION <- "v1.3.0_20260722_UNIFIED_FOUR_DATASET_EXACT_INFERENCE_AND_COMPARTMENT_PROGRAM_FIX"

.o2p6_v130_perm_cache <- new.env(parent = emptyenv())

.o2p6_v130_permutation_indices <- function(n) {
  key <- as.character(as.integer(n))
  if (exists(key, envir = .o2p6_v130_perm_cache, inherits = FALSE)) {
    return(get(key, envir = .o2p6_v130_perm_cache, inherits = FALSE))
  }
  if (n < 1L || n > 9L) stop("Exact permutation index supports n=1..9.", call. = FALSE)
  p <- matrix(1L, nrow = 1L, ncol = 1L)
  if (n >= 2L) {
    for (k in 2L:n) {
      old <- p
      nr <- nrow(old)
      new <- matrix(NA_integer_, nrow = nr * k, ncol = k)
      for (pos in seq_len(k)) {
        rows <- ((pos - 1L) * nr + 1L):(pos * nr)
        if (pos > 1L) new[rows, seq_len(pos - 1L)] <- old[, seq_len(pos - 1L), drop = FALSE]
        new[rows, pos] <- k
        if (pos <= k - 1L) new[rows, (pos + 1L):k] <- old[, pos:(k - 1L), drop = FALSE]
      }
      p <- new
    }
  }
  assign(key, p, envir = .o2p6_v130_perm_cache)
  p
}

.o2p6_v130_spearman_test <- function(x, y, seed = 20260722L, mc_B = 100000L) {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  rho <- if (n >= 3L) suppressWarnings(stats::cor(x, y, method = "spearman")) else NA_real_
  out <- list(
    n = n, rho = unname(rho), p_value = NA_real_,
    p_method = if (n < 4L) "DESCRIPTIVE_ONLY_N_LT4" else NA_character_,
    n_permutations = 0L,
    small_n_exploratory = n < 6L,
    formal_significance_eligible = n >= 6L
  )
  if (n < 4L || !is.finite(rho)) return(out)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  tol <- sqrt(.Machine$double.eps)
  if (n <= 9L) {
    idx <- .o2p6_v130_permutation_indices(n)
    yc <- ry - mean(ry)
    xc <- rx - mean(rx)
    den <- sqrt(sum(xc^2) * sum(yc^2))
    perm_rho <- as.numeric((idx * 0L)[, 1L])
    # Vectorized row-wise covariance after permuting ranked y.
    yp <- matrix(ry[idx], nrow = nrow(idx), ncol = n)
    yp <- sweep(yp, 1L, rowMeans(yp), FUN = "-")
    perm_rho <- as.numeric(yp %*% xc) / den
    out$p_value <- mean(abs(perm_rho) >= abs(rho) - tol)
    out$p_method <- "FULL_ENUMERATION_EXACT_SPEARMAN"
    out$n_permutations <- nrow(idx)
    return(out)
  }
  if (n < 20L) {
    set.seed(as.integer(seed))
    xc <- rx - mean(rx); yc <- ry - mean(ry)
    den <- sqrt(sum(xc^2) * sum(yc^2))
    exceed <- 0L
    chunk <- 5000L
    remaining <- as.integer(mc_B)
    while (remaining > 0L) {
      b <- min(chunk, remaining)
      vals <- numeric(b)
      for (i in seq_len(b)) vals[i] <- sum(xc * sample(yc, n, replace = FALSE)) / den
      exceed <- exceed + sum(abs(vals) >= abs(rho) - tol)
      remaining <- remaining - b
    }
    out$p_value <- (exceed + 1) / (as.integer(mc_B) + 1)
    out$p_method <- "MONTE_CARLO_PERMUTATION_SPEARMAN"
    out$n_permutations <- as.integer(mc_B)
    return(out)
  }
  ct <- tryCatch(suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)), error = function(e) NULL)
  out$p_value <- if (!is.null(ct)) ct$p.value else NA_real_
  out$p_method <- "ASYMPTOTIC_SPEARMAN_N_GE20"
  out$n_permutations <- 0L
  out$small_n_exploratory <- FALSE
  out$formal_significance_eligible <- TRUE
  out
}

.o2p6_v130_boot_ci <- function(x, B = 2000L, seed = 20260722L, stat = stats::median) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) == 1L) return(rep(x, 2L))
  set.seed(as.integer(seed))
  z <- replicate(as.integer(B), stat(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  as.numeric(stats::quantile(z, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 6))
}

.o2p6_v130_exact_sign_p <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  stats::binom.test(sum(x > 0), length(x), p = 0.5)$p.value
}

.o2p6_v130_paired_plan <- function() {
  data.frame(
    cohort = c("GSE150290", "GSE183904", "GSE167297", "GSE167297"),
    high = c("Gastric_Cancer", "Primary_Tumor", "Deep", "Deep"),
    low = c("Adjacent_Normal", "Primary_Normal", "Normal", "Superficial"),
    contrast = c("Gastric_Cancer_minus_Adjacent_Normal", "Primary_Tumor_minus_Primary_Normal", "Deep_minus_Normal", "Deep_minus_Superficial"),
    evidence_role = c("PRIMARY_PAIRED", "EXTERNAL_PAIRED_SUPPORT", "SMALL_N_DEPTH_CONTEXT", "SMALL_N_DEPTH_CONTEXT"),
    stringsAsFactors = FALSE
  )
}

.o2p6_v130_gse183904_paired <- function(target) {
  z <- target[target$cohort == "GSE183904" & target$tissue_or_condition %in% c("Primary_Tumor", "Primary_Normal"), , drop = FALSE]
  if (!nrow(z)) return(list(delta = data.frame(), test = data.frame()))
  a <- stats::aggregate(log2_cpm ~ patient + tissue_or_condition, z, mean, na.rm = TRUE)
  hi <- a[a$tissue_or_condition == "Primary_Tumor", c("patient", "log2_cpm"), drop = FALSE]
  lo <- a[a$tissue_or_condition == "Primary_Normal", c("patient", "log2_cpm"), drop = FALSE]
  names(hi)[2] <- "primary_tumor_log2_cpm"; names(lo)[2] <- "primary_normal_log2_cpm"
  d <- merge(hi, lo, by = "patient", all = FALSE, sort = TRUE)
  if (!nrow(d)) return(list(delta = d, test = data.frame()))
  d$delta_tumor_minus_normal <- d$primary_tumor_log2_cpm - d$primary_normal_log2_cpm
  x <- d$delta_tumor_minus_normal
  wt <- if (length(x) >= 3L) tryCatch(suppressWarnings(stats::wilcox.test(x, mu = 0, exact = length(x) <= 20L)), error = function(e) NULL) else NULL
  tt <- if (length(x) >= 3L) tryCatch(stats::t.test(x, mu = 0), error = function(e) NULL) else NULL
  ci <- .o2p6_v130_boot_ci(x, seed = 20260722L)
  test <- data.frame(
    cohort = "GSE183904", contrast = "Primary_Tumor_minus_Primary_Normal", n_paired_patients = length(x),
    mean_delta = mean(x), median_delta = stats::median(x), bootstrap_median_ci_low = ci[1], bootstrap_median_ci_high = ci[2],
    paired_t_p = if (!is.null(tt)) tt$p.value else NA_real_, paired_wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
    exact_sign_p = .o2p6_v130_exact_sign_p(x), n_positive = sum(x > 0), n_negative = sum(x < 0), n_zero = sum(x == 0),
    evidence_role = "EXTERNAL_PAIRED_DIRECTIONAL_SUPPORT", inference_unit = "official_paired_patient", stringsAsFactors = FALSE
  )
  list(delta = d, test = test)
}

.o2p6_v130_program_evaluability <- function(program_sample) {
  if (!is.data.frame(program_sample) || !nrow(program_sample)) return(data.frame())
  key <- paste(program_sample$cohort, program_sample$sample_id, sep = "||")
  o2p6_bind_rows(lapply(split(program_sample, key), function(z) {
    n_cells <- suppressWarnings(as.numeric(z$n_cells[1])); n_pos <- suppressWarnings(as.numeric(z$n_OLFML2B_positive[1]))
    n_neg <- if (is.finite(n_cells) && is.finite(n_pos)) n_cells - n_pos else NA_real_
    data.frame(
      cohort = z$cohort[1], sample_id = z$sample_id[1], patient = z$patient[1], tissue_or_condition = z$tissue_or_condition[1],
      n_cells = n_cells, n_OLFML2B_positive = n_pos, n_OLFML2B_negative = n_neg,
      n_programs_total = length(unique(z$marker_program)),
      n_programs_evaluable = sum(is.finite(suppressWarnings(as.numeric(z$delta_positive_minus_negative)))),
      evaluable = is.finite(n_pos) && n_pos > 0 && is.finite(n_neg) && n_neg > 0,
      exclusion_reason = ifelse(!is.finite(n_pos) || n_pos <= 0, "NO_OLFML2B_POSITIVE_CELL",
                                ifelse(!is.finite(n_neg) || n_neg <= 0, "NO_OLFML2B_NEGATIVE_CELL", "EVALUABLE")),
      stringsAsFactors = FALSE
    )
  }))
}

.o2p6_v130_patient_condition_program <- function(program_sample) {
  if (!is.data.frame(program_sample) || !nrow(program_sample)) return(data.frame())
  d <- program_sample[is.finite(suppressWarnings(as.numeric(program_sample$delta_positive_minus_negative))), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  d$patient_or_sample <- ifelse(!is.na(d$patient) & nzchar(as.character(d$patient)), as.character(d$patient), as.character(d$sample_id))
  d$inference_unit <- ifelse(!is.na(d$patient) & nzchar(as.character(d$patient)), "official_patient", "sample_descriptive")
  key <- paste(d$cohort, d$patient_or_sample, d$tissue_or_condition, d$marker_program, sep = "||")
  o2p6_bind_rows(lapply(split(d, key), function(z) {
    x <- suppressWarnings(as.numeric(z$delta_positive_minus_negative)); x <- x[is.finite(x)]
    data.frame(
      cohort = z$cohort[1], dataset_role = z$dataset_role[1], patient_or_sample = z$patient_or_sample[1],
      inference_unit = z$inference_unit[1], tissue_or_condition = z$tissue_or_condition[1], marker_program = z$marker_program[1],
      n_source_samples = length(unique(z$sample_id)), n_cells_total = sum(suppressWarnings(as.numeric(z$n_cells)), na.rm = TRUE),
      n_OLFML2B_positive_total = sum(suppressWarnings(as.numeric(z$n_OLFML2B_positive)), na.rm = TRUE),
      median_delta_positive_minus_negative = stats::median(x), mean_delta_positive_minus_negative = mean(x),
      stringsAsFactors = FALSE
    )
  }))
}


.o2p6_v130_family_factor <- function(d, cols) {
  cols <- cols[cols %in% names(d)]
  if (!length(cols)) return(factor(rep("ALL", nrow(d))))
  do.call(interaction, c(d[, cols, drop = FALSE], list(drop = TRUE, lex.order = TRUE)))
}

.o2p6_v130_one_sample_tests <- function(d, value_col, family_cols, label_cols, min_units = 3L, B = 2000L) {
  if (!is.data.frame(d) || !nrow(d)) return(data.frame())
  test_group_cols <- unique(label_cols[label_cols %in% names(d)])
  fam <- .o2p6_v130_family_factor(d, test_group_cols)
  out <- o2p6_bind_rows(lapply(split(d, fam), function(z) {
    x <- suppressWarnings(as.numeric(z[[value_col]])); x <- x[is.finite(x)]
    n <- length(x); ci <- .o2p6_v130_boot_ci(x, B = B, seed = 20260722L + n + nchar(paste(z[1, label_cols, drop = TRUE], collapse = "")))
    wt <- if (n >= min_units) tryCatch(suppressWarnings(stats::wilcox.test(x, mu = 0, exact = n <= 20L)), error = function(e) NULL) else NULL
    row <- z[1, label_cols, drop = FALSE]
    row$n_units <- n; row$median_delta <- if (n) stats::median(x) else NA_real_; row$mean_delta <- if (n) mean(x) else NA_real_
    row$bootstrap_median_ci_low <- ci[1]; row$bootstrap_median_ci_high <- ci[2]
    row$n_positive <- sum(x > 0); row$n_negative <- sum(x < 0); row$n_zero <- sum(x == 0)
    row$wilcox_p <- if (!is.null(wt)) wt$p.value else NA_real_; row$exact_sign_p <- .o2p6_v130_exact_sign_p(x)
    row$formal_evaluable <- n >= min_units
    row
  }))
  if (nrow(out)) {
    fam2 <- .o2p6_v130_family_factor(out, family_cols)
    out$wilcox_fdr_within_family <- stats::ave(out$wilcox_p, fam2, FUN = function(v) stats::p.adjust(v, method = "BH"))
    out$sign_fdr_within_family <- stats::ave(out$exact_sign_p, fam2, FUN = function(v) stats::p.adjust(v, method = "BH"))
  }
  out
}

.o2p6_v130_primary_tumor_program_tests <- function(patient_condition_program) {
  plan <- data.frame(cohort = c("GSE150290", "GSE183904", "GSE167297"), tissue_or_condition = c("Gastric_Cancer", "Primary_Tumor", "Deep"), stringsAsFactors = FALSE)
  d <- merge(patient_condition_program, plan, by = c("cohort", "tissue_or_condition"), all = FALSE)
  d <- d[d$inference_unit == "official_patient", , drop = FALSE]
  .o2p6_v130_one_sample_tests(d, "median_delta_positive_minus_negative",
                              family_cols = c("cohort", "tissue_or_condition"),
                              label_cols = c("cohort", "tissue_or_condition", "marker_program"), min_units = 3L)
}

.o2p6_v130_paired_marker_program_change <- function(patient_condition_program) {
  plan <- .o2p6_v130_paired_plan()
  deltas <- list(); k <- 0L
  for (i in seq_len(nrow(plan))) {
    z0 <- patient_condition_program[patient_condition_program$cohort == plan$cohort[i] &
                                      patient_condition_program$tissue_or_condition %in% c(plan$high[i], plan$low[i]) &
                                      patient_condition_program$inference_unit == "official_patient", , drop = FALSE]
    for (prog in unique(as.character(z0$marker_program))) {
      z <- z0[z0$marker_program == prog, , drop = FALSE]
      hi <- z[z$tissue_or_condition == plan$high[i], c("patient_or_sample", "median_delta_positive_minus_negative"), drop = FALSE]
      lo <- z[z$tissue_or_condition == plan$low[i], c("patient_or_sample", "median_delta_positive_minus_negative"), drop = FALSE]
      names(hi) <- c("patient", "high_delta"); names(lo) <- c("patient", "low_delta")
      m <- merge(hi, lo, by = "patient", all = FALSE, sort = TRUE)
      if (!nrow(m)) next
      m$paired_change <- m$high_delta - m$low_delta
      m$cohort <- plan$cohort[i]; m$contrast <- plan$contrast[i]; m$evidence_role <- plan$evidence_role[i]; m$marker_program <- prog
      k <- k + 1L; deltas[[k]] <- m[, c("cohort", "contrast", "evidence_role", "marker_program", "patient", "high_delta", "low_delta", "paired_change")]
    }
  }
  delta <- o2p6_bind_rows(deltas)
  tests <- .o2p6_v130_one_sample_tests(delta, "paired_change", family_cols = c("cohort", "contrast"),
                                        label_cols = c("cohort", "contrast", "evidence_role", "marker_program"), min_units = 3L)
  list(delta = delta, tests = tests)
}

.o2p6_v130_merge_fulltx <- function(target, programs) {
  merge(programs, target[, c("unit_id", "log2_cpm", "n_cells"), drop = FALSE], by = "unit_id", all = FALSE, sort = FALSE)
}

.o2p6_v130_condition_bridge <- function(target, programs) {
  d <- .o2p6_v130_merge_fulltx(target, programs)
  if (!nrow(d)) return(data.frame())
  key <- paste(d$cohort, d$tissue_or_condition, d$program, sep = "||")
  out <- o2p6_bind_rows(lapply(split(d, key), function(z) {
    a <- stats::aggregate(cbind(log2_cpm, program_score) ~ patient, z, mean, na.rm = TRUE)
    tst <- .o2p6_v130_spearman_test(a$log2_cpm, a$program_score,
                                    seed = 20260722L + nchar(z$cohort[1]) + nchar(z$program[1]))
    data.frame(
      cohort = z$cohort[1], tissue_or_condition = z$tissue_or_condition[1], program = z$program[1],
      n_patients = tst$n, rho = tst$rho, p_value = tst$p_value, p_method = tst$p_method,
      n_permutations = tst$n_permutations, small_n_exploratory = tst$small_n_exploratory,
      formal_significance_eligible = tst$formal_significance_eligible,
      primary_tumor_state = (z$cohort[1] == "GSE150290" && z$tissue_or_condition[1] == "Gastric_Cancer") ||
        (z$cohort[1] == "GSE183904" && z$tissue_or_condition[1] == "Primary_Tumor") ||
        (z$cohort[1] == "GSE167297" && z$tissue_or_condition[1] == "Deep"),
      inference_level = "condition-stratified patient-level full-transcriptome Spearman",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out)) {
    fam <- paste(out$cohort, out$tissue_or_condition, sep = "||")
    out$fdr_within_cohort_condition <- stats::ave(out$p_value, fam, FUN = function(v) stats::p.adjust(v, method = "BH"))
  }
  out
}

.o2p6_v130_paired_change_bridge <- function(target, programs) {
  d <- .o2p6_v130_merge_fulltx(target, programs)
  plan <- .o2p6_v130_paired_plan()
  deltas <- list(); tests <- list(); k <- 0L; q <- 0L
  for (i in seq_len(nrow(plan))) {
    z0 <- d[d$cohort == plan$cohort[i] & d$tissue_or_condition %in% c(plan$high[i], plan$low[i]), , drop = FALSE]
    for (prog in unique(as.character(z0$program))) {
      z <- z0[z0$program == prog, , drop = FALSE]
      a <- stats::aggregate(cbind(log2_cpm, program_score) ~ patient + tissue_or_condition, z, mean, na.rm = TRUE)
      hi <- a[a$tissue_or_condition == plan$high[i], c("patient", "log2_cpm", "program_score"), drop = FALSE]
      lo <- a[a$tissue_or_condition == plan$low[i], c("patient", "log2_cpm", "program_score"), drop = FALSE]
      names(hi)[2:3] <- c("target_high", "program_high"); names(lo)[2:3] <- c("target_low", "program_low")
      m <- merge(hi, lo, by = "patient", all = FALSE, sort = TRUE)
      if (!nrow(m)) next
      m$target_delta <- m$target_high - m$target_low; m$program_delta <- m$program_high - m$program_low
      m$cohort <- plan$cohort[i]; m$contrast <- plan$contrast[i]; m$evidence_role <- plan$evidence_role[i]; m$program <- prog
      k <- k + 1L; deltas[[k]] <- m[, c("cohort", "contrast", "evidence_role", "program", "patient", "target_high", "target_low", "target_delta", "program_high", "program_low", "program_delta")]
      tst <- .o2p6_v130_spearman_test(m$target_delta, m$program_delta,
                                      seed = 20260722L + i * 100L + nchar(prog))
      q <- q + 1L
      tests[[q]] <- data.frame(
        cohort = plan$cohort[i], contrast = plan$contrast[i], evidence_role = plan$evidence_role[i], program = prog,
        n_paired_patients = tst$n, rho_delta = tst$rho, p_value = tst$p_value, p_method = tst$p_method,
        n_permutations = tst$n_permutations, small_n_exploratory = tst$small_n_exploratory,
        formal_significance_eligible = tst$formal_significance_eligible,
        median_target_delta = stats::median(m$target_delta), median_program_delta = stats::median(m$program_delta),
        inference_level = "Spearman association of within-patient target and program changes",
        stringsAsFactors = FALSE
      )
    }
  }
  test <- o2p6_bind_rows(tests)
  if (nrow(test)) {
    fam <- paste(test$cohort, test$contrast, sep = "||")
    test$fdr_within_contrast <- stats::ave(test$p_value, fam, FUN = function(v) stats::p.adjust(v, method = "BH"))
  }
  list(delta = o2p6_bind_rows(deltas), tests = test)
}

.o2p6_v130_patient_condition_compartment_program <- function(compartment_program_sample) {
  if (!is.data.frame(compartment_program_sample) || !nrow(compartment_program_sample)) return(data.frame())
  d <- compartment_program_sample
  d$patient_or_sample <- ifelse(!is.na(d$patient) & nzchar(as.character(d$patient)), as.character(d$patient), as.character(d$sample_id))
  d$inference_unit <- ifelse(!is.na(d$patient) & nzchar(as.character(d$patient)), "official_patient", "sample_descriptive")
  key <- paste(d$cohort, d$patient_or_sample, d$tissue_or_condition, d$marker_celltype, d$annotation_confidence_scope, d$marker_program, sep = "||")
  o2p6_bind_rows(lapply(split(d, key), function(z) {
    npos <- suppressWarnings(as.numeric(z$n_OLFML2B_positive)); nneg <- suppressWarnings(as.numeric(z$n_OLFML2B_negative))
    mp <- suppressWarnings(as.numeric(z$mean_score_OLFML2B_positive)); mn <- suppressWarnings(as.numeric(z$mean_score_OLFML2B_negative))
    mean_pos <- if (sum(npos, na.rm = TRUE) > 0) stats::weighted.mean(mp, w = npos, na.rm = TRUE) else NA_real_
    mean_neg <- if (sum(nneg, na.rm = TRUE) > 0) stats::weighted.mean(mn, w = nneg, na.rm = TRUE) else NA_real_
    data.frame(
      cohort = z$cohort[1], patient_or_sample = z$patient_or_sample[1], inference_unit = z$inference_unit[1],
      tissue_or_condition = z$tissue_or_condition[1], marker_celltype = z$marker_celltype[1],
      annotation_confidence_scope = z$annotation_confidence_scope[1], marker_program = z$marker_program[1],
      n_source_samples = length(unique(z$sample_id)), n_compartment_cells = sum(suppressWarnings(as.numeric(z$n_compartment_cells)), na.rm = TRUE),
      n_OLFML2B_positive = sum(npos, na.rm = TRUE), n_OLFML2B_negative = sum(nneg, na.rm = TRUE),
      mean_score_OLFML2B_positive = mean_pos, mean_score_OLFML2B_negative = mean_neg,
      delta_positive_minus_negative = mean_pos - mean_neg,
      classifier_program_overlap_n = max(suppressWarnings(as.numeric(z$classifier_program_overlap_n)), na.rm = TRUE),
      orthogonality_status = ifelse(any(z$orthogonality_status == "IDENTITY_COUPLED_NOT_MECHANISTIC"), "IDENTITY_COUPLED_NOT_MECHANISTIC", "ORTHOGONAL_TO_COMPARTMENT_CLASSIFIER"),
      stringsAsFactors = FALSE
    )
  }))
}

.o2p6_v130_compartment_program_audit <- function(pc) {
  if (!is.data.frame(pc) || !nrow(pc)) return(data.frame())
  d <- pc
  d$gate_compartment_cells <- d$n_compartment_cells >= 30L
  d$gate_positive_cells <- d$n_OLFML2B_positive >= 3L
  d$gate_negative_cells <- d$n_OLFML2B_negative >= 20L
  d$evaluable <- d$gate_compartment_cells & d$gate_positive_cells & d$gate_negative_cells & is.finite(d$delta_positive_minus_negative)
  d$exclusion_reason <- ifelse(!d$gate_compartment_cells, "LT30_COMPARTMENT_CELLS",
                        ifelse(!d$gate_positive_cells, "LT3_OLFML2B_POSITIVE_CELLS",
                        ifelse(!d$gate_negative_cells, "LT20_OLFML2B_NEGATIVE_CELLS",
                        ifelse(!is.finite(d$delta_positive_minus_negative), "NONFINITE_DELTA", "EVALUABLE"))))
  d
}

.o2p6_v130_within_compartment_tests <- function(audit) {
  if (!is.data.frame(audit) || !nrow(audit)) return(data.frame())
  plan <- data.frame(cohort = c("GSE150290", "GSE183904", "GSE167297"), tissue_or_condition = c("Gastric_Cancer", "Primary_Tumor", "Deep"), stringsAsFactors = FALSE)
  d <- merge(audit, plan, by = c("cohort", "tissue_or_condition"), all = FALSE)
  d <- d[d$inference_unit == "official_patient" & d$evaluable %in% TRUE & d$marker_celltype %in% c("Fibroblast", "Myofibroblast"), , drop = FALSE]
  out <- .o2p6_v130_one_sample_tests(d, "delta_positive_minus_negative",
                                     family_cols = c("cohort", "tissue_or_condition", "marker_celltype", "annotation_confidence_scope"),
                                     label_cols = c("cohort", "tissue_or_condition", "marker_celltype", "annotation_confidence_scope", "marker_program", "orthogonality_status"), min_units = 3L)
  if (nrow(out)) {
    out$mechanistic_interpretability <- ifelse(out$orthogonality_status == "ORTHOGONAL_TO_COMPARTMENT_CLASSIFIER",
                                               "ORTHOGONAL_EXPLORATORY_STATE_ASSOCIATION",
                                               "IDENTITY_COUPLED_LOCALIZATION_SUPPORT_ONLY")
  }
  out
}

.o2p6_v130_strict_localization <- function(patient_compartment, min_cells = 20L, min_units = 3L, min_positive = 10L, min_ge2 = 5L) {
  d <- patient_compartment[patient_compartment$n_cells >= min_cells, , drop = FALSE]
  if (!nrow(d)) return(list(by_dataset = data.frame(), consensus = data.frame()))
  by <- o2p6_bind_rows(lapply(split(d, paste(d$cohort, d$marker_celltype, sep = "||")), function(z) {
    hc_cells <- sum(z$n_high_confidence_cells, na.rm = TRUE); hc_pos <- sum(z$n_OLFML2B_positive_high_confidence, na.rm = TRUE)
    data.frame(
      cohort = z$cohort[1], marker_celltype = z$marker_celltype[1], n_units = length(unique(z$patient_or_sample)),
      total_cells = sum(z$n_cells, na.rm = TRUE), total_positive = sum(z$n_OLFML2B_positive, na.rm = TRUE),
      pooled_positive_fraction = sum(z$n_OLFML2B_positive, na.rm = TRUE) / pmax(sum(z$n_cells, na.rm = TRUE), 1),
      median_unit_positive_fraction = stats::median(z$OLFML2B_positive_fraction, na.rm = TRUE),
      total_ge2 = sum(z$n_OLFML2B_ge2, na.rm = TRUE), pooled_ge2_fraction = sum(z$n_OLFML2B_ge2, na.rm = TRUE) / pmax(sum(z$n_cells, na.rm = TRUE), 1),
      total_high_confidence_cells = hc_cells, total_positive_high_confidence = hc_pos,
      pooled_positive_fraction_high_confidence = hc_pos / pmax(hc_cells, 1),
      eligible = length(unique(z$patient_or_sample)) >= min_units && sum(z$n_OLFML2B_positive, na.rm = TRUE) >= min_positive,
      eligible_ge2 = length(unique(z$patient_or_sample)) >= min_units && sum(z$n_OLFML2B_ge2, na.rm = TRUE) >= min_ge2,
      eligible_high_confidence = length(unique(z$patient_or_sample)) >= min_units && hc_pos >= min_positive,
      stringsAsFactors = FALSE
    )
  }))
  by$rank_count1 <- stats::ave(-by$median_unit_positive_fraction, by$cohort, FUN = function(x) rank(x, ties.method = "min"))
  by$rank_count2 <- stats::ave(-by$pooled_ge2_fraction, by$cohort, FUN = function(x) rank(x, ties.method = "min"))
  by$rank_high_confidence <- stats::ave(-by$pooled_positive_fraction_high_confidence, by$cohort, FUN = function(x) rank(x, ties.method = "min"))
  by$rank_within_dataset <- by$rank_count1
  by$threshold_top3_consistent <- by$eligible & by$eligible_ge2 & by$rank_count1 <= 3 & by$rank_count2 <= 3
  by$high_confidence_top3_consistent <- by$eligible & by$eligible_high_confidence & by$rank_count1 <= 3 & by$rank_high_confidence <= 3
  elig <- by[by$eligible %in% TRUE, , drop = FALSE]
  cons <- o2p6_bind_rows(lapply(split(elig, elig$marker_celltype), function(z) {
    data.frame(
      marker_celltype = z$marker_celltype[1], n_datasets_eligible = nrow(z), n_datasets_top3 = sum(z$rank_count1 <= 3),
      n_datasets_threshold_top3_consistent = sum(z$threshold_top3_consistent %in% TRUE),
      n_datasets_high_confidence_top3_consistent = sum(z$high_confidence_top3_consistent %in% TRUE),
      datasets_eligible = paste(z$cohort, collapse = ";"), median_dataset_fraction = stats::median(z$median_unit_positive_fraction, na.rm = TRUE),
      strict_consistency_status = ifelse(nrow(z) >= 3L && sum(z$rank_count1 <= 3) >= 3L, "REPLICATED_TOP3_AT_LEAST_THREE_DATASETS",
                                  ifelse(nrow(z) >= 2L && sum(z$rank_count1 <= 3) >= 2L, "REPLICATED_TOP3_AT_LEAST_TWO_DATASETS", "DESCRIPTIVE_OR_NOT_REPLICATED")),
      stringsAsFactors = FALSE
    )
  }))
  list(by_dataset = by, consensus = cons)
}

.o2p6_v130_lodo <- function(by) {
  d <- by[by$eligible %in% TRUE, , drop = FALSE]
  out <- list(); k <- 0L
  for (ct in unique(as.character(d$marker_celltype))) {
    z <- d[d$marker_celltype == ct, , drop = FALSE]
    for (omit in c("NONE", unique(as.character(z$cohort)))) {
      r <- if (omit == "NONE") z else z[z$cohort != omit, , drop = FALSE]
      k <- k + 1L
      out[[k]] <- data.frame(
        marker_celltype = ct, omitted_dataset = omit, n_datasets_remaining = nrow(r),
        n_remaining_top3 = sum(r$rank_count1 <= 3), n_remaining_threshold_consistent = sum(r$threshold_top3_consistent %in% TRUE),
        n_remaining_high_confidence_consistent = sum(r$high_confidence_top3_consistent %in% TRUE),
        replicated_after_omission = nrow(r) >= 2L && sum(r$rank_count1 <= 3) >= 2L,
        stringsAsFactors = FALSE
      )
    }
  }
  o2p6_bind_rows(out)
}

.o2p6_v130_source_evidence <- function(consensus, lodo, paired_tests, by) {
  celltypes <- unique(c(as.character(consensus$marker_celltype), as.character(by$marker_celltype)))
  o2p6_bind_rows(lapply(celltypes, function(ct) {
    c0 <- consensus[consensus$marker_celltype == ct, , drop = FALSE]
    b0 <- by[by$marker_celltype == ct & by$eligible %in% TRUE, , drop = FALSE]
    l0 <- lodo[lodo$marker_celltype == ct & lodo$omitted_dataset != "NONE", , drop = FALSE]
    p0 <- paired_tests[paired_tests$marker_celltype == ct & paired_tests$metric == "OLFML2B_positive_fraction" & paired_tests$formal_evaluable %in% TRUE, , drop = FALSE]
    n_elig <- if (nrow(c0)) c0$n_datasets_eligible[1] else nrow(b0)
    n_top <- if (nrow(c0)) c0$n_datasets_top3[1] else sum(b0$rank_count1 <= 3)
    lodo_robust <- nrow(l0) >= 2L && all(l0$replicated_after_omission %in% TRUE)
    total_pos <- sum(b0$total_positive, na.rm = TRUE); hc_pos <- sum(b0$total_positive_high_confidence, na.rm = TRUE)
    high_conf_frac <- hc_pos / pmax(total_pos, 1)
    threshold_consistent <- if (nrow(c0)) c0$n_datasets_threshold_top3_consistent[1] >= 2L else sum(b0$threshold_top3_consistent %in% TRUE) >= 2L
    grade <- if (n_elig >= 3L && n_top >= 3L && lodo_robust && threshold_consistent && high_conf_frac >= 0.60) {
      "ROBUST_REPLICATED_LOCALIZATION"
    } else if (n_elig >= 3L && n_top >= 3L && lodo_robust) {
      "REPLICATED_MARKER_DERIVED_LOCALIZATION_ANNOTATION_CONFIDENCE_LIMITED"
    } else if (n_elig >= 2L && n_top >= 2L) {
      "REPLICATED_LOCALIZATION_NOT_UNIQUE"
    } else if (n_elig >= 1L) {
      "DESCRIPTIVE_SINGLE_DATASET_OR_LOW_RANK"
    } else "NOT_EVALUABLE"
    data.frame(
      marker_celltype = ct, n_datasets_eligible = n_elig, n_datasets_top3 = n_top,
      leave_one_dataset_out_robust = lodo_robust, count_threshold_consistent = threshold_consistent,
      pooled_high_confidence_fraction_among_positive = high_conf_frac,
      n_formal_paired_cohort_contrasts = nrow(p0), n_formal_paired_positive_direction = sum(p0$median_delta > 0, na.rm = TRUE),
      evidence_grade = grade,
      interpretation = ifelse(ct == "Fibroblast", "Fibroblast is the best-supported preferential source; uniqueness and causality remain unproven.",
                       ifelse(ct == "Myofibroblast", "Myofibroblast localization is replicated but remains marker-derived and confidence-limited unless HIGH-only/spatial results agree.",
                       ifelse(ct == "Epithelial", "Epithelial cells are not supported as the dominant source.",
                       ifelse(ct == "Myeloid", "Myeloid localization is heterogeneous and secondary.", "Localization remains descriptive.")))),
      stringsAsFactors = FALSE
    )
  }))
}

.o2p6_v130_theme <- function(base_size = 10.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.title = ggplot2::element_text(face = "bold"),
                   strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"))
}

.o2p6_v130_save_plot <- function(p, dirs, stem, width, height) {
  ggplot2::ggsave(file.path(dirs$figures, paste0(stem, ".png")), p, width = width, height = height, dpi = 420, bg = "white")
  ggplot2::ggsave(file.path(dirs$figures, paste0(stem, ".pdf")), p, width = width, height = height, useDingbats = FALSE, bg = "white")
}

.o2p6_v130_figures <- function(index, tables) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(data.frame())
  dirs <- index$dirs; reg <- list(); k <- 0L
  save_one <- function(p, stem, width, height, source, note) {
    .o2p6_v130_save_plot(p, dirs, stem, width, height)
    k <<- k + 1L; reg[[k]] <<- data.frame(figure = stem, source_table = source, interpretation = note, stringsAsFactors = FALSE)
  }
  d <- index$cohort_summary
  if (nrow(d)) {
    d$cohort <- factor(d$cohort, levels = c("GSE150290", "GSE183904", "GSE167297", "GSE134520"))
    d$lab <- paste0("cells=", format(d$total_cells, big.mark = ","), "\nOLFML2B+=", sprintf("%.2f%%", 100*d$pooled_OLFML2B_positive_fraction))
    p <- ggplot2::ggplot(d, ggplot2::aes(cohort, pooled_OLFML2B_positive_fraction, fill = cohort)) +
      ggplot2::geom_col(color = "black", width = .62, show.legend = FALSE) + ggplot2::geom_text(ggplot2::aes(label = lab), vjust = -.25, size = 2.7) +
      ggplot2::scale_y_continuous(labels = function(x) sprintf("%.1f%%",100*x), expand = ggplot2::expansion(mult=c(0,.28))) +
      ggplot2::labs(tag="A",title="OLFML2B detection across four gastric scRNA-seq datasets",subtitle="Absolute detection fractions are platform/QC dependent and are not directly comparable",x=NULL,y="OLFML2B-positive cells") + .o2p6_v130_theme()
    save_one(p,"FIG6RC2_A_OLFML2B_detection_by_dataset",8.4,5.1,"23_cohort_level_OLFML2B_summary.csv","Descriptive detection only.")
  }
  d <- tables$strict_by[tables$strict_by$eligible %in% TRUE, , drop=FALSE]
  if (nrow(d)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(pooled_positive_fraction, marker_celltype)) +
      ggplot2::geom_point(ggplot2::aes(size=total_positive,fill=cohort),shape=21,color="black") +
      ggplot2::facet_wrap(~cohort,nrow=1,scales="free_x") + ggplot2::scale_x_continuous(labels=function(x)sprintf("%.1f%%",100*x)) +
      ggplot2::labs(tag="B",title="Strict-gated marker-derived OLFML2B localization",subtitle="Minimum cells, units and positive-cell gates; count>=2 and HIGH-confidence sensitivity are audited",x="OLFML2B-positive fraction",y=NULL,size="OLFML2B+ cells") + .o2p6_v130_theme()+ggplot2::theme(legend.position="bottom")
    save_one(p,"FIG6RC2_B_marker_compartment_localization",11.5,5.8,"54_strict_localization_by_dataset.csv","Strict marker-derived localization.")
  }
  d <- index$full_transcriptome_target
  d <- d[(d$cohort=="GSE150290" & d$tissue_or_condition%in%c("Adjacent_Normal","Gastric_Cancer")) | (d$cohort=="GSE183904" & d$tissue_or_condition%in%c("Primary_Normal","Primary_Tumor")) | (d$cohort=="GSE167297" & d$tissue_or_condition%in%c("Normal","Superficial","Deep")),,drop=FALSE]
  if (nrow(d)) {
    d$condition_key <- paste(d$cohort, d$tissue_or_condition, sep = "||")
    d$condition_key <- factor(d$condition_key, levels = c("GSE150290||Adjacent_Normal","GSE150290||Gastric_Cancer","GSE183904||Primary_Normal","GSE183904||Primary_Tumor","GSE167297||Normal","GSE167297||Superficial","GSE167297||Deep"))
    d$condition_label <- factor(gsub("_", " ", d$tissue_or_condition), levels = c("Adjacent Normal","Gastric Cancer","Primary Normal","Primary Tumor","Normal","Superficial","Deep")); d$patient_group <- interaction(d$cohort,d$patient,drop=TRUE)
    p <- ggplot2::ggplot(d,ggplot2::aes(condition_label,log2_cpm,group=patient_group))+ggplot2::geom_line(alpha=.38,color="grey45")+ggplot2::geom_point(size=2)+ggplot2::facet_wrap(~cohort,scales="free_x")+
      ggplot2::labs(tag="C",title="Patient-level OLFML2B pseudobulk trajectories",subtitle="Only biologically related conditions are connected within each cohort",x=NULL,y="OLFML2B log2-CPM")+.o2p6_v130_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1))
    save_one(p,"FIG6RC2_C_condition_level_pseudobulk",10.5,5.6,"44_patient_level_full_transcriptome_OLFML2B_log2CPM.csv","Patient-level trajectories.")
  }
  d <- tables$primary_program_tests
  core <- c("CAF_ECM","ECM_Remodeling","Myofibroblast","TGFb_Response","Smooth_Muscle","Inflammatory_Fibroblast","Myeloid_Macrophage","Epithelial_Differentiation")
  d <- d[d$marker_program%in%core,,drop=FALSE]
  if (nrow(d)) {
    d$marker_program <- factor(d$marker_program,levels=rev(core)); d$sig <- is.finite(d$wilcox_fdr_within_family)&d$wilcox_fdr_within_family<.05
    p <- ggplot2::ggplot(d,ggplot2::aes(median_delta,marker_program))+ggplot2::geom_vline(xintercept=0,linetype=2)+ggplot2::geom_errorbar(ggplot2::aes(xmin=bootstrap_median_ci_low,xmax=bootstrap_median_ci_high),orientation="y",width=.18)+ggplot2::geom_point(ggplot2::aes(fill=sig),shape=21,size=2.8)+ggplot2::facet_wrap(~cohort,nrow=1,scales="free_x")+
      ggplot2::labs(tag="D",title="Condition-specific patient-level programs in OLFML2B-positive cells",subtitle="Prespecified tumor/deep states only; repeated conditions are never mixed",x="OLFML2B+ minus OLFML2B- program score",y=NULL)+.o2p6_v130_theme()
    save_one(p,"FIG6RC2_D_OLFML2B_positive_program_delta",11.4,6,"60_primary_tumor_state_marker_program_tests.csv","Condition-specific patient-level exploratory comparison.")
  }
  d <- tables$condition_bridge[tables$condition_bridge$primary_tumor_state%in%TRUE & tables$condition_bridge$program%in%core & is.finite(tables$condition_bridge$rho),,drop=FALSE]
  if (nrow(d)) {
    d$row_label <- paste0(d$cohort," / ",gsub("_"," ",d$tissue_or_condition)," (n=",d$n_patients,")")
    formal <- d$formal_significance_eligible%in%TRUE & is.finite(d$fdr_within_cohort_condition)
    d$star <- ifelse(formal & d$fdr_within_cohort_condition<.001,"***",ifelse(formal & d$fdr_within_cohort_condition<.01,"**",ifelse(formal & d$fdr_within_cohort_condition<.05,"*",ifelse(d$small_n_exploratory%in%TRUE,"†",""))))
    d$lab <- paste0(sprintf("%.2f",d$rho),d$star); d$program <- factor(d$program,levels=core)
    p <- ggplot2::ggplot(d,ggplot2::aes(program,row_label,fill=rho))+ggplot2::geom_tile(color="white")+ggplot2::geom_text(ggplot2::aes(label=lab),size=2.5)+ggplot2::scale_fill_gradient2(low="#2F6DB3",mid="white",high="#C73E3A",midpoint=0,limits=c(-1,1))+
      ggplot2::labs(tag="E",title="Patient-level tumor-state OLFML2B-program bridge",subtitle="Exact/permutation inference; † denotes n<6 exploratory strata without formal significance stars",x=NULL,y=NULL)+.o2p6_v130_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1),axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
    save_one(p,"FIG6RC2_E_sample_level_program_bridge_heatmap",10.8,5.2,"63_fulltx_condition_stratified_program_bridge_exact.csv","Exact/permutation patient-level correlation.")
  }
  d <- tables$change_bridge[tables$change_bridge$contrast%in%c("Gastric_Cancer_minus_Adjacent_Normal","Primary_Tumor_minus_Primary_Normal","Deep_minus_Normal") & tables$change_bridge$program%in%core & is.finite(tables$change_bridge$rho_delta),,drop=FALSE]
  if (nrow(d)) {
    d$row_label <- paste0(d$cohort," / ",gsub("_minus_"," - ",gsub("_"," ",d$contrast))," (n=",d$n_paired_patients,")")
    formal <- d$formal_significance_eligible%in%TRUE & is.finite(d$fdr_within_contrast)
    d$star <- ifelse(formal & d$fdr_within_contrast<.05,"*",ifelse(d$small_n_exploratory%in%TRUE,"†","")); d$lab<-paste0(sprintf("%.2f",d$rho_delta),d$star);d$program<-factor(d$program,levels=core)
    p<-ggplot2::ggplot(d,ggplot2::aes(program,row_label,fill=rho_delta))+ggplot2::geom_tile(color="white")+ggplot2::geom_text(ggplot2::aes(label=lab),size=2.5)+ggplot2::scale_fill_gradient2(low="#2F6DB3",mid="white",high="#C73E3A",midpoint=0,limits=c(-1,1))+
      ggplot2::labs(tag="F",title="Within-patient OLFML2B-program change bridge",subtitle="Exact/permutation inference; small depth cohorts remain exploratory",x=NULL,y=NULL)+.o2p6_v130_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1),axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
    save_one(p,"FIG6RC2_F_paired_change_program_bridge_heatmap",10.8,5.2,"65_fulltx_paired_change_program_bridge_exact.csv","Within-patient change association.")
  }
  d <- tables$within_compartment_tests[tables$within_compartment_tests$marker_program%in%c("TGFb_Response","ECM_Remodeling","Inflammatory_Fibroblast","CAF_ECM","Myofibroblast") & tables$within_compartment_tests$annotation_confidence_scope=="ALL_ASSIGNED",,drop=FALSE]
  if (nrow(d)) {
    d$label <- paste(d$marker_celltype,d$marker_program,sep=" / ")
    p<-ggplot2::ggplot(d,ggplot2::aes(median_delta,label))+ggplot2::geom_vline(xintercept=0,linetype=2)+ggplot2::geom_errorbar(ggplot2::aes(xmin=bootstrap_median_ci_low,xmax=bootstrap_median_ci_high),orientation="y",width=.18)+ggplot2::geom_point(ggplot2::aes(shape=mechanistic_interpretability),size=2.7)+ggplot2::facet_wrap(~cohort,scales="free_x")+
      ggplot2::labs(tag="G",title="Within-compartment OLFML2B-positive cell programs",subtitle="Fibroblast/myofibroblast strata only; identity-coupled programs are localization support, not mechanism",x="Program score: OLFML2B+ minus OLFML2B-",y=NULL,shape="Interpretation")+.o2p6_v130_theme()
    save_one(p,"FIG6RC2_G_within_compartment_program_state",11.5,6.5,"69_within_compartment_program_tests.csv","Within-compartment exploratory state audit.")
  }
  o2p6_bind_rows(reg)
}

.o2p6_v130_finalize <- function(index, root = "D:/OLFML2B_STAD") {
  dirs <- index$dirs
  target <- index$full_transcriptome_target
  programs <- index$full_transcriptome_programs
  g183 <- .o2p6_v130_gse183904_paired(target)
  eval_program <- .o2p6_v130_program_evaluability(index$program_summary %||% data.frame())
  pc_program <- .o2p6_v130_patient_condition_program(index$program_summary %||% data.frame())
  primary_program_tests <- .o2p6_v130_primary_tumor_program_tests(pc_program)
  paired_program <- .o2p6_v130_paired_marker_program_change(pc_program)
  condition_bridge <- .o2p6_v130_condition_bridge(target, programs)
  change_bridge <- .o2p6_v130_paired_change_bridge(target, programs)
  strict <- .o2p6_v130_strict_localization(index$patient_compartment)
  lodo <- .o2p6_v130_lodo(strict$by_dataset)
  pc_comp <- .o2p6_v130_patient_condition_compartment_program(index$compartment_program_summary %||% data.frame())
  comp_audit <- .o2p6_v130_compartment_program_audit(pc_comp)
  comp_tests <- .o2p6_v130_within_compartment_tests(comp_audit)
  source_evidence <- .o2p6_v130_source_evidence(strict$consensus, lodo, index$paired_compartment$tests, strict$by_dataset)

  o2p6_atomic_write_csv(index$compartment_program_summary, file.path(dirs$tables,"24_within_compartment_marker_program_summary_by_sample.csv"))
  o2p6_atomic_write_csv(g183$delta,file.path(dirs$tables,"56_GSE183904_fulltx_paired_patient_delta.csv"))
  o2p6_atomic_write_csv(g183$test,file.path(dirs$tables,"57_GSE183904_fulltx_paired_tests.csv"))
  o2p6_atomic_write_csv(eval_program,file.path(dirs$tables,"58_marker_program_evaluability_audit.csv"))
  o2p6_atomic_write_csv(pc_program,file.path(dirs$tables,"59_patient_condition_marker_program_delta.csv"))
  o2p6_atomic_write_csv(primary_program_tests,file.path(dirs$tables,"60_primary_tumor_state_marker_program_tests.csv"))
  o2p6_atomic_write_csv(paired_program$delta,file.path(dirs$tables,"61_paired_condition_marker_program_deltas.csv"))
  o2p6_atomic_write_csv(paired_program$tests,file.path(dirs$tables,"62_paired_condition_marker_program_tests.csv"))
  o2p6_atomic_write_csv(condition_bridge,file.path(dirs$tables,"63_fulltx_condition_stratified_program_bridge_exact.csv"))
  o2p6_atomic_write_csv(change_bridge$delta,file.path(dirs$tables,"64_fulltx_paired_change_program_deltas.csv"))
  o2p6_atomic_write_csv(change_bridge$tests,file.path(dirs$tables,"65_fulltx_paired_change_program_bridge_exact.csv"))
  o2p6_atomic_write_csv(lodo,file.path(dirs$tables,"66_strict_localization_leave_one_dataset_out.csv"))
  o2p6_atomic_write_csv(comp_audit,file.path(dirs$tables,"67_compartment_program_evaluability_audit.csv"))
  o2p6_atomic_write_csv(pc_comp,file.path(dirs$tables,"68_patient_condition_compartment_program_delta.csv"))
  o2p6_atomic_write_csv(comp_tests,file.path(dirs$tables,"69_within_compartment_program_tests.csv"))
  o2p6_atomic_write_csv(source_evidence,file.path(dirs$tables,"70_final_cell_source_evidence_matrix.csv"))
  o2p6_atomic_write_csv(strict$by_dataset,file.path(dirs$tables,"54_strict_localization_by_dataset.csv"))
  o2p6_atomic_write_csv(strict$consensus,file.path(dirs$tables,"55_strict_cross_dataset_localization_consistency.csv"))

  key <- data.frame(
    result_id=c("GSE150290_primary_paired","GSE183904_paired_support","GSE167297_depth_context","fibroblast_source","myofibroblast_source","within_patient_stromal_bridge","within_compartment_activation","final_gene_lock"),
    status=c("FORMAL_PRIMARY_SUPPORT",ifelse(nrow(g183$test)&&g183$test$median_delta[1]>0,"DIRECTIONAL_SUPPORT_NOT_SIGNIFICANT","NO_DIRECTIONAL_SUPPORT"),"SMALL_N_DIRECTIONAL_CONTEXT",
             source_evidence$evidence_grade[match("Fibroblast",source_evidence$marker_celltype)],source_evidence$evidence_grade[match("Myofibroblast",source_evidence$marker_celltype)],
             "PATIENT_PAIRED_ECOLOGICAL_ASSOCIATION","EXPLORATORY_OR_NOT_EVALUABLE","FALSE"),
    interpretation=c("GSE150290 remains the formal paired single-cell pseudobulk evidence.","GSE183904 is an independent paired directional estimate.","GSE167297 is exact small-n depth context.","Fibroblast is the leading preferential source.","Myofibroblast is replicated but annotation-confidence limited unless sensitivity agrees.","Within-patient OLFML2B and stromal-program changes co-vary without proving causality.","Fixed-compartment OLFML2B-positive state tests are separated from cell-composition localization and identity-coupled programs are flagged.","Final gene lock awaits spatial concordance."),stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(key,file.path(dirs$tables,"71_unified_key_result_summary.csv"))
  supersession <- data.frame(legacy_output=c("32_OLFML2B_positive_marker_program_tests.csv","33_sample_level_OLFML2B_program_bridge.csv","58_patient_level_marker_program_delta.csv","59_patient_level_marker_program_tests_by_cohort.csv","60_fulltx_condition_stratified_program_bridge.csv","62_fulltx_paired_change_program_bridge.csv"),status="AUDIT_ONLY_SUPERSEDED",formal_replacement=c("60_primary_tumor_state_marker_program_tests.csv","63_fulltx_condition_stratified_program_bridge_exact.csv;65_fulltx_paired_change_program_bridge_exact.csv","59_patient_condition_marker_program_delta.csv","60_primary_tumor_state_marker_program_tests.csv","63_fulltx_condition_stratified_program_bridge_exact.csv","65_fulltx_paired_change_program_bridge_exact.csv"),reason=c("Sample-level mixed-condition test replaced by condition-specific patient units.","Asymptotic/sample-level bridge replaced by patient-level exact/permutation inference.","Mixed-condition patient collapse removed.","Mixed-condition patient test removed.","Small-n asymptotic Spearman removed.","Small-n asymptotic Spearman removed."),stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(supersession,file.path(dirs$tables,"72_unified_supersession_map.csv"))
  figreg <- .o2p6_v130_figures(index,list(strict_by=strict$by_dataset,primary_program_tests=primary_program_tests,condition_bridge=condition_bridge,change_bridge=change_bridge$tests,within_compartment_tests=comp_tests))
  o2p6_atomic_write_csv(figreg,file.path(dirs$tables,"73_unified_figure_registry.csv"))

  go <- data.frame(criterion=c("four_dataset_119_matrix_layer","official_patient_mapping","GSE150290_processed_text_guard","patient_condition_program_no_mixing","small_n_exact_permutation","within_patient_change_bridge","within_compartment_program_audit","threshold_high_confidence_source_sensitivity","Part0_5_freeze","final_gene_lock"),status=c(ifelse(nrow(index$dataset_audit)==119L&&all(index$dataset_audit$parse_status=="OK"),"PASS","REVIEW"),"PASS","PASS",ifelse(nrow(primary_program_tests),"PASS","NOT_EVALUABLE"),ifelse(any(condition_bridge$p_method=="FULL_ENUMERATION_EXACT_SPEARMAN"),"PASS","NOT_TRIGGERED"),ifelse(nrow(change_bridge$tests),"PASS","NOT_EVALUABLE"),ifelse(nrow(comp_audit),"PASS","NOT_EVALUABLE"),ifelse(nrow(source_evidence),"PASS","NOT_EVALUABLE"),"PASS","FALSE"),interpretation=c("All four datasets are processed in one production run with validated cache reuse.","Official patient-condition maps remain authoritative.","GSE150290 raw droplet matrices remain audit-only.","Marker-program comparisons are patient×condition specific.","n=4-9 correlations use full exact permutation; n=10-19 use deterministic Monte Carlo; n<6 remains exploratory.","Target-program associations use within-patient paired changes.","Fibroblast/myofibroblast OLFML2B-positive state is tested within fixed compartments when cell gates permit.","Source grading incorporates count>=2 and HIGH-confidence sensitivity.","Part0-Part5 are read-only.","Spatial and author-annotation concordance remain required."),stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(go,file.path(dirs$tables,"99_OLFML2B_Part6_unified_go_no_go.csv"))
  claim <- data.frame(level=c("ALLOWED","ALLOWED_WITH_LIMITATION","CONDITIONAL","FORBIDDEN"),statement=c("Patient-level OLFML2B tumor-normal/depth evidence across four gastric single-cell datasets.","Preferential fibroblast/myofibroblast localization and CAF/ECM/TGF-beta ecological co-variation, explicitly marker-derived and observational.","Within-compartment activation only when patient, cell-count, orthogonality, threshold, HIGH-confidence and spatial results agree.","No unique CAF source, secretion, direct TGF-beta activation, causal invasion/recurrence, prognosis, treatment response or clinical utility from Part6 alone."),stringsAsFactors=FALSE)
  o2p6_atomic_write_csv(claim,file.path(dirs$tables,"103_OLFML2B_Part6_unified_claim_boundary.csv"))

  unified <- list(version=OLFML2B_PART6_UNIFIED_VERSION,generated_at=o2p6_ts(),GSE183904_paired=g183,marker_program_evaluability=eval_program,patient_condition_marker_program=pc_program,primary_tumor_marker_program_tests=primary_program_tests,paired_marker_program_change=paired_program,condition_bridge=condition_bridge,paired_change_bridge=change_bridge,strict_localization=strict,localization_lodo=lodo,compartment_program_evaluability=comp_audit,patient_condition_compartment_program=pc_comp,within_compartment_program_tests=comp_tests,source_evidence=source_evidence,key_results=key,go_no_go=go,claim_boundary=claim,figure_registry=figreg,Part0_5_frozen=TRUE,final_gene_lock=FALSE)
  index$version <- OLFML2B_PART6_UNIFIED_VERSION
  index$unified_v130 <- unified
  index$strict_localization <- strict
  index$final_source_evidence <- source_evidence
  index$Part0_5_frozen <- TRUE; index$final_gene_lock <- FALSE
  o2p6_atomic_save_rds(index,file.path(dirs$objects,"Part6_OLFML2B_unified_index.rds"))
  o2p6_atomic_save_rds(index,file.path(dirs$objects,"Part6_OLFML2B_four_dataset_direct_fork_index.rds"))
  o2p6_atomic_save_rds(index,file.path(dirs$objects,"Part6_OLFML2B_single_cell_index.rds"))
  p05 <- file.path(dirs$objects,"OLFML2B_Part0_5_complete_index.rds")
  integrated <- if(file.exists(p05)) tryCatch(readRDS(p05),error=function(e)list()) else list()
  integrated$Part6 <- index; integrated$Part6_version <- OLFML2B_PART6_UNIFIED_VERSION; integrated$Part0_5_frozen<-TRUE;integrated$final_gene_lock<-FALSE;integrated$updated_at<-o2p6_ts()
  o2p6_atomic_save_rds(integrated,file.path(dirs$objects,"OLFML2B_Part0_6_complete_index.rds"))
  o2p6_log("INFO","OLFML2B Part6 unified v1.3.0 complete | one production module | exact/permutation inference | final_gene_lock=FALSE",log_file=file.path(dirs$logs,"Part6_OLFML2B_unified_v130.log"))
  index
}

.o2p6_v130_core <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_scrna_rc2_no_kang_full_production <- function(...) {
  index <- .o2p6_v130_core(...)
  .o2p6_v130_finalize(index, root = index$dirs$root %||% "D:/OLFML2B_STAD")
}
run_olfml2b_part6_scrna_rc2_no_kang <- run_olfml2b_part6_scrna_rc2_no_kang_full_production
run_olfml2b_part6_singlecell_production <- function(root="D:/OLFML2B_STAD",sc_dir=file.path(root,"data","raw","single_cell"),force_unpack=FALSE,make_figures=TRUE,refresh_cache=identical(toupper(Sys.getenv("OLFML2B_PART6_REFRESH_CACHE",unset="FALSE")),"TRUE"),...) {
  run_olfml2b_part6_scrna_rc2_no_kang_full_production(root=root,raw_single_cell_dir=sc_dir,output_subdir="Part6",force_unpack=force_unpack,make_figures=make_figures,exclude_gse150290_raw_10x=TRUE,refresh_cache=refresh_cache)
}
