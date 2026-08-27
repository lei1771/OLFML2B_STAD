# ============================================================================
# OLFML2B-STAD Part7: spatial transcriptomic competing-source and stromal-niche evidence layer
# Version: v1.0.5_20260722_SELF_EXCLUDED_NORMALIZED_DISTANCE_AND_SPATIAL_PALETTE_REPAIR
# ============================================================================
# Purpose:
#   Rebuild Part7 around the loader that was proven to read all 10 Visium
#   samples in the independent diagnostic:
#       OK_SCORE_AND_SPATIAL_POSSIBLE = 10/10
#
#   This file deliberately avoids Part10 and directly generates native Part7
#   score tables, spatial maps, spatial co-localization summaries, neighborhood
#   enrichment, niche classification and CNS-style figures.
#
# Run:
#   setwd("D:/OLFML2B_STAD")
#   source("R/07_OLFML2B_PART7_SPATIAL_TRANSCRIPTOMICS.R")
#   p7 <- run_olfml2b_part7_spatial_transcriptomics(
#       root = "D:/OLFML2B_STAD",
#       spatial_dir = "D:/OLFML2B_STAD/data/raw/spatial/GSE251950",
#       output_subdir = "Part7",
#       allow_download = FALSE,
#       force_download = FALSE,
#       force_unpack = FALSE,
#       make_figures = TRUE
#   )
# ============================================================================

options(stringsAsFactors = FALSE)
OLFML2B_PART7_VERSION <- "v1.0.5_20260722_SELF_EXCLUDED_NORMALIZED_DISTANCE_AND_SPATIAL_PALETTE_REPAIR"

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

o2b_p7_official_gse251950_sample_map <- function() {
  data.frame(
    sample_id = c(
      "20_00331_LI_SING", "21_00731_LI_SING", "21_00732_LI_SING", "21_00733_LI_SING", "21_00734_LI_SING",
      "21_01251_LI_SING", "21_01252_LI_SING", "21_01253_LI_SING", "21_01254_LI_SING", "21_01675_LI_SING"
    ),
    geo_accession = paste0("GSM", 7990473:7990482),
    geo_title = c("GC1", "GC2", "GC3", "GC4", "GC5", "GC6", "GC6-PM", "GC7", "GC8", "GC9"),
    patient_id = c("GC1", "GC2", "GC3", "GC4", "GC5", "GC6", "GC6", "GC7", "GC8", "GC9"),
    tissue_role = c(rep("primary_gastric_cancer", 6L), "paired_metastasis", rep("primary_gastric_cancer", 3L)),
    mapping_source = "GSE251950 official GEO sample titles/descriptions",
    stringsAsFactors = FALSE
  )
}

o2b_p7_normalize_sample_id <- function(x) {
  x <- basename(trimws(as.character(x)))
  x <- sub("^GSM[0-9]+_", "", x, ignore.case = TRUE)
  x
}

o2b_p7_extract_geo_accession <- function(x) {
  x <- basename(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  hit <- grepl("^GSM[0-9]+", x, ignore.case = TRUE)
  out[hit] <- toupper(sub("^(GSM[0-9]+).*$", "\\1", x[hit], ignore.case = TRUE))
  out
}

o2b_p7_attach_official_sample_metadata <- function(d) {
  if (!is.data.frame(d) || !nrow(d) || !"sample_id" %in% names(d)) return(d)

  raw_id <- if ("sample_id_raw" %in% names(d)) as.character(d$sample_id_raw) else as.character(d$sample_id)
  clean_id <- o2b_p7_normalize_sample_id(d$sample_id)
  geo_from_raw <- o2b_p7_extract_geo_accession(raw_id)
  geo_existing <- if ("geo_accession" %in% names(d)) toupper(trimws(as.character(d$geo_accession))) else rep(NA_character_, nrow(d))
  geo_existing[is.na(geo_existing) | !nzchar(geo_existing)] <- NA_character_
  geo_candidate <- ifelse(!is.na(geo_existing), geo_existing, geo_from_raw)

  map <- o2b_p7_official_gse251950_sample_map()
  idx_sample <- match(clean_id, map$sample_id)
  idx_geo <- match(geo_candidate, map$geo_accession)
  idx <- idx_sample
  idx[is.na(idx)] <- idx_geo[is.na(idx)]

  d$sample_id_raw <- raw_id
  d$sample_id <- clean_id
  for (nm in setdiff(names(map), "sample_id")) d[[nm]] <- map[[nm]][idx]
  d$patient_mapping_status <- ifelse(!is.na(d$patient_id) & nzchar(d$patient_id), "OFFICIAL_GEO_MAP", "UNMAPPED")
  d
}

.o2b_part7_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.o2b_part7_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.o2b_part7_env_valid <- nzchar(.o2b_part7_env_root) && dir.exists(.o2b_part7_env_root)
.o2b_part7_code_root <- if (.o2b_part7_env_valid) {
  normalizePath(.o2b_part7_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.o2b_part7_entry) && file.exists(.o2b_part7_entry)) {
  dirname(normalizePath(.o2b_part7_entry, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .o2b_part7_code_root)

o2b_p7_ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

o2b_p7_log <- function(level = "INFO", ..., log_file = NULL, echo = TRUE) {
  line <- sprintf("[%s] [%s] [OLFML2B-P7] %s", o2b_p7_ts(), toupper(level), paste0(..., collapse = ""))
  if (isTRUE(echo)) message(line)
  if (!is.null(log_file) && nzchar(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
  invisible(line)
}

o2b_p7_stop <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

o2b_p7_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) o2b_p7_stop("Cannot create directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

o2b_p7_dirs <- function(root, output_subdir = "Part7") {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  part <- as.character(output_subdir)[1]
  list(
    root = root,
    tables = o2b_p7_dir(file.path(root, "output", "tables", part)),
    figures = o2b_p7_dir(file.path(root, "output", "figures", part)),
    reports = o2b_p7_dir(file.path(root, "output", "reports", part)),
    qc = o2b_p7_dir(file.path(root, "output", "qc", part)),
    logs = o2b_p7_dir(file.path(root, "logs", "runtime", part)),
    objects = o2b_p7_dir(file.path(root, "output", "objects"))
  )
}

o2b_p7_bind_rows <- function(xs) {
  if (is.null(xs)) return(data.frame())
  if (is.data.frame(xs)) return(xs)
  xs <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  aligned <- lapply(xs, function(z) {
    for (m in setdiff(cols, names(z))) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

o2b_p7_write_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  if (is.atomic(x) && is.null(dim(x))) x <- data.frame(value = x, stringsAsFactors = FALSE)
  if (!is.data.frame(x) && !is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, tmp, row.names = row.names, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) o2b_p7_stop("Failed to write CSV: ", path)
  invisible(path)
}

o2b_p7_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, tmp, compress = "xz")
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) o2b_p7_stop("Failed to save RDS: ", path)
  invisible(path)
}

o2b_p7_pkg_versions <- function(pkgs) {
  data.frame(
    package = pkgs,
    available = vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    version = vapply(pkgs, function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

o2b_p7_palette <- function() {
  c(
    red = "#B2182B",
    blue = "#2166AC",
    dark = "#1F2937",
    grey = "#9CA3AF",
    light_grey = "#F3F4F6",
    green = "#1B9E77",
    purple = "#7570B3",
    orange = "#D95F02",
    pale_red = "#F4A3A8",
    pale_blue = "#9ECAE1"
  )
}

o2b_p7_theme <- function(base_size = 11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "#1F2937", linewidth = 0.32),
      axis.ticks = ggplot2::element_line(color = "#1F2937", linewidth = 0.28),
      axis.text = ggplot2::element_text(color = "#1F2937", size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", color = "#111827"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, color = "#111827", size = base_size + 2, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(hjust = 0, color = "#4B5563", size = base_size - 1, margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(hjust = 0, color = "#6B7280", size = base_size - 2),
      legend.title = ggplot2::element_text(face = "bold", color = "#111827"),
      legend.text = ggplot2::element_text(color = "#111827"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "#1F2937", linewidth = 0.28),
      strip.text = ggplot2::element_text(face = "bold", color = "#111827", lineheight = 0.95),
      plot.margin = ggplot2::margin(10, 14, 10, 14)
    )
}

o2b_p7_fmt_p <- function(p, digits = 3) {
  p <- suppressWarnings(as.numeric(p))
  if (!is.finite(p)) return("P=NA")
  if (p < 1e-4) return("P<1e-4")
  paste0("P=", formatC(p, format = "f", digits = digits))
}

o2b_p7_fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

o2b_p7_save_plot <- function(plot, stem, fig_dir, width = 7, height = 5, dpi = 600) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(character()))
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(fig_dir, paste0(stem, c(".png", ".pdf", ".tiff")))
  formats <- c("png", "pdf", "tiff")
  saved <- character()
  errors <- list()
  save_one <- function(i) {
    tryCatch({
      if (i == 1L) {
        ggplot2::ggsave(files[i], plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
      } else if (i == 2L) {
        ggplot2::ggsave(files[i], plot, width = width, height = height, bg = "white", limitsize = FALSE, useDingbats = FALSE)
      } else {
        ggplot2::ggsave(files[i], plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE, compression = "lzw")
      }
      saved <<- c(saved, files[i])
      TRUE
    }, error = function(e) {
      errors[[length(errors) + 1L]] <<- data.frame(
        figure_stem = stem,
        format = formats[i],
        output_file = files[i],
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      FALSE
    })
  }
  invisible(vapply(seq_along(files), save_one, logical(1)))
  if (length(errors)) {
    audit_path <- file.path(fig_dir, "Part7_figure_export_errors.csv")
    err <- o2b_p7_bind_rows(errors)
    if (file.exists(audit_path)) {
      old <- tryCatch(utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
      err <- o2b_p7_bind_rows(list(old, err))
    }
    utils::write.csv(err, audit_path, row.names = FALSE, na = "")
    warning("Part7 figure export warning for ", stem, ": ", paste(unique(err$error_message), collapse = " | "), call. = FALSE)
  }
  invisible(saved)
}

o2b_p7_file_size <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

o2b_p7_is_gzip <- function(path) {
  if (!file.exists(path)) return(NA)
  con <- file(path, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  magic <- tryCatch(readBin(con, what = "raw", n = 2L), error = function(e) raw(0))
  if (length(magic) < 2L) return(FALSE)
  identical(as.integer(magic), c(0x1f, 0x8b))
}

o2b_p7_open_text <- function(path) {
  if (isTRUE(o2b_p7_is_gzip(path))) gzfile(path, "rt") else file(path, "rt")
}

o2b_p7_find_first <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) normalizePath(hit[1], winslash = "/", mustWork = FALSE) else NA_character_
}

o2b_p7_asset_patterns <- function() {
  list(
    matrix = "(^|[_-])matrix\\.mtx(\\.gz)?$",
    features = "(^|[_-])(features|genes)\\.tsv(\\.gz)?$",
    barcodes = "(^|[_-])barcodes\\.tsv(\\.gz)?$",
    positions = "(^|[_-])tissue_positions(_list)?\\.csv(\\.gz)?$",
    h5 = "(^|[_-])filtered_feature_bc_matrix\\.h5$"
  )
}

o2b_p7_match_asset <- function(files, pattern, preferred_dir = NULL) {
  files <- unique(as.character(files))
  files <- files[!is.na(files) & file.exists(files) & !dir.exists(files)]
  if (!length(files)) return(NA_character_)
  hit <- files[grepl(pattern, basename(files), ignore.case = TRUE, perl = TRUE)]
  if (!length(hit)) return(NA_character_)
  preferred <- rep(FALSE, length(hit))
  if (!is.null(preferred_dir) && length(preferred_dir) && nzchar(preferred_dir[1])) {
    preferred_dir <- normalizePath(preferred_dir[1], winslash = "/", mustWork = FALSE)
    preferred <- normalizePath(dirname(hit), winslash = "/", mustWork = FALSE) == preferred_dir
  }
  hit <- hit[order(!preferred, nchar(hit), hit)]
  normalizePath(hit[1], winslash = "/", mustWork = FALSE)
}

o2b_p7_resolve_visium_files <- function(sample_dir) {
  sample_dir <- normalizePath(sample_dir, winslash = "/", mustWork = FALSE)
  pats <- o2b_p7_asset_patterns()
  direct_files <- list.files(sample_dir, recursive = FALSE, full.names = TRUE, all.files = FALSE)
  all_files <- list.files(sample_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)

  matrix_file <- o2b_p7_match_asset(direct_files, pats$matrix)
  if (is.na(matrix_file)) matrix_file <- o2b_p7_match_asset(all_files, pats$matrix)
  matrix_dir <- if (!is.na(matrix_file) && file.exists(matrix_file)) dirname(matrix_file) else sample_dir
  matrix_dir_files <- list.files(matrix_dir, recursive = FALSE, full.names = TRUE, all.files = FALSE)

  features_file <- o2b_p7_match_asset(matrix_dir_files, pats$features, preferred_dir = matrix_dir)
  if (is.na(features_file)) features_file <- o2b_p7_match_asset(all_files, pats$features, preferred_dir = matrix_dir)
  barcodes_file <- o2b_p7_match_asset(matrix_dir_files, pats$barcodes, preferred_dir = matrix_dir)
  if (is.na(barcodes_file)) barcodes_file <- o2b_p7_match_asset(all_files, pats$barcodes, preferred_dir = matrix_dir)
  positions_file <- o2b_p7_match_asset(all_files, pats$positions)
  h5_file <- o2b_p7_match_asset(all_files, pats$h5)

  list(
    matrix = matrix_file,
    features = features_file,
    barcodes = barcodes_file,
    positions = positions_file,
    h5 = h5_file,
    matrix_base_dir = matrix_dir
  )
}

o2b_p7_find_sample_dirs <- function(standardized_dir, raw_dir = NULL) {
  roots <- character()
  if (!is.null(standardized_dir) && length(standardized_dir) && !is.na(standardized_dir[1]) && dir.exists(standardized_dir[1])) roots <- c(roots, standardized_dir[1])
  if (!is.null(raw_dir) && length(raw_dir) && !is.na(raw_dir[1]) && dir.exists(raw_dir[1])) roots <- c(roots, raw_dir[1])
  roots <- unique(normalizePath(roots, winslash = "/", mustWork = FALSE))
  if (!length(roots)) return(character())

  candidates <- character()
  source_priority <- integer()
  for (i in seq_along(roots)) {
    one <- unique(c(roots[i], list.dirs(roots[i], recursive = TRUE, full.names = TRUE)))
    one <- one[dir.exists(one)]
    candidates <- c(candidates, one)
    source_priority <- c(source_priority, rep.int(i, length(one)))
  }
  normalized_candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  keep_unique <- !duplicated(normalized_candidates)
  candidates <- candidates[keep_unique]
  source_priority <- source_priority[keep_unique]

  resolved <- lapply(candidates, o2b_p7_resolve_visium_files)
  complete <- vapply(resolved, function(f) {
    all(file.exists(unlist(f[c("matrix", "features", "barcodes")])))
  }, logical(1))
  direct_complete <- vapply(seq_along(candidates), function(i) {
    f <- resolved[[i]]
    if (!complete[i]) return(FALSE)
    asset_dirs <- normalizePath(dirname(unlist(f[c("matrix", "features", "barcodes")])), winslash = "/", mustWork = FALSE)
    candidate_dir <- normalizePath(candidates[i], winslash = "/", mustWork = FALSE)
    all(asset_dirs == candidate_dir)
  }, logical(1))

  # A dataset root can recursively see one section's assets and appear complete.
  # Only the ten official section identities are eligible as sample directories.
  map <- o2b_p7_official_gse251950_sample_map()
  sid_self <- o2b_p7_normalize_sample_id(basename(candidates))
  sid_parent <- o2b_p7_normalize_sample_id(basename(dirname(candidates)))
  geo_self <- o2b_p7_extract_geo_accession(basename(candidates))
  geo_parent <- o2b_p7_extract_geo_accession(basename(dirname(candidates)))

  idx_sid_self <- match(sid_self, map$sample_id)
  idx_sid_parent <- match(sid_parent, map$sample_id)
  idx_geo_self <- match(geo_self, map$geo_accession)
  idx_geo_parent <- match(geo_parent, map$geo_accession)
  idx_official <- idx_sid_self
  idx_official[is.na(idx_official)] <- idx_sid_parent[is.na(idx_official)]
  idx_official[is.na(idx_official)] <- idx_geo_self[is.na(idx_official)]
  idx_official[is.na(idx_official)] <- idx_geo_parent[is.na(idx_official)]
  official_candidate <- !is.na(idx_official)

  keep <- complete & official_candidate
  candidates <- candidates[keep]
  source_priority <- source_priority[keep]
  direct_complete <- direct_complete[keep]
  idx_official <- idx_official[keep]
  if (!length(candidates)) return(character())

  depth <- lengths(strsplit(normalizePath(candidates, winslash = "/", mustWork = FALSE), "/", fixed = TRUE))
  ord <- order(source_priority, !direct_complete, -depth, idx_official, candidates)
  candidates <- candidates[ord]
  idx_official <- idx_official[ord]
  keep_section <- !duplicated(idx_official)
  candidates <- candidates[keep_section]
  idx_official <- idx_official[keep_section]
  candidates <- candidates[order(idx_official, candidates)]
  unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
}

o2b_p7_input_file_inventory <- function(root_dir) {
  if (is.null(root_dir) || !length(root_dir) || is.na(root_dir[1]) || !dir.exists(root_dir[1])) return(data.frame())
  files <- list.files(root_dir[1], recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) return(data.frame())
  pats <- o2b_p7_asset_patterns()
  kind <- rep("other", length(files))
  for (nm in names(pats)) kind[grepl(pats[[nm]], basename(files), ignore.case = TRUE, perl = TRUE)] <- nm
  data.frame(
    root = normalizePath(root_dir[1], winslash = "/", mustWork = FALSE),
    path = normalizePath(files, winslash = "/", mustWork = FALSE),
    basename = basename(files),
    size_bytes = as.numeric(file.info(files)$size),
    asset_type = kind,
    geo_accession = o2b_p7_extract_geo_accession(basename(files)),
    stringsAsFactors = FALSE
  )
}

o2b_p7_unpack_local_archives <- function(raw_dir, force_unpack = FALSE, raw_tar_name = "GSE251950_RAW.tar") {
  rows <- list()
  if (!isTRUE(force_unpack) || is.null(raw_dir) || !length(raw_dir) || is.na(raw_dir[1])) return(data.frame())
  if (!dir.exists(raw_dir[1])) dir.create(raw_dir[1], recursive = TRUE, showWarnings = FALSE)
  outer <- file.path(raw_dir[1], raw_tar_name)
  sample_archives <- list.files(raw_dir[1], pattern = "^GSM[0-9]+.*\\.tar(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (!length(sample_archives) && file.exists(outer)) {
    err <- tryCatch({ utils::untar(outer, exdir = raw_dir[1]); NA_character_ }, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <- data.frame(archive = outer, destination = raw_dir[1], status = ifelse(is.na(err), "UNPACKED_OUTER", "FAIL_OUTER"), error = err, stringsAsFactors = FALSE)
    sample_archives <- list.files(raw_dir[1], pattern = "^GSM[0-9]+.*\\.tar(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  }
  for (archive in sample_archives) {
    dest <- file.path(raw_dir[1], sub("\\.tar(\\.gz)?$", "", basename(archive), ignore.case = TRUE))
    if (dir.exists(dest) && length(list.files(dest, recursive = TRUE, all.files = FALSE)) > 0L) {
      rows[[length(rows) + 1L]] <- data.frame(archive = archive, destination = dest, status = "SKIPPED_EXISTING", error = NA_character_, stringsAsFactors = FALSE)
      next
    }
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    err <- tryCatch({ utils::untar(archive, exdir = dest); NA_character_ }, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <- data.frame(archive = archive, destination = dest, status = ifelse(is.na(err), "UNPACKED_SAMPLE", "FAIL_SAMPLE"), error = err, stringsAsFactors = FALSE)
  }
  o2b_p7_bind_rows(rows)
}

o2b_p7_read_first_lines <- function(path, n = 5L) {
  if (!file.exists(path)) return(character())
  con <- o2b_p7_open_text(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  tryCatch(readLines(con, n = n, warn = FALSE), error = function(e) paste0("READ_ERROR: ", conditionMessage(e)))
}

o2b_p7_parse_mtx_header <- function(path) {
  lines <- o2b_p7_read_first_lines(path, n = 50L)
  if (!length(lines) || (length(lines) == 1L && grepl("^READ_ERROR:", lines[1]))) {
    return(data.frame(header_ok = FALSE, header_nrow = NA_integer_, header_ncol = NA_integer_, header_nnz = NA_integer_, header_error = paste(lines, collapse = ";"), stringsAsFactors = FALSE))
  }
  hdr <- lines[1]
  body <- lines[!grepl("^%", lines)]
  err <- NA_character_
  dims <- rep(NA_integer_, 3)
  if (length(body)) {
    vals <- suppressWarnings(as.integer(strsplit(trimws(body[1]), "\\s+")[[1]]))
    if (length(vals) >= 3 && all(is.finite(vals[1:3]))) dims <- vals[1:3] else err <- paste0("Could not parse dimension line: ", body[1])
  } else {
    err <- "No non-comment dimension line found"
  }
  data.frame(
    header_ok = isTRUE(grepl("^%%MatrixMarket", hdr)) && all(is.finite(dims)),
    header_nrow = dims[1],
    header_ncol = dims[2],
    header_nnz = dims[3],
    header_error = err,
    stringsAsFactors = FALSE
  )
}

o2b_p7_read_mtx <- function(path) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    return(list(ok = FALSE, object = NULL, error = "Matrix package not available", dim = c(NA_integer_, NA_integer_), nnzero = NA_integer_))
  }
  con <- o2b_p7_open_text(path)
  mat <- tryCatch(Matrix::readMM(con), error = function(e) e)
  try(close(con), silent = TRUE)
  if (inherits(mat, "error")) {
    mat2 <- tryCatch(Matrix::readMM(path), error = function(e) e)
    if (!inherits(mat2, "error")) mat <- mat2
  }
  if (inherits(mat, "error")) {
    return(list(ok = FALSE, object = NULL, error = conditionMessage(mat), dim = c(NA_integer_, NA_integer_), nnzero = NA_integer_))
  }
  mat <- tryCatch(methods::as(mat, "dgCMatrix"), error = function(e) Matrix::Matrix(mat, sparse = TRUE))
  list(ok = TRUE, object = mat, error = NA_character_, dim = dim(mat), nnzero = tryCatch(Matrix::nnzero(mat), error = function(e) NA_integer_))
}

o2b_p7_read_features <- function(path) {
  con <- o2b_p7_open_text(path)
  feat <- tryCatch(utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", sep = "\t", fill = TRUE, check.names = FALSE), error = function(e) e)
  try(close(con), silent = TRUE)
  if (inherits(feat, "error")) return(list(ok = FALSE, data = data.frame(), error = conditionMessage(feat), nrow = NA_integer_, ncol = NA_integer_))
  list(ok = TRUE, data = feat, error = NA_character_, nrow = nrow(feat), ncol = ncol(feat))
}

o2b_p7_read_barcodes <- function(path) {
  con <- o2b_p7_open_text(path)
  bc <- tryCatch(readLines(con, warn = FALSE), error = function(e) e)
  try(close(con), silent = TRUE)
  if (inherits(bc, "error")) return(list(ok = FALSE, data = character(), error = conditionMessage(bc), n = NA_integer_))
  bc <- trimws(as.character(bc))
  bc <- bc[nzchar(bc)]
  list(ok = TRUE, data = bc, error = NA_character_, n = length(bc))
}

o2b_p7_read_positions <- function(path) {
  if (!file.exists(path)) return(list(ok = FALSE, data = data.frame(), error = "positions missing", nrow = NA_integer_, ncol = NA_integer_, barcode_col = NA_character_))
  first <- o2b_p7_read_first_lines(path, n = 1L)
  has_header <- length(first) && grepl("barcode|in_tissue|array_row|pxl", first[1], ignore.case = TRUE)
  pos <- tryCatch(utils::read.csv(path, header = has_header, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) e)
  if (inherits(pos, "error")) return(list(ok = FALSE, data = data.frame(), error = conditionMessage(pos), nrow = NA_integer_, ncol = NA_integer_, barcode_col = NA_character_))
  if (!has_header) {
    cn <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    names(pos)[seq_len(min(ncol(pos), length(cn)))] <- cn[seq_len(min(ncol(pos), length(cn)))]
  }
  bc_col <- grep("^barcode$|barcode", names(pos), ignore.case = TRUE, value = TRUE)[1] %||% NA_character_
  list(ok = TRUE, data = pos, error = NA_character_, nrow = nrow(pos), ncol = ncol(pos), barcode_col = bc_col)
}

o2b_p7_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^[\"']+|[\"']+$", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("\\..*$", "", x)
  x
}

o2b_p7_marker_sets <- function() {
  list(
    OLFML2B = c("OLFML2B"),
    OLFML2A = c("OLFML2A"),
    OLFM2 = c("OLFM2"),
    Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "FAP", "PDGFRA", "SPARC"),
    Myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1", "COL1A1", "COL3A1", "POSTN"),
    Smooth_Muscle = c("ACTA2", "TAGLN", "MYH11", "CNN1", "DES", "MYL9", "TPM2"),
    Pericyte = c("RGS5", "CSPG4", "MCAM", "RBP1", "NOTCH3", "ABCC9", "KCNJ8", "PDGFRB"),
    Myeloid_Macrophage = c("LYZ", "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "MRC1"),
    Endothelial = c("PECAM1", "EMCN", "VWF", "KDR", "ESAM", "RAMP2", "PLVAP"),
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MSLN", "MUC1"),
    CD8_Cytotoxic = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "CCL5", "GZMB", "PRF1"),
    CAF_ECM = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN", "LUM", "FAP", "PDGFRA", "POSTN", "SPARC", "FN1"),
    ECM_Remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "POSTN", "THBS2", "MMP2", "MMP11", "SPARC", "DCN", "LUM"),
    TGFb_Response = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF", "INHBA", "PMEPA1"),
    Inflammatory_Fibroblast = c("CXCL12", "IL6", "CXCL14", "CFD", "C7", "APOD", "PI16", "COL14A1"),
    Myofibroblast_State = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1", "POSTN", "COL1A1", "COL3A1"),
    Pericyte_State = c("RGS5", "CSPG4", "MCAM", "RBP1", "NOTCH3", "ABCC9", "KCNJ8", "PDGFRB"),
    Smooth_Muscle_State = c("ACTA2", "TAGLN", "MYH11", "CNN1", "DES", "MYL9", "TPM2"),
    Myeloid_Macrophage_State = c("LYZ", "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "MRC1"),
    Endothelial_Angiogenic = c("PECAM1", "VWF", "KDR", "ESAM", "RAMP2", "PLVAP", "ANGPT2", "KDR"),
    Epithelial_Differentiation = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2", "KRT20")
  )
}

o2b_p7_source_features <- function() {
  c("Fibroblast", "Myofibroblast", "Smooth_Muscle", "Pericyte", "Myeloid_Macrophage", "Endothelial", "Epithelial")
}

o2b_p7_state_features <- function() {
  c("CAF_ECM", "ECM_Remodeling", "TGFb_Response", "Inflammatory_Fibroblast",
    "Myofibroblast_State", "Pericyte_State", "Smooth_Muscle_State",
    "Myeloid_Macrophage_State", "Endothelial_Angiogenic",
    "Epithelial_Differentiation", "CD8_Cytotoxic", "Immune_Exclusion_Index")
}

o2b_p7_primary_features <- function() {
  c("Fibroblast", "Myofibroblast", "CAF_ECM", "ECM_Remodeling", "TGFb_Response")
}

o2b_p7_all_context_features <- function() {
  unique(c(o2b_p7_source_features(), o2b_p7_state_features()))
}


o2b_p7_feature_col_diag <- function(feat) {
  ms <- o2b_p7_marker_sets()
  all_markers <- unique(unlist(ms, use.names = FALSE))
  cols <- seq_len(min(6L, ncol(feat)))
  rows <- lapply(cols, function(j) {
    g <- o2b_p7_clean_gene(feat[[j]])
    data.frame(
      feature_col_index = j,
      n_nonmissing = sum(nzchar(g)),
      contains_OLFML2B_exact = any(g == "OLFML2B"),
      contains_OLFML2A_exact = any(g == "OLFML2A"),
      contains_OLFM2_exact = any(g == "OLFM2"),
      marker_overlap_total = sum(unique(g) %in% all_markers),
      fibroblast_marker_overlap = sum(unique(g) %in% ms$Fibroblast),
      myofibroblast_marker_overlap = sum(unique(g) %in% ms$Myofibroblast),
      pericyte_marker_overlap = sum(unique(g) %in% ms$Pericyte),
      myeloid_marker_overlap = sum(unique(g) %in% ms$Myeloid_Macrophage),
      epithelial_marker_overlap = sum(unique(g) %in% ms$Epithelial),
      examples = paste(head(unique(g[nzchar(g)]), 5), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  o2b_p7_bind_rows(rows)
}


o2b_p7_pick_gene_col <- function(feat_diag) {
  if (!nrow(feat_diag)) return(NA_integer_)
  score <- 1000 * as.integer(feat_diag$contains_OLFML2B_exact %in% TRUE) +
    100 * feat_diag$marker_overlap_total +
    feat_diag$n_nonmissing / 1000000
  feat_diag$feature_col_index[which.max(score)]
}

o2b_p7_strip_suffix <- function(x) sub("-[0-9]+$", "", trimws(as.character(x)))

o2b_p7_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - m) / s
}

o2b_p7_gene_vector <- function(logmat, gene_vec, gene) {
  gene <- toupper(gene)
  idx <- which(gene_vec == gene)
  if (!length(idx)) return(rep(NA_real_, ncol(logmat)))
  if (length(idx) == 1L) return(as.numeric(logmat[idx, ]))
  as.numeric(Matrix::colMeans(logmat[idx, , drop = FALSE], na.rm = TRUE))
}

o2b_p7_signature_score <- function(logmat, gene_vec, genes) {
  genes <- unique(toupper(genes))
  idx <- which(gene_vec %in% genes)
  if (!length(idx)) return(rep(NA_real_, ncol(logmat)))
  as.numeric(Matrix::colMeans(logmat[idx, , drop = FALSE], na.rm = TRUE))
}

o2b_p7_load_one_sample <- function(sample_dir, log_file = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) o2b_p7_stop("Matrix package is required for Part7.")
  sample_id <- basename(sample_dir)
  f <- o2b_p7_resolve_visium_files(sample_dir)

  hdr <- if (file.exists(f$matrix)) o2b_p7_parse_mtx_header(f$matrix) else data.frame(header_ok = FALSE, header_nrow = NA_integer_, header_ncol = NA_integer_, header_nnz = NA_integer_, header_error = "matrix missing")
  mat_res <- if (file.exists(f$matrix)) o2b_p7_read_mtx(f$matrix) else list(ok = FALSE, object = NULL, error = "matrix missing", dim = c(NA_integer_, NA_integer_), nnzero = NA_integer_)
  feat_res <- if (file.exists(f$features)) o2b_p7_read_features(f$features) else list(ok = FALSE, data = data.frame(), error = "features missing", nrow = NA_integer_, ncol = NA_integer_)
  bc_res <- if (file.exists(f$barcodes)) o2b_p7_read_barcodes(f$barcodes) else list(ok = FALSE, data = character(), error = "barcodes missing", n = NA_integer_)
  pos_res <- if (file.exists(f$positions)) o2b_p7_read_positions(f$positions) else list(ok = FALSE, data = data.frame(), error = "positions missing", nrow = NA_integer_, ncol = NA_integer_, barcode_col = NA_character_)

  feat_diag <- if (isTRUE(feat_res$ok)) o2b_p7_feature_col_diag(feat_res$data) else data.frame()
  gene_col <- o2b_p7_pick_gene_col(feat_diag)
  gene_vec <- if (isTRUE(feat_res$ok) && is.finite(gene_col)) o2b_p7_clean_gene(feat_res$data[[gene_col]]) else character()

  orientation <- "NOT_EVALUATED"
  if (isTRUE(mat_res$ok) && isTRUE(feat_res$ok) && isTRUE(bc_res$ok)) {
    dd <- mat_res$dim
    if (identical(as.integer(dd[1]), as.integer(feat_res$nrow)) && identical(as.integer(dd[2]), as.integer(bc_res$n))) {
      orientation <- "FEATURES_BY_BARCODES_OK"
    } else if (identical(as.integer(dd[2]), as.integer(feat_res$nrow)) && identical(as.integer(dd[1]), as.integer(bc_res$n))) {
      orientation <- "TRANSPOSE_NEEDED"
    } else {
      orientation <- "NO_DIM_MATCH"
    }
  }

  pos_bc <- character()
  if (isTRUE(pos_res$ok) && !is.na(pos_res$barcode_col) && pos_res$barcode_col %in% names(pos_res$data)) {
    pos_bc <- trimws(as.character(pos_res$data[[pos_res$barcode_col]]))
  }
  bc <- bc_res$data
  exact_overlap <- if (length(bc) && length(pos_bc)) length(intersect(bc, pos_bc)) else NA_integer_
  stripped_overlap <- if (length(bc) && length(pos_bc)) length(intersect(o2b_p7_strip_suffix(bc), o2b_p7_strip_suffix(pos_bc))) else NA_integer_

  status <- "OK_DIRECT_MTX_DIAGNOSTIC_LOADER"
  if (!file.exists(f$matrix)) status <- "FAIL_MATRIX_MISSING"
  else if (!file.exists(f$features)) status <- "FAIL_FEATURES_MISSING"
  else if (!file.exists(f$barcodes)) status <- "FAIL_BARCODES_MISSING"
  else if (!isTRUE(mat_res$ok)) status <- "FAIL_READMM"
  else if (!isTRUE(feat_res$ok)) status <- "FAIL_FEATURES_READ"
  else if (!isTRUE(bc_res$ok)) status <- "FAIL_BARCODES_READ"
  else if (identical(orientation, "NO_DIM_MATCH")) status <- "FAIL_MATRIX_FEATURE_BARCODE_DIM_MISMATCH"
  else if (!any(gene_vec == "OLFML2B")) status <- "FAIL_NO_OLFML2B_IN_FEATURES"
  else if (!(exact_overlap > 0 || stripped_overlap > 0)) status <- "WARN_SCORE_OK_NO_POSITION_OVERLAP"

  qc <- data.frame(
    sample_id = sample_id,
    sample_dir = sample_dir,
    matrix_path = f$matrix %||% NA_character_,
    features_path = f$features %||% NA_character_,
    barcodes_path = f$barcodes %||% NA_character_,
    positions_path = f$positions %||% NA_character_,
    has_matrix_mtx = file.exists(f$matrix),
    has_features = file.exists(f$features),
    has_barcodes = file.exists(f$barcodes),
    has_spatial_positions = file.exists(f$positions),
    matrix_is_gzip_magic = o2b_p7_is_gzip(f$matrix),
    features_is_gzip_magic = o2b_p7_is_gzip(f$features),
    barcodes_is_gzip_magic = o2b_p7_is_gzip(f$barcodes),
    matrix_header_ok = hdr$header_ok[1],
    matrix_header_nrow = hdr$header_nrow[1],
    matrix_header_ncol = hdr$header_ncol[1],
    matrix_header_nnz = hdr$header_nnz[1],
    readmm_ok = isTRUE(mat_res$ok),
    readmm_error = mat_res$error %||% NA_character_,
    mat_nrow = as.integer(mat_res$dim[1]),
    mat_ncol = as.integer(mat_res$dim[2]),
    mat_nnz = as.integer(mat_res$nnzero),
    features_read_ok = isTRUE(feat_res$ok),
    features_nrow = as.integer(feat_res$nrow),
    features_ncol = as.integer(feat_res$ncol),
    barcodes_read_ok = isTRUE(bc_res$ok),
    barcodes_n = as.integer(bc_res$n),
    positions_read_ok = isTRUE(pos_res$ok),
    positions_nrow = as.integer(pos_res$nrow),
    positions_ncol = as.integer(pos_res$ncol),
    barcode_position_overlap = as.integer(exact_overlap),
    barcode_position_overlap_suffix_stripped = as.integer(stripped_overlap),
    gene_col_index_selected = as.integer(gene_col),
    target_OLFML2B_found = any(gene_vec == "OLFML2B"),
    marker_overlap_selected_gene_col = if (length(gene_vec)) sum(unique(gene_vec) %in% unique(unlist(o2b_p7_marker_sets(), use.names = FALSE))) else 0L,
    orientation_call = orientation,
    parser_status = status,
    loader_notes = "Diagnostic loader integrated into production Part7 v1.9.2",
    stringsAsFactors = FALSE
  )

  if (!startsWith(status, "OK") && !startsWith(status, "WARN_SCORE_OK")) {
    return(list(ok = FALSE, scores = data.frame(), qc = qc, feature_diag = transform(feat_diag, sample_id = sample_id), coverage = data.frame()))
  }

  mat <- mat_res$object
  if (identical(orientation, "TRANSPOSE_NEEDED")) mat <- Matrix::t(mat)
  mat <- tryCatch(methods::as(mat, "dgCMatrix"), error = function(e) Matrix::Matrix(mat, sparse = TRUE))

  lib <- Matrix::colSums(mat)
  lib[!is.finite(lib) | lib <= 0] <- NA_real_
  norm <- Matrix::t(Matrix::t(mat) / lib * 10000)
  norm@x[!is.finite(norm@x)] <- 0
  logmat <- log1p(norm)

  ms <- o2b_p7_marker_sets()
  n_features_spot <- Matrix::colSums(mat > 0)
  score <- data.frame(
    spot_id = bc,
    total_counts = as.numeric(lib),
    n_features = as.numeric(n_features_spot),
    OLFML2B = o2b_p7_gene_vector(logmat, gene_vec, "OLFML2B"),
    OLFML2A = o2b_p7_gene_vector(logmat, gene_vec, "OLFML2A"),
    OLFM2 = o2b_p7_gene_vector(logmat, gene_vec, "OLFM2"),
    Fibroblast = o2b_p7_signature_score(logmat, gene_vec, ms$Fibroblast),
    Myofibroblast = o2b_p7_signature_score(logmat, gene_vec, ms$Myofibroblast),
    Smooth_Muscle = o2b_p7_signature_score(logmat, gene_vec, ms$Smooth_Muscle),
    Pericyte = o2b_p7_signature_score(logmat, gene_vec, ms$Pericyte),
    Myeloid_Macrophage = o2b_p7_signature_score(logmat, gene_vec, ms$Myeloid_Macrophage),
    Endothelial = o2b_p7_signature_score(logmat, gene_vec, ms$Endothelial),
    Epithelial = o2b_p7_signature_score(logmat, gene_vec, ms$Epithelial),
    CD8_Cytotoxic = o2b_p7_signature_score(logmat, gene_vec, ms$CD8_Cytotoxic),
    CAF_ECM = o2b_p7_signature_score(logmat, gene_vec, ms$CAF_ECM),
    ECM_Remodeling = o2b_p7_signature_score(logmat, gene_vec, ms$ECM_Remodeling),
    TGFb_Response = o2b_p7_signature_score(logmat, gene_vec, ms$TGFb_Response),
    Inflammatory_Fibroblast = o2b_p7_signature_score(logmat, gene_vec, ms$Inflammatory_Fibroblast),
    Myofibroblast_State = o2b_p7_signature_score(logmat, gene_vec, ms$Myofibroblast_State),
    Pericyte_State = o2b_p7_signature_score(logmat, gene_vec, ms$Pericyte_State),
    Smooth_Muscle_State = o2b_p7_signature_score(logmat, gene_vec, ms$Smooth_Muscle_State),
    Myeloid_Macrophage_State = o2b_p7_signature_score(logmat, gene_vec, ms$Myeloid_Macrophage_State),
    Endothelial_Angiogenic = o2b_p7_signature_score(logmat, gene_vec, ms$Endothelial_Angiogenic),
    Epithelial_Differentiation = o2b_p7_signature_score(logmat, gene_vec, ms$Epithelial_Differentiation),
    stringsAsFactors = FALSE
  )

  # Per-section derived axes. OLFML2B is excluded from every context signature.
  score$Immune_Exclusion_Index <- o2b_p7_z(score$CAF_ECM) + o2b_p7_z(score$ECM_Remodeling) +
    o2b_p7_z(score$TGFb_Response) - o2b_p7_z(score$CD8_Cytotoxic)
  score$TME_exclusion_score <- score$Immune_Exclusion_Index
  score$OLFML2B_Stromal_composite_score <- o2b_p7_z(score$OLFML2B) +
    o2b_p7_z(score$Fibroblast) + o2b_p7_z(score$Myofibroblast) +
    o2b_p7_z(score$ECM_Remodeling) - o2b_p7_z(score$Epithelial)
  score$normalization <- "log1p_CPM_10000_full_matrix"
  score$sample_id <- sample_id

  coverage <- o2b_p7_bind_rows(lapply(names(ms), function(nm) {
    genes <- unique(toupper(ms[[nm]]))
    measured <- intersect(genes, unique(gene_vec))
    data.frame(
      sample_id = sample_id,
      signature = nm,
      n_genes_defined = length(genes),
      n_genes_measured = length(measured),
      coverage_fraction = if (length(genes)) length(measured) / length(genes) else NA_real_,
      measured_genes = paste(measured, collapse = ";"),
      status = ifelse(nm %in% c("OLFML2B", "OLFML2A", "OLFM2"),
                      ifelse(length(measured) >= 1L, "MEASURED", "NOT_MEASURED"),
                      ifelse(length(measured) >= 3L && length(measured) / length(genes) >= 0.4,
                             "PASS_STATE_PROGRAM_COVERAGE",
                             ifelse(length(measured) >= 2L && length(measured) / length(genes) >= 0.25,
                                    "PASS_SOURCE_MARKER_COVERAGE", "LOW_COVERAGE"))),
      stringsAsFactors = FALSE
    )
  }))

  # Attach positions.
  score$barcode_match_key <- score$spot_id
  if (isTRUE(pos_res$ok) && !is.na(pos_res$barcode_col) && pos_res$barcode_col %in% names(pos_res$data)) {
    pos <- pos_res$data
    names(pos) <- make.names(names(pos), unique = TRUE)
    bc_col <- make.names(pos_res$barcode_col)
    pos$barcode_match_key <- trimws(as.character(pos[[bc_col]]))
    if (exact_overlap <= 0 && stripped_overlap > 0) {
      score$barcode_match_key <- o2b_p7_strip_suffix(score$spot_id)
      pos$barcode_match_key <- o2b_p7_strip_suffix(pos$barcode_match_key)
    }

    pxl_row_col <- grep("pxl.*row|row.*fullres", names(pos), ignore.case = TRUE, value = TRUE)[1]
    pxl_col_col <- grep("pxl.*col|col.*fullres", names(pos), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(pxl_row_col)) pxl_row_col <- grep("array.*row", names(pos), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(pxl_col_col)) pxl_col_col <- grep("array.*col", names(pos), ignore.case = TRUE, value = TRUE)[1]
    keep <- unique(c("barcode_match_key", pxl_row_col, pxl_col_col, grep("in.tissue|in_tissue", names(pos), ignore.case = TRUE, value = TRUE)[1]))
    keep <- keep[!is.na(keep) & keep %in% names(pos)]
    pos2 <- pos[, keep, drop = FALSE]
    if (!"pxl_row" %in% names(pos2)) names(pos2)[names(pos2) == pxl_row_col] <- "pxl_row"
    if (!"pxl_col" %in% names(pos2)) names(pos2)[names(pos2) == pxl_col_col] <- "pxl_col"
    score <- merge(score, pos2, by = "barcode_match_key", all.x = TRUE, sort = FALSE)
  }
  if (!"pxl_row" %in% names(score)) score$pxl_row <- NA_real_
  if (!"pxl_col" %in% names(score)) score$pxl_col <- NA_real_

  qc$n_scores = nrow(score)
  qc$n_spots_with_coordinates = sum(is.finite(suppressWarnings(as.numeric(score$pxl_row))) & is.finite(suppressWarnings(as.numeric(score$pxl_col))))
  qc$olfml2b_detected_spots = sum(score$OLFML2B > 0, na.rm = TRUE)
  qc$olfml2b_detected_fraction = mean(score$OLFML2B > 0, na.rm = TRUE)
  qc$fibroblast_marker_available = sum(unique(gene_vec) %in% ms$Fibroblast)
  qc$myofibroblast_marker_available = sum(unique(gene_vec) %in% ms$Myofibroblast)
  qc$cd8_marker_available = sum(unique(gene_vec) %in% ms$CD8_Cytotoxic)

  score$barcode_match_key <- NULL
  list(ok = TRUE, scores = score, qc = qc, feature_diag = transform(feat_diag, sample_id = sample_id), coverage = coverage)
}

o2b_p7_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) {
    return(c(n = sum(ok), rho = NA_real_, p_value = NA_real_))
  }
  tt <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(tt$estimate), p_value = tt$p.value)
}


o2b_p7_define_olfml2b_groups <- function(x, min_positive_spots = 30L, high_quantile = 0.75, low_quantile = 0.25) {
  x <- suppressWarnings(as.numeric(x))
  n <- sum(is.finite(x))
  pos <- is.finite(x) & x > 0
  n_pos <- sum(pos)
  det_frac <- if (n > 0) n_pos / n else NA_real_

  out <- rep(NA_character_, length(x))
  method <- "NOT_EVALUABLE"
  high_cutoff <- NA_real_
  low_cutoff <- NA_real_

  if (n < 20L || n_pos < min_positive_spots) {
    return(list(
      group = out,
      method = paste0("NOT_EVALUABLE_POSITIVE_SPOTS_LT_", min_positive_spots),
      n_total = n,
      n_positive = n_pos,
      detected_fraction = det_frac,
      high_cutoff = high_cutoff,
      low_cutoff = low_cutoff
    ))
  }

  # Standard case: enough positive spots and the 75th percentile is above zero.
  q_all <- stats::quantile(x[is.finite(x)], c(low_quantile, high_quantile), na.rm = TRUE, names = FALSE)
  if (is.finite(q_all[2]) && q_all[2] > 0) {
    out[is.finite(x) & x >= q_all[2]] <- "OLFML2B_high"
    out[is.finite(x) & x <= q_all[1]] <- "OLFML2B_low"
    method <- "ALL_SPOTS_TOP_BOTTOM_QUARTILE"
    high_cutoff <- q_all[2]
    low_cutoff <- q_all[1]
  } else {
    # Sparse case: top-quartile over all spots would be zero and would misclassify
    # many undetected spots as OLFML2B-high. Define high among positive spots only,
    # while low is OLFML2B-undetected.
    q_pos <- stats::quantile(x[pos], high_quantile, na.rm = TRUE, names = FALSE)
    if (!is.finite(q_pos) || q_pos <= 0) q_pos <- min(x[pos], na.rm = TRUE)
    out[is.finite(x) & x >= q_pos & x > 0] <- "OLFML2B_high"
    out[is.finite(x) & x <= 0] <- "OLFML2B_low"
    method <- "SPARSE_POSITIVE_ONLY_TOP_QUARTILE_VS_UNDETECTED"
    high_cutoff <- q_pos
    low_cutoff <- 0
  }

  list(
    group = out,
    method = method,
    n_total = n,
    n_positive = n_pos,
    detected_fraction = det_frac,
    high_cutoff = high_cutoff,
    low_cutoff = low_cutoff
  )
}

o2b_p7_define_binary_high <- function(x, min_positive_spots = 30L) {
  g <- o2b_p7_define_olfml2b_groups(x, min_positive_spots = min_positive_spots)
  list(
    high = g$group %in% "OLFML2B_high",
    method = g$method,
    high_cutoff = g$high_cutoff,
    detected_fraction = g$detected_fraction,
    n_positive = g$n_positive
  )
}


o2b_p7_corr_tables <- function(scores) {
  features <- intersect(c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Epithelial_Tumor", "Myeloid", "Checkpoint_Exhaustion", "Immune_Exclusion_Index", "TME_exclusion_score", "OLFML2B_TME_composite_score"), names(scores))
  rows <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    for (f in features) {
      cc <- o2b_p7_spearman(d$OLFML2B, d[[f]])
      rows[[length(rows) + 1L]] <- data.frame(sample_id = sid, feature = f, n = as.integer(cc["n"]), rho = as.numeric(cc["rho"]), p_value = as.numeric(cc["p_value"]), status = ifelse(is.finite(cc["rho"]), "OK", "NOT_EVALUABLE"), analysis_role = "spot_level_continuous_spearman", normalization = "log1p_CPM_10000", stringsAsFactors = FALSE)
    }
  }
  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) out$fdr <- stats::ave(out$p_value, out$feature, FUN = function(x) stats::p.adjust(x, method = "BH"))
  out
}

o2b_p7_high_low_tables <- function(scores) {
  features <- intersect(c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Epithelial_Tumor", "Myeloid", "Immune_Exclusion_Index"), names(scores))
  rows <- list()
  threshold_rows <- list()

  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    grp <- o2b_p7_define_olfml2b_groups(d$OLFML2B, min_positive_spots = 30L)
    d$spot_group <- grp$group

    threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
      sample_id = sid,
      grouping_method = grp$method,
      n_total = grp$n_total,
      n_positive = grp$n_positive,
      detected_fraction = grp$detected_fraction,
      high_cutoff = grp$high_cutoff,
      low_cutoff = grp$low_cutoff,
      n_high = sum(d$spot_group == "OLFML2B_high", na.rm = TRUE),
      n_low = sum(d$spot_group == "OLFML2B_low", na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    dd <- d[d$spot_group %in% c("OLFML2B_high", "OLFML2B_low"), , drop = FALSE]
    for (f in features) {
      a <- dd[[f]][dd$spot_group == "OLFML2B_high"]
      b <- dd[[f]][dd$spot_group == "OLFML2B_low"]
      pv <- if (sum(is.finite(a)) >= 10 && sum(is.finite(b)) >= 10) {
        tryCatch(stats::wilcox.test(a, b)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sid,
        grouping = grp$method,
        analysis_role = "spatial_spot_level_tme_shift_adaptive_olfml2b_high",
        feature = f,
        group_a = "OLFML2B_high",
        group_b = "OLFML2B_low",
        n_a = sum(is.finite(a)),
        n_b = sum(is.finite(b)),
        mean_a = mean(a, na.rm = TRUE),
        mean_b = mean(b, na.rm = TRUE),
        delta_a_minus_b = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
        p_value = pv,
        status = ifelse(is.finite(pv), "OK", "NOT_EVALUABLE"),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) out$fdr <- stats::ave(out$p_value, out$feature, FUN = function(x) stats::p.adjust(x, method = "BH"))
  attr(out, "thresholds") <- o2b_p7_bind_rows(threshold_rows)
  out
}


o2b_p7_sign_test_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x != 0]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  k <- sum(x > 0)
  tryCatch(stats::binom.test(k, n, p = 0.5)$p.value, error = function(e) NA_real_)
}

o2b_p7_wilcox_zero_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 3L || stats::sd(x, na.rm = TRUE) <= 0) return(NA_real_)
  tryCatch(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
}

o2b_p7_loo_direction_stability <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 4L) {
    return(list(
      loo_min_median = NA_real_,
      loo_max_median = NA_real_,
      loo_direction_stability = NA_real_
    ))
  }
  full_med <- stats::median(x, na.rm = TRUE)
  if (!is.finite(full_med) || full_med == 0) full_sign <- sign(mean(x, na.rm = TRUE)) else full_sign <- sign(full_med)
  loo <- vapply(seq_along(x), function(i) stats::median(x[-i], na.rm = TRUE), numeric(1))
  list(
    loo_min_median = min(loo, na.rm = TRUE),
    loo_max_median = max(loo, na.rm = TRUE),
    loo_direction_stability = mean(sign(loo) == full_sign, na.rm = TRUE)
  )
}

o2b_p7_summarise_section_effects <- function(effects, effect_col = "effect", group_cols = c("effect_source", "feature")) {
  if (!is.data.frame(effects) || !nrow(effects)) return(data.frame())
  effects <- o2b_p7_attach_official_sample_metadata(effects)
  effects[[effect_col]] <- suppressWarnings(as.numeric(effects[[effect_col]]))
  effects <- effects[is.finite(effects[[effect_col]]), , drop = FALSE]
  if (!nrow(effects)) return(data.frame())

  key <- do.call(paste, c(effects[, group_cols, drop = FALSE], sep = "\r"))
  rows <- lapply(split(seq_len(nrow(effects)), key), function(ii) {
    d <- effects[ii, , drop = FALSE]
    section_x <- d[[effect_col]]
    patient_id <- if ("patient_id" %in% names(d)) as.character(d$patient_id) else as.character(d$sample_id)
    patient_id[is.na(patient_id) | !nzchar(patient_id)] <- as.character(d$sample_id[is.na(patient_id) | !nzchar(patient_id)])
    patient_effect <- stats::aggregate(section_x, by = list(patient_id = patient_id), FUN = function(z) stats::median(z, na.rm = TRUE))
    x <- patient_effect$x
    loo <- o2b_p7_loo_direction_stability(x)
    data.frame(
      effect_source = if ("effect_source" %in% names(d)) d$effect_source[1] else NA_character_,
      threshold_method = if ("threshold_method" %in% names(d)) d$threshold_method[1] else NA_character_,
      feature = if ("feature" %in% names(d)) d$feature[1] else NA_character_,
      n_sections = length(unique(d$sample_id)),
      n_patients = length(unique(patient_id)),
      median_effect = stats::median(x, na.rm = TRUE),
      q1_effect = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
      q3_effect = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
      positive_direction_count = sum(x > 0, na.rm = TRUE),
      negative_direction_count = sum(x < 0, na.rm = TRUE),
      zero_direction_count = sum(x == 0, na.rm = TRUE),
      positive_direction_fraction = mean(x > 0, na.rm = TRUE),
      sign_test_p = o2b_p7_sign_test_p(x),
      wilcoxon_p_against_zero = o2b_p7_wilcox_zero_p(x),
      loo_min_median = loo$loo_min_median,
      loo_max_median = loo$loo_max_median,
      loo_direction_stability = loo$loo_direction_stability,
      inference_unit = "patient_level_after_within_patient_section_aggregation",
      interpretation_boundary = "Spot effects are summarized by section and then collapsed within patient; GC6 primary/metastasis sections are not treated as independent patients.",
      stringsAsFactors = FALSE
    )
  })
  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) {
    out$sign_test_fdr <- stats::p.adjust(out$sign_test_p, method = "BH")
    out$wilcoxon_fdr_against_zero <- stats::p.adjust(out$wilcoxon_p_against_zero, method = "BH")
  }
  out
}

o2b_p7_continuous_sample_level_effects <- function(corr, spatial) {
  rows <- list()

  if (is.data.frame(corr) && nrow(corr)) {
    d <- corr
    d$effect_source <- "spot_level_spearman"
    d$effect <- suppressWarnings(as.numeric(d$rho))
    d$n_observations <- suppressWarnings(as.integer(d$n))
    d$primary_role <- "continuous_primary"
    rows[[length(rows) + 1L]] <- d[, intersect(c("sample_id", "feature", "effect_source", "effect", "p_value", "fdr", "n_observations", "primary_role", "analysis_role"), names(d)), drop = FALSE]
  }

  if (is.list(spatial) && is.data.frame(spatial$bivar) && nrow(spatial$bivar)) {
    d <- spatial$bivar
    d$effect_source <- "neighbor_mean_colocalization"
    d$feature <- d$neighbor_feature
    d$effect <- suppressWarnings(as.numeric(d$rho_target_vs_neighbor_mean))
    d$n_observations <- suppressWarnings(as.integer(d$n))
    d$primary_role <- "continuous_primary"
    rows[[length(rows) + 1L]] <- d[, intersect(c("sample_id", "feature", "effect_source", "effect", "p_value", "fdr", "n_observations", "primary_role"), names(d)), drop = FALSE]
  }

  if (is.list(spatial) && is.data.frame(spatial$auto) && nrow(spatial$auto)) {
    d <- spatial$auto
    d$effect_source <- "spatial_autocorrelation_morans_i"
    d$effect <- suppressWarnings(as.numeric(d$morans_i))
    d$p_value <- NA_real_
    d$fdr <- NA_real_
    d$n_observations <- suppressWarnings(as.integer(d$n_spots))
    d$primary_role <- "continuous_primary"
    rows[[length(rows) + 1L]] <- d[, intersect(c("sample_id", "feature", "effect_source", "effect", "p_value", "fdr", "n_observations", "primary_role"), names(d)), drop = FALSE]
  }

  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) {
    out <- out[is.finite(suppressWarnings(as.numeric(out$effect))), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

o2b_p7_threshold_groups <- function(x, method, min_positive_spots = 30L) {
  x <- suppressWarnings(as.numeric(x))
  pos <- is.finite(x) & x > 0
  out <- rep(NA_character_, length(x))
  high_cutoff <- NA_real_
  low_cutoff <- NA_real_
  n_pos <- sum(pos)
  n_total <- sum(is.finite(x))

  if (method == "detected_vs_undetected") {
    out[is.finite(x) & x > 0] <- "high"
    out[is.finite(x) & x <= 0] <- "low"
    high_cutoff <- .Machine$double.eps
    low_cutoff <- 0
  } else if (method == "positive_only_median") {
    if (n_pos >= min_positive_spots) {
      cut <- stats::median(x[pos], na.rm = TRUE)
      out[pos & x > cut] <- "high"
      out[pos & x <= cut] <- "low"
      high_cutoff <- cut
      low_cutoff <- cut
    }
  } else if (method == "positive_only_upper_quartile") {
    if (n_pos >= min_positive_spots) {
      qs <- stats::quantile(x[pos], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      out[pos & x >= qs[2]] <- "high"
      out[pos & x <= qs[1]] <- "low"
      high_cutoff <- qs[2]
      low_cutoff <- qs[1]
    }
  } else if (method == "adaptive_sparse_aware") {
    g <- o2b_p7_define_olfml2b_groups(x, min_positive_spots = min_positive_spots)
    out[g$group == "OLFML2B_high"] <- "high"
    out[g$group == "OLFML2B_low"] <- "low"
    high_cutoff <- g$high_cutoff
    low_cutoff <- g$low_cutoff
  }

  list(
    group = out,
    method = method,
    n_total = n_total,
    n_positive = n_pos,
    detected_fraction = if (n_total > 0) n_pos / n_total else NA_real_,
    high_cutoff = high_cutoff,
    low_cutoff = low_cutoff,
    n_high = sum(out == "high", na.rm = TRUE),
    n_low = sum(out == "low", na.rm = TRUE)
  )
}

o2b_p7_threshold_sensitivity <- function(scores, min_group_spots = 10L) {
  features <- intersect(c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Epithelial_Tumor", "Myeloid", "Immune_Exclusion_Index"), names(scores))
  methods <- c("detected_vs_undetected", "positive_only_median", "positive_only_upper_quartile", "adaptive_sparse_aware")
  rows <- list()
  audit_rows <- list()

  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    for (m in methods) {
      grp <- o2b_p7_threshold_groups(d$OLFML2B, method = m, min_positive_spots = 30L)
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        sample_id = sid,
        threshold_method = m,
        n_total = grp$n_total,
        n_positive = grp$n_positive,
        detected_fraction = grp$detected_fraction,
        high_cutoff = grp$high_cutoff,
        low_cutoff = grp$low_cutoff,
        n_high = grp$n_high,
        n_low = grp$n_low,
        stringsAsFactors = FALSE
      )

      if (grp$n_high < min_group_spots || grp$n_low < min_group_spots) {
        for (f in features) {
          rows[[length(rows) + 1L]] <- data.frame(
            sample_id = sid, threshold_method = m, feature = f,
            n_high = grp$n_high, n_low = grp$n_low,
            mean_high = NA_real_, mean_low = NA_real_, delta_high_minus_low = NA_real_,
            p_value = NA_real_, status = "NOT_EVALUABLE_INSUFFICIENT_GROUP_SIZE",
            inference_role = "threshold_sensitivity_only",
            stringsAsFactors = FALSE
          )
        }
        next
      }

      for (f in features) {
        a <- suppressWarnings(as.numeric(d[[f]][grp$group == "high"]))
        b <- suppressWarnings(as.numeric(d[[f]][grp$group == "low"]))
        pv <- if (sum(is.finite(a)) >= min_group_spots && sum(is.finite(b)) >= min_group_spots) {
          tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
        } else {
          NA_real_
        }
        rows[[length(rows) + 1L]] <- data.frame(
          sample_id = sid,
          threshold_method = m,
          feature = f,
          n_high = sum(is.finite(a)),
          n_low = sum(is.finite(b)),
          mean_high = mean(a, na.rm = TRUE),
          mean_low = mean(b, na.rm = TRUE),
          delta_high_minus_low = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
          p_value = pv,
          status = ifelse(is.finite(pv), "OK", "NOT_EVALUABLE"),
          inference_role = "threshold_sensitivity_only",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  by_sample <- o2b_p7_bind_rows(rows)
  audit <- o2b_p7_bind_rows(audit_rows)
  if (nrow(by_sample)) {
    by_sample$fdr_within_table <- stats::p.adjust(by_sample$p_value, method = "BH")
    tmp <- by_sample
    tmp$effect_source <- "threshold_sensitivity_delta"
    tmp$effect <- tmp$delta_high_minus_low
    summary <- o2b_p7_summarise_section_effects(
      tmp[tmp$status == "OK", , drop = FALSE],
      effect_col = "effect",
      group_cols = c("threshold_method", "feature")
    )
    if (nrow(summary)) {
      summary$primary_role <- "sensitivity_not_primary"
      summary$direction_consistency_label <- paste0(summary$positive_direction_count, "/", summary$n_sections, " positive")
    }
  } else {
    summary <- data.frame()
  }

  list(audit = audit, by_sample = by_sample, summary = summary)
}

o2b_p7_knn_index <- function(x, y, k = 6L, max_n = 6000L) {
  n <- length(x)
  if (n < 3L || n > max_n || any(!is.finite(x) | !is.finite(y))) return(NULL)
  coords <- cbind(as.numeric(x), as.numeric(y))
  dm <- as.matrix(stats::dist(coords))
  diag(dm) <- Inf
  t(apply(dm, 1L, function(z) order(z)[seq_len(min(k, length(z) - 1L))]))
}

o2b_p7_moran_i <- function(x, nbr) {
  x <- as.numeric(x)
  if (is.null(nbr) || length(x) != nrow(nbr)) return(NA_real_)
  ok <- is.finite(x)
  if (sum(ok) < 20L || stats::sd(x[ok]) <= 0) return(NA_real_)
  if (!all(ok)) return(NA_real_)
  z <- x - mean(x, na.rm = TRUE)
  n <- length(z)
  num <- 0
  w <- 0
  for (i in seq_len(nrow(nbr))) {
    js <- nbr[i, ]
    js <- js[is.finite(js)]
    if (!length(js)) next
    num <- num + sum(z[i] * z[js], na.rm = TRUE)
    w <- w + length(js)
  }
  den <- sum(z^2, na.rm = TRUE)
  if (!is.finite(den) || den <= 0 || w <= 0) return(NA_real_)
  (n / w) * (num / den)
}

o2b_p7_moran_permutation <- function(x, nbr, observed = NULL, B = 999L, seed = 92912L) {
  x <- as.numeric(x)
  if (is.null(observed)) observed <- o2b_p7_moran_i(x, nbr)
  if (!is.finite(observed) || is.null(nbr) || length(x) != nrow(nbr) ||
      sum(is.finite(x)) < 20L || as.integer(B) < 99L) {
    return(c(p_two_sided = NA_real_, null_mean = NA_real_, null_sd = NA_real_))
  }
  set.seed(as.integer(seed))
  null <- vapply(seq_len(as.integer(B)), function(i) {
    o2b_p7_moran_i(sample(x, length(x), replace = FALSE), nbr)
  }, numeric(1))
  null <- null[is.finite(null)]
  if (!length(null)) return(c(p_two_sided = NA_real_, null_mean = NA_real_, null_sd = NA_real_))
  centered_observed <- abs(observed - mean(null))
  centered_null <- abs(null - mean(null))
  c(
    p_two_sided = (1 + sum(centered_null >= centered_observed)) / (1 + length(null)),
    null_mean = mean(null),
    null_sd = stats::sd(null)
  )
}

o2b_p7_spatial_stats <- function(scores, distance_permutation_B = 999L,
                                moran_permutation_B = 999L) {
  features <- intersect(c("OLFML2B", "CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index", "OLFML2B_TME_composite_score"), names(scores))
  context_features <- intersect(c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Epithelial_Tumor", "Myeloid", "Immune_Exclusion_Index"), names(scores))
  auto_rows <- list()
  bivar_rows <- list()
  neigh_rows <- list()
  nn_rows <- list()
  niche_rows <- list()
  enrich_rows <- list()

  feature_high <- function(x, sparse_positive = TRUE) {
    x <- suppressWarnings(as.numeric(x))
    out <- rep(FALSE, length(x))
    ok <- is.finite(x)
    if (sum(ok) < 20L) return(out)
    q75 <- suppressWarnings(stats::quantile(x[ok], 0.75, na.rm = TRUE, names = FALSE))
    if (is.finite(q75) && q75 > 0) {
      out[ok & x >= q75] <- TRUE
    } else if (isTRUE(sparse_positive)) {
      out[ok & x > 0] <- TRUE
    }
    # Prevent non-informative high flags that include all or no spots.
    if (sum(out, na.rm = TRUE) <= 0L || sum(out, na.rm = TRUE) >= sum(ok)) out[] <- FALSE
    out
  }

  corrected_or <- function(exposure, outcome) {
    exposure <- exposure %in% TRUE
    outcome <- outcome %in% TRUE
    a <- sum(exposure & outcome, na.rm = TRUE)
    b <- sum(exposure & !outcome, na.rm = TRUE)
    c <- sum(!exposure & outcome, na.rm = TRUE)
    d <- sum(!exposure & !outcome, na.rm = TRUE)
    or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    pv <- tryCatch(stats::fisher.test(tab)$p.value, error = function(e) NA_real_)
    list(a = a, b = b, c = c, d = d, odds_ratio = or, log2_or = log2(or), p_value = pv)
  }

  nearest_dist <- function(from, to, x, y) {
    from <- from %in% TRUE
    to <- to %in% TRUE
    valid_xy <- is.finite(suppressWarnings(as.numeric(x))) & is.finite(suppressWarnings(as.numeric(y)))
    from <- from & valid_xy
    to <- to & valid_xy
    n_from <- sum(from, na.rm = TRUE)
    n_to <- sum(to, na.rm = TRUE)
    if (n_from < 1L || n_to < 1L) return(rep(NA_real_, n_from))
    a <- cbind(x[from], y[from])
    b <- cbind(x[to], y[to])
    apply(a, 1L, function(z) {
      dd <- sqrt((z[1] - b[, 1])^2 + (z[2] - b[, 2])^2)
      if (any(is.finite(dd))) min(dd, na.rm = TRUE) else NA_real_
    })
  }

  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    d$pxl_row <- suppressWarnings(as.numeric(d$pxl_row))
    d$pxl_col <- suppressWarnings(as.numeric(d$pxl_col))
    d <- d[is.finite(d$pxl_row) & is.finite(d$pxl_col), , drop = FALSE]
    if (nrow(d) < 50L) next

    nbr <- o2b_p7_knn_index(d$pxl_col, d$pxl_row, k = 6L)
    if (is.null(nbr)) next

    for (f in features) {
      observed_moran <- o2b_p7_moran_i(d[[f]], nbr)
      moran_perm <- o2b_p7_moran_permutation(
        d[[f]], nbr, observed = observed_moran,
        B = moran_permutation_B,
        seed = 92912L + match(sid, unique(scores$sample_id)) * 1000L + match(f, features)
      )
      auto_rows[[length(auto_rows) + 1L]] <- data.frame(
        sample_id = sid, feature = f, n_spots = nrow(d),
        morans_i = observed_moran,
        moran_permutation_p_two_sided = as.numeric(moran_perm["p_two_sided"]),
        moran_null_mean = as.numeric(moran_perm["null_mean"]),
        moran_null_sd = as.numeric(moran_perm["null_sd"]),
        moran_permutation_B = as.integer(moran_permutation_B),
        k_neighbors = 6L,
        status = ifelse(is.finite(observed_moran) && is.finite(moran_perm["p_two_sided"]),
                        "OK_PERMUTATION_CALIBRATED", "NOT_EVALUABLE"),
        inference_role = "continuous_primary_section_effect",
        stringsAsFactors = FALSE
      )
    }

    for (f in setdiff(features, "OLFML2B")) {
      neigh_mean <- rowMeans(matrix(d[[f]][as.vector(nbr)], nrow = nrow(nbr)), na.rm = TRUE)
      cc <- o2b_p7_spearman(d$OLFML2B, neigh_mean)
      bivar_rows[[length(bivar_rows) + 1L]] <- data.frame(
        sample_id = sid, target = "OLFML2B", neighbor_feature = f,
        n = as.integer(cc["n"]), rho_target_vs_neighbor_mean = as.numeric(cc["rho"]),
        p_value = as.numeric(cc["p_value"]), k_neighbors = 6L,
        status = ifelse(is.finite(cc["rho"]), "OK", "NOT_EVALUABLE"),
        inference_role = "continuous_primary_neighbor_colocalization",
        stringsAsFactors = FALSE
      )

      # Supportive threshold-only neighborhood comparison. The primary inference
      # remains continuous and section-level.
      grp2 <- o2b_p7_define_olfml2b_groups(d$OLFML2B, min_positive_spots = 30L)
      g <- grp2$group
      a <- neigh_mean[g == "OLFML2B_high"]
      b <- neigh_mean[g == "OLFML2B_low"]
      pv <- if (sum(is.finite(a)) >= 10 && sum(is.finite(b)) >= 10) {
        tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else NA_real_
      neigh_rows[[length(neigh_rows) + 1L]] <- data.frame(
        sample_id = sid, neighbor_feature = f,
        n_high = sum(is.finite(a)), n_low = sum(is.finite(b)),
        mean_neighbor_high = mean(a, na.rm = TRUE), mean_neighbor_low = mean(b, na.rm = TRUE),
        delta_high_minus_low = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
        p_value = pv, k_neighbors = 6L,
        status = ifelse(is.finite(pv), "OK", "NOT_EVALUABLE"),
        inference_role = "threshold_supportive_not_primary",
        stringsAsFactors = FALSE
      )
    }

    olfml2b_binary <- o2b_p7_define_binary_high(d$OLFML2B, min_positive_spots = 30L)
    d$OLFML2B_high <- olfml2b_binary$high %in% TRUE
    d$OLFML2B_high_method <- olfml2b_binary$method
    d$OLFML2B_high_cutoff <- olfml2b_binary$high_cutoff
    d$CAF_high <- feature_high(d$CAF_TGFb_axis)
    d$ECM_high <- feature_high(d$ECM_Remodeling)
    d$TGFb_high <- feature_high(d$TGFb_Response)
    d$CD8_high <- feature_high(d$CD8_Cytotoxic)
    d$EPI_high <- feature_high(d$Epithelial_Tumor)
    d$MYELOID_high <- feature_high(d$Myeloid)
    d$EXCLUSION_high <- feature_high(d$Immune_Exclusion_Index)

    # Independent context class for visualization only. It does NOT include
    # OLFML2B in the class definition, avoiding circular enrichment tests.
    d$niche <- "Other"
    d$niche[d$CD8_high %in% TRUE] <- "CD8_cytotoxic_high"
    d$niche[d$EPI_high %in% TRUE] <- "Epithelial_like"
    d$niche[(d$CAF_high %in% TRUE) | (d$ECM_high %in% TRUE)] <- "CAF_ECM_rich"
    d$niche[(d$EXCLUSION_high %in% TRUE)] <- "Immune_exclusion_high"

    niche_rows[[length(niche_rows) + 1L]] <- d[, intersect(c(
      "sample_id", "spot_id", "pxl_row", "pxl_col", "OLFML2B",
      "CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic",
      "Epithelial_Tumor", "Myeloid", "Immune_Exclusion_Index",
      "OLFML2B_high", "OLFML2B_high_method", "OLFML2B_high_cutoff",
      "CAF_high", "ECM_high", "TGFb_high", "CD8_high", "EPI_high", "MYELOID_high", "EXCLUSION_high", "niche"
    ), names(d)), drop = FALSE]

    flags <- list(
      CAF_TGFb_axis_high = d$CAF_high,
      ECM_Remodeling_high = d$ECM_high,
      TGFb_Response_high = d$TGFb_high,
      CD8_Cytotoxic_high = d$CD8_high,
      Epithelial_Tumor_high = d$EPI_high,
      Myeloid_high = d$MYELOID_high,
      Immune_Exclusion_Index_high = d$EXCLUSION_high
    )

    for (nm in names(flags)) {
      ft <- corrected_or(d$OLFML2B_high, flags[[nm]])
      enrich_rows[[length(enrich_rows) + 1L]] <- data.frame(
        sample_id = sid,
        context_feature = nm,
        niche = nm,
        n_spots = nrow(d),
        n_OLFML2B_high = sum(d$OLFML2B_high, na.rm = TRUE),
        n_context_high = sum(flags[[nm]], na.rm = TRUE),
        n_overlap = ft$a,
        odds_ratio = ft$odds_ratio,
        log2_odds_ratio = ft$log2_or,
        log2_odds_ratio_capped = max(min(ft$log2_or, 5), -5),
        p_value = ft$p_value,
        status = ifelse(is.finite(ft$odds_ratio) && sum(flags[[nm]], na.rm = TRUE) > 0, "OK", "NOT_EVALUABLE"),
        inference_role = "independent_context_overlap_sensitivity",
        correction = "Haldane_Anscombe_0.5_for_OR; Fisher_exact_for_P",
        stringsAsFactors = FALSE
      )
    }

    # Permutation-calibrated distance sensitivity. This is supportive only and
    # not used as the main Part7 inference.
    from <- d$OLFML2B_high %in% TRUE
    for (dest in c("CAF_high", "ECM_high", "TGFb_high", "CD8_high", "EXCLUSION_high")) {
      to <- d[[dest]] %in% TRUE
      vals <- nearest_dist(from, to, d$pxl_col, d$pxl_row)
      obs <- stats::median(vals, na.rm = TRUE)
      random_medians <- rep(NA_real_, as.integer(distance_permutation_B))
      if (sum(from, na.rm = TRUE) >= 10 && sum(to, na.rm = TRUE) >= 10 && sum(to, na.rm = TRUE) < nrow(d)) {
        set.seed(2025 + match(sid, unique(scores$sample_id)) + match(dest, c("CAF_high", "ECM_high", "TGFb_high", "CD8_high", "EXCLUSION_high")) * 100L)
        for (bb in seq_len(as.integer(distance_permutation_B))) {
          rnd <- rep(FALSE, nrow(d))
          rnd[sample.int(nrow(d), size = sum(from, na.rm = TRUE), replace = FALSE)] <- TRUE
          random_medians[bb] <- stats::median(nearest_dist(rnd, to, d$pxl_col, d$pxl_row), na.rm = TRUE)
        }
      }
      rand_med <- stats::median(random_medians, na.rm = TRUE)
      perm_p <- if (is.finite(obs) && any(is.finite(random_medians))) {
        (1 + sum(random_medians <= obs, na.rm = TRUE)) / (1 + sum(is.finite(random_medians)))
      } else NA_real_
      nn_rows[[length(nn_rows) + 1L]] <- data.frame(
        sample_id = sid, from_group = "OLFML2B_high", to_group = dest,
        n_from = sum(from, na.rm = TRUE), n_to = sum(to, na.rm = TRUE),
        median_distance_observed = obs,
        median_distance_random = rand_med,
        delta_observed_minus_random = obs - rand_med,
        permutation_p_closer_than_random = perm_p,
        permutation_B = as.integer(distance_permutation_B),
        status = ifelse(is.finite(obs) && is.finite(rand_med), "OK", "NOT_EVALUABLE"),
        inference_role = "distance_sensitivity_permutation_supportive_only",
        stringsAsFactors = FALSE
      )
    }
  }

  auto <- o2b_p7_bind_rows(auto_rows)
  bivar <- o2b_p7_bind_rows(bivar_rows)
  neigh <- o2b_p7_bind_rows(neigh_rows)
  niche <- o2b_p7_bind_rows(niche_rows)
  enrich <- o2b_p7_bind_rows(enrich_rows)
  nn <- o2b_p7_bind_rows(nn_rows)

  if (nrow(bivar)) bivar$fdr <- stats::p.adjust(bivar$p_value, method = "BH")
  if (nrow(auto)) auto$moran_fdr <- stats::p.adjust(auto$moran_permutation_p_two_sided, method = "BH")
  if (nrow(neigh)) neigh$fdr <- stats::p.adjust(neigh$p_value, method = "BH")
  if (nrow(enrich)) enrich$fdr <- stats::p.adjust(enrich$p_value, method = "BH")
  if (nrow(nn)) nn$fdr <- stats::p.adjust(nn$permutation_p_closer_than_random, method = "BH")

  list(auto = auto, bivar = bivar, neigh = neigh, niche = niche, enrich = enrich, nn = nn)
}

o2b_p7_select_representative_sample <- function(scores) {
  d <- aggregate(cbind(n_spots = scores$spot_id, detected = scores$OLFML2B > 0) ~ sample_id, scores, function(x) length(x))
  # aggregate above names are awkward; recompute robustly
  rows <- lapply(split(scores, scores$sample_id), function(z) {
    data.frame(sample_id = z$sample_id[1], n_spots = nrow(z), olfml2b_detected_fraction = mean(z$OLFML2B > 0, na.rm = TRUE), coord_fraction = mean(is.finite(as.numeric(z$pxl_row)) & is.finite(as.numeric(z$pxl_col))), stringsAsFactors = FALSE)
  })
  x <- o2b_p7_bind_rows(rows)
  x <- x[order(-x$coord_fraction, -x$olfml2b_detected_fraction, -x$n_spots), , drop = FALSE]
  x$sample_id[1]
}

o2b_p7_make_figures <- function(scores, corr, shift, spatial, dirs, continuous_summary = NULL, threshold_summary = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  pal <- o2b_p7_palette()
  fig_dir <- dirs$figures

  score <- scores
  score$pxl_row <- suppressWarnings(as.numeric(score$pxl_row))
  score$pxl_col <- suppressWarnings(as.numeric(score$pxl_col))
  score <- score[is.finite(score$pxl_row) & is.finite(score$pxl_col), , drop = FALSE]
  if (!nrow(score)) return(invisible(FALSE))

  sid <- o2b_p7_select_representative_sample(score)
  d <- score[score$sample_id == sid, , drop = FALSE]
  d$plot_y <- -d$pxl_row

  stat <- paste0("sample=", sid, "; spots=", nrow(d), "; OLFML2B+=", o2b_p7_fmt_num(mean(d$OLFML2B > 0, na.rm = TRUE) * 100, 1), "%")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = pxl_col, y = plot_y, color = OLFML2B)) +
    ggplot2::geom_point(size = 0.72, alpha = 0.95) +
    ggplot2::scale_color_gradient(low = "white", high = pal["red"], name = "OLFML2B") +
    ggplot2::coord_fixed() +
    o2b_p7_theme(11) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank()) +
    ggplot2::labs(title = "Spatial OLFML2B expression map", subtitle = stat)
  o2b_p7_save_plot(p, "FIG7A_spatial_OLFML2B_map_CNS", fig_dir, width = 6.2, height = 5.5)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = pxl_col, y = plot_y, color = CAF_TGFb_axis)) +
    ggplot2::geom_point(size = 0.72, alpha = 0.95) +
    ggplot2::scale_color_gradient(low = "white", high = pal["red"], name = "CAF/TGFb") +
    ggplot2::coord_fixed() +
    o2b_p7_theme(11) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank()) +
    ggplot2::labs(title = "Spatial CAF/TGFβ axis map", subtitle = paste0("same section as Figure 7A; spots=", nrow(d)))
  o2b_p7_save_plot(p, "FIG7B_spatial_CAF_TGFb_map_CNS", fig_dir, width = 6.2, height = 5.5)

  ub_grp <- o2b_p7_define_binary_high(d$OLFML2B, min_positive_spots = 30L)
  q2 <- stats::quantile(d$CAF_TGFb_axis, 0.75, na.rm = TRUE)
  d$OLFML2B_high <- ub_grp$high %in% TRUE
  d$CAF_high <- d$CAF_TGFb_axis >= q2
  d$bivariate_class <- "low / low"
  d$bivariate_class[d$OLFML2B_high & !d$CAF_high] <- "OLFML2B-high only"
  d$bivariate_class[!d$OLFML2B_high & d$CAF_high] <- "CAF-high only"
  d$bivariate_class[d$OLFML2B_high & d$CAF_high] <- "double-high"
  d$bivariate_class <- as.character(d$bivariate_class)
  d$bivariate_class[is.na(d$bivariate_class) | !nzchar(d$bivariate_class)] <- "low / low"
  biv_levels <- c("low / low", "OLFML2B-high only", "CAF-high only", "double-high")
  biv_cols <- c("#D1D5DB", "#2166AC", "#D95F02", "#B2182B")
  names(biv_cols) <- biv_levels
  d$bivariate_class <- factor(d$bivariate_class, levels = biv_levels)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = pxl_col, y = plot_y, color = bivariate_class)) +
    ggplot2::geom_point(size = 0.72, alpha = 0.95, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = biv_cols, breaks = biv_levels, drop = FALSE, na.translate = FALSE, name = NULL) +
    ggplot2::coord_fixed() +
    o2b_p7_theme(11) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(title = "Bivariate OLFML2B-high and CAF/TGFβ-high niche map", subtitle = paste0("OLFML2B-high=", ub_grp$method, "; double-high spots=", sum(d$bivariate_class == "double-high", na.rm = TRUE)))
  o2b_p7_save_plot(p, "FIG7C_spatial_bivariate_OLFML2B_CAF_map_CNS", fig_dir, width = 7.0, height = 5.5)

  if (nrow(corr)) {
    cc <- corr[corr$feature %in% c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index"), , drop = FALSE]
    cc$label <- ifelse(is.finite(cc$rho), sprintf("%.2f%s", cc$rho, ifelse(cc$fdr < 0.05, "*", "")), "")
    p <- ggplot2::ggplot(cc, ggplot2::aes(x = feature, y = sample_id, fill = rho)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
      ggplot2::scale_fill_gradient2(low = pal["blue"], mid = "white", high = pal["red"], midpoint = 0, name = "rho", na.value = "#F3F4F6") +
      o2b_p7_theme(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
      ggplot2::labs(title = "Spot-level OLFML2B spatial correlation summary", subtitle = "Cell values show Spearman rho; * indicates FDR<0.05", x = NULL, y = NULL)
    o2b_p7_save_plot(p, "FIG7D_spatial_correlation_heatmap_CNS", fig_dir, width = 8.3, height = max(4.2, 0.28 * length(unique(cc$sample_id)) + 2))
  }

  if (nrow(spatial$auto)) {
    aa <- spatial$auto[spatial$auto$feature %in% c("OLFML2B", "CAF_TGFb_axis", "ECM_Remodeling", "CD8_Cytotoxic", "Immune_Exclusion_Index"), , drop = FALSE]
    aa <- aa[is.finite(suppressWarnings(as.numeric(aa$morans_i))) & !is.na(aa$feature), , drop = FALSE]
    if (nrow(aa)) {
    p <- ggplot2::ggplot(aa, ggplot2::aes(x = morans_i, y = feature)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = pal["grey"]) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F3F4F6", color = pal["dark"], width = 0.55) +
      ggplot2::geom_jitter(height = 0.08, width = 0, size = 1.7, alpha = 0.75, color = pal["red"]) +
      o2b_p7_theme(10) +
      ggplot2::labs(title = "Spatial autocorrelation across Visium samples", subtitle = paste0("Moran's I; samples=", length(unique(aa$sample_id))), x = "Moran's I", y = NULL)
    o2b_p7_save_plot(p, "FIG7E_spatial_autocorrelation_summary_CNS", fig_dir, width = 7.0, height = 4.8)
    }
  }

  if (nrow(spatial$bivar)) {
    bb <- spatial$bivar[spatial$bivar$neighbor_feature %in% c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index"), , drop = FALSE]
    bb <- bb[is.finite(suppressWarnings(as.numeric(bb$rho_target_vs_neighbor_mean))) & !is.na(bb$neighbor_feature), , drop = FALSE]
    if (nrow(bb)) {
    bb$significance <- ifelse(is.finite(bb$fdr) & bb$fdr < 0.05, "FDR<0.05", "NS")
    p <- ggplot2::ggplot(bb, ggplot2::aes(x = rho_target_vs_neighbor_mean, y = neighbor_feature)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = pal["grey"]) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F3F4F6", color = pal["dark"], width = 0.55) +
      ggplot2::geom_jitter(ggplot2::aes(color = significance), height = 0.08, width = 0, size = 1.7, alpha = 0.8) +
      ggplot2::scale_color_manual(values = c("FDR<0.05" = "#B2182B", "NS" = "#9CA3AF"), name = NULL) +
      o2b_p7_theme(10) +
      ggplot2::labs(title = "OLFML2B neighbor-mean spatial co-localization", subtitle = "Spearman correlation with KNN neighbor-averaged feature scores; not bivariate Moran's I", x = "Spearman rho", y = NULL)
    o2b_p7_save_plot(p, "FIG7F_spatial_neighbor_mean_colocalization_CNS", fig_dir, width = 7.3, height = 4.8)
    }
  }

  if (nrow(spatial$neigh)) {
    nn <- spatial$neigh[spatial$neigh$neighbor_feature %in% c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index"), , drop = FALSE]
    nn <- nn[is.finite(suppressWarnings(as.numeric(nn$delta_high_minus_low))) & !is.na(nn$neighbor_feature), , drop = FALSE]
    if (nrow(nn)) {
    nn$significance <- ifelse(is.finite(nn$fdr) & nn$fdr < 0.05, "FDR<0.05", "NS")
    p <- ggplot2::ggplot(nn, ggplot2::aes(x = delta_high_minus_low, y = neighbor_feature)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = pal["grey"]) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F3F4F6", color = pal["dark"], width = 0.55) +
      ggplot2::geom_jitter(ggplot2::aes(color = significance), height = 0.08, width = 0, size = 1.7, alpha = 0.8) +
      ggplot2::scale_color_manual(values = c("FDR<0.05" = "#B2182B", "NS" = "#9CA3AF"), name = NULL) +
      o2b_p7_theme(10) +
      ggplot2::labs(title = "Neighborhood enrichment around OLFML2B-high spots", subtitle = "Difference in neighbor-averaged feature scores: OLFML2B-high minus OLFML2B-low", x = "Delta neighbor score", y = NULL)
    o2b_p7_save_plot(p, "FIG7G_spatial_neighborhood_enrichment_CNS", fig_dir, width = 7.5, height = 4.8)
    }
  }

  if (nrow(spatial$enrich)) {
    ee <- spatial$enrich
    value_col <- if ("log2_odds_ratio_capped" %in% names(ee)) "log2_odds_ratio_capped" else "log2_odds_ratio"
    feature_col <- if ("context_feature" %in% names(ee)) "context_feature" else "niche"
    ee[[value_col]] <- suppressWarnings(as.numeric(ee[[value_col]]))
    ee <- ee[is.finite(ee[[value_col]]) & !is.na(ee[[feature_col]]), , drop = FALSE]
    if (nrow(ee)) {
      ee$label <- ifelse(is.finite(ee$odds_ratio), sprintf("%.1f%s", pmin(pmax(log2(pmax(ee$odds_ratio, .Machine$double.xmin)), -5), 5), ifelse(is.finite(ee$fdr) & ee$fdr < 0.05, "*", "")), "")
      p <- ggplot2::ggplot(ee, ggplot2::aes(x = .data[[feature_col]], y = sample_id, fill = .data[[value_col]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.25) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 2.4) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-5, 5), na.value = "#F3F4F6", name = "log2(OR) cap") +
        o2b_p7_theme(9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(title = "OLFML2B-high overlap with independent spatial context features", subtitle = "Supportive sensitivity analysis; OR uses Haldane-Anscombe correction and is not the primary Part7 inference", x = NULL, y = NULL)
      o2b_p7_save_plot(p, "FIG7H_spatial_niche_enrichment_heatmap_CNS", fig_dir, width = 9.4, height = max(4.2, 0.28 * length(unique(ee$sample_id)) + 2))
    }
  }

  # Supplementary all-sample maps.
  all_map <- function(feature, stem, title, high = pal["red"]) {
    dd <- score
    dd$plot_y <- -dd$pxl_row
    dd <- dd[is.finite(dd$pxl_col) & is.finite(dd$plot_y) & is.finite(suppressWarnings(as.numeric(dd[[feature]]))), , drop = FALSE]
    if (!nrow(dd)) return(invisible(FALSE))
    p <- ggplot2::ggplot(dd, ggplot2::aes(x = pxl_col, y = plot_y, color = .data[[feature]])) +
      ggplot2::geom_point(size = 0.28, alpha = 0.9) +
      ggplot2::scale_color_gradient(low = "white", high = high, name = feature) +
      ggplot2::coord_fixed() +
      ggplot2::facet_wrap(~sample_id) +
      o2b_p7_theme(8) +
      ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), legend.position = "right") +
      ggplot2::labs(title = title, subtitle = paste0("all samples; total spots=", nrow(dd)))
    o2b_p7_save_plot(p, stem, fig_dir, width = 10.5, height = 7.5)
  }
  all_map("OLFML2B", "SFIG7A_all_sample_OLFML2B_spatial_maps_CNS", "All-sample OLFML2B spatial maps", pal["red"])
  all_map("CAF_TGFb_axis", "SFIG7B_all_sample_CAF_TGFb_spatial_maps_CNS", "All-sample CAF/TGFβ spatial maps", pal["red"])
  all_map("CD8_Cytotoxic", "SFIG7C_all_sample_CD8_Cytotoxic_spatial_maps_CNS", "All-sample CD8/cytotoxic spatial maps", pal["blue"])


  # Best-practice primary inference figure: section-level continuous effects.
  if (is.data.frame(continuous_summary) && nrow(continuous_summary)) {
    cs <- continuous_summary[
      continuous_summary$effect_source %in% c("spot_level_spearman", "neighbor_mean_colocalization", "spatial_autocorrelation_morans_i") &
        continuous_summary$feature %in% c("OLFML2B", "CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index"),
      ,
      drop = FALSE
    ]
    cs <- cs[is.finite(suppressWarnings(as.numeric(cs$median_effect))), , drop = FALSE]
    if (nrow(cs)) {
      cs$effect_source <- factor(
        cs$effect_source,
        levels = c("spot_level_spearman", "neighbor_mean_colocalization", "spatial_autocorrelation_morans_i"),
        labels = c("Spot-level continuous association", "Neighbor co-localization", "Spatial autocorrelation")
      )
      cs$feature <- factor(cs$feature, levels = rev(unique(cs$feature)))
      p <- ggplot2::ggplot(cs, ggplot2::aes(x = median_effect, y = feature)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = pal["grey"], linewidth = 0.35) +
        ggplot2::geom_segment(ggplot2::aes(x = q1_effect, xend = q3_effect, yend = feature), linewidth = 0.75, color = pal["dark"]) +
        ggplot2::geom_point(size = 2.3, color = pal["red"]) +
        ggplot2::facet_wrap(~effect_source, scales = "free_x") +
        o2b_p7_theme(10) +
        ggplot2::labs(
          title = "Section-level continuous spatial evidence",
          subtitle = "Median effect and IQR across patients after within-patient section aggregation; GC6 primary/metastasis are clustered",
          x = "Section-level effect size",
          y = NULL
        )
      o2b_p7_save_plot(p, "FIG7I_continuous_section_level_summary_CNS", fig_dir, width = 10.0, height = 5.2)
    }
  }

  # Threshold sensitivity figure: supportive only, not the primary evidence.
  if (is.data.frame(threshold_summary) && nrow(threshold_summary)) {
    ts <- threshold_summary[
      threshold_summary$feature %in% c("CAF_TGFb_axis", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic", "Immune_Exclusion_Index"),
      ,
      drop = FALSE
    ]
    ts <- ts[is.finite(suppressWarnings(as.numeric(ts$median_effect))), , drop = FALSE]
    if (nrow(ts)) {
      ts$threshold_method <- factor(
        ts$threshold_method,
        levels = c("detected_vs_undetected", "positive_only_median", "positive_only_upper_quartile", "adaptive_sparse_aware")
      )
      ts$label <- paste0(ts$positive_direction_count, "/", ts$n_patients)
      p <- ggplot2::ggplot(ts, ggplot2::aes(x = threshold_method, y = feature, fill = median_effect)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.25) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8, color = "#111827") +
        ggplot2::scale_fill_gradient2(low = pal["blue"], mid = "white", high = pal["red"], midpoint = 0, name = "Median delta", na.value = "#F3F4F6") +
        o2b_p7_theme(9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(
          title = "Threshold sensitivity, supportive only",
          subtitle = "Cell text shows positive-direction patients / evaluable patients; thresholds remain sensitivity-only",
          x = NULL,
          y = NULL
        )
      o2b_p7_save_plot(p, "FIG7J_threshold_sensitivity_summary_CNS", fig_dir, width = 9.0, height = 5.0)
    }
  }

  invisible(TRUE)
}

o2b_p7_loader_diagnostic_figure <- function(qc, dirs) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(qc)) return(invisible(FALSE))
  ss <- as.data.frame(table(qc$parser_status), stringsAsFactors = FALSE)
  names(ss) <- c("parser_status", "n_samples")
  pal <- o2b_p7_palette()
  p <- ggplot2::ggplot(ss, ggplot2::aes(x = stats::reorder(parser_status, n_samples), y = n_samples)) +
    ggplot2::geom_col(fill = pal["grey"], width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = n_samples), hjust = -0.1, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    o2b_p7_theme(11) +
    ggplot2::labs(title = "Part7 spatial loader diagnostic", subtitle = "Biological spatial maps require a non-empty score table", x = NULL, y = "Samples")
  o2b_p7_save_plot(p, "FIG7Z_loader_diagnostic_CNS", dirs$figures, width = 8.5, height = 4.8)
  invisible(TRUE)
}


# ============================================================================
# OLFML2B-specific inferential overrides and competing-source modules
# ============================================================================

o2b_p7_exact_signflip_p <- function(x, statistic = c("median", "mean"), B = 100000L, seed = 20260722L) {
  statistic <- match.arg(statistic)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) return(list(p_value = NA_real_, method = "NOT_EVALUABLE_N_LT_3", n = n))
  stat_fun <- if (statistic == "median") stats::median else mean
  obs <- abs(stat_fun(x, na.rm = TRUE))
  if (!is.finite(obs)) return(list(p_value = NA_real_, method = "NOT_EVALUABLE", n = n))
  if (n <= 15L) {
    ids <- 0:(2^n - 1L)
    null <- vapply(ids, function(id) {
      bits <- as.integer(intToBits(id))[seq_len(n)]
      s <- ifelse(bits == 1L, 1, -1)
      abs(stat_fun(x * s, na.rm = TRUE))
    }, numeric(1))
    p <- sum(null >= obs - 1e-12) / length(null)
    return(list(p_value = p, method = paste0("EXACT_SIGN_FLIP_2^", n), n = n))
  }
  set.seed(as.integer(seed))
  null <- vapply(seq_len(as.integer(B)), function(i) {
    s <- sample(c(-1, 1), n, replace = TRUE)
    abs(stat_fun(x * s, na.rm = TRUE))
  }, numeric(1))
  list(p_value = (1 + sum(null >= obs - 1e-12)) / (1 + length(null)),
       method = paste0("MONTE_CARLO_SIGN_FLIP_B", as.integer(B)), n = n)
}

o2b_p7_bootstrap_ci <- function(x, statistic = c("median", "mean"), B = 2000L, seed = 20260722L) {
  statistic <- match.arg(statistic)
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) < 3L) return(c(low = NA_real_, high = NA_real_))
  stat_fun <- if (statistic == "median") stats::median else mean
  set.seed(as.integer(seed))
  vals <- replicate(as.integer(B), stat_fun(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  q <- stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  stats::setNames(q, c("low", "high"))
}

o2b_p7_summarise_section_effects <- function(effects, effect_col = "effect",
                                             group_cols = c("effect_source", "feature"),
                                             bootstrap_B = 2000L) {
  if (!is.data.frame(effects) || !nrow(effects)) return(data.frame())
  effects <- o2b_p7_attach_official_sample_metadata(effects)
  effects[[effect_col]] <- suppressWarnings(as.numeric(effects[[effect_col]]))
  effects <- effects[is.finite(effects[[effect_col]]), , drop = FALSE]
  if (!nrow(effects)) return(data.frame())
  key <- do.call(paste, c(effects[, group_cols, drop = FALSE], sep = "\r"))
  rows <- lapply(split(seq_len(nrow(effects)), key), function(ii) {
    d <- effects[ii, , drop = FALSE]
    pid <- as.character(d$patient_id)
    pid[is.na(pid) | !nzchar(pid)] <- as.character(d$sample_id[is.na(pid) | !nzchar(pid)])
    patient <- stats::aggregate(d[[effect_col]], by = list(patient_id = pid),
                                FUN = function(z) stats::median(z, na.rm = TRUE))
    x <- suppressWarnings(as.numeric(patient$x)); x <- x[is.finite(x)]
    loo <- o2b_p7_loo_direction_stability(x)
    sf <- o2b_p7_exact_signflip_p(x, statistic = "median")
    seed_key <- paste(vapply(group_cols, function(g) as.character(d[[g]][1]), character(1)), collapse = "|")
    ci <- o2b_p7_bootstrap_ci(x, statistic = "median", B = bootstrap_B,
                              seed = 20260722L + sum(utf8ToInt(seed_key)))
    row <- as.list(d[1, group_cols, drop = FALSE])
    cbind(as.data.frame(row, stringsAsFactors = FALSE), data.frame(
      n_sections = length(unique(d$sample_id)), n_patients = length(x),
      median_effect = if (length(x)) stats::median(x) else NA_real_,
      bootstrap_ci_low = unname(ci["low"]), bootstrap_ci_high = unname(ci["high"]),
      q1_effect = if (length(x)) as.numeric(stats::quantile(x, 0.25, names = FALSE)) else NA_real_,
      q3_effect = if (length(x)) as.numeric(stats::quantile(x, 0.75, names = FALSE)) else NA_real_,
      positive_direction_count = sum(x > 0), negative_direction_count = sum(x < 0),
      zero_direction_count = sum(x == 0),
      positive_direction_fraction = if (length(x)) mean(x > 0) else NA_real_,
      exact_signflip_p = sf$p_value, exact_signflip_method = sf$method,
      wilcoxon_p_against_zero = o2b_p7_wilcox_zero_p(x),
      loo_min_median = loo$loo_min_median, loo_max_median = loo$loo_max_median,
      loo_direction_stability = loo$loo_direction_stability,
      inference_unit = "patient_after_within_patient_section_median",
      small_n_exploratory = length(x) < 6L,
      formal_significance_eligible = length(x) >= 6L,
      interpretation_boundary = "Spot effects are summarized within section and then collapsed within patient; spot-level P values are not biological inference.",
      stringsAsFactors = FALSE
    ))
  })
  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) {
    fam <- if ("effect_source" %in% names(out)) out$effect_source else rep("all", nrow(out))
    out$exact_signflip_fdr <- stats::ave(out$exact_signflip_p, fam,
                                        FUN = function(v) stats::p.adjust(v, method = "BH"))
    out$wilcoxon_fdr <- stats::ave(out$wilcoxon_p_against_zero, fam,
                                   FUN = function(v) stats::p.adjust(v, method = "BH"))
  }
  out
}

o2b_p7_corr_tables <- function(scores) {
  features <- intersect(o2b_p7_all_context_features(), names(scores))
  rows <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    for (f in features) {
      cc <- o2b_p7_spearman(d$OLFML2B, d[[f]])
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sid, feature = f, n = as.integer(cc["n"]),
        rho = as.numeric(cc["rho"]), p_value = as.numeric(cc["p_value"]),
        status = ifelse(is.finite(cc["rho"]), "OK", "NOT_EVALUABLE"),
        analysis_role = "section_level_spot_continuous_effect_only",
        normalization = "log1p_CPM_10000_full_matrix", stringsAsFactors = FALSE)
    }
  }
  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) out$fdr_descriptive <- stats::ave(out$p_value, out$feature,
                                                    FUN = function(v) stats::p.adjust(v, method = "BH"))
  out
}

o2b_p7_high_low_tables <- function(scores) {
  features <- intersect(o2b_p7_all_context_features(), names(scores))
  rows <- list(); thresholds <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    grp <- o2b_p7_define_olfml2b_groups(d$OLFML2B, min_positive_spots = 20L)
    d$spot_group <- grp$group
    thresholds[[length(thresholds)+1L]] <- data.frame(
      sample_id=sid, grouping_method=grp$method, n_total=grp$n_total,
      n_positive=grp$n_positive, detected_fraction=grp$detected_fraction,
      high_cutoff=grp$high_cutoff, low_cutoff=grp$low_cutoff,
      n_high=sum(d$spot_group=="OLFML2B_high",na.rm=TRUE),
      n_low=sum(d$spot_group=="OLFML2B_low",na.rm=TRUE), stringsAsFactors=FALSE)
    for (f in features) {
      a <- d[[f]][d$spot_group=="OLFML2B_high"]
      b <- d[[f]][d$spot_group=="OLFML2B_low"]
      pv <- if (sum(is.finite(a))>=10L && sum(is.finite(b))>=10L)
        tryCatch(stats::wilcox.test(a,b,exact=FALSE)$p.value,error=function(e) NA_real_) else NA_real_
      rows[[length(rows)+1L]] <- data.frame(
        sample_id=sid, grouping=grp$method, feature=f,
        n_high=sum(is.finite(a)), n_low=sum(is.finite(b)),
        mean_high=mean(a,na.rm=TRUE), mean_low=mean(b,na.rm=TRUE),
        delta_high_minus_low=mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE),
        p_value=pv, status=ifelse(is.finite(pv),"OK","NOT_EVALUABLE"),
        inference_role="threshold_sensitivity_only", stringsAsFactors=FALSE)
    }
  }
  out <- o2b_p7_bind_rows(rows)
  if (nrow(out)) out$fdr_descriptive <- stats::ave(out$p_value,out$feature,
                                                    FUN=function(v) stats::p.adjust(v,method="BH"))
  attr(out,"thresholds") <- o2b_p7_bind_rows(thresholds)
  out
}

o2b_p7_threshold_sensitivity <- function(scores, min_group_spots = 10L) {
  features <- intersect(o2b_p7_all_context_features(), names(scores))
  methods <- c("detected_vs_undetected", "positive_only_median",
               "positive_only_upper_quartile", "adaptive_sparse_aware")
  rows <- list(); audit <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id==sid,,drop=FALSE]
    for (m in methods) {
      grp <- o2b_p7_threshold_groups(d$OLFML2B, m, min_positive_spots=20L)
      audit[[length(audit)+1L]] <- data.frame(sample_id=sid,threshold_method=m,
        n_total=grp$n_total,n_positive=grp$n_positive,detected_fraction=grp$detected_fraction,
        high_cutoff=grp$high_cutoff,low_cutoff=grp$low_cutoff,n_high=grp$n_high,n_low=grp$n_low,
        stringsAsFactors=FALSE)
      for (f in features) {
        a <- d[[f]][grp$group=="high"]; b <- d[[f]][grp$group=="low"]
        ok <- sum(is.finite(a))>=min_group_spots && sum(is.finite(b))>=min_group_spots
        pv <- if (ok) tryCatch(stats::wilcox.test(a,b,exact=FALSE)$p.value,error=function(e) NA_real_) else NA_real_
        rows[[length(rows)+1L]] <- data.frame(sample_id=sid,threshold_method=m,feature=f,
          n_high=sum(is.finite(a)),n_low=sum(is.finite(b)),mean_high=mean(a,na.rm=TRUE),mean_low=mean(b,na.rm=TRUE),
          delta_high_minus_low=mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE),p_value=pv,
          status=ifelse(ok && is.finite(pv),"OK","NOT_EVALUABLE"),
          inference_role="threshold_sensitivity_only",stringsAsFactors=FALSE)
      }
    }
  }
  by_sample <- o2b_p7_bind_rows(rows); au <- o2b_p7_bind_rows(audit)
  if (nrow(by_sample)) by_sample$fdr_descriptive <- stats::p.adjust(by_sample$p_value,method="BH")
  tmp <- by_sample[by_sample$status=="OK",,drop=FALSE]
  summary <- data.frame()
  if (nrow(tmp)) {
    tmp$effect_source <- "threshold_sensitivity_delta"; tmp$effect <- tmp$delta_high_minus_low
    summary <- o2b_p7_summarise_section_effects(tmp,"effect",c("threshold_method","feature"))
    if (nrow(summary)) summary$primary_role <- "sensitivity_not_primary"
  }
  list(audit=au,by_sample=by_sample,summary=summary)
}

o2b_p7_residualize <- function(y, total_counts, n_features) {
  y <- suppressWarnings(as.numeric(y)); tc <- log1p(suppressWarnings(as.numeric(total_counts)))
  nf <- log1p(suppressWarnings(as.numeric(n_features)))
  ok <- is.finite(y) & is.finite(tc) & is.finite(nf)
  out <- rep(NA_real_,length(y))
  if (sum(ok)<20L || stats::sd(y[ok])<=0) return(out)
  fit <- tryCatch(stats::lm(y[ok] ~ tc[ok] + nf[ok]),error=function(e) NULL)
  if (!is.null(fit)) out[ok] <- stats::residuals(fit)
  out
}

o2b_p7_technical_residual_sensitivity <- function(scores) {
  features <- intersect(o2b_p7_all_context_features(),names(scores)); rows <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id==sid,,drop=FALSE]
    rt <- o2b_p7_residualize(d$OLFML2B,d$total_counts,d$n_features)
    for (f in features) {
      rf <- o2b_p7_residualize(d[[f]],d$total_counts,d$n_features)
      cc <- o2b_p7_spearman(rt,rf)
      rows[[length(rows)+1L]] <- data.frame(sample_id=sid,feature=f,n=as.integer(cc["n"]),
        residual_rho=as.numeric(cc["rho"]),p_value=as.numeric(cc["p_value"]),
        status=ifelse(is.finite(cc["rho"]),"OK","NOT_EVALUABLE"),
        adjustment="within_section_residuals_from_log1p_total_counts_and_log1p_n_features",
        stringsAsFactors=FALSE)
    }
  }
  section <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(rows))
  tmp <- section; tmp$effect_source <- "technical_residual_spearman"; tmp$effect <- tmp$residual_rho
  patient <- o2b_p7_summarise_section_effects(tmp[tmp$status=="OK",,drop=FALSE],"effect",c("effect_source","feature"))
  list(section=section,patient=patient)
}

o2b_p7_ridge_fit <- function(y, x, seed = 20260722L) {
  y <- suppressWarnings(as.numeric(y)); x <- as.matrix(x)
  ok <- is.finite(y) & apply(x,1,function(z) all(is.finite(z)))
  y <- y[ok]; x <- x[ok,,drop=FALSE]
  if (length(y)<100L || stats::sd(y)<=0 || ncol(x)<2L) return(NULL)
  xs <- scale(x); ys <- as.numeric(scale(y))
  keep <- apply(xs,2,function(z) all(is.finite(z)) && stats::sd(z)>0)
  xs <- xs[,keep,drop=FALSE]
  if (ncol(xs)<2L) return(NULL)
  lambdas <- c(0.001,0.01,0.1,1,10,100)
  set.seed(as.integer(seed)); folds <- sample(rep(seq_len(5L),length.out=nrow(xs)))
  mse <- vapply(lambdas,function(lam){
    mean(vapply(seq_len(5L),function(k){
      tr <- folds!=k; te <- folds==k
      xtx <- crossprod(xs[tr,,drop=FALSE]) + diag(lam,ncol(xs))
      beta <- tryCatch(solve(xtx,crossprod(xs[tr,,drop=FALSE],ys[tr])),error=function(e) rep(NA_real_,ncol(xs)))
      mean((ys[te]-as.numeric(xs[te,,drop=FALSE]%*%beta))^2,na.rm=TRUE)
    },numeric(1)),na.rm=TRUE)
  },numeric(1))
  lam <- lambdas[which.min(mse)]
  beta <- solve(crossprod(xs)+diag(lam,ncol(xs)),crossprod(xs,ys))
  list(beta=as.numeric(beta),names=colnames(xs),lambda=lam,cv_mse=min(mse),n=nrow(xs))
}

o2b_p7_source_competition <- function(scores) {
  source_features <- intersect(o2b_p7_source_features(),names(scores)); rows <- list()
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id==sid,,drop=FALSE]
    x <- d[,source_features,drop=FALSE]
    x$log_total_counts <- log1p(d$total_counts); x$log_n_features <- log1p(d$n_features)
    fit <- o2b_p7_ridge_fit(d$OLFML2B,x,seed=20260722L+match(sid,unique(scores$sample_id)))
    if (is.null(fit)) next
    for (j in seq_along(fit$names)) rows[[length(rows)+1L]] <- data.frame(
      sample_id=sid,predictor=fit$names[j],standardized_ridge_coefficient=fit$beta[j],
      lambda=fit$lambda,cv_mse=fit$cv_mse,n_spots=fit$n,
      inference_role="source_disambiguation_sensitivity_not_cell_deconvolution",stringsAsFactors=FALSE)
  }
  section <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(rows))
  biological <- section[section$predictor %in% source_features,,drop=FALSE]
  tmp <- biological; tmp$effect_source <- "ridge_competing_source"; tmp$feature <- tmp$predictor; tmp$effect <- tmp$standardized_ridge_coefficient
  patient <- o2b_p7_summarise_section_effects(tmp,"effect",c("effect_source","feature"))
  list(section=section,patient=patient)
}

o2b_p7_classifier_overlap_audit <- function() {
  ms <- o2b_p7_marker_sets(); sources <- o2b_p7_source_features()
  states <- setdiff(o2b_p7_state_features(),"Immune_Exclusion_Index")
  rows <- list()
  for (s in sources) for (p in states) {
    ov <- intersect(ms[[s]],ms[[p]])
    rows[[length(rows)+1L]] <- data.frame(source_classifier=s,state_program=p,
      n_source_genes=length(ms[[s]]),n_program_genes=length(ms[[p]]),n_overlap=length(ov),
      overlap_genes=paste(ov,collapse=";"),
      interpretation=ifelse(length(ov)>0,"IDENTITY_COUPLED_OR_PARTLY_COUPLED","RELATIVELY_ORTHOGONAL"),
      stringsAsFactors=FALSE)
  }
  o2b_p7_bind_rows(rows)
}


o2b_p7_pretty_feature_labels <- function(x) {
  dict <- c(
    Fibroblast = "Fibroblast",
    Myofibroblast = "Myofibroblast",
    Smooth_Muscle = "Smooth muscle",
    Pericyte = "Pericyte",
    Myeloid_Macrophage = "Myeloid/macrophage",
    Endothelial = "Endothelial",
    Epithelial = "Epithelial",
    CAF_ECM = "CAF/ECM",
    ECM_Remodeling = "ECM remodeling",
    TGFb_Response = "TGF-β response",
    CD8_Cytotoxic = "CD8/cytotoxic"
  )
  x <- as.character(x)
  out <- unname(dict[x])
  out[is.na(out)] <- gsub("_", " ", x[is.na(out)], fixed = TRUE)
  out
}

o2b_p7_distance_feature_colors <- function() {
  c(
    Fibroblast = "#B2182B",
    Myofibroblast = "#D6604D",
    Pericyte = "#7B3294",
    Smooth_Muscle = "#C2A5CF",
    Myeloid_Macrophage = "#1B9E77",
    CD8_Cytotoxic = "#2166AC"
  )
}

o2b_p7_distance_stats <- function(scores, distance_permutation_B = 999L) {
  rows <- list()
  distance_features <- intersect(
    c("Fibroblast", "Myofibroblast", "Pericyte", "Smooth_Muscle", "Myeloid_Macrophage", "CD8_Cytotoxic"),
    names(scores)
  )
  high_flag <- function(x) {
    x <- suppressWarnings(as.numeric(x)); ok <- is.finite(x); out <- rep(FALSE, length(x))
    if (sum(ok) < 20L) return(out)
    q <- stats::quantile(x[ok], 0.75, names = FALSE, na.rm = TRUE)
    if (is.finite(q) && q > 0) out[ok & x >= q] <- TRUE else out[ok & x > 0] <- TRUE
    if (sum(out, na.rm = TRUE) == 0L || sum(out, na.rm = TRUE) >= sum(ok, na.rm = TRUE)) out[] <- FALSE
    out
  }
  qcut <- function(x) {
    x <- suppressWarnings(as.numeric(x)); r <- rank(x, ties.method = "average", na.last = "keep")
    out <- rep(NA_integer_, length(x)); ok <- is.finite(r)
    if (sum(ok)) out[ok] <- pmin(4L, pmax(1L, ceiling(4 * r[ok] / sum(ok))))
    out
  }
  nearest_other_for_all <- function(context_flag, x, y) {
    context_flag <- context_flag %in% TRUE
    x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
    valid <- is.finite(x) & is.finite(y)
    context_idx <- which(context_flag & valid)
    out <- rep(NA_real_, length(x))
    if (!length(context_idx)) return(out)
    for (ii in which(valid)) {
      cand <- context_idx[context_idx != ii]
      if (!length(cand)) next
      dd <- sqrt((x[ii] - x[cand])^2 + (y[ii] - y[cand])^2)
      if (any(is.finite(dd))) out[ii] <- min(dd[is.finite(dd)])
    }
    out
  }
  median_knn_spacing <- function(nbr, x, y) {
    if (is.null(nbr) || !length(nbr)) return(NA_real_)
    x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
    rr <- row(nbr); cc <- nbr
    dd <- sqrt((x[rr] - x[cc])^2 + (y[rr] - y[cc])^2)
    dd <- dd[is.finite(dd) & dd > 0]
    if (length(dd)) stats::median(dd, na.rm = TRUE) else NA_real_
  }
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    d$pxl_row <- suppressWarnings(as.numeric(d$pxl_row)); d$pxl_col <- suppressWarnings(as.numeric(d$pxl_col)
    )
    d <- d[is.finite(d$pxl_row) & is.finite(d$pxl_col), , drop = FALSE]
    if (nrow(d) < 50L) next
    nbr <- o2b_p7_knn_index(d$pxl_col, d$pxl_row, k = 6L)
    spacing_raw <- median_knn_spacing(nbr, d$pxl_col, d$pxl_row)
    target_group <- o2b_p7_define_olfml2b_groups(d$OLFML2B, min_positive_spots = 20L)
    from <- target_group$group %in% "OLFML2B_high"; from[is.na(from)] <- FALSE
    n_from <- sum(from, na.rm = TRUE)
    strata <- interaction(qcut(log1p(d$total_counts)), qcut(log1p(d$n_features)), drop = TRUE, lex.order = TRUE)
    matched_indices <- function(target_flag) {
      idx <- integer()
      for (st in unique(strata[!is.na(strata)])) {
        need <- sum(target_flag & strata == st, na.rm = TRUE)
        pool <- which(strata == st & !target_flag)
        if (need > 0L && length(pool) > 0L) idx <- c(idx, sample(pool, size = min(need, length(pool)), replace = FALSE))
      }
      idx
    }
    for (f in distance_features) {
      to <- high_flag(d[[f]]); n_to <- sum(to, na.rm = TRUE)
      nearest_raw <- nearest_other_for_all(to, d$pxl_col, d$pxl_row)
      obs_values_raw <- nearest_raw[from & is.finite(nearest_raw)]
      obs_raw <- if (length(obs_values_raw)) stats::median(obs_values_raw, na.rm = TRUE) else NA_real_
      obs_spot <- if (is.finite(obs_raw) && is.finite(spacing_raw) && spacing_raw > 0) obs_raw / spacing_raw else NA_real_
      n_overlap <- sum(from & to, na.rm = TRUE)
      overlap_fraction <- if (n_from > 0L) n_overlap / n_from else NA_real_
      null_raw <- rep(NA_real_, as.integer(distance_permutation_B))
      achieved <- integer(as.integer(distance_permutation_B))
      if (n_from >= 10L && n_to >= 10L && n_to < nrow(d)) {
        set.seed(20260722L + match(sid, unique(scores$sample_id)) * 100L + match(f, distance_features))
        for (b in seq_len(as.integer(distance_permutation_B))) {
          idx <- matched_indices(from); achieved[b] <- length(idx)
          vals <- nearest_raw[idx]; vals <- vals[is.finite(vals)]
          if (length(idx) >= max(5L, ceiling(0.8 * n_from)) && length(vals) >= max(5L, ceiling(0.8 * n_from))) {
            null_raw[b] <- stats::median(vals, na.rm = TRUE)
          }
        }
      }
      null_raw_med <- if (any(is.finite(null_raw))) stats::median(null_raw[is.finite(null_raw)], na.rm = TRUE) else NA_real_
      null_spot <- if (is.finite(null_raw_med) && is.finite(spacing_raw) && spacing_raw > 0) null_raw_med / spacing_raw else NA_real_
      delta_spot <- obs_spot - null_spot
      p <- if (is.finite(obs_spot) && any(is.finite(null_raw))) {
        null_spot_values <- null_raw[is.finite(null_raw)] / spacing_raw
        (1 + sum(null_spot_values <= obs_spot, na.rm = TRUE)) / (1 + length(null_spot_values))
      } else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sid, to_feature = f, n_target_high = n_from, n_context_high = n_to,
        n_target_context_overlap = n_overlap, target_context_overlap_fraction = overlap_fraction,
        n_target_with_finite_nearest_other = length(obs_values_raw),
        median_knn6_spacing_raw = spacing_raw,
        observed_median_nearest_other_distance_raw = obs_raw,
        null_median_nearest_other_distance_raw = null_raw_med,
        delta_nearest_other_distance_raw = obs_raw - null_raw_med,
        observed_median_distance_spot_units = obs_spot,
        null_median_distance_spot_units = null_spot,
        delta_observed_minus_null_spot_units = delta_spot,
        observed_median_distance = obs_spot,
        null_median_distance = null_spot,
        delta_observed_minus_null = delta_spot,
        empirical_p_closer_than_null = p,
        permutation_B = as.integer(distance_permutation_B),
        median_matched_control_spots = if (length(achieved)) stats::median(achieved, na.rm = TRUE) else NA_real_,
        target_grouping_method = target_group$method,
        target_detected_fraction = target_group$detected_fraction,
        target_positive_spots = target_group$n_positive,
        matching = "within-section quartiles of log1p total counts and log1p detected features",
        distance_definition = "nearest OTHER context-high spot; the query spot itself is excluded",
        distance_unit = "section-specific median KNN6 spot spacing",
        overlap_definition = "same-spot target-high and context-high overlap reported separately from distance",
        status = ifelse(is.finite(obs_spot) && is.finite(null_spot), "OK", "NOT_EVALUABLE"),
        inference_role = "threshold_distance_sensitivity_self_excluded_normalized_supportive_only",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(rows))
  if (nrow(out)) out$fdr_descriptive <- stats::p.adjust(out$empirical_p_closer_than_null, method = "BH")
  out
}

o2b_p7_spatial_stats <- function(scores, distance_permutation_B = 999L, moran_permutation_B = 999L) {
  auto_rows <- list(); neighbor_rows <- list(); graph_rows <- list()
  auto_features <- intersect(c("OLFML2B", "Fibroblast", "Myofibroblast", "Pericyte", "Smooth_Muscle", "Myeloid_Macrophage", "Epithelial", "CAF_ECM", "ECM_Remodeling", "TGFb_Response", "CD8_Cytotoxic"), names(scores))
  neighbor_features <- intersect(o2b_p7_all_context_features(), names(scores))
  for (sid in unique(scores$sample_id)) {
    d <- scores[scores$sample_id == sid, , drop = FALSE]
    d$pxl_row <- suppressWarnings(as.numeric(d$pxl_row)); d$pxl_col <- suppressWarnings(as.numeric(d$pxl_col))
    d <- d[is.finite(d$pxl_row) & is.finite(d$pxl_col), , drop = FALSE]
    if (nrow(d) < 50L) next
    nbr <- o2b_p7_knn_index(d$pxl_col, d$pxl_row, k = 6L)
    graph_rows[[length(graph_rows) + 1L]] <- data.frame(
      sample_id = sid, n_spots = nrow(d), k_neighbors = 6L,
      graph_status = ifelse(is.null(nbr), "NOT_EVALUABLE", "KNN6_WITHIN_SECTION"), stringsAsFactors = FALSE
    )
    if (is.null(nbr)) next
    for (f in auto_features) {
      obs <- o2b_p7_moran_i(d[[f]], nbr)
      pm <- o2b_p7_moran_permutation(d[[f]], nbr, obs, moran_permutation_B,
        seed = 92912L + match(sid, unique(scores$sample_id)) * 1000L + match(f, auto_features))
      auto_rows[[length(auto_rows) + 1L]] <- data.frame(
        sample_id = sid, feature = f, n_spots = nrow(d), morans_i = obs,
        permutation_p = as.numeric(pm["p_two_sided"]), null_mean = as.numeric(pm["null_mean"]),
        null_sd = as.numeric(pm["null_sd"]), permutation_B = as.integer(moran_permutation_B),
        status = ifelse(is.finite(obs), "OK", "NOT_EVALUABLE"), stringsAsFactors = FALSE
      )
    }
    for (f in neighbor_features) {
      nm <- rowMeans(matrix(d[[f]][as.vector(nbr)], nrow = nrow(nbr)), na.rm = TRUE)
      cc <- o2b_p7_spearman(d$OLFML2B, nm)
      neighbor_rows[[length(neighbor_rows) + 1L]] <- data.frame(
        sample_id = sid, neighbor_feature = f, n = as.integer(cc["n"]),
        rho_target_vs_neighbor_mean = as.numeric(cc["rho"]), p_value = as.numeric(cc["p_value"]),
        k_neighbors = 6L, status = ifelse(is.finite(cc["rho"]), "OK", "NOT_EVALUABLE"),
        statistic_name = "Spearman correlation with KNN neighbor-mean feature score",
        is_bivariate_morans_i = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  auto <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(auto_rows))
  neighbor <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(neighbor_rows))
  graph <- o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(graph_rows))
  distance <- o2b_p7_distance_stats(scores, distance_permutation_B = distance_permutation_B)
  if (nrow(auto)) auto$moran_fdr_descriptive <- stats::p.adjust(auto$permutation_p, method = "BH")
  if (nrow(neighbor)) neighbor$fdr_descriptive <- stats::p.adjust(neighbor$p_value, method = "BH")
  list(auto = auto, bivar = neighbor, nn = distance, graph = graph)
}


o2b_p7_continuous_sample_level_effects <- function(corr,spatial){
  rows<-list()
  if(is.data.frame(corr)&&nrow(corr)) rows[[1]]<-data.frame(sample_id=corr$sample_id,feature=corr$feature,
    effect_source="same_spot_spearman",effect=corr$rho,n_observations=corr$n,stringsAsFactors=FALSE)
  if(is.list(spatial)&&is.data.frame(spatial$bivar)&&nrow(spatial$bivar)) rows[[2]]<-data.frame(sample_id=spatial$bivar$sample_id,
    feature=spatial$bivar$neighbor_feature,effect_source="neighbor_mean_spearman",effect=spatial$bivar$rho_target_vs_neighbor_mean,
    n_observations=spatial$bivar$n,stringsAsFactors=FALSE)
  if(is.list(spatial)&&is.data.frame(spatial$auto)&&nrow(spatial$auto)) rows[[3]]<-data.frame(sample_id=spatial$auto$sample_id,
    feature=spatial$auto$feature,effect_source="global_morans_i",effect=spatial$auto$morans_i,
    n_observations=spatial$auto$n_spots,stringsAsFactors=FALSE)
  o2b_p7_bind_rows(rows)
}

o2b_p7_normalize_spatial_facet_coordinates <- function(data, sample_col = "sample_id", x_col = "pxl_col", y_col = "pxl_row") {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  needed <- c(sample_col, x_col, y_col)
  if (!all(needed %in% names(data))) stop("Missing spatial coordinate columns: ", paste(setdiff(needed, names(data)), collapse = "; "), call. = FALSE)
  rows <- lapply(split(data, as.character(data[[sample_col]]), drop = TRUE), function(z) {
    x <- suppressWarnings(as.numeric(z[[x_col]]))
    y <- suppressWarnings(as.numeric(z[[y_col]]))
    ok <- is.finite(x) & is.finite(y)
    z <- z[ok, , drop = FALSE]
    x <- x[ok]; y <- y[ok]
    if (!nrow(z)) return(z)
    xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
    common_scale <- max(diff(xr), diff(yr), 1)
    z$facet_x <- (x - mean(xr)) / common_scale
    z$facet_y <- -(y - mean(yr)) / common_scale
    z$facet_common_scale <- common_scale
    z
  })
  out <- o2b_p7_bind_rows(rows)
  rownames(out) <- NULL
  out
}

o2b_p7_part6_part7_concordance <- function(root,patient_summary,ridge_patient){
  p6_candidates<-file.path(root,"output","tables","Part6",c("70_final_cell_source_evidence_matrix.csv","64_final_cell_source_evidence_matrix.csv"))
  p6_path<-p6_candidates[file.exists(p6_candidates)][1]
  expected<-data.frame(feature=c("Fibroblast","Myofibroblast","Pericyte","Smooth_Muscle","Myeloid_Macrophage","Epithelial"),
    part6_expected_role=c("robust_primary","replicated_confidence_limited","competing_source","competing_source","heterogeneous_secondary","negative_control_source"),stringsAsFactors=FALSE)
  same<-patient_summary[patient_summary$effect_source=="same_spot_spearman" & patient_summary$feature %in% expected$feature,,drop=FALSE]
  neigh<-patient_summary[patient_summary$effect_source=="neighbor_mean_spearman" & patient_summary$feature %in% expected$feature,,drop=FALSE]
  ridge<-ridge_patient[ridge_patient$feature %in% expected$feature,,drop=FALSE]
  out<-expected
  out$same_spot_median_effect<-same$median_effect[match(out$feature,same$feature)]
  out$same_spot_positive_patients<-same$positive_direction_count[match(out$feature,same$feature)]
  out$neighbor_median_effect<-neigh$median_effect[match(out$feature,neigh$feature)]
  out$ridge_median_coefficient<-ridge$median_effect[match(out$feature,ridge$feature)]
  out$part6_file_used<-if(length(p6_path)&&!is.na(p6_path)) p6_path else "PART6_EVIDENCE_FILE_NOT_FOUND; EXPECTED_ROLE_ONLY"
  out$overall_spatial_support<-ifelse(out$same_spot_median_effect>0 & out$neighbor_median_effect>0,"SAME_SPOT_AND_NEIGHBOR_POSITIVE",
    ifelse(out$same_spot_median_effect>0,"SAME_SPOT_ONLY_OR_NEIGHBOR_UNCLEAR","NOT_SUPPORTED_OR_HETEROGENEOUS"))
  out
}

o2b_p7_make_figures <- function(scores, corr, shift, spatial, dirs, continuous_summary = NULL, threshold_summary = NULL,
                                ridge_patient = NULL, technical_patient = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  reg <- list(); pal <- o2b_p7_palette()
  pretty <- o2b_p7_pretty_feature_labels
  # A: detection by section, with the paired metastasis explicitly identified.
  det <- stats::aggregate(scores$OLFML2B > 0,
    by = list(sample_id = scores$sample_id, patient_id = scores$patient_id, geo_title = scores$geo_title, tissue_role = scores$tissue_role), FUN = mean)
  names(det)[5] <- "detected_fraction"
  det$section_label <- ifelse(nzchar(as.character(det$geo_title)), as.character(det$geo_title), as.character(det$sample_id))
  det$tissue_label <- ifelse(det$tissue_role == "primary_gastric_cancer", "Primary gastric cancer", "Paired metastasis sensitivity")
  p <- ggplot2::ggplot(det, ggplot2::aes(x = section_label, y = detected_fraction, fill = tissue_label)) +
    ggplot2::geom_col(width = 0.76) +
    ggplot2::scale_fill_manual(values = c("Primary gastric cancer" = "#355C7D", "Paired metastasis sensitivity" = "#E07A5F"), name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = ggplot2::expansion(mult = c(0, 0.08))) +
    o2b_p7_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "top") +
    ggplot2::labs(title = "OLFML2B spatial detectability across GSE251950 sections", subtitle = "GC6-PM is retained as paired-metastasis sensitivity and excluded from formal primary-patient inference", x = NULL, y = "OLFML2B-positive spots")
  f <- o2b_p7_save_plot(p, "FIG7A_OLFML2B_detection_by_section", dirs$figures, 8.5, 4.7)
  reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7A", source_table = "08_spatial_target_detection_by_patient.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)

  # B: grey undetected background plus a warm positive-only sequential scale.
  m <- scores[is.finite(suppressWarnings(as.numeric(scores$pxl_row))) & is.finite(suppressWarnings(as.numeric(scores$pxl_col))), , drop = FALSE]
  m <- o2b_p7_normalize_spatial_facet_coordinates(m)
  if (nrow(m)) {
    ord <- o2b_p7_official_gse251950_sample_map()
    title_map <- setNames(as.character(ord$geo_title), as.character(ord$sample_id))
    m$facet_label <- paste0(ifelse(is.na(title_map[m$sample_id]), m$patient_id, title_map[m$sample_id]), "  |  ", m$sample_id)
    facet_levels <- paste0(title_map[ord$sample_id], "  |  ", ord$sample_id)
    m$facet_label <- factor(m$facet_label, levels = facet_levels)
    pos <- m[is.finite(m$OLFML2B) & m$OLFML2B > 0, , drop = FALSE]
    cap <- if (nrow(pos)) as.numeric(stats::quantile(pos$OLFML2B, 0.99, na.rm = TRUE, names = FALSE)) else 1
    if (!is.finite(cap) || cap <= 0) cap <- 1
    if (nrow(pos)) pos$OLFML2B_positive_capped <- pmin(pos$OLFML2B, cap)
    p <- ggplot2::ggplot(m, ggplot2::aes(x = facet_x, y = facet_y)) +
      ggplot2::geom_point(color = "#D7DAE0", size = 0.28, alpha = 0.72) +
      {if (nrow(pos)) ggplot2::geom_point(data = pos, ggplot2::aes(color = OLFML2B_positive_capped), size = 0.42, alpha = 0.96) else NULL} +
      ggplot2::facet_wrap(~facet_label, scales = "fixed", ncol = 5) +
      ggplot2::coord_equal(xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.55), expand = FALSE) +
      ggplot2::scale_color_gradientn(colors = c("#FEE8C8", "#FDBB84", "#FC8D59", "#D7301F", "#7F0000"), limits = c(0, cap), name = "Positive OLFML2B\nlog1p CPM10K", guide = ggplot2::guide_colorbar(title.position = "top", barwidth = 10, barheight = 0.55)) +
      o2b_p7_theme(8) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), legend.position = "bottom", strip.text = ggplot2::element_text(size = 7.3), panel.spacing = grid::unit(1.5, "mm")) +
      ggplot2::labs(title = "OLFML2B-positive spatial spots across all ten Visium sections", subtitle = paste0("Grey = undetected; warm colors = positive-only expression; values capped at the global positive 99th percentile (", formatC(cap, digits = 2, format = "f"), "). Geometry is preserved within each section."), x = NULL, y = NULL)
    f <- o2b_p7_save_plot(p, "FIG7B_OLFML2B_all_section_spatial_maps", dirs$figures, 11.4, 7.2)
    reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7B", source_table = "20_spot_level_continuous_scores.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)
  }

  # C: patient-level same-spot effects with non-deprecated horizontal intervals.
  if (is.data.frame(continuous_summary) && nrow(continuous_summary)) {
    d <- continuous_summary[continuous_summary$effect_source == "same_spot_spearman" & continuous_summary$feature %in% c(o2b_p7_primary_features(), "Pericyte", "Smooth_Muscle", "Myeloid_Macrophage", "Epithelial", "CD8_Cytotoxic"), , drop = FALSE]
    d <- d[is.finite(d$median_effect) & is.finite(d$bootstrap_ci_low) & is.finite(d$bootstrap_ci_high), , drop = FALSE]
    if (nrow(d)) {
      d$feature_label <- pretty(d$feature)
      d$feature_label <- factor(d$feature_label, levels = rev(d$feature_label[order(d$median_effect)]))
      d$direction <- ifelse(d$median_effect >= 0, "Positive", "Negative")
      p <- ggplot2::ggplot(d, ggplot2::aes(x = median_effect, y = feature_label)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#6B7280") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), orientation = "y", width = 0.16, linewidth = 0.55) +
        ggplot2::geom_point(ggplot2::aes(fill = direction), shape = 21, size = 3, color = "#111827", stroke = 0.35) +
        ggplot2::scale_fill_manual(values = c(Positive = "#B2182B", Negative = "#2166AC"), name = NULL) +
        o2b_p7_theme(10) + ggplot2::theme(legend.position = "top") +
        ggplot2::labs(title = "Patient-level continuous spatial associations", subtitle = "Median same-spot Spearman effect after within-patient section aggregation", x = "Median Spearman rho", y = NULL)
      f <- o2b_p7_save_plot(p, "FIG7C_patient_level_same_spot_associations", dirs$figures, 7.6, 6.3)
      reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7C", source_table = "22_patient_level_continuous_associations.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)
    }
  }

  # D: all seven prespecified competing sources, explicit factor levels and labels.
  if (is.data.frame(ridge_patient) && nrow(ridge_patient)) {
    source_levels <- c("Epithelial", "Smooth_Muscle", "Endothelial", "Myeloid_Macrophage", "Pericyte", "Myofibroblast", "Fibroblast")
    d <- ridge_patient[as.character(ridge_patient$feature) %in% source_levels, , drop = FALSE]
    d$feature <- as.character(d$feature)
    d <- d[is.finite(d$median_effect) & is.finite(d$bootstrap_ci_low) & is.finite(d$bootstrap_ci_high), , drop = FALSE]
    if (nrow(d)) {
      d$feature_plot <- factor(d$feature, levels = source_levels, labels = pretty(source_levels))
      d$direction <- ifelse(d$median_effect >= 0, "Positive coefficient", "Negative coefficient")
      p <- ggplot2::ggplot(d, ggplot2::aes(x = median_effect, y = feature_plot)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#6B7280") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), orientation = "y", width = 0.16, linewidth = 0.6) +
        ggplot2::geom_point(ggplot2::aes(fill = direction), shape = 21, size = 3.2, color = "#111827", stroke = 0.4) +
        ggplot2::scale_fill_manual(values = c("Positive coefficient" = "#B2182B", "Negative coefficient" = "#2166AC"), name = NULL) +
        o2b_p7_theme(10) + ggplot2::theme(legend.position = "top") +
        ggplot2::labs(title = "Competing-source spatial model", subtitle = "All seven prespecified source programs; standardized ridge coefficients summarized at patient level", x = "Median standardized coefficient", y = NULL)
      f <- o2b_p7_save_plot(p, "FIG7D_competing_source_ridge", dirs$figures, 7.4, 5.8)
      reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7D", source_table = "41_multivariable_source_disambiguation_patient.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)
    }
  }

  # E: same-spot versus neighbor-mean effects.
  if (is.data.frame(continuous_summary) && nrow(continuous_summary)) {
    d <- continuous_summary[continuous_summary$effect_source %in% c("same_spot_spearman", "neighbor_mean_spearman") & continuous_summary$feature %in% c(o2b_p7_primary_features(), "Pericyte", "Smooth_Muscle", "Myeloid_Macrophage", "Epithelial"), , drop = FALSE]
    d <- d[is.finite(d$median_effect), , drop = FALSE]
    if (nrow(d)) {
      d$feature_label <- pretty(d$feature)
      d$effect_label <- ifelse(d$effect_source == "same_spot_spearman", "Same spot", "KNN6 neighbor mean")
      p <- ggplot2::ggplot(d, ggplot2::aes(x = feature_label, y = median_effect, color = effect_label, shape = effect_label, group = feature_label)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#6B7280") + ggplot2::geom_line(color = "#BFC5CE", linewidth = 0.45) +
        ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.16)) +
        ggplot2::scale_color_manual(values = c("Same spot" = "#B2182B", "KNN6 neighbor mean" = "#2166AC"), name = NULL) +
        ggplot2::scale_shape_manual(values = c("Same spot" = 16, "KNN6 neighbor mean" = 17), name = NULL) +
        o2b_p7_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 38, hjust = 1), legend.position = "top") +
        ggplot2::labs(title = "Same-spot and neighboring spatial effects", x = NULL, y = "Patient-level median rho")
      f <- o2b_p7_save_plot(p, "FIG7E_same_spot_vs_neighbor", dirs$figures, 8.8, 5.2)
      reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7E", source_table = "22_patient_level_continuous_associations.csv;33_patient_level_neighbor_associations.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)
    }
  }

  # F: Moran's I.
  if (is.list(spatial) && is.data.frame(spatial$auto) && nrow(spatial$auto)) {
    d <- spatial$auto[spatial$auto$feature %in% c("OLFML2B", "Fibroblast", "Myofibroblast", "ECM_Remodeling", "TGFb_Response"), , drop = FALSE]
    d$feature_label <- pretty(d$feature)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = feature_label, y = morans_i, group = sample_id)) + ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#6B7280") +
      ggplot2::geom_line(alpha = 0.24, color = "#9CA3AF") + ggplot2::geom_point(size = 1.8, color = "#355C7D", alpha = 0.82) +
      ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 23, size = 3.2, fill = "#F4A261", color = "#111827") +
      o2b_p7_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
      ggplot2::labs(title = "Spatial autocorrelation across sections", subtitle = "Blue dots = sections; orange diamonds = medians", x = NULL, y = "Global Moran's I")
    f <- o2b_p7_save_plot(p, "FIG7F_global_morans_i", dirs$figures, 7.7, 5.1)
    reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7F", source_table = "50_global_moransI_by_section.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)
  }

  # G: self-excluded, spot-spacing-normalized distance sensitivity.
  if (is.list(spatial) && is.data.frame(spatial$nn) && nrow(spatial$nn)) {
    d <- spatial$nn[spatial$nn$status == "OK" & spatial$nn$tissue_role == "primary_gastric_cancer", , drop = FALSE]
    if (nrow(d)) {
      feature_levels <- c("Fibroblast", "Myofibroblast", "Pericyte", "Smooth_Muscle", "Myeloid_Macrophage", "CD8_Cytotoxic")
      d$feature_plot <- factor(as.character(d$to_feature), levels = feature_levels, labels = pretty(feature_levels))
      cols <- o2b_p7_distance_feature_colors(); names(cols) <- pretty(names(cols))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = feature_plot, y = delta_observed_minus_null_spot_units, color = feature_plot)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#6B7280") +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F3F4F6", color = "#6B7280", width = 0.58, linewidth = 0.45) +
        ggplot2::geom_jitter(width = 0.12, size = 2.1, alpha = 0.86) +
        ggplot2::scale_color_manual(values = cols, guide = "none") +
        o2b_p7_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
        ggplot2::labs(title = "Nearest-other-context distance sensitivity", subtitle = "Self spots excluded; distances normalized to each section's median KNN6 spacing. Negative values indicate closer than matched controls.", x = NULL, y = "Observed − matched-null distance (spot-spacing units)", caption = "Same-spot overlap is reported separately and is not counted as zero spatial distance.")
      f <- o2b_p7_save_plot(p, "FIG7G_spatial_distance_sensitivity", dirs$figures, 8.5, 5.3)
      reg[[length(reg) + 1L]] <- data.frame(figure_id = "FIG7G", source_table = "70_spatial_distance_permutation_by_section.csv", output = paste(basename(f), collapse = ";"), stringsAsFactors = FALSE)

      ov <- d[is.finite(d$target_context_overlap_fraction), , drop = FALSE]
      if (nrow(ov)) {
        p2 <- ggplot2::ggplot(ov, ggplot2::aes(x = feature_plot, y = target_context_overlap_fraction, color = feature_plot)) +
          ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F3F4F6", color = "#6B7280", width = 0.58, linewidth = 0.45) +
          ggplot2::geom_jitter(width = 0.12, size = 2.1, alpha = 0.86) +
          ggplot2::scale_color_manual(values = cols, guide = "none") +
          ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
          o2b_p7_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
          ggplot2::labs(title = "Same-spot overlap of OLFML2B-high and context-high regions", subtitle = "Descriptive threshold sensitivity separated from nearest-other-context distance", x = NULL, y = "Target-high spots also context-high")
        f2 <- o2b_p7_save_plot(p2, "SFIG7A_OLFML2B_context_same_spot_overlap", dirs$figures, 8.5, 5.1)
        reg[[length(reg) + 1L]] <- data.frame(figure_id = "SFIG7A", source_table = "72_spatial_target_context_overlap_by_section.csv", output = paste(basename(f2), collapse = ";"), stringsAsFactors = FALSE)
      }
    }
  }
  o2b_p7_write_csv(o2b_p7_bind_rows(reg), file.path(dirs$tables, "98_part7_figure_registry.csv"))
  invisible(TRUE)
}




o2b_p7_existing_repair_ready <- function(dirs) {
  required <- file.path(dirs$tables, c(
    "20_spot_level_continuous_scores.csv", "21_section_level_continuous_associations.csv",
    "22_patient_level_continuous_associations.csv", "27_technical_burden_residual_sensitivity_by_patient.csv",
    "30_spatial_neighbor_graph_audit.csv", "32_section_level_neighbor_associations.csv",
    "41_multivariable_source_disambiguation_patient.csv", "50_global_moransI_by_section.csv",
    "62_threshold_sensitivity_summary.csv", "02_spatial_sample_discovery_audit.csv",
    "04_spatial_loader_and_coordinate_audit.csv", "09_signature_gene_coverage_audit.csv"
  ))
  all(file.exists(required))
}

o2b_p7_read_csv_existing <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

o2b_p7_write_distance_outputs <- function(distance, dirs, bootstrap_B = 2000L) {
  o2b_p7_write_csv(distance, file.path(dirs$tables, "70_spatial_distance_permutation_by_section.csv"))
  tmp <- distance
  tmp$effect_source <- "nearest_other_context_distance_spot_units"
  tmp$feature <- tmp$to_feature
  tmp$effect <- tmp$delta_observed_minus_null_spot_units
  patient <- o2b_p7_summarise_section_effects(tmp[tmp$tissue_role == "primary_gastric_cancer" & tmp$status == "OK", , drop = FALSE], "effect", c("effect_source", "feature"), bootstrap_B)
  if (nrow(patient)) patient$distance_definition <- "nearest OTHER context-high spot; self excluded; normalized by section median KNN6 spacing"
  o2b_p7_write_csv(patient, file.path(dirs$tables, "71_spatial_distance_patient_summary.csv"))
  overlap <- distance[, intersect(c("sample_id", "patient_id", "geo_title", "tissue_role", "to_feature", "n_target_high", "n_context_high", "n_target_context_overlap", "target_context_overlap_fraction", "target_grouping_method", "overlap_definition", "status"), names(distance)), drop = FALSE]
  o2b_p7_write_csv(overlap, file.path(dirs$tables, "72_spatial_target_context_overlap_by_section.csv"))
  ovtmp <- distance; ovtmp$effect_source <- "same_spot_target_context_overlap"; ovtmp$feature <- ovtmp$to_feature; ovtmp$effect <- ovtmp$target_context_overlap_fraction
  ovpatient <- o2b_p7_summarise_section_effects(ovtmp[ovtmp$tissue_role == "primary_gastric_cancer" & is.finite(ovtmp$effect), , drop = FALSE], "effect", c("effect_source", "feature"), bootstrap_B)
  if (nrow(ovpatient)) ovpatient$interpretation_boundary <- "Descriptive threshold overlap; not cell identity, interaction, or causal evidence."
  o2b_p7_write_csv(ovpatient, file.path(dirs$tables, "73_spatial_target_context_overlap_patient_summary.csv"))
  list(patient = patient, overlap = overlap, overlap_patient = ovpatient)
}

o2b_p7_update_v105_boundaries <- function(dirs) {
  claim_path <- file.path(dirs$tables, "100_Part7_spatial_claim_boundary.csv")
  claims <- if (file.exists(claim_path)) o2b_p7_read_csv_existing(claim_path) else data.frame(domain = character(), decision = character(), stringsAsFactors = FALSE)
  if (nrow(claims) && "domain" %in% names(claims)) {
    claims$decision[claims$domain == "distance"] <- "Nearest-other-context distance excludes the query spot itself and is normalized by section median KNN6 spacing; supportive sensitivity only."
  }
  if (!"same_spot_overlap" %in% claims$domain) claims <- rbind(claims, data.frame(domain = "same_spot_overlap", decision = "Target-high/context-high overlap is reported separately from distance and is descriptive threshold sensitivity only.", stringsAsFactors = FALSE))
  o2b_p7_write_csv(claims, claim_path)
  key_path <- file.path(dirs$tables, "82_part7_key_result_summary.csv")
  key <- if (file.exists(key_path)) o2b_p7_read_csv_existing(key_path) else data.frame(item = character(), value = character(), stringsAsFactors = FALSE)
  if (nrow(key) && "item" %in% names(key)) key$value[key$item == "version"] <- OLFML2B_PART7_VERSION
  add <- data.frame(item = c("distance_definition", "spatial_map_palette"), value = c("nearest-other; self-excluded; normalized to median KNN6 spot spacing", "grey undetected background plus warm positive-only 99th-percentile-capped scale"), stringsAsFactors = FALSE)
  key <- key[!key$item %in% add$item, , drop = FALSE]; key <- rbind(key, add)
  o2b_p7_write_csv(key, key_path)
  go_path <- file.path(dirs$tables, "99_Part7_spatial_go_no_go.csv")
  if (file.exists(go_path)) {
    go <- o2b_p7_read_csv_existing(go_path); go$distance_layer_status <- "REPAIRED_SELF_EXCLUDED_AND_SPOT_SPACING_NORMALIZED"; go$final_gene_lock <- FALSE
    o2b_p7_write_csv(go, go_path)
  }
  supersession <- data.frame(
    superseded_output = c("70_spatial_distance_permutation_by_section.csv v1.0.4", "71_spatial_distance_patient_summary.csv v1.0.4", "FIG7D v1.0.4", "FIG7G v1.0.4", "FIG7B v1.0.4"),
    replacement = c("v1.0.5 self-excluded normalized distance", "v1.0.5 patient summary in spot-spacing units", "v1.0.5 complete seven-source forest", "v1.0.5 nearest-other normalized distance", "v1.0.5 grey-undetected/warm-positive spatial map"),
    reason = c("self-distance could create zero and extreme effects", "derived from invalid unnormalized/self-inclusive distance", "factor/label conversion omitted source rows", "uses superseded distance values", "low expression was visually dominant and positive spots lacked contrast"),
    stringsAsFactors = FALSE
  )
  o2b_p7_write_csv(supersession, file.path(dirs$tables, "74_v105_supersession_audit.csv"))
  claims
}

o2b_p7_repair_from_existing <- function(root, dirs, make_figures = TRUE, bootstrap_B = 2000L, distance_permutation_B = 999L, log_file = NULL) {
  o2b_p7_log("INFO", "Validated continuation mode: reusing existing 32,291-spot score table and frozen core Part7 statistics; raw Visium matrices will not be reread.", log_file = log_file)
  scores <- o2b_p7_read_csv_existing(file.path(dirs$tables, "20_spot_level_continuous_scores.csv"))
  if (!all(c("sample_id", "patient_id", "OLFML2B", "pxl_row", "pxl_col") %in% names(scores))) stop("Existing Part7 score table lacks required columns.", call. = FALSE)
  if (length(unique(scores$sample_id)) != 10L || length(unique(stats::na.omit(scores$patient_id))) != 9L) stop("Existing score table does not satisfy the 10-section/9-patient contract.", call. = FALSE)
  corr <- o2b_p7_read_csv_existing(file.path(dirs$tables, "21_section_level_continuous_associations.csv"))
  patient_summary <- o2b_p7_read_csv_existing(file.path(dirs$tables, "22_patient_level_continuous_associations.csv"))
  technical_patient <- o2b_p7_read_csv_existing(file.path(dirs$tables, "27_technical_burden_residual_sensitivity_by_patient.csv"))
  ridge_patient <- o2b_p7_read_csv_existing(file.path(dirs$tables, "41_multivariable_source_disambiguation_patient.csv"))
  threshold_summary <- o2b_p7_read_csv_existing(file.path(dirs$tables, "62_threshold_sensitivity_summary.csv"))
  auto <- o2b_p7_read_csv_existing(file.path(dirs$tables, "50_global_moransI_by_section.csv"))
  bivar <- o2b_p7_read_csv_existing(file.path(dirs$tables, "32_section_level_neighbor_associations.csv"))
  graph <- o2b_p7_read_csv_existing(file.path(dirs$tables, "30_spatial_neighbor_graph_audit.csv"))
  distance <- o2b_p7_distance_stats(scores, distance_permutation_B = distance_permutation_B)
  spatial <- list(auto = auto, bivar = bivar, graph = graph, nn = distance)
  distance_outputs <- o2b_p7_write_distance_outputs(distance, dirs, bootstrap_B = bootstrap_B)
  claims <- o2b_p7_update_v105_boundaries(dirs)
  figure_audit <- data.frame(requested = isTRUE(make_figures), status = ifelse(isTRUE(make_figures), "PENDING", "SKIPPED_BY_USER"), error_message = "", stringsAsFactors = FALSE)
  if (isTRUE(make_figures)) {
    tryCatch({
      o2b_p7_make_figures(scores, corr, data.frame(), spatial, dirs, patient_summary, threshold_summary, ridge_patient, technical_patient)
      figure_audit$status <- if (file.exists(file.path(dirs$figures, "Part7_figure_export_errors.csv"))) "COMPLETED_WITH_FORMAT_WARNINGS" else "PASS"
    }, error = function(e) {
      figure_audit$status <- "FAILED_NONFATAL"; figure_audit$error_message <- conditionMessage(e)
      o2b_p7_log("WARN", "Part7 v1.0.5 figure regeneration failed non-fatally: ", conditionMessage(e), log_file = log_file)
    })
  }
  o2b_p7_write_csv(figure_audit, file.path(dirs$tables, "98b_part7_figure_export_audit.csv"))
  object_path <- file.path(dirs$objects, "Part7_OLFML2B_spatial_index.rds")
  index <- if (file.exists(object_path)) tryCatch(readRDS(object_path), error = function(e) list()) else list()
  index$version <- OLFML2B_PART7_VERSION; index$generated_at <- o2b_p7_ts(); index$root <- root; index$dirs <- dirs
  index$scores <- scores; index$corr <- corr; index$patient_summary <- patient_summary; index$spatial <- spatial
  index$distance_patient <- distance_outputs$patient; index$distance_overlap <- distance_outputs$overlap
  index$distance_overlap_patient <- distance_outputs$overlap_patient; index$claim_boundary <- claims
  index$figure_audit <- figure_audit; index$final_gene_lock <- FALSE
  o2b_p7_save_rds(index, object_path)
  p06 <- file.path(root, "output", "objects", "OLFML2B_Part0_6_complete_index.rds")
  integrated <- list(part0_6 = if (file.exists(p06)) readRDS(p06) else NULL, part7 = index, version = OLFML2B_PART7_VERSION, final_gene_lock = FALSE)
  o2b_p7_save_rds(integrated, file.path(dirs$objects, "OLFML2B_Part0_7_complete_index.rds"))
  o2b_p7_log("INFO", "OLFML2B Part7 v1.0.5 repair complete | ", nrow(scores), " spots reused | nearest-other normalized distance rebuilt | figures optimized | final_gene_lock=FALSE", log_file = log_file)
  invisible(index)
}

run_olfml2b_part7_spatial_transcriptomics <- function(
  root = "D:/OLFML2B_STAD",
  spatial_dir = file.path(root, "data", "raw", "spatial", "GSE251950"),
  standardized_dir = file.path(root, "data", "processed", "spatial", "GSE251950_standardized"),
  output_subdir = "Part7", make_figures = TRUE,
  bootstrap_B = 2000L, moran_permutation_B = 999L,
  distance_permutation_B = 999L, reuse_existing_validated_results = TRUE, ...
) {
  root <- normalizePath(root,winslash="/",mustWork=FALSE)
  dirs <- o2b_p7_dirs(root,output_subdir)
  log_file <- file.path(dirs$logs,"Part7_spatial_transcriptomics.log")
  o2b_p7_log("INFO","Starting OLFML2B Part7 | ",OLFML2B_PART7_VERSION,log_file=log_file)
  o2b_p7_log("INFO","Root: ",root,log_file=log_file)
  o2b_p7_log("INFO","Raw spatial dir: ",spatial_dir,log_file=log_file)
  o2b_p7_log("INFO","Standardized dir: ",standardized_dir,log_file=log_file)
  params <- data.frame(parameter=c("version","dataset","target","root","spatial_dir","standardized_dir","inference_unit","normalization","bootstrap_B","moran_permutation_B","distance_permutation_B","reuse_existing_validated_results","final_gene_lock"),
    value=c(OLFML2B_PART7_VERSION,"GSE251950_10x_Visium_10_sections_9_patients","OLFML2B",root,spatial_dir,standardized_dir,
      "patient_after_within_patient_section_aggregation","full_matrix_log1p_CPM10000",bootstrap_B,moran_permutation_B,distance_permutation_B,reuse_existing_validated_results,FALSE),stringsAsFactors=FALSE)
  o2b_p7_write_csv(params,file.path(dirs$tables,"00_part7_run_parameters.csv"))
  o2b_p7_write_csv(o2b_p7_pkg_versions(c("Matrix","ggplot2")),file.path(dirs$reports,"Part7_package_versions.csv"))
  if (isTRUE(reuse_existing_validated_results) && o2b_p7_existing_repair_ready(dirs)) {
    return(o2b_p7_repair_from_existing(root, dirs, make_figures = make_figures, bootstrap_B = bootstrap_B, distance_permutation_B = distance_permutation_B, log_file = log_file))
  }
  sample_dirs <- o2b_p7_find_sample_dirs(standardized_dir=standardized_dir,raw_dir=spatial_dir)
  manifest <- o2b_p7_attach_official_sample_metadata(data.frame(sample_id=basename(sample_dirs),sample_dir=sample_dirs,stringsAsFactors=FALSE))
  o2b_p7_write_csv(manifest,file.path(dirs$tables,"02_spatial_sample_discovery_audit.csv"))
  official <- o2b_p7_official_gse251950_sample_map()
  o2b_p7_write_csv(official,file.path(dirs$tables,"03_official_section_patient_map.csv"))
  map_audit <- data.frame(n_sections=nrow(manifest),n_mapped=sum(manifest$patient_mapping_status=="OFFICIAL_GEO_MAP",na.rm=TRUE),
    n_patients=length(unique(stats::na.omit(manifest$patient_id))),
    status=ifelse(nrow(manifest)==10L&&length(unique(stats::na.omit(manifest$patient_id)))==9L&&all(manifest$patient_mapping_status=="OFFICIAL_GEO_MAP"),"PASS_10_SECTIONS_9_PATIENTS","REVIEW_MAPPING"),stringsAsFactors=FALSE)
  o2b_p7_write_csv(map_audit,file.path(dirs$tables,"03b_section_patient_mapping_audit.csv"))
  score_rows<-list();qc_rows<-list();diag_rows<-list();coverage_rows<-list()
  for(i in seq_along(sample_dirs)){
    sid<-basename(sample_dirs[i]);o2b_p7_log("INFO","Loading sample ",i,"/",length(sample_dirs),": ",sid,log_file=log_file)
    one<-tryCatch(o2b_p7_load_one_sample(sample_dirs[i],log_file=log_file),error=function(e) list(ok=FALSE,scores=data.frame(),
      qc=data.frame(sample_id=sid,sample_dir=sample_dirs[i],parser_status="FAIL_RUNTIME_ERROR",loader_notes=conditionMessage(e),stringsAsFactors=FALSE),feature_diag=data.frame(),coverage=data.frame()))
    score_rows[[length(score_rows)+1L]]<-one$scores;qc_rows[[length(qc_rows)+1L]]<-one$qc
    diag_rows[[length(diag_rows)+1L]]<-one$feature_diag;coverage_rows[[length(coverage_rows)+1L]]<-one$coverage
  }
  scores<-o2b_p7_attach_official_sample_metadata(o2b_p7_bind_rows(score_rows));qc<-o2b_p7_bind_rows(qc_rows)
  diag<-o2b_p7_bind_rows(diag_rows);coverage<-o2b_p7_bind_rows(coverage_rows)
  o2b_p7_write_csv(qc,file.path(dirs$tables,"04_spatial_loader_and_coordinate_audit.csv"))
  o2b_p7_write_csv(diag,file.path(dirs$tables,"05_spatial_feature_column_diagnostic.csv"))
  o2b_p7_write_csv(coverage,file.path(dirs$tables,"09_signature_gene_coverage_audit.csv"))
  identity<-data.frame(sample_id=qc$sample_id,exact_OLFML2B_found=qc$target_OLFML2B_found,
    OLFML2A_detected=coverage$n_genes_measured[match(paste(qc$sample_id,"OLFML2A"),paste(coverage$sample_id,coverage$signature))]>0,
    OLFM2_detected=coverage$n_genes_measured[match(paste(qc$sample_id,"OLFM2"),paste(coverage$sample_id,coverage$signature))]>0,
    matching_policy="exact_symbol_only; OLFML2A and OLFM2 excluded from target",stringsAsFactors=FALSE)
  o2b_p7_write_csv(identity,file.path(dirs$tables,"06_exact_OLFML2B_identity_audit.csv"))
  o2b_p7_write_csv(scores,file.path(dirs$tables,"20_spot_level_continuous_scores.csv"))
  if(!nrow(scores)){
    go<-data.frame(status="NO_GO_NO_SCORE_TABLE",final_gene_lock=FALSE,stringsAsFactors=FALSE)
    o2b_p7_write_csv(go,file.path(dirs$tables,"99_Part7_spatial_go_no_go.csv"));return(invisible(list(status="NO_GO_NO_SCORE_TABLE")))
  }
  det_section<-o2b_p7_bind_rows(lapply(split(scores,scores$sample_id),function(d) data.frame(sample_id=d$sample_id[1],patient_id=d$patient_id[1],
    n_spots=nrow(d),n_positive=sum(d$OLFML2B>0,na.rm=TRUE),positive_fraction=mean(d$OLFML2B>0,na.rm=TRUE),
    median_positive_expression=if(sum(d$OLFML2B>0,na.rm=TRUE)) stats::median(d$OLFML2B[d$OLFML2B>0],na.rm=TRUE) else NA_real_,stringsAsFactors=FALSE)))
  det_patient<-o2b_p7_bind_rows(lapply(split(det_section,det_section$patient_id),function(d) data.frame(patient_id=d$patient_id[1],n_sections=nrow(d),
    total_spots=sum(d$n_spots),total_positive=sum(d$n_positive),positive_fraction=sum(d$n_positive)/sum(d$n_spots),
    evaluable=any(d$n_positive>=20L|d$positive_fraction>=0.005),stringsAsFactors=FALSE)))
  o2b_p7_write_csv(det_section,file.path(dirs$tables,"07_spatial_target_detection_by_section.csv"))
  o2b_p7_write_csv(det_patient,file.path(dirs$tables,"08_spatial_target_detection_by_patient.csv"))
  overlap<-o2b_p7_classifier_overlap_audit();o2b_p7_write_csv(overlap,file.path(dirs$tables,"10_signature_classifier_overlap_audit.csv"))
  corr<-o2b_p7_attach_official_sample_metadata(o2b_p7_corr_tables(scores))
  shift<-o2b_p7_high_low_tables(scores);threshold_audit<-attr(shift,"thresholds");shift<-o2b_p7_attach_official_sample_metadata(shift)
  spatial<-o2b_p7_spatial_stats(scores,distance_permutation_B,moran_permutation_B)
  effects<-o2b_p7_continuous_sample_level_effects(corr,spatial)
  effects_meta<-o2b_p7_attach_official_sample_metadata(effects)
  patient_summary_all_sections<-o2b_p7_summarise_section_effects(effects_meta,"effect",c("effect_source","feature"),bootstrap_B)
  primary_effects<-effects_meta[effects_meta$tissue_role=="primary_gastric_cancer",,drop=FALSE]
  patient_summary<-o2b_p7_summarise_section_effects(primary_effects,"effect",c("effect_source","feature"),bootstrap_B)
  if(nrow(patient_summary)) patient_summary$primary_tissue_policy<-"one_primary_section_per_patient; paired_metastasis excluded from formal main inference"
  technical<-o2b_p7_technical_residual_sensitivity(scores)
  technical_primary_section<-technical$section[technical$section$tissue_role=="primary_gastric_cancer",,drop=FALSE]
  technical_primary_tmp<-technical_primary_section; technical_primary_tmp$effect_source<-"technical_residual_spearman"; technical_primary_tmp$effect<-technical_primary_tmp$residual_rho
  technical$patient_primary<-o2b_p7_summarise_section_effects(technical_primary_tmp[technical_primary_tmp$status=="OK",,drop=FALSE],"effect",c("effect_source","feature"),bootstrap_B)
  ridge<-o2b_p7_source_competition(scores)
  ridge_primary_section<-ridge$section[ridge$section$tissue_role=="primary_gastric_cancer",,drop=FALSE]
  ridge_primary_tmp<-ridge_primary_section[ridge_primary_section$predictor %in% o2b_p7_source_features(),,drop=FALSE]
  ridge_primary_tmp$effect_source<-"ridge_competing_source";ridge_primary_tmp$feature<-ridge_primary_tmp$predictor;ridge_primary_tmp$effect<-ridge_primary_tmp$standardized_ridge_coefficient
  ridge$patient_primary<-o2b_p7_summarise_section_effects(ridge_primary_tmp,"effect",c("effect_source","feature"),bootstrap_B)
  threshold<-o2b_p7_threshold_sensitivity(scores,10L)
  concordance<-o2b_p7_part6_part7_concordance(root,patient_summary,ridge$patient_primary)
  o2b_p7_write_csv(corr,file.path(dirs$tables,"21_section_level_continuous_associations.csv"))
  o2b_p7_write_csv(patient_summary,file.path(dirs$tables,"22_patient_level_continuous_associations.csv"))
  o2b_p7_write_csv(patient_summary_all_sections,file.path(dirs$tables,"22b_all_section_within_patient_sensitivity.csv"))
  o2b_p7_write_csv(patient_summary[patient_summary$feature %in% o2b_p7_primary_features(),,drop=FALSE],file.path(dirs$tables,"23_primary_family_patient_inference.csv"))
  o2b_p7_write_csv(patient_summary[patient_summary$feature %in% c("Pericyte","Smooth_Muscle","Myeloid_Macrophage","Endothelial","Epithelial"),,drop=FALSE],file.path(dirs$tables,"24_competing_source_patient_inference.csv"))
  o2b_p7_write_csv(technical$section,file.path(dirs$tables,"26_technical_burden_residual_sensitivity_by_section.csv"))
  o2b_p7_write_csv(technical$patient_primary,file.path(dirs$tables,"27_technical_burden_residual_sensitivity_by_patient.csv"))
  o2b_p7_write_csv(technical$patient,file.path(dirs$tables,"27b_technical_residual_all_section_sensitivity.csv"))
  o2b_p7_write_csv(spatial$graph,file.path(dirs$tables,"30_spatial_neighbor_graph_audit.csv"))
  o2b_p7_write_csv(spatial$bivar,file.path(dirs$tables,"32_section_level_neighbor_associations.csv"))
  neigh_tmp<-data.frame();if(nrow(spatial$bivar)){neigh_tmp<-spatial$bivar;neigh_tmp$effect_source<-"neighbor_mean_spearman";neigh_tmp$feature<-neigh_tmp$neighbor_feature;neigh_tmp$effect<-neigh_tmp$rho_target_vs_neighbor_mean}
  neigh_tmp<-o2b_p7_attach_official_sample_metadata(neigh_tmp)
  neigh_patient<-o2b_p7_summarise_section_effects(neigh_tmp[neigh_tmp$tissue_role=="primary_gastric_cancer",,drop=FALSE],"effect",c("effect_source","feature"),bootstrap_B)
  o2b_p7_write_csv(neigh_patient,file.path(dirs$tables,"33_patient_level_neighbor_associations.csv"))
  o2b_p7_write_csv(ridge$section,file.path(dirs$tables,"40_multivariable_source_disambiguation_section.csv"))
  o2b_p7_write_csv(ridge$patient_primary,file.path(dirs$tables,"41_multivariable_source_disambiguation_patient.csv"))
  o2b_p7_write_csv(ridge$patient,file.path(dirs$tables,"41b_ridge_all_section_sensitivity.csv"))
  o2b_p7_write_csv(spatial$auto,file.path(dirs$tables,"50_global_moransI_by_section.csv"))
  moran_tmp<-spatial$auto;moran_tmp$effect_source<-"global_morans_i";moran_tmp$effect<-moran_tmp$morans_i
  moran_tmp<-o2b_p7_attach_official_sample_metadata(moran_tmp)
  moran_patient<-o2b_p7_summarise_section_effects(moran_tmp[moran_tmp$tissue_role=="primary_gastric_cancer",,drop=FALSE],"effect",c("effect_source","feature"),bootstrap_B)
  o2b_p7_write_csv(moran_patient,file.path(dirs$tables,"51_global_moransI_patient_summary.csv"))
  o2b_p7_write_csv(threshold$audit,file.path(dirs$tables,"60_threshold_definition_audit.csv"))
  o2b_p7_write_csv(threshold$by_sample,file.path(dirs$tables,"61_threshold_sensitivity_by_section.csv"))
  o2b_p7_write_csv(threshold$summary,file.path(dirs$tables,"62_threshold_sensitivity_summary.csv"))
  distance_outputs <- o2b_p7_write_distance_outputs(spatial$nn, dirs, bootstrap_B = bootstrap_B)
  distance_patient <- distance_outputs$patient
  o2b_p7_write_csv(concordance,file.path(dirs$tables,"80_Part6_Part7_concordance_matrix.csv"))
  evaluable_patients<-sum(det_patient$evaluable %in% TRUE)
  data_gate<-nrow(manifest)==10L && length(unique(stats::na.omit(manifest$patient_id)))==9L && all(startsWith(qc$parser_status,"OK")) && all(qc$n_spots_with_coordinates>0)
  primary<-patient_summary[patient_summary$effect_source=="same_spot_spearman" & patient_summary$feature %in% o2b_p7_primary_features(),,drop=FALSE]
  strong<-nrow(primary)>0 && any(primary$n_patients>=8L & primary$positive_direction_count>=7L & primary$exact_signflip_fdr<0.05 & primary$loo_direction_stability>=0.875,na.rm=TRUE)
  supportive<-nrow(primary)>0 && any(primary$n_patients>=6L & primary$positive_direction_count>=6L & primary$median_effect>0,na.rm=TRUE)
  grade<-if(!data_gate||evaluable_patients<6L) "TARGET_NOT_SPATIALLY_EVALUABLE" else if(strong) "STRONG_PATIENT_LEVEL_SPATIAL_SUPPORT" else if(supportive) "SUPPORTIVE_DIRECTIONAL_SPATIAL_CONTEXT" else "HETEROGENEOUS_OR_INCONCLUSIVE"
  go<-data.frame(status=ifelse(data_gate&&evaluable_patients>=6L,"GO_PATIENT_LEVEL_SPATIAL_CONTEXT_AVAILABLE","REVIEW_SPATIAL_INFERENCE_GATE"),
    spatial_evidence_grade=grade,n_sections=nrow(manifest),n_patients=length(unique(stats::na.omit(manifest$patient_id))),n_spots=nrow(scores),
    n_evaluable_patients=evaluable_patients,data_gate_pass=data_gate,final_gene_lock=FALSE,
    interpretation="Orthogonal spatial tissue-context support only; not causal, prognostic, treatment-selection, cell-cell interaction, or CAF-specific proof.",stringsAsFactors=FALSE)
  o2b_p7_write_csv(go,file.path(dirs$tables,"99_Part7_spatial_go_no_go.csv"))
  claims<-data.frame(domain=c("analysis_unit","primary_inference","source_claim","state_claim","thresholds","distance","same_spot_overlap","causality","prognosis","final_gene_lock"),
    decision=c("Patient; sections and spots are nested observations","Primary gastric-cancer section per patient is formal main analysis; paired metastasis is sensitivity/context only","Preferential spatial association only; not unique cell source or deconvolution","ECM/TGFb co-variation only; not pathway activation by OLFML2B","Sensitivity and visualization only","Nearest-other-context distance excludes the query spot itself and is normalized by section median KNN6 spacing; supportive sensitivity only","Target-high/context-high overlap is reported separately from distance and is descriptive threshold sensitivity only","Not supported","Not supported by Part7","FALSE"),stringsAsFactors=FALSE)
  o2b_p7_write_csv(claims,file.path(dirs$tables,"100_Part7_spatial_claim_boundary.csv"))
  key<-data.frame(item=c("version","dataset","sections","patients","spots","evaluable_patients","spatial_grade","primary_features","competing_sources","distance_definition","spatial_map_palette","final_gene_lock"),
    value=c(OLFML2B_PART7_VERSION,"GSE251950",nrow(manifest),length(unique(stats::na.omit(manifest$patient_id))),nrow(scores),evaluable_patients,grade,
      paste(o2b_p7_primary_features(),collapse=";"),paste(c("Pericyte","Smooth_Muscle","Myeloid_Macrophage","Endothelial","Epithelial"),collapse=";"),"nearest-other; self-excluded; normalized to median KNN6 spot spacing","grey undetected plus warm positive-only 99th-percentile-capped scale",FALSE),stringsAsFactors=FALSE)
  o2b_p7_write_csv(key,file.path(dirs$tables,"82_part7_key_result_summary.csv"))
  figure_audit <- data.frame(
    requested = isTRUE(make_figures),
    status = ifelse(isTRUE(make_figures), "PENDING", "SKIPPED_BY_USER"),
    error_message = "",
    stringsAsFactors = FALSE
  )
  if(isTRUE(make_figures)) {
    tryCatch({
      o2b_p7_make_figures(scores,corr,shift,spatial,dirs,patient_summary,threshold$summary,ridge$patient_primary,technical$patient_primary)
      figure_audit$status <- if (file.exists(file.path(dirs$figures, "Part7_figure_export_errors.csv"))) "COMPLETED_WITH_FORMAT_WARNINGS" else "PASS"
    }, error=function(e){
      figure_audit$status <- "FAILED_NONFATAL"
      figure_audit$error_message <- conditionMessage(e)
      o2b_p7_log("WARN","Part7 figure generation failed non-fatally: ",conditionMessage(e),log_file=log_file)
    })
  }
  o2b_p7_write_csv(figure_audit,file.path(dirs$tables,"98b_part7_figure_export_audit.csv"))
  index<-list(version=OLFML2B_PART7_VERSION,generated_at=o2b_p7_ts(),root=root,spatial_dir=spatial_dir,dirs=dirs,manifest=manifest,qc=qc,coverage=coverage,
    scores=scores,corr=corr,patient_summary=patient_summary,patient_summary_all_sections=patient_summary_all_sections,technical=technical,ridge=ridge,spatial=spatial,threshold=threshold,concordance=concordance,distance_patient=distance_outputs$patient,distance_overlap=distance_outputs$overlap,distance_overlap_patient=distance_outputs$overlap_patient,figure_audit=figure_audit,go_no_go=go,claim_boundary=claims,final_gene_lock=FALSE)
  o2b_p7_save_rds(index,file.path(dirs$objects,"Part7_OLFML2B_spatial_index.rds"))
  p06<-file.path(root,"output","objects","OLFML2B_Part0_6_complete_index.rds")
  integrated<-list(part0_6=if(file.exists(p06)) readRDS(p06) else NULL,part7=index,version=OLFML2B_PART7_VERSION,final_gene_lock=FALSE)
  o2b_p7_save_rds(integrated,file.path(dirs$objects,"OLFML2B_Part0_7_complete_index.rds"))
  o2b_p7_log("INFO","OLFML2B Part7 complete | ",nrow(scores)," spots | ",nrow(manifest)," sections | ",length(unique(stats::na.omit(manifest$patient_id)))," patients | grade=",grade," | final_gene_lock=FALSE",log_file=log_file)
  invisible(index)
}

run_part7 <- run_olfml2b_part7_spatial_transcriptomics
