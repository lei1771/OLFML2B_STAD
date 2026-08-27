# ==============================================================================
# File 01/06: TCGA helpers and Part1 TCGA-STAD acquisition/processing.
# ==============================================================================
.compact_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.compact_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.compact_env_valid <- nzchar(.compact_env_root) && dir.exists(.compact_env_root) &&
  file.exists(file.path(.compact_env_root, "00_OLFML2B_PART0_CONFIG_CORE.R")) &&
  file.exists(file.path(.compact_env_root, "05_OLFML2B_PART5_PDC_PRODUCTION.R"))
.compact_code_root <- if (.compact_env_valid) {
  normalizePath(.compact_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.compact_entry) && file.exists(.compact_entry)) {
  dirname(normalizePath(.compact_entry, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .compact_code_root)
if (!exists("run_part0", mode = "function")) {
  sys.source(file.path(.compact_code_root, "00_OLFML2B_PART0_CONFIG_CORE.R"), envir = environment(), chdir = FALSE)
}

# ==============================================================================
# TCGA-STAD acquisition and harmonisation helpers for the OLFML2B project.
# Uses official GDC through TCGAbiolinks with retry, TLS transport bridging, query caching and full audit.
# ==============================================================================

olfml2b_gdc_status <- function(log_file = NULL) {
  olfml2b_require_packages(c("jsonlite", "curl"), log_file)
  url <- "https://api.gdc.cancer.gov/status"
  profiles <- c("default", "ipv4")
  rows <- vector("list", length(profiles))

  for (i in seq_along(profiles)) {
    profile <- profiles[i]
    handle <- curl::new_handle(
      connecttimeout = 30,
      timeout = 60,
      followlocation = TRUE,
      useragent = paste0("OLFML2B-STAD/", olfml2b_version())
    )
    if (identical(profile, "ipv4")) {
      curl::handle_setopt(handle, ipresolve = 1L)
    }

    started <- Sys.time()
    fetched <- tryCatch(
      curl::curl_fetch_memory(url, handle = handle),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

    if (inherits(fetched, "error")) {
      rows[[i]] <- data.frame(
        endpoint = url,
        network_profile = profile,
        reachable = FALSE,
        http_status = NA_integer_,
        api_status = NA_character_,
        data_release = NA_character_,
        elapsed_seconds = elapsed,
        error = conditionMessage(fetched),
        timestamp = olfml2b_timestamp(),
        stringsAsFactors = FALSE
      )
      next
    }

    txt <- tryCatch(rawToChar(fetched$content), error = function(e) "")
    parsed <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
    api_status <- if (is.list(parsed) && !is.null(parsed$status)) as.character(parsed$status) else NA_character_
    data_release <- if (is.list(parsed) && !is.null(parsed$data_release)) as.character(parsed$data_release) else NA_character_
    ok <- isTRUE(fetched$status_code >= 200L && fetched$status_code < 300L)
    rows[[i]] <- data.frame(
      endpoint = url,
      network_profile = profile,
      reachable = ok,
      http_status = as.integer(fetched$status_code),
      api_status = api_status,
      data_release = data_release,
      elapsed_seconds = elapsed,
      error = if (ok) NA_character_ else paste0("HTTP ", fetched$status_code),
      timestamp = olfml2b_timestamp(),
      stringsAsFactors = FALSE
    )
  }

  out <- olfml2b_bind_rows_safe(rows)
  if (any(out$reachable %in% TRUE)) {
    good <- out[out$reachable %in% TRUE, , drop = FALSE]
    olfml2b_log(
      "INFO",
      "GDC status endpoint reachable via ", paste(unique(good$network_profile), collapse = ", "),
      " | API status=", paste(unique(stats::na.omit(good$api_status)), collapse = ", "),
      log_file = log_file
    )
  } else {
    olfml2b_log(
      "WARN",
      "GDC status endpoint was not reachable by either default or IPv4 curl. ",
      "The RNA query will still be retried because the status endpoint can fail transiently.",
      log_file = log_file
    )
  }
  out
}

olfml2b_install_tcgabiolinks_status_bridge <- function(log_file = NULL) {
  ns <- asNamespace("TCGAbiolinks")
  if (!exists("getGDCInfo", envir = ns, inherits = FALSE)) {
    olfml2b_log(
      "WARN",
      "TCGAbiolinks::getGDCInfo was not found; no status transport bridge was installed.",
      log_file = log_file
    )
    return(list(applied = FALSE, original = NULL))
  }

  original <- get("getGDCInfo", envir = ns, inherits = FALSE)
  bridge <- function() {
    url <- "https://api.gdc.cancer.gov/status"
    errors <- character()
    for (profile in c("default", "ipv4")) {
      handle <- curl::new_handle(
        connecttimeout = 30,
        timeout = 60,
        followlocation = TRUE,
        useragent = "OLFML2B-STAD GDC status bridge"
      )
      if (identical(profile, "ipv4")) curl::handle_setopt(handle, ipresolve = 1L)
      fetched <- tryCatch(
        curl::curl_fetch_memory(url, handle = handle),
        error = function(e) e
      )
      if (inherits(fetched, "error")) {
        errors <- c(errors, paste0(profile, ": ", conditionMessage(fetched)))
        next
      }
      if (!isTRUE(fetched$status_code >= 200L && fetched$status_code < 300L)) {
        errors <- c(errors, paste0(profile, ": HTTP ", fetched$status_code))
        next
      }
      txt <- rawToChar(fetched$content)
      parsed <- tryCatch(
        jsonlite::fromJSON(txt, simplifyDataFrame = TRUE),
        error = function(e) e
      )
      if (!inherits(parsed, "error")) return(parsed)
      errors <- c(errors, paste0(profile, ": JSON parse: ", conditionMessage(parsed)))
    }
    stop(
      "Unable to retrieve GDC status through the explicit curl bridge: ",
      paste(errors, collapse = " | "),
      call. = FALSE
    )
  }

  applied <- tryCatch({
    utils::assignInNamespace("getGDCInfo", bridge, ns = "TCGAbiolinks")
    TRUE
  }, error = function(e) {
    olfml2b_log(
      "WARN",
      "Unable to install temporary TCGAbiolinks GDC status bridge: ", conditionMessage(e),
      log_file = log_file
    )
    FALSE
  })
  if (applied) {
    olfml2b_log(
      "INFO",
      "Installed temporary TCGAbiolinks GDC status bridge using explicit curl with IPv4 fallback.",
      log_file = log_file
    )
  }
  list(applied = applied, original = original)
}

olfml2b_restore_tcgabiolinks_status_bridge <- function(bridge_state, log_file = NULL) {
  if (is.null(bridge_state) || !isTRUE(bridge_state$applied) || !is.function(bridge_state$original)) {
    return(invisible(FALSE))
  }
  ok <- tryCatch({
    utils::assignInNamespace("getGDCInfo", bridge_state$original, ns = "TCGAbiolinks")
    TRUE
  }, error = function(e) {
    olfml2b_log(
      "WARN",
      "Unable to restore original TCGAbiolinks::getGDCInfo: ", conditionMessage(e),
      log_file = log_file
    )
    FALSE
  })
  invisible(ok)
}

olfml2b_gdc_query_cache_valid <- function(path, log_file = NULL) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  query <- tryCatch(readRDS(path), error = function(e) {
    olfml2b_log("WARN", "Unable to read cached GDC query: ", conditionMessage(e), log_file = log_file)
    NULL
  })
  if (is.null(query)) return(NULL)
  n <- tryCatch(nrow(TCGAbiolinks::getResults(query)), error = function(e) 0L)
  if (!is.finite(n) || n < 1L) {
    olfml2b_log("WARN", "Cached GDC query is invalid or empty: ", path, log_file = log_file)
    return(NULL)
  }
  attr(query, "workflow_used") <- attr(query, "workflow_used") %||% "STAR - Counts"
  olfml2b_log("INFO", "Using cached GDC STAR query: ", path, " | files=", n, log_file = log_file)
  query
}

olfml2b_tcga_query_rna <- function(
  project = "TCGA-STAD",
  sample_types = c("Primary Tumor", "Solid Tissue Normal"),
  query_cache = NULL,
  log_file = NULL
) {
  olfml2b_require_packages(c("TCGAbiolinks", "curl", "jsonlite"), log_file)

  old_url_method <- getOption("url.method")
  options(url.method = "libcurl")
  on.exit(options(url.method = old_url_method), add = TRUE)

  bridge_state <- olfml2b_install_tcgabiolinks_status_bridge(log_file = log_file)
  on.exit(olfml2b_restore_tcgabiolinks_status_bridge(bridge_state, log_file = log_file), add = TRUE)

  cached <- olfml2b_gdc_query_cache_valid(query_cache, log_file = log_file)
  if (!is.null(cached)) return(cached)

  retries <- suppressWarnings(as.integer(Sys.getenv("OLFML2B_GDC_QUERY_RETRIES", unset = "4")))
  if (!is.finite(retries) || retries < 1L) retries <- 4L
  retries <- min(retries, 8L)
  base_wait <- suppressWarnings(as.numeric(Sys.getenv("OLFML2B_GDC_RETRY_BASE_SECONDS", unset = "10")))
  if (!is.finite(base_wait) || base_wait < 0) base_wait <- 10

  profiles <- rep(c("default", "ipv4"), length.out = retries)
  error_rows <- vector("list", retries)

  for (attempt in seq_len(retries)) {
    profile <- profiles[attempt]
    if (attempt > 1L) {
      wait <- min(base_wait * (attempt - 1L), 60)
      olfml2b_log(
        "WARN",
        "Retrying GDC STAR query after ", wait, " seconds | attempt ", attempt, "/", retries,
        " | network_profile=", profile,
        log_file = log_file
      )
      Sys.sleep(wait)
    } else {
      olfml2b_log(
        "INFO",
        "GDC STAR query attempt ", attempt, "/", retries,
        " | network_profile=", profile,
        log_file = log_file
      )
    }

    started <- Sys.time()
    query_call <- function() {
      TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type = sample_types
      )
    }

    q <- tryCatch(
      {
        if (requireNamespace("httr", quietly = TRUE)) {
          cfg <- httr::config(
            connecttimeout = 60,
            timeout = 600,
            ipresolve = if (identical(profile, "ipv4")) 1L else 0L
          )
          httr::with_config(cfg, query_call())
        } else {
          query_call()
        }
      },
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

    if (!inherits(q, "error")) {
      n <- tryCatch(nrow(TCGAbiolinks::getResults(q)), error = function(e) 0L)
      if (is.finite(n) && n > 0L) {
        attr(q, "workflow_used") <- "STAR - Counts"
        attr(q, "query_network_profile") <- profile
        attr(q, "query_attempt") <- attempt
        if (!is.null(query_cache) && nzchar(query_cache)) {
          olfml2b_atomic_save_rds(q, query_cache)
        }
        olfml2b_log(
          "INFO",
          "TCGA RNA query resolved with workflow: STAR - Counts | files=", n,
          " | attempt=", attempt, " | network_profile=", profile,
          log_file = log_file
        )
        return(q)
      }
      err <- "GDCquery returned no files"
    } else {
      err <- conditionMessage(q)
    }

    error_rows[[attempt]] <- data.frame(
      attempt = attempt,
      network_profile = profile,
      elapsed_seconds = elapsed,
      error = err,
      timestamp = olfml2b_timestamp(),
      stringsAsFactors = FALSE
    )
    olfml2b_log(
      "WARN",
      "GDC STAR query attempt ", attempt, " failed | network_profile=", profile,
      " | ", err,
      log_file = log_file
    )
  }

  errors <- olfml2b_bind_rows_safe(error_rows)
  if (!is.null(query_cache) && nzchar(query_cache)) {
    error_file <- sub("\\.rds$", "_errors.csv", query_cache, ignore.case = TRUE)
    if (identical(error_file, query_cache)) error_file <- paste0(query_cache, "_errors.csv")
    olfml2b_atomic_write_csv(errors, error_file)
  }
  detail <- paste0(
    "TCGA-STAD STAR - Counts query failed after ", retries, " attempts. ",
    "The obsolete HTSeq workflow was not attempted because current GDC RNA quantification accepts STAR - Counts. ",
    "This is usually a transient GDC/TLS/proxy connectivity failure rather than a cohort or analysis error. ",
    "Rerun the same command later; a successful query will be cached automatically. Last error: ",
    if (nrow(errors)) tail(errors$error, 1L) else "unknown"
  )
  olfml2b_abort(detail)
}

olfml2b_archive_gdc_work_artifacts <- function(work_dir, chunks_dir, manifests_dir, tools_dir, log_file = NULL) {
  if (!dir.exists(work_dir)) return(invisible(data.frame()))
  files <- list.files(work_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(invisible(data.frame()))
  rows <- list()
  for (src in files) {
    bn <- basename(src)
    destination_dir <- if (grepl("\\.tar\\.gz$", bn, ignore.case = TRUE)) {
      chunks_dir
    } else if (grepl("manifest|MANIFEST", bn, ignore.case = TRUE)) {
      manifests_dir
    } else if (grepl("^gdc[-_]client|gdc-client_configuration", bn, ignore.case = TRUE)) {
      tools_dir
    } else {
      file.path(work_dir, "unclassified")
    }
    olfml2b_safe_dir_create(destination_dir)
    dest <- file.path(destination_dir, bn)
    final <- tryCatch(
      olfml2b_move_file_safe(src, dest, overwrite = FALSE),
      error = function(e) {
        olfml2b_log("WARN", "Unable to archive GDC work artifact ", src, ": ", conditionMessage(e), log_file = log_file)
        NA_character_
      }
    )
    rows[[length(rows) + 1L]] <- data.frame(
      source = src,
      destination = final,
      size_bytes = if (!is.na(final) && file.exists(final)) file.info(final)$size else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  olfml2b_bind_rows_safe(rows)
}


olfml2b_gdc_local_files_available <- function(directory, min_files = 10L) {
  if (!dir.exists(directory)) return(FALSE)
  files <- list.files(directory, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (!length(files)) return(FALSE)
  info <- file.info(files)
  files <- files[info$isdir %in% FALSE]
  files <- files[!grepl("\\.partial$|\\.tmp$|\\.log$|\\.errors\\.txt$", files, ignore.case = TRUE)]
  files <- files[file.info(files)$size > 0]
  length(files) >= min_files
}

olfml2b_gdc_download_prepare <- function(query, directory, cache_rds, overwrite = FALSE,
                                       files_per_chunk = 20L, log_file = NULL,
                                       work_dir = NULL, chunks_dir = NULL,
                                       manifests_dir = NULL, tools_dir = NULL) {
  olfml2b_require_packages(c("TCGAbiolinks", "SummarizedExperiment"), log_file)

  ## True local-first rule:
  ## 1) if prepared RDS exists, skip query/download/prepare here and read it.
  ## 2) if prepared RDS is absent but GDC raw files already exist in `directory`,
  ##    skip GDCdownload and run GDCprepare directly from local files.
  if (file.exists(cache_rds) && !overwrite) {
    olfml2b_log("INFO", "Using cached GDC prepared object; skipping GDCdownload: ", cache_rds, log_file = log_file)
    return(readRDS(cache_rds))
  }

  olfml2b_safe_dir_create(directory)
  work_dir <- olfml2b_safe_dir_create(work_dir %||% file.path(dirname(directory), "gdc_work"))
  chunks_dir <- olfml2b_safe_dir_create(chunks_dir %||% file.path(dirname(directory), "gdc_chunks"))
  manifests_dir <- olfml2b_safe_dir_create(manifests_dir %||% file.path(dirname(directory), "manifests"))
  tools_dir <- olfml2b_safe_dir_create(tools_dir %||% file.path(dirname(directory), "tools"))

  local_gdc_available <- !overwrite && olfml2b_gdc_local_files_available(directory)
  if (isTRUE(local_gdc_available)) {
    olfml2b_log(
      "INFO",
      "Local GDC files detected; skipping GDCdownload and running GDCprepare directly from: ",
      directory,
      log_file = log_file
    )
  }

  download_errors <- character()
  methods <- if (isTRUE(local_gdc_available)) "local_existing" else c("api", "client")
  prepared <- NULL

  for (method in methods) {
    ok <- tryCatch({
      if (!identical(method, "local_existing")) {
        olfml2b_with_dir(work_dir, {
          if (identical(method, "api")) {
            TCGAbiolinks::GDCdownload(
              query = query, method = method, directory = directory,
              files.per.chunk = files_per_chunk
            )
          } else {
            TCGAbiolinks::GDCdownload(
              query = query, method = method, directory = directory
            )
          }
        })
        archive <- olfml2b_archive_gdc_work_artifacts(
          work_dir = work_dir,
          chunks_dir = chunks_dir,
          manifests_dir = manifests_dir,
          tools_dir = tools_dir,
          log_file = log_file
        )
        if (nrow(archive)) {
          olfml2b_atomic_write_csv(archive, file.path(manifests_dir, paste0("GDC_work_artifact_archive_", method, ".csv")))
        }
      }
      TRUE
    }, error = function(e) {
      download_errors <<- c(download_errors, paste(method, conditionMessage(e), sep = ": "))
      FALSE
    })
    if (!ok) next

    prepared <- tryCatch(
      TCGAbiolinks::GDCprepare(query = query, directory = directory, summarizedExperiment = TRUE),
      error = function(e) {
        download_errors <<- c(download_errors, paste("GDCprepare after", method, conditionMessage(e), sep = ": "))
        NULL
      }
    )
    if (!is.null(prepared)) break
  }

  if (is.null(prepared)) {
    olfml2b_atomic_write_lines(download_errors, paste0(cache_rds, ".errors.txt"))
    olfml2b_abort("GDC local/remote prepare failed. See ", paste0(cache_rds, ".errors.txt"))
  }
  olfml2b_atomic_save_rds(prepared, cache_rds)
  prepared
}

olfml2b_tcga_assay_candidates <- function(se) {
  an <- SummarizedExperiment::assayNames(se)
  data.frame(
    assay = an,
    normalized_name = olfml2b_clean_names(an),
    stringsAsFactors = FALSE
  )
}

olfml2b_tcga_pick_assay <- function(se, type = c("counts", "tpm", "fpkm_uq", "fpkm")) {
  type <- match.arg(type)
  an <- SummarizedExperiment::assayNames(se)
  clean <- olfml2b_clean_names(an)
  patterns <- switch(type,
    counts = c("unstranded", "counts", "htseq_counts"),
    tpm = c("tpm_unstrand", "tpm_unstranded", "tpm"),
    fpkm_uq = c("fpkm_uq_unstrand", "fpkm_uq", "fpkm_uq_unstranded"),
    fpkm = c("fpkm_unstrand", "fpkm_unstranded", "fpkm")
  )
  idx <- integer()
  for (p in patterns) {
    idx <- c(idx, which(clean == p), which(grepl(p, clean, fixed = TRUE)))
  }
  idx <- unique(idx)
  if (!length(idx)) return(NULL)
  # Avoid selecting fpkm_uq when plain fpkm was requested.
  if (type == "fpkm") idx <- idx[!grepl("uq", clean[idx])]
  if (!length(idx)) return(NULL)
  SummarizedExperiment::assay(se, an[idx[1]])
}

olfml2b_tcga_row_annotation <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  if (!nrow(rd)) rd <- data.frame(row_id = rownames(se), stringsAsFactors = FALSE)
  nms <- names(rd)
  symbol_col <- olfml2b_case_when_column(rd, c("gene_name", "external_gene_name", "gene_symbol", "symbol"))
  ensembl_col <- olfml2b_case_when_column(rd, c("gene_id", "ensembl_gene_id", "ensembl"))
  entrez_col <- olfml2b_case_when_column(rd, c("entrez_gene_id", "entrezid", "entrez"))
  symbol <- if (!is.na(symbol_col)) as.character(rd[[symbol_col]]) else rep(NA_character_, nrow(rd))
  ensembl <- if (!is.na(ensembl_col)) as.character(rd[[ensembl_col]]) else rownames(se)
  ensembl <- sub("\\.[0-9]+$", "", ensembl)
  if (sum(nzchar(symbol) & !is.na(symbol)) < 0.7 * nrow(rd)) {
    olfml2b_require_packages(c("AnnotationDbi", "org.Hs.eg.db"))
    missing <- is.na(symbol) | !nzchar(symbol)
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ensembl[missing],
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    symbol[missing] <- unname(mapped)
  }
  data.frame(
    row_index = seq_len(nrow(rd)),
    ensembl_id = ensembl,
    symbol = toupper(trimws(symbol)),
    entrez_id = if (!is.na(entrez_col)) as.character(rd[[entrez_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

olfml2b_tcga_collapse_genes <- function(mat, annotation, summary = c("median", "sum")) {
  summary <- match.arg(summary)
  olfml2b_assert(nrow(mat) == nrow(annotation), "Matrix/annotation row mismatch")
  ok <- !is.na(annotation$symbol) & nzchar(annotation$symbol)
  mat <- mat[ok, , drop = FALSE]
  annotation <- annotation[ok, , drop = FALSE]
  split_idx <- split(seq_len(nrow(mat)), annotation$symbol)
  if (summary == "sum") {
    out <- vapply(split_idx, function(idx) colSums(mat[idx, , drop = FALSE], na.rm = TRUE), numeric(ncol(mat)))
  } else {
    out <- vapply(split_idx, function(idx) {
      if (length(idx) == 1L) as.numeric(mat[idx, ]) else matrixStats::colMedians(mat[idx, , drop = FALSE], na.rm = TRUE)
    }, numeric(ncol(mat)))
  }
  out <- t(out)
  colnames(out) <- colnames(mat)
  storage.mode(out) <- "numeric"
  out
}

olfml2b_tcga_sample_metadata_from_se <- function(se) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  barcode <- colnames(se)
  if (is.null(barcode)) {
    barcode_col <- olfml2b_case_when_column(cd, c("barcode", "sample", "submitter_id"))
    olfml2b_assert(!is.na(barcode_col), "Cannot identify TCGA sample barcode")
    barcode <- as.character(cd[[barcode_col]])
  }
  data.frame(
    sample_id = barcode,
    patient_id = olfml2b_tcga_patient(barcode),
    sample_type_code = olfml2b_tcga_sample_type_code(barcode),
    sample_type = ifelse(olfml2b_tcga_is_primary(barcode), "Primary Tumor",
                         ifelse(olfml2b_tcga_is_normal(barcode), "Solid Tissue Normal", "Other")),
    stringsAsFactors = FALSE
  )
}

olfml2b_tcga_select_unique_samples <- function(expr_list, sample_meta) {
  olfml2b_assert(length(expr_list) > 0L, "expr_list is empty")
  all_ids <- colnames(expr_list[[1]])
  for (m in expr_list) olfml2b_assert(identical(colnames(m), all_ids), "TCGA assay sample order mismatch")
  sample_meta <- sample_meta[match(all_ids, sample_meta$sample_id), , drop = FALSE]
  sample_meta$aliquot_priority <- ifelse(sample_meta$sample_type_code == 1L, 1L,
                                         ifelse(sample_meta$sample_type_code == 11L, 2L, 9L))
  sample_meta$vial <- substr(sample_meta$sample_id, 16L, 16L)
  ord <- order(sample_meta$patient_id, sample_meta$aliquot_priority, sample_meta$vial, sample_meta$sample_id)
  key <- paste(sample_meta$patient_id, sample_meta$sample_type_code, sep = "|")
  keep_ord <- !duplicated(key[ord])
  keep <- sort(ord[keep_ord])
  expr_list <- lapply(expr_list, function(m) m[, keep, drop = FALSE])
  sample_meta <- sample_meta[keep, , drop = FALSE]
  list(expr = expr_list, meta = sample_meta)
}

olfml2b_tcga_clinical <- function(project = "TCGA-STAD", cache_rds = NULL,
                                  overwrite = FALSE, log_file = NULL) {
  olfml2b_require_packages("TCGAbiolinks", log_file)
  if (!is.null(cache_rds) && file.exists(cache_rds) && !overwrite) {
    cached <- tryCatch(readRDS(cache_rds), error = function(e) NULL)
    if (is.list(cached) && all(c("harmonized", "raw") %in% names(cached)) &&
        is.data.frame(cached$harmonized) && nrow(cached$harmonized) > 0L) {
      olfml2b_log("INFO", "Using cached TCGA clinical object: ", cache_rds, log_file = log_file)
      return(cached)
    }
    olfml2b_log("WARN", "Ignoring unreadable/incomplete TCGA clinical cache: ", cache_rds, log_file = log_file)
  }
  clinical <- TCGAbiolinks::GDCquery_clinic(project = project, type = "clinical")
  clinical <- as.data.frame(clinical, stringsAsFactors = FALSE, check.names = FALSE)
  names(clinical) <- make.unique(olfml2b_clean_names(names(clinical)), sep = "_")
  id_col <- olfml2b_case_when_column(clinical, c("submitter_id", "bcr_patient_barcode", "case_submitter_id", "patient"))
  olfml2b_assert(!is.na(id_col), "Unable to identify patient ID in GDC clinical table")
  patient_id <- substr(as.character(clinical[[id_col]]), 1L, 12L)
  age_col <- olfml2b_case_when_column(clinical, c("age_at_diagnosis", "age"))
  sex_col <- olfml2b_case_when_column(clinical, c("gender", "sex"))
  stage_col <- olfml2b_case_when_column(clinical, c("ajcc_pathologic_stage", "tumor_stage", "pathologic_stage", "clinical_stage"))
  vital_col <- olfml2b_case_when_column(clinical, c("vital_status", "status"))
  death_col <- olfml2b_case_when_column(clinical, c("days_to_death"))
  follow_col <- olfml2b_case_when_column(clinical, c("days_to_last_follow_up", "days_to_last_known_alive", "days_to_follow_up"))
  progression_col <- olfml2b_case_when_column(clinical, c("progression_or_recurrence", "disease_progression"))
  age_raw <- if (!is.na(age_col)) olfml2b_numeric(clinical[[age_col]]) else rep(NA_real_, nrow(clinical))
  age_years <- ifelse(age_raw > 200, age_raw / 365.25, age_raw)
  death_days <- if (!is.na(death_col)) olfml2b_numeric(clinical[[death_col]]) else rep(NA_real_, nrow(clinical))
  follow_days <- if (!is.na(follow_col)) olfml2b_numeric(clinical[[follow_col]]) else rep(NA_real_, nrow(clinical))
  vital <- if (!is.na(vital_col)) olfml2b_binary_event(clinical[[vital_col]]) else ifelse(is.finite(death_days), 1L, 0L)
  os_time <- ifelse(vital == 1L & is.finite(death_days), death_days, follow_days)
  out <- data.frame(
    patient_id = patient_id,
    age = age_years,
    sex = if (!is.na(sex_col)) olfml2b_sex(clinical[[sex_col]]) else factor(NA_character_, levels = c("Female", "Male")),
    stage_raw = if (!is.na(stage_col)) as.character(clinical[[stage_col]]) else NA_character_,
    stage = if (!is.na(stage_col)) olfml2b_stage_group(clinical[[stage_col]]) else factor(NA_character_, levels = c("Stage 0","Stage I","Stage II","Stage III","Stage IV"), ordered = TRUE),
    os_time_days = os_time,
    os_event = vital,
    progression_or_recurrence = if (!is.na(progression_col)) as.character(clinical[[progression_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out$.follow_rank <- ifelse(is.finite(out$os_time_days), out$os_time_days, -Inf)
  out <- out[order(out$patient_id, -out$.follow_rank), , drop = FALSE]
  out <- out[!duplicated(out$patient_id), , drop = FALSE]
  out$.follow_rank <- NULL
  result <- list(harmonized = out, raw = clinical)
  if (!is.null(cache_rds)) olfml2b_atomic_save_rds(result, cache_rds)
  result
}

olfml2b_tcga_subtypes <- function(tumor = "stad", cache_rds = NULL,
                                  overwrite = FALSE, log_file = NULL) {
  olfml2b_require_packages("TCGAbiolinks", log_file)
  if (!is.null(cache_rds) && file.exists(cache_rds) && !overwrite) {
    cached <- tryCatch(readRDS(cache_rds), error = function(e) NULL)
    if (is.list(cached) && all(c("harmonized", "raw") %in% names(cached)) &&
        is.data.frame(cached$harmonized)) {
      olfml2b_log("INFO", "Using cached TCGA molecular subtype object: ", cache_rds, log_file = log_file)
      return(cached)
    }
    olfml2b_log("WARN", "Ignoring unreadable/incomplete TCGA subtype cache: ", cache_rds, log_file = log_file)
  }
  raw <- tryCatch(
    TCGAbiolinks::TCGAquery_subtype(tumor = tumor),
    error = function(e) {
      olfml2b_log("WARN", "TCGAquery_subtype failed: ", conditionMessage(e), log_file = log_file)
      NULL
    }
  )
  if (is.null(raw)) return(list(harmonized = data.frame(patient_id = character()), raw = data.frame()))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- make.unique(olfml2b_clean_names(names(raw)), sep = "_")
  id_col <- olfml2b_case_when_column(raw, c("patient", "bcr_patient_barcode", "barcode", "submitter_id"))
  olfml2b_assert(!is.na(id_col), "Unable to identify patient ID in TCGA subtype table")
  subtype_candidates <- grep("subtype|molecular|tcga|cluster", names(raw), ignore.case = TRUE, value = TRUE)
  score <- vapply(subtype_candidates, function(col) {
    vals <- unique(olfml2b_clean_text(raw[[col]]))
    vals <- vals[nzchar(vals)]
    n <- length(vals)
    (grepl("subtype", col, ignore.case = TRUE) * 10) + (n >= 3 && n <= 8) * 5 - abs(n - 4)
  }, numeric(1))
  subtype_col <- if (length(score)) subtype_candidates[which.max(score)] else NA_character_
  purity_candidates <- grep("purity|cpe|absolute|estimate", names(raw), ignore.case = TRUE, value = TRUE)
  purity_col <- NA_character_
  if (length(purity_candidates)) {
    pscore <- vapply(purity_candidates, function(col) sum(is.finite(olfml2b_numeric(raw[[col]]))), numeric(1))
    purity_col <- purity_candidates[which.max(pscore)]
  }
  out <- data.frame(
    patient_id = substr(as.character(raw[[id_col]]), 1L, 12L),
    molecular_subtype = if (!is.na(subtype_col)) as.character(raw[[subtype_col]]) else NA_character_,
    purity = if (!is.na(purity_col)) olfml2b_numeric(raw[[purity_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$patient_id), , drop = FALSE]
  result <- list(harmonized = out, raw = raw, subtype_column = subtype_col, purity_column = purity_col)
  if (!is.null(cache_rds)) olfml2b_atomic_save_rds(result, cache_rds)
  result
}

olfml2b_tcga_merge_metadata <- function(sample_meta, clinical, subtype) {
  meta <- merge(sample_meta, clinical, by = "patient_id", all.x = TRUE, sort = FALSE)
  if (nrow(subtype)) meta <- merge(meta, subtype, by = "patient_id", all.x = TRUE, sort = FALSE)
  meta <- meta[match(sample_meta$sample_id, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  meta
}

olfml2b_tcga_expression_qc <- function(expr, meta, target_gene = "OLFML2B") {
  olfml2b_assert(target_gene %in% rownames(expr), "Target gene missing from TCGA expression matrix")
  sample_stats <- data.frame(
    sample_id = colnames(expr),
    patient_id = meta$patient_id[match(colnames(expr), meta$sample_id)],
    sample_type = meta$sample_type[match(colnames(expr), meta$sample_id)],
    backtransformed_expression_sum = colSums(pmax(2^expr - 1, 0), na.rm = TRUE),
    median_expression = matrixStats::colMedians(expr, na.rm = TRUE),
    detected_genes = colSums(expr > 0, na.rm = TRUE),
    olfml2b_expression = as.numeric(expr[target_gene, ]),
    stringsAsFactors = FALSE
  )
  gene_stats <- data.frame(
    gene = rownames(expr),
    mean = rowMeans(expr, na.rm = TRUE),
    sd = matrixStats::rowSds(expr, na.rm = TRUE),
    detected_fraction = rowMeans(expr > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(sample = sample_stats, gene = gene_stats)
}

olfml2b_tcga_query_gene_cnv <- function(project = "TCGA-STAD", log_file = NULL) {
  olfml2b_require_packages("TCGAbiolinks", log_file)
  attempts <- list(
    list(data.category = "Copy Number Variation", data.type = "Gene Level Copy Number"),
    list(data.category = "Copy Number Variation", data.type = "Gene Level Copy Number", workflow.type = "ASCAT2"),
    list(data.category = "Copy Number Variation", data.type = "Masked Copy Number Segment")
  )
  errors <- character()
  for (args in attempts) {
    q <- tryCatch(do.call(TCGAbiolinks::GDCquery, c(list(project = project), args)), error = function(e) {
      errors <<- c(errors, conditionMessage(e)); NULL
    })
    if (!is.null(q)) {
      n <- tryCatch(nrow(TCGAbiolinks::getResults(q)), error = function(e) 0L)
      if (n > 0L) {
        attr(q, "query_args") <- args
        return(q)
      }
    }
  }
  olfml2b_warn("No GDC CNV query resolved: ", paste(errors, collapse = " | "))
  NULL
}

olfml2b_tcga_query_methylation <- function(project = "TCGA-STAD", log_file = NULL) {
  olfml2b_require_packages("TCGAbiolinks", log_file)
  attempts <- list(
    list(data.category = "DNA Methylation", data.type = "Methylation Beta Value", platform = "Illumina Human Methylation 450"),
    list(data.category = "DNA Methylation", data.type = "Methylation Beta Value"),
    list(data.category = "DNA Methylation", platform = "Illumina Human Methylation 450")
  )
  errors <- character()
  for (args in attempts) {
    q <- tryCatch(do.call(TCGAbiolinks::GDCquery, c(list(project = project), args)), error = function(e) {
      errors <<- c(errors, conditionMessage(e)); NULL
    })
    if (!is.null(q)) {
      n <- tryCatch(nrow(TCGAbiolinks::getResults(q)), error = function(e) 0L)
      if (n > 0L) {
        attr(q, "query_args") <- args
        return(q)
      }
    }
  }
  olfml2b_warn("No GDC methylation query resolved: ", paste(errors, collapse = " | "))
  NULL
}

olfml2b_extract_olfml2b_from_gene_level_cnv <- function(obj, target_gene = "OLFML2B") {
  if (inherits(obj, "SummarizedExperiment")) {
    rd <- as.data.frame(SummarizedExperiment::rowData(obj))
    sym_col <- olfml2b_case_when_column(rd, c("gene_name", "symbol", "external_gene_name"))
    if (!is.na(sym_col)) {
      idx <- which(toupper(as.character(rd[[sym_col]])) == target_gene)
      if (length(idx)) {
        mat <- SummarizedExperiment::assay(obj)
        return(data.frame(sample_id = colnames(mat), olfml2b_cnv = as.numeric(mat[idx[1], ]), stringsAsFactors = FALSE))
      }
    }
  }
  if (is.data.frame(obj)) {
    nms <- names(obj)
    gene_col <- olfml2b_case_when_column(obj, c("gene_name", "symbol", "gene"))
    sample_col <- olfml2b_case_when_column(obj, c("sample", "sample_id", "submitter_id"))
    value_col <- olfml2b_case_when_column(obj, c("copy_number", "segment_mean", "value"))
    if (!anyNA(c(gene_col, sample_col, value_col))) {
      x <- obj[toupper(as.character(obj[[gene_col]])) == target_gene, , drop = FALSE]
      return(data.frame(sample_id = as.character(x[[sample_col]]), olfml2b_cnv = olfml2b_numeric(x[[value_col]]), stringsAsFactors = FALSE))
    }
  }
  data.frame()
}



olfml2b_tcga_matrix_integrity_audit <- function(expression, counts = NULL, tpm = NULL, fpkm_uq = NULL,
                                              expression_method = NA_character_) {
  mats <- list(expression = expression, counts = counts, tpm = tpm, fpkm_uq = fpkm_uq)
  mats <- mats[!vapply(mats, is.null, logical(1))]
  summary_rows <- lapply(names(mats), function(nm) {
    x <- mats[[nm]]
    finite <- is.finite(x)
    data.frame(
      record_type = "matrix_summary",
      object = nm,
      comparison = "",
      n_genes = nrow(x),
      n_samples = ncol(x),
      finite_fraction = mean(finite),
      min = suppressWarnings(min(x[finite], na.rm = TRUE)),
      median = suppressWarnings(stats::median(x[finite], na.rm = TRUE)),
      max = suppressWarnings(max(x[finite], na.rm = TRUE)),
      identical = NA,
      max_abs_difference = NA_real_,
      expression_method = expression_method,
      stringsAsFactors = FALSE
    )
  })
  pairs <- combn(names(mats), 2L, simplify = FALSE)
  comparison_rows <- lapply(pairs, function(pair) {
    a <- mats[[pair[1]]]; b <- mats[[pair[2]]]
    same_dim <- identical(dim(a), dim(b))
    same_names <- identical(rownames(a), rownames(b)) && identical(colnames(a), colnames(b))
    nr <- min(nrow(a), 250L); nc <- min(ncol(a), 100L)
    ri <- unique(round(seq(1, nrow(a), length.out = nr)))
    ci <- unique(round(seq(1, ncol(a), length.out = nc)))
    diff <- if (same_dim) suppressWarnings(max(abs(a[ri, ci, drop = FALSE] - b[ri, ci, drop = FALSE]), na.rm = TRUE)) else NA_real_
    data.frame(
      record_type = "pairwise_comparison",
      object = NA_character_,
      comparison = paste(pair, collapse = " vs "),
      n_genes = NA_integer_, n_samples = NA_integer_, finite_fraction = NA_real_,
      min = NA_real_, median = NA_real_, max = NA_real_,
      identical = same_dim && same_names && identical(a[ri, ci, drop = FALSE], b[ri, ci, drop = FALSE]),
      max_abs_difference = diff,
      expression_method = expression_method,
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, c(summary_rows, comparison_rows))
  if (!is.null(tpm) && identical(dim(expression), dim(tpm))) {
    nr <- min(nrow(tpm), 500L); nc <- min(ncol(tpm), 150L)
    ri <- unique(round(seq(1, nrow(tpm), length.out = nr)))
    ci <- unique(round(seq(1, ncol(tpm), length.out = nc)))
    expected <- log2(pmax(tpm[ri, ci, drop = FALSE], 0) + 1)
    observed <- expression[ri, ci, drop = FALSE]
    audit <- rbind(audit, data.frame(
      record_type = "expression_contract",
      object = "expression",
      comparison = "expression vs log2(TPM+1)",
      n_genes = length(ri), n_samples = length(ci), finite_fraction = mean(is.finite(observed)),
      min = NA_real_, median = NA_real_, max = NA_real_, identical = NA,
      max_abs_difference = suppressWarnings(max(abs(observed - expected), na.rm = TRUE)),
      expression_method = expression_method,
      stringsAsFactors = FALSE
    ))
  }
  audit
}

# ==============================================================================
# Part1 - TCGA-STAD official GDC download and processing
# Outputs a fully harmonized tumor/normal gene-expression and clinical object.
# ==============================================================================

run_part1 <- function(ctx = NULL) {
  if (is.null(ctx)) {
    code_root <- olfml2b_find_code_root()
    ctx <- olfml2b_load_context()
  }
  ctx$part_paths <- olfml2b_part_paths(ctx, "Part1")
  olfml2b_require_packages(c(
    "TCGAbiolinks", "SummarizedExperiment", "S4Vectors", "matrixStats",
    "edgeR", "ggplot2", "data.table", "AnnotationDbi", "org.Hs.eg.db"
  ))
  log_file <- file.path(ctx$part_paths$logs, "Part1_TCGA_STAD_Acquire_Process.log")
  error_file <- file.path(ctx$part_paths$reports, "pipeline_errors.csv")
  out_rds <- file.path(ctx$dirs$objects, "Part1_TCGA_STAD.rds")
  if (file.exists(out_rds) && !ctx$overwrite_results) {
    olfml2b_log("INFO", "Part1 output exists; loading: ", out_rds, log_file = log_file)
    return(readRDS(out_rds))
  }
  olfml2b_log("INFO", "Starting Part1 | ", ctx$version, log_file = log_file)
  se_cache <- file.path(ctx$dirs$cache_tcga, "TCGA_STAD_RNA_SummarizedExperiment.rds")
  query_cache <- file.path(ctx$dirs$cache_tcga, "TCGA_STAD_STAR_Counts_query.rds")
  clinical_cache <- file.path(ctx$dirs$cache_tcga, "TCGA_STAD_clinical_harmonized.rds")
  subtype_cache <- file.path(ctx$dirs$cache_tcga, "TCGA_STAD_molecular_subtypes.rds")
  query_report_file <- file.path(ctx$part_paths$reports, "Part1_TCGA_RNA_query_results.csv")
  prepared_cache_available <- file.exists(se_cache) && !ctx$overwrite_downloads

  if (prepared_cache_available) {
    olfml2b_log(
      "INFO",
      "Prepared TCGA RNA cache is available; skipping GDC status and query network calls: ", se_cache,
      log_file = log_file
    )
    gdc_status <- data.frame(
      endpoint = "https://api.gdc.cancer.gov/status",
      network_profile = "not_checked",
      reachable = NA,
      http_status = NA_integer_,
      api_status = "skipped_prepared_cache_available",
      data_release = NA_character_,
      elapsed_seconds = 0,
      error = NA_character_,
      timestamp = olfml2b_timestamp(),
      stringsAsFactors = FALSE
    )
  } else {
    gdc_status <- olfml2b_run_stage(
      "GDC status",
      olfml2b_gdc_status(log_file),
      log_file, error_file, required = FALSE
    )
  }
  if (is.data.frame(gdc_status)) {
    olfml2b_atomic_write_csv(gdc_status, file.path(ctx$part_paths$reports, "Part1_GDC_status.csv"))
  }

  if (prepared_cache_available) {
    query <- NULL
    workflow_used <- "STAR - Counts (prepared cache)"
    query_results <- if (file.exists(query_report_file)) {
      tryCatch(utils::read.csv(query_report_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
    } else {
      data.frame()
    }
  } else {
    query <- olfml2b_run_stage(
      "GDC RNA query",
      olfml2b_tcga_query_rna(
        ctx$contract$tcga_project,
        query_cache = query_cache,
        log_file = log_file
      ),
      log_file, error_file, required = TRUE
    )
    query_results <- as.data.frame(TCGAbiolinks::getResults(query))
    olfml2b_atomic_write_csv(query_results, query_report_file)
    workflow_used <- attr(query, "workflow_used") %||% "STAR - Counts"
  }

  legacy_gdc_project_dir <- file.path(ctx$dirs$raw_tcga, ctx$contract$tcga_project)
  gdc_data_dir <- if (dir.exists(legacy_gdc_project_dir)) ctx$dirs$raw_tcga else ctx$dirs$raw_tcga_gdc
  olfml2b_log("INFO", "Resolved GDC data directory: ", gdc_data_dir, log_file = log_file)
  se <- olfml2b_run_stage(
    "GDC RNA download and prepare",
    olfml2b_gdc_download_prepare(
      query = query,
      directory = gdc_data_dir,
      cache_rds = se_cache,
      overwrite = ctx$overwrite_downloads,
      files_per_chunk = 20L,
      log_file = log_file,
      work_dir = file.path(ctx$dirs$temp_gdc, "Part1_RNA"),
      chunks_dir = ctx$dirs$raw_tcga_gdc_chunks,
      manifests_dir = ctx$dirs$raw_tcga_manifests,
      tools_dir = ctx$dirs$raw_tcga_tools
    ),
    log_file, error_file, required = TRUE
  )
  assay_audit <- olfml2b_tcga_assay_candidates(se)
  olfml2b_atomic_write_csv(assay_audit, file.path(ctx$part_paths$reports, "Part1_TCGA_assay_audit.csv"))

  annotation <- olfml2b_run_stage(
    "TCGA row annotation",
    olfml2b_tcga_row_annotation(se),
    log_file, error_file, required = TRUE
  )
  olfml2b_atomic_write_csv(annotation, file.path(ctx$part_paths$reports, "Part1_TCGA_gene_annotation.csv"))

  counts_raw <- olfml2b_tcga_pick_assay(se, "counts")
  tpm_raw <- olfml2b_tcga_pick_assay(se, "tpm")
  fpkm_uq_raw <- olfml2b_tcga_pick_assay(se, "fpkm_uq")
  olfml2b_assert(!is.null(counts_raw) || !is.null(tpm_raw), "Neither counts nor TPM assay was available in TCGA object")
  if (!is.null(counts_raw)) {
    counts_raw <- as.matrix(counts_raw)
    storage.mode(counts_raw) <- "numeric"
  }
  if (!is.null(tpm_raw)) {
    tpm_raw <- as.matrix(tpm_raw)
    storage.mode(tpm_raw) <- "numeric"
  }
  if (!is.null(fpkm_uq_raw)) {
    fpkm_uq_raw <- as.matrix(fpkm_uq_raw)
    storage.mode(fpkm_uq_raw) <- "numeric"
  }

  counts_gene <- if (!is.null(counts_raw)) olfml2b_tcga_collapse_genes(counts_raw, annotation, summary = "sum") else NULL
  tpm_gene <- if (!is.null(tpm_raw)) olfml2b_tcga_collapse_genes(tpm_raw, annotation, summary = "median") else NULL
  fpkm_uq_gene <- if (!is.null(fpkm_uq_raw)) olfml2b_tcga_collapse_genes(fpkm_uq_raw, annotation, summary = "median") else NULL
  if (!is.null(tpm_gene)) {
    expression <- log2(pmax(tpm_gene, 0) + 1)
    expression_method <- "log2(TPM+1)_from_GDC_STAR"
  } else {
    olfml2b_assert(!is.null(counts_gene), "Counts are required when TPM is absent")
    dge <- edgeR::DGEList(counts = round(pmax(counts_gene, 0)))
    keep <- edgeR::filterByExpr(dge, group = ifelse(olfml2b_tcga_is_primary(colnames(counts_gene)), "Tumor", "Normal"))
    dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
    expression <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    expression_method <- "edgeR_TMM_log2CPM_due_to_missing_TPM"
  }
  olfml2b_assert("OLFML2B" %in% rownames(expression), "OLFML2B is absent after TCGA gene aggregation")

  sample_meta <- olfml2b_tcga_sample_metadata_from_se(se)
  expr_list <- list(expression = expression)
  if (!is.null(counts_gene)) expr_list$counts <- counts_gene
  if (!is.null(tpm_gene)) expr_list$tpm <- tpm_gene
  if (!is.null(fpkm_uq_gene)) expr_list$fpkm_uq <- fpkm_uq_gene
  selected <- olfml2b_tcga_select_unique_samples(expr_list, sample_meta)
  expr_list <- selected$expr
  sample_meta <- selected$meta
  expression <- expr_list$expression
  matrix_integrity_audit <- olfml2b_tcga_matrix_integrity_audit(
    expression = expression,
    counts = expr_list$counts %||% NULL,
    tpm = expr_list$tpm %||% NULL,
    fpkm_uq = expr_list$fpkm_uq %||% NULL,
    expression_method = expression_method
  )
  olfml2b_atomic_write_csv(
    matrix_integrity_audit,
    file.path(ctx$part_paths$reports, "Part1_TCGA_matrix_integrity_audit.csv")
  )
  contract_idx <- which(
    !is.na(matrix_integrity_audit$comparison) &
      matrix_integrity_audit$comparison == "expression vs log2(TPM+1)"
  )
  olfml2b_assert(
    length(contract_idx) == 1L,
    "TCGA expression-contract audit did not produce exactly one contract row; rows=",
    length(contract_idx)
  )
  contract_row <- matrix_integrity_audit[contract_idx, , drop = FALSE]
  contract_diff <- suppressWarnings(as.numeric(contract_row$max_abs_difference[1L]))
  olfml2b_assert(
    is.finite(contract_diff) && contract_diff < 1e-8,
    "TCGA primary expression matrix is not log2(TPM+1); max sampled difference=",
    if (is.finite(contract_diff)) signif(contract_diff, 6) else as.character(contract_diff),
    ". See Part1_TCGA_matrix_integrity_audit.csv for finite-value and assay diagnostics."
  )

  clinical_obj <- olfml2b_run_stage(
    "GDC clinical harmonisation",
    olfml2b_tcga_clinical(ctx$contract$tcga_project, clinical_cache, ctx$overwrite_downloads, log_file),
    log_file, error_file, required = TRUE
  )
  subtype_obj <- olfml2b_run_stage(
    "TCGA molecular subtype retrieval",
    olfml2b_tcga_subtypes("stad", subtype_cache, ctx$overwrite_downloads, log_file),
    log_file, error_file, required = FALSE
  )
  subtype_harmonized <- if (olfml2b_is_failure(subtype_obj)) data.frame(patient_id = character()) else subtype_obj$harmonized
  metadata <- olfml2b_tcga_merge_metadata(sample_meta, clinical_obj$harmonized, subtype_harmonized)
  metadata$olfml2b_expression <- as.numeric(expression["OLFML2B", metadata$sample_id])
  metadata$olfml2b_z <- olfml2b_z(metadata$olfml2b_expression)
  metadata$is_tumor <- metadata$sample_type == "Primary Tumor"
  metadata$is_normal <- metadata$sample_type == "Solid Tissue Normal"
  metadata$paired_available <- metadata$patient_id %in% intersect(
    metadata$patient_id[metadata$is_tumor], metadata$patient_id[metadata$is_normal]
  )
  rownames(metadata) <- metadata$sample_id

  # Reusable derived inputs are separated from raw downloads and final results.
  olfml2b_atomic_save_rds(expression, file.path(ctx$dirs$derived_expression, "TCGA_STAD_log_expression_gene_symbol.rds"))
  if (!is.null(expr_list$counts)) {
    olfml2b_atomic_save_rds(expr_list$counts, file.path(ctx$dirs$derived_expression, "TCGA_STAD_gene_counts.rds"), compress = "xz")
  }
  olfml2b_atomic_write_csv(metadata, file.path(ctx$dirs$derived_clinical, "TCGA_STAD_sample_metadata_harmonized.csv"))

  qc <- olfml2b_tcga_expression_qc(expression, metadata, "OLFML2B")
  olfml2b_atomic_write_csv(metadata, file.path(ctx$part_paths$tables, "Part1_TCGA_sample_metadata.csv"))
  olfml2b_atomic_write_csv(clinical_obj$raw, file.path(ctx$part_paths$reports, "Part1_TCGA_clinical_raw.csv"))
  olfml2b_atomic_write_csv(clinical_obj$harmonized, file.path(ctx$part_paths$reports, "Part1_TCGA_clinical_harmonized.csv"))
  if (!olfml2b_is_failure(subtype_obj)) {
    olfml2b_atomic_write_csv(subtype_obj$raw, file.path(ctx$part_paths$reports, "Part1_TCGA_subtypes_raw.csv"))
    olfml2b_atomic_write_csv(subtype_obj$harmonized, file.path(ctx$part_paths$reports, "Part1_TCGA_subtypes_harmonized.csv"))
  }
  olfml2b_atomic_write_csv(qc$sample, file.path(ctx$part_paths$tables, "Part1_TCGA_sample_QC.csv"))
  olfml2b_atomic_write_csv(qc$gene, file.path(ctx$part_paths$tables, "Part1_TCGA_gene_QC.csv"))

  sample_flow <- olfml2b_make_sample_flow("TCGA_STAD", list(
    GDC_files_queried = if (nrow(query_results)) nrow(query_results) else NA_integer_,
    prepared_aliquots = ncol(se),
    unique_patient_sampletype_aliquots = ncol(expression),
    primary_tumors = sum(metadata$is_tumor),
    solid_tissue_normals = sum(metadata$is_normal),
    paired_patients = length(unique(metadata$patient_id[metadata$paired_available]))
  ))
  olfml2b_atomic_write_csv(sample_flow, file.path(ctx$part_paths$tables, "Part1_TCGA_sample_flow.csv"))

  cohort_audit <- data.frame(
    cohort = "TCGA_STAD",
    workflow_used = workflow_used,
    expression_method = expression_method,
    n_genes = nrow(expression),
    n_samples = ncol(expression),
    n_tumor = sum(metadata$is_tumor),
    n_normal = sum(metadata$is_normal),
    n_pairs = length(intersect(metadata$patient_id[metadata$is_tumor], metadata$patient_id[metadata$is_normal])),
    n_os_complete_tumor = sum(metadata$is_tumor & is.finite(metadata$os_time_days) & metadata$os_event %in% c(0, 1)),
    os_events_tumor = sum(metadata$is_tumor & metadata$os_event == 1L, na.rm = TRUE),
    n_stage_tumor = sum(metadata$is_tumor & !is.na(metadata$stage)),
    n_subtype_tumor = sum(metadata$is_tumor & !is.na(metadata$molecular_subtype) & nzchar(olfml2b_clean_text(metadata$molecular_subtype))),
    n_purity_tumor = sum(metadata$is_tumor & is.finite(metadata$purity)),
    olfml2b_mean_tumor = mean(metadata$olfml2b_expression[metadata$is_tumor], na.rm = TRUE),
    olfml2b_mean_normal = mean(metadata$olfml2b_expression[metadata$is_normal], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  olfml2b_atomic_write_csv(cohort_audit, file.path(ctx$part_paths$tables, "Part1_TCGA_cohort_audit.csv"))

  # Core QC figures
  plot_df <- metadata
  p_olfml2b <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_type, y = olfml2b_expression)) +
    ggplot2::geom_violin(trim = FALSE, fill = "grey90", color = "grey35") +
    ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
    ggplot2::labs(title = "TCGA-STAD OLFML2B expression", subtitle = expression_method,
                  x = NULL, y = "OLFML2B expression") + olfml2b_base_theme()
  olfml2b_save_plot(p_olfml2b, file.path(ctx$part_paths$figures, "Part1_TCGA_OLFML2B_tumor_normal.png"), 6.5, 5.2)

  sample_subset <- qc$sample[order(qc$sample$detected_genes), ]
  p_detect <- ggplot2::ggplot(sample_subset, ggplot2::aes(x = seq_along(detected_genes), y = detected_genes, shape = sample_type)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = "TCGA-STAD detected-gene distribution", x = "Samples ordered by detected genes", y = "Detected genes") +
    olfml2b_base_theme()
  olfml2b_save_plot(p_detect, file.path(ctx$part_paths$figures, "Part1_TCGA_detected_genes_QC.png"), 7.5, 5.2)

  p_density <- ggplot2::ggplot(metadata, ggplot2::aes(x = olfml2b_expression, fill = sample_type)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::labs(title = "TCGA-STAD OLFML2B expression density", x = "OLFML2B expression", y = "Density", fill = "Sample type") +
    olfml2b_base_theme()
  olfml2b_save_plot(p_density, file.path(ctx$part_paths$figures, "Part1_TCGA_OLFML2B_density.png"), 7, 5)

  result <- list(
    version = ctx$version,
    generated_at = olfml2b_timestamp(),
    cohort = "TCGA_STAD",
    project = ctx$contract$tcga_project,
    workflow_used = workflow_used,
    expression_method = expression_method,
    expression = expression,
    counts = expr_list$counts %||% NULL,
    tpm = expr_list$tpm %||% NULL,
    fpkm_uq = expr_list$fpkm_uq %||% NULL,
    gene_annotation = annotation,
    sample_metadata = metadata,
    clinical_raw = clinical_obj$raw,
    subtype_raw = if (!olfml2b_is_failure(subtype_obj)) subtype_obj$raw else NULL,
    qc = qc,
    matrix_integrity_audit = matrix_integrity_audit,
    sample_flow = sample_flow,
    cohort_audit = cohort_audit,
    query_results = as.data.frame(query_results),
    source_files_manifest = olfml2b_file_manifest(ctx$dirs$raw_tcga)
  )
  olfml2b_atomic_save_rds(result, out_rds)
  olfml2b_atomic_write_csv(result$source_files_manifest, file.path(ctx$part_paths$reports, "Part1_TCGA_raw_file_manifest.csv"))
  olfml2b_capture_session(file.path(ctx$part_paths$reports, "Part1_sessionInfo.txt"))
  olfml2b_log("INFO", "Part1 complete: ", out_rds, log_file = log_file)
  invisible(result)
}


# ==============================================================================
# OLFML2B explicit Part1 entry
# ==============================================================================

run_olfml2b_part1_tcga <- function(ctx = NULL) {
  if (is.null(ctx)) ctx <- olfml2b_load_context()
  ctx$contract <- olfml2b_default_contract()
  res <- run_part1(ctx)
  # Add explicit fact-checked target identifiers to audit output.
  pdir <- olfml2b_part_paths(ctx, "Part1")
  audit_path <- file.path(pdir$tables, "Part1_TCGA_cohort_audit.csv")
  if (file.exists(audit_path)) {
    audit <- utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
    audit$target_ensembl <- ctx$contract$target_ensembl
    audit$target_entrez <- ctx$contract$target_entrez
    audit$target_hgnc <- ctx$contract$target_hgnc
    audit$target_locus_grch38 <- paste0(ctx$contract$target_chromosome, ":", ctx$contract$target_grch38_start, "-", ctx$contract$target_grch38_end, ":", ctx$contract$target_strand)
    olfml2b_atomic_write_csv(audit, audit_path)
  }
  invisible(res)
}
