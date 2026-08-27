# ==============================================================================
# OLFML2B-STAD adapted compact engine
# File 00/06: core utilities, embedded contracts, statistics and Part0.
# ==============================================================================
.compact_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.compact_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.compact_valid_root <- function(path) {
  nzchar(path) && dir.exists(path) &&
    file.exists(file.path(path, "00_OLFML2B_PART0_CONFIG_CORE.R")) &&
    file.exists(file.path(path, "03_OLFML2B_PART3_BULK_SURVIVAL.R"))
}
.compact_code_root <- if (.compact_valid_root(.compact_env_root)) {
  # The launcher has already resolved and hash-verified this directory.  Keep
  # it: nested sys.source() on Windows does not reliably expose $ofile.
  normalizePath(.compact_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.compact_entry) && file.exists(.compact_entry)) {
  dirname(normalizePath(.compact_entry, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .compact_code_root)

# ==============================================================================
# OLFML2B-STAD Part0-5 production release
# Core utilities: paths, logging, atomic IO, package bootstrap, parsing, auditing
# Self-contained production code. No user-supplied config/manifest is required.
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

olfml2b_abort <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

olfml2b_warn <- function(...) {
  warning(paste0(..., collapse = ""), call. = FALSE, immediate. = TRUE)
}

olfml2b_assert <- function(condition, ...) {
  if (!isTRUE(condition)) olfml2b_abort(...)
  invisible(TRUE)
}

olfml2b_timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

olfml2b_safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    olfml2b_assert(ok || dir.exists(path), "Cannot create directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

olfml2b_find_code_root <- function(start = getwd()) {
  env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
  if (!nzchar(env_root)) env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
  required <- c("00_OLFML2B_PART0_CONFIG_CORE.R", "03_OLFML2B_PART3_BULK_SURVIVAL.R")
  valid_root <- function(path) dir.exists(path) && all(file.exists(file.path(path, required)))
  if (nzchar(env_root) && valid_root(env_root)) {
    root <- normalizePath(env_root, winslash = "/", mustWork = TRUE)
    Sys.setenv(OLFML2B_STAD_CODE_ROOT = root, OLFML2B_STAD_CODE_ROOT = root)
    return(root)
  }
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)
  candidates <- unique(c(current, file.path(current, "R"), dirname(current), file.path(dirname(current), "R")))
  for (cand in candidates) {
    if (valid_root(cand)) {
      root <- normalizePath(cand, winslash = "/", mustWork = TRUE)
      Sys.setenv(OLFML2B_STAD_CODE_ROOT = root, OLFML2B_STAD_CODE_ROOT = root)
      return(root)
    }
  }
  for (i in 0:12) {
    if (valid_root(current)) {
      root <- normalizePath(current, winslash = "/", mustWork = TRUE)
      Sys.setenv(OLFML2B_STAD_CODE_ROOT = root, OLFML2B_STAD_CODE_ROOT = root)
      return(root)
    }
    rdir <- file.path(current, "R")
    if (valid_root(rdir)) {
      root <- normalizePath(rdir, winslash = "/", mustWork = TRUE)
      Sys.setenv(OLFML2B_STAD_CODE_ROOT = root, OLFML2B_STAD_CODE_ROOT = root)
      return(root)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  olfml2b_abort("Cannot locate OLFML2B code root. Put files 00-05 in <project>/R or define OLFML2B_STAD_CODE_ROOT.")
}


olfml2b_resolve_project_root <- function(code_root = NULL) {
  env_root <- Sys.getenv("OLFML2B_STAD_ROOT", unset = "")
  if (!nzchar(env_root)) env_root <- Sys.getenv("OLFML2B_STAD_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  if (is.null(code_root)) {
    code_root <- tryCatch(olfml2b_find_code_root(), error = function(e) getwd())
  }
  code_root <- normalizePath(code_root, winslash = "/", mustWork = FALSE)
  # Recommended layout: <project_root>/R contains all seven source files.
  # When the code folder is literally named R, infer the project root as its parent.
  if (identical(tolower(basename(code_root)), "r")) {
    return(normalizePath(dirname(code_root), winslash = "/", mustWork = FALSE))
  }
  normalizePath(code_root, winslash = "/", mustWork = FALSE)
}

olfml2b_with_dir <- function(path, expr) {
  path <- olfml2b_safe_dir_create(path)
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  force(expr)
}

olfml2b_part_paths <- function(ctx, part) {
  part <- as.character(part)[1L]
  olfml2b_assert(grepl("^Part[0-9]$", part), "Invalid part identifier: ", part)
  list(
    tables = olfml2b_safe_dir_create(file.path(ctx$dirs$tables_root, part)),
    figures = olfml2b_safe_dir_create(file.path(ctx$dirs$figures_root, part)),
    reports = olfml2b_safe_dir_create(file.path(ctx$dirs$reports_root, part)),
    logs = olfml2b_safe_dir_create(file.path(ctx$dirs$logs_runtime, part)),
    qc = olfml2b_safe_dir_create(file.path(ctx$dirs$qc_root, part))
  )
}


# ------------------------------------------------------------------------------
# Publication compatibility interface for downstream Part3-Part5
# ------------------------------------------------------------------------------
# The core utilities predate the target-specific Part3-Part5 evidence layers.
# These downstream layers expect a stable directory contract named
# olfml2b_build_dirs(). Keep original Part0 semantics and add this missing public
# interface here, instead of patching startup scripts or downstream analysis code.
olfml2b_build_dirs <- function(root = olfml2b_resolve_project_root()) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  safe_dir <- function(path) olfml2b_safe_dir_create(file.path(path))

  project_data_root <- safe_dir(file.path(root, "data"))
  shared_data_env <- Sys.getenv("OLFML2B_STAD_SHARED_DATA_ROOT", unset = "")
  shared_data_ok <- nzchar(shared_data_env) && dir.exists(shared_data_env)
  source_data_root <- if (shared_data_ok) {
    normalizePath(shared_data_env, winslash = "/", mustWork = TRUE)
  } else {
    project_data_root
  }

  # Raw downloads and immutable caches may be shared with the prior STAD project.
  # Target-derived objects and all outputs remain isolated under OLFML2B_STAD.
  raw_root <- safe_dir(file.path(source_data_root, "raw"))
  cache_root <- safe_dir(file.path(source_data_root, "cache"))
  derived_root <- safe_dir(file.path(project_data_root, "derived"))

  output_root <- safe_dir(file.path(root, "output"))
  tables_root <- safe_dir(file.path(output_root, "tables"))
  figures_root <- safe_dir(file.path(output_root, "figures"))
  reports_root <- safe_dir(file.path(output_root, "reports"))
  qc_root <- safe_dir(file.path(output_root, "qc"))
  objects_root <- safe_dir(file.path(output_root, "objects"))
  logs_root <- safe_dir(file.path(root, "logs"))
  logs_runtime <- safe_dir(file.path(logs_root, "runtime"))
  logs_install <- safe_dir(file.path(logs_root, "install"))
  temp_root <- safe_dir(file.path(root, "temp"))

  list(
    root = root,
    code_root = normalizePath(Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = file.path(root, "R")), winslash = "/", mustWork = FALSE),
    data = project_data_root,
    source_data_root = source_data_root,
    shared_data_enabled = shared_data_ok,
    raw = raw_root,
    cache = cache_root,
    derived = derived_root,
    output = output_root,
    tables_root = tables_root,
    figures_root = figures_root,
    reports_root = reports_root,
    qc_root = qc_root,
    objects = objects_root,
    logs = logs_root,
    logs_runtime = logs_runtime,
    logs_install = logs_install,
    raw_tcga = safe_dir(file.path(raw_root, "TCGA_STAD")),
    raw_tcga_gdc = safe_dir(file.path(raw_root, "TCGA_STAD", "GDC")),
    raw_tcga_gdc_chunks = safe_dir(file.path(raw_root, "TCGA_STAD", "GDC_chunks")),
    raw_tcga_manifests = safe_dir(file.path(raw_root, "TCGA_STAD", "manifests")),
    raw_tcga_tools = safe_dir(file.path(raw_root, "TCGA_STAD", "tools")),
    raw_geo = safe_dir(file.path(raw_root, "GEO")),
    cache_tcga = safe_dir(file.path(cache_root, "TCGA_STAD")),
    cache_geo = safe_dir(file.path(cache_root, "GEO")),
    temp = temp_root,
    temp_gdc = safe_dir(file.path(temp_root, "GDC")),
    derived_expression = safe_dir(file.path(derived_root, "expression")),
    derived_clinical = safe_dir(file.path(derived_root, "clinical"))
  )
}

# Stable aliases used by repaired visualization layers. Part3 may define legacy
# wrappers later, but Part0 should expose these common interfaces first.
if (!exists("olfml2b_write_csv", mode = "function")) {
  olfml2b_write_csv <- function(x, path) olfml2b_atomic_write_csv(x, path)
}
if (!exists("olfml2b_bind_rows", mode = "function")) {
  olfml2b_bind_rows <- function(xs) olfml2b_bind_rows_safe(xs)
}

olfml2b_code_root_generated_patterns <- function() {
  c(
    "^\\.Rhistory$", "^\\.RData$", "^df\\.rds$", "^results\\.rds$",
    "^PACKAGE_INSTALL\\.log$", "^PACKAGE_INSTALL_STATUS\\.csv$",
    "^PACKAGE_REQUESTED\\.csv$", "^Fri_.*\\.tar\\.gz$",
    "^gdc[-_]manifest.*\\.txt$", "^MANIFEST\\.txt$",
    "^gdc-client.*\\.(zip|tar\\.gz|exe)$", "^gdc-client_configuration\\..*$"
  )
}

olfml2b_audit_code_root_cleanliness <- function(code_root) {
  files <- list.files(code_root, all.files = TRUE, no.. = TRUE, full.names = FALSE)
  patterns <- olfml2b_code_root_generated_patterns()
  contaminated <- files[vapply(files, function(x) any(vapply(patterns, grepl, logical(1), x = x, perl = TRUE)), logical(1))]
  data.frame(
    file = contaminated,
    issue = rep("runtime_generated_file_in_code_root", length(contaminated)),
    stringsAsFactors = FALSE
  )
}

olfml2b_move_file_safe <- function(source, destination, overwrite = FALSE) {
  olfml2b_assert(file.exists(source), "Source file does not exist: ", source)
  olfml2b_safe_dir_create(dirname(destination))
  if (file.exists(destination) && !overwrite) {
    stem <- tools::file_path_sans_ext(basename(destination))
    ext <- tools::file_ext(destination)
    suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
    destination <- file.path(dirname(destination), paste0(stem, "_", suffix, if (nzchar(ext)) paste0(".", ext) else ""))
  }
  ok <- file.rename(source, destination)
  if (!ok) {
    ok <- file.copy(source, destination, overwrite = overwrite)
    if (ok) unlink(source, recursive = TRUE, force = TRUE)
  }
  olfml2b_assert(ok && file.exists(destination), "Failed to move ", source, " to ", destination)
  normalizePath(destination, winslash = "/", mustWork = TRUE)
}

olfml2b_log <- function(level = "INFO", ..., log_file = NULL, echo = TRUE) {
  line <- sprintf("[%s] [%s] %s", olfml2b_timestamp(), toupper(level), paste0(..., collapse = ""))
  if (isTRUE(echo)) message(line)
  if (!is.null(log_file)) {
    olfml2b_safe_dir_create(dirname(log_file))
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
  invisible(line)
}

olfml2b_with_log <- function(expr, stage, log_file) {
  started <- Sys.time()
  olfml2b_log("INFO", "START ", stage, log_file = log_file)
  ans <- tryCatch(
    force(expr),
    error = function(e) {
      olfml2b_log("ERROR", stage, " failed: ", conditionMessage(e), log_file = log_file)
      stop(e)
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  olfml2b_log("INFO", sprintf("END %s | %.1f sec", stage, elapsed), log_file = log_file)
  ans
}

olfml2b_atomic_save_rds <- function(object, path, compress = "xz") {
  olfml2b_safe_dir_create(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  olfml2b_assert(file.exists(tmp) && file.info(tmp)$size > 0, "Failed to write temporary RDS: ", tmp)
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  olfml2b_assert(ok && file.exists(path), "Failed to commit RDS: ", path)
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

olfml2b_atomic_write_lines <- function(lines, path, use_bytes = TRUE) {
  olfml2b_safe_dir_create(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  writeLines(enc2utf8(as.character(lines)), tmp, useBytes = use_bytes)
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  olfml2b_assert(ok && file.exists(path), "Failed to write text file: ", path)
  invisible(path)
}



olfml2b_bind_rows_safe <- function(xs, template = NULL) {
  if (is.null(xs)) xs <- list()
  if (!is.list(xs)) xs <- list(xs)
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) {
    if (is.data.frame(template)) return(template[0, , drop = FALSE])
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  aligned <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

olfml2b_atomic_write_csv <- function(x, path, row.names = FALSE, na = "") {
  olfml2b_safe_dir_create(dirname(path))
  # Empty optional outputs are materialised as valid zero-row CSV files rather
  # than failing on NULL. Lists must be converted explicitly by the caller.
  if (is.null(x)) x <- data.frame()
  if (is.atomic(x) && is.null(dim(x))) x <- data.frame(value = x, stringsAsFactors = FALSE)
  olfml2b_assert(is.data.frame(x) || is.matrix(x), "CSV output must be a data.frame or matrix: ", path)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  utils::write.csv(x, tmp, row.names = row.names, na = na, fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  olfml2b_assert(ok && file.exists(path), "Failed to write CSV: ", path)
  invisible(path)
}

olfml2b_write_tsv <- function(x, path, row.names = FALSE, na = "") {
  olfml2b_safe_dir_create(dirname(path))
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = row.names,
    col.names = TRUE, na = na, fileEncoding = "UTF-8"
  )
  invisible(path)
}

olfml2b_read_csv <- function(path, ...) {
  olfml2b_assert(file.exists(path), "Missing CSV: ", path)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE, na.strings = c("", "NA", "N/A", "null"), ...))
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

olfml2b_file_manifest <- function(root, recursive = TRUE) {
  if (!dir.exists(root)) return(data.frame())
  files <- list.files(root, recursive = recursive, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  # Interrupted downloads are not scientific inputs and must never appear in
  # a frozen input manifest as if they were complete source files.
  files <- files[!grepl("\\.partial$|\\.part$|\\.tmp$", files, ignore.case = TRUE)]
  if (!length(files)) return(data.frame())
  info <- file.info(files)
  data.frame(
    relative_path = substring(normalizePath(files, winslash = "/", mustWork = TRUE),
                              nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L),
    size_bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
}

olfml2b_hash_object <- function(x) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = FALSE)
  unname(tools::md5sum(tmp))
}

olfml2b_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

olfml2b_require_packages <- function(packages, log_file = NULL) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    olfml2b_log("ERROR", "Missing packages: ", paste(missing, collapse = ", "), log_file = log_file)
    olfml2b_abort("Missing required packages. Run INSTALL_PACKAGES.R first: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

olfml2b_repo_candidates <- function() {
  # Official CRAN first, then two widely used mirrors as network fallbacks.
  c(
    "https://cloud.r-project.org",
    "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
    "https://mirrors.ustc.edu.cn/CRAN"
  )
}

olfml2b_probe_repository <- function(repo, timeout_sec = 45L) {
  probe <- paste0(sub("/$", "", repo), "/src/contrib/PACKAGES.gz")
  dest <- tempfile(fileext = ".gz")
  on.exit(unlink(dest, force = TRUE), add = TRUE)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(as.numeric(old_timeout), timeout_sec))
  ok <- tryCatch({
    suppressWarnings(utils::download.file(
      probe, destfile = dest, mode = "wb", quiet = TRUE,
      method = getOption("download.file.method", "libcurl")
    ))
    file.exists(dest) && file.info(dest)$size > 100L
  }, error = function(e) FALSE, warning = function(w) FALSE)
  isTRUE(ok)
}

olfml2b_select_cran_repo <- function(candidates = olfml2b_repo_candidates(), log_file = NULL) {
  for (repo in unique(candidates)) {
    olfml2b_log("INFO", "Testing CRAN repository: ", repo, log_file = log_file)
    if (olfml2b_probe_repository(repo)) {
      olfml2b_log("INFO", "Selected CRAN repository: ", repo, log_file = log_file)
      return(repo)
    }
    olfml2b_log("WARN", "Repository probe failed: ", repo, log_file = log_file)
  }
  # Return official CRAN so the subsequent installation error remains explicit.
  olfml2b_warn("No CRAN mirror probe succeeded; installation will still try cloud.r-project.org")
  "https://cloud.r-project.org"
}

olfml2b_remove_stale_locks <- function(lib = .libPaths()[1L], log_file = NULL) {
  locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
  if (!length(locks)) return(invisible(character()))
  for (x in locks) {
    olfml2b_log("WARN", "Removing stale package-install lock: ", x, log_file = log_file)
    unlink(x, recursive = TRUE, force = TRUE)
  }
  invisible(locks)
}

olfml2b_install_one_cran <- function(pkg, repos, lib, retries = 3L,
                                   dependencies = c("Depends", "Imports", "LinkingTo"),
                                   log_file = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  repos <- unique(repos)
  for (attempt in seq_len(max(1L, retries))) {
    repo <- repos[((attempt - 1L) %% length(repos)) + 1L]
    olfml2b_log("INFO", sprintf("CRAN install %s | attempt %d/%d | %s", pkg, attempt, retries, repo),
              log_file = log_file)
    options(repos = c(CRAN = repo))
    ok <- tryCatch({
      olfml2b_abort("Runtime CRAN installation is disabled; restore the frozen environment before analysis.")
    }, error = function(e) {
      olfml2b_log("ERROR", "CRAN installation error for ", pkg, ": ", conditionMessage(e),
                log_file = log_file)
      FALSE
    })
    if (isTRUE(ok)) return(TRUE)
    Sys.sleep(min(10, 2 * attempt))
  }
  FALSE
}

olfml2b_install_one_bioc <- function(pkg, lib, retries = 3L,
                                   dependencies = c("Depends", "Imports", "LinkingTo"),
                                   log_file = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  for (attempt in seq_len(max(1L, retries))) {
    olfml2b_log("INFO", sprintf("Bioconductor install %s | attempt %d/%d", pkg, attempt, retries),
              log_file = log_file)
    ok <- tryCatch({
      olfml2b_abort("Runtime Bioconductor installation is disabled; restore the frozen environment before analysis.")
    }, error = function(e) {
      olfml2b_log("ERROR", "Bioconductor installation error for ", pkg, ": ", conditionMessage(e),
                log_file = log_file)
      FALSE
    })
    if (isTRUE(ok)) return(TRUE)
    Sys.sleep(min(10, 2 * attempt))
  }
  FALSE
}

olfml2b_install_packages <- function(cran = character(), bioc = character(), github = character(),
                                   update = FALSE, ncpus = 1L, timeout_sec = 3600L,
                                   retries = 3L, log_file = NULL) {
  # ncpus/update are retained for backward compatibility; sequential installs
  # are intentional on Windows because large parallel batches are more prone to
  # RStudio download timeouts and partial ZIP extraction.
  invisible(update)
  invisible(ncpus)
  lib <- .libPaths()[1L]
  if (!dir.exists(lib)) dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(lib, .libPaths())))

  old <- options(
    timeout = max(as.numeric(getOption("timeout", 60)), as.numeric(timeout_sec)),
    download.file.method = "libcurl",
    pkgType = if (.Platform$OS.type == "windows") "win.binary" else getOption("pkgType")
  )
  on.exit(options(old), add = TRUE)
  olfml2b_remove_stale_locks(lib, log_file = log_file)

  candidates <- olfml2b_repo_candidates()
  selected <- olfml2b_select_cran_repo(candidates, log_file = log_file)
  candidates <- unique(c(selected, candidates))
  options(repos = c(CRAN = selected))

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    ok <- olfml2b_install_one_cran(
      "BiocManager", repos = candidates, lib = lib, retries = retries,
      dependencies = c("Depends", "Imports", "LinkingTo"), log_file = log_file
    )
    olfml2b_assert(ok, "Unable to install BiocManager after retries.")
  }

  cran <- unique(as.character(cran[nzchar(cran)]))
  bioc <- unique(as.character(bioc[nzchar(bioc)]))
  github <- unique(as.character(github[nzchar(github)]))
  status <- list()

  for (pkg in cran) {
    ok <- olfml2b_install_one_cran(
      pkg, repos = candidates, lib = lib, retries = retries,
      dependencies = c("Depends", "Imports", "LinkingTo"), log_file = log_file
    )
    status[[length(status) + 1L]] <- data.frame(
      package = pkg, repository = "CRAN", installed = ok,
      version = olfml2b_package_version(pkg), stringsAsFactors = FALSE
    )
  }

  # Let BiocManager select the release compatible with the installed R version.
  # Do not force the newest Bioconductor release.
  bioc_version <- tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_)
  olfml2b_log("INFO", "Bioconductor version selected for this R: ", bioc_version, log_file = log_file)
  options(repos = BiocManager::repositories())
  for (pkg in bioc) {
    ok <- olfml2b_install_one_bioc(
      pkg, lib = lib, retries = retries,
      dependencies = c("Depends", "Imports", "LinkingTo"), log_file = log_file
    )
    status[[length(status) + 1L]] <- data.frame(
      package = pkg, repository = "Bioconductor", installed = ok,
      version = olfml2b_package_version(pkg), stringsAsFactors = FALSE
    )
  }

  if (length(github)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      olfml2b_install_one_cran("remotes", candidates, lib, retries, log_file = log_file)
    }
    for (repo in github) {
      pkg <- sub(".*/", "", repo)
      ok <- requireNamespace(pkg, quietly = TRUE)
      if (!ok) {
        ok <- FALSE
        olfml2b_log("ERROR", "Runtime GitHub installation is disabled for frozen analysis: ", repo,
                   log_file = log_file)
      }
      status[[length(status) + 1L]] <- data.frame(
        package = pkg, repository = "GitHub", installed = ok,
        version = olfml2b_package_version(pkg), stringsAsFactors = FALSE
      )
    }
  }

  result <- olfml2b_bind_rows_safe(status)
  attr(result, "cran_repository") <- selected
  attr(result, "bioconductor_version") <- bioc_version
  result
}

olfml2b_download <- function(urls, destfile, overwrite = FALSE, min_bytes = 1L,
                           retries = 4L, timeout_sec = 7200L, log_file = NULL) {
  urls <- unique(as.character(urls[nzchar(urls)]))
  olfml2b_assert(length(urls) > 0L, "No download URL supplied for ", destfile)
  if (file.exists(destfile) && !overwrite && file.info(destfile)$size >= min_bytes) {
    olfml2b_log("INFO", "Local file exists; skipping download: ", destfile, log_file = log_file)
    return(normalizePath(destfile, winslash = "/", mustWork = TRUE))
  }
  olfml2b_safe_dir_create(dirname(destfile))
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout %||% 60, timeout_sec))
  methods <- unique(c("libcurl", "auto", if (.Platform$OS.type == "windows") "wininet" else character()))
  errors <- character()
  for (url in urls) {
    for (attempt in seq_len(retries)) {
      for (method in methods) {
        tmp <- paste0(destfile, ".partial")
        if (file.exists(tmp)) unlink(tmp, force = TRUE)
        olfml2b_log("INFO", "Download attempt ", attempt, "/", retries,
                  " | method=", method, " | ", url, log_file = log_file)
        ok <- tryCatch({
          status <- utils::download.file(url, tmp, mode = "wb", quiet = FALSE, method = method)
          identical(status, 0L) && file.exists(tmp) && file.info(tmp)$size >= min_bytes
        }, error = function(e) {
          errors <<- c(errors, paste(url, method, attempt, conditionMessage(e), sep = " | "))
          FALSE
        }, warning = function(w) {
          errors <<- c(errors, paste(url, method, attempt, conditionMessage(w), sep = " | "))
          invokeRestart("muffleWarning")
        })
        if (isTRUE(ok)) {
          if (file.exists(destfile)) unlink(destfile, force = TRUE)
          moved <- file.rename(tmp, destfile)
          if (!moved) {
            moved <- file.copy(tmp, destfile, overwrite = TRUE)
            unlink(tmp, force = TRUE)
          }
          if (moved && file.exists(destfile) && file.info(destfile)$size >= min_bytes) {
            olfml2b_log("INFO", "Downloaded ", basename(destfile), " | ",
                      format(file.info(destfile)$size, big.mark = ","), " bytes", log_file = log_file)
            return(normalizePath(destfile, winslash = "/", mustWork = TRUE))
          }
        }
        if (file.exists(tmp)) unlink(tmp, force = TRUE)
      }
      Sys.sleep(min(2^attempt, 15))
    }
  }
  err_file <- paste0(destfile, ".download_errors.txt")
  olfml2b_atomic_write_lines(errors, err_file)
  olfml2b_abort("Download failed after all official routes: ", paste(urls, collapse = "; "),
              ". See ", err_file)
}

olfml2b_clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

olfml2b_clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

olfml2b_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- olfml2b_clean_text(x)
  y <- gsub(",", "", y, fixed = TRUE)
  # Extract the first valid decimal/scientific numeric token.  The previous
  # implementation placed an escaped hyphen inside a TRE character class;
  # on Windows/R 4.5 this can be interpreted as an invalid range and abort
  # optional subtype harmonisation.  Keeping the sign outside a character
  # class is portable across TRE and PCRE, and perl=TRUE is explicit.
  pattern <- "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
  hit <- regexpr(pattern, y, perl = TRUE)
  start <- as.integer(hit)
  width <- attr(hit, "match.length")
  token <- rep(NA_character_, length(y))
  ok <- is.finite(start) & start > 0L & is.finite(width) & width > 0L
  token[ok] <- substring(y[ok], start[ok], start[ok] + width[ok] - 1L)
  suppressWarnings(as.numeric(token))
}

olfml2b_binary_event <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  raw <- tolower(trimws(olfml2b_clean_text(x)))
  out <- rep(NA_integer_, length(raw))
  num <- suppressWarnings(as.numeric(raw))
  numeric_ok <- is.finite(num) & num %in% c(0, 1)
  out[numeric_ok] <- as.integer(num[numeric_ok])
  # Censoring phrases are evaluated before event phrases because strings such
  # as "no recurrence" and "disease free" contain event-related substrings.
  censor_patterns <- c(
    "^alive$", "^living$", "^censored$", "^no$", "^false$", "^negative$",
    "no[ _-]*recurrence", "non[ _-]*recur", "disease[ _-]*free",
    "without[ _-]*(recurrence|relapse|progression)", "event[ _-]*free",
    "^0$"
  )
  event_patterns <- c(
    "^dead$", "^deceased$", "^died$", "^death$", "^yes$", "^true$", "^positive$",
    "recurred", "recurrence", "relapse", "progressed", "progression",
    "death|deceased|died", "^event$", "^1$"
  )
  for (pat in censor_patterns) out[is.na(out) & grepl(pat, raw, perl = TRUE)] <- 0L
  for (pat in event_patterns) out[is.na(out) & grepl(pat, raw, perl = TRUE)] <- 1L
  out
}

olfml2b_time_to_days <- function(x, source_name = "", unit_override = NULL) {
  value <- olfml2b_numeric(x)
  n <- tolower(source_name)
  unit <- if (!is.null(unit_override) && length(unit_override) &&
              !is.na(unit_override[1]) && nzchar(as.character(unit_override)[1])) {
    tolower(as.character(unit_override)[1])
  } else if (grepl("month|mos|mth", n)) {
    "months"
  } else if (grepl("year|yrs", n)) {
    "years"
  } else if (grepl("day", n)) {
    "days"
  } else {
    NA_character_
  }
  if (is.na(unit)) return(rep(NA_real_, length(value)))
  multiplier <- switch(unit, day = 1, days = 1, month = 30.4375, months = 30.4375,
                       year = 365.25, years = 365.25, NA_real_)
  if (!is.finite(multiplier)) return(rep(NA_real_, length(value)))
  value * multiplier
}

olfml2b_stage_group <- function(x) {
  raw <- toupper(olfml2b_clean_text(x))
  raw <- gsub("PATHOLOGIC|PATHOLOGICAL|CLINICAL|AJCC|STAGE|PT|\"|\\'", "", raw)
  raw <- gsub("[^A-Z0-9]", "", raw)
  out <- rep(NA_character_, length(raw))
  out[grepl("^0", raw)] <- "Stage 0"
  out[grepl("^I($|A|B|C|1)", raw) | grepl("^1", raw)] <- "Stage I"
  out[grepl("^II($|A|B|C|2)", raw) | grepl("^2", raw)] <- "Stage II"
  out[grepl("^III($|A|B|C|3)", raw) | grepl("^3", raw)] <- "Stage III"
  out[grepl("^IV($|A|B|C|4)", raw) | grepl("^4", raw)] <- "Stage IV"
  factor(out, levels = c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"), ordered = TRUE)
}

olfml2b_stage_numeric <- function(x) {
  f <- olfml2b_stage_group(x)
  as.numeric(f) - 1L
}

olfml2b_sex <- function(x) {
  raw <- tolower(olfml2b_clean_text(x))
  out <- rep(NA_character_, length(raw))
  out[grepl("female|woman|^f$", raw)] <- "Female"
  out[grepl("male|man|^m$", raw) & !grepl("female", raw)] <- "Male"
  factor(out, levels = c("Female", "Male"))
}

olfml2b_z <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2L || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}


# Standardized-score helper used throughout the OLFML2B-STAD pipeline.
olfml2b_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

olfml2b_rank_inverse_normal <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  n <- sum(ok)
  if (n > 1L) out[ok] <- stats::qnorm((rank(x[ok], ties.method = "average") - 0.5) / n)
  out
}

olfml2b_tcga_patient <- function(barcode) substr(as.character(barcode), 1L, 12L)

olfml2b_tcga_sample_type_code <- function(barcode) {
  b <- as.character(barcode)
  suppressWarnings(as.integer(substr(b, 14L, 15L)))
}

olfml2b_tcga_is_primary <- function(barcode) olfml2b_tcga_sample_type_code(barcode) == 1L
olfml2b_tcga_is_normal <- function(barcode) olfml2b_tcga_sample_type_code(barcode) == 11L

olfml2b_select_one_tcga_aliquot <- function(barcodes) {
  barcodes <- as.character(barcodes)
  patient <- olfml2b_tcga_patient(barcodes)
  type_code <- olfml2b_tcga_sample_type_code(barcodes)
  priority <- ifelse(type_code == 1L, 1L, ifelse(type_code == 11L, 2L, 9L))
  vial <- substr(barcodes, 16L, 16L)
  order_idx <- order(patient, priority, vial, barcodes)
  keep <- !duplicated(patient[order_idx])
  order_idx[keep]
}

olfml2b_case_when_column <- function(df, candidates, min_nonmissing = 1L) {
  nms <- names(df)
  clean <- olfml2b_clean_names(nms)
  candidates <- olfml2b_clean_names(candidates)
  for (cand in candidates) {
    exact <- which(clean == cand)
    contains <- which(grepl(cand, clean, fixed = TRUE))
    idx <- unique(c(exact, contains))
    for (j in idx) {
      nonmissing <- sum(!is.na(df[[j]]) & nzchar(olfml2b_clean_text(df[[j]])))
      if (nonmissing >= min_nonmissing) return(nms[[j]])
    }
  }
  NA_character_
}

olfml2b_first_nonmissing <- function(...) {
  xs <- list(...)
  if (!length(xs)) return(NULL)
  n <- max(vapply(xs, length, integer(1)))
  out <- rep(NA, n)
  for (x in xs) {
    x <- rep_len(x, n)
    take <- is.na(out) | !nzchar(olfml2b_clean_text(out))
    valid <- !is.na(x) & nzchar(olfml2b_clean_text(x))
    out[take & valid] <- x[take & valid]
  }
  out
}

olfml2b_safe_model_matrix <- function(formula, data) {
  tryCatch(stats::model.matrix(formula, data = data), error = function(e) NULL)
}

olfml2b_count_complete <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (!length(cols)) return(0L)
  sum(stats::complete.cases(df[, cols, drop = FALSE]))
}

olfml2b_save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300, bg = "white") {
  olfml2b_assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required to save plots")
  olfml2b_safe_dir_create(dirname(filename))
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE)
  invisible(filename)
}

olfml2b_base_theme <- function(base_size = 11) {
  olfml2b_assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required")
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.38),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.32),
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, color = "black"),
      plot.subtitle = ggplot2::element_text(color = "#374151"),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", color = "black")
    )
}

olfml2b_pub_palette <- function() {
  c(red = "#C73E3A", blue = "#2F6DB3", grey = "#8C8C8C", dark = "#111827",
    light_red = "#F2B8B5", light_blue = "#BFD2EE", purple = "#6A4C93", orange = "#D97706")
}

olfml2b_fmt_p <- function(p, digits = 3) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(!is.finite(p), "P = NA",
         ifelse(p < 1e-4, "P < 1e-4", paste0("P = ", formatC(p, format = "f", digits = digits))))
}

olfml2b_fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.finite(x), "NA", formatC(x, format = "f", digits = digits))
}

olfml2b_fmt_hr_label <- function(hr, lcl, ucl, p, n = NA, events = NA) {
  paste0(
    "HR = ", olfml2b_fmt_num(hr, 2), " (95% CI ", olfml2b_fmt_num(lcl, 2), "–", olfml2b_fmt_num(ucl, 2), ")\n",
    olfml2b_fmt_p(p, 3),
    ifelse(is.finite(suppressWarnings(as.numeric(n))), paste0(", n = ", n), ""),
    ifelse(is.finite(suppressWarnings(as.numeric(events))), paste0(", events = ", events), "")
  )
}

olfml2b_fmt_rho_label <- function(rho, p, n = NA, fdr = NA) {
  paste0(
    "rho = ", olfml2b_fmt_num(rho, 2), ", ", olfml2b_fmt_p(p, 3),
    ifelse(is.finite(suppressWarnings(as.numeric(fdr))), paste0(", FDR = ", olfml2b_fmt_num(fdr, 3)), ""),
    ifelse(is.finite(suppressWarnings(as.numeric(n))), paste0(", n = ", n), "")
  )
}

olfml2b_save_pub_plot <- function(plot, stem, width = 7, height = 5, dpi = 300, bg = "white") {
  olfml2b_assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required to save plots")
  olfml2b_safe_dir_create(dirname(stem))
  out <- c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), tiff = paste0(stem, ".tiff"))
  ggplot2::ggsave(out[["pdf"]], plot = plot, width = width, height = height, bg = bg, limitsize = FALSE, useDingbats = FALSE)
  ggplot2::ggsave(out[["png"]], plot = plot, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE)
  ggplot2::ggsave(out[["tiff"]], plot = plot, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE, compression = "lzw")
  invisible(out)
}

olfml2b_write_figure_registry <- function(rows, path) {
  if (is.null(rows)) rows <- list()
  if (is.data.frame(rows)) rows <- list(rows)
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  reg <- if (length(rows)) olfml2b_bind_rows_safe(rows) else data.frame()
  if (nrow(reg)) {
    reg$generated_at <- olfml2b_timestamp()
  }
  olfml2b_atomic_write_csv(reg, path)
  invisible(reg)
}

olfml2b_figure_registry_row <- function(figure_id, title, source_table = NA_character_, stat_method = NA_character_, n = NA_character_, caption = NA_character_, stem = NA_character_) {
  data.frame(
    figure_id = figure_id,
    title = title,
    source_table = source_table,
    stat_method = stat_method,
    n = as.character(n),
    output_pdf = ifelse(is.na(stem), NA_character_, paste0(stem, ".pdf")),
    output_png = ifelse(is.na(stem), NA_character_, paste0(stem, ".png")),
    output_tiff = ifelse(is.na(stem), NA_character_, paste0(stem, ".tiff")),
    caption = caption,
    stringsAsFactors = FALSE
  )
}


olfml2b_capture_session <- function(path, extra = list()) {
  info <- utils::capture.output(sessionInfo())
  if (length(extra)) {
    info <- c(info, "", "Additional context:", utils::capture.output(str(extra, max.level = 2)))
  }
  olfml2b_atomic_write_lines(info, path)
}

olfml2b_write_error_record <- function(stage, error, path) {
  record <- data.frame(
    timestamp = olfml2b_timestamp(),
    stage = stage,
    error_class = paste(class(error), collapse = ";"),
    message = conditionMessage(error),
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    old <- tryCatch(olfml2b_read_csv(path), error = function(e) data.frame())
    record <- rbind(old, record)
  }
  olfml2b_atomic_write_csv(record, path)
}

olfml2b_run_stage <- function(stage, expr, log_file, error_file, required = TRUE) {
  olfml2b_log("INFO", "===== ", stage, " =====", log_file = log_file)
  result <- tryCatch(
    force(expr),
    error = function(e) {
      olfml2b_write_error_record(stage, e, error_file)
      olfml2b_log(if (required) "ERROR" else "WARN", stage, ": ", conditionMessage(e), log_file = log_file)
      if (required) stop(e)
      structure(list(stage = stage, error = conditionMessage(e)), class = "olfml2b_optional_failure")
    }
  )
  result
}

olfml2b_is_failure <- function(x) inherits(x, "olfml2b_optional_failure")

olfml2b_object_summary <- function(x, object_name = deparse(substitute(x))) {
  data.frame(
    object = object_name,
    class = paste(class(x), collapse = ";"),
    dimensions = paste(dim(x) %||% length(x), collapse = " x "),
    size_bytes = as.numeric(object.size(x)),
    hash_md5 = olfml2b_hash_object(x),
    stringsAsFactors = FALSE
  )
}

olfml2b_assert_unique <- function(x, label) {
  dup <- unique(x[duplicated(x) & !is.na(x)])
  olfml2b_assert(!length(dup), label, " contains duplicated identifiers: ", paste(utils::head(dup, 20), collapse = ", "))
  invisible(TRUE)
}

olfml2b_make_sample_flow <- function(cohort, steps) {
  olfml2b_assert(is.list(steps) && length(steps), "steps must be a named list")
  data.frame(
    cohort = cohort,
    step_order = seq_along(steps),
    step = names(steps),
    n = as.integer(unlist(steps, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

olfml2b_safe_quantile <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1)) {
  if (!any(is.finite(x))) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

olfml2b_summary_numeric <- function(x, name = "value") {
  q <- olfml2b_safe_quantile(x)
  data.frame(
    variable = name,
    n_total = length(x),
    n_nonmissing = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5],
    stringsAsFactors = FALSE
  )
}

olfml2b_dir_size <- function(path) {
  if (!dir.exists(path)) return(0)
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  sum(file.info(files)$size, na.rm = TRUE)
}

olfml2b_format_bytes <- function(n) {
  units <- c("B", "KB", "MB", "GB", "TB")
  if (!is.finite(n) || n < 0) return(NA_character_)
  idx <- min(floor(log(max(n, 1), 1024)) + 1L, length(units))
  sprintf("%.2f %s", n / 1024^(idx - 1L), units[idx])
}

# ==============================================================================
# Embedded study contract and frozen biological resources.
# Nothing in this file is read from external user-supplied CSV/RDS files.
# Public cohort data are downloaded by the pipeline itself.
# ==============================================================================

olfml2b_version <- function() "v2.3.1_20260710_EMPIRICAL_DIRECTION_AND_ROBUST_META_CONTRACT"

olfml2b_required_cran <- function(mode = "core") {
  mode <- tolower(mode)
  if (!mode %in% c("acquisition", "core", "scrna", "methylation", "full")) {
    stop("Unknown installation mode: ", mode, call. = FALSE)
  }

  # Keep this list intentionally minimal. Only packages directly used by the
  # selected pipeline profile are declared here. Their mandatory Depends,
  # Imports and LinkingTo dependencies are installed automatically. Suggested
  # packages are deliberately excluded because they are not required to run
  # the pipeline and can pull very large Bayesian/modeling stacks.
  acquisition <- c(
    "data.table", "ggplot2", "matrixStats", "readxl", "R.utils",
    "jsonlite", "curl", "httr", "yaml"
  )
  core_extra <- c(
    "survival", "metafor", "sandwich", "lmtest", "mice", "digest"
  )
  # The RC2 production implementation streams processed matrices and uses
  # sparse Matrix operations; it does not require a Seurat object workflow.
  scrna_extra <- c("Matrix")

  x <- acquisition
  if (mode %in% c("core", "methylation", "full")) x <- c(x, core_extra)
  if (mode %in% c("scrna", "full")) x <- c(x, scrna_extra)
  unique(x)
}

olfml2b_required_bioc <- function(mode = "core") {
  mode <- tolower(mode)
  if (!mode %in% c("acquisition", "core", "scrna", "methylation", "full")) {
    stop("Unknown installation mode: ", mode, call. = FALSE)
  }

  acquisition <- c(
    "TCGAbiolinks", "SummarizedExperiment", "S4Vectors", "Biobase",
    "GEOquery", "AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db", "edgeR"
  )
  core_extra <- c("GSVA")
  methylation_extra <- c(
    "IlluminaHumanMethylation450kanno.ilmn12.hg19", "minfi"
  )

  x <- acquisition
  if (mode %in% c("core", "scrna", "methylation", "full")) x <- c(x, core_extra)
  if (mode %in% c("methylation", "full")) x <- c(x, methylation_extra)
  unique(x)
}


olfml2b_known_platform_packages <- function() {
  # Public platform-annotation contract used by original Part2. Keep this
  # deliberately conservative: only map platforms with established
  # Bioconductor annotation packages. Platforms not listed fall back to the
  # GEO GPL/AnnotGPL parser inside Part2 instead of forcing package installs.
  c(
    GPL570 = "hgu133plus2.db",
    GPL96 = "hgu133a.db",
    GPL97 = "hgu133b.db",
    GPL571 = "hgu133a2.db",
    GPL6244 = "hugene10sttranscriptcluster.db",
    GPL6947 = "illuminaHumanv3.db",
    GPL10558 = "illuminaHumanv4.db",
    GPL6884 = "illuminaHumanv4.db"
  )
}

olfml2b_gene_annotation <- function() {
  data.frame(
    official_symbol = "OLFML2B",
    approved_name = "olfactomedin like 2B",
    target_ensembl = "ENSG00000162745",
    target_entrez = "25903",
    target_hgnc = "HGNC:24558",
    target_uniprot = "Q68BL8",
    target_synonyms = "MGC51337;photomedin-2",
    excluded_symbols = "OLFML2A;OLFM2",
    chromosome = "chr1",
    cytoband = "1q23.3",
    genome_build_primary = "GRCh38.p14",
    grch38_start_ensembl = 161983192L,
    grch38_end_ensembl = 162023869L,
    strand = "-",
    coordinate_note = "NCBI Gene 25903 / GRCh38.p14; coordinates are audited metadata and are not used to select results.",
    stringsAsFactors = FALSE
  )
}

olfml2b_default_contract <- function() {
  list(
    project = "OLFML2B_STAD",
    target_gene = "OLFML2B",
    target_aliases = c("OLFML2B", "MGC51337"),
    target_entrez = "25903",
    target_ensembl = "ENSG00000162745",
    target_hgnc = "HGNC:24558",
    target_uniprot = "Q68BL8",
    target_chromosome = "chr1",
    target_cytoband = "1q23.3",
    target_grch38_start = 161983192L,
    target_grch38_end = 162023869L,
    target_strand = "-",
    excluded_symbols = c("OLFML2A", "OLFM2"),
    tcga_project = "TCGA-STAD",
    seed = 20260721L,
    expression_primary = "log2_TPM_plus_1",
    exposure_primary = "within_cohort_z_score",
    primary_survival_endpoint = "RECURRENCE_FAMILY",
    secondary_survival_endpoints = c("OS"),
    recurrence_endpoint_family = "recurrence_survival",
    recurrence_endpoint_policy = "retain_dataset_label_and_pool_only_by_family",
    geo_characteristic_parser_version = "colon_first_v2_4_0",
    target_mapping_audit_policy = "exact_OLFML2B_symbol_or_Entrez_25903; no substring or paralog merging",
    adjustment_primary = c("age", "sex", "stage"),
    no_best_cutpoint_primary = TRUE,
    multiple_testing = "BH",
    alpha = 0.05,
    min_gene_set_overlap = 10L,
    min_survival_events_univariable = 20L,
    min_events_per_parameter = 10L,
    meta_method = "REML_modified_Hartung_Knapp",
    geo_cohorts = c("GSE62254", "GSE15459", "GSE26253", "GSE84437"),
    stage_order = c("Stage I", "Stage II", "Stage III", "Stage IV"),
    expected_direction = "empirical_two_sided",
    resolved_primary_direction = "adverse_high_for_recurrence_family",
    direction_lock = "Direction was frozen after the completed genome-wide discovery/external-validation screen; it is not selected separately within a cohort.",
    allow_negative_results = TRUE,
    main_hypothesis = paste(
      "Higher OLFML2B is tested primarily as an adverse recurrence-family marker and as a marker of CAF/ECM/TGF-beta-rich spatial niches in gastric cancer.",
      "OS is secondary. Bulk and spatial associations do not establish malignant-cell origin, mediation or causality."
    ),
    recurrence = c("recurrence", "recur", "relapse", "disease_recurrence", "recurrence_status"),
    os_time = c("os_time", "overall_survival_time", "overall_survival", "survival_time", "follow_up_time", "followup_time", "days_to_death", "days_to_last_follow_up"),
    os_event = c("os_event", "overall_survival_event", "vital_status", "death_event", "dead", "status"),
    dfs_time = c("dfs_time", "disease_free_survival_time", "disease_free_time", "time_to_recurrence", "recurrence_free_survival_time", "recurrence_free_survival_time_month", "rfs_time", "relapse_free_survival_time"),
    dfs_event = c("dfs_event", "disease_free_survival_event", "recurrence_event", "rfs_event", "relapse_event", "recurrence_status", "status_0_non_recurrence_1_recurrence"),
    treatment = c("treatment", "adjuvant_treatment", "chemotherapy", "chemoradiotherapy", "therapy"),
    tissue = c("tissue", "tissue_type", "sample_type", "source_name_ch1", "tumor_normal", "disease_state")
  )
}


# ==============================================================================
# Embedded public-data and clinical parsing contracts.
# These functions are intentionally defined in Part0 because both Part0 reports
# and Part2 acquisition/harmonisation depend on them. They were present in the
# source UBE2Q2 engine and must not be lost during target-gene adaptation.
# ==============================================================================

olfml2b_embedded_geo_manifest <- function() {
  data.frame(
    accession = c("GSE62254", "GSE15459", "GSE26253", "GSE84437"),
    role = c(
      "external_bulk_primary",
      "external_bulk_replication",
      "external_bulk_recurrence_replication",
      "external_bulk_os_formal_subset"
    ),
    organism = rep("Homo sapiens", 4L),
    platform = c("GPL570", "GPL570", "GPL8432", "GPL6947"),
    expected_samples = c(300L, 192L, 432L, 433L),
    tissue_scope = c(
      "gastric primary tumors",
      "gastric primary tumors after eight official exclusions",
      "FFPE gastric tumors after surgery and adjuvant chemoradiotherapy",
      "GSE84426/GSE84433 formal OS samples; GSE147163 is context-only"
    ),
    endpoint_priority = c("OS;DFS", "OS", "RFS", "OS"),
    matrix_url = c(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/matrix/GSE62254_series_matrix.txt.gz",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15459/matrix/GSE15459_series_matrix.txt.gz",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26253/matrix/GSE26253_series_matrix.txt.gz",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84437/matrix/GSE84437_series_matrix.txt.gz"
    ),
    raw_url = c(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/suppl/GSE62254_RAW.tar",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15459/suppl/GSE15459_RAW.tar",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26253/suppl/GSE26253_RAW.tar",
      NA_character_
    ),
    official_page = paste0(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
      c("GSE62254", "GSE15459", "GSE26253", "GSE84437")
    ),
    stringsAsFactors = FALSE
  )
}

olfml2b_gse15459_exclusions <- function() {
  data.frame(
    sample_title = c(
      "GC-011LGE-T", "GC-021LAH-T", "GC-035PCC-T", "GC-038LYC-T",
      "GC-026-GJK-T", "GC-039-TSC-T", "GC-2000619T", "GC-980327T"
    ),
    reason = c(
      rep("failed_quality_control_in_GEO_record", 4L),
      rep("not_gastric_adenocarcinoma_in_GEO_record", 4L)
    ),
    stringsAsFactors = FALSE
  )
}

olfml2b_geo_prefix <- function(accession) {
  accession <- toupper(trimws(as.character(accession)[1]))
  if (!grepl("^GSE[0-9]+$", accession)) {
    olfml2b_abort("Invalid GEO accession: ", accession)
  }
  num <- suppressWarnings(as.integer(sub("^GSE", "", accession)))
  olfml2b_assert(is.finite(num), "Cannot parse GEO accession: ", accession)
  sprintf("GSE%dnnn", floor(num / 1000L))
}

olfml2b_geo_urls <- function(accession) {
  accession <- toupper(trimws(as.character(accession)[1]))
  prefix <- olfml2b_geo_prefix(accession)
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s", prefix, accession)
  list(
    matrix = sprintf("%s/matrix/%s_series_matrix.txt.gz", base, accession),
    family_soft = sprintf("%s/soft/%s_family.soft.gz", base, accession),
    family_miniml = sprintf("%s/miniml/%s_family.xml.tgz", base, accession),
    raw_tar = sprintf("%s/suppl/%s_RAW.tar", base, accession),
    suppl_base = sprintf("%s/suppl", base),
    page = sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s", accession)
  )
}

olfml2b_public_clinical_routes <- function(accession) {
  accession <- toupper(trimws(as.character(accession)[1]))
  empty <- data.frame(
    source = character(), url = character(), filename = character(),
    stringsAsFactors = FALSE
  )
  routes <- list(
    GSE62254 = data.frame(
      source = "Nature_Medicine_original_article_supplement",
      url = "https://static-content.springer.com/esm/art%3A10.1038%2Fnm.3850/MediaObjects/41591_2015_BFnm3850_MOESM34_ESM.xls",
      filename = "GSE62254_ACRG_Clinical_Information.xls",
      stringsAsFactors = FALSE
    ),
    GSE15459 = data.frame(
      source = "NCBI_GEO_official_supplement",
      url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15459/suppl/GSE15459_outcome.xls",
      filename = "GSE15459_outcome.xls",
      stringsAsFactors = FALSE
    ),
    GSE26253 = empty,
    GSE84437 = empty
  )
  routes[[accession]] %||% empty
}

olfml2b_public_data_source_citations <- function() {
  data.frame(
    resource = c(
      "TCGA-STAD", "GSE62254 expression", "GSE62254 clinical", "GSE15459",
      "GSE26253", "GSE84437 SuperSeries", "GSE84426 formal subset",
      "GSE84433 formal subset", "GSE147163 context-only subset"
    ),
    provider = c(
      "NCI Genomic Data Commons", "NCBI GEO", "Nature Medicine article supplement",
      rep("NCBI GEO", 6L)
    ),
    accession_or_doi = c(
      "TCGA-STAD", "GSE62254", "10.1038/nm.3850 Supplementary Data 1",
      "GSE15459", "GSE26253", "GSE84437", "GSE84426", "GSE84433", "GSE147163"
    ),
    retrieval = c(
      "TCGAbiolinks/GDC API",
      "GEO Series Matrix",
      "direct public supplement download",
      "GEO Series Matrix + outcome supplement",
      "GEO Series Matrix sample characteristics (RFS/status/stage)",
      "GEO SuperSeries matrix partitioned by GSM provenance before transformation",
      "formal 2016 OS subset",
      "formal 2016 OS subset",
      "2020 neqc context-only subset; excluded from formal OS models"
    ),
    stringsAsFactors = FALSE
  )
}

olfml2b_clinical_synonyms <- function() {
  list(
    sample_id = c("geo_accession", "sample_id", "sample", "gsm", "array_id", "accession"),
    patient_id = c("patient_id", "patient", "case_id", "subject_id", "subject", "manuscript_id", "id", "sample_title", "title"),
    title = c("title", "sample_title", "source_name_ch1"),
    age = c("age", "age_at_diagnosis", "age_years", "age_yrs", "patient_age"),
    sex = c("sex", "gender", "patient_sex"),
    stage = c("stage", "pathologic_stage", "pathological_stage", "ajcc_stage", "tnm_stage", "tumor_stage"),
    t_stage = c("t_stage", "pathologic_t", "pt", "tumor_t_stage", "t category"),
    n_stage = c("n_stage", "pathologic_n", "pn", "nodal_stage", "n category"),
    m_stage = c("m_stage", "pathologic_m", "pm", "metastasis_stage", "m category"),
    grade = c("grade", "histologic_grade", "tumor_grade"),
    lauren = c("lauren", "lauren_classification", "histology", "histological_type", "diffuse_intestinal"),
    subtype = c("acr_type", "acrg_subtype", "molecular_subtype", "subtype", "tcga_subtype"),
    recurrence = c("recurrence", "recur", "relapse", "disease_recurrence", "recurrence_status"),
    os_time = c("os_time", "overall_survival_time", "overall_survival", "survival_time", "follow_up_time", "followup_time", "days_to_death", "days_to_last_follow_up"),
    os_event = c("os_event", "overall_survival_event", "vital_status", "death_event", "dead", "status"),
    dfs_time = c("dfs_time", "disease_free_survival_time", "disease_free_time", "time_to_recurrence", "recurrence_free_survival_time", "recurrence_free_survival_time_month", "rfs_time", "relapse_free_survival_time"),
    dfs_event = c("dfs_event", "disease_free_survival_event", "recurrence_event", "rfs_event", "relapse_event", "recurrence_status", "status_0_non_recurrence_1_recurrence"),
    treatment = c("treatment", "adjuvant_treatment", "chemotherapy", "chemoradiotherapy", "therapy"),
    tissue = c("tissue", "tissue_type", "sample_type", "source_name_ch1", "tumor_normal", "disease_state")
  )
}

olfml2b_event_terms <- function() {
  list(
    event = c("dead", "deceased", "died", "death", "recurred", "recurrence", "relapse", "progressed", "progression", "event", "yes", "true", "1"),
    censored = c("alive", "living", "no recurrence", "non-recurred", "disease free", "censored", "no", "false", "0")
  )
}

olfml2b_hallmark_hypoxia <- function() {
  # Frozen human-symbol version used by this project. OLFML2B is removed if present.
  # The list is embedded to prevent runtime drift and dependence on msigdbr/MSigDB login.
  unique(c(
    "ADM","ADORA2B","AK4","AKAP12","ALDOA","ALDOB","ALDOC","AMPD3","ANGPTL4","ANKZF1",
    "ANXA2","ATF3","ATP7A","B3GALT6","BACH1","BCAN","BCL2","BGN","BHLHE40","BNIP3",
    "BTG1","CA9","CASP6","CAV1","CAVIN1","CCNG2","CDKN1A","CITED2","COL5A1","CP",
    "CSRP2","CXCR4","DCN","DDIT3","DDIT4","DDAH1","DPYSL4","DTNA","DUSP1","EDN2",
    "EFNA1","EFNA3","EGFR","EGR1","EMP1","ENO1","ENO2","ENO3","ERO1A","ERRFI1",
    "ETS1","F3","FAM162A","FOS","FOSL2","GAA","GALK1","GAPDH","GBE1","GCK",
    "GPC1","GPI","GRHPR","GYS1","HAS1","HDLBP","HEXA","HK1","HK2","HMOX1",
    "HOXB9","HS3ST1","HSPA5","HSD17B11","IDS","IGFBP1","IGFBP3","IL6","INHA","IRS2",
    "ISG20","JMJD6","JUN","KDM3A","KDM4B","KLF6","KLF7","KLHL24","KRT18","LDHA",
    "LOX","MAP3K1","MIF","MT1E","MT2A","MXI1","MYH9","NAGK","NCAN","NDRG1",
    "NDUFA4L2","NR4A1","P4HA1","P4HA2","PAM","PDK1","PFKFB3","PFKL","PGAM1","PGF",
    "PGK1","PGM1","PIM1","PKLR","PLOD1","PLOD2","PNRC1","PPARGC1A","PPP1R15A","PRDX5",
    "PRKCA","PYGM","RORA","RRAGD","S100A4","SAP30","SCARB1","SDC2","SERPINE1","SIAH2",
    "SLC2A1","SLC2A3","SLC2A14","SLC6A6","SRPX","STC1","STC2","TGFB3","TGFBI","TGM2",
    "TIPARP","TKTL1","TMEM45A","TNFAIP3","TPBG","TPD52","TPI1","VEGFA","VHL","WSB1",
    "XPNPEP1","ZFP36","ZNF292","AHR","ARL4A","BHLHE41","BRS3","C3","CBX3","CDKN1B",
    "CHST2","CITED1","CLEC11A","CX3CR1","ENO4","FAM117B","FOSB","GBE1","HILPDA","HMGCS1",
    "IER3","KLF10","LONP1","MAFF","MALAT1","MMP9","NDRG2","NFIL3","PCK1","PFKP",
    "PLAUR","PLIN2","PPFIA4","RBPJ","RGS2","RNF24","SLC16A3","SLC25A1","SLC25A13","SLC25A30",
    "SLC37A4","SOD2","SPAG4","TES","TFRC","TMEM64","UGP2","VLDLR","WDR45B","ZNF395"
  ))
}

olfml2b_hypoxia_core15 <- function() {
  c("ADM", "ALDOA", "BNIP3", "CA9", "ENO1", "HK2", "LDHA", "MIF", "NDRG1",
    "P4HA1", "PDK1", "PGAM1", "SLC2A1", "TPI1", "VEGFA")
}

olfml2b_hif_targets <- function() {
  unique(c(
    "CA9", "SLC2A1", "SLC2A3", "VEGFA", "LDHA", "BNIP3", "NDRG1", "PDK1",
    "P4HA1", "P4HA2", "PLOD1", "PLOD2", "LOX", "ADM", "ENO1", "HK2", "PGK1",
    "PFKFB3", "MCT4", "SLC16A3", "HILPDA", "ANGPTL4", "ERO1A", "IGFBP3", "STC1", "STC2"
  ))
}

olfml2b_celltype_markers <- function() {
  list(
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1"),
    T_NK = c("CD3D", "CD3E", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1"),
    B = c("CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA"),
    Plasma = c("MZB1", "JCHAIN", "SDC1", "IGHG1", "IGKC"),
    Myeloid = c("LYZ", "LST1", "FCER1G", "TYROBP", "CTSS", "AIF1"),
    Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA"),
    Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2"),
    Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"),
    Smooth_muscle = c("ACTA2", "TAGLN", "MYL9", "RGS5", "MCAM")
  )
}

olfml2b_output_contract <- function() {
  # v2.3.0: Part8 is the frozen PRJEB25780/TIGER anti-PD-1 molecular-context
  # boundary analysis, and Part9 is the integrated evidence registry/runner.
  # The previous integrated-evidence Part8 has been promoted to Part9 to avoid
  # naming conflict with the formal immunotherapy layer.
  data.frame(
    part = c("Part0", "Part1", "Part2", "Part3", "Part4", "Part5", "Part6", "Part7", "Part8", "Part9"),
    required_object = c(
      "Part0_context.rds",
      "Part1_TCGA_STAD.rds",
      "Part2_GEO_bulk_index.rds",
      "Part3_OLFML2B_specialized_bioinformatics_index.rds",
      "Part4_immune_TME_production_index.rds",
      "Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds",
      "Part6_single_cell_evidence_layer_index.rds",
      "Part7_spatial_transcriptomics_index.rds",
      "Part8_Immunotherapy_TIGER_only_molecular_context.rds",
      "Part9_integrated_pipeline_index.rds"
    ),
    evidence_role = c(
      "context_and_directory_contract",
      "TCGA_STAD_expression_clinical_acquisition",
      "GEO_STAD_bulk_cohort_harmonisation",
      "bulk_survival_recurrence_and_module_screen",
      "bulk_TME_association_and_attenuation_analysis",
      "single_PDC000614_within_plex_case_paired_protein_support",
      "single_cell_state_evidence_layer_text_count_plus_atlas",
      "spatial_transcriptomic_niche_evidence_layer",
      "PRJEB25780_TIGER_anti_PD1_molecular_context_boundary",
      "one_click_runner_and_integrated_evidence_registry"
    ),
    evidence_boundary = c(
      "reproducibility infrastructure only",
      "source data harmonisation only",
      "source data harmonisation only",
      "association and prognostic modelling; not causal inference",
      "TME attenuation / ecological state; not formal mediation proof",
      "single-cohort orthogonal support; not multi-proteome validation",
      "cell-state and patient/sample-level ecological support; not direct causal evidence",
      "orthogonal tissue-context validation layer; not population-level prognostic cohort",
      "immunotherapy extrapolation boundary; not standalone anti-PD-1 response biomarker validation",
      "runner, provenance and evidence registry; no primary causal proof"
    ),
    stringsAsFactors = FALSE
  )
}

olfml2b_directory_contract <- function(ctx) {
  values <- unlist(ctx$dirs, use.names = TRUE)
  data.frame(
    key = names(values),
    path = unname(values),
    exists = vapply(values, dir.exists, logical(1)),
    role = ifelse(
      grepl("^raw", names(values)), "downloaded_public_raw_data",
      ifelse(grepl("^cache", names(values)), "reusable_processed_cache",
      ifelse(grepl("^derived|processed", names(values)), "derived_analysis_input",
      ifelse(grepl("^logs", names(values)), "runtime_and_audit_logs",
      ifelse(grepl("^temp", names(values)), "temporary_work_area",
      ifelse(grepl("^release", names(values)), "frozen_release_material",
      ifelse(grepl("tables|figures|reports|objects|final|output|qc", names(values)), "analysis_output", "project_root"))))))),
    stringsAsFactors = FALSE
  )
}

olfml2b_write_embedded_contract <- function(ctx, part_paths) {
  contract_flat <- data.frame(
    parameter = names(ctx$contract),
    value = vapply(ctx$contract, function(x) paste(as.character(x), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  olfml2b_atomic_write_csv(contract_flat, file.path(part_paths$reports, "Part0_embedded_analysis_contract.csv"))
  olfml2b_atomic_write_csv(olfml2b_embedded_geo_manifest(), file.path(part_paths$reports, "Part0_embedded_geo_manifest.csv"))
  olfml2b_atomic_write_csv(olfml2b_gse15459_exclusions(), file.path(part_paths$reports, "Part0_embedded_GSE15459_exclusions.csv"))
  gene_sets <- rbind(
    data.frame(gene_set = "HALLMARK_HYPOXIA_FROZEN", gene = olfml2b_hallmark_hypoxia(), stringsAsFactors = FALSE),
    data.frame(gene_set = "HYPOXIA_CORE_15", gene = olfml2b_hypoxia_core15(), stringsAsFactors = FALSE),
    data.frame(gene_set = "HIF_TARGETS", gene = olfml2b_hif_targets(), stringsAsFactors = FALSE)
  )
  olfml2b_atomic_write_csv(gene_sets, file.path(part_paths$reports, "Part0_embedded_gene_sets.csv"))
  olfml2b_atomic_write_csv(olfml2b_output_contract(), file.path(part_paths$reports, "Part0_output_contract.csv"))
  olfml2b_atomic_write_csv(olfml2b_directory_contract(ctx), file.path(part_paths$reports, "Part0_directory_contract.csv"))
}

run_part0 <- function(
  root = olfml2b_resolve_project_root(),
  allow_internet = TRUE,
  overwrite_downloads = FALSE,
  overwrite_results = TRUE,
  run_scrna = FALSE,
  run_cnv = TRUE,
  run_methylation = FALSE,
  auto_install = FALSE,
  save_large_objects = FALSE,
  cores = max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L)),
  seed = 20260702L
) {
  code_root <- olfml2b_find_code_root()
  dirs <- olfml2b_build_dirs(root)
  provisional_ctx <- list(dirs = dirs)
  part_paths <- olfml2b_part_paths(provisional_ctx, "Part0")
  log_file <- file.path(part_paths$logs, "Part0_Project_Init.log")

  olfml2b_log("INFO", "Initializing ", olfml2b_version(), log_file = log_file)
  olfml2b_log("INFO", "Code root: ", code_root, log_file = log_file)
  olfml2b_log("INFO", "Project root: ", dirs$root, log_file = log_file)

  set.seed(seed)
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
  options(timeout = max(getOption("timeout", 60), 7200))
  options(download.file.method.GEOquery = "auto")

  package_profile <- if (isTRUE(run_scrna) && isTRUE(run_methylation)) {
    "full"
  } else if (isTRUE(run_methylation)) {
    "methylation"
  } else if (isTRUE(run_scrna)) {
    "scrna"
  } else {
    "core"
  }
  requested_packages <- unique(c(
    olfml2b_required_cran(package_profile),
    olfml2b_required_bioc(package_profile)
  ))

  if (isTRUE(auto_install)) {
    olfml2b_install_packages(
      cran = olfml2b_required_cran(package_profile),
      bioc = olfml2b_required_bioc(package_profile),
      update = FALSE,
      ncpus = 1L,
      timeout_sec = 3600L,
      retries = 3L,
      log_file = file.path(dirs$logs_install, "Part0_package_install.log")
    )
  }

  ctx <- list(
    version = olfml2b_version(),
    io_contract_version = "2.5.0",
    created_at = olfml2b_timestamp(),
    code_root = code_root,
    project_root = dirs$root,
    dirs = dirs,
    contract = olfml2b_default_contract(),
    allow_internet = isTRUE(allow_internet),
    overwrite_downloads = isTRUE(overwrite_downloads),
    overwrite_results = isTRUE(overwrite_results),
    run_scrna = isTRUE(run_scrna),
    run_cnv = isTRUE(run_cnv),
    run_methylation = isTRUE(run_methylation),
    save_large_objects = isTRUE(save_large_objects),
    cores = as.integer(cores),
    seed = as.integer(seed),
    package_install_profile = package_profile,
    package_versions_at_init = setNames(
      vapply(requested_packages, olfml2b_package_version, character(1)),
      requested_packages
    )
  )

  olfml2b_write_embedded_contract(ctx, part_paths)

  contamination <- olfml2b_audit_code_root_cleanliness(code_root)
  olfml2b_atomic_write_csv(contamination, file.path(part_paths$reports, "Part0_code_root_cleanliness_audit.csv"))
  if (nrow(contamination)) {
    olfml2b_log(
      "WARN",
      "Working root contains loose runtime files: ", paste(contamination$file, collapse = ", "),
      ". Remove these loose files before continuing; generated outputs should live in subdirectories.",
      log_file = log_file
    )
  }

  context_path <- file.path(dirs$objects, "Part0_context.rds")
  olfml2b_atomic_save_rds(ctx, context_path)
  olfml2b_capture_session(
    file.path(part_paths$reports, "Part0_sessionInfo.txt"),
    extra = list(
      version = ctx$version,
      io_contract_version = ctx$io_contract_version,
      root = ctx$project_root,
      flags = ctx[c("run_scrna", "run_cnv", "run_methylation", "save_large_objects")]
    )
  )
  olfml2b_atomic_write_csv(
    data.frame(
      item = c(
        "version", "io_contract_version", "code_root", "project_root", "allow_internet",
        "overwrite_downloads", "run_scrna", "run_cnv", "run_methylation",
        "save_large_objects", "package_install_profile", "cores", "seed"
      ),
      value = c(
        ctx$version, ctx$io_contract_version, ctx$code_root, ctx$project_root,
        ctx$allow_internet, ctx$overwrite_downloads, ctx$run_scrna, ctx$run_cnv,
        ctx$run_methylation, ctx$save_large_objects, ctx$package_install_profile,
        ctx$cores, ctx$seed
      ),
      stringsAsFactors = FALSE
    ),
    file.path(part_paths$reports, "Part0_run_context.csv")
  )
  olfml2b_log("INFO", "Part0 complete: ", context_path, log_file = log_file)
  invisible(ctx)
}

olfml2b_load_context <- function(root = olfml2b_resolve_project_root()) {
  path <- file.path(normalizePath(root, winslash = "/", mustWork = FALSE), "output", "objects", "Part0_context.rds")
  olfml2b_assert(file.exists(path), "Part0 context is missing. Run run_part0() first: ", path)
  ctx <- readRDS(path)
  olfml2b_assert(!is.null(ctx$io_contract_version), "Legacy Part0 context detected. Rerun Part0 with v2.2.6.")
  olfml2b_assert(identical(as.character(ctx$version), olfml2b_version()),
               "Part0 context version mismatch: cached=", as.character(ctx$version),
               " current=", olfml2b_version(), ". Rerun Part0.")
  ctx
}


# ==============================================================================
# OLFML2B-adapted wrappers and fact-check audit
# ==============================================================================

run_olfml2b_part0 <- function(...) {
  ctx <- run_part0(...)
  ctx$contract <- olfml2b_default_contract()
  olfml2b_atomic_save_rds(ctx, file.path(ctx$dirs$objects, "Part0_context.rds"))
  ann <- olfml2b_gene_annotation()
  olfml2b_atomic_write_csv(ann, file.path(ctx$dirs$tables_root, "Part0", "Part0_OLFML2B_gene_annotation_factcheck.csv"))
  olfml2b_atomic_write_csv(
    data.frame(
      exact_target_symbol = "OLFML2B",
      allowed_synonyms = "MGC51337;photomedin-2",
      forbidden_exact_symbols = "OLFML2A;OLFM2",
      extraction_rule = "Use exact approved gene symbol == 'OLFML2B' or audited Entrez 25903; never use substring matching or merge OLFML2A/OLFM2.",
      stringsAsFactors = FALSE
    ),
    file.path(ctx$dirs$tables_root, "Part0", "Part0_OLFML2B_symbol_disambiguation_guard.csv")
  )
  invisible(ctx)
}

olfml2b_load_context <- function(root = NULL) {
  if (!is.null(root)) Sys.setenv(OLFML2B_STAD_ROOT = normalizePath(root, winslash = "/", mustWork = FALSE))
  path <- file.path(olfml2b_resolve_project_root(), "output", "objects", "Part0_context.rds")
  if (file.exists(path)) readRDS(path) else run_olfml2b_part0(root = olfml2b_resolve_project_root())
}


# ==============================================================================
# CNS-style visual contract extension | v20260709_CNS_VISUAL_CONTRACT
# Purpose: project-wide theme, palette, safe labels, dynamic export, and statistics
# labels for Part1-Part9 native figures. This replaces Part10-style post hoc
# figure assembly with per-Part publication-quality native visualization.
# ==============================================================================

olfml2b_wrap_label <- function(x, width = 55) {
  x <- as.character(x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
}
olfml2b_wrap_axis <- function(x, width = 18) olfml2b_wrap_label(x, width)
olfml2b_wrap_title <- function(x, width = 62) olfml2b_wrap_label(x, width)
olfml2b_wrap_subtitle <- function(x, width = 88) olfml2b_wrap_label(x, width)

olfml2b_safe_factor <- function(x, levels = NULL, reverse = FALSE) {
  x_chr <- as.character(x)
  if (is.null(levels)) levels <- unique(x_chr)
  levels <- unique(as.character(levels))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  if (isTRUE(reverse)) levels <- rev(levels)
  factor(x_chr, levels = levels)
}

olfml2b_pub_palette <- function() {
  c(
    red = "#B2182B",
    blue = "#2166AC",
    dark = "#1F2937",
    grey = "#9CA3AF",
    light_grey = "#F3F4F6",
    green = "#1B9E77",
    purple = "#7570B3",
    orange = "#D95F02",
    light_red = "#F4A6A6",
    light_blue = "#A6BDD7"
  )
}

olfml2b_base_theme <- function(base_size = 11) {
  olfml2b_assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required")
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.36),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.28),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", color = "black", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, color = "black", size = base_size + 2, lineheight = 0.98),
      plot.subtitle = ggplot2::element_text(color = "#374151", hjust = 0, size = base_size - 1, lineheight = 0.98),
      plot.caption = ggplot2::element_text(color = "#6B7280", hjust = 0, size = base_size - 2, lineheight = 0.95),
      plot.margin = ggplot2::margin(8, 14, 8, 10),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.30),
      strip.text = ggplot2::element_text(face = "bold", color = "black", size = base_size - 1, lineheight = 0.95)
    )
}

olfml2b_plot_height <- function(n = NULL, base = 4.6, per_row = 0.30, min_height = 4.2, max_height = 13) {
  n <- suppressWarnings(as.numeric(n))[1]
  if (!is.finite(n) || n <= 0) return(base)
  max(min_height, min(max_height, base + per_row * max(0, n - 6)))
}
olfml2b_plot_width <- function(n = NULL, base = 6.8, per_col = 0.28, min_width = 6.2, max_width = 14) {
  n <- suppressWarnings(as.numeric(n))[1]
  if (!is.finite(n) || n <= 0) return(base)
  max(min_width, min(max_width, base + per_col * max(0, n - 4)))
}

olfml2b_fmt_p <- function(p, digits = 3) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(!is.finite(p), "P=NA",
         ifelse(p < 1e-4, "P<1e-4", paste0("P=", formatC(p, format = "f", digits = digits))))
}
olfml2b_fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.finite(x), "NA", formatC(x, format = "f", digits = digits))
}
olfml2b_fmt_fdr <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.finite(x), "FDR=NA", ifelse(x < 1e-4, "FDR<1e-4", paste0("FDR=", formatC(x, format = "f", digits = digits))))
}
olfml2b_fmt_rho_label <- function(rho, p, n = NA, fdr = NA) {
  paste0("rho=", olfml2b_fmt_num(rho, 2), "; ", olfml2b_fmt_p(p, 3),
         ifelse(is.finite(suppressWarnings(as.numeric(fdr))), paste0("; ", olfml2b_fmt_fdr(fdr, 3)), ""),
         ifelse(is.finite(suppressWarnings(as.numeric(n))), paste0("; n=", n), ""))
}
olfml2b_fmt_hr_label <- function(hr, lcl, ucl, p, n = NA, events = NA) {
  paste0("HR=", olfml2b_fmt_num(hr, 2), " (95% CI ", olfml2b_fmt_num(lcl, 2), "–", olfml2b_fmt_num(ucl, 2), ")\n",
         olfml2b_fmt_p(p, 3),
         ifelse(is.finite(suppressWarnings(as.numeric(n))), paste0("; n=", n), ""),
         ifelse(is.finite(suppressWarnings(as.numeric(events))), paste0("; events=", events), ""))
}

olfml2b_read_csv_if_exists <- function(path, check.names = FALSE) {
  if (length(path) == 0L || is.na(path[1]) || !file.exists(path[1]) || file.info(path[1])$size <= 0) return(data.frame())
  tryCatch(utils::read.csv(path[1], stringsAsFactors = FALSE, check.names = check.names), error = function(e) data.frame())
}

olfml2b_save_pub_plot <- function(plot, stem, width = 7, height = 5, dpi = 600, bg = "white") {
  olfml2b_assert(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 is required to save plots")
  olfml2b_safe_dir_create(dirname(stem))
  out <- c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), tiff = paste0(stem, ".tiff"))
  ggplot2::ggsave(out[["pdf"]], plot = plot, width = width, height = height, bg = bg, limitsize = FALSE, useDingbats = FALSE)
  ggplot2::ggsave(out[["png"]], plot = plot, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE)
  ggplot2::ggsave(out[["tiff"]], plot = plot, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE, compression = "lzw")
  invisible(out)
}

olfml2b_save_cns_plot <- function(plot, stem, kind = "general", n_rows = NULL, n_cols = NULL, width = NULL, height = NULL, dpi = 600) {
  if (is.null(width)) {
    width <- switch(as.character(kind)[1],
                    forest = 8.2, heatmap = olfml2b_plot_width(n_cols, base = 7.4), dotplot = 7.2,
                    spatial = 6.2, waterfall = 8.0, scatter = 6.6, paired = 6.2, 7.0)
  }
  if (is.null(height)) {
    height <- switch(as.character(kind)[1],
                     forest = olfml2b_plot_height(n_rows, base = 4.8, per_row = 0.36),
                     heatmap = olfml2b_plot_height(n_rows, base = 4.6, per_row = 0.28),
                     dotplot = olfml2b_plot_height(n_rows, base = 4.8, per_row = 0.30),
                     spatial = 5.8, waterfall = 4.8, scatter = 5.0, paired = 4.8, 5.0)
  }
  olfml2b_save_pub_plot(plot, stem, width = width, height = height, dpi = dpi)
}

# ==============================================================================
# IF 7-8 locked methodology contract | v2.0.0 | RC2-primary full rebuild
# This final definition block is intentionally last: R resolves these public
# functions at call time, so the locked contract supersedes legacy defaults
# while preserving the complete acquisition/processing implementation above.
# ==============================================================================

OLFML2B_IF78_RELEASE <- "v1.0.3_20260721_SELF_CONTAINED_INTEGRITY_AND_CONTRACT_FIX"

.olfml2b_if78_legacy_contract <- olfml2b_default_contract
olfml2b_default_contract <- function() {
  x <- .olfml2b_if78_legacy_contract()
  x$release <- OLFML2B_IF78_RELEASE
  x$seed <- 20260721L
  x$expression_primary <- "log2_TPM_plus_1_or_platform_appropriate_log2_expression"
  x$exposure_primary <- "within_cohort_z_score_continuous"
  x$primary_survival_endpoint <- "RECURRENCE_FAMILY"
  x$secondary_survival_endpoints <- c("OS")
  x$adjustment_primary <- c("age10", "sex", "stage")
  x$adjustment_policy <- "M1-common is fixed a priori as OLFML2B_z + age10 + sex + stage; age10 is centered age per 10 years; no outcome- or availability-driven covariate selection"
  x$missing_data_policy <- "complete-case only when core-covariate missingness is <=5%; otherwise cohort-level MICE with outcome, event and Nelson-Aalen auxiliary information, m=20, seed=20260721"
  x$min_primary_n <- 80L
  x$min_primary_events <- 40L
  x$min_events_per_parameter <- 10L
  x$meta_method <- "REML_modified_Hartung_Knapp"
  x$meta_primary_min_k <- 3L
  x$multiple_testing <- "BH_within_prespecified_analysis_family"
  x$scrna_cohort <- "RC2_NO_KANG"
  x$scrna_primary_cohorts <- c("GSE150290", "GSE183904")
  x$scrna_depth_context <- "GSE167297"
  x$scrna_supportive_cohort <- "GSE134520"
  x$scrna_excluded_from_formal_layer <- "GSE206785_Kang"
  x$scrna_inference_unit <- "officially_mapped_patient_or_sample; never cells as independent replicates"
  x$tme_model_policy <- "M0/M1/M2 nested models use the same complete-case patient set within cohort-endpoint"
  x$quality_gate_policy <- "structural and diagnostic criteria only; effect direction, magnitude, confidence interval and P value never determine pipeline success"
  x$claim_ceiling <- "observational_multilayer_association_and_context; no causality, clinical utility, treatment selection or standalone ICI prediction"
  x$main_hypothesis <- "OLFML2B is evaluated as a prespecified continuous exposure for prognosis and CAF/TGFb/ECM-rich ecological context in gastric cancer; null and discordant results remain reportable."
  x$expected_direction <- "two_sided_no_result_gate"
  x$resolved_primary_direction <- NA_character_
  x$direction_lock <- "No direction is used for QC, model eligibility or pipeline continuation."
  x
}

.olfml2b_if78_legacy_required_cran <- olfml2b_required_cran
olfml2b_required_cran <- function(mode = "core") {
  unique(c(.olfml2b_if78_legacy_required_cran(mode), "digest", "mice"))
}

.olfml2b_if78_legacy_output_contract <- olfml2b_output_contract
olfml2b_output_contract <- function() {
  x <- .olfml2b_if78_legacy_output_contract()
  x <- x[x$part %in% paste0("Part", 0:5), , drop = FALSE]
  i <- match("Part5", x$part)
  if (!is.na(i)) {
    x$required_object[i] <- "Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds"
    x$evidence_role[i] <- "PDC000614 exact-target within-plex case-paired direction audit"
    x$evidence_boundary[i] <- "Target non-quantification is a coverage result; no protein validation claim without auditable mapping and sufficient matched units"
  }
  rownames(x) <- NULL
  x
}

# Runtime dependency mutation is prohibited in a reproducible analysis. These
# compatibility entry points now fail with an actionable message instead of
# modifying the user's R library during a scientific run.
olfml2b_install_one_cran <- function(...) {
  olfml2b_abort("Runtime package installation is disabled in ", OLFML2B_IF78_RELEASE,
               ". Restore the frozen environment before starting the pipeline.")
}
olfml2b_install_one_bioc <- olfml2b_install_one_cran
olfml2b_install_packages <- olfml2b_install_one_cran

olfml2b_assert_dependencies <- function(run_scrna = FALSE, run_methylation = FALSE) {
  profiles <- "core"
  if (isTRUE(run_scrna)) profiles <- c(profiles, "scrna")
  if (isTRUE(run_methylation)) profiles <- c(profiles, "methylation")
  pkgs <- unique(unlist(lapply(profiles, function(z) c(olfml2b_required_cran(z), olfml2b_required_bioc(z))), use.names = FALSE))
  present <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  audit <- data.frame(
    package = pkgs,
    available = present,
    version = vapply(pkgs, olfml2b_package_version, character(1)),
    stringsAsFactors = FALSE
  )
  if (any(!present)) {
    olfml2b_abort(
      "Frozen-environment preflight failed. Missing packages: ",
      paste(pkgs[!present], collapse = ", "),
      ". Install/restore them before analysis; the pipeline will not install packages at runtime."
    )
  }
  audit
}

.olfml2b_if78_run_part0_core <- run_part0
run_part0 <- function(..., auto_install = FALSE) {
  if (isTRUE(auto_install)) {
    olfml2b_abort("auto_install=TRUE is prohibited by the frozen IF7-8 reproducibility contract.")
  }
  args <- list(...)
  dep_audit <- olfml2b_assert_dependencies(
    run_scrna = isTRUE(args$run_scrna %||% FALSE),
    run_methylation = isTRUE(args$run_methylation %||% FALSE)
  )
  args$auto_install <- FALSE
  ctx <- do.call(.olfml2b_if78_run_part0_core, args)
  ctx$version <- OLFML2B_IF78_RELEASE
  ctx$io_contract_version <- "IF78-2.0"
  ctx$contract <- olfml2b_default_contract()
  ctx$dependency_audit <- dep_audit
  olfml2b_atomic_write_csv(dep_audit, file.path(ctx$dirs$tables_root, "Part0", "Part0_frozen_dependency_preflight.csv"))
  olfml2b_atomic_save_rds(ctx, file.path(ctx$dirs$objects, "Part0_context.rds"))
  invisible(ctx)
}

olfml2b_version <- function() OLFML2B_IF78_RELEASE

# ==============================================================================
# OLFML2B-STAD Part0-5 methodology repair release | v1.2.1
# Synchronizes the active PDC000614 protein scope, same-patient Part4 nested
# models, and output contracts. This final block intentionally supersedes older
# release constants without changing upstream acquisition implementations.
# ==============================================================================
OLFML2B_PART0_5_RELEASE <- "v1.2.1_20260722_PART4_COX_BOOTSTRAP_STABILITY_REPAIR"
OLFML2B_IF78_RELEASE <- OLFML2B_PART0_5_RELEASE

.olfml2b_v120_contract_core <- olfml2b_default_contract
olfml2b_default_contract <- function() {
  x <- .olfml2b_v120_contract_core()
  x$release <- OLFML2B_PART0_5_RELEASE
  x$primary_survival_endpoint <- "RECURRENCE_FAMILY"
  x$secondary_survival_endpoints <- c("OS")
  x$tme_model_policy <- paste(
    "Within each cohort-endpoint-axis family, M0 clinical only, M1 clinical plus OLFML2B,",
    "and M2 clinical plus OLFML2B plus the frozen TME axis use the identical patient set;",
    "OLFML2B attenuation is M1 versus M2 and is bootstrapped at patient level."
  )
  x$tme_clinical_policy <- paste(
    "STRICT_COMMON uses age10, sex and stage when all are structurally available;",
    "otherwise a prespecified AVAILABLE_CLINICAL tier uses the available subset without",
    "outcome- or significance-driven covariate selection."
  )
  x$protein_primary_cohort <- "PDC000614_standardized"
  x$protein_primary_matrix <- "PDC000614 TMT18 sample/reference log-ratio matrix"
  x$protein_primary_estimand <- "within-case, within-plex Tumor-minus-Normal OLFML2B protein log-ratio"
  x$protein_primary_value_type <- "Log_Ratio"
  x$protein_sensitivity_value_type <- "Unshared_Log_Ratio"
  x$protein_mapping_threshold <- 0.90
  x$protein_min_pairs_supportive <- 10L
  x$protein_min_pairs_formal <- 20L
  x$protein_pairing_policy <- paste(
    "Tumor and Normal values must be paired within the same analytical TMT plex before",
    "case-level collapse; cross-plex subtraction is prohibited."
  )
  x$protein_claim_ceiling <- "single_cohort_case_paired_orthogonal_protein_support_only"
  x$quality_gate_policy <- paste(
    "Structural availability, mapping, inferential-unit and diagnostic criteria determine",
    "pipeline evaluability; effect direction, magnitude and P value never determine continuation."
  )
  x
}

.olfml2b_v120_output_contract_core <- olfml2b_output_contract
olfml2b_output_contract <- function() {
  x <- .olfml2b_v120_output_contract_core()
  x <- x[x$part %in% paste0("Part", 0:5), , drop = FALSE]
  i4 <- match("Part4", x$part)
  if (!is.na(i4)) {
    x$required_object[i4] <- "Part4_immune_TME_production_index.rds"
    x$evidence_role[i4] <- "bulk CAF/ECM/TGFb association plus same-patient nested attenuation"
    x$evidence_boundary[i4] <- "Ecological association and coefficient attenuation only; not mediation or causality"
  }
  i5 <- match("Part5", x$part)
  if (!is.na(i5)) {
    x$required_object[i5] <- "Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds"
    x$evidence_role[i5] <- "PDC000614 exact-target within-plex case-paired protein direction audit"
    x$evidence_boundary[i5] <- "Single-cohort orthogonal protein support only; biological evidence status is separated from structural evaluability"
  }
  rownames(x) <- NULL
  x
}

olfml2b_version <- function() OLFML2B_PART0_5_RELEASE

# ==============================================================================
# OLFML2B-STAD Part0-5 complete contract repair | v1.3.0
# Freezes GEO tumor-only roles, source endpoint labels, formal tumor-only Part4
# ecology, cohort-specific clinical adjustment, multiplicity/collinearity/PH
# audits, and PDC000614 analytical-plex robustness diagnostics.
# ==============================================================================
OLFML2B_PART0_5_RELEASE <- "v1.3.0_20260722_GEO_TUMOR_ONLY_TME_ENDPOINT_PLEX_COMPLETE_REPAIR"
OLFML2B_IF78_RELEASE <- OLFML2B_PART0_5_RELEASE

.olfml2b_v130_contract_core <- olfml2b_default_contract
olfml2b_default_contract <- function() {
  x <- .olfml2b_v130_contract_core()
  x$release <- OLFML2B_PART0_5_RELEASE
  x$geo_sample_context_policy <- paste(
    "GSE62254, GSE15459, GSE26253 and the 433-sample formal GSE84437 partition are",
    "tumor-only prognostic cohorts; GSE147163 is tumor-only molecular context.",
    "No GEO cohort is used for tumor-normal expression or diagnostic ROC analysis."
  )
  x$geo_sample_context_levels <- c("Tumor", "Normal", "Unknown")
  x$geo_normal_comparison_allowed <- FALSE
  x$geo_source_endpoint_labels <- c(
    GSE62254 = "DFS", GSE15459 = "OS", GSE26253 = "RFS",
    GSE84437 = "OS", GSE147163 = "NONE_CONTEXT_ONLY"
  )
  x$tme_formal_sample_policy <- paste(
    "All formal Part4 ecological correlations, high-low contrasts, checkpoint analyses,",
    "subtype analyses, quadrant summaries and survival models are tumor-only.",
    "TCGA normal tissues are retained only in a separate descriptive context audit."
  )
  x$tme_core_axes <- c("CAF_Core", "ECM_Remodeling", "TGFb_Response")
  x$tme_exploratory_axes <- c("CAF_TGFb_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
  x$tme_duplicate_axis_policy <- "Cytotoxic_axis is an internal alias of CD8_Cytotoxic and is not tested as an additional hypothesis."
  x$tme_clinical_policy <- paste(
    "TCGA-STAD, GSE62254 and GSE15459 use age10, sex and source-reported overall stage when structurally usable;",
    "GSE26253 uses source-reported overall stage because age/sex are structurally absent;",
    "GSE84437 uses age10, sex, source pT and source pN with source-subseries-stratified baseline hazards."
  )
  x$tme_lrt_policy <- "M0-to-M1 and M1-to-M2 likelihood-ratio P values must be finite for stable nested model families."
  x$tme_multiplicity_policy <- paste(
    "M1 reference models are not repeatedly counted in BH adjustment; M2 core axes and",
    "M2 exploratory axes form separate within-cohort-endpoint multiplicity families."
  )
  x$tme_collinearity_policy <- "Spearman correlation, VIF, condition index and sign reversal are audited for every M2 family."
  x$tme_ph_policy <- "OLFML2B PH review is supplemented by a prespecified OLFML2B-by-log-time interaction sensitivity model."
  x$protein_plex_completeness_policy <- paste(
    "OLFML2B measurement completeness is audited by analytical TMT plex; an entire-plex",
    "non-quantification event is not interpreted as tissue-specific biological absence."
  )
  x$protein_plex_robustness_policy <- paste(
    "Primary paired direction is supplemented by analytical-plex-level direction and",
    "leave-one-plex-out sensitivity; plex sensitivity cannot upgrade the evidence claim."
  )
  x
}

.olfml2b_v130_output_contract_core <- olfml2b_output_contract
olfml2b_output_contract <- function() {
  x <- .olfml2b_v130_output_contract_core()
  i2 <- match("Part2", x$part)
  if (!is.na(i2)) {
    x$evidence_role[i2] <- "four tumor-only GEO prognostic cohorts with source endpoint labels"
    x$evidence_boundary[i2] <- "No GEO normal comparator; GSE26253 is displayed as RFS and GSE62254 as DFS"
  }
  i4 <- match("Part4", x$part)
  if (!is.na(i4)) {
    x$evidence_role[i4] <- "tumor-only CAF/ECM/TGFb ecology plus cohort-specific same-patient nested attenuation"
    x$evidence_boundary[i4] <- "Core and exploratory multiplicity families are separate; attenuation is not mediation"
  }
  i5 <- match("Part5", x$part)
  if (!is.na(i5)) {
    x$evidence_role[i5] <- "PDC000614 within-plex case-paired direction plus plex-completeness and leave-one-plex-out audits"
    x$evidence_boundary[i5] <- "Single-cohort direction remains uncertain when analytical-plex sensitivity is present"
  }
  rownames(x) <- NULL
  x
}

olfml2b_version <- function() OLFML2B_PART0_5_RELEASE
