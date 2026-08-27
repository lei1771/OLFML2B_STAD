# ==============================================================================
# OLFML2B-STAD Part5: PDC000614 case-paired TMT18 protein validation
# Version: v1.1.0_20260721_PDC000614_CASE_PAIRED_TMT18_PRODUCTION
# ==============================================================================
# Active protein cohort: PDC000614_standardized only.
# Primary estimand: within-case Tumor minus Normal OLFML2B protein log-ratio
# reconstructed from sample/reference TMT18 values after audited channel ->
# sample -> biospecimen -> case/tissue mapping.
#
# Frozen analysis hierarchy:
#   * Primary value type: Log_Ratio (all quantified peptides).
#   * Sensitivity value type: Unshared_Log_Ratio.
#   * Primary inferential unit: case-level paired Tumor-Normal delta.
#   * Technical/sample repeats are collapsed by median within case/tissue.
#   * Mapping coverage >=90% is required for formal inference.
#   * >=10 paired cases permits supportive direction analysis.
#   * >=20 paired cases permits formal single-cohort orthogonal support.
#   * Direction/significance never determines pipeline completion.
#   * No survival or clinical-prediction claim is made from this protein cohort.
# ==============================================================================

options(stringsAsFactors = FALSE)
OLFML2B_PART5_VERSION <- "v1.1.0_20260721_PDC000614_CASE_PAIRED_TMT18_PRODUCTION"
OLFML2B_PART5_COHORT <- "PDC000614_standardized"
OLFML2B_PART5_TARGET <- "OLFML2B"
OLFML2B_PART5_PRIMARY_VALUE_TYPE <- "Log_Ratio"
OLFML2B_PART5_SENSITIVITY_VALUE_TYPE <- "Unshared_Log_Ratio"
OLFML2B_PART5_MAPPING_COVERAGE_MIN <- 0.90
OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE <- 10L
OLFML2B_PART5_MIN_PAIRS_FORMAL <- 20L
OLFML2B_PART5_BOOTSTRAP_B <- 2000L
OLFML2B_PART5_SEED <- 20260721L

.ol_p5_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.ol_p5_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.ol_p5_code_root <- if (nzchar(.ol_p5_env_root) && dir.exists(.ol_p5_env_root)) {
  normalizePath(.ol_p5_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.ol_p5_entry) && file.exists(.ol_p5_entry)) {
  dirname(normalizePath(.ol_p5_entry, winslash = "/", mustWork = TRUE))
} else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .ol_p5_code_root)

if (!exists("olfml2b_build_dirs", mode = "function")) {
  p0 <- file.path(.ol_p5_code_root, "00_OLFML2B_PART0_CONFIG_CORE.R")
  if (!file.exists(p0)) stop("Missing Part0 core script.", call. = FALSE)
  sys.source(p0, envir = parent.frame(), chdir = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

ol_p5_abort <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
ol_p5_log <- function(level = "INFO", ..., log_file = NULL) {
  line <- sprintf("[%s] [%s] [OLFML2B-P5] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  message(line)
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
  invisible(line)
}

ol_p5_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  invisible(path)
}

ol_p5_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(x, tmp)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }
  invisible(path)
}

ol_p5_bind_rows <- function(xs, schema = NULL) {
  xs <- Filter(is.data.frame, xs)
  if (!length(xs)) {
    if (is.null(schema)) return(data.frame())
    return(as.data.frame(setNames(lapply(schema, function(z) z[0]), names(schema)), stringsAsFactors = FALSE))
  }
  cols <- unique(c(names(schema %||% list()), unlist(lapply(xs, names), use.names = FALSE)))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}

ol_p5_clean_text <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^[\\\"']+|[\\\"']+$", "", x, perl = TRUE)
  x
}

ol_p5_norm_id <- function(x) {
  x <- toupper(ol_p5_clean_text(x))
  gsub("[^A-Z0-9]+", "", x)
}

ol_p5_gene_tokens <- function(x) {
  x <- toupper(ol_p5_clean_text(x))
  unique(Filter(nzchar, unlist(strsplit(x, "[;,|/[:space:]]+", perl = TRUE))))
}

ol_p5_has_exact_gene <- function(x, gene) {
  gene <- toupper(gene)
  vapply(x, function(z) gene %in% ol_p5_gene_tokens(z), logical(1))
}

ol_p5_safe_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", ol_p5_clean_text(x), fixed = TRUE)))

ol_p5_make_dirs <- function(root, output_subdir = "Part5") {
  dirs <- list(
    tables = file.path(root, "output", "tables", output_subdir),
    figures = file.path(root, "output", "figures", output_subdir),
    reports = file.path(root, "output", "reports", output_subdir),
    qc = file.path(root, "output", "qc", output_subdir),
    logs = file.path(root, "logs", "runtime", output_subdir),
    objects = file.path(root, "output", "objects")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

ol_p5_reset_outputs <- function(dirs) {
  for (nm in c("tables", "figures", "reports", "qc", "logs")) {
    d <- dirs[[nm]]
    if (dir.exists(d)) unlink(list.files(d, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE, force = TRUE)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  old_objects <- list.files(dirs$objects, pattern = "^Part5_PDC.*\\.rds$|^Part5_PDC_production_index\\.rds$", full.names = TRUE)
  if (length(old_objects)) unlink(old_objects, force = TRUE)
  invisible(TRUE)
}

ol_p5_find_one <- function(directory, pattern, label, required = TRUE) {
  hit <- list.files(directory, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  hit <- hit[file.info(hit)$isdir %in% FALSE]
  if (!length(hit)) {
    if (required) ol_p5_abort("Missing required ", label, " under: ", directory)
    return(NA_character_)
  }
  if (length(hit) > 1L) {
    hit <- hit[order(file.info(hit)$size, decreasing = TRUE)]
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

ol_p5_target_file_specs <- function(pdc614_dir) {
  protein_dir <- file.path(pdc614_dir, "protein_assembly")
  metadata_dir <- file.path(pdc614_dir, "metadata")
  if (!dir.exists(protein_dir)) ol_p5_abort("PDC000614 protein_assembly directory is missing: ", protein_dir)
  if (!dir.exists(metadata_dir)) ol_p5_abort("PDC000614 metadata directory is missing: ", metadata_dir)
  files <- c(
    tmt18 = ol_p5_find_one(protein_dir, "tmt18\\.tsv$", "TMT18 matrix"),
    summary = ol_p5_find_one(protein_dir, "summary\\.tsv$", "protein summary", required = FALSE),
    primary_candidate = file.path(protein_dir, "PDC000614_proteome_primary_candidate.tsv"),
    sample_map = ol_p5_find_one(protein_dir, "sample\\.txt$", "sample map"),
    label_map = ol_p5_find_one(protein_dir, "label\\.txt$", "label map", required = FALSE),
    biospecimen = file.path(metadata_dir, "PDC000614_biospecimen_latest.tsv")
  )
  present <- !is.na(files) & nzchar(files) & file.exists(files)
  size_mb <- rep(NA_real_, length(files))
  md5 <- rep(NA_character_, length(files))
  size_mb[present] <- round(file.info(files[present])$size / 1024^2, 3)
  md5[present] <- unname(tools::md5sum(files[present]))
  data.frame(
    component = names(files),
    file = unname(files),
    present = present,
    required = names(files) %in% c("tmt18", "sample_map", "biospecimen"),
    size_mb = size_mb,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

ol_p5_read_tabular <- function(path) {
  if (!file.exists(path)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- tryCatch(data.table::fread(path, data.table = FALSE, showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.null(x)) return(x)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

ol_p5_header_tokens <- function(path) {
  line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (!length(line)) ol_p5_abort("Empty TMT18 file: ", path)
  ol_p5_clean_text(strsplit(line, "\t", fixed = TRUE)[[1]])
}


ol_p5_exact_gene_detection_first_field <- function(path, gene = OLFML2B_PART5_TARGET, chunk_size = 5000L) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(data.frame(file = path, gene = gene, exact_rows = 0L, status = "FILE_MISSING", stringsAsFactors = FALSE))
  }
  con <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  readLines(con, n = 1L, warn = FALSE)
  n <- 0L
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE, encoding = "UTF-8")
    if (!length(lines)) break
    first <- sub("\\t.*$", "", lines)
    n <- n + sum(ol_p5_has_exact_gene(first, gene))
  }
  data.frame(file = normalizePath(path, winslash = "/", mustWork = FALSE), gene = gene, exact_rows = n, status = ifelse(n > 0L, "EXACT_TARGET_DETECTED", "TARGET_NOT_DETECTED"), stringsAsFactors = FALSE)
}

ol_p5_quant_map_from_header <- function(header_tokens) {
  h <- ol_p5_clean_text(header_tokens)
  unshared <- grepl("(\\.|[ _])Unshared(\\.|[ _])Log(\\.|[ _])Ratio$", h, ignore.case = TRUE, perl = TRUE)
  all_ratio <- grepl("(\\.|[ _])Log(\\.|[ _])Ratio$", h, ignore.case = TRUE, perl = TRUE)
  pos <- which(unshared | all_ratio)
  if (!length(pos)) ol_p5_abort("No TMT18 Log Ratio columns were detected in the raw header.")
  value_type <- ifelse(unshared[pos], "Unshared_Log_Ratio", "Log_Ratio")
  sample_id <- h[pos]
  sample_id <- sub("(\\.|[ _])Unshared(\\.|[ _])Log(\\.|[ _])Ratio$", "", sample_id, ignore.case = TRUE, perl = TRUE)
  sample_id <- sub("(\\.|[ _])Log(\\.|[ _])Ratio$", "", sample_id, ignore.case = TRUE, perl = TRUE)
  data.frame(
    position = pos,
    quant_column = h[pos],
    sample_id = ol_p5_clean_text(sample_id),
    sample_id_norm = ol_p5_norm_id(sample_id),
    value_type = value_type,
    stringsAsFactors = FALSE
  )
}

ol_p5_scan_selected_rows <- function(path, genes, quant_map, chunk_size = 4000L) {
  genes <- unique(toupper(genes))
  header <- ol_p5_header_tokens(path)
  con <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  readLines(con, n = 1L, warn = FALSE)
  rows <- list(); line_no <- 1L; mismatch <- list()
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE, encoding = "UTF-8")
    if (!length(lines)) break
    for (line in lines) {
      line_no <- line_no + 1L
      first_field <- sub("\t.*$", "", line)
      tokens_gene <- ol_p5_gene_tokens(first_field)
      matched <- intersect(genes, tokens_gene)
      if (!length(matched)) next
      tokens <- ol_p5_clean_text(strsplit(line, "\t", fixed = TRUE)[[1]])
      if (length(tokens) != length(header)) {
        mismatch[[length(mismatch) + 1L]] <- data.frame(
          line_number = line_no, first_field = first_field,
          n_header_tokens = length(header), n_row_tokens = length(tokens),
          stringsAsFactors = FALSE
        )
        next
      }
      for (g in matched) {
        vals <- ol_p5_safe_num(tokens[quant_map$position])
        rows[[length(rows) + 1L]] <- data.frame(
          gene = g,
          source_gene_field = first_field,
          source_line_number = line_no,
          quant_column = quant_map$quant_column,
          sample_id = quant_map$sample_id,
          sample_id_norm = quant_map$sample_id_norm,
          value_type = quant_map$value_type,
          log_ratio = vals,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    values = ol_p5_bind_rows(rows),
    mismatch = ol_p5_bind_rows(mismatch),
    header = header
  )
}

ol_p5_classify_sample_id <- function(x) {
  xx <- toupper(ol_p5_clean_text(x))
  out <- rep("SAMPLE", length(xx))
  out[is.na(xx) | !nzchar(xx)] <- "MISSING"
  out[grepl("POOL|REFERENCE|(^|[^A-Z])REF([^A-Z]|$)|126C", xx, perl = TRUE)] <- "POOL_REFERENCE"
  out[grepl("^NCI|NCI7", xx)] <- "NCI_INTERNAL"
  out[grepl("^T[0-9]+$", xx)] <- "INTERNAL_T_SAMPLE_CODE"
  out
}

ol_p5_flatten_sample_map <- function(sample_map) {
  if (!nrow(sample_map)) return(data.frame())
  channel_cols <- names(sample_map)[grepl("^X?[0-9]{3}[NC]$|^X?13[0-9][NC]$", names(sample_map), ignore.case = TRUE)]
  if (!length(channel_cols)) channel_cols <- names(sample_map)[grepl("^X[0-9]", names(sample_map), ignore.case = TRUE)]
  analytical <- names(sample_map)[grepl("AnalyticalSample|Analytical.Sample|Experiment|Plex|TMT", names(sample_map), ignore.case = TRUE)]
  if (!length(analytical)) analytical <- names(sample_map)[1]
  if (!length(channel_cols)) return(data.frame())
  out <- lapply(seq_len(nrow(sample_map)), function(i) data.frame(
    analytical_sample = as.character(sample_map[[analytical[1]]][i]),
    channel_raw = channel_cols,
    channel = sub("^X", "", channel_cols),
    sample_id = as.character(unlist(sample_map[i, channel_cols, drop = TRUE])),
    stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, out)
  out$sample_id <- ol_p5_clean_text(out$sample_id)
  out$sample_id_norm <- ol_p5_norm_id(out$sample_id)
  out$sample_id_class <- ol_p5_classify_sample_id(out$sample_id)
  rownames(out) <- NULL
  out
}

ol_p5_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) NA_character_ else x[1]
}

ol_p5_classify_tissue <- function(x) {
  xx <- toupper(paste(as.character(x), collapse = " | "))
  if (grepl("NORMAL|ADJACENT|NAT|NON[-_ ]?TUMOU?R|NON[-_ ]?TUMOR|NOT TUMOU?R|NOT TUMOR", xx, perl = TRUE)) return("Normal")
  if (grepl("TUMOU?R|TUMOR|PRIMARY|MALIGNANT|CARCINOMA|CANCER", xx, perl = TRUE)) return("Tumor")
  if (grepl("POOL|REFERENCE|(^|[^A-Z])REF([^A-Z]|$)", xx, perl = TRUE)) return("Reference")
  NA_character_
}

ol_p5_biospecimen_index <- function(biospec) {
  if (!nrow(biospec)) return(list(data = biospec, index = data.frame()))
  b <- as.data.frame(lapply(biospec, as.character), stringsAsFactors = FALSE, check.names = FALSE)
  idx <- list()
  for (r in seq_len(nrow(b))) {
    for (cc in names(b)) {
      v <- ol_p5_clean_text(b[[cc]][r])
      if (is.na(v) || !nzchar(v)) next
      idx[[length(idx) + 1L]] <- data.frame(value_norm = ol_p5_norm_id(v), row = r, column = cc, stringsAsFactors = FALSE)
    }
  }
  list(data = b, index = unique(ol_p5_bind_rows(idx)))
}

ol_p5_case_from_biospecimen <- function(row) {
  nms <- names(row)
  pri <- nms[grepl("case.*submitter|case.*id|participant|patient|subject|donor", nms, ignore.case = TRUE)]
  for (cc in pri) {
    v <- ol_p5_first_nonempty(row[[cc]])
    if (!is.na(v)) return(v)
  }
  NA_character_
}

ol_p5_tissue_from_biospecimen <- function(row) {
  nms <- names(row)
  pri <- nms[grepl("sample.*type|tissue.*type|specimen|diagnosis|descriptor|tumou?r|normal", nms, ignore.case = TRUE)]
  txt <- if (length(pri)) unlist(row[1, pri, drop = TRUE]) else unlist(row[1, , drop = TRUE])
  ol_p5_classify_tissue(txt)
}

ol_p5_map_quant_to_biospecimen <- function(quant_map, sample_long, biospec) {
  qm <- quant_map
  qm$sample_id_class <- ol_p5_classify_sample_id(qm$sample_id)
  qm$analytical_sample <- NA_character_
  qm$channel <- NA_character_
  if (nrow(sample_long)) {
    m <- match(qm$sample_id_norm, sample_long$sample_id_norm)
    qm$analytical_sample <- sample_long$analytical_sample[m]
    qm$channel <- sample_long$channel[m]
  }
  bi <- ol_p5_biospecimen_index(biospec)
  rec <- lapply(seq_len(nrow(qm)), function(i) {
    sid <- qm$sample_id_norm[i]
    hit <- unique(bi$index$row[bi$index$value_norm == sid])
    status <- if (!length(hit)) "NO_BIOSPECIMEN_MATCH" else if (length(hit) == 1L) "UNIQUE_EXACT_NORMALIZED_MATCH" else "AMBIGUOUS_MULTIPLE_MATCHES"
    row <- if (length(hit) == 1L) bi$data[hit, , drop = FALSE] else NULL
    data.frame(
      quant_column = qm$quant_column[i],
      sample_id = qm$sample_id[i],
      biospecimen_match_status = status,
      biospecimen_row = if (length(hit) == 1L) hit else NA_integer_,
      case_id = if (!is.null(row)) ol_p5_case_from_biospecimen(row) else NA_character_,
      tissue_class = if (!is.null(row)) ol_p5_tissue_from_biospecimen(row) else if (qm$sample_id_class[i] == "POOL_REFERENCE") "Reference" else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  rec <- ol_p5_bind_rows(rec)
  out <- merge(qm, rec, by = c("quant_column", "sample_id"), all.x = TRUE, sort = FALSE)
  out$usable_for_tumor_normal <- out$sample_id_class == "SAMPLE" &
    out$biospecimen_match_status == "UNIQUE_EXACT_NORMALIZED_MATCH" &
    out$tissue_class %in% c("Tumor", "Normal") &
    !is.na(out$case_id) & nzchar(out$case_id)
  out
}

ol_p5_mapping_summary <- function(mapped) {
  primary <- mapped[mapped$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  eligible <- primary$sample_id_class == "SAMPLE"
  denom <- sum(eligible)
  usable <- sum(primary$usable_for_tumor_normal & eligible, na.rm = TRUE)
  data.frame(
    metric = c(
      "n_quant_columns_all_value_types", "n_primary_quant_columns", "n_primary_sample_columns",
      "n_primary_usable_tumor_normal_columns", "primary_mapping_coverage",
      "n_primary_tumor_columns", "n_primary_normal_columns", "n_unique_mapped_cases"
    ),
    value = c(
      nrow(mapped), nrow(primary), denom, usable,
      if (denom > 0) usable / denom else NA_real_,
      sum(primary$usable_for_tumor_normal & primary$tissue_class == "Tumor", na.rm = TRUE),
      sum(primary$usable_for_tumor_normal & primary$tissue_class == "Normal", na.rm = TRUE),
      length(unique(primary$case_id[primary$usable_for_tumor_normal]))
    ),
    stringsAsFactors = FALSE
  )
}

ol_p5_case_tissue_aggregate <- function(values) {
  d <- values[values$usable_for_tumor_normal & is.finite(values$log_ratio), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  keys <- interaction(d$gene, d$value_type, d$case_id, d$tissue_class, drop = TRUE, lex.order = TRUE)
  out <- lapply(split(d, keys), function(z) data.frame(
    gene = z$gene[1], value_type = z$value_type[1], case_id = z$case_id[1], tissue_class = z$tissue_class[1],
    n_channels = nrow(z), median_log_ratio = stats::median(z$log_ratio, na.rm = TRUE),
    mean_log_ratio = mean(z$log_ratio, na.rm = TRUE),
    iqr_log_ratio = stats::IQR(z$log_ratio, na.rm = TRUE), stringsAsFactors = FALSE
  ))
  ol_p5_bind_rows(out)
}

ol_p5_pair_case_tissue <- function(case_tissue) {
  if (!nrow(case_tissue)) return(data.frame())
  tum <- case_tissue[case_tissue$tissue_class == "Tumor", , drop = FALSE]
  nor <- case_tissue[case_tissue$tissue_class == "Normal", , drop = FALSE]
  if (!nrow(tum) || !nrow(nor)) return(data.frame())
  names(tum)[names(tum) %in% c("n_channels", "median_log_ratio", "mean_log_ratio", "iqr_log_ratio")] <- paste0(names(tum)[names(tum) %in% c("n_channels", "median_log_ratio", "mean_log_ratio", "iqr_log_ratio")], "_tumor")
  names(nor)[names(nor) %in% c("n_channels", "median_log_ratio", "mean_log_ratio", "iqr_log_ratio")] <- paste0(names(nor)[names(nor) %in% c("n_channels", "median_log_ratio", "mean_log_ratio", "iqr_log_ratio")], "_normal")
  out <- merge(tum, nor, by = c("gene", "value_type", "case_id"), all = FALSE, sort = FALSE)
  out$paired_delta_tumor_minus_normal <- out$median_log_ratio_tumor - out$median_log_ratio_normal
  out
}

ol_p5_bootstrap_median_ci <- function(x, B = OLFML2B_PART5_BOOTSTRAP_B, seed = OLFML2B_PART5_SEED) {
  x <- x[is.finite(x)]
  if (length(x) < 5L || B < 100L) return(c(low = NA_real_, high = NA_real_))
  set.seed(seed)
  z <- replicate(B, stats::median(sample(x, length(x), replace = TRUE)))
  stats::quantile(z, c(0.025, 0.975), na.rm = TRUE, names = FALSE) |>
    setNames(c("low", "high"))
}

ol_p5_direction_test <- function(paired, value_type, mapping_pass, B = OLFML2B_PART5_BOOTSTRAP_B) {
  d <- paired[paired$gene == OLFML2B_PART5_TARGET & paired$value_type == value_type, , drop = FALSE]
  x <- d$paired_delta_tumor_minus_normal
  x <- x[is.finite(x)]
  n <- length(x)
  ci <- if (mapping_pass) ol_p5_bootstrap_median_ci(x, B = B, seed = OLFML2B_PART5_SEED + ifelse(value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, 0L, 1L)) else c(low = NA_real_, high = NA_real_)
  wp <- if (mapping_pass && n >= OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE) tryCatch(stats::wilcox.test(x, mu = 0, alternative = "two.sided", exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
  nz <- x[x != 0]
  sp <- if (mapping_pass && length(nz) >= OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE) stats::binom.test(sum(nz > 0), length(nz), p = 0.5, alternative = "two.sided")$p.value else NA_real_
  role <- if (!mapping_pass) "INFERENCE_BLOCKED_MAPPING" else if (value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE && n >= OLFML2B_PART5_MIN_PAIRS_FORMAL) "FORMAL_PRIMARY_SINGLE_COHORT" else if (n >= OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE) "SUPPORTIVE_SENSITIVITY" else "DESCRIPTIVE_INSUFFICIENT_PAIRS"
  direction <- if (!n) "NOT_EVALUABLE" else if (is.finite(ci[1]) && ci[1] > 0) "POSITIVE_CI_EXCLUDES_ZERO" else if (is.finite(ci[2]) && ci[2] < 0) "NEGATIVE_CI_EXCLUDES_ZERO" else if (stats::median(x) > 0) "POSITIVE_UNCERTAIN" else if (stats::median(x) < 0) "NEGATIVE_UNCERTAIN" else "NULL_MEDIAN"
  data.frame(
    gene = OLFML2B_PART5_TARGET, cohort = OLFML2B_PART5_COHORT, value_type = value_type,
    n_paired_cases = n, median_delta = if (n) stats::median(x) else NA_real_, mean_delta = if (n) mean(x) else NA_real_,
    bootstrap_median_ci_low = unname(ci[1]), bootstrap_median_ci_high = unname(ci[2]),
    positive_fraction = if (n) mean(x > 0) else NA_real_, rank_biserial_sign = if (length(nz)) (sum(nz > 0) - sum(nz < 0)) / length(nz) else NA_real_,
    wilcox_signed_rank_p = wp, exact_sign_p = sp, inference_role = role, evidence_direction = direction,
    claim_boundary = "Case-paired single-cohort protein direction; no causality, prognosis, or clinical utility claim.",
    stringsAsFactors = FALSE
  )
}

ol_p5_catalog <- function() {
  list(
    CAF_Core = c("FAP","ACTA2","PDGFRA","PDGFRB","TAGLN","THY1","DCN","LUM","POSTN"),
    ECM_Remodeling = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","SPARC","MMP2","MMP14","LOX","PLOD2"),
    TGFb_Response = c("TGFBI","SERPINE1","SMAD3","SMAD7","CTGF","THBS1","ITGA5","PMEPA1","INHBA"),
    Epithelial_Control = c("EPCAM","KRT8","KRT18","KRT19"),
    Proliferation_Control = c("MKI67","PCNA","MCM2","MCM5"),
    Cytotoxic_Control = c("CD8A","CD8B","GZMA","GZMB","PRF1","NKG7","GNLY","IFNG")
  )
}

ol_p5_sentinel_catalog <- function() {
  data.frame(
    gene = c("EPCAM","KRT8","KRT18","KRT19","MKI67","FAP","ACTA2","COL1A1","COL1A2","COL3A1","SPARC","TGFBI","SERPINE1","CTGF","INHBA","PGC","PGA3","TFF1","TFF2","CD8A","CD8B","GZMB","NKG7"),
    group = c(rep("Tumor_Epithelial",5), rep("CAF_ECM",6), rep("TGFb",4), rep("Gastric_Normal",4), rep("Immune_Control",4)),
    expected_direction = c(rep("POSITIVE",15), rep("NEGATIVE",4), rep("NO_FIXED_EXPECTATION",4)),
    stringsAsFactors = FALSE
  )
}

ol_p5_sentinel_summary <- function(paired, sentinel_catalog) {
  empty <- data.frame(gene = character(), n_pairs = integer(), median_delta = numeric(), mean_delta = numeric(), positive_fraction = numeric(), wilcox_p = numeric(), stringsAsFactors = FALSE)
  if (!nrow(paired) || !all(c("value_type","gene","paired_delta_tumor_minus_normal") %in% names(paired))) {
    out <- merge(sentinel_catalog, empty, by = "gene", all.x = TRUE, sort = FALSE)
  } else {
    d <- paired[paired$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & paired$gene %in% sentinel_catalog$gene, , drop = FALSE]
    rows <- lapply(split(d, d$gene), function(z) {
      x <- z$paired_delta_tumor_minus_normal[is.finite(z$paired_delta_tumor_minus_normal)]
      wp <- if (length(x) >= 5L) tryCatch(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
      data.frame(gene = z$gene[1], n_pairs = length(x), median_delta = if (length(x)) stats::median(x) else NA_real_, mean_delta = if (length(x)) mean(x) else NA_real_, positive_fraction = if (length(x)) mean(x > 0) else NA_real_, wilcox_p = wp, stringsAsFactors = FALSE)
    })
    out <- merge(sentinel_catalog, if (length(rows)) ol_p5_bind_rows(rows) else empty, by = "gene", all.x = TRUE, sort = FALSE)
  }
  out$direction_concordant <- with(out, ifelse(expected_direction == "POSITIVE" & is.finite(median_delta), median_delta > 0, ifelse(expected_direction == "NEGATIVE" & is.finite(median_delta), median_delta < 0, NA)))
  out$wilcox_fdr <- stats::p.adjust(out$wilcox_p, method = "BH")
  out
}

ol_p5_orientation_qc <- function(sentinel_summary) {
  pos <- sentinel_summary[sentinel_summary$expected_direction == "POSITIVE" & is.finite(sentinel_summary$n_pairs) & sentinel_summary$n_pairs >= 5L & !is.na(sentinel_summary$direction_concordant), , drop = FALSE]
  neg <- sentinel_summary[sentinel_summary$expected_direction == "NEGATIVE" & is.finite(sentinel_summary$n_pairs) & sentinel_summary$n_pairs >= 5L & !is.na(sentinel_summary$direction_concordant), , drop = FALSE]
  pcon <- if (nrow(pos)) mean(pos$direction_concordant) else NA_real_
  ncon <- if (nrow(neg)) mean(neg$direction_concordant) else NA_real_
  status <- if (nrow(pos) >= 3L && nrow(neg) >= 2L && pcon <= 0.25 && ncon <= 0.25) "FAIL_POSSIBLE_GLOBAL_LABEL_INVERSION" else if (nrow(pos) + nrow(neg) < 5L) "REVIEW_INSUFFICIENT_SENTINELS" else "PASS_NO_GLOBAL_INVERSION_PATTERN"
  data.frame(n_positive_expected_evaluable = nrow(pos), positive_expected_concordance = pcon, n_negative_expected_evaluable = nrow(neg), negative_expected_concordance = ncon, status = status, interpretation = "Sentinel directions are a label-orientation QC, not a biological success gate.", stringsAsFactors = FALSE)
}

ol_p5_module_scores <- function(paired, modules, min_genes = 5L) {
  d <- paired[paired$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  coverage <- list(); scores <- list()
  for (m in names(modules)) {
    genes <- modules[[m]]
    measured <- intersect(genes, unique(d$gene))
    coverage[[m]] <- data.frame(module = m, n_catalog_genes = length(genes), n_measured_genes = length(measured), coverage_fraction = length(measured)/length(genes), min_genes_required = min_genes, status = ifelse(length(measured) >= min_genes, "PASS", "INSUFFICIENT_GENE_COVERAGE"), measured_genes = paste(measured, collapse = ";"), stringsAsFactors = FALSE)
    if (length(measured) < min_genes) next
    dm <- d[d$gene %in% measured, c("case_id","gene","paired_delta_tumor_minus_normal"), drop = FALSE]
    dm$gene_z <- ave(dm$paired_delta_tumor_minus_normal, dm$gene, FUN = function(x) {s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else as.numeric(scale(x))})
    by_case <- split(dm, dm$case_id)
    scores[[m]] <- ol_p5_bind_rows(lapply(by_case, function(z) data.frame(module = m, case_id = z$case_id[1], n_genes_available = sum(is.finite(z$gene_z)), module_score = if (sum(is.finite(z$gene_z)) >= min(3L, min_genes)) mean(z$gene_z, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)))
  }
  list(coverage = ol_p5_bind_rows(coverage), scores = ol_p5_bind_rows(scores))
}

ol_p5_target_module_cor <- function(paired, module_scores, primary_modules = c("CAF_Core","ECM_Remodeling","TGFb_Response")) {
  empty <- data.frame(module=character(), family=character(), n_cases=integer(), spearman_rho=numeric(), p_value=numeric(), analysis_role=character(), fdr_within_family=numeric(), stringsAsFactors=FALSE)
  if (!nrow(paired) || !nrow(module_scores) || !all(c("module","case_id","module_score") %in% names(module_scores))) return(empty)
  target <- paired[paired$gene == OLFML2B_PART5_TARGET & paired$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, c("case_id","paired_delta_tumor_minus_normal"), drop = FALSE]
  if (!nrow(target)) return(empty)
  names(target)[2] <- "OLFML2B_paired_delta"
  rows <- lapply(split(module_scores, module_scores$module), function(z) {
    d <- merge(target, z[, c("case_id","module_score")], by = "case_id", all = FALSE)
    d <- d[is.finite(d$OLFML2B_paired_delta) & is.finite(d$module_score), , drop = FALSE]
    ct <- if (nrow(d) >= 10L) tryCatch(suppressWarnings(stats::cor.test(d$OLFML2B_paired_delta, d$module_score, method = "spearman", exact = FALSE)), error = function(e) NULL) else NULL
    data.frame(module = z$module[1], family = ifelse(z$module[1] %in% primary_modules, "PRIMARY_CAF_ECM_TGFB", "CONTROL"), n_cases = nrow(d), spearman_rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_, p_value = if (!is.null(ct)) ct$p.value else NA_real_, analysis_role = "Exploratory patient-level protein co-variation", stringsAsFactors = FALSE)
  })
  out <- ol_p5_bind_rows(rows)
  out$fdr_within_family <- ave(out$p_value, out$family, FUN = function(x) stats::p.adjust(x, method = "BH"))
  out
}

ol_p5_save_figures <- function(dirs, paired, direction_tests, sentinel_summary, module_scores, module_cor) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(data.frame())
  figs <- list()
  d <- paired[paired$gene == OLFML2B_PART5_TARGET & paired$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  if (nrow(d)) {
    long <- rbind(
      data.frame(case_id=d$case_id, tissue="Normal", value=d$median_log_ratio_normal),
      data.frame(case_id=d$case_id, tissue="Tumor", value=d$median_log_ratio_tumor)
    )
    long$tissue <- factor(long$tissue, levels=c("Normal","Tumor"))
    p1 <- ggplot2::ggplot(long, ggplot2::aes(tissue, value, group=case_id)) + ggplot2::geom_line(alpha=.35) + ggplot2::geom_point(size=1.7) + ggplot2::theme_classic() + ggplot2::labs(title="OLFML2B protein abundance by paired case", x=NULL, y="Sample/reference log ratio")
    stem <- file.path(dirs$figures, "Figure_P5A_OLFML2B_case_paired_abundance")
    ggplot2::ggsave(paste0(stem,".pdf"), p1, width=6.2, height=4.8, useDingbats=FALSE)
    ggplot2::ggsave(paste0(stem,".png"), p1, width=6.2, height=4.8, dpi=600)
    figs[[length(figs)+1L]] <- data.frame(figure="P5A", stem=stem, status="PASS", stringsAsFactors=FALSE)
    d <- d[order(d$paired_delta_tumor_minus_normal), , drop=FALSE]
    d$case_id <- factor(d$case_id, levels=d$case_id)
    p2 <- ggplot2::ggplot(d, ggplot2::aes(case_id, paired_delta_tumor_minus_normal)) + ggplot2::geom_col() + ggplot2::geom_hline(yintercept=0, linewidth=.3) + ggplot2::coord_flip() + ggplot2::theme_classic() + ggplot2::labs(title="OLFML2B case-paired protein delta", x="Case", y="Tumor - Normal log ratio")
    stem <- file.path(dirs$figures, "Figure_P5B_OLFML2B_paired_delta_waterfall")
    ggplot2::ggsave(paste0(stem,".pdf"), p2, width=6.8, height=max(4.8, .18*nrow(d)+1.5), useDingbats=FALSE, limitsize=FALSE)
    ggplot2::ggsave(paste0(stem,".png"), p2, width=6.8, height=max(4.8, .18*nrow(d)+1.5), dpi=600, limitsize=FALSE)
    figs[[length(figs)+1L]] <- data.frame(figure="P5B", stem=stem, status="PASS", stringsAsFactors=FALSE)
  }
  if (nrow(sentinel_summary)) {
    s <- sentinel_summary[is.finite(sentinel_summary$median_delta), , drop=FALSE]
    if (nrow(s)) {
      s$gene <- factor(s$gene, levels=s$gene[order(s$median_delta)])
      p <- ggplot2::ggplot(s, ggplot2::aes(gene, median_delta)) + ggplot2::geom_col() + ggplot2::geom_hline(yintercept=0, linewidth=.3) + ggplot2::coord_flip() + ggplot2::facet_wrap(~group, scales="free_y") + ggplot2::theme_classic() + ggplot2::labs(title="Sentinel protein direction audit", x=NULL, y="Median paired delta")
      stem <- file.path(dirs$figures, "Figure_P5C_sentinel_direction_audit")
      ggplot2::ggsave(paste0(stem,".pdf"), p, width=8.2, height=7.0, useDingbats=FALSE)
      ggplot2::ggsave(paste0(stem,".png"), p, width=8.2, height=7.0, dpi=600)
      figs[[length(figs)+1L]] <- data.frame(figure="P5C", stem=stem, status="PASS", stringsAsFactors=FALSE)
    }
  }
  ol_p5_bind_rows(figs)
}

ol_p5_output_manifest <- function(dirs) {
  roots <- unlist(dirs[c("tables","figures","reports","qc","logs")], use.names=FALSE)
  files <- unique(unlist(lapply(roots, function(d) if (dir.exists(d)) list.files(d, recursive=TRUE, full.names=TRUE) else character()), use.names=FALSE))
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(data.frame())
  data.frame(path=normalizePath(files, winslash="/", mustWork=FALSE), size_bytes=file.info(files)$size, md5=unname(tools::md5sum(files)), stringsAsFactors=FALSE)
}

run_olfml2b_part5_pdc_production <- function(
  root = "D:/OLFML2B_STAD",
  pdc614_dir = file.path(root, "data", "PDC_STAD", "PDC000614_standardized"),
  output_subdir = "Part5",
  mapping_coverage_min = OLFML2B_PART5_MAPPING_COVERAGE_MIN,
  min_pairs_supportive = OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE,
  min_pairs_formal = OLFML2B_PART5_MIN_PAIRS_FORMAL,
  min_genes_per_module_score = 5L,
  bootstrap_B = OLFML2B_PART5_BOOTSTRAP_B,
  make_figures = TRUE,
  reset_part5_outputs = TRUE
) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  pdc614_dir <- normalizePath(pdc614_dir, winslash = "/", mustWork = FALSE)
  dirs <- ol_p5_make_dirs(root, output_subdir)
  if (isTRUE(reset_part5_outputs)) ol_p5_reset_outputs(dirs)
  log_file <- file.path(dirs$logs, "Part5_PDC000614_case_paired_production.log")
  ol_p5_log("INFO", "Starting Part5 | ", OLFML2B_PART5_VERSION, log_file=log_file)
  ol_p5_log("INFO", "Active protein cohort: ", OLFML2B_PART5_COHORT, log_file=log_file)
  ol_p5_log("INFO", "PDC000614 directory: ", pdc614_dir, log_file=log_file)

  params <- data.frame(
    parameter = c("version","root","pdc614_dir","target","primary_value_type","sensitivity_value_type","mapping_coverage_min","min_pairs_supportive","min_pairs_formal","bootstrap_B","inference_unit","claim_ceiling"),
    value = c(OLFML2B_PART5_VERSION,root,pdc614_dir,OLFML2B_PART5_TARGET,OLFML2B_PART5_PRIMARY_VALUE_TYPE,OLFML2B_PART5_SENSITIVITY_VALUE_TYPE,mapping_coverage_min,min_pairs_supportive,min_pairs_formal,bootstrap_B,"case-level paired Tumor-Normal delta after within-case/tissue median collapse","single-cohort orthogonal paired-protein support only"),
    stringsAsFactors=FALSE
  )
  ol_p5_write_csv(params, file.path(dirs$tables,"00_run_parameters.csv"))

  specs <- ol_p5_target_file_specs(pdc614_dir)
  ol_p5_write_csv(specs, file.path(dirs$tables,"01_pdc614_file_manifest.csv"))
  if (any(specs$required & !specs$present)) ol_p5_abort("Required PDC000614 files are missing. See 01_pdc614_file_manifest.csv")

  f <- setNames(specs$file, specs$component)
  sample_map <- ol_p5_read_tabular(f[["sample_map"]])
  label_map <- if (isTRUE(specs$present[specs$component=="label_map"])) ol_p5_read_tabular(f[["label_map"]]) else data.frame()
  biospec <- ol_p5_read_tabular(f[["biospecimen"]])
  header <- ol_p5_header_tokens(f[["tmt18"]])
  quant_map <- ol_p5_quant_map_from_header(header)
  schema <- data.frame(source=c("TMT18_header","sample_map","label_map","biospecimen"), n_rows=c(NA_integer_,nrow(sample_map),nrow(label_map),nrow(biospec)), n_cols=c(length(header),ncol(sample_map),ncol(label_map),ncol(biospec)), columns=c(paste(header,collapse=";"),paste(names(sample_map),collapse=";"),paste(names(label_map),collapse=";"),paste(names(biospec),collapse=";")), stringsAsFactors=FALSE)
  ol_p5_write_csv(schema, file.path(dirs$tables,"02_input_schema_audit.csv"))

  processed_detection <- ol_p5_bind_rows(lapply(c("tmt18","summary","primary_candidate"), function(nm) {
    z <- ol_p5_exact_gene_detection_first_field(f[[nm]], OLFML2B_PART5_TARGET)
    z$matrix_role <- nm
    z
  }))
  ol_p5_write_csv(processed_detection, file.path(dirs$tables,"03a_OLFML2B_processed_matrix_detection.csv"))

  modules <- ol_p5_catalog(); sentinels <- ol_p5_sentinel_catalog()
  genes <- unique(c(OLFML2B_PART5_TARGET, unlist(modules, use.names=FALSE), sentinels$gene))
  scan <- ol_p5_scan_selected_rows(f[["tmt18"]], genes, quant_map)
  ol_p5_write_csv(scan$mismatch, file.path(dirs$qc,"03b_selected_row_token_mismatch_audit.csv"))
  target_row_mismatch <- nrow(scan$mismatch) > 0L && any(vapply(scan$mismatch$first_field, function(z) OLFML2B_PART5_TARGET %in% ol_p5_gene_tokens(z), logical(1)))
  if (target_row_mismatch) ol_p5_abort("OLFML2B target row token count does not match the TMT18 header.")
  values <- scan$values
  if (!nrow(values)) {
    values <- data.frame(gene=character(), source_gene_field=character(), source_line_number=integer(), quant_column=character(), sample_id=character(), sample_id_norm=character(), value_type=character(), log_ratio=numeric(), stringsAsFactors=FALSE)
  }
  detected_genes <- sort(unique(values$gene))
  target_detected <- OLFML2B_PART5_TARGET %in% detected_genes
  detection <- data.frame(gene=genes, detected_in_tmt18=genes %in% detected_genes, n_source_rows=vapply(genes, function(g) length(unique(values$source_line_number[values$gene==g])), integer(1)), n_finite_values=vapply(genes, function(g) sum(is.finite(values$log_ratio[values$gene==g])), integer(1)), stringsAsFactors=FALSE)
  ol_p5_write_csv(detection, file.path(dirs$tables,"03c_target_and_module_detection_audit.csv"))
  if (!target_detected) ol_p5_log("WARN", "OLFML2B is not quantified in the PDC000614 TMT18 matrix.", log_file=log_file)

  sample_long <- ol_p5_flatten_sample_map(sample_map)
  mapped <- ol_p5_map_quant_to_biospecimen(quant_map, sample_long, biospec)
  mapping_summary <- ol_p5_mapping_summary(mapped)
  ol_p5_write_csv(mapped, file.path(dirs$tables,"04_quant_column_sample_biospecimen_map.csv"))
  ol_p5_write_csv(sample_long, file.path(dirs$tables,"05_sample_map_long.csv"))
  ol_p5_write_csv(mapping_summary, file.path(dirs$tables,"06_mapping_coverage_summary.csv"))
  map_lookup <- mapped[, c("quant_column","sample_id","sample_id_norm","value_type","sample_id_class","analytical_sample","channel","biospecimen_match_status","case_id","tissue_class","usable_for_tumor_normal"), drop=FALSE]
  values <- merge(values, map_lookup, by=c("quant_column","sample_id","sample_id_norm","value_type"), all.x=TRUE, sort=FALSE)
  ol_p5_write_csv(values[values$gene==OLFML2B_PART5_TARGET,,drop=FALSE], file.path(dirs$tables,"07_OLFML2B_reference_relative_values.csv"))

  case_tissue <- ol_p5_case_tissue_aggregate(values)
  paired <- ol_p5_pair_case_tissue(case_tissue)
  target_case <- case_tissue[case_tissue$gene==OLFML2B_PART5_TARGET,,drop=FALSE]
  target_paired <- paired[paired$gene==OLFML2B_PART5_TARGET,,drop=FALSE]
  ol_p5_write_csv(target_case, file.path(dirs$tables,"08_OLFML2B_case_tissue_aggregates.csv"))
  ol_p5_write_csv(target_paired, file.path(dirs$tables,"09_OLFML2B_case_paired_deltas.csv"))

  coverage <- mapping_summary$value[mapping_summary$metric=="primary_mapping_coverage"]
  mapping_pass <- length(coverage)==1L && is.finite(coverage) && coverage >= mapping_coverage_min
  tests <- ol_p5_bind_rows(lapply(c(OLFML2B_PART5_PRIMARY_VALUE_TYPE, OLFML2B_PART5_SENSITIVITY_VALUE_TYPE), function(v) ol_p5_direction_test(paired,v,mapping_pass,B=bootstrap_B)))
  tests$fdr_target_value_type_family <- stats::p.adjust(tests$wilcox_signed_rank_p, method="BH")
  ol_p5_write_csv(tests, file.path(dirs$tables,"10_OLFML2B_case_paired_direction_tests.csv"))

  sentinel_paired <- paired[paired$gene %in% sentinels$gene,,drop=FALSE]
  sentinel_summary <- ol_p5_sentinel_summary(paired,sentinels)
  orientation <- ol_p5_orientation_qc(sentinel_summary)
  ol_p5_write_csv(sentinel_paired, file.path(dirs$tables,"11_sentinel_case_paired_deltas.csv"))
  ol_p5_write_csv(sentinel_summary, file.path(dirs$tables,"12_sentinel_direction_summary.csv"))
  ol_p5_write_csv(orientation, file.path(dirs$tables,"13_label_orientation_qc.csv"))

  mod <- ol_p5_module_scores(paired,modules,min_genes=min_genes_per_module_score)
  modcor <- ol_p5_target_module_cor(paired,mod$scores)
  ol_p5_write_csv(mod$coverage, file.path(dirs$tables,"14_module_gene_coverage.csv"))
  ol_p5_write_csv(mod$scores, file.path(dirs$tables,"15_module_case_scores.csv"))
  ol_p5_write_csv(modcor, file.path(dirs$tables,"16_OLFML2B_module_correlations.csv"))

  primary_test <- tests[tests$value_type==OLFML2B_PART5_PRIMARY_VALUE_TYPE,,drop=FALSE]
  n_pairs <- if (nrow(primary_test)) primary_test$n_paired_cases[1] else 0L
  orientation_fail <- nrow(orientation) && orientation$status[1]=="FAIL_POSSIBLE_GLOBAL_LABEL_INVERSION"
  overall <- if (!target_detected) "NO_GO_TARGET_NOT_QUANTIFIED" else if (!mapping_pass) "HOLD_MAPPING_COVERAGE_BELOW_90_PERCENT" else if (orientation_fail) "HOLD_POSSIBLE_GLOBAL_LABEL_INVERSION" else if (n_pairs < min_pairs_supportive) "HOLD_INSUFFICIENT_CASE_PAIRS" else if (n_pairs < min_pairs_formal) "GO_SUPPORTIVE_CASE_PAIRED_PROTEIN_EVIDENCE" else "GO_FORMAL_SINGLE_COHORT_CASE_PAIRED_PROTEIN_SUPPORT"
  direction <- if (nrow(primary_test)) primary_test$evidence_direction[1] else "NOT_EVALUABLE"
  go <- data.frame(
    criterion=c("active_protein_cohort","required_files","exact_target_detection","primary_mapping_coverage","case_paired_inferential_units","primary_direction_test","sentinel_label_orientation","module_context","pipeline_continuation","claim_ceiling","overall_status"),
    status=c(OLFML2B_PART5_COHORT,ifelse(all(specs$present[specs$required]),"PASS","FAIL"),ifelse(target_detected,"PASS","NO_GO"),ifelse(mapping_pass,"PASS",paste0("FAIL_",round(100*coverage,1),"_PERCENT")),ifelse(n_pairs>=min_pairs_formal,"PASS_FORMAL",ifelse(n_pairs>=min_pairs_supportive,"PASS_SUPPORTIVE","FAIL_INSUFFICIENT")),ifelse(n_pairs>=min_pairs_supportive,paste0("EVALUABLE_",direction),"NOT_EVALUABLE"),orientation$status[1],ifelse(any(mod$coverage$status=="PASS"),"EVALUABLE","NOT_EVALUABLE"),"PASS_RESULT_INDEPENDENT","SINGLE_COHORT_ORTHOGONAL_CASE_PAIRED_PROTEIN_SUPPORT_ONLY",overall),
    boundary=c("Only the project-local standardized PDC000614 TMT18 cohort is active.","TMT18, sample map and biospecimen metadata are mandatory; label map is an audit cross-check.","OLFML2B is exact-token matched; OLFML2A and OLFM2 are never merged.","At least 90% of primary sample columns must map to case and Tumor/Normal tissue.","Technical/sample repeats are collapsed by median before case-level pairing.","Direction and P value are reported but never used to decide pipeline completion.","Sentinels detect a global label inversion pattern; they are not a target-success gate.","CAF/ECM/TGFb protein co-variation is exploratory and patient-level.","Null or discordant protein results remain reportable.","No causal, survival, diagnostic, treatment-selection or multi-cohort protein-validation claim.","Overall status is structural/evaluability based and independent of effect sign or significance."),
    stringsAsFactors=FALSE
  )
  ol_p5_write_csv(go, file.path(dirs$tables,"20_part5_pdc614_go_no_go_summary.csv"))

  amendment <- data.frame(
    amendment_date="2026-07-21", active_release=OLFML2B_PART5_VERSION,
    active_protein_scope="PDC000614_standardized only",
    amendment_reason="The active protein methodology was rebuilt around a target-quantified TMT18 cohort with auditable sample, biospecimen, tissue and case pairing.",
    transparency_rule="This is a documented protocol amendment made before interpreting the PDC000614 case-paired OLFML2B direction; all thresholds are frozen in code.",
    stringsAsFactors=FALSE
  )
  ol_p5_write_csv(amendment, file.path(dirs$tables,"21_protein_scope_method_amendment.csv"))

  figure_registry <- if (isTRUE(make_figures)) tryCatch(ol_p5_save_figures(dirs,paired,tests,sentinel_summary,mod$scores,modcor), error=function(e) data.frame(figure="ALL",stem=NA_character_,status=paste0("WARN: ",conditionMessage(e)),stringsAsFactors=FALSE)) else data.frame()
  ol_p5_write_csv(figure_registry, file.path(dirs$tables,"98_figure_registry.csv"))
  ol_p5_write_csv(ol_p5_output_manifest(dirs), file.path(dirs$tables,"99_output_file_manifest.csv"))

  index <- list(
    version=OLFML2B_PART5_VERSION, generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
    active_cohort=OLFML2B_PART5_COHORT, root=root, pdc614_dir=pdc614_dir, dirs=dirs,
    params=params, file_manifest=specs, input_schema=schema, detection=detection,
    quant_mapping=mapped, mapping_summary=mapping_summary,
    target_reference_relative_values=values[values$gene==OLFML2B_PART5_TARGET,,drop=FALSE],
    target_case_tissue=target_case, target_paired=target_paired, direction_tests=tests,
    sentinel_summary=sentinel_summary, label_orientation_qc=orientation,
    module_coverage=mod$coverage, module_scores=mod$scores, module_correlations=modcor,
    go_no_go=go, overall_status=overall, final_gene_lock=FALSE,
    claim_boundary="PDC000614 provides single-cohort case-paired orthogonal protein evidence only."
  )
  ol_p5_save_rds(index,file.path(dirs$objects,"Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds"))
  ol_p5_save_rds(index,file.path(dirs$objects,"Part5_PDC_production_index.rds"))
  ol_p5_log("INFO","Part5 completed | overall_status=",overall," | paired_cases=",n_pairs," | direction=",direction,log_file=log_file)
  invisible(index)
}

# ==============================================================================
# Part5 audit hardening v1.2.0
# Adds dual-denominator mapping coverage, label-reagent cross-check, target
# measurement/attrition audit, mandatory within-plex Tumor-Normal pairing,
# module-specific coverage thresholds, and separate structural versus biological
# evidence states.
# ==============================================================================
OLFML2B_PART5_VERSION <- "v1.2.0_20260722_PDC614_WITHIN_PLEX_AUDIT_HARDENING"

# Use the thresholds passed to the production runner rather than hidden global
# constants, so the audit table and the inferential role cannot drift apart.
ol_p5_direction_test <- function(
  paired,
  value_type,
  mapping_pass,
  B = OLFML2B_PART5_BOOTSTRAP_B,
  min_pairs_supportive = OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE,
  min_pairs_formal = OLFML2B_PART5_MIN_PAIRS_FORMAL
) {
  d <- paired[paired$gene == OLFML2B_PART5_TARGET & paired$value_type == value_type, , drop = FALSE]
  x <- d$paired_delta_tumor_minus_normal
  x <- x[is.finite(x)]
  n <- length(x)
  ci <- if (mapping_pass) {
    ol_p5_bootstrap_median_ci(
      x,
      B = B,
      seed = OLFML2B_PART5_SEED + ifelse(value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, 0L, 1L)
    )
  } else {
    c(low = NA_real_, high = NA_real_)
  }
  wp <- if (mapping_pass && n >= min_pairs_supportive) {
    tryCatch(stats::wilcox.test(x, mu = 0, alternative = "two.sided", exact = FALSE)$p.value, error = function(e) NA_real_)
  } else NA_real_
  nz <- x[x != 0]
  sp <- if (mapping_pass && length(nz) >= min_pairs_supportive) {
    stats::binom.test(sum(nz > 0), length(nz), p = 0.5, alternative = "two.sided")$p.value
  } else NA_real_
  role <- if (!mapping_pass) {
    "INFERENCE_BLOCKED_MAPPING"
  } else if (value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE && n >= min_pairs_formal) {
    "FORMAL_PRIMARY_SINGLE_COHORT"
  } else if (n >= min_pairs_supportive) {
    "SUPPORTIVE_SENSITIVITY"
  } else {
    "DESCRIPTIVE_INSUFFICIENT_PAIRS"
  }
  direction <- if (!n) {
    "NOT_EVALUABLE"
  } else if (is.finite(ci[1L]) && ci[1L] > 0) {
    "POSITIVE_CI_EXCLUDES_ZERO"
  } else if (is.finite(ci[2L]) && ci[2L] < 0) {
    "NEGATIVE_CI_EXCLUDES_ZERO"
  } else if (stats::median(x) > 0) {
    "POSITIVE_UNCERTAIN"
  } else if (stats::median(x) < 0) {
    "NEGATIVE_UNCERTAIN"
  } else {
    "NULL_MEDIAN"
  }
  data.frame(
    gene = OLFML2B_PART5_TARGET, cohort = OLFML2B_PART5_COHORT,
    value_type = value_type, n_paired_cases = n,
    median_delta = if (n) stats::median(x) else NA_real_,
    mean_delta = if (n) mean(x) else NA_real_,
    bootstrap_median_ci_low = unname(ci[1L]),
    bootstrap_median_ci_high = unname(ci[2L]),
    positive_fraction = if (n) mean(x > 0) else NA_real_,
    rank_biserial_sign = if (length(nz)) (sum(nz > 0) - sum(nz < 0)) / length(nz) else NA_real_,
    wilcox_signed_rank_p = wp, exact_sign_p = sp,
    min_pairs_supportive = min_pairs_supportive,
    min_pairs_formal = min_pairs_formal,
    inference_role = role, evidence_direction = direction,
    claim_boundary = "Case-paired single-cohort protein direction; no causality, prognosis, or clinical utility claim.",
    stringsAsFactors = FALSE
  )
}

ol_p5_mapping_summary <- function(mapped) {
  primary <- mapped[mapped$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  eligible <- primary$sample_id_class == "SAMPLE"
  n_all <- nrow(primary)
  n_eligible <- sum(eligible, na.rm = TRUE)
  n_usable <- sum(primary$usable_for_tumor_normal & eligible, na.rm = TRUE)
  data.frame(
    metric = c(
      "n_quant_columns_all_value_types", "n_primary_quant_columns_all", "n_primary_eligible_sample_columns",
      "n_primary_usable_tumor_normal_columns", "eligible_sample_mapping_coverage",
      "all_primary_column_usable_fraction", "primary_mapping_coverage",
      "n_primary_excluded_reference_or_internal_columns", "n_primary_tumor_columns",
      "n_primary_normal_columns", "n_unique_mapped_cases"
    ),
    value = c(
      nrow(mapped), n_all, n_eligible, n_usable,
      if (n_eligible > 0) n_usable / n_eligible else NA_real_,
      if (n_all > 0) n_usable / n_all else NA_real_,
      if (n_eligible > 0) n_usable / n_eligible else NA_real_,
      n_all - n_eligible,
      sum(primary$usable_for_tumor_normal & primary$tissue_class == "Tumor", na.rm = TRUE),
      sum(primary$usable_for_tumor_normal & primary$tissue_class == "Normal", na.rm = TRUE),
      length(unique(primary$case_id[primary$usable_for_tumor_normal]))
    ),
    denominator = c(
      "all value-type columns", "all primary Log_Ratio columns", "eligible biological SAMPLE columns",
      "eligible biological SAMPLE columns", "eligible biological SAMPLE columns",
      "all primary Log_Ratio columns", "eligible biological SAMPLE columns",
      "all primary Log_Ratio columns", "usable primary columns", "usable primary columns", "usable primary columns"
    ),
    stringsAsFactors = FALSE
  )
}

ol_p5_label_reagent_crosscheck <- function(sample_map, label_map) {
  empty <- data.frame(label_reagent = character(), normalized_reagent = character(), linked_to_label_definition = logical(), match_location = character(), status = character(), stringsAsFactors = FALSE)
  if (!nrow(sample_map)) return(transform(empty, status = character()))
  rc <- names(sample_map)[grepl("label.*reagent|reagent.*label", names(sample_map), ignore.case = TRUE)]
  if (!length(rc)) return(data.frame(label_reagent = NA_character_, normalized_reagent = NA_character_, linked_to_label_definition = NA, match_location = NA_character_, status = "NOT_EVALUABLE_SAMPLE_MAP_HAS_NO_LABEL_REAGENT_COLUMN", stringsAsFactors = FALSE))
  reagents <- unique(ol_p5_clean_text(sample_map[[rc[1L]]]))
  reagents <- reagents[!is.na(reagents) & nzchar(reagents)]
  if (!length(reagents)) return(data.frame(label_reagent = NA_character_, normalized_reagent = NA_character_, linked_to_label_definition = NA, match_location = NA_character_, status = "NOT_EVALUABLE_NO_LABEL_REAGENT_VALUE", stringsAsFactors = FALSE))
  if (!nrow(label_map)) return(data.frame(label_reagent = reagents, normalized_reagent = ol_p5_norm_id(reagents), linked_to_label_definition = FALSE, match_location = NA_character_, status = "LABEL_MAP_NOT_AVAILABLE", stringsAsFactors = FALSE))
  header_norm <- ol_p5_norm_id(names(label_map))
  value_vec <- unlist(lapply(label_map, as.character), use.names = FALSE)
  value_norm <- ol_p5_norm_id(value_vec)
  rows <- lapply(reagents, function(r) {
    rn <- ol_p5_norm_id(r)
    hm <- which(header_norm == rn | (nzchar(rn) & grepl(rn, header_norm, fixed = TRUE)))
    vm <- which(value_norm == rn | (nzchar(rn) & grepl(rn, value_norm, fixed = TRUE)))
    linked <- length(hm) > 0L || length(vm) > 0L
    data.frame(label_reagent = r, normalized_reagent = rn, linked_to_label_definition = linked, match_location = if (length(hm)) "label_map_header" else if (length(vm)) "label_map_value" else NA_character_, status = ifelse(linked, "PASS_REAGENT_DEFINITION_LINKED", "REVIEW_REAGENT_DEFINITION_NOT_LINKED"), stringsAsFactors = FALSE)
  })
  ol_p5_bind_rows(rows)
}

ol_p5_case_plex_tissue_aggregate <- function(values) {
  d <- values[values$usable_for_tumor_normal & is.finite(values$log_ratio) & !is.na(values$analytical_sample) & nzchar(values$analytical_sample), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  keys <- interaction(d$gene, d$value_type, d$case_id, d$analytical_sample, d$tissue_class, drop = TRUE, lex.order = TRUE)
  ol_p5_bind_rows(lapply(split(d, keys), function(z) data.frame(
    gene = z$gene[1L], value_type = z$value_type[1L], case_id = z$case_id[1L], analytical_sample = z$analytical_sample[1L], tissue_class = z$tissue_class[1L],
    n_channels = nrow(z), median_log_ratio = stats::median(z$log_ratio, na.rm = TRUE), mean_log_ratio = mean(z$log_ratio, na.rm = TRUE), iqr_log_ratio = stats::IQR(z$log_ratio, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
}

ol_p5_pair_within_plex <- function(plex_tissue) {
  if (!nrow(plex_tissue)) return(data.frame())
  tum <- plex_tissue[plex_tissue$tissue_class == "Tumor", , drop = FALSE]
  nor <- plex_tissue[plex_tissue$tissue_class == "Normal", , drop = FALSE]
  if (!nrow(tum) || !nrow(nor)) return(data.frame())
  metrics <- c("n_channels", "median_log_ratio", "mean_log_ratio", "iqr_log_ratio")
  names(tum)[names(tum) %in% metrics] <- paste0(names(tum)[names(tum) %in% metrics], "_tumor")
  names(nor)[names(nor) %in% metrics] <- paste0(names(nor)[names(nor) %in% metrics], "_normal")
  out <- merge(tum, nor, by = c("gene", "value_type", "case_id", "analytical_sample"), all = FALSE, sort = FALSE)
  out$paired_delta_tumor_minus_normal <- out$median_log_ratio_tumor - out$median_log_ratio_normal
  out$pairing_scope <- "WITHIN_ANALYTICAL_TMT_PLEX"
  out
}

ol_p5_collapse_case_pairs <- function(plex_pairs) {
  if (!nrow(plex_pairs)) return(data.frame())
  keys <- interaction(plex_pairs$gene, plex_pairs$value_type, plex_pairs$case_id, drop = TRUE, lex.order = TRUE)
  ol_p5_bind_rows(lapply(split(plex_pairs, keys), function(z) data.frame(
    gene = z$gene[1L], value_type = z$value_type[1L], case_id = z$case_id[1L],
    n_same_plex_pairs = nrow(z), analytical_samples = paste(sort(unique(z$analytical_sample)), collapse = ";"),
    n_channels_tumor = sum(z$n_channels_tumor, na.rm = TRUE), n_channels_normal = sum(z$n_channels_normal, na.rm = TRUE),
    median_log_ratio_tumor = stats::median(z$median_log_ratio_tumor, na.rm = TRUE),
    median_log_ratio_normal = stats::median(z$median_log_ratio_normal, na.rm = TRUE),
    paired_delta_tumor_minus_normal = stats::median(z$paired_delta_tumor_minus_normal, na.rm = TRUE),
    mean_paired_delta_across_plexes = mean(z$paired_delta_tumor_minus_normal, na.rm = TRUE),
    pairing_scope = "CASE_LEVEL_COLLAPSE_OF_WITHIN_PLEX_DELTAS",
    stringsAsFactors = FALSE
  )))
}

ol_p5_target_attrition_audit <- function(mapped, values, plex_pairs) {
  pm <- mapped[mapped$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & mapped$usable_for_tumor_normal, , drop = FALSE]
  if (!nrow(pm)) return(list(attrition = data.frame(), tissue = data.frame(), same_plex = data.frame()))
  potential <- unique(pm[, c("case_id", "tissue_class"), drop = FALSE])
  cases <- sort(unique(potential$case_id))
  tv <- values[values$gene == OLFML2B_PART5_TARGET & values$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & values$usable_for_tumor_normal & is.finite(values$log_ratio), , drop = FALSE]
  attrition <- ol_p5_bind_rows(lapply(cases, function(cid) {
    has_tumor_map <- any(potential$case_id == cid & potential$tissue_class == "Tumor")
    has_normal_map <- any(potential$case_id == cid & potential$tissue_class == "Normal")
    tumor_finite <- any(tv$case_id == cid & tv$tissue_class == "Tumor")
    normal_finite <- any(tv$case_id == cid & tv$tissue_class == "Normal")
    same_plex <- any(plex_pairs$case_id == cid & plex_pairs$gene == OLFML2B_PART5_TARGET & plex_pairs$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE)
    data.frame(case_id = cid, mapped_tumor = has_tumor_map, mapped_normal = has_normal_map, target_finite_tumor = tumor_finite, target_finite_normal = normal_finite, same_plex_target_pair = same_plex,
      analysis_status = if (!has_tumor_map || !has_normal_map) "NOT_A_MAPPED_TUMOR_NORMAL_CASE" else if (!tumor_finite && !normal_finite) "TARGET_MISSING_BOTH_TISSUES" else if (!tumor_finite) "TARGET_MISSING_TUMOR" else if (!normal_finite) "TARGET_MISSING_NORMAL" else if (!same_plex) "TARGET_BOTH_TISSUES_BUT_NO_SAME_PLEX_PAIR" else "INCLUDED_SAME_PLEX_CASE_PAIR", stringsAsFactors = FALSE)
  }))
  tissue <- ol_p5_bind_rows(lapply(c("Tumor", "Normal"), function(tt) {
    den <- unique(pm$quant_column[pm$tissue_class == tt])
    num <- unique(tv$quant_column[tv$tissue_class == tt])
    data.frame(tissue_class = tt, n_mapped_primary_columns = length(den), n_finite_OLFML2B_columns = length(num), measurement_completeness = if (length(den)) length(num) / length(den) else NA_real_, interpretation = "Missing target quantification is not imputed as zero and is not biological absence.", stringsAsFactors = FALSE)
  }))
  same_plex <- plex_pairs[plex_pairs$gene == OLFML2B_PART5_TARGET & plex_pairs$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, c("case_id", "analytical_sample", "n_channels_tumor", "n_channels_normal", "paired_delta_tumor_minus_normal", "pairing_scope"), drop = FALSE]
  summary <- data.frame(
    step = c(
      "cases_with_any_mapped_primary_sample", "cases_with_mapped_tumor",
      "cases_with_mapped_normal", "mapped_tumor_normal_cases",
      "cases_with_finite_target_tumor", "cases_with_finite_target_normal",
      "cases_with_finite_target_both_tissues", "included_same_plex_target_pairs"
    ),
    n_cases = c(
      length(cases),
      sum(vapply(cases, function(cid) any(potential$case_id == cid & potential$tissue_class == "Tumor"), logical(1))),
      sum(vapply(cases, function(cid) any(potential$case_id == cid & potential$tissue_class == "Normal"), logical(1))),
      sum(attrition$mapped_tumor & attrition$mapped_normal),
      sum(attrition$target_finite_tumor),
      sum(attrition$target_finite_normal),
      sum(attrition$target_finite_tumor & attrition$target_finite_normal),
      sum(attrition$same_plex_target_pair)
    ),
    interpretation = c(
      "All cases represented by at least one eligible mapped primary sample column.",
      "Cases with at least one mapped Tumor sample.",
      "Cases with at least one mapped Normal sample.",
      "Potential case-paired denominator before target-specific missingness.",
      "Cases with at least one finite OLFML2B Tumor value.",
      "Cases with at least one finite OLFML2B Normal value.",
      "Cases with finite target values in both tissues, irrespective of plex.",
      "Final inferential cases with Tumor and Normal paired within the same analytical TMT plex."
    ),
    stringsAsFactors = FALSE
  )
  list(attrition = attrition, attrition_summary = summary, tissue = tissue, same_plex = same_plex)
}

ol_p5_sentinel_catalog <- function() {
  data.frame(
    gene = c("EPCAM","KRT8","KRT18","KRT19","MKI67","FAP","ACTA2","COL1A1","COL1A2","COL3A1","SPARC","TGFBI","SERPINE1","CTGF","INHBA","PGC","PGA3","TFF1","TFF2","CD8A","CD8B","GZMB","NKG7"),
    group = c(rep("Tumor_Epithelial",5), rep("CAF_ECM",6), rep("TGFb",4), rep("Gastric_Normal",4), rep("Immune_Control",4)),
    expected_direction = c(rep("POSITIVE",6), "NO_FIXED_EXPECTATION", rep("POSITIVE",8), rep("NEGATIVE",4), rep("NO_FIXED_EXPECTATION",4)),
    rationale = c(rep("Expected tumor-positive orientation sentinel",6), "ACTA2 is context-dependent because normal gastric smooth muscle can be abundant", rep("Expected tumor-positive orientation sentinel",8), rep("Expected normal-gastric-positive orientation sentinel",4), rep("No fixed direction; contextual immune control",4)),
    stringsAsFactors = FALSE
  )
}

ol_p5_module_scores <- function(paired, modules, min_genes = 3L, min_fraction = 0.60) {
  d <- paired[paired$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  coverage <- list(); scores <- list()
  for (m in names(modules)) {
    genes <- unique(modules[[m]])
    measured <- intersect(genes, unique(d$gene))
    required <- min(length(genes), max(as.integer(min_genes), ceiling(min_fraction * length(genes))))
    status <- ifelse(length(measured) >= required, "PASS", "INSUFFICIENT_GENE_COVERAGE")
    coverage[[m]] <- data.frame(module = m, n_catalog_genes = length(genes), n_measured_genes = length(measured), coverage_fraction = length(measured) / length(genes), min_genes_required = required, min_fraction_required = min_fraction, status = status, measured_genes = paste(measured, collapse = ";"), stringsAsFactors = FALSE)
    if (length(measured) < required) next
    dm <- d[d$gene %in% measured, c("case_id", "gene", "paired_delta_tumor_minus_normal"), drop = FALSE]
    dm$gene_z <- ave(dm$paired_delta_tumor_minus_normal, dm$gene, FUN = function(x) { s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else as.numeric(scale(x)) })
    by_case <- split(dm, dm$case_id)
    case_required <- max(3L, ceiling(min_fraction * length(measured)))
    scores[[m]] <- ol_p5_bind_rows(lapply(by_case, function(z) data.frame(module = m, case_id = z$case_id[1L], n_genes_available = sum(is.finite(z$gene_z)), min_case_genes_required = case_required, module_score = if (sum(is.finite(z$gene_z)) >= case_required) mean(z$gene_z, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)))
  }
  list(coverage = ol_p5_bind_rows(coverage), scores = ol_p5_bind_rows(scores))
}

run_olfml2b_part5_pdc_production <- function(
  root = "D:/OLFML2B_STAD",
  pdc614_dir = file.path(root, "data", "PDC_STAD", "PDC000614_standardized"),
  output_subdir = "Part5",
  mapping_coverage_min = OLFML2B_PART5_MAPPING_COVERAGE_MIN,
  min_pairs_supportive = OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE,
  min_pairs_formal = OLFML2B_PART5_MIN_PAIRS_FORMAL,
  min_genes_per_module_score = 3L,
  bootstrap_B = OLFML2B_PART5_BOOTSTRAP_B,
  make_figures = TRUE,
  reset_part5_outputs = TRUE
) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  pdc614_dir <- normalizePath(pdc614_dir, winslash = "/", mustWork = FALSE)
  dirs <- ol_p5_make_dirs(root, output_subdir)
  if (isTRUE(reset_part5_outputs)) ol_p5_reset_outputs(dirs)
  log_file <- file.path(dirs$logs, "Part5_PDC000614_case_paired_production.log")
  ol_p5_log("INFO", "Starting Part5 | ", OLFML2B_PART5_VERSION, log_file = log_file)
  ol_p5_log("INFO", "Active protein cohort: ", OLFML2B_PART5_COHORT, log_file = log_file)
  ol_p5_log("INFO", "PDC000614 directory: ", pdc614_dir, log_file = log_file)

  params <- data.frame(parameter = c("version","root","pdc614_dir","target","primary_value_type","sensitivity_value_type","mapping_coverage_min","min_pairs_supportive","min_pairs_formal","bootstrap_B","inference_unit","pairing_rule","claim_ceiling"), value = c(OLFML2B_PART5_VERSION, root, pdc614_dir, OLFML2B_PART5_TARGET, OLFML2B_PART5_PRIMARY_VALUE_TYPE, OLFML2B_PART5_SENSITIVITY_VALUE_TYPE, mapping_coverage_min, min_pairs_supportive, min_pairs_formal, bootstrap_B, "case-level paired Tumor-Normal delta", "Tumor and Normal must pair within the same analytical TMT plex before case collapse", "single-cohort orthogonal paired-protein support only"), stringsAsFactors = FALSE)
  ol_p5_write_csv(params, file.path(dirs$tables, "00_run_parameters.csv"))

  specs <- ol_p5_target_file_specs(pdc614_dir)
  ol_p5_write_csv(specs, file.path(dirs$tables, "01_pdc614_file_manifest.csv"))
  if (!all(specs$present[specs$required])) ol_p5_abort("One or more mandatory PDC000614 files are missing.")
  f <- setNames(specs$file, specs$component)
  header <- ol_p5_header_tokens(f[["tmt18"]])
  quant_map <- ol_p5_quant_map_from_header(header)
  sample_map <- ol_p5_read_tabular(f[["sample_map"]])
  label_map <- if (isTRUE(specs$present[specs$component == "label_map"])) ol_p5_read_tabular(f[["label_map"]]) else data.frame()
  biospec <- ol_p5_read_tabular(f[["biospecimen"]])
  schema <- data.frame(source = c("TMT18_header","sample_map","label_map","biospecimen"), n_rows = c(NA_integer_, nrow(sample_map), nrow(label_map), nrow(biospec)), n_cols = c(length(header), ncol(sample_map), ncol(label_map), ncol(biospec)), columns = c(paste(header, collapse = ";"), paste(names(sample_map), collapse = ";"), paste(names(label_map), collapse = ";"), paste(names(biospec), collapse = ";")), stringsAsFactors = FALSE)
  ol_p5_write_csv(schema, file.path(dirs$tables, "02_input_schema_audit.csv"))
  label_crosscheck <- ol_p5_label_reagent_crosscheck(sample_map, label_map)
  ol_p5_write_csv(label_crosscheck, file.path(dirs$tables, "02a_label_reagent_crosscheck.csv"))

  processed_detection <- ol_p5_bind_rows(lapply(c("tmt18","summary","primary_candidate"), function(nm) { z <- ol_p5_exact_gene_detection_first_field(f[[nm]], OLFML2B_PART5_TARGET); z$matrix_role <- nm; z }))
  ol_p5_write_csv(processed_detection, file.path(dirs$tables, "03a_OLFML2B_processed_matrix_detection.csv"))
  modules <- ol_p5_catalog(); sentinels <- ol_p5_sentinel_catalog()
  genes <- unique(c(OLFML2B_PART5_TARGET, unlist(modules, use.names = FALSE), sentinels$gene))
  scan <- ol_p5_scan_selected_rows(f[["tmt18"]], genes, quant_map)
  ol_p5_write_csv(scan$mismatch, file.path(dirs$qc, "03b_selected_row_token_mismatch_audit.csv"))
  target_row_mismatch <- nrow(scan$mismatch) > 0L && any(vapply(scan$mismatch$first_field, function(z) OLFML2B_PART5_TARGET %in% ol_p5_gene_tokens(z), logical(1)))
  if (target_row_mismatch) ol_p5_abort("OLFML2B target row token count does not match the TMT18 header.")
  values <- scan$values
  if (!nrow(values)) values <- data.frame(gene = character(), source_gene_field = character(), source_line_number = integer(), quant_column = character(), sample_id = character(), sample_id_norm = character(), value_type = character(), log_ratio = numeric(), stringsAsFactors = FALSE)
  detected_genes <- sort(unique(values$gene)); target_detected <- OLFML2B_PART5_TARGET %in% detected_genes
  detection <- data.frame(gene = genes, detected_in_tmt18 = genes %in% detected_genes, n_source_rows = vapply(genes, function(g) length(unique(values$source_line_number[values$gene == g])), integer(1)), n_finite_values = vapply(genes, function(g) sum(is.finite(values$log_ratio[values$gene == g])), integer(1)), stringsAsFactors = FALSE)
  ol_p5_write_csv(detection, file.path(dirs$tables, "03c_target_and_module_detection_audit.csv"))

  sample_long <- ol_p5_flatten_sample_map(sample_map)
  mapped <- ol_p5_map_quant_to_biospecimen(quant_map, sample_long, biospec)
  mapping_summary <- ol_p5_mapping_summary(mapped)
  ol_p5_write_csv(mapped, file.path(dirs$tables, "04_quant_column_sample_biospecimen_map.csv"))
  ol_p5_write_csv(sample_long, file.path(dirs$tables, "05_sample_map_long.csv"))
  ol_p5_write_csv(mapping_summary, file.path(dirs$tables, "06_mapping_coverage_summary.csv"))
  map_lookup <- mapped[, c("quant_column","sample_id","sample_id_norm","value_type","sample_id_class","analytical_sample","channel","biospecimen_match_status","case_id","tissue_class","usable_for_tumor_normal"), drop = FALSE]
  values <- merge(values, map_lookup, by = c("quant_column","sample_id","sample_id_norm","value_type"), all.x = TRUE, sort = FALSE)
  ol_p5_write_csv(values[values$gene == OLFML2B_PART5_TARGET, , drop = FALSE], file.path(dirs$tables, "07_OLFML2B_reference_relative_values.csv"))

  plex_tissue <- ol_p5_case_plex_tissue_aggregate(values)
  plex_pairs <- ol_p5_pair_within_plex(plex_tissue)
  paired <- ol_p5_collapse_case_pairs(plex_pairs)
  target_case <- plex_tissue[plex_tissue$gene == OLFML2B_PART5_TARGET, , drop = FALSE]
  target_plex_pairs <- plex_pairs[plex_pairs$gene == OLFML2B_PART5_TARGET, , drop = FALSE]
  target_paired <- paired[paired$gene == OLFML2B_PART5_TARGET, , drop = FALSE]
  ol_p5_write_csv(target_case, file.path(dirs$tables, "08_OLFML2B_case_plex_tissue_aggregates.csv"))
  ol_p5_write_csv(target_case, file.path(dirs$tables, "08_OLFML2B_case_tissue_aggregates.csv"))  # compatibility alias; analytical plex remains explicit
  ol_p5_write_csv(target_paired, file.path(dirs$tables, "09_OLFML2B_case_paired_deltas.csv"))
  attr <- ol_p5_target_attrition_audit(mapped, values, plex_pairs)
  ol_p5_write_csv(attr$attrition, file.path(dirs$tables, "09a_OLFML2B_pair_attrition_audit.csv"))
  ol_p5_write_csv(attr$attrition_summary, file.path(dirs$tables, "09a2_OLFML2B_pair_attrition_summary.csv"))
  ol_p5_write_csv(attr$tissue, file.path(dirs$tables, "09b_OLFML2B_tissue_measurement_completeness.csv"))
  ol_p5_write_csv(attr$same_plex, file.path(dirs$tables, "09c_OLFML2B_same_plex_pairing_audit.csv"))

  coverage <- mapping_summary$value[mapping_summary$metric == "eligible_sample_mapping_coverage"]
  all_fraction <- mapping_summary$value[mapping_summary$metric == "all_primary_column_usable_fraction"]
  mapping_pass <- length(coverage) == 1L && is.finite(coverage) && coverage >= mapping_coverage_min
  tests <- ol_p5_bind_rows(lapply(
    c(OLFML2B_PART5_PRIMARY_VALUE_TYPE, OLFML2B_PART5_SENSITIVITY_VALUE_TYPE),
    function(v) ol_p5_direction_test(
      paired, v, mapping_pass, B = bootstrap_B,
      min_pairs_supportive = min_pairs_supportive,
      min_pairs_formal = min_pairs_formal
    )
  ))
  tests$fdr_target_value_type_family <- stats::p.adjust(tests$wilcox_signed_rank_p, method = "BH")
  ol_p5_write_csv(tests, file.path(dirs$tables, "10_OLFML2B_case_paired_direction_tests.csv"))

  sentinel_paired <- paired[paired$gene %in% sentinels$gene, , drop = FALSE]
  sentinel_summary <- ol_p5_sentinel_summary(paired, sentinels)
  orientation <- ol_p5_orientation_qc(sentinel_summary)
  ol_p5_write_csv(sentinel_paired, file.path(dirs$tables, "11_sentinel_case_paired_deltas.csv"))
  ol_p5_write_csv(sentinel_summary, file.path(dirs$tables, "12_sentinel_direction_summary.csv"))
  ol_p5_write_csv(orientation, file.path(dirs$tables, "13_label_orientation_qc.csv"))

  mod <- ol_p5_module_scores(paired, modules, min_genes = min_genes_per_module_score, min_fraction = 0.60)
  modcor <- ol_p5_target_module_cor(paired, mod$scores)
  ol_p5_write_csv(mod$coverage, file.path(dirs$tables, "14_module_gene_coverage.csv"))
  ol_p5_write_csv(mod$scores, file.path(dirs$tables, "15_module_case_scores.csv"))
  ol_p5_write_csv(modcor, file.path(dirs$tables, "16_OLFML2B_module_correlations.csv"))

  primary_test <- tests[tests$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, , drop = FALSE]
  n_pairs <- if (nrow(primary_test)) primary_test$n_paired_cases[1L] else 0L
  direction <- if (nrow(primary_test)) primary_test$evidence_direction[1L] else "NOT_EVALUABLE"
  orientation_fail <- nrow(orientation) && orientation$status[1L] == "FAIL_POSSIBLE_GLOBAL_LABEL_INVERSION"
  structural_status <- if (!target_detected) "NO_GO_TARGET_NOT_QUANTIFIED" else if (!mapping_pass) "HOLD_MAPPING_COVERAGE_BELOW_THRESHOLD" else if (orientation_fail) "HOLD_POSSIBLE_GLOBAL_LABEL_INVERSION" else if (n_pairs < min_pairs_supportive) "HOLD_INSUFFICIENT_SAME_PLEX_CASE_PAIRS" else if (n_pairs < min_pairs_formal) "GO_SUPPORTIVE_CASE_PAIRED_ANALYSIS" else "GO_FORMAL_CASE_PAIRED_ANALYSIS"
  biological_status <- if (n_pairs < min_pairs_supportive) "NOT_EVALUABLE" else switch(direction, POSITIVE_CI_EXCLUDES_ZERO = "DIRECTIONALLY_POSITIVE_CI_EXCLUDES_ZERO", NEGATIVE_CI_EXCLUDES_ZERO = "DIRECTIONALLY_NEGATIVE_CI_EXCLUDES_ZERO", POSITIVE_UNCERTAIN = "DIRECTIONALLY_POSITIVE_BUT_UNCERTAIN", NEGATIVE_UNCERTAIN = "DIRECTIONALLY_NEGATIVE_BUT_UNCERTAIN", NULL_MEDIAN = "NULL_MEDIAN", "NOT_EVALUABLE")
  label_status <- if (!nrow(label_crosscheck)) "NOT_EVALUABLE" else if (all(label_crosscheck$status == "PASS_REAGENT_DEFINITION_LINKED")) "PASS" else if (any(label_crosscheck$status == "LABEL_MAP_NOT_AVAILABLE")) "NOT_AVAILABLE_OPTIONAL" else "REVIEW_OPTIONAL_CROSSCHECK"
  go <- data.frame(
    criterion = c("active_protein_cohort","required_files","exact_target_detection","eligible_sample_mapping_coverage","all_primary_column_usable_fraction","label_reagent_crosscheck","same_plex_case_paired_units","analysis_evaluability_status","biological_evidence_status","sentinel_label_orientation","module_context","pipeline_continuation","claim_ceiling"),
    status = c(OLFML2B_PART5_COHORT, ifelse(all(specs$present[specs$required]), "PASS", "FAIL"), ifelse(target_detected, "PASS", "NO_GO"), ifelse(mapping_pass, paste0("PASS_", round(100 * coverage, 1), "_PERCENT"), paste0("FAIL_", round(100 * coverage, 1), "_PERCENT")), paste0(round(100 * all_fraction, 1), "_PERCENT"), label_status, ifelse(n_pairs >= min_pairs_formal, "PASS_FORMAL", ifelse(n_pairs >= min_pairs_supportive, "PASS_SUPPORTIVE", "FAIL_INSUFFICIENT")), structural_status, biological_status, orientation$status[1L], ifelse(any(mod$coverage$status == "PASS"), "EVALUABLE", "NOT_EVALUABLE"), "PASS_RESULT_INDEPENDENT", "SINGLE_COHORT_ORTHOGONAL_CASE_PAIRED_PROTEIN_SUPPORT_ONLY"),
    boundary = c("Only project-local PDC000614 standardized TMT18 data are active.", "TMT18, sample map and biospecimen metadata are mandatory; label map is an optional reagent-definition cross-check.", "OLFML2B is exact-token matched; OLFML2A and OLFM2 are never merged.", "At least 90% of eligible biological SAMPLE columns must map to case and Tumor/Normal tissue.", "All-column fraction additionally reports excluded reference/internal channels and is not substituted for the eligible-sample denominator.", "Sample-map LabelReagent values are compared with label-map definitions; unresolved optional linkage is reported, not hidden.", "Tumor and Normal are subtracted only within the same analytical TMT plex, then collapsed to one case-level delta.", "Structural evaluability depends on files, mapping, label orientation and same-plex case count, not effect direction.", "Biological direction is reported separately from structural evaluability.", "Sentinels detect broad label inversion; ACTA2 is context-dependent and is excluded from fixed-direction voting.", "Module coverage uses max(3 genes, 60% of the catalog), capped at the catalog size.", "Null or discordant results remain reportable.", "No causal, survival, diagnostic, treatment-selection or multi-cohort protein-validation claim."),
    stringsAsFactors = FALSE
  )
  ol_p5_write_csv(go, file.path(dirs$tables, "20_part5_pdc614_go_no_go_summary.csv"))

  amendment <- data.frame(amendment_date = "2026-07-22", active_release = OLFML2B_PART5_VERSION, active_protein_scope = "PDC000614_standardized only", amendment_reason = "PDC000614 was rebuilt as the active exact-target case-paired protein cohort; v1.2.0 further freezes within-plex subtraction and separates structural evaluability from biological direction.", transparency_rule = "All thresholds and pairing rules are frozen before rerun; null and discordant outcomes remain reportable.", stringsAsFactors = FALSE)
  ol_p5_write_csv(amendment, file.path(dirs$tables, "21_protein_scope_method_amendment.csv"))

  figure_registry <- if (isTRUE(make_figures)) tryCatch(ol_p5_save_figures(dirs, paired, tests, sentinel_summary, mod$scores, modcor), error = function(e) data.frame(figure = "ALL", stem = NA_character_, status = paste0("WARN: ", conditionMessage(e)), stringsAsFactors = FALSE)) else data.frame()
  ol_p5_write_csv(figure_registry, file.path(dirs$tables, "98_figure_registry.csv"))
  ol_p5_write_csv(ol_p5_output_manifest(dirs), file.path(dirs$tables, "99_output_file_manifest.csv"))

  index <- list(version = OLFML2B_PART5_VERSION, generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), active_cohort = OLFML2B_PART5_COHORT, root = root, pdc614_dir = pdc614_dir, dirs = dirs, params = params, file_manifest = specs, input_schema = schema, label_reagent_crosscheck = label_crosscheck, detection = detection, quant_mapping = mapped, mapping_summary = mapping_summary, target_reference_relative_values = values[values$gene == OLFML2B_PART5_TARGET, , drop = FALSE], target_case_plex_tissue = target_case, target_plex_pairs = target_plex_pairs, target_paired = target_paired, target_attrition = attr$attrition, target_attrition_summary = attr$attrition_summary, target_tissue_completeness = attr$tissue, direction_tests = tests, sentinel_summary = sentinel_summary, label_orientation_qc = orientation, module_coverage = mod$coverage, module_scores = mod$scores, module_correlations = modcor, go_no_go = go, analysis_evaluability_status = structural_status, biological_evidence_status = biological_status, overall_status = structural_status, final_gene_lock = FALSE, claim_boundary = "PDC000614 provides single-cohort, within-plex, case-paired orthogonal protein support only.")
  ol_p5_save_rds(index, file.path(dirs$objects, "Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds"))
  ol_p5_save_rds(index, file.path(dirs$objects, "Part5_PDC_production_index.rds"))
  ol_p5_log("INFO", "Part5 completed | evaluability=", structural_status, " | biological_evidence=", biological_status, " | same_plex_paired_cases=", n_pairs, log_file = log_file)
  invisible(index)
}

# ==============================================================================
# Part5 v1.3.0: analytical-plex completeness and robustness audit
# ==============================================================================
OLFML2B_PART5_VERSION <- "v1.3.0_20260722_PDC614_PLEX_COMPLETENESS_AND_LOPO_FIX"

ol_p5_v130_plex_completeness <- function(mapped, target_values) {
  m <- mapped[mapped$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & mapped$sample_id_class == "SAMPLE" & mapped$usable_for_tumor_normal %in% TRUE, , drop = FALSE]
  v <- target_values[target_values$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE, c("quant_column", "log_ratio"), drop = FALSE]
  d <- merge(m, v, by = "quant_column", all.x = TRUE, sort = FALSE)
  plexes <- sort(unique(as.character(d$analytical_sample)))
  ol_p5_bind_rows(lapply(plexes, function(px) {
    z <- d[d$analytical_sample == px, , drop = FALSE]
    nt <- sum(z$tissue_class == "Tumor"); nn <- sum(z$tissue_class == "Normal")
    ft <- sum(z$tissue_class == "Tumor" & is.finite(z$log_ratio)); fn <- sum(z$tissue_class == "Normal" & is.finite(z$log_ratio))
    finite_all <- sum(is.finite(z$log_ratio)); n_all <- nrow(z)
    data.frame(
      analytical_sample = px, n_mapped_primary_columns = n_all,
      n_tumor = nt, n_normal = nn, finite_tumor = ft, finite_normal = fn,
      finite_all = finite_all, measurement_completeness = if (n_all) finite_all / n_all else NA_real_,
      target_quantified_in_plex = finite_all > 0L,
      plex_measurement_status = if (finite_all == 0L) "TARGET_NOT_QUANTIFIED_IN_ENTIRE_PLEX" else if (finite_all == n_all) "PASS_COMPLETE" else "PARTIAL_TARGET_QUANTIFICATION",
      interpretation = if (finite_all == 0L) "Entire analytical plex lacks finite OLFML2B quantification; do not interpret as tissue-specific biological absence." else "Finite values are analyzed without zero imputation.",
      stringsAsFactors = FALSE
    )
  }))
}

ol_p5_v130_leave_one_plex_out <- function(plex_pairs) {
  d <- plex_pairs[plex_pairs$gene == OLFML2B_PART5_TARGET & plex_pairs$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & is.finite(plex_pairs$paired_delta_tumor_minus_normal), , drop = FALSE]
  plexes <- sort(unique(as.character(d$analytical_sample)))
  ol_p5_bind_rows(lapply(plexes, function(px) {
    z <- d[d$analytical_sample != px, , drop = FALSE]
    x <- z$paired_delta_tumor_minus_normal
    wp <- if (length(x) >= OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE) tryCatch(stats::wilcox.test(x, mu=0, exact=FALSE)$p.value,error=function(e)NA_real_) else NA_real_
    nz <- x[x != 0]
    sp <- if (length(nz) >= OLFML2B_PART5_MIN_PAIRS_SUPPORTIVE) stats::binom.test(sum(nz>0),length(nz),p=.5)$p.value else NA_real_
    data.frame(omitted_plex=px,n_remaining_cases=length(x),median_delta=if(length(x))stats::median(x)else NA_real_,mean_delta=if(length(x))mean(x)else NA_real_,positive_fraction=if(length(x))mean(x>0)else NA_real_,wilcox_signed_rank_p=wp,exact_sign_p=sp,direction=if(!length(x))"NOT_EVALUABLE"else if(stats::median(x)>0)"POSITIVE"else if(stats::median(x)<0)"NEGATIVE"else"NULL",stringsAsFactors=FALSE)
  }))
}

ol_p5_v130_plex_direction <- function(plex_pairs) {
  d <- plex_pairs[plex_pairs$gene == OLFML2B_PART5_TARGET & plex_pairs$value_type == OLFML2B_PART5_PRIMARY_VALUE_TYPE & is.finite(plex_pairs$paired_delta_tumor_minus_normal), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  rows <- lapply(split(d, d$analytical_sample), function(z) data.frame(
    analytical_sample=z$analytical_sample[1L],n_cases=nrow(z),median_delta=stats::median(z$paired_delta_tumor_minus_normal),mean_delta=mean(z$paired_delta_tumor_minus_normal),positive_fraction=mean(z$paired_delta_tumor_minus_normal>0),plex_direction=if(stats::median(z$paired_delta_tumor_minus_normal)>0)"POSITIVE"else if(stats::median(z$paired_delta_tumor_minus_normal)<0)"NEGATIVE"else"NULL",stringsAsFactors=FALSE
  ))
  out <- ol_p5_bind_rows(rows)
  nz <- out$median_delta[out$median_delta != 0 & is.finite(out$median_delta)]
  out$plex_level_positive_count <- sum(nz > 0)
  out$plex_level_nonzero_count <- length(nz)
  out$plex_level_exact_sign_p <- if (length(nz)) stats::binom.test(sum(nz>0),length(nz),p=.5)$p.value else NA_real_
  out
}

.ol_p5_v130_core <- run_olfml2b_part5_pdc_production
run_olfml2b_part5_pdc_production <- function(...) {
  index <- .ol_p5_v130_core(...)
  index$version <- OLFML2B_PART5_VERSION
  dirs <- index$dirs
  plex_complete <- ol_p5_v130_plex_completeness(index$quant_mapping, index$target_reference_relative_values)
  lopo <- ol_p5_v130_leave_one_plex_out(index$target_plex_pairs)
  plex_dir <- ol_p5_v130_plex_direction(index$target_plex_pairs)
  whole_missing <- sum(plex_complete$plex_measurement_status == "TARGET_NOT_QUANTIFIED_IN_ENTIRE_PLEX")
  lopo_sensitive <- nrow(lopo) > 0L && any(is.finite(lopo$wilcox_signed_rank_p) & lopo$wilcox_signed_rank_p >= 0.05)
  plex_sign_p <- if (nrow(plex_dir)) unique(plex_dir$plex_level_exact_sign_p)[1L] else NA_real_
  robustness <- if (!nrow(lopo)) "NOT_EVALUABLE" else if (lopo_sensitive || (is.finite(plex_sign_p) && plex_sign_p >= 0.05)) "DIRECTION_POSITIVE_BUT_PLEX_SENSITIVE" else "DIRECTION_ROBUST_TO_LEAVE_ONE_PLEX_OUT"
  tests <- index$direction_tests %||% data.frame()
  if (nrow(tests)) {
    tests$n_analytical_plexes_with_pairs <- if (nrow(plex_dir)) nrow(plex_dir) else 0L
    tests$plex_level_exact_sign_p <- plex_sign_p
    lopo_p <- lopo$wilcox_signed_rank_p[is.finite(lopo$wilcox_signed_rank_p)]
    tests$leave_one_plex_out_wilcox_min <- if (length(lopo_p)) min(lopo_p) else NA_real_
    tests$leave_one_plex_out_wilcox_max <- if (length(lopo_p)) max(lopo_p) else NA_real_
    tests$n_entire_plex_not_quantified <- whole_missing
    tests$n_analytical_plexes_quantified <- sum(plex_complete$target_quantified_in_plex %in% TRUE)
    tests$plex_robustness_status <- robustness
    ol_p5_write_csv(tests,file.path(dirs$tables,"10_OLFML2B_case_paired_direction_tests.csv"))
    index$direction_tests <- tests
  }
  go <- index$go_no_go %||% data.frame()
  if (nrow(go)) go <- go[!go$criterion %in% c("plex_measurement_completeness","leave_one_plex_out_robustness","plex_level_direction"),,drop=FALSE]
  add <- data.frame(criterion=c("plex_measurement_completeness","plex_level_direction","leave_one_plex_out_robustness"),status=c(ifelse(whole_missing>0L,paste0("REVIEW_",whole_missing,"_ENTIRE_PLEX_NOT_QUANTIFIED"),"PASS"),ifelse(is.finite(plex_sign_p),ifelse(plex_sign_p<0.05,"DIRECTIONALLY_POSITIVE_PLEX_LEVEL","POSITIVE_UNCERTAIN_PLEX_LEVEL"),"NOT_EVALUABLE"),robustness),boundary=c("Measurement completeness is reported by analytical TMT plex; entire-plex missingness is not tissue-specific absence.","Each analytical plex contributes one median case-paired direction before plex-level sign testing.","Every paired analytical plex is omitted once; sensitivity cannot upgrade the single-cohort evidence claim."),stringsAsFactors=FALSE)
  index$go_no_go <- ol_p5_bind_rows(list(go,add))
  index$plex_measurement_completeness <- plex_complete
  index$plex_level_direction <- plex_dir
  index$leave_one_plex_out <- lopo
  index$plex_robustness_status <- robustness
  index$biological_evidence_status <- if (grepl("POSITIVE",index$biological_evidence_status %||% "") && grepl("SENSITIVE",robustness)) "DIRECTIONALLY_POSITIVE_BUT_UNCERTAIN_AND_PLEX_SENSITIVE" else index$biological_evidence_status
  index$final_gene_lock <- FALSE
  ol_p5_write_csv(plex_complete,file.path(dirs$tables,"09d_OLFML2B_plex_measurement_completeness.csv"))
  ol_p5_write_csv(lopo,file.path(dirs$tables,"10a_OLFML2B_leave_one_plex_out_sensitivity.csv"))
  ol_p5_write_csv(plex_dir,file.path(dirs$tables,"10b_OLFML2B_plex_level_direction_summary.csv"))
  ol_p5_write_csv(index$go_no_go,file.path(dirs$tables,"20_part5_pdc614_go_no_go_summary.csv"))
  ol_p5_save_rds(index,file.path(dirs$objects,"Part5_PDC000614_OLFML2B_case_paired_protein_validation_index.rds"))
  ol_p5_save_rds(index,file.path(dirs$objects,"Part5_PDC_production_index.rds"))
  invisible(index)
}
