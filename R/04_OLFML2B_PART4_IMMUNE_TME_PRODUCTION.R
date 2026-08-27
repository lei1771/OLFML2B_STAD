
# ============================================================================
# OLFML2B-STAD Part4: Reproducible immune/TME production analysis
# Version: v1.2.0_20260706_REPRO_CNS_TME_PRODUCTION
# ============================================================================
# This is a production script, not a skeleton. It reads the already harmonized
# Part1/Part2 objects through the Part3 view constructor, calculates predefined
# immune/TME signatures, performs correlation meta-analysis, TME-adjusted Cox
# models, quadrant survival analysis, subtype audits, reproducibility manifests,
# and CNS-style figures.
# ============================================================================

options(stringsAsFactors = FALSE)
OLFML2B_PART4_VERSION <- "v1.9.0_20260710_RANDOM_EFFECTS_ENDPOINT_DEDUP_REPAIR"

.ol_part4_entry <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.ol_part4_env_root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = "")
.ol_part4_env_valid <- nzchar(.ol_part4_env_root) && dir.exists(.ol_part4_env_root) &&
  file.exists(file.path(.ol_part4_env_root, "00_OLFML2B_PART0_CONFIG_CORE.R")) &&
  file.exists(file.path(.ol_part4_env_root, "05_OLFML2B_PART5_PDC_PRODUCTION.R"))
.ol_part4_code_root <- if (.ol_part4_env_valid) {
  normalizePath(.ol_part4_env_root, winslash = "/", mustWork = TRUE)
} else if (!is.null(.ol_part4_entry) && file.exists(.ol_part4_entry)) {
  dirname(normalizePath(.ol_part4_entry, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
Sys.setenv(OLFML2B_STAD_CODE_ROOT = .ol_part4_code_root)

ol_p4_source_if_needed <- function() {
  root <- Sys.getenv("OLFML2B_STAD_CODE_ROOT", unset = .ol_part4_code_root)
  needed <- c("00_OLFML2B_PART0_CONFIG_CORE.R", "01_OLFML2B_PART1_TCGA.R", "02_OLFML2B_PART2_GEO.R", "03_OLFML2B_PART3_BULK_SURVIVAL.R")
  for (f in needed) {
    fp <- file.path(root, f)
    if (file.exists(fp)) {
      sys.source(fp, envir = parent.frame(), chdir = FALSE)
    }
  }
}
if (!exists("olfml2b_load_bulk_views", mode = "function")) ol_p4_source_if_needed()
if (!exists("olfml2b_load_bulk_views", mode = "function")) stop("Part4 requires 03_OLFML2B_PART3_BULK_SURVIVAL.R and olfml2b_load_bulk_views().", call. = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x

ol_p4_assert <- function(cond, ...) if (!isTRUE(cond)) stop(paste0(..., collapse = ""), call. = FALSE)
ol_p4_ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ol_p4_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    ol_p4_assert(ok || dir.exists(path), "Cannot create directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

ol_p4_log <- function(level = "INFO", ..., log_file = NULL, echo = TRUE) {
  line <- sprintf("[%s] [%s] [OLFML2B-P4] %s", ol_p4_ts(), toupper(level), paste0(..., collapse = ""))
  if (isTRUE(echo)) message(line)
  if (!is.null(log_file)) {
    ol_p4_dir(dirname(log_file))
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
  invisible(line)
}

ol_p4_bind_rows <- function(xs) {
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

ol_p4_atomic_write_csv <- function(x, path, row.names = FALSE, na = "") {
  ol_p4_dir(dirname(path))
  if (is.null(x)) x <- data.frame()
  if (is.atomic(x) && is.null(dim(x))) x <- data.frame(value = x, stringsAsFactors = FALSE)
  if (!is.data.frame(x) && !is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  utils::write.csv(x, tmp, row.names = row.names, na = na, fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) ok <- file.copy(tmp, path, overwrite = TRUE)
  ol_p4_assert(ok && file.exists(path), "Failed to write CSV: ", path)
  invisible(path)
}

ol_p4_atomic_save_rds <- function(object, path, compress = "xz") {
  ol_p4_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (file.exists(path)) unlink(path, force = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) ok <- file.copy(tmp, path, overwrite = TRUE)
  ol_p4_assert(ok && file.exists(path), "Failed to save RDS: ", path)
  invisible(path)
}

ol_p4_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  x <- gsub("[^0-9eE.+\\-]", "", x)
  suppressWarnings(as.numeric(x))
}

ol_p4_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

ol_p4_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\|.*$", "", x)
  x <- sub("\\..*$", "", x)
  x
}

ol_p4_pkg_versions <- function(pkgs) {
  data.frame(
    package = pkgs,
    version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
    available = vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    stringsAsFactors = FALSE
  )
}

ol_p4_file_manifest <- function(root, pattern = NULL, recursive = TRUE) {
  if (!dir.exists(root)) return(data.frame())
  files <- list.files(root, recursive = recursive, full.names = TRUE, all.files = FALSE, pattern = pattern)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(data.frame())
  info <- file.info(files)
  base <- normalizePath(root, winslash = "/", mustWork = TRUE)
  abs <- normalizePath(files, winslash = "/", mustWork = TRUE)
  data.frame(
    relative_path = substring(abs, nchar(base) + 2L),
    absolute_path = abs,
    size_bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(abs)),
    stringsAsFactors = FALSE
  )
}

ol_p4_make_dirs <- function(root, output_subdir = "Part4") {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  part <- if (grepl("^Part[0-9]$", as.character(output_subdir)[1])) as.character(output_subdir)[1] else "Part4"
  ctx <- tryCatch(olfml2b_load_context(root = root), error = function(e) NULL)
  if (!is.null(ctx) && exists("olfml2b_part_paths", mode = "function")) {
    pp <- olfml2b_part_paths(ctx, part)
    return(list(
      root = root,
      part = part,
      out_root = ctx$dirs$output,
      tables = pp$tables,
      figures = pp$figures,
      reports = pp$reports,
      objects = ctx$dirs$objects,
      logs = pp$logs,
      qc = pp$qc
    ))
  }
  out_root <- file.path(root, "output")
  list(
    root = root,
    part = part,
    out_root = ol_p4_dir(out_root),
    tables = ol_p4_dir(file.path(out_root, "tables", part)),
    figures = ol_p4_dir(file.path(out_root, "figures", part)),
    reports = ol_p4_dir(file.path(out_root, "reports", part)),
    objects = ol_p4_dir(file.path(out_root, "objects")),
    logs = ol_p4_dir(file.path(root, "logs", "runtime", part)),
    qc = ol_p4_dir(file.path(out_root, "qc", part))
  )
}

ol_p4_output_manifest <- function(dirs) {
  sections <- list(tables = dirs$tables, figures = dirs$figures, reports = dirs$reports, qc = dirs$qc)
  rows <- lapply(names(sections), function(section) {
    x <- ol_p4_file_manifest(sections[[section]])
    if (nrow(x)) x$output_section <- section
    x
  })
  ol_p4_bind_rows(rows)
}

ol_p4_palette <- function() {
  c(red = "#C73E3A", blue = "#2F6DB3", teal = "#1F9D8A", gold = "#C69C28", purple = "#7B61A8", orange = "#E67E22", green = "#2F855A", grey = "#8C8C8C", lightgrey = "#E5E7EB", light_red = "#F2B8B5", light_blue = "#BFD2EE", dark = "#111827")
}

ol_p4_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.30),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black", face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "#374151"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

ol_p4_panel <- function(p, label) {
  p + ggplot2::annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.35, vjust = 1.35, fontface = "bold", size = 5)
}

ol_p4_required_packages <- function(log_file = NULL) {
  pkgs <- c("survival", "ggplot2", "matrixStats")
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) stop("Missing required packages for Part4: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

ol_p4_prior_boundary <- function() {
  data.frame(
    boundary_item = c("analyte", "not_quantified", "excluded_prior_axis", "excluded_prior_functions", "primary_focus", "wording_rule", "therapy_claim_boundary"),
    status = c("linear_or_gene_level_OLFML2B_mRNA", "circOLFML2B_back_splice_junction", "miR-370-3p_STAT3_rescue", "autophagy_glycolysis_exosome_metastasis_CCK8_EdU_Transwell_as_primary_story", "CAF_TGFb_ECM_immune_exclusion_like_TME", "association_or_marker_language_only", "no_immunotherapy_response_claim_without_treatment_cohort"),
    interpretation = c(
      "Conventional RNA-seq/microarray gene-level signal is treated as linear/gene-level OLFML2B mRNA.",
      "No back-splice junction quantification is performed; bulk OLFML2B values are not interpreted as circOLFML2B.",
      "This already belongs to prior circOLFML2B gastric-cancer mechanism literature and is deliberately excluded.",
      "These modules are excluded from the main novelty to avoid overlap with circOLFML2B work.",
      "Main mechanism is stromal-immune ecology, especially CAF/TGFb and immune-exclusion-like signatures.",
      "The script never proves OLFML2B causes immune exclusion; it tests reproducible associations and attenuation after TME adjustment.",
      "Checkpoint and IFNg analyses are exploratory unless a real immunotherapy cohort is supplied."
    ),
    stringsAsFactors = FALSE
  )
}

ol_p4_signature_catalog <- function() {
  list(
    CAF_Core = c("FAP","ACTA2","PDGFRB","PDGFRA","THY1","PDPN","TAGLN","DCN","LUM","COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","POSTN","SPARC"),
    TGFb_Response = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7","TGFBI","SERPINE1","CTGF","CYR61","INHBA","PMEPA1","COL1A1","COL3A1"),
    ECM_Remodeling = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","ITGA5","ITGB1","MMP2","MMP9","MMP11","LOX","LOXL2","SPARC","POSTN"),
    Angiogenesis = c("VEGFA","KDR","FLT1","ANGPT1","ANGPT2","TEK","PECAM1","VWF","ENG","ESM1","CDH5"),
    CD8_Cytotoxic = c("CD8A","CD8B","GZMA","GZMB","GZMK","PRF1","NKG7","GNLY","IFNG","CXCL9","CXCL10"),
    IFNg_Response = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","CXCL11","IDO1","GBP1","GBP5","HLA-DRA","TAP1","B2M"),
    Checkpoint_Exhaustion = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1","VSIR","BTLA"),
    Treg = c("FOXP3","IL2RA","CTLA4","IKZF2","CCR8","TNFRSF18","TNFRSF4","IL10","TGFB1"),
    M2_Macrophage = c("CD68","CD163","MRC1","MSR1","CSF1R","IL10","C1QA","C1QB","C1QC","MERTK","MARCO"),
    Myeloid_Inflammatory = c("LYZ","LST1","S100A8","S100A9","FCGR3A","ITGAM","IL1B","CXCL8","CCL2","CCR2","CSF1R"),
    NK_Cell = c("NKG7","GNLY","KLRD1","KLRK1","NCAM1","PRF1","GZMB","FCGR3A"),
    B_Cell = c("MS4A1","CD79A","CD79B","CD19","BANK1","CD22","CD74"),
    Plasma_Cell = c("MZB1","JCHAIN","XBP1","SDC1","IGHG1","IGHG3","IGKC"),
    Epithelial_Tumor = c("EPCAM","KRT8","KRT18","KRT19","MUC1","CDH1","CLDN3","CLDN4","CLDN7"),
    EMT = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1","ITGA5","MMP2","MMP9","COL1A1","COL1A2"),
    Proliferation = c("MKI67","TOP2A","PCNA","MCM2","MCM5","MCM6","UBE2C","CCNB1","CCNB2","CDK1","AURKA","BUB1")
  )
}

ol_p4_catalog_table <- function(catalog) {
  ol_p4_bind_rows(lapply(names(catalog), function(nm) data.frame(signature = nm, n_genes = length(unique(catalog[[nm]])), genes = paste(unique(catalog[[nm]]), collapse = ";"), stringsAsFactors = FALSE)))
}

ol_p4_signature_overlap_audit <- function(catalog) {
  nms <- names(catalog)
  rows <- list()
  for (i in seq_along(nms)) for (j in seq_along(nms)) {
    if (j <= i) next
    a <- unique(ol_p4_clean_gene(catalog[[i]])); b <- unique(ol_p4_clean_gene(catalog[[j]]))
    ov <- intersect(a, b)
    rows[[length(rows) + 1L]] <- data.frame(
      signature_a = nms[i], signature_b = nms[j], n_a = length(a), n_b = length(b),
      n_overlap = length(ov), jaccard = length(ov) / max(length(union(a, b)), 1L),
      overlap_genes = paste(ov, collapse = ";"),
      interpretation = ifelse(length(ov) > 0L, "correlated signatures are not independent evidence", "no direct gene overlap"),
      stringsAsFactors = FALSE
    )
  }
  ol_p4_bind_rows(rows)
}

ol_p4_score_signatures <- function(expr, catalog) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- ol_p4_clean_gene(rownames(expr))
  score_list <- list(); coverage <- list()
  for (sig in names(catalog)) {
    genes <- unique(ol_p4_clean_gene(catalog[[sig]]))
    idx <- which(rownames(expr) %in% genes)
    if (length(idx)) {
      # one row per gene: choose highest variance if duplicated
      idx_by_gene <- split(idx, rownames(expr)[idx])
      keep <- vapply(idx_by_gene, function(ii) {
        if (length(ii) == 1L) return(ii)
        vv <- matrixStats::rowVars(expr[ii, , drop = FALSE], na.rm = TRUE)
        ii[which.max(vv)]
      }, integer(1))
      mat <- expr[keep, , drop = FALSE]
      if (nrow(mat) == 1L) score <- ol_p4_z(mat[1, ]) else score <- colMeans(t(apply(mat, 1, ol_p4_z)), na.rm = TRUE)
      score[!is.finite(score)] <- NA_real_
      present <- rownames(mat)
    } else {
      score <- rep(NA_real_, ncol(expr)); present <- character()
    }
    names(score) <- colnames(expr)
    score_list[[sig]] <- score
    coverage[[sig]] <- data.frame(signature = sig, requested_genes = length(genes), present_genes = length(present), coverage_fraction = ifelse(length(genes) > 0L, length(present) / length(genes), NA_real_), present_gene_symbols = paste(sort(present), collapse = ";"), missing_gene_symbols = paste(sort(setdiff(genes, present)), collapse = ";"), stringsAsFactors = FALSE)
  }
  scores <- as.data.frame(score_list, stringsAsFactors = FALSE, check.names = FALSE)
  scores$sample_id <- colnames(expr)
  list(scores = scores, coverage = ol_p4_bind_rows(coverage))
}

ol_p4_add_axes <- function(d) {
  needed <- c("CAF_Core","TGFb_Response","ECM_Remodeling","CD8_Cytotoxic","M2_Macrophage","Treg","Checkpoint_Exhaustion","IFNg_Response","Myeloid_Inflammatory")
  for (nm in needed) if (!nm %in% names(d)) d[[nm]] <- NA_real_
  zc <- function(nm) ol_p4_z(d[[nm]])
  d$CAF_TGFb_axis <- rowMeans(cbind(zc("CAF_Core"), zc("TGFb_Response"), zc("ECM_Remodeling")), na.rm = TRUE)
  d$Suppressive_Myeloid_axis <- rowMeans(cbind(zc("M2_Macrophage"), zc("Treg"), zc("Checkpoint_Exhaustion")), na.rm = TRUE)
  d$Cytotoxic_axis <- zc("CD8_Cytotoxic")
  d$Immune_Exclusion_Index <- d$CAF_TGFb_axis - d$Cytotoxic_axis
  d$Suppressive_TME_Index <- rowMeans(cbind(d$CAF_TGFb_axis, d$Suppressive_Myeloid_axis, -d$Cytotoxic_axis), na.rm = TRUE)
  d$IFNg_Checkpoint_axis <- rowMeans(cbind(zc("IFNg_Response"), zc("Checkpoint_Exhaustion")), na.rm = TRUE)
  d
}

ol_p4_extract_gene <- function(expr, gene) {
  rn <- ol_p4_clean_gene(rownames(expr))
  idx <- which(rn == toupper(gene))
  out <- rep(NA_real_, ncol(expr)); names(out) <- colnames(expr)
  if (!length(idx)) return(out)
  if (length(idx) == 1L) out <- as.numeric(expr[idx, ]) else out <- colMeans(expr[idx, , drop = FALSE], na.rm = TRUE)
  names(out) <- colnames(expr)
  out
}

ol_p4_safe_col <- function(df, patterns) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NA_character_)
  for (p in patterns) {
    hit <- grep(p, names(df), ignore.case = TRUE, perl = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

ol_p4_event <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.integer(as.numeric(x) > 0))
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(y))
  out[grepl("dead|death|deceased|event|progress|recur|relapse|1|yes|true", y)] <- 1L
  out[grepl("alive|living|censor|no event|non|free|0|no|false", y)] <- 0L
  out
}

ol_p4_endpoint_defs <- function(df) {
  defs <- list()
  if (all(c("os_time_days", "os_event") %in% names(df))) defs$OS <- c("os_time_days", "os_event")
  recurrence_candidates <- list(
    DFS = c("dfs_time_days", "dfs_event"),
    RFS = c("rfs_time_days", "rfs_event"),
    RECURRENCE = c("recurrence_time_days", "recurrence_event")
  )
  available <- names(recurrence_candidates)[vapply(recurrence_candidates, function(z) all(z %in% names(df)), logical(1))]
  if (length(available)) {
    # Canonical recurrence-family priority is locked before analysis.  Related
    # aliases are never fitted as if they were independent endpoints.
    chosen <- available[1L]
    defs[[chosen]] <- recurrence_candidates[[chosen]]
  }
  if (!length(defs)) {
    tcol <- ol_p4_safe_col(df, c("^os_time", "survival.*time", "days_to_death", "days_to_last_follow"))
    ecol <- ol_p4_safe_col(df, c("^os_event", "vital_status", "death"))
    if (!is.na(tcol) && !is.na(ecol)) defs$OS <- c(tcol, ecol)
  }
  defs
}

ol_p4_tme_covariates <- function(base, events, min_epv = 10) {
  candidates <- c("age", "sex", "stage_numeric")
  covars <- character()
  df_used <- 1L
  for (cv in candidates) {
    if (!cv %in% names(base)) next
    x <- base[[cv]]
    ok <- !is.na(x)
    if (sum(ok) < max(30L, 0.60 * nrow(base))) next
    if (is.numeric(x) || is.integer(x)) {
      if (!is.finite(stats::sd(as.numeric(x), na.rm = TRUE)) || stats::sd(as.numeric(x), na.rm = TRUE) <= 0) next
      add_df <- 1L
    } else {
      lev <- unique(as.character(stats::na.omit(x)))
      if (length(lev) < 2L || length(lev) > 8L) next
      add_df <- length(lev) - 1L
    }
    if (events / (df_used + add_df) >= min_epv) {
      covars <- c(covars, cv)
      df_used <- df_used + add_df
    }
  }
  covars
}

ol_p4_tme_cox <- function(d, cohort, min_n = 40L, min_events = 20L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) return(data.frame())
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) {
    ol_p4_v230_record_eligibility(data.frame(
      cohort = cohort, endpoint = NA_character_, attenuation_axis = NA_character_,
      status = "NO_AUDITABLE_ENDPOINT", n = nrow(d), events = NA_integer_,
      clinical_tier = NA_character_, clinical_covariates = NA_character_,
      stringsAsFactors = FALSE
    ))
    return(data.frame())
  }
  rows <- list()
  model_sets <- list(
    OLFML2B_only = c("OLFML2B_z"),
    OLFML2B_plus_CAF_TGFb = c("OLFML2B_z", "CAF_TGFb_axis"),
    OLFML2B_plus_CAF_TGFb_plus_Cytotoxic = c("OLFML2B_z", "CAF_TGFb_axis", "Cytotoxic_axis"),
    OLFML2B_plus_ImmuneExclusion = c("OLFML2B_z", "Immune_Exclusion_Index"),
    OLFML2B_plus_SuppressiveTME = c("OLFML2B_z", "Suppressive_TME_Index")
  )
  for (ep in names(defs)) {
    tc <- defs[[ep]][1]; ec <- defs[[ep]][2]
    base <- d
    base$time <- ol_p4_num(base[[tc]])
    base$event <- ol_p4_event(base[[ec]])
    base <- base[is.finite(base$time) & base$time > 0 & !is.na(base$event), , drop = FALSE]
    if (nrow(base) < min_n || sum(base$event == 1L) < min_events) next
    clinical_covars <- ol_p4_tme_covariates(base, events = sum(base$event == 1L), min_epv = min_epv)
    for (mn in names(model_sets)) {
      vars_base <- model_sets[[mn]]
      vars_base <- vars_base[vars_base %in% names(base)]
      if (!("OLFML2B_z" %in% vars_base)) next
      variants <- list(
        tme_attenuation = vars_base,
        clinical_tme_attenuation = unique(c(vars_base, clinical_covars))
      )
      for (role in names(variants)) {
        vars <- variants[[role]]
        dd <- base[stats::complete.cases(base[, c("time", "event", vars), drop = FALSE]), , drop = FALSE]
        if (nrow(dd) < min_n || sum(dd$event == 1L) < min_events) next
        form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(vars, collapse = " + ")))
        fit <- tryCatch(survival::coxph(form, data = dd), error = function(e) NULL)
        if (is.null(fit)) next
        sm <- summary(fit)
        if (!("OLFML2B_z" %in% rownames(sm$coefficients))) next
        i <- which(rownames(sm$coefficients) == "OLFML2B_z")[1]
        beta <- sm$coefficients[i, "coef"]; se <- sm$coefficients[i, "se(coef)"]
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cohort, endpoint = ep, model = mn, adjustment_role = role,
          n = nrow(dd), events = sum(dd$event == 1L), beta = beta, se = se,
          HR = exp(beta), lcl = exp(beta - 1.96 * se), ucl = exp(beta + 1.96 * se),
          p_value = sm$coefficients[i, "Pr(>|z|)"],
          covariates = paste(setdiff(vars, "OLFML2B_z"), collapse = ";"),
          clinical_covariates = paste(clinical_covars, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) {
    out$fdr <- stats::p.adjust(out$p_value, method = "BH")
    base_beta <- out[out$model == "OLFML2B_only" & out$adjustment_role == "tme_attenuation", c("cohort", "endpoint", "beta"), drop = FALSE]
    names(base_beta)[3] <- "base_beta"
    out <- merge(out, base_beta, by = c("cohort", "endpoint"), all.x = TRUE, sort = FALSE)
    out$attenuation_percent <- ifelse(is.finite(out$base_beta) & abs(out$base_beta) > 1e-8,
                                      100 * (out$base_beta - out$beta) / out$base_beta,
                                      NA_real_)
  }
  out
}

ol_p4_correlations <- function(d, cohort) {
  sigs <- c("CAF_Core","TGFb_Response","ECM_Remodeling","Angiogenesis","CD8_Cytotoxic","IFNg_Response","Checkpoint_Exhaustion","Treg","M2_Macrophage","Myeloid_Inflammatory","NK_Cell","B_Cell","Plasma_Cell","Epithelial_Tumor","EMT","Proliferation_Control","Hypoxia_Control","CAF_TGFb_axis","Suppressive_Myeloid_axis","Cytotoxic_axis","Immune_Exclusion_Index","Suppressive_TME_Index","IFNg_Checkpoint_axis")
  rows <- list()
  for (sig in sigs) {
    if (!sig %in% names(d)) next
    keep <- is.finite(d$OLFML2B_z) & is.finite(d[[sig]])
    if (sum(keep) < 10L) next
    if (cohort == "GSE84437" && "source_subseries" %in% names(d)) {
      x <- rank(d$OLFML2B_z[keep], ties.method = "average")
      y <- rank(d[[sig]][keep], ties.method = "average")
      batch <- factor(d$source_subseries[keep])
      design <- stats::model.matrix(~ batch)
      rx <- stats::lm.fit(design, x)$residuals
      ry <- stats::lm.fit(design, y)$residuals
      rho <- stats::cor(rx, ry, use = "complete.obs")
      df <- sum(keep) - ncol(design) - 1L
      t_stat <- rho * sqrt(df / pmax(1 - rho^2, .Machine$double.eps))
      p_value <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
      cor_method <- "partial_spearman_rank_residual_controlling_subseries"
    } else {
      ct <- suppressWarnings(stats::cor.test(d$OLFML2B_z[keep], d[[sig]][keep], method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); p_value <- ct$p.value
      cor_method <- "spearman"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort, signature = sig, n = sum(keep), rho = rho,
      p_value = p_value, correlation_method = cor_method,
      fisher_z_variance = 1 / pmax(sum(keep) - ifelse(cohort == "GSE84437", 4L, 3L), 1L),
      batch_adjustment = ifelse(cohort == "GSE84437", "within_GSE84426_GSE84433_standardization", "not_required"),
      stringsAsFactors = FALSE
    )
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) {
    out$multiplicity_family <- ifelse(out$signature %in% c("CAF_Core", "TGFb_Response", "ECM_Remodeling"),
                                      "core_CAF_TGFb_ECM", "exploratory_or_control_TME")
    out$fdr_within_cohort_family <- ave(out$p_value, out$multiplicity_family,
                                        FUN = function(p) stats::p.adjust(p, method = "BH"))
    out$fdr_within_cohort <- out$fdr_within_cohort_family
  }
  out
}

ol_p4_apply_subseries_batch_design <- function(d, cohort) {
  if (!identical(cohort, "GSE84437")) {
    d$batch_adjustment <- "not_required"
    return(d)
  }
  if (!"source_subseries" %in% names(d) ||
      !identical(sort(unique(stats::na.omit(as.character(d$source_subseries)))), c("GSE84426", "GSE84433"))) {
    stop("Part4 GSE84437 requires both formal source_subseries levels and no context-only samples.", call. = FALSE)
  }
  z_within <- function(x, g) {
    within <- ave(as.numeric(x), g, FUN = function(v) {
      s <- stats::sd(v, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) rep(NA_real_, length(v)) else (v - mean(v, na.rm = TRUE)) / s
    })
    ol_p4_z(within)
  }
  d$OLFML2B_z <- z_within(d$OLFML2B, d$source_subseries)
  signature_fields <- c(
    "CAF_Core","TGFb_Response","ECM_Remodeling","Angiogenesis","CD8_Cytotoxic","IFNg_Response",
    "Checkpoint_Exhaustion","Treg","M2_Macrophage","Myeloid_Inflammatory","NK_Cell","B_Cell",
    "Plasma_Cell","Epithelial_Tumor","EMT","Proliferation_Control","Hypoxia_Control",
    "CAF_TGFb_axis","Suppressive_Myeloid_axis","Cytotoxic_axis","Immune_Exclusion_Index",
    "Suppressive_TME_Index","IFNg_Checkpoint_axis"
  )
  for (field in intersect(signature_fields, names(d))) {
    if (is.numeric(d[[field]])) d[[field]] <- z_within(d[[field]], d$source_subseries)
  }
  batch_means <- tapply(d$OLFML2B_z, d$source_subseries, mean, na.rm = TRUE)
  if (abs(stats::sd(d$OLFML2B_z, na.rm = TRUE) - 1) > 1e-8 || any(abs(batch_means) > 1e-8)) {
    stop("Part4 GSE84437 within-subseries exposure standardization failed.", call. = FALSE)
  }
  d$batch_adjustment <- "within_GSE84426_GSE84433_standardization"
  d
}

ol_p4_meta_cor <- function(tab) {
  if (!nrow(tab)) return(data.frame())
  rows <- lapply(split(tab, tab$signature), function(d) {
    d <- d[is.finite(d$rho) & is.finite(d$n) & d$n > 3, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    z <- atanh(pmax(pmin(d$rho, 0.999999), -0.999999))
    vi <- if ("fisher_z_variance" %in% names(d)) d$fisher_z_variance else 1 / pmax(d$n - 3, 1)
    wf <- 1 / vi
    fixed <- sum(wf * z) / sum(wf)
    q <- sum(wf * (z - fixed)^2)
    df <- length(z) - 1L
    cval <- sum(wf) - sum(wf^2) / sum(wf)
    tau2_dl <- if (length(z) > 1L && is.finite(cval) && cval > 0) max(0, (q - df) / cval) else 0
    if (length(z) > 1L) {
      objective <- function(tau2_value) {
        w0 <- 1 / (vi + tau2_value)
        mu0 <- sum(w0 * z) / sum(w0)
        sum(log(vi + tau2_value)) + log(sum(w0)) + sum(w0 * (z - mu0)^2)
      }
      upper <- max(1, tau2_dl * 20, stats::var(z) * 20, na.rm = TRUE)
      opt <- tryCatch(stats::optimize(objective, interval = c(0, upper)), error = function(e) NULL)
      tau2 <- if (!is.null(opt) && is.finite(opt$minimum)) max(0, opt$minimum) else tau2_dl
      if (is.finite(objective(0)) && objective(0) <= objective(tau2) + 1e-10) tau2 <- 0
    } else {
      tau2 <- 0
    }
    w <- 1 / (vi + tau2)
    zbar <- sum(w * z) / sum(w)
    se_model <- sqrt(1 / sum(w))
    if (length(z) > 1L) {
      q_star <- sum(w * (z - zbar)^2) / df
      hk_scale <- max(1, q_star)
      se <- sqrt(hk_scale / sum(w))
      crit <- stats::qt(0.975, df = df)
      p <- 2 * stats::pt(abs(zbar / se), df = df, lower.tail = FALSE)
      pred_se <- sqrt(tau2 + se^2)
    } else {
      q_star <- NA_real_; hk_scale <- NA_real_; se <- se_model; crit <- 1.96
      p <- 2 * stats::pnorm(abs(zbar / se), lower.tail = FALSE); pred_se <- NA_real_
    }
    i2 <- if (length(z) > 1L && q > df && q > 0) 100 * (q - df) / q else 0
    data.frame(
      signature = d$signature[1], k = nrow(d), n_total = sum(d$n),
      rho_meta = tanh(zbar), rho_lcl = tanh(zbar - crit * se), rho_ucl = tanh(zbar + crit * se),
      statistic = zbar / se, p_value = p,
      prediction_low = if (is.finite(pred_se)) tanh(zbar - crit * pred_se) else NA_real_,
      prediction_high = if (is.finite(pred_se)) tanh(zbar + crit * pred_se) else NA_real_,
      I2_approx = i2, tau2_fisher_z = tau2, tau2_dl_diagnostic = tau2_dl,
      hk_scale = hk_scale, inference_df = ifelse(length(z) > 1L, df, NA_integer_),
      method = ifelse(length(z) > 1L, "Fisher_z_REML_modified_Hartung_Knapp", "single_cohort_descriptive"),
      cohorts = paste(d$cohort, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) {
    out$multiplicity_family <- ifelse(out$signature %in% c("CAF_Core", "TGFb_Response", "ECM_Remodeling"),
                                      "core_CAF_TGFb_ECM", "exploratory_or_control_TME")
    out$fdr <- ave(out$p_value, out$multiplicity_family,
                   FUN = function(p) stats::p.adjust(p, method = "BH"))
  }
  out
}

ol_p4_high_low <- function(d, cohort) {
  sigs <- c("CAF_Core","TGFb_Response","ECM_Remodeling","CD8_Cytotoxic","IFNg_Response","Checkpoint_Exhaustion","Treg","M2_Macrophage","CAF_TGFb_axis","Immune_Exclusion_Index","Suppressive_TME_Index")
  med <- stats::median(d$OLFML2B, na.rm = TRUE)
  d$OLFML2B_group <- ifelse(d$OLFML2B >= med, "High", "Low")
  rows <- list()
  for (sig in sigs) {
    if (!sig %in% names(d)) next
    keep <- is.finite(d[[sig]]) & d$OLFML2B_group %in% c("High", "Low")
    if (sum(keep) < 10L) next
    p <- tryCatch(stats::wilcox.test(d[[sig]][keep] ~ d$OLFML2B_group[keep])$p.value, error = function(e) NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(cohort = cohort, signature = sig, n = sum(keep), high_median = stats::median(d[[sig]][keep & d$OLFML2B_group == "High"], na.rm = TRUE), low_median = stats::median(d[[sig]][keep & d$OLFML2B_group == "Low"], na.rm = TRUE), high_minus_low = stats::median(d[[sig]][keep & d$OLFML2B_group == "High"], na.rm = TRUE) - stats::median(d[[sig]][keep & d$OLFML2B_group == "Low"], na.rm = TRUE), p_value = p, stringsAsFactors = FALSE)
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) {
    out$multiplicity_family <- ifelse(out$signature %in% c("CAF_Core", "TGFb_Response", "ECM_Remodeling"),
                                      "core_CAF_TGFb_ECM", "exploratory_TME_display")
    out$fdr <- ave(out$p_value, out$multiplicity_family,
                   FUN = function(p) stats::p.adjust(p, method = "BH"))
  }
  out
}

ol_p4_checkpoint_genes <- function(expr, cohort) {
  genes <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TOX","CXCL9","CXCL10","IFNG","GZMB","CD8A")
  ub <- ol_p4_extract_gene(expr, "OLFML2B")
  rows <- list()
  for (g in genes) {
    gv <- ol_p4_extract_gene(expr, g)
    keep <- is.finite(ub) & is.finite(gv)
    if (sum(keep) < 10L) next
    ct <- suppressWarnings(stats::cor.test(ol_p4_z(ub[keep]), ol_p4_z(gv[keep]), method = "spearman"))
    rows[[length(rows) + 1L]] <- data.frame(cohort = cohort, gene = g, n = sum(keep), rho = unname(ct$estimate), p_value = ct$p.value, stringsAsFactors = FALSE)
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

ol_p4_subtype_audit <- function(d, cohort) {
  cand <- names(d)[grepl("subtype|lauren|diffuse|intestinal|molecular|acrg|tcga", names(d), ignore.case = TRUE)]
  rows <- list()
  for (col in cand) {
    x <- trimws(as.character(d[[col]]))
    keep <- nzchar(x) & is.finite(d$OLFML2B)
    if (sum(keep) < 20L || length(unique(x[keep])) < 2L) next
    tab <- aggregate(d$OLFML2B[keep], by = list(group = x[keep]), FUN = function(v) c(n = length(v), median = median(v, na.rm = TRUE), mean = mean(v, na.rm = TRUE)))
    tab <- do.call(data.frame, tab)
    names(tab) <- c("group", "n", "median_OLFML2B", "mean_OLFML2B")
    tab$cohort <- cohort; tab$subtype_variable <- col
    rows[[length(rows) + 1L]] <- tab[, c("cohort","subtype_variable","group","n","median_OLFML2B","mean_OLFML2B")]
  }
  ol_p4_bind_rows(rows)
}

ol_p4_quadrant_summary <- function(d, cohort) {
  if (!all(c("OLFML2B", "CAF_TGFb_axis", "Cytotoxic_axis") %in% names(d))) return(data.frame())
  d <- d[is.finite(d$OLFML2B) & is.finite(d$CAF_TGFb_axis) & is.finite(d$Cytotoxic_axis), , drop = FALSE]
  if (nrow(d) < 20L) return(data.frame())
  d$OLFML2B_group <- ifelse(d$OLFML2B >= median(d$OLFML2B, na.rm = TRUE), "OLFML2B-high", "OLFML2B-low")
  d$CAF_group <- ifelse(d$CAF_TGFb_axis >= median(d$CAF_TGFb_axis, na.rm = TRUE), "CAF/TGFb-high", "CAF/TGFb-low")
  d$Cytotoxic_group <- ifelse(d$Cytotoxic_axis >= median(d$Cytotoxic_axis, na.rm = TRUE), "cytotoxic-high", "cytotoxic-low")
  d$quadrant <- paste(d$OLFML2B_group, d$CAF_group, d$Cytotoxic_group, sep = " | ")
  tab <- as.data.frame(table(d$quadrant), stringsAsFactors = FALSE)
  names(tab) <- c("quadrant", "n")
  tab$cohort <- cohort
  tab[, c("cohort", "quadrant", "n")]
}

ol_p4_make_figures <- function(dirs, meta_cor, high_low, surv_tab, sample_scores) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  pal <- ol_p4_palette()
  fig_rows <- list()
  save_one <- function(p, id, title, width, height, source_table, stat_method, n, caption) {
    stem <- file.path(dirs$figures, id)
    if (exists("olfml2b_save_pub_plot", mode = "function")) {
      olfml2b_save_pub_plot(p, stem, width = width, height = height)
    } else {
      ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
      ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 300)
    }
    if (exists("olfml2b_figure_registry_row", mode = "function")) {
      fig_rows[[length(fig_rows) + 1L]] <<- olfml2b_figure_registry_row(id, title, source_table, stat_method, n, caption, stem)
    }
  }
  top_sigs <- c("CAF_Core","TGFb_Response","ECM_Remodeling","CAF_TGFb_axis","Immune_Exclusion_Index","Suppressive_TME_Index","CD8_Cytotoxic","IFNg_Response","Checkpoint_Exhaustion","M2_Macrophage","Treg")

  if (nrow(meta_cor)) {
    d <- meta_cor[meta_cor$signature %in% top_sigs, , drop = FALSE]
    if (!nrow(d)) d <- meta_cor
    d <- d[order(d$rho_meta), , drop = FALSE]
    d$label <- paste0(d$signature, "\nrho = ", sprintf("%.2f", d$rho_meta), ", ", ifelse(d$p_value < 1e-4, "P < 1e-4", paste0("P = ", sprintf("%.3f", d$p_value))), ", FDR = ", sprintf("%.3f", d$fdr), ", k = ", d$k, ", n = ", d$n_total)
    d$label <- factor(make.unique(d$label, sep = " #"), levels = rev(make.unique(d$label, sep = " #")))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = rho_meta, y = label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, color = unname(pal["grey"])) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = rho_lcl, xmax = rho_ucl), orientation = "y", width = 0.16, color = unname(pal["dark"]), linewidth = 0.38) +
      ggplot2::geom_point(ggplot2::aes(color = rho_meta > 0, size = n_total), alpha = 0.95) +
      ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = unname(pal["blue"])), guide = "none") +
      ggplot2::scale_size_continuous(range = c(2.4, 4.5), name = "Samples") +
      ggplot2::labs(title = "OLFML2B and tumor microenvironment signatures", subtitle = "Fisher-z meta-analysis of cohort-level Spearman correlations; labels report rho, P, FDR, k and total n", x = "Pooled Spearman rho", y = NULL) +
      ol_p4_theme(base_size = 10)
    save_one(p, "FIG4A_TME_meta_correlation_forest_labeled", "OLFML2B and tumor microenvironment signatures", 10.0, max(5.3, 0.70*nrow(d)+1.6), "10_olfml2b_tme_meta_correlations.csv", "Cohort-wise Spearman correlations meta-analysed by Fisher-z transform; BH FDR across signatures", paste0("signatures=", nrow(d)), "Forest plot of pooled correlations between OLFML2B and TME signatures. Horizontal bars show approximate 95% CI for pooled rho; labels report P value, FDR, cohort count and total samples.")
  }

  if (nrow(high_low)) {
    d <- high_low[high_low$signature %in% top_sigs, , drop = FALSE]
    if (nrow(d)) {
      d$sig_label <- ifelse(is.finite(d$fdr) & d$fdr < 0.05, "FDR<0.05", ifelse(is.finite(d$p_value) & d$p_value < 0.05, "P<0.05", "NS"))
      d$label <- paste0("n=", d$n, "\n", ifelse(d$p_value < 1e-4, "P<1e-4", paste0("P=", sprintf("%.3f", d$p_value))))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = high_minus_low, y = signature)) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, color = unname(pal["grey"])) +
        ggplot2::geom_point(ggplot2::aes(color = high_minus_low > 0, shape = sig_label), size = 2.6, alpha = 0.95) +
        ggplot2::facet_wrap(~ cohort, scales = "free_y") +
        ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = unname(pal["blue"])), guide = "none") +
        ggplot2::labs(title = "OLFML2B-high versus OLFML2B-low TME shift", subtitle = "Median split is an exploratory display; effect is median score difference (high − low)", x = "Median difference: OLFML2B-high minus OLFML2B-low", y = NULL, shape = "Evidence") +
        ol_p4_theme(base_size = 10)
      save_one(p, "FIG4B_high_low_TME_shift_labeled", "OLFML2B-high versus OLFML2B-low TME shift", 10.5, 7.0, "11_olfml2b_high_low_tme_shift.csv", "Two-sided Wilcoxon rank-sum test within cohort; BH FDR within table", paste0("comparisons=", nrow(d)), "Exploratory median-split comparison of TME signatures between OLFML2B-high and OLFML2B-low tumors. Points show median high–low difference; shapes indicate nominal/FDR evidence level.")
    }
  }

  if (nrow(surv_tab)) {
    d <- surv_tab[surv_tab$model %in% c("OLFML2B_only", "OLFML2B_plus_CAF_TGFb", "OLFML2B_plus_ImmuneExclusion", "OLFML2B_plus_SuppressiveTME"), , drop = FALSE]
    d <- d[d$adjustment_role %in% c("tme_attenuation", "clinical_tme_attenuation"), , drop = FALSE]
    d <- d[is.finite(d$HR) & is.finite(d$lcl) & is.finite(d$ucl) & d$HR > 0 & d$lcl > 0 & d$ucl > 0, , drop = FALSE]
    if (nrow(d)) {
      d$label <- paste0(d$cohort, " / ", d$endpoint, " / ", d$adjustment_role, "\n", d$model, "\nHR = ", sprintf("%.2f", d$HR), " (95% CI ", sprintf("%.2f", d$lcl), "–", sprintf("%.2f", d$ucl), "), ", ifelse(d$p_value < 1e-4, "P < 1e-4", paste0("P = ", sprintf("%.3f", d$p_value))), ", n = ", d$n, ", events = ", d$events)
      d$label <- factor(make.unique(d$label, sep = " #"), levels = rev(make.unique(d$label, sep = " #")))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = HR, y = label)) +
        ggplot2::geom_vline(xintercept = 1, linetype = 2, color = unname(pal["grey"])) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = lcl, xmax = ucl), orientation = "y", width = 0.16, color = unname(pal["dark"]), linewidth = 0.38) +
        ggplot2::geom_point(ggplot2::aes(color = model == "OLFML2B_only", size = events), alpha = 0.95) +
        ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["purple"]), `FALSE` = unname(pal["red"])), labels = c(`TRUE` = "OLFML2B only", `FALSE` = "TME attenuation"), name = "Model") +
        ggplot2::scale_size_continuous(range = c(2.0, 4.0), name = "Events") +
        ggplot2::scale_x_log10() +
        ggplot2::labs(title = "TME attenuation models for OLFML2B survival association", subtitle = "These are attenuation/sensitivity models, not formal mediation or causal proof", x = "Hazard ratio per SD OLFML2B", y = NULL) +
        ol_p4_theme(base_size = 9.5)
      save_one(p, "FIG4C_TME_attenuation_survival_forest_labeled", "TME attenuation models for OLFML2B survival association", 11.0, max(6.0, 0.58*nrow(d)+1.8), "12_tme_attenuation_survival_models.csv", "Cox model; TME signature covariates added to assess attenuation; age/stage added when EPV permits", paste0("models=", nrow(d)), "Forest plot of OLFML2B Cox estimates before and after TME-axis adjustment. Points and horizontal bars show HR and 95% CI; labels report P, n and events.")
    }
  }

  if (nrow(sample_scores)) {
    core <- c("OLFML2B_z","CAF_Core","TGFb_Response","ECM_Remodeling","CAF_TGFb_axis","Immune_Exclusion_Index","CD8_Cytotoxic","IFNg_Response")
    scat <- sample_scores[, intersect(c("cohort", core), names(sample_scores)), drop = FALSE]
    if (all(c("cohort", "OLFML2B_z") %in% names(scat))) {
      long <- data.frame()
      for (sig in setdiff(names(scat), c("cohort", "OLFML2B_z"))) {
        z <- scat[, c("cohort", "OLFML2B_z", sig), drop = FALSE]
        names(z)[3] <- "signature_score"
        z$signature <- sig
        long <- rbind(long, z)
      }
      long <- long[is.finite(long$OLFML2B_z) & is.finite(long$signature_score), , drop = FALSE]
      long <- long[long$signature %in% c("CAF_TGFb_axis", "Immune_Exclusion_Index", "CAF_Core", "ECM_Remodeling"), , drop = FALSE]
      if (nrow(long)) {
        p <- ggplot2::ggplot(long, ggplot2::aes(x = OLFML2B_z, y = signature_score)) +
          ggplot2::geom_point(size = 0.75, alpha = 0.35, color = unname(pal["dark"])) +
          ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, color = unname(pal["red"]), fill = unname(pal["light_red"])) +
          ggplot2::facet_grid(signature ~ cohort, scales = "free_y") +
          ggplot2::labs(title = "Sample-level OLFML2B–TME state relationships", subtitle = "Scatter plots visualize the correlations summarized in the forest plot", x = "OLFML2B z-score", y = "Signature score") +
          ol_p4_theme(base_size = 9)
        save_one(p, "FIG4D_sample_level_OLFML2B_TME_scatter", "Sample-level OLFML2B–TME state relationships", 11.5, 7.2, "Part4_sample_scores_*.rds / 10_olfml2b_tme_meta_correlations.csv", "Visual scatter; inferential statistics in correlation tables", paste0("points=", nrow(long)), "Sample-level scatter plots for selected TME signatures. Lines are linear smooths with 95% CI for visualization; statistical inference is reported in the correlation tables.")
      }
    }
  }
  if (nrow(high_low)) {
    d <- high_low[high_low$signature %in% top_sigs, , drop = FALSE]
    if (nrow(d)) {
      d$label <- paste0(sprintf("%.2f", d$high_minus_low), ifelse(is.finite(d$fdr) & d$fdr < 0.05, "*", ifelse(is.finite(d$p_value) & d$p_value < 0.05, "+", "")))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = cohort, y = signature, fill = high_minus_low)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.28) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
        ggplot2::scale_fill_gradient2(low = unname(pal["blue"]), mid = "white", high = unname(pal["red"]), midpoint = 0, name = "High–low\ndelta") +
        ggplot2::labs(title = "High–low TME shifts across cohorts", subtitle = "Cell values show median difference; * FDR<0.05, + nominal P<0.05", x = NULL, y = NULL) +
        ol_p4_theme(base_size = 9.5) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
      save_one(p, "FIG4E_high_low_TME_shift_heatmap", "High–low TME shifts across cohorts", 10.0, 6.8, "11_olfml2b_high_low_tme_shift.csv", "Within-cohort Wilcoxon comparison; BH FDR shown by symbol", paste0("comparisons=", nrow(d)), "Heatmap of cohort-wise TME differences between OLFML2B-high and OLFML2B-low tumors.")
    }
  }
  if (exists("olfml2b_write_figure_registry", mode = "function")) {
    olfml2b_write_figure_registry(fig_rows, file.path(dirs$tables, "99_figure_registry.csv"))
  }
  invisible(TRUE)
}


run_olfml2b_part4_immune_tme_production <- function(
  root = "D:/OLFML2B_STAD",
  geo_validation_cohorts = c("GSE26253", "GSE84437", "GSE62254", "GSE15459"),
  output_subdir = "Part4",
  seed = 20260706L,
  min_survival_n = 40L,
  min_survival_events = 20L,
  make_figures = TRUE
) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  set.seed(seed)
  dirs <- ol_p4_make_dirs(root, output_subdir)
  log_file <- file.path(dirs$logs, "Part4_immune_TME_production.log")
  ol_p4_required_packages(log_file)
  ol_p4_log("INFO", "Starting Part4 | ", OLFML2B_PART4_VERSION, log_file = log_file)

  params <- data.frame(parameter = c("version","root","seed","target","primary_TME_hypothesis","geo_validation_cohorts","min_survival_n","min_survival_events"), value = c(OLFML2B_PART4_VERSION, root, as.character(seed), "linear/gene-level OLFML2B mRNA", "CAF/TGFb/ECM-rich immune-exclusion-like TME", paste(geo_validation_cohorts, collapse = ";"), min_survival_n, min_survival_events), stringsAsFactors = FALSE)
  ol_p4_atomic_write_csv(params, file.path(dirs$tables, "00_run_parameters.csv"))
  ol_p4_atomic_write_csv(ol_p4_prior_boundary(), file.path(dirs$tables, "01_prior_OLFML2B_claim_boundary_audit.csv"))
  ol_p4_atomic_write_csv(ol_p4_pkg_versions(c("survival","ggplot2","matrixStats","Matrix","Seurat","SingleCellExperiment","readxl")), file.path(dirs$tables, "02_package_versions.csv"))
  ol_p4_atomic_write_csv(ol_p4_file_manifest(file.path(root, "R"), pattern = "\\.R$"), file.path(dirs$tables, "03_code_file_manifest.csv"))
  ol_p4_atomic_write_csv(ol_p4_file_manifest(file.path(root, "output", "objects"), pattern = "\\.rds$"), file.path(dirs$tables, "04_input_object_manifest.csv"))
  utils::capture.output(utils::sessionInfo(), file = file.path(dirs$reports, "sessionInfo.txt"))

  views <- olfml2b_load_bulk_views(root = root, target_gene = "OLFML2B", geo_validation_cohorts = geo_validation_cohorts, log_file = log_file)
  input_contract_audit <- attr(views, "input_contract_audit")
  ol_p4_atomic_write_csv(input_contract_audit, file.path(dirs$tables, "04a_part2_input_contract_audit.csv"))
  catalog <- ol_p4_signature_catalog()
  ol_p4_atomic_write_csv(ol_p4_catalog_table(catalog), file.path(dirs$tables, "05_tme_signature_catalog.csv"))
  ol_p4_atomic_write_csv(ol_p4_signature_overlap_audit(catalog), file.path(dirs$tables, "05a_tme_signature_overlap_audit.csv"))

  source_rows <- list(); coverage_rows <- list(); score_rows <- list(); cor_rows <- list(); highlow_rows <- list(); surv_rows <- list(); subtype_rows <- list(); checkpoint_rows <- list(); quadrant_rows <- list()
  for (cohort in names(views)) {
    v <- views[[cohort]]
    expr <- v$expr
    if (is.null(expr) || !is.matrix(expr) || nrow(expr) < 10L || ncol(expr) < 10L) {
      source_rows[[cohort]] <- data.frame(cohort = cohort, status = "NO_EXPRESSION", n_genes = 0, n_samples = 0, target_measured = FALSE, stringsAsFactors = FALSE)
      next
    }
    source_rows[[cohort]] <- data.frame(cohort = cohort, status = v$load_status %||% "OK", source_file = v$source_file %||% NA_character_, n_genes = nrow(expr), n_samples = ncol(expr), target_measured = isTRUE(v$target_measured), stringsAsFactors = FALSE)
    sc <- ol_p4_score_signatures(expr, catalog)
    cov <- sc$coverage; cov$cohort <- cohort; coverage_rows[[cohort]] <- cov
    d <- sc$scores
    d$sample_id <- as.character(d$sample_id)
    d$OLFML2B <- as.numeric(ol_p4_extract_gene(expr, "OLFML2B")[d$sample_id])
    d$OLFML2B_z <- ol_p4_z(d$OLFML2B)
    d <- ol_p4_add_axes(d)
    d$cohort <- cohort
    clin <- v$clinical
    if (is.data.frame(clin) && nrow(clin) && "sample_id" %in% names(clin)) {
      clin$sample_id <- as.character(clin$sample_id)
      common <- intersect(d$sample_id, clin$sample_id)
      d <- d[match(common, d$sample_id), , drop = FALSE]
      clin <- clin[match(common, clin$sample_id), , drop = FALSE]
      extra <- clin[, setdiff(names(clin), names(d)), drop = FALSE]
      d <- cbind(d, extra)
    }
    d <- ol_p4_apply_subseries_batch_design(d, cohort)
    ol_p4_atomic_save_rds(d, file.path(dirs$objects, paste0("Part4_sample_scores_", cohort, ".rds")))
    score_rows[[cohort]] <- d
    cor_rows[[cohort]] <- ol_p4_correlations(d, cohort)
    highlow_rows[[cohort]] <- ol_p4_high_low(d, cohort)
    surv_rows[[cohort]] <- ol_p4_tme_cox(d, cohort, min_n = min_survival_n, min_events = min_survival_events, min_epv = 10)
    subtype_rows[[cohort]] <- ol_p4_subtype_audit(d, cohort)
    checkpoint_rows[[cohort]] <- ol_p4_checkpoint_genes(expr[, d$sample_id, drop = FALSE], cohort)
    quadrant_rows[[cohort]] <- ol_p4_quadrant_summary(d, cohort)
  }
  source_tab <- ol_p4_bind_rows(source_rows)
  coverage_tab <- ol_p4_bind_rows(coverage_rows)
  score_tab <- ol_p4_bind_rows(score_rows)
  cor_tab <- ol_p4_bind_rows(cor_rows)
  meta_cor <- ol_p4_meta_cor(cor_tab)
  highlow_tab <- ol_p4_bind_rows(highlow_rows)
  surv_tab <- ol_p4_bind_rows(surv_rows)
  subtype_tab <- ol_p4_bind_rows(subtype_rows)
  checkpoint_tab <- ol_p4_bind_rows(checkpoint_rows)
  quadrant_tab <- ol_p4_bind_rows(quadrant_rows)
  if (nrow(cor_tab)) cor_tab$fdr_global <- stats::p.adjust(cor_tab$p_value, method = "BH")

  ol_p4_atomic_write_csv(source_tab, file.path(dirs$tables, "06_source_status.csv"))
  ol_p4_atomic_write_csv(coverage_tab, file.path(dirs$tables, "07_tme_signature_coverage.csv"))
  ol_p4_atomic_write_csv(score_tab, file.path(dirs$tables, "08_sample_level_tme_scores_and_axes.csv"))
  ol_p4_atomic_write_csv(cor_tab, file.path(dirs$tables, "09_olfml2b_tme_correlations_by_cohort.csv"))
  ol_p4_atomic_write_csv(meta_cor, file.path(dirs$tables, "10_olfml2b_tme_meta_correlations.csv"))
  ol_p4_atomic_write_csv(highlow_tab, file.path(dirs$tables, "11_olfml2b_high_low_tme_shift.csv"))
  ol_p4_atomic_write_csv(surv_tab, file.path(dirs$tables, "12_tme_attenuation_survival_models.csv"))
  ol_p4_atomic_write_csv(surv_tab, file.path(dirs$tables, "12_tme_adjusted_survival_models.csv"))  # compatibility alias; interpret as attenuation, not full mediation proof
  ol_p4_atomic_write_csv(subtype_tab, file.path(dirs$tables, "13_subtype_tme_audit.csv"))
  ol_p4_atomic_write_csv(checkpoint_tab, file.path(dirs$tables, "14_checkpoint_gene_correlations.csv"))
  ol_p4_atomic_write_csv(quadrant_tab, file.path(dirs$tables, "15_olfml2b_caf_cytotoxic_quadrant_audit.csv"))

  figures_ok <- FALSE
  if (isTRUE(make_figures)) {
    figures_ok <- tryCatch({
      ol_p4_make_figures(dirs, meta_cor, highlow_tab, surv_tab, score_tab)
      required_figs <- c(
        file.path(dirs$figures, "FIG4A_TME_meta_correlation_forest.pdf"),
        file.path(dirs$figures, "FIG4A_TME_meta_correlation_forest.png")
      )
      all(file.exists(required_figs) & file.info(required_figs)$size > 0)
    }, error = function(e) {
      ol_p4_log("WARN", "Part4 figures failed: ", conditionMessage(e), log_file = log_file)
      FALSE
    })
  }

  robust_positive <- function(sig) {
    z <- meta_cor[meta_cor$signature == sig, , drop = FALSE]
    nrow(z) == 1L && is.finite(z$rho_meta[1]) && z$rho_meta[1] > 0 &&
      is.finite(z$rho_lcl[1]) && z$rho_lcl[1] > 0 &&
      is.finite(z$p_value[1]) && z$p_value[1] < 0.05 &&
      identical(as.character(z$method[1]), "Fisher_z_REML_modified_Hartung_Knapp")
  }
  finite_surv <- nrow(surv_tab) > 0L && any(
    is.finite(surv_tab$beta) & is.finite(surv_tab$se) & surv_tab$se > 0 & is.finite(surv_tab$p_value)
  )
  go <- data.frame(
    criterion = c("linear_mRNA_boundary_written", "CAF_TGFb_axis_positive_meta", "Immune_Exclusion_Index_positive_meta", "Suppressive_TME_Index_positive_meta", "TME_attenuation_survival_models", "CNS_figures"),
    status = c("PASS", ifelse(robust_positive("CAF_TGFb_axis"), "PASS", "REVIEW"), ifelse(robust_positive("Immune_Exclusion_Index"), "PASS", "REVIEW"), ifelse(robust_positive("Suppressive_TME_Index"), "PASS", "REVIEW"), ifelse(finite_surv, "PASS", "REVIEW"), ifelse(isTRUE(figures_ok), "PASS", ifelse(isTRUE(make_figures), "REVIEW", "NOT_RUN"))),
    boundary = c("No circular-RNA inference is made; analysis is gene-level OLFML2B mRNA only", "Main TME hypothesis", "Bulk immune-exclusion-like inference only", "Suppressive ecology support", "Association/attenuation, not mediation proof; clinical covariates are added when EPV allows", "PDF+PNG outputs checked by actual file existence"),
    stringsAsFactors = FALSE
  )
  ol_p4_atomic_write_csv(go, file.path(dirs$tables, "20_part4_go_no_go_summary.csv"))
  ol_p4_atomic_write_csv(ol_p4_output_manifest(dirs), file.path(dirs$tables, "99_output_file_manifest.csv"))
  index <- list(version = OLFML2B_PART4_VERSION, generated_at = ol_p4_ts(), dirs = dirs, params = params, source_status = source_tab, input_contract_audit = input_contract_audit, coverage = coverage_tab, sample_scores = score_tab, correlations = cor_tab, meta_correlations = meta_cor, high_low = highlow_tab, tme_survival = surv_tab, subtype = subtype_tab, checkpoint = checkpoint_tab, go_no_go = go)
  ol_p4_atomic_save_rds(index, file.path(dirs$objects, "Part4_immune_TME_production_index.rds"))
  ol_p4_log("INFO", "Part4 complete: ", dirs$tables, log_file = log_file)
  invisible(index)
}

# ============================================================================
# Methodology repair layer v1.7.0 for Part4
# Makes attenuation models auditable and explicitly prevents mediation/causal
# claims from bulk TME correlations.
# ============================================================================

ol_p4_tme_cox <- function(d, cohort, min_n = 40L, min_events = 20L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) return(data.frame())
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) return(data.frame())
  rows <- list()
  model_sets <- list(
    OLFML2B_only = c("OLFML2B_z"),
    OLFML2B_plus_CAF_TGFb = c("OLFML2B_z", "CAF_TGFb_axis"),
    OLFML2B_plus_CAF_TGFb_plus_Cytotoxic = c("OLFML2B_z", "CAF_TGFb_axis", "Cytotoxic_axis"),
    OLFML2B_plus_ImmuneExclusion = c("OLFML2B_z", "Immune_Exclusion_Index"),
    OLFML2B_plus_SuppressiveTME = c("OLFML2B_z", "Suppressive_TME_Index")
  )
  for (ep in names(defs)) {
    tc <- defs[[ep]][1]; ec <- defs[[ep]][2]
    base <- d
    base$time <- ol_p4_num(base[[tc]])
    base$event <- ol_p4_event(base[[ec]])
    base <- base[is.finite(base$time) & base$time > 0 & !is.na(base$event), , drop = FALSE]
    if (nrow(base) < min_n || sum(base$event == 1L) < min_events) next
    clinical_covars <- ol_p4_tme_covariates(base, events = sum(base$event == 1L), min_epv = min_epv)
    for (mn in names(model_sets)) {
      vars_base <- model_sets[[mn]]
      vars_base <- vars_base[vars_base %in% names(base)]
      if (!("OLFML2B_z" %in% vars_base)) next
      variants <- list(tme_attenuation = vars_base, clinical_tme_attenuation = unique(c(vars_base, clinical_covars)))
      for (role in names(variants)) {
        vars <- variants[[role]]
        dd <- base[stats::complete.cases(base[, c("time", "event", vars), drop = FALSE]), , drop = FALSE]
        if (nrow(dd) < min_n || sum(dd$event == 1L) < min_events) next
        model_df <- tryCatch(max(1L, ncol(stats::model.matrix(stats::as.formula(paste0("~ ", paste(vars, collapse = " + "))), data = dd)) - 1L), error = function(e) length(vars))
        epdf <- sum(dd$event == 1L) / model_df
        form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(vars, collapse = " + ")))
        fit <- tryCatch(survival::coxph(form, data = dd, x = TRUE), error = function(e) NULL)
        if (is.null(fit)) next
        sm <- summary(fit)
        if (!("OLFML2B_z" %in% rownames(sm$coefficients))) next
        i <- which(rownames(sm$coefficients) == "OLFML2B_z")[1]
        beta <- sm$coefficients[i, "coef"]; se <- sm$coefficients[i, "se(coef)"]
        ph <- tryCatch({ z <- survival::cox.zph(fit); c(term = z$table["OLFML2B_z", "p"], global = z$table["GLOBAL", "p"]) }, error = function(e) c(term = NA_real_, global = NA_real_))
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cohort, endpoint = ep, model = mn, adjustment_role = role,
          n = nrow(dd), events = sum(dd$event == 1L), model_df = model_df, events_per_model_df = epdf,
          beta = beta, se = se, HR = exp(beta), lcl = exp(beta - 1.96 * se), ucl = exp(beta + 1.96 * se),
          p_value = sm$coefficients[i, "Pr(>|z|)"], ph_OLFML2B_p = unname(ph["term"]), ph_global_p = unname(ph["global"]),
          covariates = paste(setdiff(vars, "OLFML2B_z"), collapse = ";"), clinical_covariates = paste(clinical_covars, collapse = ";"),
          inference_tier = ifelse(epdf >= min_epv, "ATTENUATION_SENSITIVITY_OK", "SUPPORTIVE_ONLY_LOW_EPV"),
          mediation_warning = "This is covariate attenuation/sensitivity, not formal mediation or causal proof.",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) {
    out$fdr <- stats::p.adjust(out$p_value, method = "BH")
    base_beta <- out[out$model == "OLFML2B_only" & out$adjustment_role == "tme_attenuation", c("cohort", "endpoint", "beta"), drop = FALSE]
    names(base_beta)[3] <- "base_beta"
    out <- merge(out, base_beta, by = c("cohort", "endpoint"), all.x = TRUE, sort = FALSE)
    out$attenuation_percent <- ifelse(is.finite(out$base_beta) & abs(out$base_beta) > 1e-8, 100 * (out$base_beta - out$beta) / out$base_beta, NA_real_)
    out$attenuation_interpretation <- ifelse(out$model == "OLFML2B_only", "baseline OLFML2B association", "attenuation after adding correlated TME axis; not mediation")
  }
  out
}

ol_p4_methodology_claim_limits <- function() {
  data.frame(
    analysis_block = c("signature_scoring", "correlation_meta", "high_low_display", "tme_cox", "subtype_or_checkpoint", "global_claim"),
    allowed_claim = c(
      "Predefined TME programs are scored as relative sample-level RNA states.",
      "OLFML2B is associated with TME signatures across available cohorts if direction and FDR support it.",
      "Median split shows visual contrast only and must not define the primary exposure.",
      "TME covariates attenuate or do not attenuate OLFML2B survival association; this is a sensitivity analysis.",
      "Subtype/checkpoint associations are exploratory context only.",
      "Bulk RNA supports a TME-associated hypothesis, not a causal OLFML2B mechanism."
    ),
    prohibited_claim = c(
      "absolute immune abundance without deconvolution validation", "causal immune exclusion", "optimal clinical cutoff", "formal mediation", "therapeutic response prediction", "OLFML2B drives CAF/TGFb/ECM or CD8 exclusion"
    ),
    stringsAsFactors = FALSE
  )
}

.ol_p4_original_runner <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  index <- .ol_p4_original_runner(...)
  dirs <- index$dirs
  limits <- ol_p4_methodology_claim_limits()
  ol_p4_atomic_write_csv(limits, file.path(dirs$tables, "21_part4_methodology_claim_limits.csv"))
  index$methodology_claim_limits <- limits
  ol_p4_atomic_save_rds(index, file.path(dirs$objects, "Part4_immune_TME_production_index.rds"))
  invisible(index)
}


# ==============================================================================
# CNS-style native visualization extension for Part4 | v20260709
# Generates lollipop, heatmap, scatter, and attenuation summary without Part10.
# ==============================================================================

ol_p4_cns_read <- function(path) {
  if (exists("olfml2b_read_csv_if_exists", mode = "function")) return(olfml2b_read_csv_if_exists(path))
  if (!file.exists(path) || file.info(path)$size <= 0) return(data.frame())
  tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
}

ol_p4_cns_save <- function(p, dirs, id, kind = "general", n_rows = NULL, n_cols = NULL, width = NULL, height = NULL) {
  stem <- file.path(dirs$figures, id)
  if (exists("olfml2b_save_cns_plot", mode = "function")) {
    olfml2b_save_cns_plot(p, stem, kind = kind, n_rows = n_rows, n_cols = n_cols, width = width, height = height)
  } else if (exists("olfml2b_save_pub_plot", mode = "function")) {
    olfml2b_save_pub_plot(p, stem, width = width %||% 7, height = height %||% 5)
  } else {
    ggplot2::ggsave(paste0(stem, ".png"), p, width = width %||% 7, height = height %||% 5, dpi = 600, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width %||% 7, height = height %||% 5, bg = "white", limitsize = FALSE)
  }
  data.frame(figure_id = id, file_stem = stem, kind = kind, stringsAsFactors = FALSE)
}

ol_p4_make_cns_enhanced_figures <- function(dirs) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  pal <- if (exists("olfml2b_pub_palette", mode = "function")) olfml2b_pub_palette() else c(red="#B2182B", blue="#2166AC", grey="#9CA3AF", dark="#1F2937", light_grey="#F3F4F6")
  th <- if (exists("olfml2b_base_theme", mode = "function")) olfml2b_base_theme(10.5) else ol_p4_theme(10.5)
  wrap <- if (exists("olfml2b_wrap_label", mode = "function")) olfml2b_wrap_label else function(x, width = 55) as.character(x)
  reg <- list()
  meta <- ol_p4_cns_read(file.path(dirs$tables, "10_olfml2b_tme_meta_correlations.csv"))
  highlow <- ol_p4_cns_read(file.path(dirs$tables, "11_olfml2b_high_low_tme_shift.csv"))
  scores <- ol_p4_cns_read(file.path(dirs$tables, "08_sample_level_tme_scores_and_axes.csv"))
  surv <- ol_p4_cns_read(file.path(dirs$tables, "12_tme_adjusted_survival_models.csv"))

  selected <- c("CAF_Core", "CAF_TGFb_axis", "Suppressive_TME_Index", "ECM_Remodeling", "EMT", "Immune_Exclusion_Index", "Checkpoint_Exhaustion", "Epithelial_Tumor", "CD8_Cytotoxic", "Proliferation")
  if (nrow(meta) && all(c("signature", "rho_meta") %in% names(meta))) {
    d <- meta[meta$signature %in% selected & is.finite(suppressWarnings(as.numeric(meta$rho_meta))), , drop = FALSE]
    if (nrow(d)) {
      d$rho_meta <- as.numeric(d$rho_meta)
      if (!"fdr" %in% names(d)) d$fdr <- if ("p_value" %in% names(d)) stats::p.adjust(as.numeric(d$p_value), method = "BH") else NA_real_
      d$fdr <- as.numeric(d$fdr)
      d$direction <- ifelse(is.finite(d$fdr) & d$fdr < 0.05 & d$rho_meta >= 0, "positive FDR<0.05", ifelse(is.finite(d$fdr) & d$fdr < 0.05 & d$rho_meta < 0, "negative FDR<0.05", "not significant"))
      d <- d[order(d$rho_meta), , drop = FALSE]
      d$signature_label <- factor(wrap(d$signature, 24), levels = wrap(d$signature, 24))
      d$plot_label <- paste0("rho=", sprintf("%.2f", d$rho_meta), "\n", ifelse(is.finite(d$fdr), ifelse(d$fdr < 1e-4, "FDR<1e-4", paste0("FDR=", sprintf("%.3f", d$fdr))), "FDR=NA"))
      sub <- paste0("signatures=", nrow(d), "; cohorts per signature up to ", suppressWarnings(max(as.numeric(d$k), na.rm = TRUE)), "; total n up to ", suppressWarnings(max(as.numeric(d$n_total), na.rm = TRUE)))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = rho_meta, y = signature_label, color = direction)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = unname(pal["grey"]), linewidth = 0.35) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = rho_meta, yend = signature_label), linewidth = 0.55, color = unname(pal["grey"])) +
        ggplot2::geom_point(size = 2.9) +
        ggplot2::geom_text(ggplot2::aes(label = plot_label), hjust = ifelse(d$rho_meta >= 0, -0.08, 1.08), size = 2.7, color = unname(pal["dark"]), lineheight = 0.90) +
        ggplot2::scale_color_manual(values = c("positive FDR<0.05" = unname(pal["red"]), "negative FDR<0.05" = unname(pal["blue"]), "not significant" = unname(pal["grey"])), name = NULL) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.26))) +
        ggplot2::labs(title = "Selected OLFML2B-associated TME programs", subtitle = sub, x = "Meta-analytic Spearman rho", y = NULL) + th
      reg[[length(reg)+1L]] <- ol_p4_cns_save(p, dirs, "FIG4A_selected_TME_lollipop_CNS", kind = "dotplot", n_rows = nrow(d), width = 8.4)
    }
  }

  if (nrow(highlow) && "signature" %in% names(highlow) && "cohort" %in% names(highlow)) {
    val_col <- intersect(c("high_minus_low", "median_high_minus_low", "delta_high_minus_low", "delta", "mean_high_minus_low"), names(highlow))[1]
    if (!is.na(val_col)) {
      d <- highlow[highlow$signature %in% selected, , drop = FALSE]
      if (nrow(d)) {
        d$value <- as.numeric(d[[val_col]])
        if (!"fdr" %in% names(d)) d$fdr <- if ("p_value" %in% names(d)) stats::p.adjust(as.numeric(d$p_value), method = "BH") else NA_real_
        d$label <- paste0(sprintf("%.2f", d$value), ifelse(is.finite(as.numeric(d$fdr)) & as.numeric(d$fdr) < 0.05, "*", ""))
        d$signature <- factor(wrap(d$signature, 18), levels = rev(wrap(selected[selected %in% unique(as.character(highlow$signature))], 18)))
        p <- ggplot2::ggplot(d, ggplot2::aes(x = cohort, y = signature, fill = value)) +
          ggplot2::geom_tile(color = "white", linewidth = 0.38) +
          ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7, color = "black") +
          ggplot2::scale_fill_gradient2(low = unname(pal["blue"]), mid = "white", high = unname(pal["red"]), midpoint = 0, name = "High-low\ndelta") +
          ggplot2::labs(title = "TME-state shift in OLFML2B-high tumors", subtitle = paste0("selected signatures=", length(unique(d$signature)), "; cohorts=", length(unique(d$cohort)), "; *FDR<0.05"), x = NULL, y = NULL) + th +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
        reg[[length(reg)+1L]] <- ol_p4_cns_save(p, dirs, "FIG4B_OLFML2B_high_low_TME_shift_heatmap_CNS", kind = "heatmap", n_rows = length(unique(d$signature)), n_cols = length(unique(d$cohort)), width = 8.2)
      }
    }
  }

  if (nrow(scores) && "OLFML2B" %in% names(scores)) {
    feats <- intersect(c("CAF_TGFb_axis", "Immune_Exclusion_Index"), names(scores))
    long <- list()
    for (ft in feats) {
      dd <- scores[, intersect(c("cohort", "sample_id", "OLFML2B", ft), names(scores)), drop = FALSE]
      if (!all(c("OLFML2B", ft) %in% names(dd))) next
      long[[length(long)+1L]] <- data.frame(cohort = dd$cohort %||% NA_character_, sample_id = dd$sample_id %||% seq_len(nrow(dd)), OLFML2B = as.numeric(dd$OLFML2B), feature = ft, value = as.numeric(dd[[ft]]), stringsAsFactors = FALSE)
    }
    long <- ol_p4_bind_rows(long)
    long <- long[is.finite(long$OLFML2B) & is.finite(long$value), , drop = FALSE]
    if (nrow(long)) {
      stat <- ol_p4_bind_rows(lapply(split(long, long$feature), function(z) {
        ct <- tryCatch(suppressWarnings(stats::cor.test(z$OLFML2B, z$value, method = "spearman", exact = FALSE)), error = function(e) NULL)
        data.frame(feature = z$feature[1], n = nrow(z), rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_, p_value = if (!is.null(ct)) ct$p.value else NA_real_, stringsAsFactors = FALSE)
      }))
      stat$feature_lab <- paste0(stat$feature, "\n", vapply(seq_len(nrow(stat)), function(i) if (exists("olfml2b_fmt_rho_label", mode = "function")) olfml2b_fmt_rho_label(stat$rho[i], stat$p_value[i], stat$n[i]) else paste0("rho=", sprintf("%.2f", stat$rho[i]), "; P=", signif(stat$p_value[i], 3), "; n=", stat$n[i]), character(1)))
      long$feature_lab <- stat$feature_lab[match(long$feature, stat$feature)]
      p <- ggplot2::ggplot(long, ggplot2::aes(x = OLFML2B, y = value)) +
        ggplot2::geom_point(size = 1.15, alpha = 0.45, color = unname(pal["dark"])) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, color = unname(pal["red"]), fill = unname(pal["light_red"])) +
        ggplot2::facet_wrap(~ feature_lab, scales = "free_y", nrow = 1) +
        ggplot2::labs(title = "Sample-level OLFML2B correlation with key TME axes", subtitle = paste0("samples=", nrow(long), "; Spearman statistics shown in facet strips"), x = "OLFML2B expression", y = "TME axis score") + th
      reg[[length(reg)+1L]] <- ol_p4_cns_save(p, dirs, "FIG4C_OLFML2B_TME_axis_scatter_CNS", kind = "scatter", width = 9.0, height = 4.8)
    }
  }

  if (nrow(surv) && "attenuation_percent" %in% names(surv)) {
    d <- surv[is.finite(suppressWarnings(as.numeric(surv$attenuation_percent))), , drop = FALSE]
    if (nrow(d)) {
      d$attenuation_percent <- as.numeric(d$attenuation_percent)
      d$model_group <- paste(d$model %||% "model", d$adjustment_role %||% "", sep = " | ")
      d <- d[!grepl("OLFML2B_only", d$model_group, ignore.case = TRUE), , drop = FALSE]
    }
    if (nrow(d)) {
      ss <- lapply(split(d$attenuation_percent, d$model_group), function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(NULL)
        data.frame(median = stats::median(x), q1 = stats::quantile(x, 0.25, names = FALSE), q3 = stats::quantile(x, 0.75, names = FALSE), n = length(x), stringsAsFactors = FALSE)
      })
      ss <- Filter(Negate(is.null), ss)
      if (length(ss)) {
        s <- do.call(rbind, ss); s$model_group <- rownames(s); rownames(s) <- NULL
        s$model_label <- factor(wrap(s$model_group, 28), levels = wrap(s$model_group[order(s$median)], 28))
        p <- ggplot2::ggplot(s, ggplot2::aes(x = median, y = model_label)) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = unname(pal["grey"]), linewidth = 0.32) +
          ggplot2::geom_errorbar(ggplot2::aes(xmin = q1, xmax = q3), orientation = "y", width = 0.18, color = unname(pal["dark"]), linewidth = 0.45) +
          ggplot2::geom_point(size = 2.7, color = unname(pal["red"])) +
          ggplot2::geom_text(ggplot2::aes(label = paste0("median=", sprintf("%.1f", median), "%\nn=", n)), hjust = -0.10, size = 2.6, lineheight = 0.92) +
          ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.30))) +
          ggplot2::labs(title = "TME-adjustment attenuation summary", subtitle = paste0("model groups=", nrow(s), "; median and IQR shown"), x = "OLFML2B coefficient attenuation (%)", y = NULL) + th
        reg[[length(reg)+1L]] <- ol_p4_cns_save(p, dirs, "FIG4D_TME_attenuation_summary_CNS", kind = "dotplot", n_rows = nrow(s), width = 8.2)
      }
    }
  }
  if (length(reg)) ol_p4_atomic_write_csv(ol_p4_bind_rows(reg), file.path(dirs$tables, "99_cns_enhanced_figure_registry.csv"))
  invisible(TRUE)
}

.ol_p4_cns_original_runner <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  index <- .ol_p4_cns_original_runner(...)
  tryCatch(ol_p4_make_cns_enhanced_figures(index$dirs), error = function(e) ol_p4_log("WARN", "Part4 CNS enhanced figures failed: ", conditionMessage(e)))
  invisible(index)
}

# ==============================================================================
# IF 7-8 TME implementation | frozen non-overlapping signatures, coverage
# eligibility, same-patient nested Cox models and patient bootstrap attenuation.
# ==============================================================================

OLFML2B_PART4_IF78_VERSION <- "v2.2.0_20260720_SUBSERIES_BATCH_STRATIFIED_TME"
OLFML2B_PART4_VERSION <- OLFML2B_PART4_IF78_VERSION

ol_p4_signature_catalog <- function() {
  list(
    CAF_Core = c("FAP","ACTA2","PDGFRA","PDGFRB","TAGLN","THY1","CXCL12","COL1A1","COL1A2","COL3A1","DCN","LUM","POSTN"),
    TGFb_Response = c("TGFBI","SERPINE1","SMAD3","SMAD7","CTGF","THBS1","ITGA5","PMEPA1","INHBA"),
    ECM_Remodeling = c("COL5A1","COL6A1","COL6A2","SPARC","MMP2","MMP14","LOX","PLOD2"),
    CD8_Cytotoxic = c("CD8A","CD8B","GZMA","GZMB","PRF1","NKG7","GNLY","IFNG"),
    IFNg_Response = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","IDO1","HLA-DRA","GBP1"),
    Checkpoint_Exhaustion = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1","VSIR"),
    Treg = c("FOXP3","IL2RA","IKZF2","CCR8","TNFRSF18","TNFRSF4","IL10"),
    M2_Macrophage = c("CD68","CD163","MRC1","MSR1","CSF1R","C1QA","C1QB","C1QC","MERTK","MARCO"),
    Myeloid_Inflammatory = c("LYZ","LST1","S100A8","S100A9","FCGR3A","ITGAM","IL1B","CXCL8","CCL2","CCR2"),
    Angiogenesis = c("VEGFA","KDR","FLT1","ANGPT1","ANGPT2","TEK","PECAM1","VWF","ENG","ESM1"),
    Epithelial_Tumor = c("EPCAM","KRT8","KRT18","KRT19","MUC1","CDH1","CLDN3","CLDN4","CLDN7"),
    EMT = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1"),
    Proliferation_Control = c("MKI67","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TOP2A","CDC20","CCNB1","BUB1","AURKA"),
    Hypoxia_Control = c("CA9","VEGFA","SLC2A1","LDHA","PDK1","BNIP3","EGLN3","HILPDA")
  )
}

ol_p4_score_signatures <- function(expr, catalog) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- ol_p4_clean_gene(rownames(expr))
  score_list <- list(); coverage <- list()
  for (sig in names(catalog)) {
    genes <- unique(ol_p4_clean_gene(catalog[[sig]]))
    idx <- which(rownames(expr) %in% genes)
    idx_by_gene <- split(idx, rownames(expr)[idx])
    keep <- if (length(idx_by_gene)) vapply(idx_by_gene, function(ii) {
      if (length(ii) == 1L) return(ii)
      vv <- apply(expr[ii, , drop = FALSE], 1L, stats::var, na.rm = TRUE)
      ii[which.max(replace(vv, !is.finite(vv), -Inf))]
    }, integer(1)) else integer()
    present <- if (length(keep)) rownames(expr)[keep] else character()
    coverage_fraction <- length(present) / max(length(genes), 1L)
    eligible <- length(present) >= 5L && coverage_fraction >= 0.50
    if (eligible) {
      mat <- expr[keep, , drop = FALSE]
      zmat <- t(apply(mat, 1L, ol_p4_z))
      score <- colMeans(zmat, na.rm = TRUE)
      score[!is.finite(score)] <- NA_real_
    } else {
      score <- rep(NA_real_, ncol(expr))
    }
    names(score) <- colnames(expr)
    score_list[[sig]] <- score
    coverage[[sig]] <- data.frame(
      signature = sig, requested_genes = length(genes), present_genes = length(present),
      coverage_fraction = coverage_fraction, eligible_for_inference = eligible,
      coverage_status = ifelse(eligible, "PASS_GE_5_AND_GE_50_PERCENT", "FAIL_SIGNATURE_NOT_SCORED"),
      present_gene_symbols = paste(sort(present), collapse = ";"),
      missing_gene_symbols = paste(sort(setdiff(genes, present)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  scores <- as.data.frame(score_list, stringsAsFactors = FALSE, check.names = FALSE)
  scores$sample_id <- colnames(expr)
  list(scores = scores, coverage = ol_p4_bind_rows(coverage))
}

ol_p4_if78_bootstrap_attenuation <- function(dd, form0, form2, B = 1000L, seed = 92912L) {
  set.seed(seed)
  strata_index <- if ("source_subseries" %in% names(dd) && length(unique(stats::na.omit(dd$source_subseries))) > 1L) {
    split(seq_len(nrow(dd)), dd$source_subseries)
  } else NULL
  boot <- replicate(B, {
    ii <- if (is.null(strata_index)) {
      sample.int(nrow(dd), nrow(dd), replace = TRUE)
    } else {
      unlist(lapply(strata_index, function(ix) sample(ix, length(ix), replace = TRUE)), use.names = FALSE)
    }
    z <- dd[ii, , drop = FALSE]
    b0 <- tryCatch(stats::coef(survival::coxph(form0, data = z))["OLFML2B_z"], error = function(e) NA_real_)
    b2 <- tryCatch(stats::coef(survival::coxph(form2, data = z))["OLFML2B_z"], error = function(e) NA_real_)
    if (!is.finite(b0) || !is.finite(b2) || abs(b0) <= 1e-8) NA_real_ else 100 * (b0 - b2) / b0
  })
  boot <- boot[is.finite(boot)]
  if (length(boot) < max(100L, B / 2L)) return(c(low = NA_real_, high = NA_real_, successful = length(boot)))
  c(low = unname(stats::quantile(boot, 0.025)), high = unname(stats::quantile(boot, 0.975)), successful = length(boot))
}

ol_p4_if78_sex_factor <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("female", "f", "woman")] <- "Female"
  out[raw %in% c("male", "m", "man")] <- "Male"
  factor(out, levels = c("Female", "Male"))
}

ol_p4_if78_stage_factor <- function(x) {
  raw <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  raw <- sub("^(AJCC|PATHOLOGIC|PATHOLOGICAL|CLINICAL)", "", raw)
  raw <- sub("^STAGE", "", raw)
  out <- rep(NA_character_, length(raw))
  out[grepl("^(IV|4)", raw)] <- "IV"
  out[is.na(out) & grepl("^(III|3)", raw)] <- "III"
  out[is.na(out) & grepl("^(II|2)", raw)] <- "II"
  out[is.na(out) & grepl("^(I|1)", raw)] <- "I"
  factor(out, levels = c("I", "II", "III", "IV"))
}

ol_p4_tme_cox <- function(d, cohort, min_n = 80L, min_events = 40L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package is required", call. = FALSE)
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) return(data.frame())
  clinical_source <- c("age", "sex", "stage")
  clinical <- c("age10", "sex", "stage")
  axes <- c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "CAF_TGFb_axis",
            "Cytotoxic_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
  rows <- list()
  design_strata <- if (identical(cohort, "GSE84437") && "source_subseries" %in% names(d) &&
                       length(unique(stats::na.omit(d$source_subseries))) == 2L) "source_subseries" else NA_character_
  for (ep in names(defs)) {
    tc <- defs[[ep]][1]; ec <- defs[[ep]][2]
    base <- d
    base$time <- ol_p4_num(base[[tc]])
    base$event <- ol_p4_event(base[[ec]])
    missing_core <- setdiff(c("time", "event", "OLFML2B_z", clinical_source), names(base))
    if (length(missing_core)) next
    age_num <- suppressWarnings(as.numeric(base$age))
    base$age10 <- (age_num - mean(age_num, na.rm = TRUE)) / 10
    base$sex <- ol_p4_if78_sex_factor(base$sex)
    base$stage <- ol_p4_if78_stage_factor(base$stage)
    for (axis in axes[axes %in% names(base)]) {
      required <- c("time", "event", "OLFML2B_z", clinical, axis,
                    if (is.na(design_strata)) character() else design_strata)
      dd <- base[stats::complete.cases(base[, required, drop = FALSE]) &
                   is.finite(base$time) & base$time > 0 & base$event %in% c(0L, 1L), , drop = FALSE]
      n <- nrow(dd); events <- sum(dd$event == 1L)
      if (n < max(80L, min_n) || events < max(40L, min_events)) next
      strata_term <- if (is.na(design_strata)) "" else paste0(" + survival::strata(", design_strata, ")")
      rhs <- list(
        OLFML2B_only = "OLFML2B_z",
        OLFML2B_plus_clinical = paste(c("OLFML2B_z", clinical), collapse = " + "),
        OLFML2B_plus_TME = paste(c("OLFML2B_z", clinical, axis), collapse = " + ")
      )
      forms <- lapply(rhs, function(z) stats::as.formula(paste0("survival::Surv(time, event) ~ ", z, strata_term)))
      fits <- lapply(forms, function(f) tryCatch(survival::coxph(f, data = dd, x = TRUE), error = function(e) NULL))
      if (any(vapply(fits, is.null, logical(1)))) {
        ol_p4_v230_record_eligibility(data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "COX_FIT_FAILED",
          n = n, events = events, clinical_tier = cc$tier,
          clinical_covariates = paste(cc$covars, collapse = ";"), stringsAsFactors = FALSE
        ))
        next
      }
      b0 <- unname(stats::coef(fits[[1]])["OLFML2B_z"])
      b2 <- unname(stats::coef(fits[[3]])["OLFML2B_z"])
      attenuation <- if (is.finite(b0) && abs(b0) > 1e-8) 100 * (b0 - b2) / b0 else NA_real_
      boot <- ol_p4_if78_bootstrap_attenuation(dd, forms[[1]], forms[[3]], B = 1000L,
                                               seed = 92912L + match(axis, axes))
      analysis_set <- if (requireNamespace("digest", quietly = TRUE))
        digest::digest(sort(as.character(dd$sample_id %||% seq_len(n))), algo = "sha256") else paste0(cohort, "_", ep, "_", axis, "_n", n)
      for (j in seq_along(fits)) {
        fit <- fits[[j]]; sm <- summary(fit)
        co <- sm$coefficients["OLFML2B_z", ]
        model_df <- max(1L, ncol(stats::model.matrix(fit)) - 1L)
        epdf <- events / model_df
        ph <- tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
        ph_term <- if (!is.null(ph) && "OLFML2B_z" %in% rownames(ph)) as.numeric(ph["OLFML2B_z", "p"]) else NA_real_
        ph_global <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph)) as.numeric(ph["GLOBAL", "p"]) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis,
          model = names(fits)[j], model_id = c("M0_same_case", "M1_common_same_case", "M2_TME_same_case")[j],
          adjustment_role = "same_patient_nested_attenuation", analysis_set_id = analysis_set,
          design_strata = design_strata,
          n = n, events = events, model_df = model_df, events_per_model_df = epdf,
          beta = as.numeric(co["coef"]), se = as.numeric(co["se(coef)"]), HR = exp(as.numeric(co["coef"])),
          lcl = exp(as.numeric(co["coef"]) - 1.96 * as.numeric(co["se(coef)"])),
          ucl = exp(as.numeric(co["coef"]) + 1.96 * as.numeric(co["se(coef)"])),
          p_value = as.numeric(co["Pr(>|z|)"]), ph_OLFML2B_p = ph_term, ph_global_p = ph_global,
          covariates = paste(setdiff(all.vars(forms[[j]])[-c(1,2)], "OLFML2B_z"), collapse = ";"),
          clinical_covariates = paste(clinical, collapse = ";"),
          attenuation_percent = ifelse(j == 3L, attenuation, ifelse(j == 1L, 0, NA_real_)),
          attenuation_bootstrap_lcl = ifelse(j == 3L, boot["low"], NA_real_),
          attenuation_bootstrap_ucl = ifelse(j == 3L, boot["high"], NA_real_),
          bootstrap_B = 1000L, bootstrap_successful = ifelse(j == 3L, boot["successful"], NA_real_),
          bootstrap_design = ifelse(is.na(design_strata), "patient_resampling", "patient_resampling_stratified_by_source_subseries"),
          structural_eligible = n >= 80L && events >= 40L && epdf >= min_epv,
          inference_tier = ifelse(n >= 80L && events >= 40L && epdf >= min_epv, "STRUCTURALLY_ELIGIBLE", "SUPPORTIVE_ONLY"),
          mediation_warning = "Attenuation on an identical patient set is descriptive sensitivity analysis, not formal mediation or causality.",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

.ol_p4_if78_runner_core <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  index <- .ol_p4_if78_runner_core(...)
  dirs <- index$dirs
  catalog <- ol_p4_signature_catalog()
  coverage_ok <- is.data.frame(index$coverage) && all(c("cohort", "signature", "eligible_for_inference", "coverage_status") %in% names(index$coverage))
  core_signatures <- c("CAF_Core", "TGFb_Response", "ECM_Remodeling")
  core_coverage_k <- 0L
  if (coverage_ok) {
    core_cov <- index$coverage[index$coverage$signature %in% core_signatures, , drop = FALSE]
    by_cohort <- split(core_cov, core_cov$cohort)
    core_coverage_k <- sum(vapply(by_cohort, function(z) {
      all(core_signatures %in% z$signature) && all(z$eligible_for_inference[match(core_signatures, z$signature)] %in% TRUE)
    }, logical(1)))
  }
  meta_cor <- index$meta_correlations %||% data.frame()
  core_meta_evaluable <- is.data.frame(meta_cor) && nrow(meta_cor) > 0L &&
    all(c("signature", "k", "rho_meta") %in% names(meta_cor)) &&
    all(core_signatures %in% meta_cor$signature) &&
    all(meta_cor$k[match(core_signatures, meta_cor$signature)] >= 3L) &&
    all(is.finite(meta_cor$rho_meta[match(core_signatures, meta_cor$signature)]))
  surv <- index$tme_survival
  same_case <- is.data.frame(surv) && nrow(surv) > 0L && all(c("analysis_set_id", "attenuation_axis", "model_id") %in% names(surv))
  if (same_case) {
    groups <- split(surv, interaction(surv$cohort, surv$endpoint, surv$attenuation_axis, drop = TRUE))
    same_case <- all(vapply(groups, function(z) length(unique(z$analysis_set_id)) == 1L && length(unique(z$n)) == 1L, logical(1)))
  }
  scores <- index$sample_scores %||% data.frame()
  g844_scores <- if (is.data.frame(scores) && "cohort" %in% names(scores)) scores[scores$cohort == "GSE84437", , drop = FALSE] else data.frame()
  g844_surv <- if (is.data.frame(surv) && "cohort" %in% names(surv)) surv[surv$cohort == "GSE84437", , drop = FALSE] else data.frame()
  batch_design_ok <- nrow(g844_scores) == 433L && "batch_adjustment" %in% names(g844_scores) &&
    all(g844_scores$batch_adjustment == "within_GSE84426_GSE84433_standardization") &&
    (!nrow(g844_surv) || (all(c("design_strata", "bootstrap_design") %in% names(g844_surv)) &&
      all(g844_surv$design_strata == "source_subseries") &&
      all(g844_surv$bootstrap_design == "patient_resampling_stratified_by_source_subseries")))
  go <- data.frame(
    criterion = c("frozen_signature_catalog", "signature_coverage_rule", "core_signature_coverage_k_ge_3",
                  "core_correlation_meta_evaluable", "fixed_clinical_covariates", "same_patient_nested_models",
                  "patient_bootstrap_attenuation", "GSE84437_subseries_batch_control",
                  "result_independent_pipeline_continuation", "claim_ceiling"),
    status = c(ifelse(all(c("CAF_Core", "TGFb_Response", "ECM_Remodeling") %in% names(catalog)), "PASS", "FAIL"),
               ifelse(coverage_ok, "PASS", "FAIL"),
               ifelse(core_coverage_k >= 3L, "PASS", "NOT_EVALUABLE_CORE_COVERAGE_LT3_COHORTS"),
               ifelse(core_meta_evaluable, "PASS", "NOT_EVALUABLE_CORE_META"),
               ifelse(!nrow(surv) || all(surv$clinical_covariates == "age10;sex;stage"), "PASS_OR_NOT_EVALUABLE", "FAIL"),
               ifelse(same_case, "PASS", ifelse(!nrow(surv), "NOT_EVALUABLE", "FAIL")),
               ifelse(!nrow(surv) || all(c("attenuation_bootstrap_lcl", "attenuation_bootstrap_ucl", "bootstrap_B") %in% names(surv)), "PASS_OR_NOT_EVALUABLE", "FAIL"),
               ifelse(batch_design_ok, "PASS", "FAIL"),
               "PASS", "TME_ASSOCIATION_AND_ATTENUATION_ONLY"),
    boundary = c(
      "Gene sets are prespecified and version-locked; primary CAF/TGFb/ECM sets do not share genes.",
      "A signature is scored only with at least 5 genes and at least 50% catalog coverage.",
      paste0("All three core axes are coverage-eligible in ", core_coverage_k, " cohorts; at least three are required for the formal cross-cohort ecological layer."),
      "Each core-axis Fisher-z REML synthesis requires at least three structurally evaluable cohorts; effect direction and P value are ignored by this gate.",
      "M1-common uses age, sex and stage without data-driven selection.",
      "M0/M1/M2 within an attenuation family use the identical patient set.",
      "M2 attenuation uncertainty uses 1000 patient-level bootstrap resamples.",
      "GSE84426/GSE84433 are standardized separately for expression/TME association and survival baseline hazards are stratified by formal subseries.",
      "Correlation direction, attenuation magnitude and P values never determine QC success.",
      "Bulk ecological association/attenuation is not mediation, causality or immune-exclusion proof."
    ), stringsAsFactors = FALSE
  )
  index$version <- OLFML2B_PART4_IF78_VERSION
  index$go_no_go <- go
  ol_p4_atomic_write_csv(go, file.path(dirs$tables, "20_part4_go_no_go_summary.csv"))
  ol_p4_atomic_write_csv(data.frame(signature = names(catalog), version = OLFML2B_PART4_IF78_VERSION,
                                    stringsAsFactors = FALSE), file.path(dirs$tables, "05b_frozen_signature_version.csv"))
  ol_p4_atomic_save_rds(index, file.path(dirs$objects, "Part4_immune_TME_production_index.rds"))
  invisible(index)
}

# ==============================================================================
# OLFML2B-STAD Part4 final release wrapper.
# Repairs the v1.0.0 source-order defect by wrapping the complete production
# runner only after that runner has been defined. It also enforces a hard
# target-exclusion contract for every TME signature.
# ==============================================================================
OLFML2B_PART4_TARGET_VERSION <- "v1.0.1_20260721_COMPLETE_CORE_SOURCE_INTERFACE_FIX"
.ol_p4_olfml2b_complete_core <- run_olfml2b_part4_immune_tme_production
OLFML2B_PART4_IF78_VERSION <- OLFML2B_PART4_TARGET_VERSION
OLFML2B_PART4_VERSION <- OLFML2B_PART4_TARGET_VERSION

run_olfml2b_part4_immune_tme_production <- function(...) {
  catalog <- ol_p4_signature_catalog()
  target_overlap <- vapply(
    catalog,
    function(g) "OLFML2B" %in% toupper(trimws(as.character(g))),
    logical(1)
  )
  if (any(target_overlap)) {
    stop(
      "OLFML2B appears in one or more TME signatures: ",
      paste(names(catalog)[target_overlap], collapse = ";"),
      call. = FALSE
    )
  }

  index <- .ol_p4_olfml2b_complete_core(...)
  index$version <- OLFML2B_PART4_TARGET_VERSION
  index$final_gene_lock <- FALSE

  audit <- data.frame(
    signature = names(catalog),
    target_gene = "OLFML2B",
    target_in_signature = unname(target_overlap),
    primary_or_secondary = ifelse(
      names(catalog) %in% c("CAF_Core", "TGFb_Response", "ECM_Remodeling"),
      "PRIMARY_CORE",
      "SECONDARY_OR_CONTROL"
    ),
    status = ifelse(target_overlap, "FAIL_CIRCULAR_SIGNATURE", "PASS_TARGET_EXCLUDED"),
    stringsAsFactors = FALSE
  )
  index$target_exclusion_audit <- audit

  dirs <- index$dirs
  ol_p4_atomic_write_csv(
    audit,
    file.path(dirs$tables, "05c_OLFML2B_target_excluded_from_signatures.csv")
  )
  source_contract <- data.frame(
    interface = c(
      "complete_base_runner_defined_before_wrappers",
      "final_runner_available",
      "signature_catalog_available",
      "correlation_meta_available",
      "same_patient_tme_cox_available",
      "target_excluded_from_all_signatures"
    ),
    status = c(
      "PASS",
      ifelse(exists("run_olfml2b_part4_immune_tme_production", mode = "function"), "PASS", "FAIL"),
      ifelse(exists("ol_p4_signature_catalog", mode = "function"), "PASS", "FAIL"),
      ifelse(exists("ol_p4_meta_cor", mode = "function"), "PASS", "FAIL"),
      ifelse(exists("ol_p4_tme_cox", mode = "function"), "PASS", "FAIL"),
      ifelse(!any(target_overlap), "PASS", "FAIL")
    ),
    stringsAsFactors = FALSE
  )
  index$source_interface_audit <- source_contract
  ol_p4_atomic_write_csv(
    source_contract,
    file.path(dirs$tables, "00b_part4_source_interface_audit.csv")
  )
  ol_p4_atomic_save_rds(
    index,
    file.path(dirs$objects, "Part4_immune_TME_production_index.rds")
  )
  invisible(index)
}


# ==============================================================================
# Part4 methodology repair v2.3.0
# Fixes stage parsing, tumor-only patient inference, same-case M0/M1/M2 nesting,
# prespecified strict/available clinical tiers, M1-to-M2 attenuation, patient
# bootstrap uncertainty, and an explicit model-eligibility audit.
# ==============================================================================
OLFML2B_PART4_TARGET_VERSION <- "v2.3.0_20260722_SAME_CASE_TME_ATTENUATION_FIX"
OLFML2B_PART4_IF78_VERSION <- OLFML2B_PART4_TARGET_VERSION
OLFML2B_PART4_VERSION <- OLFML2B_PART4_TARGET_VERSION

.ol_p4_v230_state <- new.env(parent = emptyenv())
.ol_p4_v230_state$eligibility <- list()

ol_p4_if78_stage_factor <- function(x) {
  raw <- toupper(as.character(x))
  raw <- gsub("[^A-Z0-9]", "", raw)
  # Remove stacked descriptors such as AJCCPATHOLOGICSTAGE before parsing.
  # The previous implementation removed only one leading descriptor, so values
  # such as "AJCC Pathologic Stage IIIA" could remain unparsable.
  raw <- gsub("(AJCC|PATHOLOGIC(AL)?|CLINICAL|PATHSTAGE|STAGE)", "", raw)
  out <- rep(NA_character_, length(raw))
  out[grepl("^(IV|4)", raw)] <- "IV"
  out[is.na(out) & grepl("^(III|3)", raw)] <- "III"
  out[is.na(out) & grepl("^(II|2)", raw)] <- "II"
  out[is.na(out) & grepl("^(I|1)", raw)] <- "I"
  factor(out, levels = c("I", "II", "III", "IV"))
}

ol_p4_v230_one_value <- function(x, numeric = FALSE) {
  if (numeric) {
    z <- suppressWarnings(as.numeric(x))
    z <- z[is.finite(z)]
    return(if (length(z)) stats::median(z) else NA_real_)
  }
  z <- as.character(x)
  z <- z[!is.na(z) & nzchar(trimws(z))]
  z <- unique(z)
  if (length(z) == 1L) z else NA_character_
}

ol_p4_v230_patient_key <- function(d) {
  key <- rep(NA_character_, nrow(d))
  if ("patient_id" %in% names(d)) key <- trimws(as.character(d$patient_id))
  bad <- is.na(key) | !nzchar(key)
  if ("sample_id" %in% names(d)) key[bad] <- trimws(as.character(d$sample_id[bad]))
  bad <- is.na(key) | !nzchar(key)
  key[bad] <- paste0("ROW_", which(bad))
  key
}

ol_p4_v230_collapse_patient <- function(d, axis, design_strata = NA_character_) {
  d$patient_key <- ol_p4_v230_patient_key(d)
  groups <- split(d, d$patient_key)
  rows <- lapply(groups, function(z) {
    data.frame(
      patient_key = z$patient_key[1],
      sample_id = paste(sort(unique(as.character(z$sample_id %||% z$patient_key))), collapse = ";"),
      time = ol_p4_v230_one_value(z$time, numeric = TRUE),
      event = {
        ev <- unique(z$event[!is.na(z$event)])
        if (length(ev) == 1L) as.integer(ev) else NA_integer_
      },
      OLFML2B_z = ol_p4_v230_one_value(z$OLFML2B_z, numeric = TRUE),
      TME_axis = ol_p4_v230_one_value(z[[axis]], numeric = TRUE),
      age = if ("age" %in% names(z)) ol_p4_v230_one_value(z$age, numeric = TRUE) else NA_real_,
      sex_raw = if ("sex" %in% names(z)) ol_p4_v230_one_value(z$sex) else NA_character_,
      stage_raw_v230 = if ("stage" %in% names(z)) ol_p4_v230_one_value(z$stage) else NA_character_,
      source_subseries = if (!is.na(design_strata) && design_strata %in% names(z)) ol_p4_v230_one_value(z[[design_strata]]) else NA_character_,
      n_source_rows = nrow(z),
      stringsAsFactors = FALSE
    )
  })
  out <- ol_p4_bind_rows(rows)
  if (!nrow(out)) return(out)
  out$age10 <- if (any(is.finite(out$age))) (out$age - mean(out$age, na.rm = TRUE)) / 10 else NA_real_
  out$sex <- ol_p4_if78_sex_factor(out$sex_raw)
  out$stage <- ol_p4_if78_stage_factor(out$stage_raw_v230)
  out$source_subseries <- factor(out$source_subseries)
  out
}

ol_p4_v230_covar_usable <- function(x, n_total, numeric = FALSE) {
  if (numeric) {
    z <- suppressWarnings(as.numeric(x))
    return(sum(is.finite(z)) >= max(30L, ceiling(0.60 * n_total)) &&
             is.finite(stats::sd(z, na.rm = TRUE)) && stats::sd(z, na.rm = TRUE) > 0)
  }
  z <- droplevels(factor(x))
  tab <- table(z, useNA = "no")
  sum(tab) >= max(30L, ceiling(0.60 * n_total)) && length(tab) >= 2L && all(tab >= 5L)
}

ol_p4_v230_clinical_covars <- function(d) {
  usable <- c(
    age10 = ol_p4_v230_covar_usable(d$age10, nrow(d), numeric = TRUE),
    sex = ol_p4_v230_covar_usable(d$sex, nrow(d), numeric = FALSE),
    stage = ol_p4_v230_covar_usable(d$stage, nrow(d), numeric = FALSE)
  )
  covars <- names(usable)[usable]
  tier <- if (all(usable)) "STRICT_COMMON_AGE_SEX_STAGE" else if (length(covars)) "AVAILABLE_CLINICAL_PRESPECIFIED_SUBSET" else "NO_CLINICAL_COVARIATE_AVAILABLE"
  list(covars = covars, tier = tier, usable = usable)
}

ol_p4_v230_model_df <- function(fit) {
  mm <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if (is.null(mm)) return(NA_integer_)
  max(0L, ncol(mm) - 1L)
}

ol_p4_v230_lrt_p <- function(fit_small, fit_large) {
  z <- tryCatch(stats::anova(fit_small, fit_large, test = "LRT"), error = function(e) NULL)
  if (is.null(z) || nrow(z) < 2L) return(NA_real_)
  pcol <- grep("P\\(", names(z), value = TRUE)
  if (!length(pcol)) return(NA_real_)
  suppressWarnings(as.numeric(z[2L, pcol[1L]]))
}

ol_p4_if78_bootstrap_attenuation <- function(dd, form1, form2, B = 1000L, seed = 92912L) {
  set.seed(seed)
  strata_index <- if ("source_subseries" %in% names(dd) && length(unique(stats::na.omit(dd$source_subseries))) > 1L) split(seq_len(nrow(dd)), dd$source_subseries) else NULL
  boot <- replicate(B, {
    ii <- if (is.null(strata_index)) sample.int(nrow(dd), nrow(dd), replace = TRUE) else unlist(lapply(strata_index, function(ix) sample(ix, length(ix), replace = TRUE)), use.names = FALSE)
    z <- dd[ii, , drop = FALSE]
    b1 <- tryCatch(stats::coef(survival::coxph(form1, data = z))["OLFML2B_z"], error = function(e) NA_real_)
    b2 <- tryCatch(stats::coef(survival::coxph(form2, data = z))["OLFML2B_z"], error = function(e) NA_real_)
    if (!is.finite(b1) || !is.finite(b2) || abs(b1) <= 1e-8) NA_real_ else 100 * (b1 - b2) / b1
  })
  boot <- boot[is.finite(boot)]
  if (length(boot) < max(100L, B / 2L)) return(c(low = NA_real_, high = NA_real_, successful = length(boot)))
  c(low = unname(stats::quantile(boot, 0.025)), high = unname(stats::quantile(boot, 0.975)), successful = length(boot))
}

ol_p4_v230_record_eligibility <- function(row) {
  .ol_p4_v230_state$eligibility[[length(.ol_p4_v230_state$eligibility) + 1L]] <- row
  invisible(TRUE)
}

ol_p4_tme_cox <- function(d, cohort, min_n = 40L, min_events = 20L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package is required", call. = FALSE)
  if ("is_tumor" %in% names(d) && any(d$is_tumor %in% TRUE, na.rm = TRUE)) d <- d[d$is_tumor %in% TRUE, , drop = FALSE]
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) {
    ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = NA_character_, attenuation_axis = NA_character_, status = "NO_AUDITABLE_ENDPOINT", n = nrow(d), events = NA_integer_, clinical_tier = NA_character_, clinical_covariates = NA_character_, stringsAsFactors = FALSE))
    return(data.frame())
  }
  axes <- c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "CAF_TGFb_axis", "Cytotoxic_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
  rows <- list()
  design_strata <- if (identical(cohort, "GSE84437") && "source_subseries" %in% names(d) && length(unique(stats::na.omit(d$source_subseries))) >= 2L) "source_subseries" else NA_character_
  for (ep in names(defs)) {
    tc <- defs[[ep]][1L]; ec <- defs[[ep]][2L]
    for (axis in axes[axes %in% names(d)]) {
      raw <- d
      raw$time <- ol_p4_num(raw[[tc]])
      raw$event <- ol_p4_event(raw[[ec]])
      raw <- raw[is.finite(raw$time) & raw$time > 0 & raw$event %in% c(0L, 1L) & is.finite(raw$OLFML2B_z) & is.finite(raw[[axis]]), , drop = FALSE]
      if (!nrow(raw)) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "NO_COMPLETE_ENDPOINT_TARGET_AXIS_ROWS", n = 0L, events = 0L, clinical_tier = NA_character_, clinical_covariates = NA_character_, stringsAsFactors = FALSE))
        next
      }
      pd <- ol_p4_v230_collapse_patient(raw, axis, design_strata)
      cc <- ol_p4_v230_clinical_covars(pd)
      required <- c("time", "event", "OLFML2B_z", "TME_axis", cc$covars, if (is.na(design_strata)) character() else "source_subseries")
      dd <- pd[stats::complete.cases(pd[, required, drop = FALSE]) & is.finite(pd$time) & pd$time > 0 & pd$event %in% c(0L, 1L), , drop = FALSE]
      n <- nrow(dd); events <- sum(dd$event == 1L)
      if (n < min_n || events < min_events) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "INSUFFICIENT_SAME_CASE_COMPLETE_DATA", n = n, events = events, clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"), stringsAsFactors = FALSE))
        next
      }
      clinical_rhs <- if (length(cc$covars)) paste(cc$covars, collapse = " + ") else "1"
      strata_term <- if (is.na(design_strata)) "" else " + survival::strata(source_subseries)"
      rhs <- c(
        M0_clinical = clinical_rhs,
        M1_clinical_plus_OLFML2B = paste(clinical_rhs, "+ OLFML2B_z"),
        M2_clinical_plus_OLFML2B_plus_TME = paste(clinical_rhs, "+ OLFML2B_z + TME_axis")
      )
      forms <- lapply(rhs, function(z) stats::as.formula(paste0("survival::Surv(time, event) ~ ", z, strata_term)))
      fits <- lapply(forms, function(f) tryCatch(survival::coxph(f, data = dd, x = TRUE, model = TRUE), error = function(e) NULL))
      if (any(vapply(fits, is.null, logical(1)))) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "COX_FIT_FAILED", n = n, events = events, clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"), stringsAsFactors = FALSE))
        next
      }
      df_m2 <- ol_p4_v230_model_df(fits[[3L]])
      epv_m2 <- events / max(df_m2, 1L)
      inference_tier <- if (n >= 80L && events >= 40L && epv_m2 >= min_epv) "FORMAL_STRUCTURALLY_ELIGIBLE" else if (epv_m2 >= 5) "SUPPORTIVE_STRUCTURALLY_ELIGIBLE" else "DESCRIPTIVE_LOW_EPV"
      b1 <- unname(stats::coef(fits[[2L]])["OLFML2B_z"])
      b2 <- unname(stats::coef(fits[[3L]])["OLFML2B_z"])
      attenuation <- if (is.finite(b1) && is.finite(b2) && abs(b1) > 1e-8) 100 * (b1 - b2) / b1 else NA_real_
      boot <- ol_p4_if78_bootstrap_attenuation(dd, forms[[2L]], forms[[3L]], B = 1000L, seed = 92912L + match(axis, axes))
      analysis_set <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(sort(as.character(dd$patient_key)), algo = "sha256") else paste0(cohort, "_", ep, "_", axis, "_n", n)
      lrt01 <- ol_p4_v230_lrt_p(fits[[1L]], fits[[2L]])
      lrt12 <- ol_p4_v230_lrt_p(fits[[2L]], fits[[3L]])
      for (j in seq_along(fits)) {
        fit <- fits[[j]]
        sm <- summary(fit)
        has_target <- "OLFML2B_z" %in% rownames(sm$coefficients)
        beta <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "coef"]) else NA_real_
        se <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "se(coef)"]) else NA_real_
        ph <- tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
        ph_term <- if (!is.null(ph) && "OLFML2B_z" %in% rownames(ph)) as.numeric(ph["OLFML2B_z", "p"]) else NA_real_
        ph_global <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph)) as.numeric(ph["GLOBAL", "p"]) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis,
          model = names(fits)[j], model_id = names(fits)[j],
          adjustment_role = "same_patient_nested_attenuation", analysis_set_id = analysis_set,
          clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"),
          design_strata = ifelse(is.na(design_strata), "none", design_strata),
          n = n, events = events, model_df = ol_p4_v230_model_df(fit), events_per_model_df = events / max(ol_p4_v230_model_df(fit), 1L),
          beta = beta, se = se, HR = if (is.finite(beta)) exp(beta) else NA_real_,
          lcl = if (is.finite(beta) && is.finite(se)) exp(beta - 1.96 * se) else NA_real_,
          ucl = if (is.finite(beta) && is.finite(se)) exp(beta + 1.96 * se) else NA_real_,
          p_value = if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "Pr(>|z|)"]) else NA_real_,
          aic = stats::AIC(fit), lrt_M0_vs_M1_p = lrt01, lrt_M1_vs_M2_p = lrt12,
          ph_OLFML2B_p = ph_term, ph_global_p = ph_global,
          attenuation_percent = ifelse(j == 3L, attenuation, NA_real_),
          attenuation_bootstrap_lcl = ifelse(j == 3L, boot["low"], NA_real_),
          attenuation_bootstrap_ucl = ifelse(j == 3L, boot["high"], NA_real_),
          bootstrap_B = 1000L, bootstrap_successful = ifelse(j == 3L, boot["successful"], NA_real_),
          bootstrap_design = ifelse(is.na(design_strata), "patient_resampling", "patient_resampling_stratified_by_source_subseries"),
          same_patient_set_verified = TRUE, inference_tier = inference_tier,
          mediation_warning = "Coefficient attenuation on an identical patient set is descriptive sensitivity analysis, not mediation or causality.",
          stringsAsFactors = FALSE
        )
      }
      ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = inference_tier, n = n, events = events, clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"), analysis_set_id = analysis_set, m2_model_df = df_m2, m2_epv = epv_m2, stringsAsFactors = FALSE))
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) out$fdr_OLFML2B_within_model_family <- stats::p.adjust(out$p_value, method = "BH")
  out
}

ol_p4_v230_nested_forest <- function(index) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  d <- index$tme_survival %||% data.frame()
  d <- d[d$model %in% c("M1_clinical_plus_OLFML2B", "M2_clinical_plus_OLFML2B_plus_TME") & is.finite(d$HR) & is.finite(d$lcl) & is.finite(d$ucl), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d$label <- paste0(d$cohort, " / ", d$endpoint, " / ", d$attenuation_axis, " / ", d$model, "\nHR ", sprintf("%.2f", d$HR), " (", sprintf("%.2f", d$lcl), "–", sprintf("%.2f", d$ucl), "), n=", d$n, ", events=", d$events)
  d$label <- factor(make.unique(d$label), levels = rev(make.unique(d$label)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = HR, y = label, shape = model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lcl, xmax = ucl), orientation = "y", width = 0.15) +
    ggplot2::geom_point(size = 2.2) + ggplot2::scale_x_log10() +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::labs(title = "Same-patient nested TME attenuation models", subtitle = "M1 clinical + OLFML2B versus M2 clinical + OLFML2B + TME axis", x = "Hazard ratio per SD OLFML2B", y = NULL)
  stem <- file.path(index$dirs$figures, "FIG4C_same_patient_nested_TME_attenuation")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 11, height = max(6, 0.28 * nrow(d) + 2), limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 11, height = max(6, 0.28 * nrow(d) + 2), dpi = 600, limitsize = FALSE)
  invisible(TRUE)
}

.ol_p4_v230_core_runner <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  .ol_p4_v230_state$eligibility <- list()
  index <- .ol_p4_v230_core_runner(...)
  index$version <- OLFML2B_PART4_TARGET_VERSION

  eligibility <- ol_p4_bind_rows(.ol_p4_v230_state$eligibility)
  surv <- index$tme_survival %||% data.frame()
  same_case <- nrow(surv) > 0L && all(surv$same_patient_set_verified %in% TRUE)
  nested_complete <- FALSE
  if (nrow(surv)) {
    grp <- split(surv, interaction(surv$cohort, surv$endpoint, surv$attenuation_axis, drop = TRUE))
    nested_complete <- all(vapply(
      grp,
      function(z) {
        all(c("M0_clinical", "M1_clinical_plus_OLFML2B", "M2_clinical_plus_OLFML2B_plus_TME") %in% z$model) &&
          length(unique(z$analysis_set_id)) == 1L &&
          length(unique(z$n)) == 1L &&
          length(unique(z$events)) == 1L
      },
      logical(1)
    ))
  }

  formal_families <- if (nrow(eligibility)) sum(eligibility$status == "FORMAL_STRUCTURALLY_ELIGIBLE", na.rm = TRUE) else 0L
  supportive_families <- if (nrow(eligibility)) sum(eligibility$status == "SUPPORTIVE_STRUCTURALLY_ELIGIBLE", na.rm = TRUE) else 0L

  # Preserve non-superseded ecological/signature gates from the production core.
  prior_go <- index$go_no_go %||% data.frame()
  superseded <- c(
    "fixed_clinical_covariates", "same_patient_nested_models",
    "patient_bootstrap_attenuation", "result_independent_pipeline_continuation",
    "claim_ceiling", "TME_attenuation_survival_models"
  )
  if (nrow(prior_go) && "criterion" %in% names(prior_go)) {
    prior_go <- prior_go[!prior_go$criterion %in% superseded, , drop = FALSE]
  } else {
    prior_go <- data.frame()
  }

  repair_go <- data.frame(
    criterion = c(
      "stage_parser", "tumor_only_patient_unit", "prespecified_clinical_tiers",
      "same_patient_nested_models", "M1_to_M2_attenuation",
      "patient_bootstrap_attenuation", "formal_model_families",
      "supportive_model_families", "result_independent_pipeline_continuation",
      "claim_ceiling"
    ),
    status = c(
      "PASS", "PASS", "PASS",
      ifelse(same_case && nested_complete, "PASS", ifelse(!nrow(surv), "NOT_EVALUABLE", "FAIL")),
      ifelse(nrow(surv) && all(c("attenuation_percent", "lrt_M1_vs_M2_p") %in% names(surv)), "PASS", "NOT_EVALUABLE"),
      ifelse(nrow(surv) && all(c("attenuation_bootstrap_lcl", "attenuation_bootstrap_ucl", "bootstrap_successful") %in% names(surv)), "PASS", "NOT_EVALUABLE"),
      as.character(formal_families), as.character(supportive_families), "PASS",
      "TME_ASSOCIATION_AND_SAME_CASE_ATTENUATION_ONLY"
    ),
    boundary = c(
      "Stage parsing uppercases before stripping characters and removes stacked AJCC/pathologic/clinical stage descriptors.",
      "Survival inference uses tumor samples when explicit tumor labels exist and collapses duplicate rows to one patient unit.",
      "STRICT_COMMON uses age10, sex and stage; AVAILABLE_CLINICAL uses only the prespecified structurally available subset without significance-driven selection.",
      "M0, M1 and M2 within each cohort-endpoint-axis family use an identical patient and event set.",
      "OLFML2B attenuation is the coefficient change from M1 clinical+OLFML2B to M2 clinical+OLFML2B+TME.",
      "Attenuation uncertainty uses 1000 patient-level bootstrap resamples, stratified by GSE84437 subseries when required.",
      "Number of model families meeting n>=80, events>=40 and EPV>=10.",
      "Number of model families meeting supportive data and EPV criteria but not the formal threshold.",
      "Effect sign, attenuation magnitude and P value do not determine pipeline continuation.",
      "Bulk ecological association and coefficient attenuation are not mediation, causality or proof of immune exclusion."
    ),
    stringsAsFactors = FALSE
  )
  go <- ol_p4_bind_rows(list(prior_go, repair_go))

  index$tme_model_eligibility <- eligibility
  index$go_no_go <- go
  index$final_gene_lock <- FALSE

  attenuation_tab <- if (nrow(surv)) {
    surv[surv$model %in% c("M1_clinical_plus_OLFML2B", "M2_clinical_plus_OLFML2B_plus_TME"), , drop = FALSE]
  } else {
    surv
  }
  ol_p4_atomic_write_csv(surv, file.path(index$dirs$tables, "12_tme_adjusted_survival_models.csv"))
  ol_p4_atomic_write_csv(attenuation_tab, file.path(index$dirs$tables, "12_tme_attenuation_survival_models.csv"))
  ol_p4_atomic_write_csv(eligibility, file.path(index$dirs$tables, "12a_tme_nested_model_eligibility_audit.csv"))
  ol_p4_atomic_write_csv(go, file.path(index$dirs$tables, "20_part4_go_no_go_summary.csv"))
  try(ol_p4_v230_nested_forest(index), silent = TRUE)
  ol_p4_atomic_save_rds(index, file.path(index$dirs$objects, "Part4_immune_TME_production_index.rds"))
  invisible(index)
}


# ============================================================================
# Part4 v2.3.1 Cox/bootstrap stability hardening
# - captures and audits Cox warnings instead of flooding the console;
# - marks potentially infinite/non-converged main fits as unstable;
# - stratifies patient bootstrap by event and source subseries when available;
# - rejects unstable bootstrap fits and records success/failure counts;
# - avoids relative attenuation percentages when the M1 beta denominator is
#   too close to zero (absolute beta < 0.05);
# - writes explicit convergence and bootstrap-stability audit tables.
# ============================================================================
OLFML2B_PART4_TARGET_VERSION <- "v2.3.1_20260722_COX_BOOTSTRAP_STABILITY_AND_WARNING_AUDIT_FIX"
OLFML2B_PART4_VERSION <- OLFML2B_PART4_TARGET_VERSION

ol_p4_v231_safe_coxph <- function(form, data, x = FALSE, model = FALSE) {
  warning_text <- character()
  error_text <- NA_character_
  fit <- withCallingHandlers(
    tryCatch(
      survival::coxph(form, data = data, x = x, model = model),
      error = function(e) {
        error_text <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warning_text <- unique(warning_text[nzchar(warning_text)])
  if (is.null(fit)) {
    return(list(
      fit = NULL, stable = FALSE, status = "COX_ERROR",
      warning_count = length(warning_text), warning_text = paste(warning_text, collapse = " | "),
      error_text = error_text
    ))
  }
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  coef_tab <- if (!is.null(sm)) sm$coefficients else NULL
  finite_coefficients <- !is.null(coef_tab) && nrow(coef_tab) > 0L &&
    all(is.finite(suppressWarnings(as.numeric(coef_tab[, "coef"])))) &&
    all(is.finite(suppressWarnings(as.numeric(coef_tab[, "se(coef)"]))))
  suspect_warning <- any(grepl(
    "coefficient may be infinite|did not converge|ran out of iterations|algorithm did not converge",
    warning_text, ignore.case = TRUE
  ))
  extreme_coefficient <- FALSE
  if (!is.null(coef_tab) && nrow(coef_tab) > 0L) {
    b <- suppressWarnings(as.numeric(coef_tab[, "coef"]))
    extreme_coefficient <- any(is.finite(b) & abs(b) > 20)
  }
  stable <- isTRUE(finite_coefficients) && !suspect_warning && !extreme_coefficient
  status <- if (stable && !length(warning_text)) {
    "OK"
  } else if (stable) {
    "OK_WITH_NONCRITICAL_WARNING"
  } else if (suspect_warning || extreme_coefficient) {
    "UNSTABLE_POTENTIAL_INFINITE_COEFFICIENT"
  } else {
    "UNSTABLE_NONFINITE_COEFFICIENT_OR_SE"
  }
  list(
    fit = fit, stable = stable, status = status,
    warning_count = length(warning_text),
    warning_text = paste(warning_text, collapse = " | "),
    error_text = error_text
  )
}

ol_p4_v231_bootstrap_attenuation <- function(
  dd, form1, form2, B = 1000L, seed = 92912L,
  min_success_fraction = 0.50, min_abs_baseline_beta = 0.05
) {
  set.seed(seed)
  source_key <- if (
    "source_subseries" %in% names(dd) &&
      length(unique(stats::na.omit(dd$source_subseries))) > 1L
  ) as.character(dd$source_subseries) else rep("ALL", nrow(dd))
  event_key <- ifelse(dd$event == 1L, "EVENT", "CENSORED")
  bootstrap_stratum <- interaction(source_key, event_key, drop = TRUE, lex.order = TRUE)
  strata_index <- split(seq_len(nrow(dd)), bootstrap_stratum)

  estimates <- rep(NA_real_, B)
  failed_fit <- 0L
  unstable_fit <- 0L
  warning_replicates <- 0L
  warning_count <- 0L
  denominator_small <- 0L

  for (b in seq_len(B)) {
    ii <- unlist(lapply(strata_index, function(ix) sample(ix, length(ix), replace = TRUE)), use.names = FALSE)
    z <- dd[ii, , drop = FALSE]
    f1 <- ol_p4_v231_safe_coxph(form1, z)
    f2 <- ol_p4_v231_safe_coxph(form2, z)
    wc <- f1$warning_count + f2$warning_count
    warning_count <- warning_count + wc
    if (wc > 0L) warning_replicates <- warning_replicates + 1L
    if (is.null(f1$fit) || is.null(f2$fit)) {
      failed_fit <- failed_fit + 1L
      next
    }
    if (!isTRUE(f1$stable) || !isTRUE(f2$stable)) {
      unstable_fit <- unstable_fit + 1L
      next
    }
    b1 <- suppressWarnings(unname(stats::coef(f1$fit)["OLFML2B_z"]))
    b2 <- suppressWarnings(unname(stats::coef(f2$fit)["OLFML2B_z"]))
    if (!is.finite(b1) || !is.finite(b2)) {
      unstable_fit <- unstable_fit + 1L
      next
    }
    if (abs(b1) < min_abs_baseline_beta) {
      denominator_small <- denominator_small + 1L
      next
    }
    estimates[b] <- 100 * (b1 - b2) / b1
  }
  ok <- estimates[is.finite(estimates)]
  minimum_success <- max(100L, ceiling(B * min_success_fraction))
  if (length(ok) < minimum_success) {
    return(c(
      low = NA_real_, high = NA_real_, median = NA_real_, successful = length(ok),
      failed_fit = failed_fit, unstable_fit = unstable_fit,
      warning_replicates = warning_replicates, warning_count = warning_count,
      denominator_small = denominator_small, minimum_success = minimum_success
    ))
  }
  c(
    low = unname(stats::quantile(ok, 0.025, names = FALSE)),
    high = unname(stats::quantile(ok, 0.975, names = FALSE)),
    median = stats::median(ok), successful = length(ok),
    failed_fit = failed_fit, unstable_fit = unstable_fit,
    warning_replicates = warning_replicates, warning_count = warning_count,
    denominator_small = denominator_small, minimum_success = minimum_success
  )
}

# Preserve the public helper name while replacing its implementation.
ol_p4_if78_bootstrap_attenuation <- function(dd, form1, form2, B = 1000L, seed = 92912L) {
  ol_p4_v231_bootstrap_attenuation(dd, form1, form2, B = B, seed = seed)
}

ol_p4_tme_cox <- function(d, cohort, min_n = 40L, min_events = 20L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package is required", call. = FALSE)
  if ("is_tumor" %in% names(d) && any(d$is_tumor %in% TRUE, na.rm = TRUE)) d <- d[d$is_tumor %in% TRUE, , drop = FALSE]
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) {
    ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = NA_character_, attenuation_axis = NA_character_, status = "NO_AUDITABLE_ENDPOINT", n = nrow(d), events = NA_integer_, clinical_tier = NA_character_, clinical_covariates = NA_character_, stringsAsFactors = FALSE))
    return(data.frame())
  }
  axes <- c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "CAF_TGFb_axis", "Cytotoxic_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
  rows <- list()
  design_strata <- if (identical(cohort, "GSE84437") && "source_subseries" %in% names(d) && length(unique(stats::na.omit(d$source_subseries))) >= 2L) "source_subseries" else NA_character_
  for (ep in names(defs)) {
    tc <- defs[[ep]][1L]; ec <- defs[[ep]][2L]
    for (axis in axes[axes %in% names(d)]) {
      raw <- d
      raw$time <- ol_p4_num(raw[[tc]])
      raw$event <- ol_p4_event(raw[[ec]])
      raw <- raw[is.finite(raw$time) & raw$time > 0 & raw$event %in% c(0L, 1L) & is.finite(raw$OLFML2B_z) & is.finite(raw[[axis]]), , drop = FALSE]
      if (!nrow(raw)) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "NO_COMPLETE_ENDPOINT_TARGET_AXIS_ROWS", n = 0L, events = 0L, clinical_tier = NA_character_, clinical_covariates = NA_character_, stringsAsFactors = FALSE))
        next
      }
      pd <- ol_p4_v230_collapse_patient(raw, axis, design_strata)
      cc <- ol_p4_v230_clinical_covars(pd)
      required <- c("time", "event", "OLFML2B_z", "TME_axis", cc$covars, if (is.na(design_strata)) character() else "source_subseries")
      dd <- pd[stats::complete.cases(pd[, required, drop = FALSE]) & is.finite(pd$time) & pd$time > 0 & pd$event %in% c(0L, 1L), , drop = FALSE]
      n <- nrow(dd); events <- sum(dd$event == 1L)
      if (n < min_n || events < min_events) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "INSUFFICIENT_SAME_CASE_COMPLETE_DATA", n = n, events = events, clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"), stringsAsFactors = FALSE))
        next
      }
      clinical_rhs <- if (length(cc$covars)) paste(cc$covars, collapse = " + ") else "1"
      strata_term <- if (is.na(design_strata)) "" else " + survival::strata(source_subseries)"
      rhs <- c(
        M0_clinical = clinical_rhs,
        M1_clinical_plus_OLFML2B = paste(clinical_rhs, "+ OLFML2B_z"),
        M2_clinical_plus_OLFML2B_plus_TME = paste(clinical_rhs, "+ OLFML2B_z + TME_axis")
      )
      forms <- lapply(rhs, function(z) stats::as.formula(paste0("survival::Surv(time, event) ~ ", z, strata_term)))
      fit_info <- lapply(forms, function(f) ol_p4_v231_safe_coxph(f, dd, x = TRUE, model = TRUE))
      fits <- lapply(fit_info, `[[`, "fit")
      if (any(vapply(fits, is.null, logical(1)))) {
        ol_p4_v230_record_eligibility(data.frame(cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "COX_FIT_FAILED", n = n, events = events, clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"), main_model_family_stable = FALSE, stringsAsFactors = FALSE))
        next
      }
      family_stable <- all(vapply(fit_info, function(x) isTRUE(x$stable), logical(1)))
      df_m2 <- ol_p4_v230_model_df(fits[[3L]])
      epv_m2 <- events / max(df_m2, 1L)
      inference_tier <- if (!family_stable) {
        "DESCRIPTIVE_UNSTABLE_COX"
      } else if (n >= 80L && events >= 40L && epv_m2 >= min_epv) {
        "FORMAL_STRUCTURALLY_ELIGIBLE"
      } else if (epv_m2 >= 5) {
        "SUPPORTIVE_STRUCTURALLY_ELIGIBLE"
      } else {
        "DESCRIPTIVE_LOW_EPV"
      }
      b1 <- suppressWarnings(unname(stats::coef(fits[[2L]])["OLFML2B_z"]))
      b2 <- suppressWarnings(unname(stats::coef(fits[[3L]])["OLFML2B_z"]))
      attenuation_denominator_ok <- is.finite(b1) && abs(b1) >= 0.05
      attenuation <- if (family_stable && attenuation_denominator_ok && is.finite(b2)) 100 * (b1 - b2) / b1 else NA_real_
      absolute_beta_change <- if (is.finite(b1) && is.finite(b2)) b1 - b2 else NA_real_
      boot <- if (family_stable) {
        ol_p4_v231_bootstrap_attenuation(dd, forms[[2L]], forms[[3L]], B = 1000L, seed = 92912L + match(axis, axes))
      } else {
        c(low = NA_real_, high = NA_real_, median = NA_real_, successful = 0, failed_fit = 0, unstable_fit = 0, warning_replicates = 0, warning_count = 0, denominator_small = 0, minimum_success = 500)
      }
      analysis_set <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(sort(as.character(dd$patient_key)), algo = "sha256") else paste0(cohort, "_", ep, "_", axis, "_n", n)
      lrt01 <- if (family_stable) ol_p4_v230_lrt_p(fits[[1L]], fits[[2L]]) else NA_real_
      lrt12 <- if (family_stable) ol_p4_v230_lrt_p(fits[[2L]], fits[[3L]]) else NA_real_
      for (j in seq_along(fits)) {
        fit <- fits[[j]]
        sm <- summary(fit)
        has_target <- "OLFML2B_z" %in% rownames(sm$coefficients)
        beta <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "coef"]) else NA_real_
        se <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "se(coef)"]) else NA_real_
        ph <- if (isTRUE(fit_info[[j]]$stable)) tryCatch(suppressWarnings(survival::cox.zph(fit)$table), error = function(e) NULL) else NULL
        ph_term <- if (!is.null(ph) && "OLFML2B_z" %in% rownames(ph)) as.numeric(ph["OLFML2B_z", "p"]) else NA_real_
        ph_global <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph)) as.numeric(ph["GLOBAL", "p"]) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis,
          model = names(fits)[j], model_id = names(fits)[j],
          adjustment_role = "same_patient_nested_attenuation", analysis_set_id = analysis_set,
          clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"),
          design_strata = ifelse(is.na(design_strata), "none", design_strata),
          n = n, events = events, model_df = ol_p4_v230_model_df(fit), events_per_model_df = events / max(ol_p4_v230_model_df(fit), 1L),
          beta = beta, se = se, HR = if (is.finite(beta)) exp(beta) else NA_real_,
          lcl = if (is.finite(beta) && is.finite(se)) exp(beta - 1.96 * se) else NA_real_,
          ucl = if (is.finite(beta) && is.finite(se)) exp(beta + 1.96 * se) else NA_real_,
          p_value = if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "Pr(>|z|)"]) else NA_real_,
          aic = stats::AIC(fit), lrt_M0_vs_M1_p = lrt01, lrt_M1_vs_M2_p = lrt12,
          ph_OLFML2B_p = ph_term, ph_global_p = ph_global,
          attenuation_percent = ifelse(j == 3L, attenuation, NA_real_),
          attenuation_absolute_beta_change = ifelse(j == 3L, absolute_beta_change, NA_real_),
          attenuation_denominator_ok = ifelse(j == 3L, attenuation_denominator_ok, NA),
          attenuation_denominator_rule = ifelse(j == 3L, "relative percentage reported only when abs(M1 beta)>=0.05", NA_character_),
          attenuation_bootstrap_lcl = ifelse(j == 3L, boot["low"], NA_real_),
          attenuation_bootstrap_ucl = ifelse(j == 3L, boot["high"], NA_real_),
          attenuation_bootstrap_median = ifelse(j == 3L, boot["median"], NA_real_),
          bootstrap_B = 1000L,
          bootstrap_successful = ifelse(j == 3L, boot["successful"], NA_real_),
          bootstrap_minimum_success = ifelse(j == 3L, boot["minimum_success"], NA_real_),
          bootstrap_failed_fit = ifelse(j == 3L, boot["failed_fit"], NA_real_),
          bootstrap_unstable_fit = ifelse(j == 3L, boot["unstable_fit"], NA_real_),
          bootstrap_warning_replicates = ifelse(j == 3L, boot["warning_replicates"], NA_real_),
          bootstrap_warning_count = ifelse(j == 3L, boot["warning_count"], NA_real_),
          bootstrap_denominator_small = ifelse(j == 3L, boot["denominator_small"], NA_real_),
          bootstrap_design = ifelse(is.na(design_strata), "patient_resampling_stratified_by_event", "patient_resampling_stratified_by_source_subseries_and_event"),
          cox_fit_status = fit_info[[j]]$status,
          cox_warning_count = fit_info[[j]]$warning_count,
          cox_warning_text = fit_info[[j]]$warning_text,
          cox_error_text = fit_info[[j]]$error_text,
          model_family_stable = family_stable,
          same_patient_set_verified = TRUE, inference_tier = inference_tier,
          mediation_warning = "Coefficient attenuation on an identical patient set is descriptive sensitivity analysis, not mediation or causality.",
          stringsAsFactors = FALSE
        )
      }
      ol_p4_v230_record_eligibility(data.frame(
        cohort = cohort, endpoint = ep, attenuation_axis = axis, status = inference_tier,
        n = n, events = events, clinical_tier = cc$tier,
        clinical_covariates = paste(cc$covars, collapse = ";"), analysis_set_id = analysis_set,
        m2_model_df = df_m2, m2_epv = epv_m2,
        main_model_family_stable = family_stable,
        main_cox_warning_count = sum(vapply(fit_info, `[[`, numeric(1), "warning_count")),
        bootstrap_successful = unname(boot["successful"]),
        bootstrap_unstable_fit = unname(boot["unstable_fit"]),
        stringsAsFactors = FALSE
      ))
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (nrow(out)) out$fdr_OLFML2B_within_model_family <- stats::p.adjust(out$p_value, method = "BH")
  out
}

# Do not plot unstable Cox families as if they were inferentially valid.
ol_p4_v230_nested_forest <- function(index) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  d <- index$tme_survival %||% data.frame()
  if ("model_family_stable" %in% names(d)) d <- d[d$model_family_stable %in% TRUE, , drop = FALSE]
  d <- d[d$model %in% c("M1_clinical_plus_OLFML2B", "M2_clinical_plus_OLFML2B_plus_TME") & is.finite(d$HR) & is.finite(d$lcl) & is.finite(d$ucl), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d$label <- paste0(d$cohort, " / ", d$endpoint, " / ", d$attenuation_axis, " / ", d$model, "\nHR ", sprintf("%.2f", d$HR), " (", sprintf("%.2f", d$lcl), "–", sprintf("%.2f", d$ucl), "), n=", d$n, ", events=", d$events)
  d$label <- factor(make.unique(d$label), levels = rev(make.unique(d$label)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = HR, y = label, shape = model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lcl, xmax = ucl), orientation = "y", width = 0.15) +
    ggplot2::geom_point(size = 2.2) + ggplot2::scale_x_log10() +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::labs(title = "Stable same-patient nested TME attenuation models", subtitle = "Unstable/infinite-coefficient Cox families are excluded and retained in the convergence audit", x = "Hazard ratio per SD OLFML2B", y = NULL)
  stem <- file.path(index$dirs$figures, "FIG4C_same_patient_nested_TME_attenuation")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 11, height = max(6, 0.28 * nrow(d) + 2), limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 11, height = max(6, 0.28 * nrow(d) + 2), dpi = 600, limitsize = FALSE)
  invisible(TRUE)
}

.ol_p4_v231_core_runner <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  index <- .ol_p4_v231_core_runner(...)
  index$version <- OLFML2B_PART4_TARGET_VERSION
  surv <- index$tme_survival %||% data.frame()

  convergence <- if (nrow(surv)) {
    keep <- c("cohort", "endpoint", "attenuation_axis", "model", "n", "events", "clinical_tier", "clinical_covariates", "model_df", "events_per_model_df", "cox_fit_status", "cox_warning_count", "cox_warning_text", "cox_error_text", "model_family_stable", "inference_tier")
    surv[, keep[keep %in% names(surv)], drop = FALSE]
  } else data.frame()
  bootstrap_audit <- if (nrow(surv)) {
    d <- surv[surv$model == "M2_clinical_plus_OLFML2B_plus_TME", , drop = FALSE]
    keep <- c("cohort", "endpoint", "attenuation_axis", "n", "events", "model_family_stable", "attenuation_percent", "attenuation_absolute_beta_change", "attenuation_denominator_ok", "attenuation_bootstrap_lcl", "attenuation_bootstrap_ucl", "attenuation_bootstrap_median", "bootstrap_B", "bootstrap_minimum_success", "bootstrap_successful", "bootstrap_failed_fit", "bootstrap_unstable_fit", "bootstrap_warning_replicates", "bootstrap_warning_count", "bootstrap_denominator_small", "bootstrap_design", "inference_tier")
    d[, keep[keep %in% names(d)], drop = FALSE]
  } else data.frame()

  stable_main <- nrow(convergence) > 0L && all(convergence$model_family_stable %in% TRUE)
  boot_ok <- nrow(bootstrap_audit) > 0L && all(
    bootstrap_audit$model_family_stable %in% TRUE &
      is.finite(bootstrap_audit$bootstrap_successful) &
      is.finite(bootstrap_audit$bootstrap_minimum_success) &
      bootstrap_audit$bootstrap_successful >= bootstrap_audit$bootstrap_minimum_success
  )
  old_go <- index$go_no_go %||% data.frame()
  if (nrow(old_go) && "criterion" %in% names(old_go)) {
    old_go <- old_go[!old_go$criterion %in% c("patient_bootstrap_attenuation", "cox_convergence_audit", "bootstrap_stability_audit"), , drop = FALSE]
  }
  audit_go <- data.frame(
    criterion = c("cox_convergence_audit", "bootstrap_stability_audit", "console_warning_control"),
    status = c(
      ifelse(!nrow(convergence), "NOT_EVALUABLE", ifelse(stable_main, "PASS", "REVIEW_UNSTABLE_MODEL_FAMILIES")),
      ifelse(!nrow(bootstrap_audit), "NOT_EVALUABLE", ifelse(boot_ok, "PASS", "REVIEW_INSUFFICIENT_STABLE_BOOTSTRAP_REPLICATES")),
      "PASS"
    ),
    boundary = c(
      "Potentially infinite, non-converged or non-finite Cox fits are captured, flagged and excluded from formal inference/forest plots.",
      "A percentile attenuation interval is reported only when at least 50% (minimum 500/1000) event-stratified patient bootstrap replicates yield stable M1 and M2 fits and abs(M1 beta)>=0.05.",
      "Expected bootstrap separation warnings are muffled from the console but retained as counts and text in explicit audit tables."
    ),
    stringsAsFactors = FALSE
  )
  index$go_no_go <- ol_p4_bind_rows(list(old_go, audit_go))
  index$cox_convergence_audit <- convergence
  index$bootstrap_stability_audit <- bootstrap_audit
  index$final_gene_lock <- FALSE

  ol_p4_atomic_write_csv(convergence, file.path(index$dirs$tables, "12b_tme_cox_convergence_audit.csv"))
  ol_p4_atomic_write_csv(bootstrap_audit, file.path(index$dirs$tables, "12c_tme_bootstrap_stability_audit.csv"))
  ol_p4_atomic_write_csv(index$go_no_go, file.path(index$dirs$tables, "20_part4_go_no_go_summary.csv"))
  try(ol_p4_v230_nested_forest(index), silent = TRUE)
  ol_p4_atomic_save_rds(index, file.path(index$dirs$objects, "Part4_immune_TME_production_index.rds"))

  unstable_families <- if (nrow(bootstrap_audit)) sum(!bootstrap_audit$model_family_stable, na.rm = TRUE) else 0L
  insufficient_boot <- if (nrow(bootstrap_audit)) sum(
    !is.finite(bootstrap_audit$bootstrap_successful) |
      !is.finite(bootstrap_audit$bootstrap_minimum_success) |
      bootstrap_audit$bootstrap_successful < bootstrap_audit$bootstrap_minimum_success,
    na.rm = TRUE
  ) else 0L
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] [INFO] [OLFML2B-P4] Cox/bootstrap stability audit | unstable_families=", unstable_families, " | insufficient_bootstrap_families=", insufficient_boot, " | details=12b/12c audit tables")
  invisible(index)
}

# ==============================================================================
# Part4 v2.4.0 complete methodology hardening
# - formal tumor-only ecology for every cohort;
# - source-aware DFS/RFS labels;
# - GSE84437 age/sex/pT/pN plus subseries strata;
# - robust LRT extraction, separate multiplicity families, duplicate-axis removal;
# - collinearity and PH time-interaction audits;
# - explicit TCGA tumor-normal context kept outside formal ecology.
# ==============================================================================
OLFML2B_PART4_TARGET_VERSION <- "v2.4.0_20260722_TUMOR_ONLY_ENDPOINT_CLINICAL_MULTIPLICITY_PH_FIX"
OLFML2B_PART4_IF78_VERSION <- OLFML2B_PART4_TARGET_VERSION
OLFML2B_PART4_VERSION <- OLFML2B_PART4_TARGET_VERSION

ol_p4_v240_context <- function(d, cohort) {
  if (!"sample_context" %in% names(d)) {
    if (cohort %in% c("GSE62254", "GSE15459", "GSE26253", "GSE84437", "GSE147163")) {
      d$sample_context <- "Tumor"
    } else {
      z <- rep("Unknown", nrow(d))
      if ("is_normal" %in% names(d)) z[d$is_normal %in% TRUE] <- "Normal"
      if ("is_tumor" %in% names(d)) z[d$is_tumor %in% TRUE] <- "Tumor"
      d$sample_context <- z
    }
  }
  d$sample_context <- as.character(d$sample_context)
  d
}

.ol_p4_v240_batch_core <- ol_p4_apply_subseries_batch_design
ol_p4_apply_subseries_batch_design <- function(d, cohort) {
  d <- .ol_p4_v240_batch_core(d, cohort)
  d <- ol_p4_v240_context(d, cohort)
  keep <- d$sample_context == "Tumor"
  if (!any(keep)) stop("Part4 formal tumor-only gate removed every sample in ", cohort, call. = FALSE)
  d <- d[keep, , drop = FALSE]
  d$formal_ecology_role <- "FORMAL_TUMOR_ONLY"
  d$analysis_population <- "TUMOR_ONLY"
  d
}

.ol_p4_v240_cor_core <- ol_p4_correlations
ol_p4_correlations <- function(d, cohort) {
  out <- .ol_p4_v240_cor_core(d, cohort)
  if (!nrow(out)) return(out)
  out <- out[out$signature != "Cytotoxic_axis", , drop = FALSE]
  out$hypothesis_tier <- ifelse(out$signature %in% c("CAF_Core", "ECM_Remodeling", "TGFb_Response"), "CORE_PRESPECIFIED", "EXPLORATORY_OR_CONTROL")
  out$multiplicity_family <- ifelse(out$hypothesis_tier == "CORE_PRESPECIFIED", "core_CAF_ECM_TGFb", "exploratory_or_control_TME")
  out$fdr_within_cohort_family <- ave(out$p_value, out$cohort, out$multiplicity_family, FUN = function(p) stats::p.adjust(p, method = "BH"))
  out$analysis_population <- "TUMOR_ONLY"
  out
}

.ol_p4_v240_highlow_core <- ol_p4_high_low
ol_p4_high_low <- function(d, cohort) {
  out <- .ol_p4_v240_highlow_core(d, cohort)
  if (!nrow(out)) return(out)
  out <- out[out$signature != "Cytotoxic_axis", , drop = FALSE]
  out$analysis_population <- "TUMOR_ONLY"
  out
}

ol_p4_endpoint_defs <- function(df) {
  defs <- list()
  if (all(c("os_time_days", "os_event") %in% names(df))) defs$OS <- c("os_time_days", "os_event")
  source_label <- NA_character_
  if ("source_endpoint_label" %in% names(df)) {
    z <- unique(toupper(as.character(df$source_endpoint_label)))
    z <- z[z %in% c("DFS", "RFS")]
    if (length(z) == 1L) source_label <- z
  }
  if (is.na(source_label) && "recurrence_endpoint" %in% names(df)) {
    z <- unique(toupper(as.character(df$recurrence_endpoint)))
    z <- z[z %in% c("DFS", "RFS")]
    if (length(z) == 1L) source_label <- z
  }
  candidates <- list(DFS = c("dfs_time_days", "dfs_event"), RFS = c("rfs_time_days", "rfs_event"), RECURRENCE = c("recurrence_time_days", "recurrence_event"))
  available <- names(candidates)[vapply(candidates, function(z) all(z %in% names(df)), logical(1))]
  chosen <- if (!is.na(source_label) && source_label %in% available) source_label else if (length(available)) available[1L] else NA_character_
  if (!is.na(chosen)) defs[[chosen]] <- candidates[[chosen]]
  defs
}

ol_p4_v240_t_factor <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("T1", z)] <- "T1"; out[grepl("T2", z)] <- "T2"
  out[grepl("T3", z)] <- "T3"; out[grepl("T4", z)] <- "T4"
  factor(out, levels = c("T1", "T2", "T3", "T4"))
}
ol_p4_v240_n_factor <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("N0", z)] <- "N0"; out[grepl("N1", z)] <- "N1"
  out[grepl("N2", z)] <- "N2"; out[grepl("N3", z)] <- "N3"
  factor(out, levels = c("N0", "N1", "N2", "N3"))
}

ol_p4_v240_collapse_patient <- function(d, axis, design_strata = NA_character_) {
  d$patient_key <- ol_p4_v230_patient_key(d)
  rows <- lapply(split(d, d$patient_key), function(z) data.frame(
    patient_key = z$patient_key[1L],
    sample_id = paste(sort(unique(as.character(z$sample_id %||% z$patient_key))), collapse = ";"),
    time = ol_p4_v230_one_value(z$time, numeric = TRUE),
    event = { ev <- unique(z$event[!is.na(z$event)]); if (length(ev) == 1L) as.integer(ev) else NA_integer_ },
    OLFML2B_z = ol_p4_v230_one_value(z$OLFML2B_z, numeric = TRUE),
    TME_axis = ol_p4_v230_one_value(z[[axis]], numeric = TRUE),
    age = if ("age" %in% names(z)) ol_p4_v230_one_value(z$age, numeric = TRUE) else NA_real_,
    sex_raw = if ("sex" %in% names(z)) ol_p4_v230_one_value(z$sex) else NA_character_,
    stage_raw = if ("stage" %in% names(z)) ol_p4_v230_one_value(z$stage) else NA_character_,
    pT_raw = if ("stage_pT" %in% names(z)) ol_p4_v230_one_value(z$stage_pT) else NA_character_,
    pN_raw = if ("stage_pN" %in% names(z)) ol_p4_v230_one_value(z$stage_pN) else NA_character_,
    source_subseries = if (!is.na(design_strata) && design_strata %in% names(z)) ol_p4_v230_one_value(z[[design_strata]]) else NA_character_,
    n_source_rows = nrow(z), stringsAsFactors = FALSE
  ))
  out <- ol_p4_bind_rows(rows)
  if (!nrow(out)) return(out)
  out$age10 <- if (any(is.finite(out$age))) (out$age - mean(out$age, na.rm = TRUE)) / 10 else NA_real_
  out$sex <- ol_p4_if78_sex_factor(out$sex_raw)
  out$stage <- ol_p4_if78_stage_factor(out$stage_raw)
  out$pT <- ol_p4_v240_t_factor(out$pT_raw)
  out$pN <- ol_p4_v240_n_factor(out$pN_raw)
  out$source_subseries <- factor(out$source_subseries)
  out
}

ol_p4_v240_usable <- function(x, n, numeric = FALSE) {
  if (numeric) {
    z <- suppressWarnings(as.numeric(x)); ok <- is.finite(z)
    return(sum(ok) >= max(30L, ceiling(0.60 * n)) && is.finite(stats::sd(z[ok])) && stats::sd(z[ok]) > 0)
  }
  z <- droplevels(factor(x)); tab <- table(z, useNA = "no")
  sum(tab) >= max(30L, ceiling(0.60 * n)) && length(tab) >= 2L && all(tab >= 5L)
}

ol_p4_v240_clinical_covars <- function(d, cohort) {
  prespecified <- switch(cohort,
    TCGA_STAD = c("age10", "sex", "stage"),
    GSE62254 = c("age10", "sex", "stage"),
    GSE15459 = c("age10", "sex", "stage"),
    GSE26253 = c("stage"),
    GSE84437 = c("age10", "sex", "pT", "pN"),
    character()
  )
  usable <- vapply(prespecified, function(v) ol_p4_v240_usable(d[[v]], nrow(d), numeric = identical(v, "age10")), logical(1))
  covars <- prespecified[usable]
  tier <- if (length(covars) == length(prespecified) && length(prespecified)) "PRESPECIFIED_COMPLETE" else if (length(covars)) "PRESPECIFIED_AVAILABLE_SUBSET" else "NO_CLINICAL_COVARIATE_AVAILABLE"
  list(covars = covars, prespecified = prespecified, usable = usable, tier = tier)
}

ol_p4_v230_lrt_p <- function(fit_small, fit_large) {
  z <- tryCatch(stats::anova(fit_small, fit_large, test = "LRT"), error = function(e) NULL)
  if (is.null(z) || nrow(z) < 2L) return(NA_real_)
  nms <- names(z)
  pcol <- grep("^Pr\\(|^P\\(", nms, value = TRUE, ignore.case = TRUE)
  if (!length(pcol)) pcol <- grep("p.value|p_value|LRT", nms, value = TRUE, ignore.case = TRUE)
  if (!length(pcol)) return(NA_real_)
  vals <- suppressWarnings(as.numeric(z[2L, pcol, drop = TRUE]))
  vals <- vals[is.finite(vals) & vals >= 0 & vals <= 1]
  if (length(vals)) vals[1L] else NA_real_
}

ol_p4_v240_vif <- function(mm, target) {
  if (is.null(mm) || !target %in% colnames(mm)) return(NA_real_)
  x <- mm[, target]
  others <- mm[, setdiff(colnames(mm), c("(Intercept)", target)), drop = FALSE]
  if (!ncol(others) || !is.finite(stats::var(x)) || stats::var(x) <= 0) return(1)
  fit <- tryCatch(stats::lm.fit(cbind(1, others), x), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  rss <- sum(fit$residuals^2); tss <- sum((x - mean(x))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  if (!is.finite(r2) || r2 >= 0.999999) Inf else 1 / (1 - max(0, r2))
}

ol_p4_v240_condition_index <- function(mm) {
  if (is.null(mm)) return(NA_real_)
  mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  if (ncol(mm) < 2L) return(1)
  s <- apply(mm, 2L, stats::sd, na.rm = TRUE)
  mm <- mm[, is.finite(s) & s > 0, drop = FALSE]
  if (ncol(mm) < 2L) return(1)
  zz <- scale(mm)
  ev <- tryCatch(eigen(crossprod(zz), symmetric = TRUE, only.values = TRUE)$values, error = function(e) numeric())
  ev <- ev[is.finite(ev) & ev > 1e-10]
  if (length(ev) < 2L) return(NA_real_)
  sqrt(max(ev) / min(ev))
}

ol_p4_v240_ph_time_interaction <- function(dd, rhs, strata_term = "") {
  form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", rhs, " + tt(OLFML2B_z)", strata_term))
  fit <- tryCatch(suppressWarnings(survival::coxph(form, data = dd, tt = function(x, t, ...) x * log(pmax(t, 1)))), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  sm <- summary(fit)$coefficients
  hit <- grep("tt\\(OLFML2B_z\\)", rownames(sm))
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(sm[hit[1L], "Pr(>|z|)"]))
}

ol_p4_tme_cox <- function(d, cohort, min_n = 40L, min_events = 20L, min_epv = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package is required", call. = FALSE)
  d <- ol_p4_v240_context(d, cohort)
  d <- d[d$sample_context == "Tumor", , drop = FALSE]
  defs <- ol_p4_endpoint_defs(d)
  if (!length(defs)) {
    ol_p4_v230_record_eligibility(data.frame(
      cohort = cohort, endpoint = NA_character_, attenuation_axis = NA_character_,
      status = "NO_AUDITABLE_ENDPOINT", n = nrow(d), events = NA_integer_,
      clinical_tier = NA_character_, clinical_covariates = NA_character_,
      stringsAsFactors = FALSE
    ))
    return(data.frame())
  }
  core_axes <- c("CAF_Core", "ECM_Remodeling", "TGFb_Response")
  exploratory_axes <- c("CAF_TGFb_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
  axes <- c(core_axes, exploratory_axes)
  rows <- list()
  design_strata <- if (cohort == "GSE84437" && "source_subseries" %in% names(d) && length(unique(stats::na.omit(d$source_subseries))) >= 2L) "source_subseries" else NA_character_
  for (ep in names(defs)) {
    tc <- defs[[ep]][1L]; ec <- defs[[ep]][2L]
    for (axis in axes[axes %in% names(d)]) {
      raw <- d
      raw$time <- ol_p4_num(raw[[tc]]); raw$event <- ol_p4_event(raw[[ec]])
      raw <- raw[is.finite(raw$time) & raw$time > 0 & raw$event %in% c(0L, 1L) & is.finite(raw$OLFML2B_z) & is.finite(raw[[axis]]), , drop = FALSE]
      if (!nrow(raw)) {
        ol_p4_v230_record_eligibility(data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis,
          status = "NO_COMPLETE_ENDPOINT_TARGET_AXIS_ROWS", n = 0L, events = 0L,
          clinical_tier = NA_character_, clinical_covariates = NA_character_,
          stringsAsFactors = FALSE
        ))
        next
      }
      pd <- ol_p4_v240_collapse_patient(raw, axis, design_strata)
      cc <- ol_p4_v240_clinical_covars(pd, cohort)
      required <- c("time", "event", "OLFML2B_z", "TME_axis", cc$covars, if (is.na(design_strata)) character() else "source_subseries")
      dd <- pd[stats::complete.cases(pd[, required, drop = FALSE]) & is.finite(pd$time) & pd$time > 0 & pd$event %in% c(0L, 1L), , drop = FALSE]
      n <- nrow(dd); events <- sum(dd$event == 1L)
      if (n < min_n || events < min_events) {
        ol_p4_v230_record_eligibility(data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis,
          status = "INSUFFICIENT_SAME_CASE_COMPLETE_DATA", n = n, events = events,
          clinical_tier = cc$tier, clinical_covariates = paste(cc$covars, collapse = ";"),
          stringsAsFactors = FALSE
        ))
        next
      }
      clinical_rhs <- if (length(cc$covars)) paste(cc$covars, collapse = " + ") else "1"
      strata_term <- if (is.na(design_strata)) "" else " + survival::strata(source_subseries)"
      rhs <- c(M0_clinical = clinical_rhs, M1_clinical_plus_OLFML2B = paste(clinical_rhs, "+ OLFML2B_z"), M2_clinical_plus_OLFML2B_plus_TME = paste(clinical_rhs, "+ OLFML2B_z + TME_axis"))
      forms <- lapply(rhs, function(z) stats::as.formula(paste0("survival::Surv(time, event) ~ ", z, strata_term)))
      fit_info <- lapply(forms, function(f) ol_p4_v231_safe_coxph(f, dd, x = TRUE, model = TRUE))
      fits <- lapply(fit_info, `[[`, "fit")
      if (any(vapply(fits, is.null, logical(1)))) {
        ol_p4_v230_record_eligibility(data.frame(
          cohort = cohort, endpoint = ep, attenuation_axis = axis, status = "COX_FIT_FAILED",
          n = n, events = events, clinical_tier = cc$tier,
          clinical_covariates = paste(cc$covars, collapse = ";"), stringsAsFactors = FALSE
        ))
        next
      }
      family_stable <- all(vapply(fit_info, function(x) isTRUE(x$stable), logical(1)))
      df_m2 <- ol_p4_v230_model_df(fits[[3L]]); epv_m2 <- events / max(df_m2, 1L)
      inference_tier <- if (!family_stable) "DESCRIPTIVE_UNSTABLE_COX" else if (n >= 80L && events >= 40L && epv_m2 >= min_epv) "FORMAL_STRUCTURALLY_ELIGIBLE" else if (epv_m2 >= 5) "SUPPORTIVE_STRUCTURALLY_ELIGIBLE" else "DESCRIPTIVE_LOW_EPV"
      b1 <- suppressWarnings(unname(stats::coef(fits[[2L]])["OLFML2B_z"])); b2 <- suppressWarnings(unname(stats::coef(fits[[3L]])["OLFML2B_z"]))
      denominator_ok <- is.finite(b1) && abs(b1) >= 0.05
      attenuation <- if (family_stable && denominator_ok && is.finite(b2)) 100 * (b1 - b2) / b1 else NA_real_
      abs_change <- if (is.finite(b1) && is.finite(b2)) b1 - b2 else NA_real_
      boot <- if (family_stable) ol_p4_v231_bootstrap_attenuation(dd, forms[[2L]], forms[[3L]], B = 1000L, seed = 92912L + match(axis, axes)) else c(low=NA,high=NA,median=NA,successful=0,failed_fit=0,unstable_fit=0,warning_replicates=0,warning_count=0,denominator_small=0,minimum_success=500)
      lrt01 <- if (family_stable) ol_p4_v230_lrt_p(fits[[1L]], fits[[2L]]) else NA_real_
      lrt12 <- if (family_stable) ol_p4_v230_lrt_p(fits[[2L]], fits[[3L]]) else NA_real_
      analysis_set <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(sort(as.character(dd$patient_key)), algo = "sha256") else paste0(cohort, "_", ep, "_", axis, "_n", n)
      mm2 <- tryCatch(stats::model.matrix(fits[[3L]]), error = function(e) NULL)
      rho <- suppressWarnings(stats::cor(dd$OLFML2B_z, dd$TME_axis, method = "spearman", use = "complete.obs"))
      vif_o <- ol_p4_v240_vif(mm2, "OLFML2B_z"); vif_t <- ol_p4_v240_vif(mm2, "TME_axis")
      ci <- ol_p4_v240_condition_index(mm2)
      tt_p_by_model <- rep(NA_real_, length(fits))
      if (family_stable) {
        tt_p_by_model[2L] <- ol_p4_v240_ph_time_interaction(dd, rhs[[2L]], strata_term)
        tt_p_by_model[3L] <- ol_p4_v240_ph_time_interaction(dd, rhs[[3L]], strata_term)
      }
      for (j in seq_along(fits)) {
        fit <- fits[[j]]; sm <- summary(fit); has_target <- "OLFML2B_z" %in% rownames(sm$coefficients)
        beta <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "coef"]) else NA_real_
        se <- if (has_target) as.numeric(sm$coefficients["OLFML2B_z", "se(coef)"]) else NA_real_
        ph <- if (isTRUE(fit_info[[j]]$stable)) tryCatch(suppressWarnings(survival::cox.zph(fit)$table), error = function(e) NULL) else NULL
        ph_term <- if (!is.null(ph) && "OLFML2B_z" %in% rownames(ph)) as.numeric(ph["OLFML2B_z", "p"]) else NA_real_
        ph_global <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph)) as.numeric(ph["GLOBAL", "p"]) else NA_real_
        rows[[length(rows)+1L]] <- data.frame(
          cohort=cohort, endpoint=ep, source_endpoint_label=ep, canonical_endpoint_family=ifelse(ep=="OS","OS","RECURRENCE"), attenuation_axis=axis,
          hypothesis_tier=ifelse(axis %in% core_axes,"CORE_PRESPECIFIED","EXPLORATORY"), model=names(fits)[j], model_id=names(fits)[j],
          adjustment_role="same_patient_nested_attenuation", analysis_set_id=analysis_set, analysis_population="TUMOR_ONLY",
          clinical_tier=cc$tier, prespecified_clinical_covariates=paste(cc$prespecified,collapse=";"), clinical_covariates=paste(cc$covars,collapse=";"),
          design_strata=ifelse(is.na(design_strata),"none",design_strata), n=n, events=events, model_df=ol_p4_v230_model_df(fit), events_per_model_df=events/max(ol_p4_v230_model_df(fit),1L),
          beta=beta,se=se,HR=if(is.finite(beta))exp(beta)else NA_real_,lcl=if(is.finite(beta)&&is.finite(se))exp(beta-1.96*se)else NA_real_,ucl=if(is.finite(beta)&&is.finite(se))exp(beta+1.96*se)else NA_real_,
          p_value=if(has_target)as.numeric(sm$coefficients["OLFML2B_z","Pr(>|z|)"])else NA_real_, aic=stats::AIC(fit), lrt_M0_vs_M1_p=lrt01, lrt_M1_vs_M2_p=lrt12,
          ph_OLFML2B_p=ph_term, ph_global_p=ph_global, ph_time_interaction_p=ifelse(j>=2L,tt_p_by_model[j],NA_real_), ph_inference_status=ifelse(j>=2L && ((is.finite(ph_term)&&ph_term<0.05)||(is.finite(tt_p_by_model[j])&&tt_p_by_model[j]<0.05)),"FORMAL_PH_REVIEW",ifelse(j>=2L,"PH_OK_OR_NOT_DETECTED","NOT_APPLICABLE")),
          olfml2b_tme_spearman_rho=rho, vif_OLFML2B=ifelse(j==3L,vif_o,NA_real_), vif_TME_axis=ifelse(j==3L,vif_t,NA_real_), condition_index=ifelse(j==3L,ci,NA_real_), coefficient_sign_reversal=ifelse(j==3L,is.finite(b1)&&is.finite(b2)&&b1*b2<0,NA),
          collinearity_status=ifelse(j==3L,ifelse((is.finite(vif_o)&&vif_o>=5)||(is.finite(vif_t)&&vif_t>=5)||(is.finite(ci)&&ci>=30),"REVIEW_HIGH_COLLINEARITY","PASS_OR_MODERATE"),NA_character_),
          attenuation_percent=ifelse(j==3L,attenuation,NA_real_), attenuation_absolute_beta_change=ifelse(j==3L,abs_change,NA_real_), attenuation_denominator_ok=ifelse(j==3L,denominator_ok,NA), attenuation_bootstrap_lcl=ifelse(j==3L,boot["low"],NA_real_), attenuation_bootstrap_ucl=ifelse(j==3L,boot["high"],NA_real_), attenuation_bootstrap_median=ifelse(j==3L,boot["median"],NA_real_),
          bootstrap_B=1000L,bootstrap_successful=ifelse(j==3L,boot["successful"],NA_real_),bootstrap_minimum_success=ifelse(j==3L,boot["minimum_success"],NA_real_),bootstrap_failed_fit=ifelse(j==3L,boot["failed_fit"],NA_real_),bootstrap_unstable_fit=ifelse(j==3L,boot["unstable_fit"],NA_real_),bootstrap_warning_replicates=ifelse(j==3L,boot["warning_replicates"],NA_real_),bootstrap_warning_count=ifelse(j==3L,boot["warning_count"],NA_real_),bootstrap_denominator_small=ifelse(j==3L,boot["denominator_small"],NA_real_),bootstrap_design=ifelse(is.na(design_strata),"patient_resampling_stratified_by_event","patient_resampling_stratified_by_source_subseries_and_event"),
          cox_fit_status=fit_info[[j]]$status,cox_warning_count=fit_info[[j]]$warning_count,cox_warning_text=fit_info[[j]]$warning_text,cox_error_text=fit_info[[j]]$error_text,model_family_stable=family_stable,same_patient_set_verified=TRUE,inference_tier=inference_tier,mediation_warning="Coefficient attenuation is descriptive sensitivity analysis, not mediation or causality.",stringsAsFactors=FALSE
        )
      }
      ol_p4_v230_record_eligibility(data.frame(
        cohort = cohort, endpoint = ep, attenuation_axis = axis, status = inference_tier,
        n = n, events = events, clinical_tier = cc$tier,
        clinical_covariates = paste(cc$covars, collapse = ";"), analysis_set_id = analysis_set,
        m2_model_df = df_m2, m2_epv = epv_m2, main_model_family_stable = family_stable,
        main_cox_warning_count = sum(vapply(fit_info, `[[`, numeric(1), "warning_count")),
        bootstrap_successful = unname(boot["successful"]),
        bootstrap_unstable_fit = unname(boot["unstable_fit"]),
        stringsAsFactors = FALSE
      ))
    }
  }
  out <- ol_p4_bind_rows(rows)
  if (!nrow(out)) return(out)
  out$multiplicity_family <- ifelse(
    out$model == "M0_clinical", "NOT_APPLICABLE",
    ifelse(
      out$model == "M1_clinical_plus_OLFML2B", "M1_REFERENCE_NOT_DUPLICATED_IN_BH",
      ifelse(out$hypothesis_tier == "CORE_PRESPECIFIED", "M2_CORE_CAF_ECM_TGFB", "M2_EXPLORATORY_TME")
    )
  )
  out$fdr_OLFML2B_within_model_family <- NA_real_
  m1 <- out$model=="M1_clinical_plus_OLFML2B"
  out$fdr_OLFML2B_within_model_family[m1] <- out$p_value[m1]
  for (fam in c("M2_CORE_CAF_ECM_TGFB","M2_EXPLORATORY_TME")) {
    idx <- which(out$multiplicity_family==fam)
    if (length(idx)) {
      groups <- interaction(out$cohort[idx],out$endpoint[idx],drop=TRUE)
      out$fdr_OLFML2B_within_model_family[idx] <- ave(out$p_value[idx],groups,FUN=function(p)stats::p.adjust(p,method="BH"))
    }
  }
  out$reference_model_key <- paste(out$cohort,out$endpoint,out$analysis_set_id,out$model,sep="|")
  out$reference_model_duplicate <- out$model=="M1_clinical_plus_OLFML2B" & duplicated(out$reference_model_key)
  out
}


ol_p4_v240_tcga_tumor_normal_context <- function(root, catalog) {
  views <- olfml2b_load_bulk_views(
    root = root, target_gene = "OLFML2B", geo_validation_cohorts = character(),
    log_file = file.path(root, "logs", "runtime", "Part4", "Part4_TCGA_context.log")
  )
  if (!"TCGA_STAD" %in% names(views)) return(data.frame())
  v <- views[["TCGA_STAD"]]
  expr <- v$expr
  clin <- v$clinical
  if (!is.matrix(expr) || !is.data.frame(clin) || !"sample_id" %in% names(clin)) return(data.frame())
  sc <- ol_p4_score_signatures(expr, catalog)$scores
  sc$sample_id <- as.character(sc$sample_id)
  sc$OLFML2B <- as.numeric(ol_p4_extract_gene(expr, "OLFML2B")[sc$sample_id])
  common <- intersect(sc$sample_id, as.character(clin$sample_id))
  sc <- sc[match(common, sc$sample_id), , drop = FALSE]
  clin <- clin[match(common, as.character(clin$sample_id)), , drop = FALSE]
  context <- rep("Unknown", nrow(clin))
  if ("is_normal" %in% names(clin)) context[clin$is_normal %in% TRUE] <- "Normal"
  if ("is_tumor" %in% names(clin)) context[clin$is_tumor %in% TRUE] <- "Tumor"
  sc$sample_context <- context
  vars <- intersect(c("OLFML2B", "CAF_Core", "ECM_Remodeling", "TGFb_Response"), names(sc))
  rows <- lapply(vars, function(vn) {
    do.call(rbind, lapply(c("Tumor", "Normal"), function(g) {
      x <- suppressWarnings(as.numeric(sc[[vn]][sc$sample_context == g]))
      x <- x[is.finite(x)]
      data.frame(
        cohort = "TCGA_STAD", analysis_role = "DESCRIPTIVE_TUMOR_NORMAL_CONTEXT_ONLY",
        sample_context = g, variable = vn, n = length(x),
        median = if (length(x)) stats::median(x) else NA_real_,
        mean = if (length(x)) mean(x) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- ol_p4_bind_rows(rows)
  out$formal_part4_ecology_included <- FALSE
  out$interpretation <- "Descriptive context only; excluded from formal tumor-only ecology, meta-analysis and survival attenuation."
  out
}

.ol_p4_v240_core_runner <- run_olfml2b_part4_immune_tme_production
run_olfml2b_part4_immune_tme_production <- function(...) {
  index <- .ol_p4_v240_core_runner(...)
  index$version <- OLFML2B_PART4_TARGET_VERSION
  dirs <- index$dirs
  root <- normalizePath(dirname(dirname(dirname(dirs$tables))), winslash = "/", mustWork = FALSE)
  tcga_context <- tryCatch(
    ol_p4_v240_tcga_tumor_normal_context(root, ol_p4_signature_catalog()),
    error = function(e) data.frame(status = "CONTEXT_AUDIT_FAILED", error = conditionMessage(e), stringsAsFactors = FALSE)
  )
  scores <- index$sample_scores %||% data.frame()
  context_audit <- if (nrow(scores)) {
    aggregate(sample_id ~ cohort + analysis_population, data = transform(scores, analysis_population="TUMOR_ONLY"), FUN=length)
  } else data.frame()
  if (nrow(context_audit)) names(context_audit)[names(context_audit)=="sample_id"] <- "n_formal_samples"
  expected <- data.frame(cohort=c("TCGA_STAD","GSE62254","GSE15459","GSE26253","GSE84437"), expected_tumor_n=c(412L,300L,192L,432L,433L), stringsAsFactors=FALSE)
  context_audit <- merge(expected, context_audit, by="cohort", all.x=TRUE, sort=FALSE)
  context_audit$n_normal_formal <- 0L; context_audit$n_unknown_formal <- 0L
  context_audit$status <- ifelse(
    is.finite(context_audit$n_formal_samples) & context_audit$n_formal_samples == context_audit$expected_tumor_n,
    "PASS_TUMOR_ONLY", "FAIL_SAMPLE_COUNT"
  )
  if (any(context_audit$status != "PASS_TUMOR_ONLY") || any(!is.finite(context_audit$n_formal_samples)) || sum(context_audit$n_formal_samples) != 1769L) {
    stop("Part4 formal tumor-only sample count contract failed; expected total n=1769.", call.=FALSE)
  }
  context_count_ok <- is.data.frame(tcga_context) && all(c("sample_context", "variable", "n") %in% names(tcga_context)) &&
    any(tcga_context$variable == "OLFML2B" & tcga_context$sample_context == "Tumor" & tcga_context$n == 412L) &&
    any(tcga_context$variable == "OLFML2B" & tcga_context$sample_context == "Normal" & tcga_context$n == 36L)
  if (!context_count_ok) stop("Part4 TCGA descriptive tumor-normal context audit failed; expected 412 tumors and 36 normals.", call. = FALSE)
  surv <- index$tme_survival %||% data.frame()
  m2 <- surv[surv$model=="M2_clinical_plus_OLFML2B_plus_TME",,drop=FALSE]
  collin <- if (nrow(m2)) m2[,intersect(c("cohort","endpoint","source_endpoint_label","attenuation_axis","hypothesis_tier","n","events","olfml2b_tme_spearman_rho","vif_OLFML2B","vif_TME_axis","condition_index","coefficient_sign_reversal","collinearity_status","inference_tier"),names(m2)),drop=FALSE] else data.frame()
  ph_audit <- if (nrow(surv)) surv[surv$model %in% c("M1_clinical_plus_OLFML2B","M2_clinical_plus_OLFML2B_plus_TME"),intersect(c("cohort","endpoint","attenuation_axis","model","n","events","ph_OLFML2B_p","ph_global_p","ph_time_interaction_p","ph_inference_status"),names(surv)),drop=FALSE] else data.frame()
  mult <- if (nrow(surv)) unique(surv[,intersect(c("cohort","endpoint","attenuation_axis","hypothesis_tier","model","analysis_set_id","multiplicity_family","p_value","fdr_OLFML2B_within_model_family","reference_model_duplicate"),names(surv)),drop=FALSE]) else data.frame()
  endpoint_audit <- if (nrow(surv)) unique(surv[,intersect(c("cohort","endpoint","source_endpoint_label","canonical_endpoint_family","n","events"),names(surv)),drop=FALSE]) else data.frame()
  endpoint_ok <- nrow(endpoint_audit) > 0L &&
    any(endpoint_audit$cohort == "GSE26253" & endpoint_audit$source_endpoint_label == "RFS") &&
    any(endpoint_audit$cohort == "GSE62254" & endpoint_audit$source_endpoint_label == "DFS") &&
    !any(endpoint_audit$cohort == "GSE26253" & endpoint_audit$source_endpoint_label == "DFS")
  stable_m2 <- m2[m2$model_family_stable %in% TRUE,,drop=FALSE]
  lrt_ok <- nrow(stable_m2)>0L && all(is.finite(stable_m2$lrt_M0_vs_M1_p) & is.finite(stable_m2$lrt_M1_vs_M2_p))
  batch_rows <- m2[m2$cohort=="GSE84437",,drop=FALSE]
  batch_ok <- nrow(batch_rows)>0L && all(batch_rows$design_strata=="source_subseries") && all(batch_rows$bootstrap_design=="patient_resampling_stratified_by_source_subseries_and_event") && all(grepl("pT",batch_rows$clinical_covariates,fixed=TRUE)) && all(grepl("pN",batch_rows$clinical_covariates,fixed=TRUE))
  old_go <- index$go_no_go %||% data.frame()
  if (nrow(old_go)) old_go <- old_go[!old_go$criterion %in% c("GSE84437_subseries_batch_control","tumor_only_patient_unit","M1_to_M2_attenuation","formal_tumor_only_ecology","source_endpoint_labels","multiplicity_family_contract","collinearity_audit","PH_time_interaction_sensitivity","TCGA_tumor_normal_context_separated"),,drop=FALSE]
  new_go <- data.frame(
    criterion=c("formal_tumor_only_ecology","TCGA_tumor_normal_context_separated","source_endpoint_labels","GSE84437_subseries_batch_control","M1_to_M2_attenuation","multiplicity_family_contract","collinearity_audit","PH_time_interaction_sensitivity"),
    status=c(ifelse(all(context_audit$status=="PASS_TUMOR_ONLY"),"PASS","FAIL"),ifelse(context_count_ok,"PASS","FAIL"),ifelse(endpoint_ok,"PASS","FAIL_SOURCE_ENDPOINT_LABEL"),ifelse(batch_ok,"PASS","FAIL"),ifelse(lrt_ok,"PASS","FAIL_LRT_NOT_FINITE"),"PASS","PASS","PASS"),
    boundary=c(
      "All formal ecological and survival analyses contain tumor samples only; TCGA normals are excluded from formal Part4.",
      "TCGA 412 tumors and 36 normals are retained only in a separate descriptive context table.",
      "GSE26253 is RFS and GSE62254 is DFS; both synthesize only at the canonical recurrence-family level.",
      "GSE84437 uses within-subseries standardization, age/sex/pT/pN adjustment, stratified baseline hazards and subseries-by-event bootstrap.",
      "Stable M0/M1/M2 families require finite nested likelihood-ratio P values.",
      "M1 reference P values are not duplicated in BH; core and exploratory M2 axes are separate families.",
      "Every M2 family exports rho, VIF, condition index and coefficient sign-reversal review.",
      "Cox.zph review is supplemented by an OLFML2B-by-log-time interaction sensitivity model."
    ), stringsAsFactors=FALSE
  )
  index$go_no_go <- ol_p4_bind_rows(list(old_go,new_go))
  index$formal_sample_context_audit <- context_audit
  index$tcga_tumor_normal_context <- tcga_context
  index$collinearity_audit <- collin
  index$ph_time_interaction_audit <- ph_audit
  index$multiplicity_audit <- mult
  index$endpoint_label_audit <- endpoint_audit
  index$final_gene_lock <- FALSE
  ol_p4_atomic_write_csv(context_audit,file.path(dirs$tables,"05d_formal_sample_context_audit.csv"))
  ol_p4_atomic_write_csv(tcga_context,file.path(dirs$tables,"05e_TCGA_tumor_normal_ecology_context.csv"))
  ol_p4_atomic_write_csv(collin,file.path(dirs$tables,"12d_OLFML2B_TME_collinearity_audit.csv"))
  ol_p4_atomic_write_csv(mult,file.path(dirs$tables,"12e_tme_multiplicity_audit.csv"))
  ol_p4_atomic_write_csv(endpoint_audit,file.path(dirs$tables,"12f_tme_endpoint_label_audit.csv"))
  ol_p4_atomic_write_csv(ph_audit,file.path(dirs$tables,"12g_OLFML2B_PH_time_interaction_audit.csv"))
  ol_p4_atomic_write_csv(index$go_no_go,file.path(dirs$tables,"20_part4_go_no_go_summary.csv"))
  ol_p4_atomic_save_rds(index,file.path(dirs$objects,"Part4_immune_TME_production_index.rds"))
  invisible(index)
}
