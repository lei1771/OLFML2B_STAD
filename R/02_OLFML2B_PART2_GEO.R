# ==============================================================================
# File 02/06: GEO helpers and Part2 bulk cohort acquisition/processing.
# Patched: v2.3.5_20260711_VISIBLE_PATCH_AND_FINAL_ENDPOINT_GATE
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
OLFML2B_PART2_PATCH_VERSION <- "v1.0.0_20260721_OLFML2B_ROW_MATCHED_TNM_CONTRACT"
OLFML2B_PART2_CLINICAL_CONTRACT_VERSION <- "v1.0.0_20260721_OLFML2B_STRUCTURAL_MISSINGNESS_CONTRACT"
if (!exists("run_part0", mode = "function")) {
  sys.source(file.path(.compact_code_root, "00_OLFML2B_PART0_CONFIG_CORE.R"), envir = environment(), chdir = FALSE)
}

# ==============================================================================
# GEO acquisition, platform annotation, clinical harmonisation and audit helpers.
# Official GEO Series Matrix / SOFT / supplementary files are downloaded by code.
# No local input file is expected from the user.
# ==============================================================================

olfml2b_geo_cache_paths <- function(accession, dirs) {
  list(
    cohort_dir = olfml2b_safe_dir_create(file.path(dirs$raw_geo, accession)),
    matrix_gz = file.path(dirs$raw_geo, accession, paste0(accession, "_series_matrix.txt.gz")),
    eset_rds = file.path(dirs$cache_geo, paste0(accession, "_ExpressionSet.rds")),
    family_rds = file.path(dirs$cache_geo, paste0(accession, "_family.rds")),
    gpl_dir = olfml2b_safe_dir_create(file.path(dirs$raw_geo, "GPL")),
    suppl_dir = olfml2b_safe_dir_create(file.path(dirs$raw_geo, accession, "supplementary"))
  )
}

olfml2b_geo_download_series_matrix <- function(accession, dirs, overwrite = FALSE, log_file = NULL) {
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  urls <- olfml2b_geo_urls(accession)
  olfml2b_download(
    urls = c(urls$matrix,
             sub("https://ftp.ncbi.nlm.nih.gov", "https://www.ncbi.nlm.nih.gov", urls$matrix, fixed = TRUE)),
    destfile = paths$matrix_gz,
    overwrite = overwrite,
    min_bytes = 1000L,
    log_file = log_file
  )
}

olfml2b_geo_pick_eset <- function(obj, accession, expected_platform = NULL) {
  if (inherits(obj, "ExpressionSet")) return(obj)
  if (is.list(obj) && length(obj)) {
    candidates <- obj[vapply(obj, inherits, logical(1), what = "ExpressionSet")]
    if (!length(candidates)) olfml2b_abort(accession, ": GEOquery returned no ExpressionSet")
    if (!is.null(expected_platform)) {
      platform <- vapply(candidates, function(es) {
        ann <- tryCatch(Biobase::annotation(es), error = function(e) "")
        if (!nzchar(ann)) {
          pd <- Biobase::pData(es)
          col <- grep("platform", names(pd), ignore.case = TRUE, value = TRUE)[1]
          if (!is.na(col)) as.character(pd[[col]][1]) else ""
        } else ann
      }, character(1))
      hit <- which(toupper(platform) == toupper(expected_platform))
      if (length(hit)) return(candidates[[hit[1]]])
    }
    sizes <- vapply(candidates, function(es) ncol(Biobase::exprs(es)), integer(1))
    return(candidates[[which.max(sizes)]])
  }
  olfml2b_abort(accession, ": unsupported GEO object class: ", paste(class(obj), collapse = ";"))
}

olfml2b_geo_load_eset <- function(accession, dirs, expected_platform = NULL,
                                overwrite = FALSE, log_file = NULL) {
  olfml2b_require_packages(c("GEOquery", "Biobase"), log_file)
  paths <- olfml2b_geo_cache_paths(accession, dirs)

  if (file.exists(paths$eset_rds) && !overwrite) {
    sz <- round(file.info(paths$eset_rds)$size / 1024^2, 2)
    olfml2b_log("INFO", accession, ": ExpressionSet cache exists; skipping Series Matrix parsing: ",
              paths$eset_rds, " (", sz, " MB)", log_file = log_file)
    return(readRDS(paths$eset_rds))
  }

  matrix_file <- tryCatch(
    olfml2b_geo_download_series_matrix(accession, dirs, overwrite, log_file),
    error = function(e) {
      olfml2b_log("WARN", accession, ": direct matrix download/local check failed: ", conditionMessage(e), log_file = log_file)
      NA_character_
    }
  )

  obj <- NULL
  if (!is.na(matrix_file) && file.exists(matrix_file)) {
    sz <- round(file.info(matrix_file)$size / 1024^2, 2)
    olfml2b_log(
      "INFO",
      accession, ": local Series Matrix is available; now parsing locally with GEOquery::getGEO(filename=...). ",
      "This is not a download. Large cohorts such as GSE62254 can take 10-30 minutes on first parse. file=",
      matrix_file, " (", sz, " MB)",
      log_file = log_file
    )
    t0 <- Sys.time()
    obj <- tryCatch(
      GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, AnnotGPL = FALSE),
      error = function(e) {
        olfml2b_log("WARN", accession, ": local series-matrix parsing failed: ", conditionMessage(e), log_file = log_file)
        NULL
      }
    )
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    olfml2b_log("INFO", accession, ": local Series Matrix parse finished in ", elapsed, " min", log_file = log_file)
  }

  if (is.null(obj)) {
    olfml2b_log("INFO", accession, ": using GEOquery accession fallback because local parsing was unavailable/failed", log_file = log_file)
    t0 <- Sys.time()
    obj <- GEOquery::getGEO(
      accession, GSEMatrix = TRUE, AnnotGPL = FALSE,
      destdir = paths$cohort_dir, getGPL = FALSE
    )
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    olfml2b_log("INFO", accession, ": GEOquery accession fallback finished in ", elapsed, " min", log_file = log_file)
  }

  olfml2b_log("INFO", accession, ": selecting ExpressionSet and writing cache: ", paths$eset_rds, log_file = log_file)
  eset <- olfml2b_geo_pick_eset(obj, accession, expected_platform)
  olfml2b_atomic_save_rds(eset, paths$eset_rds)
  olfml2b_log("INFO", accession, ": ExpressionSet cached. samples=", ncol(Biobase::exprs(eset)),
            " probes=", nrow(Biobase::exprs(eset)), log_file = log_file)
  eset
}

olfml2b_geo_get_family <- function(accession, dirs, overwrite = FALSE, log_file = NULL) {
  olfml2b_require_packages("GEOquery", log_file)
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  if (file.exists(paths$family_rds) && !overwrite) return(readRDS(paths$family_rds))
  obj <- tryCatch(
    GEOquery::getGEO(accession, GSEMatrix = FALSE, destdir = paths$cohort_dir, getGPL = FALSE),
    error = function(e) {
      olfml2b_log("WARN", accession, ": family SOFT retrieval failed: ", conditionMessage(e), log_file = log_file)
      NULL
    }
  )
  if (!is.null(obj)) olfml2b_atomic_save_rds(obj, paths$family_rds)
  obj
}

olfml2b_parse_characteristic_cell <- function(value) {
  value <- olfml2b_clean_text(value)
  if (!nzchar(value)) return(NULL)
  # GEO characteristic labels may contain '=' inside the label itself, e.g.
  # 'status (0=non-recurrence, 1=recurrence): 1'.  Therefore the field/value
  # separator must be parsed by colon first; treating '=' as an equal-priority
  # separator truncates the label to 'status (0' and silently loses RFS events.
  pos <- regexpr("\\s*:\\s*", value, perl = TRUE)
  if (pos[1] > 0L) {
    width <- attr(pos, "match.length")[1]
    key <- substr(value, 1L, pos[1] - 1L)
    val <- substr(value, pos[1] + width, nchar(value))
  } else {
    pos <- regexpr("\\s*=\\s*", value, perl = TRUE)
    if (pos[1] <= 0L) return(NULL)
    width <- attr(pos, "match.length")[1]
    key <- substr(value, 1L, pos[1] - 1L)
    val <- substr(value, pos[1] + width, nchar(value))
  }
  key <- olfml2b_clean_names(key)
  val <- olfml2b_clean_text(val)
  if (!nzchar(key)) return(NULL)
  list(key = key, value = val)
}

olfml2b_expand_characteristics <- function(pdata) {
  x <- as.data.frame(pdata, stringsAsFactors = FALSE, check.names = FALSE)
  original_names <- names(x)
  names(x) <- make.unique(olfml2b_clean_names(names(x)), sep = "_")
  char_cols <- grep("characteristics", names(x), ignore.case = TRUE, value = TRUE)
  extracted <- list()
  add_value <- function(key, i, value) {
    if (is.null(extracted[[key]])) extracted[[key]] <<- rep(NA_character_, nrow(x))
    current <- extracted[[key]][i]
    if (is.na(current) || !nzchar(olfml2b_clean_text(current))) extracted[[key]][i] <<- value
  }
  if (length(char_cols)) {
    for (col in char_cols) {
      values <- olfml2b_clean_text(x[[col]])
      for (i in seq_along(values)) {
        val <- values[[i]]
        if (!nzchar(val)) next
        parsed <- olfml2b_parse_characteristic_cell(val)
        if (!is.null(parsed)) add_value(parsed$key, i, parsed$value)
        # Some matrices concatenate multiple key:value records with semicolons.
        # Split only when more than one colon-delimited record is apparent.
        if (grepl(";", val, fixed = TRUE) && length(gregexpr(":", val, fixed = TRUE)[[1]]) > 1L) {
          pieces <- strsplit(val, ";", fixed = TRUE)[[1]]
          for (piece in pieces) {
            parsed_piece <- olfml2b_parse_characteristic_cell(piece)
            if (!is.null(parsed_piece)) add_value(parsed_piece$key, i, parsed_piece$value)
          }
        }
      }
    }
  }
  if (length(extracted)) {
    extra_df <- as.data.frame(extracted, stringsAsFactors = FALSE, check.names = FALSE)
    duplicates <- intersect(names(extra_df), names(x))
    if (length(duplicates)) {
      names(extra_df)[match(duplicates, names(extra_df))] <- paste0(duplicates, "_characteristic")
    }
    x <- cbind(x, extra_df)
  }
  attr(x, "original_names") <- original_names
  x
}

olfml2b_geo_platform_from_eset <- function(eset) {
  ann <- tryCatch(Biobase::annotation(eset), error = function(e) "")
  if (nzchar(ann) && grepl("GPL[0-9]+", ann, ignore.case = TRUE)) {
    return(toupper(regmatches(ann, regexpr("GPL[0-9]+", ann, ignore.case = TRUE))))
  }
  pd <- Biobase::pData(eset)
  vals <- unlist(pd[, grep("platform", names(pd), ignore.case = TRUE), drop = FALSE], use.names = FALSE)
  hit <- regmatches(vals, regexpr("GPL[0-9]+", vals, ignore.case = TRUE))
  hit <- hit[nzchar(hit)]
  if (length(hit)) toupper(hit[1]) else NA_character_
}

olfml2b_annotation_db_mapping <- function(probes, pkg, log_file = NULL) {
  olfml2b_require_packages(c("AnnotationDbi", pkg), log_file)
  db <- get(pkg, envir = asNamespace(pkg))
  keytypes <- AnnotationDbi::keytypes(db)
  keytype <- if ("PROBEID" %in% keytypes) "PROBEID" else if ("PROBE" %in% keytypes) "PROBE" else keytypes[1]
  cols <- intersect(c("SYMBOL", "ENTREZID", "GENENAME"), AnnotationDbi::columns(db))
  tab <- AnnotationDbi::select(db, keys = unique(probes), keytype = keytype, columns = cols)
  names(tab)[names(tab) == keytype] <- "probe_id"
  names(tab) <- olfml2b_clean_names(names(tab))
  symbol_col <- olfml2b_case_when_column(tab, c("symbol", "gene_symbol"))
  entrez_col <- olfml2b_case_when_column(tab, c("entrezid", "entrez_gene_id"))
  out <- data.frame(
    probe_id = as.character(tab$probe_id),
    symbol = if (!is.na(symbol_col)) as.character(tab[[symbol_col]]) else NA_character_,
    entrez_id = if (!is.na(entrez_col)) as.character(tab[[entrez_col]]) else NA_character_,
    mapping_source = paste0("Bioconductor_", pkg),
    stringsAsFactors = FALSE
  )
  out$symbol <- toupper(trimws(out$symbol))
  out <- out[!is.na(out$probe_id) & nzchar(out$probe_id), , drop = FALSE]
  out
}

olfml2b_gpl_candidate_symbol_column <- function(tab) {
  nms <- names(tab)
  clean <- olfml2b_clean_names(nms)
  preferred <- c("gene_symbol", "symbol", "gene_assignment", "gene_symbol_from_gene", "gene")
  for (p in preferred) {
    idx <- which(clean == p)
    if (length(idx)) return(nms[idx[1]])
  }
  idx <- grep("symbol|gene_assignment|gene symbol", nms, ignore.case = TRUE)
  if (length(idx)) nms[idx[1]] else NA_character_
}

olfml2b_gpl_candidate_entrez_column <- function(tab) {
  nms <- names(tab)
  idx <- grep("entrez|gene_id|gene id", nms, ignore.case = TRUE)
  if (length(idx)) nms[idx[1]] else NA_character_
}

olfml2b_parse_gpl_symbol <- function(x) {
  x <- olfml2b_clean_text(x)
  # Affymetrix / Illumina annotation cells may contain 'GENE /// description',
  # 'symbol // entrez', or multiple assignments. Keep all plausible HGNC symbols;
  # ambiguous probes are later excluded from the target-gene primary analysis.
  x <- gsub("///.*$", "", x)
  x <- gsub("//.*$", "", x)
  x <- gsub(";.*$", "", x)
  x <- trimws(x)
  x[!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", x)] <- NA_character_
  toupper(x)
}

olfml2b_download_gpl_table <- function(platform, dirs, overwrite = FALSE, log_file = NULL) {
  olfml2b_require_packages("GEOquery", log_file)
  rds <- file.path(dirs$cache_geo, paste0(platform, "_AnnotGPL.rds"))
  if (file.exists(rds) && !overwrite) return(readRDS(rds))
  obj <- GEOquery::getGEO(platform, AnnotGPL = TRUE, destdir = file.path(dirs$raw_geo, "GPL"))
  tab <- GEOquery::Table(obj)
  olfml2b_atomic_save_rds(tab, rds)
  tab
}

olfml2b_target_entrez_id <- function() {
  ann <- tryCatch(olfml2b_gene_annotation(), error = function(e) NULL)
  target_entrez <- if (is.list(ann) && !is.null(ann$target_entrez)) {
    as.character(ann$target_entrez[1])
  } else {
    "25903"
  }
  if (!identical(target_entrez, "25903")) {
    olfml2b_abort(
      "OLFML2B Entrez guard failed. Expected 25903 for OLFML2B, observed ",
      target_entrez,
      ". This prevents accidental carry-over from unrelated projects."
    )
  }
  target_entrez
}

olfml2b_assert_target_mapping_clean <- function(mapping_audits, target_gene = "OLFML2B") {
  if (is.null(mapping_audits) || !is.data.frame(mapping_audits) || !nrow(mapping_audits)) {
    return(invisible(TRUE))
  }
  target_gene <- toupper(target_gene)
  if ("symbol" %in% names(mapping_audits)) {
    symbols <- toupper(trimws(as.character(mapping_audits$symbol)))
    bad_symbol <- !is.na(symbols) & nzchar(symbols) & symbols != target_gene
    if (any(bad_symbol)) {
      bad <- unique(symbols[bad_symbol])
      olfml2b_abort(
        "Part2 target mapping audit contains non-OLFML2B symbols: ",
        paste(bad, collapse = ";"),
        ". Stop and inspect platform mapping before survival analysis."
      )
    }
  }
  if ("entrez_id" %in% names(mapping_audits)) {
    entrez <- trimws(as.character(mapping_audits$entrez_id))
    entrez <- entrez[!is.na(entrez) & nzchar(entrez)]
    bad_entrez <- setdiff(unique(entrez), olfml2b_target_entrez_id())
    if (length(bad_entrez)) {
      olfml2b_abort(
        "Part2 target mapping audit contains non-OLFML2B Entrez IDs: ",
        paste(bad_entrez, collapse = ";"),
        ". Expected only Entrez 25903 for OLFML2B."
      )
    }
  }
  invisible(TRUE)
}

olfml2b_target_probe_audit <- function(unique_pairs, probes, target_gene = "OLFML2B",
                                     target_entrez = "25903", ambiguous = character()) {
  if (!identical(toupper(target_gene), "OLFML2B") || !identical(as.character(target_entrez), olfml2b_target_entrez_id())) {
    olfml2b_abort("Target-gene guard failed in olfml2b_target_probe_audit: target_gene=", target_gene,
                "; target_entrez=", target_entrez, ". Expected exact OLFML2B / Entrez 25903 only; OLFML2A/OLFM2 and other confusable paralogs must not enter the target audit.")
  }
  symbol <- toupper(trimws(as.character(unique_pairs$symbol)))
  entrez <- trimws(as.character(unique_pairs$entrez_id))
  symbol[is.na(unique_pairs$symbol)] <- NA_character_
  entrez[is.na(unique_pairs$entrez_id)] <- NA_character_
  support <- (!is.na(symbol) & symbol == toupper(target_gene)) |
    (!is.na(entrez) & entrez == target_entrez)
  evidence <- unique_pairs[which(support), , drop = FALSE]
  if (!nrow(evidence)) {
    empty <- data.frame(
      probe_id = character(), symbol = character(), entrez_id = character(),
      mapping_source = character(), symbol_annotations = character(),
      entrez_annotations = character(), target_match_reason = character(),
      ambiguous_symbol = logical(), in_expression = logical(),
      primary_eligible = logical(), stringsAsFactors = FALSE
    )
    return(list(summary = empty, evidence = evidence))
  }
  probe_ids <- unique(as.character(evidence$probe_id))
  rows <- lapply(probe_ids, function(probe_id) {
    all_probe <- unique_pairs[unique_pairs$probe_id == probe_id, , drop = FALSE]
    syms <- sort(unique(toupper(trimws(as.character(all_probe$symbol)))))
    syms <- syms[!is.na(syms) & nzchar(syms)]
    ents <- sort(unique(trimws(as.character(all_probe$entrez_id))))
    ents <- ents[!is.na(ents) & nzchar(ents)]
    sources <- sort(unique(olfml2b_clean_text(all_probe$mapping_source)))
    sources <- sources[nzchar(sources)]
    symbol_support <- toupper(target_gene) %in% syms
    entrez_support <- target_entrez %in% ents
    conflict_symbols <- setdiff(syms, toupper(target_gene))
    is_ambiguous <- probe_id %in% ambiguous || length(syms) > 1L
    reason <- paste(c(if (symbol_support) "symbol" else NULL,
                      if (entrez_support) "entrez" else NULL), collapse = "+")
    data.frame(
      probe_id = probe_id,
      symbol = if (symbol_support || entrez_support) toupper(target_gene) else if (length(syms)) syms[1] else NA_character_,
      entrez_id = if (entrez_support) target_entrez else if (length(ents)) ents[1] else NA_character_,
      mapping_source = paste(sources, collapse = ";"),
      symbol_annotations = paste(syms, collapse = ";"),
      entrez_annotations = paste(ents, collapse = ";"),
      target_match_reason = reason,
      ambiguous_symbol = is_ambiguous,
      in_expression = probe_id %in% probes,
      primary_eligible = probe_id %in% probes && !is_ambiguous && (symbol_support || entrez_support) && !length(conflict_symbols),
      stringsAsFactors = FALSE
    )
  })
  list(summary = do.call(rbind, rows), evidence = evidence)
}

olfml2b_platform_mapping <- function(platform, probes, dirs, overwrite = FALSE, log_file = NULL) {
  known <- olfml2b_known_platform_packages()
  mapping_parts <- list()
  if (platform %in% names(known)) {
    pkg <- unname(known[[platform]])
    part <- tryCatch(
      olfml2b_annotation_db_mapping(probes, pkg, log_file),
      error = function(e) {
        olfml2b_log("WARN", platform, ": annotation DB mapping failed: ", conditionMessage(e), log_file = log_file)
        NULL
      }
    )
    if (!is.null(part)) mapping_parts[[length(mapping_parts) + 1L]] <- part
  }
  mapped_fraction <- if (length(mapping_parts)) {
    length(intersect(probes, mapping_parts[[1]]$probe_id[!is.na(mapping_parts[[1]]$symbol)])) / length(unique(probes))
  } else 0
  if (mapped_fraction < 0.85 || platform == "GPL8432") {
    gpl <- olfml2b_download_gpl_table(platform, dirs, overwrite, log_file)
    id_col <- olfml2b_case_when_column(gpl, c("id", "id_ref", "probe_id"))
    if (is.na(id_col)) id_col <- names(gpl)[1]
    sym_col <- olfml2b_gpl_candidate_symbol_column(gpl)
    entrez_col <- olfml2b_gpl_candidate_entrez_column(gpl)
    olfml2b_assert(!is.na(sym_col), platform, ": no gene-symbol column found in GPL annotation")
    part <- data.frame(
      probe_id = as.character(gpl[[id_col]]),
      symbol = olfml2b_parse_gpl_symbol(gpl[[sym_col]]),
      entrez_id = if (!is.na(entrez_col)) as.character(gpl[[entrez_col]]) else NA_character_,
      mapping_source = paste0("GEO_AnnotGPL_", platform),
      stringsAsFactors = FALSE
    )
    mapping_parts[[length(mapping_parts) + 1L]] <- part
  }
  olfml2b_assert(length(mapping_parts) > 0L, platform, ": no platform annotation source was available")
  mapping <- do.call(rbind, mapping_parts)
  mapping <- mapping[mapping$probe_id %in% probes, , drop = FALSE]
  mapping$symbol <- toupper(trimws(mapping$symbol))
  unique_pairs <- unique(mapping[, c("probe_id", "symbol", "entrez_id", "mapping_source")])
  valid_symbol <- !is.na(unique_pairs$symbol) & nzchar(unique_pairs$symbol)
  symbol_count <- if (any(valid_symbol)) {
    aggregate(symbol ~ probe_id, unique_pairs[valid_symbol, , drop = FALSE],
              FUN = function(z) length(unique(z)))
  } else data.frame(probe_id = character(), symbol = integer())
  ambiguous <- symbol_count$probe_id[symbol_count$symbol > 1L]
  usable <- unique_pairs[
    !unique_pairs$probe_id %in% ambiguous & valid_symbol,
    , drop = FALSE
  ]
  usable <- usable[!duplicated(usable$probe_id), , drop = FALSE]
  target_audit <- olfml2b_target_probe_audit(
    unique_pairs = unique_pairs,
    probes = probes,
    target_gene = "OLFML2B",
    target_entrez = "25903",
    ambiguous = ambiguous
  )
  target_rows <- target_audit$summary
  target_probe_status <- target_rows[, intersect(
    c("probe_id", "symbol", "entrez_id", "target_match_reason", "ambiguous_symbol", "primary_eligible"),
    names(target_rows)
  ), drop = FALSE]
  list(
    usable = usable,
    all = unique_pairs,
    ambiguous_probes = ambiguous,
    target_rows = target_rows,
    target_annotations = target_audit$evidence,
    target_probe_status = target_probe_status,
    mapped_fraction = length(unique(usable$probe_id)) / length(unique(probes))
  )
}

olfml2b_transform_expression_if_needed <- function(expr) {
  q <- stats::quantile(expr[is.finite(expr)], probs = c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE)
  needs_log <- q[5] > 100 || (q[5] - q[1] > 50 && q[2] >= 0)
  transformed <- expr
  if (needs_log) transformed <- log2(pmax(expr, 0) + 1)
  list(matrix = transformed, was_log2_transformed = needs_log, quantiles_before = q,
       quantiles_after = stats::quantile(transformed[is.finite(transformed)], probs = c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE))
}

olfml2b_aggregate_probes_to_genes <- function(expr_probe, mapping, target_gene = "OLFML2B") {
  olfml2b_assert(is.matrix(expr_probe), "expr_probe must be a matrix")
  usable <- mapping$usable
  usable <- usable[usable$probe_id %in% rownames(expr_probe), , drop = FALSE]
  olfml2b_assert(nrow(usable) > 100L, "Too few unambiguous probe mappings: ", nrow(usable))
  index <- split(match(usable$probe_id, rownames(expr_probe)), usable$symbol)
  gene_mat <- vapply(index, function(idx) {
    if (length(idx) == 1L) as.numeric(expr_probe[idx, ]) else matrixStats::colMedians(expr_probe[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(expr_probe)))
  gene_mat <- t(gene_mat)
  colnames(gene_mat) <- colnames(expr_probe)
  target_probes <- unique(usable$probe_id[toupper(usable$symbol) == toupper(target_gene)])
  olfml2b_assert(length(target_probes) >= 1L, "No unambiguous OLFML2B probe mapped on this platform")
  target_probe_matrix <- expr_probe[target_probes, , drop = FALSE]
  target_probe_stats <- data.frame(
    probe_id = target_probes,
    mean = rowMeans(target_probe_matrix, na.rm = TRUE),
    sd = matrixStats::rowSds(target_probe_matrix, na.rm = TRUE),
    iqr = matrixStats::rowIQRs(target_probe_matrix, na.rm = TRUE),
    missing_fraction = rowMeans(!is.finite(target_probe_matrix)),
    stringsAsFactors = FALSE
  )
  list(matrix = gene_mat, target_probes = target_probes,
       target_probe_matrix = target_probe_matrix, target_probe_stats = target_probe_stats)
}

olfml2b_geo_supplementary_listing <- function(accession, dirs, log_file = NULL) {
  olfml2b_require_packages("GEOquery", log_file)
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  listing <- tryCatch(
    GEOquery::getGEOSuppFiles(accession, makeDirectory = FALSE, baseDir = paths$suppl_dir, fetch_files = FALSE),
    error = function(e) {
      olfml2b_log("WARN", accession, ": supplementary listing failed: ", conditionMessage(e), log_file = log_file)
      NULL
    }
  )
  if (is.null(listing) || !nrow(listing)) return(data.frame())
  out <- data.frame(
    url = rownames(listing),
    size = if ("size" %in% names(listing)) as.numeric(listing$size) else NA_real_,
    filename = basename(rownames(listing)),
    stringsAsFactors = FALSE
  )
  out
}


olfml2b_local_supplementary_files <- function(accession, dirs, max_bytes = 100 * 1024^2) {
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  if (!dir.exists(paths$suppl_dir)) return(data.frame())
  files <- list.files(paths$suppl_dir, pattern = "\\.(xls|xlsx|csv|tsv|txt|soft|json)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  info <- file.info(files)
  keep <- !is.na(info$size) & info$size > 0 & info$size <= max_bytes
  files <- files[keep]
  info <- info[keep, , drop = FALSE]
  if (!length(files)) return(data.frame())
  data.frame(
    url = NA_character_,
    size = as.numeric(info$size),
    filename = basename(files),
    local_path = normalizePath(files, winslash = "/", mustWork = TRUE),
    status = "local_existing",
    source = "local_existing_supplementary",
    stringsAsFactors = FALSE
  )
}

olfml2b_download_small_supplementary <- function(accession, dirs, overwrite = FALSE,
                                               max_bytes = 100 * 1024^2, log_file = NULL) {
  ## True local-first rule: if tabular supplementary files already exist locally,
  ## return them and do not query the GEO supplementary listing.
  local <- if (!overwrite) olfml2b_local_supplementary_files(accession, dirs, max_bytes = max_bytes) else data.frame()
  if (nrow(local)) {
    olfml2b_log("INFO", accession, ": using existing supplementary files; skipping GEO listing/download.", log_file = log_file)
    return(local)
  }

  listing <- olfml2b_geo_supplementary_listing(accession, dirs, log_file)
  if (!nrow(listing)) return(data.frame())
  is_tabular <- grepl("\\.(xls|xlsx|csv|tsv|txt|soft|json)(\\.gz)?$", listing$filename, ignore.case = TRUE)
  keep <- is_tabular & (is.na(listing$size) | listing$size <= max_bytes)
  listing <- listing[keep, , drop = FALSE]
  if (!nrow(listing)) return(data.frame())
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  listing$local_path <- NA_character_
  listing$status <- "not_attempted"
  for (i in seq_len(nrow(listing))) {
    dest <- file.path(paths$suppl_dir, listing$filename[i])
    if (file.exists(dest) && !overwrite && file.info(dest)$size > 0) {
      listing$local_path[i] <- normalizePath(dest, winslash = "/", mustWork = TRUE)
      listing$status[i] <- "local_existing"
      next
    }
    result <- tryCatch(
      olfml2b_download(listing$url[i], dest, overwrite = overwrite, min_bytes = 1L, log_file = log_file),
      error = function(e) NA_character_
    )
    listing$local_path[i] <- result
    listing$status[i] <- if (!is.na(result) && file.exists(result)) "downloaded" else "failed"
  }
  listing
}

olfml2b_download_public_clinical_files <- function(accession, dirs, overwrite = FALSE,
                                                 log_file = NULL) {
  routes <- olfml2b_public_clinical_routes(accession)
  if (!nrow(routes)) return(data.frame())
  paths <- olfml2b_geo_cache_paths(accession, dirs)
  out <- routes
  out$size <- NA_real_
  out$local_path <- NA_character_
  out$status <- "not_attempted"
  for (i in seq_len(nrow(out))) {
    dest <- file.path(paths$suppl_dir, out$filename[i])
    if (file.exists(dest) && !overwrite && file.info(dest)$size > 0) {
      out$local_path[i] <- normalizePath(dest, winslash = "/", mustWork = TRUE)
      out$status[i] <- "local_existing"
      out$size[i] <- file.info(dest)$size
      olfml2b_log("INFO", accession, ": using existing public clinical file; skipping download: ", dest, log_file = log_file)
      next
    }
    result <- tryCatch(
      olfml2b_download(out$url[i], dest, overwrite = overwrite, min_bytes = 100L, log_file = log_file),
      error = function(e) {
        olfml2b_log("WARN", accession, ": public clinical route failed [", out$source[i], "]: ",
                  conditionMessage(e), log_file = log_file)
        NA_character_
      }
    )
    out$local_path[i] <- result
    out$status[i] <- if (!is.na(result) && file.exists(result)) "downloaded" else "failed"
    if (!is.na(result) && file.exists(result)) out$size[i] <- file.info(result)$size
  }
  out[, c("url", "size", "filename", "local_path", "status", "source"), drop = FALSE]
}

olfml2b_bind_download_listings <- function(...) {
  xs <- list(...)
  xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(xs)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::rbindlist(xs, fill = TRUE, use.names = TRUE)))
  }
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}

olfml2b_read_tabular_file <- function(path) {
  olfml2b_assert(file.exists(path), "Missing tabular file: ", path)
  lower <- tolower(path)
  actual <- path
  temp_unzipped <- NULL
  if (grepl("\\.gz$", lower)) {
    temp_unzipped <- tempfile(fileext = sub(".*(\\.[a-z0-9]+)\\.gz$", "\\1", lower))
    R.utils::gunzip(path, destname = temp_unzipped, overwrite = TRUE, remove = FALSE)
    actual <- temp_unzipped
    lower <- sub("\\.gz$", "", lower)
    on.exit(unlink(temp_unzipped), add = TRUE)
  }
  if (grepl("\\.xlsx?$", lower)) {
    sheets <- readxl::excel_sheets(actual)
    tabs <- list()
    for (s in sheets) {
      default <- tryCatch(
        as.data.frame(readxl::read_excel(actual, sheet = s), stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) data.frame()
      )
      if (nrow(default)) tabs[[s]] <- default

      # Several legacy GEO outcome workbooks place a title or merged-cell note
      # above the real header.  Read a raw copy and promote the most plausible
      # clinical header instead of silently treating the title as column names.
      raw <- tryCatch(
        as.data.frame(readxl::read_excel(actual, sheet = s, col_names = FALSE), stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) data.frame()
      )
      if (nrow(raw) >= 2L) {
        scan_n <- min(30L, nrow(raw))
        header_score <- vapply(seq_len(scan_n), function(i) {
          cells <- tolower(olfml2b_clean_text(unlist(raw[i, , drop = TRUE], use.names = FALSE)))
          cells <- cells[nzchar(cells) & cells != "na"]
          if (!length(cells)) return(-Inf)
          clinical_tokens <- sum(grepl(
            "sample|patient|tumou?r|geo|accession|overall|survival|follow|month|day|status|vital|death|dead|age|sex|stage|lauren",
            cells, perl = TRUE
          ))
          20 * clinical_tokens + min(length(unique(cells)), 15L)
        }, numeric(1))
        header_row <- which.max(header_score)
        if (length(header_row) && header_row < nrow(raw) && is.finite(header_score[header_row]) && header_score[header_row] >= 25) {
          hdr <- olfml2b_clean_text(unlist(raw[header_row, , drop = TRUE], use.names = FALSE))
          hdr[!nzchar(hdr) | is.na(hdr)] <- paste0("unnamed_", which(!nzchar(hdr) | is.na(hdr)))
          promoted <- raw[(header_row + 1L):nrow(raw), , drop = FALSE]
          names(promoted) <- make.unique(hdr, sep = "_")
          promoted <- promoted[rowSums(!is.na(promoted) & nzchar(trimws(as.matrix(promoted)))) > 0L, , drop = FALSE]
          if (nrow(promoted)) tabs[[paste0(s, "__detected_header_row_", header_row)]] <- promoted
        }
      }
    }
    return(tabs)
  }
  if (grepl("\\.csv$", lower)) {
    return(list(data = data.table::fread(actual, data.table = FALSE, fill = TRUE)))
  }
  if (grepl("\\.(tsv|txt|soft)$", lower)) {
    tab <- tryCatch(data.table::fread(actual, data.table = FALSE, fill = TRUE), error = function(e) NULL)
    if (is.null(tab)) tab <- utils::read.delim(actual, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
    return(list(data = tab))
  }
  list()
}

olfml2b_keyvalue_metadata <- function(df) {
  x <- olfml2b_expand_characteristics(df)
  # Generate a searchable concatenation per row for fallback extraction.
  x$.all_metadata <- apply(x, 1, function(row) paste(olfml2b_clean_text(row), collapse = " | "))
  x
}

olfml2b_choose_clinical_column <- function(df, concept, cohort = NULL) {
  synonyms <- olfml2b_clinical_synonyms()[[concept]] %||% concept
  nms <- names(df)
  clean <- olfml2b_clean_names(nms)
  score <- rep(0, length(nms))
  for (syn in olfml2b_clean_names(synonyms)) {
    score <- pmax(score,
                  ifelse(clean == syn, 100, 0),
                  ifelse(startsWith(clean, paste0(syn, "_")), 80, 0),
                  ifelse(grepl(syn, clean, fixed = TRUE), 50, 0))
  }
  # Dataset-specific deterministic preferences, embedded rather than supplied by CSV.
  if (!is.null(cohort)) {
    preferred <- switch(cohort,
      GSE62254 = list(
        subtype = c("acr_subtype", "acrg_subtype", "molecular_subtype"),
        os_time = c("overall_survival_months", "os_months", "os_time"),
        os_event = c("overall_survival_status", "os_event", "vital_status"),
        dfs_time = c("recurrence_free_survival_months", "rfs_months", "time_to_recurrence"),
        dfs_event = c("recurrence_status", "rfs_event", "recurrence")
      ),
      GSE15459 = list(
        os_time = c("overall_survival_months", "survival_months", "os_time"),
        os_event = c("vital_status", "survival_status", "os_event")
      ),
      GSE26253 = list(
        dfs_time = c("recurrence_free_survival_time_month", "recurrence_free_survival_time_months", "rfs_months", "recurrence_free_survival_months", "time_to_recurrence"),
        dfs_event = c("status_0_non_recurrence_1_recurrence", "recurrence_status", "rfs_event", "status"),
        os_time = c("overall_survival_time_month", "os_months", "overall_survival_months"),
        os_event = c("overall_survival_status", "os_event", "vital_status")
      ),
      list()
    )
    if (concept %in% names(preferred)) {
      for (i in seq_along(preferred[[concept]])) {
        p <- olfml2b_clean_names(preferred[[concept]][i])
        score[clean == p] <- score[clean == p] + 200 - i
      }
    }
    if (identical(cohort, "GSE26253")) {
      # In GSE26253 the generic field named "status" is explicitly recurrence
      # status, not vital status. Prevent accidental duplication into OS.
      if (concept %in% c("os_event", "os_time")) {
        bad <- grepl("recurrence|relapse|rfs|disease_free", clean) | clean == "status"
        score[bad] <- -Inf
      }
      if (identical(concept, "dfs_event")) {
        score[clean == "status"] <- pmax(score[clean == "status"], 450)
      }
    }
  }
  nonmissing <- vapply(df, function(z) sum(!is.na(z) & nzchar(olfml2b_clean_text(z))), integer(1))
  score <- score + pmin(nonmissing / max(1, nrow(df)) * 20, 20)
  idx <- which.max(score)
  if (!length(idx) || score[idx] < 20) return(NA_character_)
  nms[idx]
}


olfml2b_canonical_join_id <- function(x, numeric_only = FALSE) {
  raw <- toupper(olfml2b_clean_text(x))
  raw <- gsub("\\.0$", "", raw)
  raw <- gsub("^PATIENT[ _-]*", "", raw)
  raw <- gsub("^T(?=[0-9]+$)", "", raw, perl = TRUE)
  raw <- gsub("^TR[_-]*(?=[0-9]+$)", "", raw, perl = TRUE)
  if (numeric_only) {
    m <- regexpr("[0-9]+", raw, perl = TRUE)
    digits <- rep(NA_character_, length(raw))
    hit <- !is.na(m) & m > 0L
    digits[hit] <- regmatches(raw, m)[hit]
    digits <- sub("^0+", "", digits)
    digits[!nzchar(digits)] <- NA_character_
    return(digits)
  }
  gsub("[^A-Z0-9]", "", raw)
}

olfml2b_parse_gse62254_clinical <- function(df) {
  x <- as.data.frame(df, stringsAsFactors = FALSE)
  names(x) <- make.unique(olfml2b_clean_names(names(x)), sep = "_")
  pick <- function(patterns, required = FALSE) {
    for (pat in patterns) {
      hit <- which(names(x) == pat)
      if (length(hit)) return(names(x)[hit[1]])
    }
    for (pat in patterns) {
      hit <- grep(pat, names(x), ignore.case = TRUE, perl = TRUE)
      if (length(hit)) return(names(x)[hit[1]])
    }
    if (required) olfml2b_abort("GSE62254 clinical supplement lacks required field matching: ", paste(patterns, collapse = ";"))
    NA_character_
  }
  get <- function(col, default = NA_character_) {
    if (is.na(col)) rep(default, nrow(x)) else x[[col]]
  }
  cols <- list(
    patient = pick(c("tumor_id", "^scri_no$"), TRUE),
    sex = pick(c("^sex$")),
    age = pick(c("^age$")),
    stage = pick(c("^pstage$", "^stage_tnm$")),
    t = pick(c("^t$")), n = pick(c("^n$")), m = pick(c("^m$")),
    lauren = pick(c("^lauren$", "lauren_1_intestinal")),
    recurrence = pick(c("^documented_recurrence_no_0_yes_1_unknown_2$", "documented_recurrence"), TRUE),
    fu_status = pick(c("^fu_status0_alive_without_ds", "fu_status"), TRUE),
    dfs = pick(c("^dfs_months$", "^dfs"), TRUE),
    os = pick(c("^os_months$", "^os"), TRUE),
    subtype = pick(c("^mol_subtype", "molecular_subtype"), TRUE),
    treatment = pick(c("^adj_ctx_description$", "adjuvant_ccrt"))
  )
  patient_num <- olfml2b_numeric(get(cols$patient))
  patient_id <- ifelse(is.finite(patient_num), as.character(as.integer(patient_num)), olfml2b_clean_text(get(cols$patient)))
  recurrence_raw <- olfml2b_numeric(get(cols$recurrence))
  dfs_event <- ifelse(recurrence_raw == 0, 0L, ifelse(recurrence_raw == 1, 1L, NA_integer_))
  fu <- olfml2b_numeric(get(cols$fu_status))
  os_event <- ifelse(fu %in% c(0, 1, 5), 0L, ifelse(fu %in% c(2, 3, 4), 1L, NA_integer_))
  subtype_num <- olfml2b_numeric(get(cols$subtype))
  subtype <- c("MSS_TP53_inactive", "MSS_TP53_active", "MSI", "EMT")[match(subtype_num, 0:3)]
  stage_raw <- olfml2b_clean_text(get(cols$stage))
  out <- data.frame(
    sample_id = paste0("T", patient_id),
    sample_title = paste0("T", patient_id),
    patient_id = patient_id,
    age = olfml2b_numeric(get(cols$age)),
    sex = olfml2b_sex(get(cols$sex)),
    stage_raw = stage_raw,
    stage = olfml2b_stage_group(stage_raw),
    t_stage = olfml2b_clean_text(get(cols$t)),
    n_stage = olfml2b_clean_text(get(cols$n)),
    m_stage = olfml2b_clean_text(get(cols$m)),
    grade = NA_character_,
    lauren = olfml2b_clean_text(get(cols$lauren)),
    molecular_subtype = subtype,
    treatment = olfml2b_clean_text(get(cols$treatment)),
    tissue = "gastric_primary_tumor",
    os_time_days = olfml2b_numeric(get(cols$os)) * 30.4375,
    os_event = as.integer(os_event),
    dfs_time_days = olfml2b_numeric(get(cols$dfs)) * 30.4375,
    dfs_event = as.integer(dfs_event),
    stringsAsFactors = FALSE
  )
  audit <- data.frame(
    cohort = "GSE62254",
    concept = c("patient_id", "age", "sex", "stage", "t_stage", "n_stage", "m_stage", "lauren", "subtype", "treatment", "os_time", "os_event", "dfs_time", "dfs_event"),
    selected_column = c(cols$patient, cols$age, cols$sex, cols$stage, cols$t, cols$n, cols$m,
                        cols$lauren, cols$subtype, cols$treatment, cols$os, cols$fu_status,
                        cols$dfs, cols$recurrence),
    stringsAsFactors = FALSE
  )
  audit$nonmissing <- vapply(audit$selected_column, function(col) {
    if (is.na(col) || !col %in% names(x)) 0L else sum(!is.na(x[[col]]) & nzchar(olfml2b_clean_text(x[[col]])))
  }, integer(1))
  list(clinical = out, expanded = x, audit = audit, parser = "embedded_exact_GSE62254_ACRG")
}


# ==============================================================================
# GSE62254/ACRG clinical hardening helpers.
# The GEO Series Matrix uses titles such as T107, whereas the Nature Medicine
# supplement uses Tumor ID 107. Generic clinical parsing can accidentally treat
# the numeric part of T107 as age or survival time. These helpers detect that
# leakage and ensure the exact ACRG supplement overwrites pData-derived fields.
# ==============================================================================

olfml2b_numeric_patient_key <- function(x) {
  key <- olfml2b_canonical_join_id(x, numeric_only = TRUE)
  suppressWarnings(as.integer(key))
}

olfml2b_fraction_equal_numeric_patient_key <- function(value, patient_id, tolerance = 1e-8) {
  v <- suppressWarnings(as.numeric(value))
  p <- olfml2b_numeric_patient_key(patient_id)

  # Generic cohorts can contain clinical vectors derived from supplemental
  # files before sample-level merging.  In that situation value and patient_id
  # may legitimately have different lengths.  Leakage auditing is only valid
  # for paired sample-level vectors, so avoid R vector recycling warnings and
  # return zero rather than inventing a suspicious fraction.
  if (!length(v) || !length(p) || length(v) != length(p)) {
    return(0)
  }

  ok <- is.finite(v) & is.finite(p)
  if (!any(ok)) return(0)
  # GEO sample titles such as GC-017TPC-T may be converted by a generic numeric
  # parser into -17.  Treat both signed and absolute matches as identity leakage.
  mean((abs(v[ok] - p[ok]) <= tolerance) | (abs(abs(v[ok]) - p[ok]) <= tolerance))
}

olfml2b_clinical_patient_id_leakage_audit <- function(clinical, cohort) {
  if (!is.data.frame(clinical) || !nrow(clinical)) {
    return(data.frame(
      cohort = cohort,
      variable = character(),
      fraction_equal_numeric_patient_id = numeric(),
      suspicious = logical(),
      stringsAsFactors = FALSE
    ))
  }
  id_source <- if ("patient_id" %in% names(clinical)) clinical$patient_id else if ("sample_title" %in% names(clinical)) clinical$sample_title else clinical$sample_id
  vars <- intersect(c("age", "os_time_days", "dfs_time_days", "rfs_time_days"), names(clinical))
  rows <- lapply(vars, function(v) {
    frac <- olfml2b_fraction_equal_numeric_patient_key(clinical[[v]], id_source)
    data.frame(
      cohort = cohort,
      variable = v,
      fraction_equal_numeric_patient_id = frac,
      suspicious = is.finite(frac) && frac >= 0.60,
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

olfml2b_clinical_plausibility_audit <- function(clinical, cohort) {
  n <- if (is.data.frame(clinical)) nrow(clinical) else 0L
  leak <- olfml2b_clinical_patient_id_leakage_audit(clinical, cohort)
  leak_age <- if (nrow(leak) && any(leak$variable == "age")) leak$fraction_equal_numeric_patient_id[match("age", leak$variable)] else 0
  leak_os <- if (nrow(leak) && any(leak$variable == "os_time_days")) leak$fraction_equal_numeric_patient_id[match("os_time_days", leak$variable)] else 0
  leak_dfs <- if (nrow(leak) && any(leak$variable == "dfs_time_days")) leak$fraction_equal_numeric_patient_id[match("dfs_time_days", leak$variable)] else 0
  data.frame(
    cohort = cohort,
    n = n,
    n_age = if ("age" %in% names(clinical)) sum(is.finite(clinical$age)) else 0L,
    median_age = if ("age" %in% names(clinical) && any(is.finite(clinical$age))) stats::median(clinical$age, na.rm = TRUE) else NA_real_,
    n_os_complete = if (all(c("os_time_days", "os_event") %in% names(clinical))) sum(is.finite(clinical$os_time_days) & clinical$os_time_days > 0 & clinical$os_event %in% c(0L, 1L)) else 0L,
    os_events = if ("os_event" %in% names(clinical)) sum(clinical$os_event == 1L, na.rm = TRUE) else 0L,
    median_os_days = if ("os_time_days" %in% names(clinical) && any(is.finite(clinical$os_time_days) & clinical$os_time_days > 0)) stats::median(clinical$os_time_days[is.finite(clinical$os_time_days) & clinical$os_time_days > 0], na.rm = TRUE) else NA_real_,
    n_dfs_complete = if (all(c("dfs_time_days", "dfs_event") %in% names(clinical))) sum(is.finite(clinical$dfs_time_days) & clinical$dfs_time_days > 0 & clinical$dfs_event %in% c(0L, 1L)) else 0L,
    dfs_events = if ("dfs_event" %in% names(clinical)) sum(clinical$dfs_event == 1L, na.rm = TRUE) else 0L,
    median_dfs_days = if ("dfs_time_days" %in% names(clinical) && any(is.finite(clinical$dfs_time_days) & clinical$dfs_time_days > 0)) stats::median(clinical$dfs_time_days[is.finite(clinical$dfs_time_days) & clinical$dfs_time_days > 0], na.rm = TRUE) else NA_real_,
    patient_id_leak_age = leak_age,
    patient_id_leak_os_time = leak_os,
    patient_id_leak_dfs_time = leak_dfs,
    leakage_fail = any(leak$suspicious %in% TRUE),
    stringsAsFactors = FALSE
  )
}

olfml2b_is_exact_gse62254_acrg_candidate <- function(candidate) {
  is.list(candidate) && is.list(candidate$parsed) && identical(candidate$parsed$parser %||% NA_character_, "embedded_exact_GSE62254_ACRG")
}

olfml2b_gse62254_clinical_fields_to_overwrite <- function() {
  c(
    "age", "sex", "stage_raw", "stage", "t_stage", "n_stage", "m_stage",
    "grade", "lauren", "molecular_subtype", "treatment", "tissue",
    "os_time_days", "os_event", "dfs_time_days", "dfs_event"
  )
}


olfml2b_valid_clinical_category_mask <- function(value, field) {
  x <- toupper(trimws(olfml2b_clean_text(value)))
  x_compact <- gsub("[[:space:]_.-]+", "", x)
  out <- rep(FALSE, length(x))
  if (identical(field, "stage_raw")) {
    out <- grepl("^(AJCC[[:space:]_.-]*)?(PATHOLOGIC(AL)?[[:space:]_.-]*)?(STAGE[[:space:]_.-]*)?(0|I|II|III|IV|[1-4])([ABC])?$", x, perl = TRUE)
  } else if (identical(field, "t_stage")) {
    out <- grepl("^(P?T)?(IS|[0-4][ABC]?)$", x_compact, perl = TRUE)
  } else if (identical(field, "n_stage")) {
    out <- grepl("^(P?N)?[0-3][ABC]?$", x_compact, perl = TRUE)
  } else if (identical(field, "m_stage")) {
    out <- grepl("^(P?M)?[01][ABC]?$", x_compact, perl = TRUE)
  } else if (identical(field, "grade")) {
    out <- grepl("^(G|GRADE)?[1-4]$|WELL|MODERATE|POOR|UNDIFFERENTIATED", x_compact, perl = TRUE)
  } else if (identical(field, "lauren")) {
    out <- grepl("INTESTINAL|DIFFUSE|MIXED|INDETERMINATE", x, perl = TRUE)
  } else if (identical(field, "molecular_subtype")) {
    out <- grepl("MSI|EMT|TP53|MSS|EBV|CIN|GENOMICALLY[[:space:]_.-]*STABLE|CHROMOSOMAL[[:space:]_.-]*INSTABILITY", x, perl = TRUE)
  }
  out[is.na(out)] <- FALSE
  out
}

olfml2b_identity_like_mask <- function(value, clinical, field = NA_character_) {
  if (!length(value) || !is.data.frame(clinical) || !nrow(clinical)) return(rep(FALSE, length(value)))
  if (length(value) != nrow(clinical)) return(rep(FALSE, length(value)))
  value_raw <- toupper(olfml2b_clean_text(value))
  id_cols <- intersect(c("sample_id", "sample_title", "patient_id"), names(clinical))
  if (!length(id_cols)) return(rep(FALSE, length(value)))
  key_value <- olfml2b_canonical_join_id(value_raw)
  num_value <- olfml2b_canonical_join_id(value_raw, numeric_only = TRUE)
  direct <- canonical <- numeric <- rep(FALSE, length(value_raw))

  # Identity leakage is a row-level property.  Comparing a clinical value with
  # every patient identifier in the cohort deletes legitimate low-cardinality
  # variables (for example all T3 values merely because patient 3 exists).
  # Compare only against the identifiers belonging to the same sample.
  for (id_col in id_cols) {
    id_raw <- toupper(olfml2b_clean_text(clinical[[id_col]]))
    key_id <- olfml2b_canonical_join_id(id_raw)
    num_id <- olfml2b_canonical_join_id(id_raw, numeric_only = TRUE)
    direct <- direct | (nzchar(value_raw) & nzchar(id_raw) & value_raw == id_raw)
    canonical <- canonical |
      (!is.na(key_value) & nzchar(key_value) & !is.na(key_id) & nzchar(key_id) & key_value == key_id)
    # Numeric fallback is reserved for ID-like strings with at least two
    # digits.  One-digit TNM/grade codes are never treated as identities.
    numeric <- numeric |
      (!is.na(num_value) & nchar(num_value) >= 2L &
         !is.na(num_id) & nchar(num_id) >= 2L & num_value == num_id)
  }
  plausible_category <- if (!is.na(field) && nzchar(field)) {
    olfml2b_valid_clinical_category_mask(value, field)
  } else {
    rep(FALSE, length(value_raw))
  }
  (direct | canonical | numeric) & !plausible_category
}

olfml2b_sanitize_identity_leakage_fields <- function(clinical, cohort) {
  if (!is.data.frame(clinical) || !nrow(clinical)) return(list(clinical = clinical, audit = data.frame()))
  fields <- intersect(
    c("stage_raw", "t_stage", "n_stage", "m_stage", "grade", "lauren", "molecular_subtype", "treatment", "tissue"),
    names(clinical)
  )
  rows <- list()
  for (field in fields) {
    protected <- olfml2b_valid_clinical_category_mask(clinical[[field]], field)
    mask <- olfml2b_identity_like_mask(clinical[[field]], clinical, field = field)
    n_removed <- sum(mask, na.rm = TRUE)
    if (n_removed > 0L) {
      clinical[[field]][mask] <- NA
    }
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort,
      field = field,
      n_identity_like_values_removed = n_removed,
      n_valid_categories_protected = sum(protected, na.rm = TRUE),
      identity_comparison_scope = "ROW_MATCHED_SAMPLE_PATIENT_IDS_ONLY",
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) rows <- list(data.frame())
  list(clinical = clinical, audit = do.call(rbind, rows))
}

olfml2b_sanitize_gse15459_survival_and_clinical <- function(clinical) {
  if (!is.data.frame(clinical) || !nrow(clinical)) return(list(clinical = clinical, audit = data.frame()))
  original <- clinical
  id_source <- if ("patient_id" %in% names(clinical)) clinical$patient_id else if ("sample_title" %in% names(clinical)) clinical$sample_title else clinical$sample_id
  age_leak <- if ("age" %in% names(clinical)) olfml2b_fraction_equal_numeric_patient_key(clinical$age, id_source) else 0
  os_leak <- if ("os_time_days" %in% names(clinical)) olfml2b_fraction_equal_numeric_patient_key(clinical$os_time_days, id_source) else 0
  dfs_leak <- if ("dfs_time_days" %in% names(clinical)) olfml2b_fraction_equal_numeric_patient_key(clinical$dfs_time_days, id_source) else 0
  implausible_age <- if ("age" %in% names(clinical)) mean(is.finite(clinical$age) & (clinical$age < 18 | clinical$age > 100), na.rm = TRUE) else 0
  implausible_os <- if ("os_time_days" %in% names(clinical)) mean(is.finite(clinical$os_time_days) & clinical$os_time_days <= 0, na.rm = TRUE) else 0
  implausible_dfs <- if ("dfs_time_days" %in% names(clinical)) mean(is.finite(clinical$dfs_time_days) & clinical$dfs_time_days <= 0, na.rm = TRUE) else 0
  hard_fail <- any(c(age_leak, os_leak, dfs_leak) >= 0.60, na.rm = TRUE) ||
    any(c(implausible_age, implausible_os, implausible_dfs) >= 0.25, na.rm = TRUE)

  # The local GSE15459 supplementary outcome table can be mis-scored by the
  # generic parser because sample labels (GC-017TPC-T) contain signed numeric
  # fragments.  If plausibility fails, keep expression and sample identifiers,
  # but explicitly remove clinical/survival endpoints so downstream code cannot
  # treat this cohort as a valid survival validation set.
  if (hard_fail) {
    for (field in intersect(c("age", "os_time_days", "dfs_time_days", "recurrence_time_days"), names(clinical))) {
      clinical[[field]] <- NA_real_
    }
    for (field in intersect(c("os_event", "dfs_event", "recurrence_event"), names(clinical))) {
      clinical[[field]] <- NA_integer_
    }
    for (field in intersect(c("stage_raw", "stage", "t_stage", "n_stage", "m_stage", "grade", "lauren", "molecular_subtype", "treatment"), names(clinical))) {
      clinical[[field]] <- NA
    }
  }
  data.frame_audit <- data.frame(
    cohort = "GSE15459",
    age_leak_fraction = age_leak,
    os_time_leak_fraction = os_leak,
    dfs_time_leak_fraction = dfs_leak,
    implausible_age_fraction = implausible_age,
    nonpositive_os_fraction = implausible_os,
    nonpositive_dfs_fraction = implausible_dfs,
    survival_fields_cleared = hard_fail,
    n_os_complete_before = if (all(c("os_time_days", "os_event") %in% names(original))) sum(is.finite(original$os_time_days) & original$os_time_days > 0 & original$os_event %in% c(0L, 1L)) else 0L,
    n_os_complete_after = if (all(c("os_time_days", "os_event") %in% names(clinical))) sum(is.finite(clinical$os_time_days) & clinical$os_time_days > 0 & clinical$os_event %in% c(0L, 1L)) else 0L,
    stringsAsFactors = FALSE
  )
  list(clinical = clinical, audit = data.frame_audit)
}

olfml2b_enforce_gse15459_authoritative_schema <- function(clinical) {
  olfml2b_assert(is.data.frame(clinical) && nrow(clinical) == 192L,
                "GSE15459 authoritative schema requires exactly 192 retained samples.")
  required <- c("sample_id", "sample_title", "patient_id", "age", "sex", "stage",
                "lauren", "molecular_subtype", "os_time_days", "os_event")
  olfml2b_assert(!length(setdiff(required, names(clinical))),
                "GSE15459 authoritative fields missing: ",
                paste(setdiff(required, names(clinical)), collapse = ";"))
  unsupported_character <- intersect(
    c("stage_raw", "t_stage", "n_stage", "m_stage", "grade", "treatment"),
    names(clinical)
  )
  unsupported_numeric <- intersect(
    c("dfs_time_days", "rfs_time_days", "recurrence_time_days"), names(clinical)
  )
  unsupported_event <- intersect(
    c("dfs_event", "rfs_event", "recurrence_event"), names(clinical)
  )
  before <- c(
    vapply(unsupported_character, function(z) sum(!is.na(clinical[[z]]) & nzchar(trimws(as.character(clinical[[z]])))), integer(1)),
    vapply(unsupported_numeric, function(z) sum(is.finite(suppressWarnings(as.numeric(clinical[[z]])))), integer(1)),
    vapply(unsupported_event, function(z) sum(clinical[[z]] %in% c(0L, 1L), na.rm = TRUE), integer(1))
  )
  for (field in unsupported_character) clinical[[field]] <- NA_character_
  for (field in unsupported_numeric) clinical[[field]] <- NA_real_
  for (field in unsupported_event) clinical[[field]] <- NA_integer_
  clinical$tissue <- "primary_gastric_adenocarcinoma"
  clinical$recurrence_endpoint <- NA_character_
  clinical$analysis_role <- "FORMAL_OS"
  clinical$clinical_schema <- "GSE15459_OFFICIAL_OUTCOME_ONLY"
  clinical$age[!is.finite(clinical$age) | clinical$age < 18 | clinical$age > 100] <- NA_real_
  clinical$os_time_days[!is.finite(clinical$os_time_days) | clinical$os_time_days <= 0] <- NA_real_
  clinical$os_event[!clinical$os_event %in% c(0L, 1L)] <- NA_integer_
  n_complete <- sum(is.finite(clinical$os_time_days) & clinical$os_time_days > 0 &
                      clinical$os_event %in% c(0L, 1L))
  n_events <- sum(clinical$os_event == 1L, na.rm = TRUE)
  olfml2b_assert(n_complete == 191L && n_events == 95L,
                "GSE15459 official OS contract failed after schema enforcement: complete=",
                n_complete, ", events=", n_events, "; expected 191 and 95.")
  audit <- data.frame(
    cohort = "GSE15459",
    field = names(before),
    n_non_authoritative_values_cleared = as.integer(before),
    authoritative_source = "GSE15459_outcome.xls_or_exact_official_merge",
    n_retained_samples = nrow(clinical),
    n_os_complete = n_complete,
    os_events = n_events,
    status = "PASS_AUTHORITATIVE_SCHEMA",
    stringsAsFactors = FALSE
  )
  if (!nrow(audit)) {
    audit <- data.frame(cohort = "GSE15459", field = NA_character_,
                        n_non_authoritative_values_cleared = 0L,
                        authoritative_source = "GSE15459_outcome.xls_or_exact_official_merge",
                        n_retained_samples = nrow(clinical), n_os_complete = n_complete,
                        os_events = n_events, status = "PASS_AUTHORITATIVE_SCHEMA",
                        stringsAsFactors = FALSE)
  }
  list(clinical = clinical, audit = audit)
}

olfml2b_extract_ids_from_metadata <- function(df) {
  sample_col <- olfml2b_choose_clinical_column(df, "sample_id")
  title_col <- olfml2b_choose_clinical_column(df, "title")
  patient_col <- olfml2b_choose_clinical_column(df, "patient_id")
  sample_id <- if (!is.na(sample_col)) as.character(df[[sample_col]]) else rownames(df)
  title <- if (!is.na(title_col)) as.character(df[[title_col]]) else sample_id
  patient_id <- if (!is.na(patient_col)) as.character(df[[patient_col]]) else title
  # Remove only obvious tissue suffixes. The original identifiers remain in
  # sample_id and sample_title for traceability and exact matching.
  patient_id <- gsub("[ _-](T|N|A|B)$", "", patient_id, ignore.case = TRUE)
  patient_id <- olfml2b_clean_text(patient_id)
  data.frame(sample_id = olfml2b_clean_text(sample_id), sample_title = olfml2b_clean_text(title),
             patient_id = patient_id, stringsAsFactors = FALSE)
}



olfml2b_extract_regex_group <- function(text, pattern) {
  hits <- regexec(pattern, text, ignore.case = TRUE, perl = TRUE)
  values <- regmatches(text, hits)
  vapply(values, function(z) if (length(z) >= 2L) z[2L] else NA_character_, character(1))
}

olfml2b_parse_gse26253_clinical <- function(df) {
  x <- olfml2b_keyvalue_metadata(df)
  ids <- olfml2b_extract_ids_from_metadata(x)
  # In this cohort the _A/_B suffix is part of the sample identity, not a
  # tumor/normal aliquot suffix. Removing it creates 44 false patient pairs
  # with discordant stage/outcome records.
  ids$patient_id <- ids$sample_title
  olfml2b_assert(!anyDuplicated(ids$patient_id),
                "GSE26253 patient identity contract failed: sample titles are not unique.")
  choose_exact <- function(candidates) {
    candidates <- olfml2b_clean_names(candidates)
    hit <- candidates[candidates %in% names(x)]
    if (length(hit)) hit[1L] else NA_character_
  }
  event_col <- choose_exact(c(
    "status_0_non_recurrence_1_recurrence",
    "status_1_recurrence_0_non_recurrence",
    "recurrence_status", "rfs_event"
  ))
  time_col <- choose_exact(c(
    "recurrence_free_survival_time_month",
    "recurrence_free_survival_time_months",
    "recurrence_free_survival_month", "rfs_months"
  ))
  stage_col <- choose_exact(c("pathological_stage", "pathologic_stage", "stage"))
  event_text <- if (!is.na(event_col)) as.character(x[[event_col]]) else rep(NA_character_, nrow(x))
  time_text <- if (!is.na(time_col)) as.character(x[[time_col]]) else rep(NA_character_, nrow(x))
  stage_text <- if (!is.na(stage_col)) as.character(x[[stage_col]]) else rep(NA_character_, nrow(x))
  all_text <- x$.all_metadata
  missing_event <- is.na(event_text) | !nzchar(olfml2b_clean_text(event_text))
  if (any(missing_event)) {
    event_text[missing_event] <- olfml2b_extract_regex_group(
      all_text[missing_event],
      "status\\s*\\([^)]*(?:non[- ]?recurrence|recurrence)[^)]*\\)\\s*[:=]\\s*([01])"
    )
  }
  missing_time <- is.na(time_text) | !nzchar(olfml2b_clean_text(time_text))
  if (any(missing_time)) {
    time_text[missing_time] <- olfml2b_extract_regex_group(
      all_text[missing_time],
      "recurrence\\s*free\\s*survival\\s*time\\s*\\(months?\\)\\s*[:=]\\s*([0-9]+(?:\\.[0-9]+)?)"
    )
  }
  missing_stage <- is.na(stage_text) | !nzchar(olfml2b_clean_text(stage_text))
  if (any(missing_stage)) {
    stage_text[missing_stage] <- olfml2b_extract_regex_group(
      all_text[missing_stage],
      "pathological\\s*stage\\s*[:=]\\s*([0-9]+|[IVX]+)"
    )
  }
  rfs_month <- olfml2b_numeric(time_text)
  rfs_event <- olfml2b_binary_event(event_text)
  n <- nrow(x)
  out <- cbind(
    ids,
    data.frame(
      age = rep(NA_real_, n),
      sex = factor(rep(NA_character_, n), levels = c("Female", "Male")),
      stage_raw = stage_text,
      stage = olfml2b_stage_group(stage_text),
      t_stage = rep(NA_character_, n), n_stage = rep(NA_character_, n), m_stage = rep(NA_character_, n),
      grade = rep(NA_character_, n), lauren = rep(NA_character_, n), molecular_subtype = rep(NA_character_, n),
      treatment = rep("surgery_plus_adjuvant_chemoradiotherapy", n),
      tissue = rep("FFPE_primary_gastric_adenocarcinoma", n),
      os_time_days = rep(NA_real_, n), os_event = rep(NA_integer_, n),
      dfs_time_days = rfs_month * 30.4375,
      dfs_event = rfs_event,
      rfs_time_months = rfs_month,
      rfs_time_days = rfs_month * 30.4375,
      rfs_event = rfs_event,
      recurrence_endpoint = rep("RFS", n),
      stringsAsFactors = FALSE
    )
  )
  audit <- data.frame(
    cohort = "GSE26253",
    concept = c("sample_id", "patient_id", "stage", "dfs_time", "dfs_event", "os_time", "os_event"),
    selected_column = c(
      olfml2b_choose_clinical_column(x, "sample_id", "GSE26253"),
      olfml2b_choose_clinical_column(x, "patient_id", "GSE26253"),
      stage_col, time_col, event_col, NA_character_, NA_character_
    ),
    nonmissing = c(
      sum(nzchar(ids$sample_id)), sum(nzchar(ids$patient_id)),
      sum(!is.na(out$stage)), sum(is.finite(out$dfs_time_days)),
      sum(out$dfs_event %in% c(0L, 1L)), 0L, 0L
    ),
    parser = "embedded_exact_GSE26253_RFS_colon_first",
    stringsAsFactors = FALSE
  )
  list(clinical = out, expanded = x, audit = audit, parser = "embedded_exact_GSE26253_RFS_colon_first")
}

olfml2b_parse_gse15459_clinical <- function(df) {
  x <- olfml2b_keyvalue_metadata(df)
  ids <- olfml2b_extract_ids_from_metadata(x)
  clean <- olfml2b_clean_names(names(x))

  pick_named <- function(candidates, value_role = c("any", "time", "event")) {
    value_role <- match.arg(value_role)
    cand <- olfml2b_clean_names(candidates)
    idx <- which(clean %in% cand)
    if (!length(idx)) {
      idx <- which(vapply(clean, function(z) any(startsWith(z, paste0(cand, "_"))), logical(1)))
    }
    if (!length(idx)) return(NA_character_)
    for (j in idx) {
      v <- x[[j]]
      if (value_role == "time") {
        num <- olfml2b_numeric(v)
        if (sum(is.finite(num)) >= 20L && length(unique(num[is.finite(num)])) >= 8L && max(num, na.rm = TRUE) > 1) return(names(x)[j])
      } else if (value_role == "event") {
        ev <- olfml2b_binary_event(v)
        if (sum(ev %in% c(0L, 1L), na.rm = TRUE) >= 20L && length(unique(stats::na.omit(ev))) >= 2L) return(names(x)[j])
      } else {
        return(names(x)[j])
      }
    }
    NA_character_
  }

  time_col <- pick_named(c(
    "os_months", "os_month", "overall_survival_months", "overall_survival_month",
    "overall_survival_time_months", "overall_survival_time", "survival_months",
    "survival_month", "survival_time", "follow_up_months", "followup_months", "os"
  ), "time")
  event_col <- pick_named(c(
    "os_event", "overall_survival_event", "overall_survival_status", "survival_status",
    "vital_status", "death_event", "death", "dead", "status",
    # Official GSE15459 supplementary workbook header:
    # "Outcome (1=dead)" -> outcome_1_dead after name cleaning.
    "outcome_1_dead", "outcome"
  ), "event")

  gsm_col <- pick_named(c("gsm_id", "geo_accession", "sample_id"))
  patient_col <- pick_named(c("id", "patient_id", "subject_id", "case_id"))
  cel_col <- pick_named(c("expression_cel_file", "cel_file", "expression_file"))
  age_col <- pick_named(c("age_at_surgery", "age"))
  sex_col <- pick_named(c("gender", "sex"))
  stage_col <- pick_named(c("stage", "pathologic_stage", "pathological_stage"))
  lauren_col <- pick_named(c("laurenclassification", "lauren_classification", "lauren"))
  subtype_col <- pick_named(c("subtype", "molecular_subtype"))

  get <- function(col, default = NA_character_) {
    if (is.na(col) || !col %in% names(x)) rep(default, nrow(x)) else as.character(x[[col]])
  }
  if (!is.na(gsm_col)) ids$sample_id <- olfml2b_clean_text(get(gsm_col))
  if (!is.na(patient_col)) ids$patient_id <- olfml2b_clean_text(get(patient_col))
  # CEL names are retained as an additional exact join key only when present;
  # GSM is still the canonical sample identifier used to merge with GEO pData.
  if (!is.na(cel_col)) ids$sample_title <- olfml2b_clean_text(sub("\\.[Cc][Ee][Ll](\\.gz)?$", "", basename(get(cel_col))))

  resolved <- !is.na(time_col) && !is.na(event_col) && !identical(time_col, event_col)
  n <- nrow(x)
  out <- cbind(
    ids,
    data.frame(
      age = olfml2b_numeric(get(age_col)),
      sex = olfml2b_sex(get(sex_col)),
      stage_raw = olfml2b_clean_text(get(stage_col)),
      stage = olfml2b_stage_group(get(stage_col)),
      t_stage = rep(NA_character_, n), n_stage = rep(NA_character_, n), m_stage = rep(NA_character_, n),
      grade = rep(NA_character_, n), lauren = olfml2b_clean_text(get(lauren_col)), molecular_subtype = olfml2b_clean_text(get(subtype_col)),
      treatment = rep(NA_character_, n), tissue = rep("primary_gastric_tumor", n),
      os_time_days = if (resolved) olfml2b_time_to_days(x[[time_col]], time_col, unit_override = "months") else rep(NA_real_, n),
      os_event = if (resolved) olfml2b_binary_event(x[[event_col]]) else rep(NA_integer_, n),
      dfs_time_days = rep(NA_real_, n), dfs_event = rep(NA_integer_, n),
      stringsAsFactors = FALSE
    )
  )

  complete <- is.finite(out$os_time_days) & out$os_time_days > 0 & out$os_event %in% c(0L, 1L)
  leak <- olfml2b_fraction_equal_numeric_patient_key(out$os_time_days, out$patient_id)
  resolved <- resolved && sum(complete) >= 30L && sum(out$os_event[complete] == 1L, na.rm = TRUE) >= 20L && leak < 0.60
  if (!resolved) {
    out$os_time_days <- NA_real_
    out$os_event <- NA_integer_
  }

  audit <- data.frame(
    cohort = "GSE15459",
    concept = c("sample_id", "patient_id", "cel_file", "age", "sex", "stage", "lauren", "subtype", "os_time", "os_event"),
    selected_column = c(
      gsm_col, patient_col, cel_col, age_col, sex_col, stage_col, lauren_col, subtype_col, time_col, event_col
    ),
    nonmissing = c(
      sum(nzchar(ids$sample_id)), sum(nzchar(ids$patient_id)), sum(nzchar(get(cel_col))),
      sum(is.finite(out$age)), sum(!is.na(out$sex)), sum(!is.na(out$stage)),
      sum(nzchar(out$lauren)), sum(nzchar(out$molecular_subtype)),
      sum(is.finite(out$os_time_days)), sum(out$os_event %in% c(0L, 1L), na.rm = TRUE)
    ),
    parser = "embedded_exact_GSE15459_official_outcome_schema",
    status = ifelse(resolved, "RESOLVED_EVALUABLE_OS", "UNRESOLVED_OR_IMPLAUSIBLE_OS"),
    stringsAsFactors = FALSE
  )
  list(clinical = out, expanded = x, audit = audit, parser = "embedded_exact_GSE15459_official_outcome_schema", resolved = resolved)
}


olfml2b_frozen_time_unit <- function(cohort, source_column) {
  n <- if (is.null(source_column) || !length(source_column) || is.na(source_column[1])) "" else tolower(as.character(source_column[1]))
  if (grepl("day", n)) return("days")
  if (grepl("month|mos|mth", n)) return("months")
  if (grepl("year|yrs", n)) return("years")
  switch(cohort,
    GSE62254 = "months", GSE15459 = "months", GSE26253 = "months",
    GSE84437 = "months", NA_character_)
}

olfml2b_parse_clinical_dataframe <- function(df, cohort) {
  if (identical(cohort, "GSE26253")) {
    return(olfml2b_parse_gse26253_clinical(df))
  }
  # The ACRG article supplement has a stable published schema. Use an exact
  # parser rather than heuristic field selection whenever that schema is seen.
  clean_input_names <- olfml2b_clean_names(names(df))
  if (identical(cohort, "GSE62254") &&
      any(clean_input_names == "tumor_id") &&
      any(grepl("^dfs", clean_input_names)) && any(grepl("^os", clean_input_names))) {
    return(olfml2b_parse_gse62254_clinical(df))
  }
  if (identical(cohort, "GSE15459")) {
    exact <- tryCatch(olfml2b_parse_gse15459_clinical(df), error = function(e) NULL)
    if (!is.null(exact) && isTRUE(exact$resolved)) return(exact)
  }
  x <- olfml2b_keyvalue_metadata(df)
  ids <- olfml2b_extract_ids_from_metadata(x)
  get_col <- function(concept) {
    col <- olfml2b_choose_clinical_column(x, concept, cohort)
    if (is.na(col)) rep(NA_character_, nrow(x)) else as.character(x[[col]])
  }
  os_time_col <- olfml2b_choose_clinical_column(x, "os_time", cohort)
  dfs_time_col <- olfml2b_choose_clinical_column(x, "dfs_time", cohort)
  out <- cbind(
    ids,
    data.frame(
      age = olfml2b_numeric(get_col("age")),
      sex = olfml2b_sex(get_col("sex")),
      stage_raw = get_col("stage"),
      stage = olfml2b_stage_group(get_col("stage")),
      t_stage = get_col("t_stage"),
      n_stage = get_col("n_stage"),
      m_stage = get_col("m_stage"),
      grade = get_col("grade"),
      lauren = get_col("lauren"),
      molecular_subtype = get_col("subtype"),
      treatment = get_col("treatment"),
      tissue = get_col("tissue"),
      os_time_days = olfml2b_time_to_days(get_col("os_time"), os_time_col %||% "",
                                         unit_override = olfml2b_frozen_time_unit(cohort, os_time_col)),
      os_event = olfml2b_binary_event(get_col("os_event")),
      dfs_time_days = olfml2b_time_to_days(get_col("dfs_time"), dfs_time_col %||% "",
                                          unit_override = olfml2b_frozen_time_unit(cohort, dfs_time_col)),
      dfs_event = olfml2b_binary_event(get_col("dfs_event")),
      stringsAsFactors = FALSE
    )
  )
  audit <- data.frame(
    cohort = cohort,
    concept = names(olfml2b_clinical_synonyms()),
    selected_column = vapply(names(olfml2b_clinical_synonyms()), function(z) olfml2b_choose_clinical_column(x, z, cohort), character(1)),
    stringsAsFactors = FALSE
  )
  audit$nonmissing <- vapply(audit$selected_column, function(col) {
    if (is.na(col) || !col %in% names(x)) 0L else sum(!is.na(x[[col]]) & nzchar(olfml2b_clean_text(x[[col]])))
  }, integer(1))
  list(clinical = out, expanded = x, audit = audit)
}

olfml2b_clinical_candidate_quality <- function(tab, cohort, source_label) {
  parsed <- olfml2b_parse_clinical_dataframe(tab, cohort)
  clin <- parsed$clinical
  plausibility <- olfml2b_clinical_plausibility_audit(clin, cohort)
  score <- 0
  score <- score + min(sum(is.finite(clin$age)), 50) / 10
  score <- score + min(sum(!is.na(clin$stage)), 100) / 10
  score <- score + min(sum(is.finite(clin$os_time_days) & clin$os_time_days > 0 & clin$os_event %in% c(0, 1)), 150) / 5
  score <- score + min(sum(is.finite(clin$dfs_time_days) & clin$dfs_time_days > 0 & clin$dfs_event %in% c(0, 1)), 150) / 5
  score <- score + min(sum(nzchar(clin$patient_id)), 100) / 20
  if (identical(cohort, "GSE62254") && identical(parsed$parser %||% NA_character_, "embedded_exact_GSE62254_ACRG")) {
    # The ACRG supplement is the authoritative clinical source for GSE62254.
    # Give it priority over pData and generic tables, provided the fields pass
    # patient-ID leakage checks.
    score <- score + 100
  }
  if (identical(cohort, "GSE15459") &&
      identical(parsed$parser %||% NA_character_, "embedded_exact_GSE15459_official_outcome_schema") &&
      isTRUE(parsed$resolved)) {
    # The official GEO outcome workbook is authoritative and must outrank
    # heuristic parsing of pData/sample labels.
    score <- score + 100
  }
  if (nrow(plausibility) && isTRUE(plausibility$leakage_fail[1])) {
    score <- score - 1e6
  }
  list(score = score, source = source_label, parsed = parsed, plausibility = plausibility)
}

olfml2b_merge_clinical_sources <- function(primary_parsed, supplementary_candidates, cohort) {
  primary <- primary_parsed$clinical
  candidates <- list(list(score = 0, source = "GEO_pData", parsed = primary_parsed))
  if (length(supplementary_candidates)) candidates <- c(candidates, supplementary_candidates)
  # Start with pData and iteratively fill missing fields from the highest-quality
  # downloaded clinical table. The merge key is chosen deterministically.
  candidates <- candidates[order(vapply(candidates, `[[`, numeric(1), "score"), decreasing = TRUE)]
  merged <- primary
  source_map <- data.frame(field = names(merged), source = "GEO_pData", stringsAsFactors = FALSE)
  for (cand in candidates) {
    sup <- cand$parsed$clinical
    if (!nrow(sup)) next
    key_options <- c("sample_id", "sample_title", "patient_id")
    best_key <- NA_character_
    best_matches <- 0L
    best_idx <- NULL
    for (key in key_options) {
      if (!key %in% names(merged) || !key %in% names(sup)) next
      direct_idx <- match(olfml2b_clean_text(merged[[key]]), olfml2b_clean_text(sup[[key]]))
      direct_matches <- sum(!is.na(direct_idx))
      if (direct_matches > best_matches) {
        best_matches <- direct_matches
        best_key <- key
        best_idx <- direct_idx
      }
      # Numeric patient-key fallback is restricted to ACRG, where GEO titles
      # are T117 and the published clinical supplement uses Tumor ID 117.
      if (identical(cohort, "GSE62254") && key %in% c("sample_title", "patient_id")) {
        norm_idx <- match(olfml2b_canonical_join_id(merged[[key]], numeric_only = TRUE),
                          olfml2b_canonical_join_id(sup[[key]], numeric_only = TRUE))
        norm_matches <- sum(!is.na(norm_idx))
        if (norm_matches > best_matches) {
          best_matches <- norm_matches
          best_key <- key
          best_idx <- norm_idx
        }
      }
    }
    if (is.na(best_key) || best_matches < max(5L, floor(0.1 * nrow(merged)))) next
    idx <- best_idx
    authoritative_gse62254 <- identical(cohort, "GSE62254") && olfml2b_is_exact_gse62254_acrg_candidate(cand)
    authoritative_gse15459 <- identical(cohort, "GSE15459") &&
      identical(cand$parsed$parser %||% NA_character_, "embedded_exact_GSE15459_official_outcome_schema") &&
      isTRUE(cand$parsed$resolved)
    overwrite_fields <- c(
      if (authoritative_gse62254) olfml2b_gse62254_clinical_fields_to_overwrite() else character(),
      if (authoritative_gse15459) c("age", "sex", "stage_raw", "stage", "lauren", "molecular_subtype", "os_time_days", "os_event") else character()
    )
    for (field in intersect(names(merged), names(sup))) {
      if (field %in% key_options) next
      current <- merged[[field]]
      incoming <- sup[[field]][idx]
      missing <- if (is.numeric(current)) !is.finite(current) else is.na(current) | !nzchar(olfml2b_clean_text(current))
      valid <- if (is.numeric(incoming)) is.finite(incoming) else !is.na(incoming) & nzchar(olfml2b_clean_text(incoming))
      # For GSE62254, pData often contains sample titles such as T107. Generic
      # parsers can convert that numeric suffix into age or survival time. The
      # exact ACRG supplement must therefore overwrite clinical fields rather
      # than only fill missing cells.
      overwrite <- field %in% overwrite_fields
      take <- (missing | overwrite) & valid
      if (any(take, na.rm = TRUE)) {
        merged[[field]][take] <- incoming[take]
        source_map$source[source_map$field == field] <- paste0(source_map$source[source_map$field == field], ";", cand$source)
      }
    }
  }
  list(clinical = merged, source_map = source_map, candidates = candidates)
}

olfml2b_build_supplementary_clinical_candidates <- function(download_listing, cohort, log_file = NULL) {
  out <- list()
  if (!nrow(download_listing)) return(out)
  for (i in seq_len(nrow(download_listing))) {
    path <- download_listing$local_path[i]
    if (is.na(path) || !file.exists(path)) next
    tabs <- tryCatch(olfml2b_read_tabular_file(path), error = function(e) list())
    if (!length(tabs)) next
    for (nm in names(tabs)) {
      tab <- tabs[[nm]]
      if (!is.data.frame(tab) || nrow(tab) < 5L || ncol(tab) < 2L) next
      candidate <- tryCatch(
        olfml2b_clinical_candidate_quality(tab, cohort, paste0(basename(path), "::", nm)),
        error = function(e) NULL
      )
      if (!is.null(candidate) && candidate$score > 1) out[[length(out) + 1L]] <- candidate
    }
  }
  out
}

olfml2b_apply_known_geo_exclusions <- function(clinical, expr, cohort) {
  removed <- data.frame()
  if (cohort != "GSE15459") return(list(clinical = clinical, expr = expr, removed = removed))
  exclusions <- olfml2b_gse15459_exclusions()
  # GEO appends "[EXCLUDED]" to these eight titles, while the frozen official
  # exclusion registry deliberately stores the original sample title.  Match
  # canonical titles so the label itself cannot defeat the exclusion gate.
  canonical_title <- function(x) {
    x <- gsub("\\s*\\[EXCLUDED\\]\\s*$", "", as.character(x), ignore.case = TRUE)
    toupper(trimws(x))
  }
  match_key <- canonical_title(clinical$sample_title)
  exclusion_key <- canonical_title(exclusions$sample_title)
  match_title <- match_key %in% exclusion_key
  if (any(match_title)) {
    removed <- clinical[match_title, c("sample_id", "sample_title"), drop = FALSE]
    removed$canonical_sample_title <- canonical_title(removed$sample_title)
    removed$reason <- exclusions$reason[match(removed$canonical_sample_title, exclusion_key)]
    keep <- !match_title
    clinical <- clinical[keep, , drop = FALSE]
    expr <- expr[, clinical$sample_id, drop = FALSE]
  }
  olfml2b_assert(nrow(removed) == 8L,
                "GSE15459 official exclusion contract failed: expected 8 removed samples, observed ",
                nrow(removed), ". Refusing to retain GEO samples flagged [EXCLUDED].")
  olfml2b_assert(nrow(clinical) == 192L && ncol(expr) == 192L,
                "GSE15459 post-exclusion sample contract failed: clinical=", nrow(clinical),
                ", expression=", ncol(expr), ", expected 192 each.")
  list(clinical = clinical, expr = expr, removed = removed)
}

olfml2b_geo_normalize_sample_order <- function(expr, clinical) {
  olfml2b_assert(!is.null(colnames(expr)), "Expression matrix lacks sample names")
  idx <- match(colnames(expr), clinical$sample_id)
  if (anyNA(idx)) {
    idx_title <- match(colnames(expr), clinical$sample_title)
    idx[is.na(idx)] <- idx_title[is.na(idx)]
  }
  keep <- !is.na(idx)
  expr <- expr[, keep, drop = FALSE]
  clinical <- clinical[idx[keep], , drop = FALSE]
  clinical$sample_id <- colnames(expr)
  rownames(clinical) <- clinical$sample_id
  list(expr = expr, clinical = clinical)
}

olfml2b_geo_cohort_audit <- function(accession, platform, expr_probe, expr_gene, clinical,
                                   mapping, target_probes, transformed, removed) {
  os_valid <- is.finite(clinical$os_time_days) & clinical$os_time_days > 0 & clinical$os_event %in% c(0L, 1L)
  dfs_valid <- is.finite(clinical$dfs_time_days) & clinical$dfs_time_days > 0 & clinical$dfs_event %in% c(0L, 1L)
  data.frame(
    cohort = accession,
    platform = platform,
    n_probe_features = nrow(expr_probe),
    n_gene_features = nrow(expr_gene),
    n_samples_expression = ncol(expr_gene),
    n_samples_clinical = nrow(clinical),
    n_olfml2b_probes = length(target_probes),
    olfml2b_probes = paste(target_probes, collapse = ";"),
    mapped_probe_fraction = mapping$mapped_fraction,
    log2_transformed_by_pipeline = transformed$was_log2_transformed,
    n_os_complete = sum(os_valid),
    os_events = sum(clinical$os_event[os_valid] == 1L, na.rm = TRUE),
    recurrence_endpoint = if (accession == "GSE62254") "DFS" else if (accession == "GSE26253") "RFS" else NA_character_,
    n_dfs_complete = sum(dfs_valid),
    dfs_events = sum(clinical$dfs_event[dfs_valid] == 1L, na.rm = TRUE),
    n_recurrence_complete = sum(dfs_valid),
    recurrence_events = sum(clinical$dfs_event[dfs_valid] == 1L, na.rm = TRUE),
    n_stage = sum(!is.na(clinical$stage)),
    n_age = sum(is.finite(clinical$age)),
    n_removed_known = nrow(removed),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Part2 - GEO bulk cohorts: GSE62254, GSE15459 and GSE26253
# Official Series Matrix and supplementary clinical files are downloaded and
# harmonized automatically. Dataset rules are embedded in code.
# ==============================================================================

olfml2b_process_one_geo_bulk <- function(accession, ctx, log_file, error_file) {
  manifest <- olfml2b_embedded_geo_manifest()
  row <- manifest[manifest$accession == accession, , drop = FALSE]
  olfml2b_assert(nrow(row) == 1L, "Unknown embedded GEO cohort: ", accession)
  platform_expected <- row$platform
  eset <- olfml2b_run_stage(
    paste0(accession, " Series Matrix"),
    olfml2b_geo_load_eset(
      accession, ctx$dirs, expected_platform = platform_expected,
      overwrite = ctx$overwrite_downloads, log_file = log_file
    ),
    log_file, error_file, required = TRUE
  )
  expr_probe <- Biobase::exprs(eset)
  storage.mode(expr_probe) <- "numeric"
  transformed <- olfml2b_transform_expression_if_needed(expr_probe)
  expr_probe <- transformed$matrix
  platform <- olfml2b_geo_platform_from_eset(eset) %||% platform_expected
  if (is.na(platform) || !nzchar(platform)) platform <- platform_expected
  olfml2b_log("INFO", accession, ": platform=", platform, " | probes=", nrow(expr_probe), " | samples=", ncol(expr_probe), log_file = log_file)

  mapping <- olfml2b_run_stage(
    paste0(accession, " platform mapping"),
    olfml2b_platform_mapping(platform, rownames(expr_probe), ctx$dirs,
                           overwrite = ctx$overwrite_downloads, log_file = log_file),
    log_file, error_file, required = TRUE
  )
  aggregated <- olfml2b_run_stage(
    paste0(accession, " probe aggregation"),
    olfml2b_aggregate_probes_to_genes(expr_probe, mapping, target_gene = ctx$contract$target_gene),
    log_file, error_file, required = TRUE
  )
  expr_gene <- aggregated$matrix

  pdata <- Biobase::pData(eset)
  primary_parsed <- olfml2b_parse_clinical_dataframe(pdata, accession)
  geo_supplementary_listing <- olfml2b_download_small_supplementary(
    accession, ctx$dirs, overwrite = ctx$overwrite_downloads,
    max_bytes = if (accession == "GSE15459") 150 * 1024^2 else 100 * 1024^2,
    log_file = log_file
  )
  public_clinical_listing <- olfml2b_download_public_clinical_files(
    accession, ctx$dirs, overwrite = ctx$overwrite_downloads, log_file = log_file
  )
  supplementary_listing <- olfml2b_bind_download_listings(
    geo_supplementary_listing, public_clinical_listing
  )
  supplementary_candidates <- olfml2b_build_supplementary_clinical_candidates(
    supplementary_listing, accession, log_file
  )
  merged_clinical <- olfml2b_merge_clinical_sources(primary_parsed, supplementary_candidates, accession)
  # The same cached workbook can legitimately be discovered through both the
  # GEO-local and article-supplement routes. Keep provenance once per source.
  if (is.data.frame(merged_clinical$source_map) && nrow(merged_clinical$source_map)) {
    merged_clinical$source_map$source <- vapply(
      strsplit(as.character(merged_clinical$source_map$source), ";", fixed = TRUE),
      function(z) paste(unique(z[nzchar(z)]), collapse = ";"),
      character(1)
    )
  }
  clinical <- merged_clinical$clinical
  clinical_plausibility_pre_order <- olfml2b_clinical_plausibility_audit(clinical, accession)
  if (identical(accession, "GSE62254")) {
    if (isTRUE(clinical_plausibility_pre_order$leakage_fail[1])) {
      olfml2b_abort(
        "GSE62254 clinical harmonisation failed patient-ID leakage audit before sample ordering. ",
        "The ACRG clinical supplement was not used to overwrite pData-derived fields. ",
        "Check Part2_GSE62254_supplementary_download_audit.csv and the Nature Medicine supplement download."
      )
    }
    olfml2b_assert(
      clinical_plausibility_pre_order$n_os_complete[1] >= 250L && clinical_plausibility_pre_order$os_events[1] >= 50L,
      "GSE62254 OS clinical repair failed: complete OS n=", clinical_plausibility_pre_order$n_os_complete[1],
      ", events=", clinical_plausibility_pre_order$os_events[1]
    )
    olfml2b_assert(
      clinical_plausibility_pre_order$n_dfs_complete[1] >= 250L && clinical_plausibility_pre_order$dfs_events[1] >= 50L,
      "GSE62254 DFS clinical repair failed: complete DFS n=", clinical_plausibility_pre_order$n_dfs_complete[1],
      ", events=", clinical_plausibility_pre_order$dfs_events[1]
    )
  }
  if (identical(accession, "GSE15459")) {
    # The official outcome workbook contains 192 outcome rows (191 with
    # positive survival time and binary outcome). Never continue silently with
    # expression-only GSE15459 when the clinical supplement was missed or an
    # older parser was installed.
    olfml2b_assert(
      clinical_plausibility_pre_order$n_os_complete[1] >= 180L &&
        clinical_plausibility_pre_order$os_events[1] >= 80L,
      "GSE15459 official outcome repair failed: complete OS n=",
      clinical_plausibility_pre_order$n_os_complete[1], ", events=",
      clinical_plausibility_pre_order$os_events[1],
      ". Expected the cached/downloaded GSE15459_outcome.xls parsed and authoritatively merged with the v2.3.5 schema."
    )
  }
  ordered <- olfml2b_geo_normalize_sample_order(expr_gene, clinical)
  expr_gene <- ordered$expr
  clinical <- ordered$clinical
  excluded <- olfml2b_apply_known_geo_exclusions(clinical, expr_gene, accession)
  clinical <- excluded$clinical
  expr_gene <- excluded$expr

  if (identical(accession, "GSE15459")) {
    # Preserve the exact official OS/age/sex/stage/Lauren/subtype merge, while
    # explicitly clearing fields that the official outcome sheet does not
    # support. This prevents sample identifiers and signed label fragments from
    # leaking into TNM, DFS, grade or treatment columns.
    gse15459_sanitized <- olfml2b_enforce_gse15459_authoritative_schema(clinical)
    clinical <- gse15459_sanitized$clinical
    identity_sanitized <- list(clinical = clinical, audit = gse15459_sanitized$audit)
  } else {
    identity_sanitized <- olfml2b_sanitize_identity_leakage_fields(clinical, accession)
    clinical <- identity_sanitized$clinical
    gse15459_sanitized <- list(clinical = clinical, audit = data.frame())
  }

  clinical_plausibility <- olfml2b_clinical_plausibility_audit(clinical, accession)
  if (identical(accession, "GSE15459")) {
    olfml2b_assert(
      clinical_plausibility$n_os_complete[1] == 191L && clinical_plausibility$os_events[1] == 95L,
      "GSE15459 valid OS was lost after sample ordering: n=", clinical_plausibility$n_os_complete[1],
      ", events=", clinical_plausibility$os_events[1]
    )
  }
  if (identical(accession, "GSE26253")) {
    olfml2b_assert(nrow(clinical) == 432L && length(unique(clinical$patient_id)) == 432L,
                  "GSE26253 identity contract failed after ordering: samples=", nrow(clinical),
                  ", unique patient_id=", length(unique(clinical$patient_id)), ".")
  }
  if (identical(accession, "GSE62254")) {
    if (isTRUE(clinical_plausibility$leakage_fail[1])) {
      olfml2b_abort("GSE62254 clinical harmonisation failed patient-ID leakage audit after sample ordering.")
    }
    # Fail immediately if categorical TNM/subtype values were mistaken for
    # patient identifiers.  The official ACRG supplement contains near-complete
    # pT/pM, complete pN and molecular-subtype annotations with the frozen level
    # sets below.  These checks deliberately test information content, not only
    # the number of rows.
    p_t_gate <- olfml2b_part2_t_component(clinical$t_stage)
    p_n_gate <- olfml2b_part2_n_component(clinical$n_stage)
    p_m_gate <- olfml2b_part2_m_component(clinical$m_stage)
    subtype_gate <- olfml2b_clean_text(clinical$molecular_subtype)
    tnm_subtype_ok <-
      sum(!is.na(p_t_gate)) >= 298L && all(c("T2", "T3", "T4") %in% unique(stats::na.omit(p_t_gate))) &&
      sum(!is.na(p_n_gate)) == 300L && all(c("N0", "N1", "N2", "N3") %in% unique(stats::na.omit(p_n_gate))) &&
      sum(!is.na(p_m_gate)) >= 297L && all(c("M0", "M1") %in% unique(stats::na.omit(p_m_gate))) &&
      sum(nzchar(subtype_gate)) == 300L &&
      all(c("MSS_TP53_inactive", "MSS_TP53_active", "MSI", "EMT") %in% unique(subtype_gate))
    olfml2b_assert(
      tnm_subtype_ok,
      "GSE62254 categorical-integrity gate failed after identity sanitization: pT n=",
      sum(!is.na(p_t_gate)), "/levels=", paste(sort(unique(stats::na.omit(p_t_gate))), collapse = "/"),
      "; pN n=", sum(!is.na(p_n_gate)), "/levels=", paste(sort(unique(stats::na.omit(p_n_gate))), collapse = "/"),
      "; pM n=", sum(!is.na(p_m_gate)), "/levels=", paste(sort(unique(stats::na.omit(p_m_gate))), collapse = "/"),
      "; subtype n=", sum(nzchar(subtype_gate)), "/levels=", paste(sort(unique(subtype_gate[nzchar(subtype_gate)])), collapse = "/"),
      ". Legitimate TNM/subtype categories may have been removed as identifiers."
    )
    olfml2b_log(
      "INFO", "GSE62254 categorical-integrity gate PASS | pT=", sum(!is.na(p_t_gate)),
      " | pN=", sum(!is.na(p_n_gate)), " | pM=", sum(!is.na(p_m_gate)),
      " | subtype=", sum(nzchar(subtype_gate)), log_file = log_file
    )
  }
  clinical$olfml2b_expression <- as.numeric(expr_gene[ctx$contract$target_gene, clinical$sample_id])
  clinical$olfml2b_z <- olfml2b_z(clinical$olfml2b_expression)
  clinical$cohort <- accession
  clinical$platform <- platform
  clinical$stage_numeric <- olfml2b_stage_numeric(clinical$stage)
  clinical$recurrence_endpoint <- if (accession == "GSE62254") "DFS" else if (accession == "GSE26253") "RFS" else NA_character_
  clinical$recurrence_time_days <- clinical$dfs_time_days
  clinical$recurrence_event <- clinical$dfs_event
  if (accession == "GSE26253") {
    clinical$rfs_time_days <- clinical$dfs_time_days
    clinical$rfs_event <- clinical$dfs_event
  }
  rownames(clinical) <- clinical$sample_id

  olfml2b_atomic_save_rds(
    expr_gene,
    file.path(ctx$dirs$derived_expression, paste0(accession, "_gene_expression.rds")),
    compress = "xz"
  )
  olfml2b_atomic_write_csv(
    clinical,
    file.path(ctx$dirs$derived_clinical, paste0(accession, "_sample_metadata_harmonized.csv"))
  )

  audit <- olfml2b_geo_cohort_audit(
    accession, platform, expr_probe, expr_gene, clinical,
    mapping, aggregated$target_probes, transformed, excluded$removed
  )
  audit$expected_samples <- row$expected_samples
  audit$sample_count_difference <- audit$n_samples_expression - audit$expected_samples
  audit$clinical_source_candidates <- length(supplementary_candidates)
  audit$clinical_source_best_score <- if (length(supplementary_candidates)) max(vapply(supplementary_candidates, `[[`, numeric(1), "score")) else 0
  audit$public_clinical_files_downloaded <- if (nrow(public_clinical_listing)) sum(public_clinical_listing$status == "downloaded") else 0L
  audit$clinical_merge_fields_filled <- sum(grepl(";", merged_clinical$source_map$source, fixed = TRUE))
  audit$olfml2b_probe_mapping_pass <- length(aggregated$target_probes) >= 1L
  audit$publication_evaluable_os <- audit$n_os_complete >= 30L && audit$os_events >= 20L
  audit$publication_evaluable_dfs <- audit$n_dfs_complete >= 30L && audit$dfs_events >= 20L
  if (identical(accession, "GSE26253")) {
    olfml2b_assert(audit$n_dfs_complete >= 400L,
                 "GSE26253 RFS extraction failed: expected near-complete GEO sample metadata, observed complete RFS n=",
                 audit$n_dfs_complete)
    olfml2b_assert(audit$dfs_events >= 20L,
                 "GSE26253 RFS event extraction failed: observed events=", audit$dfs_events)
  }

  target_mapping <- mapping$target_rows
  target_mapping$cohort <- accession
  target_mapping$platform <- platform
  target_mapping_evidence <- mapping$target_annotations
  if (nrow(target_mapping_evidence)) {
    target_mapping_evidence$cohort <- accession
    target_mapping_evidence$platform <- platform
  }
  target_probe_stats <- aggregated$target_probe_stats
  target_probe_stats$cohort <- accession
  target_probe_stats$platform <- platform

  # Report final merged provenance separately from heuristic pData candidate
  # diagnostics. The latter is useful for debugging but is not the source map
  # used by the analysis.
  source_lookup <- stats::setNames(
    as.character(merged_clinical$source_map$source),
    as.character(merged_clinical$source_map$field)
  )
  final_source <- unname(source_lookup[names(clinical)])
  final_source[is.na(final_source) | !nzchar(final_source)] <- "PIPELINE_DERIVED_OR_ENDPOINT_ALIAS"
  final_mapping_audit <- data.frame(
    cohort = accession,
    concept = names(clinical),
    selected_source = final_source,
    n_nonmissing = vapply(clinical, function(x) {
      if (is.numeric(x) || is.integer(x)) sum(is.finite(as.numeric(x)))
      else sum(!is.na(x) & nzchar(olfml2b_clean_text(x)))
    }, integer(1)),
    n_total = nrow(clinical),
    audit_role = "FINAL_MERGED_CLINICAL_PROVENANCE",
    stringsAsFactors = FALSE
  )
  primary_mapping_diagnostic <- primary_parsed$audit
  if (is.data.frame(primary_mapping_diagnostic) && nrow(primary_mapping_diagnostic)) {
    primary_mapping_diagnostic$audit_role <- "GEO_PDATA_HEURISTIC_CANDIDATE_NOT_FINAL"
    primary_mapping_diagnostic$used_as_final_mapping_audit <- FALSE
  }

  out <- list(
    version = OLFML2B_PART2_PATCH_VERSION,
    pipeline_version = ctx$version,
    release_contract_version = OLFML2B_PART2_PATCH_VERSION,
    generated_at = olfml2b_timestamp(),
    cohort = accession,
    platform = platform,
    role = row$role,
    expression = expr_gene,
    probe_expression = expr_probe,
    sample_metadata = clinical,
    clinical = clinical,
    platform_mapping = mapping$usable,
    platform_mapping_all = mapping$all,
    ambiguous_probes = mapping$ambiguous_probes,
    target_mapping = target_mapping,
    target_mapping_evidence = target_mapping_evidence,
    target_probes = aggregated$target_probes,
    target_probe_expression = aggregated$target_probe_matrix,
    target_probe_stats = target_probe_stats,
    expression_transform = transformed,
    clinical_mapping_audit = final_mapping_audit,
    primary_pdata_mapping_diagnostic = primary_mapping_diagnostic,
    clinical_source_map = merged_clinical$source_map,
    clinical_source_candidates = lapply(merged_clinical$candidates, function(z) list(score = z$score, source = z$source, parser = z$parsed$parser %||% NA_character_)),
    clinical_plausibility = clinical_plausibility,
    expanded_metadata = primary_parsed$expanded,
    supplementary_listing = supplementary_listing,
    removed_samples = excluded$removed,
    identity_leakage_sanitization_audit = identity_sanitized$audit,
    gse15459_survival_sanitization_audit = gse15459_sanitized$audit,
    cohort_audit = audit,
    source_files_manifest = olfml2b_file_manifest(file.path(ctx$dirs$raw_geo, accession))
  )

  olfml2b_atomic_save_rds(out, file.path(ctx$dirs$objects, paste0("Part2_", accession, ".rds")))
  olfml2b_atomic_write_csv(clinical, file.path(ctx$part_paths$tables, paste0("Part2_", accession, "_sample_metadata.csv")))
  olfml2b_atomic_write_csv(audit, file.path(ctx$part_paths$tables, paste0("Part2_", accession, "_cohort_audit.csv")))
  olfml2b_atomic_write_csv(target_mapping, file.path(ctx$part_paths$tables, paste0("Part2_", accession, "_OLFML2B_mapping.csv")))
  olfml2b_atomic_write_csv(target_mapping_evidence, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_OLFML2B_mapping_evidence.csv")))
  olfml2b_atomic_write_csv(target_probe_stats, file.path(ctx$part_paths$tables, paste0("Part2_", accession, "_OLFML2B_probe_stats.csv")))
  olfml2b_atomic_write_csv(final_mapping_audit, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_clinical_mapping_audit.csv")))
  olfml2b_atomic_write_csv(primary_mapping_diagnostic, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_primary_pdata_mapping_diagnostic.csv")))
  olfml2b_atomic_write_csv(merged_clinical$source_map, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_clinical_source_map.csv")))
  olfml2b_atomic_write_csv(clinical_plausibility, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_clinical_plausibility_audit.csv")))
  olfml2b_atomic_write_csv(identity_sanitized$audit, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_identity_leakage_sanitization_audit.csv")))
  if (identical(accession, "GSE15459")) {
    olfml2b_atomic_write_csv(gse15459_sanitized$audit, file.path(ctx$part_paths$reports, "Part2_GSE15459_survival_sanitization_audit.csv"))
  }
  olfml2b_atomic_write_csv(supplementary_listing, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_supplementary_download_audit.csv")))
  olfml2b_atomic_write_csv(excluded$removed, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_removed_samples.csv")))
  olfml2b_atomic_write_csv(out$source_files_manifest, file.path(ctx$part_paths$reports, paste0("Part2_", accession, "_raw_file_manifest.csv")))
  olfml2b_atomic_write_csv(olfml2b_public_data_source_citations(), file.path(ctx$part_paths$reports, "Part2_public_data_source_citations.csv"))

  p <- ggplot2::ggplot(clinical, ggplot2::aes(x = olfml2b_expression)) +
    ggplot2::geom_histogram(bins = 35, fill = "grey70", color = "white") +
    ggplot2::geom_vline(xintercept = median(clinical$olfml2b_expression, na.rm = TRUE), linetype = 2) +
    ggplot2::labs(
      title = paste0(accession, " OLFML2B expression distribution"),
      subtitle = paste0(platform, " | n=", sum(is.finite(clinical$olfml2b_expression)),
                        " | mean=", sprintf("%.3f", mean(clinical$olfml2b_expression, na.rm = TRUE)),
                        " | SD=", sprintf("%.3f", stats::sd(clinical$olfml2b_expression, na.rm = TRUE)),
                        " | probes: ", paste(aggregated$target_probes, collapse = ", ")),
      x = "OLFML2B expression", y = "Samples"
    ) + olfml2b_base_theme()
  olfml2b_save_plot(p, file.path(ctx$part_paths$figures, paste0("Part2_", accession, "_OLFML2B_distribution.png")), 7, 5)

  if (sum(!is.na(clinical$stage)) >= 20L) {
    stage_keep <- !is.na(clinical$stage) & is.finite(clinical$olfml2b_expression)
    stage_num <- as.numeric(clinical$stage[stage_keep])
    stage_cor <- suppressWarnings(stats::cor.test(
      clinical$olfml2b_expression[stage_keep], stage_num,
      method = "spearman", exact = FALSE
    ))
    stage_kw <- stats::kruskal.test(clinical$olfml2b_expression[stage_keep] ~ clinical$stage[stage_keep])
    ps <- ggplot2::ggplot(clinical[!is.na(clinical$stage), ], ggplot2::aes(x = stage, y = olfml2b_expression)) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "white") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
      ggplot2::labs(
        title = paste0(accession, " OLFML2B by stage"),
        subtitle = paste0("n=", sum(stage_keep), " | Spearman rho=", sprintf("%.3f", unname(stage_cor$estimate)),
                          ", p=", format.pval(stage_cor$p.value, digits = 3),
                          " | Kruskal-Wallis p=", format.pval(stage_kw$p.value, digits = 3)),
        x = NULL, y = "OLFML2B expression"
      ) +
      olfml2b_base_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    olfml2b_save_plot(ps, file.path(ctx$part_paths$figures, paste0("Part2_", accession, "_OLFML2B_stage.png")), 7, 5)
  }
  out
}

run_part2 <- function(ctx = NULL) {
  if (is.null(ctx)) {
    code_root <- olfml2b_find_code_root()
    ctx <- olfml2b_load_context()
  }
  ctx$part_paths <- olfml2b_part_paths(ctx, "Part2")
  olfml2b_require_packages(c(
    "GEOquery", "Biobase", "AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db",
    "matrixStats", "data.table", "readxl", "R.utils", "ggplot2"
  ))
  log_file <- file.path(ctx$part_paths$logs, "Part2_GEO_Bulk_Acquire_Process.log")
  error_file <- file.path(ctx$part_paths$reports, "pipeline_errors.csv")
  out_rds <- file.path(ctx$dirs$objects, "Part2_GEO_bulk_index.rds")
  olfml2b_log("INFO", "Starting Part2 | ", ctx$version, log_file = log_file)
  cohorts <- ctx$contract$geo_cohorts
  results <- list()
  for (accession in cohorts) {
    results[[accession]] <- olfml2b_run_stage(
      paste0("Process ", accession),
      olfml2b_process_one_geo_bulk(accession, ctx, log_file, error_file),
      log_file, error_file, required = TRUE
    )
  }
  audits <- do.call(rbind, lapply(results, `[[`, "cohort_audit"))
  mapping_audits <- do.call(rbind, lapply(results, `[[`, "target_mapping"))
  if (nrow(mapping_audits)) {
    nonempty <- apply(mapping_audits, 1L, function(z) any(!is.na(z) & nzchar(olfml2b_clean_text(z))))
    mapping_audits <- mapping_audits[nonempty, , drop = FALSE]
    olfml2b_assert(!any(is.na(mapping_audits$probe_id) | !nzchar(olfml2b_clean_text(mapping_audits$probe_id))),
                 "Part2 target mapping audit contains empty probe identifiers")
    olfml2b_assert(nrow(mapping_audits) <= 100L,
                 "Part2 target mapping audit unexpectedly expanded to ", nrow(mapping_audits), " rows")
    olfml2b_assert_target_mapping_clean(mapping_audits, target_gene = ctx$contract$target_gene)
  }
  endpoint_audit <- audits[, intersect(c(
    "cohort", "platform", "n_samples_expression", "recurrence_endpoint",
    "n_os_complete", "os_events", "n_recurrence_complete", "recurrence_events",
    "publication_evaluable_os", "publication_evaluable_dfs"
  ), names(audits)), drop = FALSE]
  index <- list(
    version = ctx$version,
    generated_at = olfml2b_timestamp(),
    cohort_files = setNames(file.path(ctx$dirs$objects, paste0("Part2_", cohorts, ".rds")), cohorts),
    cohort_audit = audits,
    endpoint_audit = endpoint_audit,
    target_mapping_audit = mapping_audits
  )
  olfml2b_atomic_save_rds(index, out_rds)
  olfml2b_atomic_write_csv(audits, file.path(ctx$part_paths$tables, "Part2_GEO_bulk_cohort_audit.csv"))
  olfml2b_atomic_write_csv(endpoint_audit, file.path(ctx$part_paths$tables, "Part2_GEO_endpoint_audit.csv"))
  olfml2b_atomic_write_csv(mapping_audits, file.path(ctx$part_paths$tables, "Part2_GEO_OLFML2B_mapping_audit.csv"))
  olfml2b_capture_session(file.path(ctx$part_paths$reports, "Part2_sessionInfo.txt"))
  olfml2b_log("INFO", "Part2 complete: ", out_rds, log_file = log_file)
  invisible(index)
}


# ==============================================================================
# GSE84437 integrated GPL6947 builder used by Part2.
# ==============================================================================
############################################################
## Integrated GSE84437 Part2-compatible builder
## Version: v1.0.1_20260705_GPL6947_ANNOTATION_PACKAGE_AND_RETRY_FIX
##
## Purpose:
## - Build D:/OLFML2B_STAD/output/objects/Part2_GSE84437.rds
## - Update D:/OLFML2B_STAD/output/objects/Part2_GEO_bulk_index.rds
## - Use local Series Matrix first.
## - Parse OS from GSE84437 phenotype fields:
##     death
##     duration overall survival
## - Aggregate probes to gene symbols using fData/GPL annotation.
##
## Run:
##   setwd("D:/OLFML2B_STAD")
##   source("R/02_OLFML2B_PART2_GEO.R")
##   gse84437_object <- olfml2b_build_gse84437_part2_object(
##       root = "D:/OLFML2B_STAD",
##       allow_gpl_download = TRUE
##   )
############################################################

options(stringsAsFactors = FALSE)

GSE84437_PART2_VERSION <- "v2.0.0_20260720_PROVENANCE_PARTITIONED_FORMAL_OS"

gse84437_partition_samples <- function(sample_ids) {
    sample_ids <- toupper(trimws(as.character(sample_ids)))
    gsm_number <- suppressWarnings(as.integer(sub("^GSM", "", sample_ids)))
    source_subseries <- rep("UNKNOWN", length(sample_ids))
    source_subseries[grepl("^GSM223", sample_ids) & gsm_number <= 2235637L] <- "GSE84426"
    source_subseries[grepl("^GSM223", sample_ids) & gsm_number >= 2235695L & gsm_number <= 2236095L] <- "GSE84433"
    source_subseries[grepl("^GSM441", sample_ids) & gsm_number >= 4419484L & gsm_number <= 4419533L] <- "GSE147163"
    analysis_role <- ifelse(source_subseries %in% c("GSE84426", "GSE84433"),
                            "FORMAL_OS", ifelse(source_subseries == "GSE147163", "CONTEXT_ONLY", "UNKNOWN"))
    processing_batch <- ifelse(
        analysis_role == "FORMAL_OS", "2016_GenomeStudio_quantile_normalized",
        ifelse(analysis_role == "CONTEXT_ONLY", "2020_lumi_neqc", "UNKNOWN")
    )
    out <- data.frame(sample_id = sample_ids, source_subseries = source_subseries,
                      analysis_role = analysis_role, processing_batch = processing_batch,
                      stringsAsFactors = FALSE)
    if (any(out$analysis_role == "UNKNOWN")) {
        stop("GSE84437 provenance contract found unknown GSM accessions: ",
             paste(head(out$sample_id[out$analysis_role == "UNKNOWN"], 20L), collapse = ";"), call. = FALSE)
    }
    counts <- table(out$source_subseries)
    expected <- c(GSE84426 = 76L, GSE84433 = 357L, GSE147163 = 50L)
    observed <- as.integer(counts[names(expected)]); observed[is.na(observed)] <- 0L
    if (!identical(unname(observed), unname(expected))) {
        stop("GSE84437 provenance contract failed. Observed ",
             paste(names(expected), observed, sep = "=", collapse = ";"),
             "; expected GSE84426=76;GSE84433=357;GSE147163=50.", call. = FALSE)
    }
    out
}

olfml2b_build_gse84437_part2_object <- function(
    root = "D:/OLFML2B_STAD",
    accession = "GSE84437",
    allow_gpl_download = TRUE,
    allow_bioc_annotation_install = FALSE,
    download_timeout_sec = 600L,
    overwrite = TRUE
) {
    msg("Version: ", GSE84437_PART2_VERSION)
    msg("Root: ", root)

    require_pkg("GEOquery")
    require_pkg("Biobase")
    require_pkg("matrixStats")
    require_pkg("ggplot2")

    suppressPackageStartupMessages({
        library(GEOquery)
        library(Biobase)
        library(matrixStats)
    })

    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    old_timeout <- getOption("timeout")
    options(timeout = max(as.integer(download_timeout_sec), old_timeout %||% 60L, 600L))
    on.exit(options(timeout = old_timeout), add = TRUE)
    raw_dir <- file.path(root, "data", "raw", "GEO", accession)
    gpl_dir <- dir_create(file.path(root, "data", "raw", "GEO", "GPL"))
    object_dir <- dir_create(file.path(root, "output", "objects"))
    table_dir <- dir_create(file.path(root, "output", "tables", "Part2"))
    report_dir <- dir_create(file.path(root, "output", "reports", "Part2"))
    derived_expression_dir <- dir_create(file.path(root, "data", "derived", "expression"))
    derived_clinical_dir <- dir_create(file.path(root, "data", "derived", "clinical"))

    sm_files <- list.files(
        raw_dir,
        pattern = "series_matrix.*\\.txt(\\.gz)?$|series_matrix.*\\.gz$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )
    if (length(sm_files) == 0L) {
        stop("No local Series Matrix found under: ", raw_dir, call. = FALSE)
    }
    sm_file <- sm_files[order(nchar(sm_files), sm_files)][1]
    msg("Using Series Matrix: ", sm_file)

    obj <- GEOquery::getGEO(filename = sm_file, getGPL = FALSE)
    eset <- pick_eset(obj)
    platform <- tryCatch(Biobase::annotation(eset), error = function(e) NA_character_)
    if (is.na(platform) || !nzchar(platform)) {
        platform <- infer_platform_from_pdata(Biobase::pData(eset))
    }
    msg("Platform: ", platform)

    expr_probe_all <- as.matrix(Biobase::exprs(eset))
    storage.mode(expr_probe_all) <- "numeric"
    partition <- gse84437_partition_samples(colnames(expr_probe_all))
    formal_ids <- partition$sample_id[partition$analysis_role == "FORMAL_OS"]
    context_ids <- partition$sample_id[partition$analysis_role == "CONTEXT_ONLY"]
    # Transform each known processing batch independently. The 2020 neqc
    # samples must never determine the transform or SD used by the 2016 OS set.
    transformed <- transform_expression_if_needed(expr_probe_all[, formal_ids, drop = FALSE])
    expr_probe <- transformed$matrix
    context_transformed <- transform_expression_if_needed(expr_probe_all[, context_ids, drop = FALSE])
    expr_probe_context <- context_transformed$matrix

    pdata_all <- Biobase::pData(eset)
    pdata_ids <- if ("geo_accession" %in% names(pdata_all)) as.character(pdata_all$geo_accession) else rownames(pdata_all)
    pdata_raw <- pdata_all[match(formal_ids, pdata_ids), , drop = FALSE]
    pdata_context <- pdata_all[match(context_ids, pdata_ids), , drop = FALSE]
    if (anyNA(match(formal_ids, pdata_ids)) || anyNA(match(context_ids, pdata_ids))) {
        stop("GSE84437 pData could not be aligned to provenance-partitioned expression columns.", call. = FALSE)
    }
    clinical_expanded <- expand_characteristics(pdata_raw)
    clinical <- parse_gse84437_clinical(clinical_expanded)
    context_expanded <- expand_characteristics(pdata_context)
    context_clinical <- parse_gse84437_clinical(context_expanded)

    mapping <- build_platform_mapping(
        fdat = tryCatch(Biobase::fData(eset), error = function(e) data.frame()),
        platform = platform,
        probe_ids = rownames(expr_probe),
        gpl_dir = gpl_dir,
        cohort_dir = raw_dir,
        allow_gpl_download = allow_gpl_download,
        allow_bioc_annotation_install = allow_bioc_annotation_install,
        download_timeout_sec = download_timeout_sec
    )

    target_info <- gse84437_extract_olfml2b_target(
        expr_probe = expr_probe,
        mapping = mapping,
        platform = platform,
        target_gene = "OLFML2B",
        target_entrez = olfml2b_target_entrez_id()
    )

    expr_gene <- aggregate_probes_to_genes(expr_probe, mapping$usable)
    expr_gene_context <- aggregate_probes_to_genes(expr_probe_context, mapping$usable)

    ordered <- normalize_sample_order(expr_gene, clinical)
    expr_gene <- ordered$expr
    clinical <- ordered$clinical
    context_ordered <- normalize_sample_order(expr_gene_context, context_clinical)
    expr_gene_context <- context_ordered$expr
    context_clinical <- context_ordered$clinical

    clinical$cohort <- accession
    clinical$platform <- platform
    formal_provenance <- partition[match(clinical$sample_id, partition$sample_id), , drop = FALSE]
    clinical$source_subseries <- formal_provenance$source_subseries
    clinical$analysis_role <- formal_provenance$analysis_role
    clinical$processing_batch <- formal_provenance$processing_batch
    clinical$clinical_schema <- "GSE84426_GSE84433_FORMAL_OS"
    if (!"stage_numeric" %in% names(clinical)) {
        clinical$stage_numeric <- if ("stage" %in% names(clinical)) stage_numeric(clinical$stage) else NA_real_
    }
    if ("OLFML2B" %in% rownames(expr_gene)) {
        clinical$olfml2b_expression <- as.numeric(expr_gene["OLFML2B", clinical$sample_id])
        clinical$olfml2b_z_within_subseries <- ave(
            clinical$olfml2b_expression, clinical$source_subseries,
            FUN = function(z) olfml2b_z(as.numeric(z))
        )
        # A final common rescaling gives one pooled SD while retaining a zero
        # mean within each formal subseries; no GSE84426/GSE84433 mean shift is
        # allowed to contribute to the exposure contrast.
        clinical$olfml2b_z <- olfml2b_z(clinical$olfml2b_z_within_subseries)
    } else {
        stop("GSE84437 gene-level matrix does not contain OLFML2B after aggregation.", call. = FALSE)
    }

    context_clinical$cohort <- "GSE147163"
    context_clinical$platform <- platform
    context_provenance <- partition[match(context_clinical$sample_id, partition$sample_id), , drop = FALSE]
    context_clinical$source_subseries <- context_provenance$source_subseries
    context_clinical$analysis_role <- context_provenance$analysis_role
    context_clinical$processing_batch <- context_provenance$processing_batch
    context_clinical$clinical_schema <- "CONTEXT_ONLY_NO_FORMAL_SURVIVAL"
    context_clinical$olfml2b_expression <- as.numeric(expr_gene_context["OLFML2B", context_clinical$sample_id])
    context_clinical$olfml2b_z <- olfml2b_z(context_clinical$olfml2b_expression)
    context_clinical$os_time_days <- NA_real_
    context_clinical$os_event <- NA_integer_

    if (nrow(clinical) != 433L || ncol(expr_gene) != 433L ||
        any(clinical$analysis_role != "FORMAL_OS") || any(grepl("^GSM441", clinical$sample_id))) {
        stop("GSE84437 formal OS contract failed after alignment; expected exactly 433 2016 formal samples and no GSE147163 rows.", call. = FALSE)
    }
    sex_tab <- table(as.character(clinical$sex), useNA = "no")
    if (!all(c("Female", "Male") %in% names(sex_tab)) || any(sex_tab[c("Female", "Male")] == 0L)) {
        stop("GSE84437 sex parsing contract failed: both Female and Male must be retained.", call. = FALSE)
    }
    formal_z_sd <- stats::sd(clinical$olfml2b_z, na.rm = TRUE)
    formal_batch_means <- tapply(clinical$olfml2b_z, clinical$source_subseries, mean, na.rm = TRUE)
    if (!is.finite(formal_z_sd) || abs(formal_z_sd - 1) > 1e-8 ||
        any(abs(formal_batch_means) > 1e-8)) {
        stop("GSE84437 within-formal-cohort z-score contract failed: SD=", formal_z_sd, call. = FALSE)
    }

    rownames(clinical) <- clinical$sample_id

    audit <- make_gse84437_audit(
        accession = accession,
        platform = platform,
        expr_probe = expr_probe,
        expr_gene = expr_gene,
        clinical = clinical,
        mapping = mapping,
        transformed = transformed,
        sm_file = sm_file
    )
    if (audit$n_os_complete[1] != 431L || audit$os_events[1] != 207L) {
        stop("GSE84437 formal endpoint contract failed: complete OS n=", audit$n_os_complete[1],
             ", events=", audit$os_events[1], "; expected 431 and 207.", call. = FALSE)
    }
    audit$n_superseries_total <- ncol(expr_probe_all)
    audit$n_formal_os <- ncol(expr_probe)
    audit$n_context_only_excluded <- ncol(expr_probe_context)
    audit$formal_z_sd <- formal_z_sd
    audit$max_abs_formal_subseries_z_mean <- max(abs(formal_batch_means))
    audit$batch_partition_status <- "PASS_PROVENANCE_PARTITIONED_BEFORE_TRANSFORM"

    target_probe <- target_info$target_probes
    target_long <- data.frame(
        sample_id = c(colnames(expr_probe), colnames(expr_probe_context)),
        olfml2b_probe_expression = c(
            colMeans(expr_probe[target_probe, , drop = FALSE], na.rm = TRUE),
            colMeans(expr_probe_context[target_probe, , drop = FALSE], na.rm = TRUE)
        ),
        stringsAsFactors = FALSE
    )
    target_long <- merge(target_long, partition, by = "sample_id", all.x = TRUE, sort = FALSE)
    batch_audit <- do.call(rbind, lapply(split(target_long, target_long$source_subseries), function(z) {
        data.frame(
            source_subseries = z$source_subseries[1], analysis_role = z$analysis_role[1],
            processing_batch = z$processing_batch[1], n = nrow(z),
            mean_olfml2b = mean(z$olfml2b_probe_expression, na.rm = TRUE),
            sd_olfml2b = stats::sd(z$olfml2b_probe_expression, na.rm = TRUE),
            median_olfml2b = stats::median(z$olfml2b_probe_expression, na.rm = TRUE),
            transform_scope = ifelse(z$analysis_role[1] == "FORMAL_OS", "formal_2016_batch_only", "context_2020_batch_only"),
            exposure_standardization = ifelse(z$analysis_role[1] == "FORMAL_OS",
                                               "within_GSE84426_or_GSE84433_then_pooled_1SD",
                                               "within_GSE147163_context_only"),
            eligible_for_formal_survival = z$analysis_role[1] == "FORMAL_OS",
            stringsAsFactors = FALSE
        )
    }))
    rownames(batch_audit) <- NULL

    out <- list(
        version = GSE84437_PART2_VERSION,
        release_contract_version = OLFML2B_PART2_PATCH_VERSION,
        generated_at = timestamp(),
        cohort = accession,
        platform = platform,
        role = "external_validation_os_formal_2016_only",
        expression = expr_gene,
        probe_expression = expr_probe,
        sample_metadata = clinical,
        platform_mapping = mapping$usable,
        platform_mapping_all = mapping$all,
        ambiguous_probes = mapping$ambiguous_probes,
        target_mapping = target_info$target_mapping,
        target_mapping_evidence = target_info$target_mapping_evidence,
        target_probes = target_info$target_probes,
        target_probe_expression = target_info$target_probe_expression,
        target_probe_stats = target_info$target_probe_stats,
        expression_transform = transformed,
        clinical_mapping_audit = attr(clinical, "clinical_mapping_audit"),
        clinical_source_map = data.frame(
            field = c("os_time_days", "os_event"),
            source = c(attr(clinical, "os_time_source"), attr(clinical, "os_event_source")),
            stringsAsFactors = FALSE
        ),
        clinical_plausibility = audit[, c(
            "cohort", "platform", "n_samples_expression", "n_samples_clinical",
            "n_os_complete", "os_events", "n_stage", "n_age"
        ), drop = FALSE],
        expanded_metadata = clinical_expanded,
        superseries_partition_audit = partition,
        batch_audit = batch_audit,
        supplementary_listing = data.frame(),
        removed_samples = transform(
            partition[partition$analysis_role == "CONTEXT_ONLY", , drop = FALSE],
            reason = "GSE147163_2020_neqc_context_only_excluded_from_formal_OS"
        ),
        identity_leakage_sanitization_audit = data.frame(),
        gse15459_survival_sanitization_audit = data.frame(),
        cohort_audit = audit,
        source_files_manifest = file_manifest(raw_dir)
    )

    context_out <- list(
        version = GSE84437_PART2_VERSION,
        release_contract_version = OLFML2B_PART2_PATCH_VERSION,
        generated_at = timestamp(), cohort = "GSE147163", parent_superseries = accession,
        platform = platform, role = "context_only_not_for_survival_or_meta_analysis",
        expression = expr_gene_context, probe_expression = expr_probe_context,
        sample_metadata = context_clinical, expression_transform = context_transformed,
        target_mapping = target_info$target_mapping,
        provenance = partition[partition$analysis_role == "CONTEXT_ONLY", , drop = FALSE],
        claim_ceiling = "molecular_context_only; not a formal prognostic validation cohort"
    )

    rds_path <- file.path(object_dir, paste0("Part2_", accession, ".rds"))
    saveRDS(out, rds_path, compress = "xz")
    context_path <- file.path(object_dir, "Part2_GSE147163_context_only.rds")
    saveRDS(context_out, context_path, compress = "xz")
    saveRDS(expr_gene, file.path(derived_expression_dir, paste0(accession, "_gene_expression.rds")), compress = "xz")
    write_csv(clinical, file.path(derived_clinical_dir, paste0(accession, "_sample_metadata_harmonized.csv")))

    write_csv(audit, file.path(table_dir, paste0("Part2_", accession, "_cohort_audit.csv")))
    write_csv(clinical, file.path(table_dir, paste0("Part2_", accession, "_sample_metadata.csv")))
    write_csv(mapping$usable, file.path(table_dir, paste0("Part2_", accession, "_platform_mapping_usable.csv")))
    write_csv(mapping$all, file.path(table_dir, paste0("Part2_", accession, "_platform_mapping_all.csv")))
    write_csv(target_info$target_mapping, file.path(table_dir, paste0("Part2_", accession, "_OLFML2B_mapping.csv")))
    write_csv(target_info$target_probe_stats, file.path(table_dir, paste0("Part2_", accession, "_OLFML2B_probe_stats.csv")))
    write_csv(target_info$target_mapping_evidence, file.path(report_dir, paste0("Part2_", accession, "_OLFML2B_mapping_evidence.csv")))
    write_csv(attr(clinical, "clinical_mapping_audit"), file.path(report_dir, paste0("Part2_", accession, "_clinical_mapping_audit.csv")))
    write_csv(partition, file.path(report_dir, "Part2_GSE84437_superseries_partition_audit.csv"))
    write_csv(batch_audit, file.path(table_dir, "Part2_GSE84437_processing_batch_audit.csv"))
    write_csv(out$removed_samples, file.path(report_dir, "Part2_GSE84437_context_only_removed_from_formal_OS.csv"))
    write_csv(context_clinical, file.path(table_dir, "Part2_GSE147163_context_only_sample_metadata.csv"))
    write_csv(file_manifest(raw_dir), file.path(report_dir, paste0("Part2_", accession, "_raw_file_manifest.csv")))

    p_dist <- ggplot2::ggplot(clinical, ggplot2::aes(x = olfml2b_z)) +
        ggplot2::geom_histogram(bins = 35, fill = "grey70", color = "white") +
        ggplot2::geom_vline(xintercept = stats::median(clinical$olfml2b_z, na.rm = TRUE), linetype = 2) +
        ggplot2::labs(
            title = "GSE84437 formal batch-aware OLFML2B distribution",
            subtitle = "GSE84426/GSE84433 centered and scaled separately, then pooled to one SD | n=433",
            x = "OLFML2B within-subseries standardized expression (pooled SD)", y = "Samples"
        ) + olfml2b_base_theme()
    olfml2b_save_plot(p_dist, file.path(root, "output", "figures", "Part2", "Part2_GSE84437_OLFML2B_distribution.png"), 7, 5)

    p_batch <- ggplot2::ggplot(target_long, ggplot2::aes(x = source_subseries, y = olfml2b_probe_expression, fill = analysis_role)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 0.7) +
        ggplot2::labs(
            title = "GSE84437 provenance and processing-batch QC",
            subtitle = "Batches transformed independently; GSE147163 is context-only and excluded from formal OS inference",
            x = NULL, y = "OLFML2B probe expression"
        ) + olfml2b_base_theme()
    olfml2b_save_plot(p_batch, file.path(root, "output", "figures", "Part2", "Part2_GSE84437_processing_batch_QC.png"), 7.5, 5)

    index_path <- file.path(object_dir, "Part2_GEO_bulk_index.rds")
    index <- update_part2_geo_index(
        index_path = index_path,
        accession = accession,
        object_path = rds_path,
        audit = audit,
        target_mapping = target_info$target_mapping
    )

    msg("Saved: ", rds_path)
    msg("Updated index: ", index_path)
    msg("OS complete n=", audit$n_os_complete[1], " | events=", audit$os_events[1])
    msg("Gene-level matrix: genes=", nrow(expr_gene), " | samples=", ncol(expr_gene))

    invisible(list(
        object = out,
        object_path = rds_path,
        index = index,
        audit = audit,
        table_dir = table_dir,
        report_dir = report_dir
    ))
}

pick_eset <- function(obj) {
    if (inherits(obj, "ExpressionSet")) return(obj)
    if (is.list(obj)) {
        esets <- obj[vapply(obj, inherits, logical(1), what = "ExpressionSet")]
        if (!length(esets)) stop("GEOquery returned no ExpressionSet.", call. = FALSE)
        sizes <- vapply(esets, function(e) ncol(Biobase::exprs(e)), integer(1))
        return(esets[[which.max(sizes)]])
    }
    stop("Unsupported GEO object class: ", paste(class(obj), collapse = ";"), call. = FALSE)
}

infer_platform_from_pdata <- function(pdata) {
    if (is.null(pdata) || !nrow(pdata)) return(NA_character_)
    vals <- unlist(pdata, use.names = FALSE)
    hit <- regmatches(vals, regexpr("GPL[0-9]+", vals, ignore.case = TRUE))
    hit <- hit[nzchar(hit)]
    if (length(hit)) toupper(hit[1]) else NA_character_
}

parse_gse84437_clinical <- function(pdata) {
    clinical <- as.data.frame(pdata, stringsAsFactors = FALSE, check.names = FALSE)
    original_names <- names(clinical)
    names(clinical) <- make.unique(clean_names(names(clinical)))

    if (!"sample_id" %in% names(clinical)) {
        if ("geo_accession" %in% names(clinical)) {
            clinical$sample_id <- as.character(clinical$geo_accession)
        } else {
            clinical$sample_id <- rownames(clinical)
        }
    }
    clinical$sample_id <- as.character(clinical$sample_id)

    if (!"patient_id" %in% names(clinical)) {
        patient_col <- find_first_col(clinical, c("^patient$", "patient_id", "case", "subject", "title"))
        clinical$patient_id <- if (!is.na(patient_col)) as.character(clinical[[patient_col]]) else clinical$sample_id
    }

    time_col <- find_first_col(
        clinical,
        c(
            "duration.*overall.*survival",
            "overall.*survival.*duration",
            "overall.*survival.*time",
            "^os.*time",
            "survival.*time",
            "follow.*up"
        )
    )
    event_col <- find_first_col(
        clinical,
        c(
            "^death$",
            "death",
            "dead",
            "vital",
            "overall.*survival.*status",
            "^os.*event"
        )
    )

    if (is.na(time_col) || is.na(event_col)) {
        audit <- data.frame(
            field = c("os_time_days", "os_event"),
            source_column = c(time_col, event_col),
            status = "FAILED_TO_RESOLVE",
            stringsAsFactors = FALSE
        )
        attr(clinical, "clinical_mapping_audit") <- audit
        attr(clinical, "os_time_source") <- time_col
        attr(clinical, "os_event_source") <- event_col
        clinical$os_time_days <- NA_real_
        clinical$os_event <- NA_integer_
    } else {
        raw_time <- clinical[[time_col]]
        raw_event <- clinical[[event_col]]
        clinical$os_time_days <- infer_survival_time_days(raw_time, time_col)
        clinical$os_event <- as_event(raw_event)
        audit <- data.frame(
            field = c("os_time_days", "os_event"),
            source_column = c(time_col, event_col),
            n_nonmissing = c(sum(is.finite(clinical$os_time_days)), sum(clinical$os_event %in% c(0L, 1L), na.rm = TRUE)),
            example_values = c(
                paste(head(unique(na.omit(as.character(raw_time))), 8), collapse = " || "),
                paste(head(unique(na.omit(as.character(raw_event))), 8), collapse = " || ")
            ),
            status = "OK",
            stringsAsFactors = FALSE
        )
        attr(clinical, "clinical_mapping_audit") <- audit
        attr(clinical, "os_time_source") <- time_col
        attr(clinical, "os_event_source") <- event_col
    }

    age_col <- find_first_col(clinical, c("^age$", "age_at", "age_year"))
    if (!is.na(age_col) && !"age" %in% names(clinical)) {
        clinical$age <- suppressWarnings(as.numeric(clinical[[age_col]]))
        clinical$age[is.finite(clinical$age) & (clinical$age < 18 | clinical$age > 120)] <- NA_real_
    } else if ("age" %in% names(clinical)) {
        clinical$age <- suppressWarnings(as.numeric(clinical$age))
        clinical$age[is.finite(clinical$age) & (clinical$age < 18 | clinical$age > 120)] <- NA_real_
    } else {
        clinical$age <- NA_real_
    }

    sex_col <- find_first_col(clinical, c("^sex$", "^gender$"))
    if (!is.na(sex_col) && !"sex" %in% names(clinical)) clinical$sex <- sex_factor(clinical[[sex_col]])
    if ("sex" %in% names(clinical)) clinical$sex <- sex_factor(clinical$sex)

    stage_col <- find_first_col(clinical, c("^stage$", "pathologic.*stage", "pathological.*stage", "ajcc.*stage", "tumor.*stage"))
    if (!is.na(stage_col) && !"stage" %in% names(clinical)) clinical$stage <- stage_factor(clinical[[stage_col]])
    if ("stage" %in% names(clinical)) clinical$stage <- stage_factor(clinical$stage)

    clinical$dfs_time_days <- NA_real_
    clinical$dfs_event <- NA_integer_
    clinical$rfs_time_days <- NA_real_
    clinical$rfs_event <- NA_integer_
    clinical$recurrence_time_days <- NA_real_
    clinical$recurrence_event <- NA_integer_

    attr(clinical, "original_names") <- original_names
    clinical
}

find_first_col <- function(df, patterns) {
    nms <- names(df)
    for (p in patterns) {
        hit <- grep(p, nms, ignore.case = TRUE, perl = TRUE, value = TRUE)
        if (length(hit)) return(hit[1])
    }
    NA_character_
}

infer_survival_time_days <- function(x, source_name) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    source_low <- tolower(source_name)
    finite <- y[is.finite(y)]
    if (!length(finite)) return(y)

    ## Units are resolved only from explicit labels or the frozen GSE84437
    ## metadata dictionary. Numeric ranges are never used to guess units.
    if (grepl("month", source_low)) return(y * 30.4375)
    if (grepl("year", source_low)) return(y * 365.25)
    if (grepl("day", source_low)) return(y)
    if (grepl("duration.*overall.*survival|overall.*survival.*duration", source_low)) return(y * 30.4375)
    rep(NA_real_, length(y))
}

as_event <- function(x) {
    if (is.null(x)) return(integer())
    if (is.numeric(x) || is.integer(x)) {
        y <- suppressWarnings(as.numeric(x))
        out <- rep(NA_integer_, length(y))
        out[is.finite(y) & y == 1] <- 1L
        out[is.finite(y) & y == 0] <- 0L
        return(out)
    }
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(y))
    out[y %in% c("1", "yes", "true", "dead", "death", "deceased", "died", "event")] <- 1L
    out[y %in% c("0", "no", "false", "alive", "living", "censored", "censor")] <- 0L
    out
}

build_platform_mapping <- function(
    fdat,
    platform,
    probe_ids,
    gpl_dir,
    cohort_dir,
    allow_gpl_download = TRUE,
    allow_bioc_annotation_install = FALSE,
    download_timeout_sec = 600L
) {
    maps <- list()

    ## 1) Try fData from Series Matrix first.
    if (!is.null(fdat) && is.data.frame(fdat) && nrow(fdat) > 0L) {
        maps[["fdata"]] <- map_probe_table(fdat, probe_ids, source_prefix = "fData")
    }

    ## 2) Try Bioconductor Illumina annotation package.
    ## GPL6947 is Illumina HumanHT-12 V3; package illuminaHumanv3.db is often
    ## more reliable than downloading the large NCBI GPL SOFT file.
    pkg_map <- map_from_bioc_annotation_package(
        platform = platform,
        probe_ids = probe_ids,
        allow_install = allow_bioc_annotation_install
    )
    if (nrow(pkg_map) > 0L) {
        maps[["bioc_annotation_package"]] <- pkg_map
    }

    ## 3) Try local GPL SOFT/annot table, then robust re-download if needed.
    gpl_file <- find_local_gpl_file(platform, gpl_dir, cohort_dir)
    if (!is.na(gpl_file) && file.exists(gpl_file)) {
        msg("Trying local GPL annotation: ", gpl_file)
        gpl_tab <- read_gpl_table(gpl_file)
        if (!is.null(gpl_tab) && nrow(gpl_tab) > 0L) {
            gpl_map <- map_probe_table(gpl_tab, probe_ids, source_prefix = basename(gpl_file))
            if (nrow(gpl_map) > 0L) {
                maps[["local_gpl"]] <- gpl_map
            } else {
                msg("Local GPL table was readable but produced no symbol mapping; ignoring it.", level = "WARN")
            }
        } else {
            msg("Local GPL file was not a valid readable platform table; ignoring it.", level = "WARN")
        }
    }

    if (isTRUE(allow_gpl_download) && !is.na(platform) && nzchar(platform)) {
        current_map_n <- nrow(bind_rows(maps))
        if (current_map_n < 1000L) {
            msg("Usable mapping is insufficient; robust-downloading platform annotation for ", platform)
            gpl_file2 <- try_download_gpl(
                platform = platform,
                gpl_dir = gpl_dir,
                timeout_sec = download_timeout_sec
            )
            if (!is.na(gpl_file2) && file.exists(gpl_file2)) {
                msg("Using downloaded GPL annotation: ", gpl_file2)
                gpl_tab2 <- read_gpl_table(gpl_file2)
                if (!is.null(gpl_tab2) && nrow(gpl_tab2) > 0L) {
                    gpl_map2 <- map_probe_table(gpl_tab2, probe_ids, source_prefix = basename(gpl_file2))
                    if (nrow(gpl_map2) > 0L) maps[["downloaded_gpl"]] <- gpl_map2
                }
            }
        }
    }

    all <- bind_rows(maps)
    if (nrow(all) == 0L) {
        usable <- data.frame(probe_id = character(), gene_symbol = character(), source = character())
        ambiguous <- data.frame()
    } else {
        all$probe_id <- as.character(all$probe_id)
        all$gene_symbol <- toupper(trimws(as.character(all$gene_symbol)))
        all <- all[nzchar(all$gene_symbol) & all$probe_id %in% probe_ids, , drop = FALSE]
        all <- all[!duplicated(paste(all$probe_id, all$gene_symbol, all$source)), , drop = FALSE]

        ## If multiple sources agree on a symbol, keep it.  Only discard probes
        ## that map to multiple distinct symbols.
        symbol_count <- tapply(all$gene_symbol, all$probe_id, function(z) length(unique(z)))
        ambiguous_ids <- names(symbol_count)[symbol_count > 1L]
        ambiguous <- all[all$probe_id %in% ambiguous_ids, , drop = FALSE]
        usable <- all[!all$probe_id %in% ambiguous_ids, , drop = FALSE]
        usable <- usable[!duplicated(paste(usable$probe_id, usable$gene_symbol)), , drop = FALSE]
    }

    if (nrow(usable) == 0L) {
        stop(
            "No usable platform mapping for ", platform, ". ",
            "Restore the frozen annotation package illuminaHumanv3.db before starting the run, or place GPL6947_family.soft.gz ",
            "under D:/OLFML2B_STAD/data/raw/GEO/GPL/.",
            call. = FALSE
        )
    }

    msg("Platform mapping usable rows: ", nrow(usable), " | unique genes: ", length(unique(usable$gene_symbol)))
    list(all = all, usable = usable, ambiguous_probes = ambiguous)
}

map_from_bioc_annotation_package <- function(platform, probe_ids, allow_install = FALSE) {
    platform <- toupper(as.character(platform %||% ""))
    pkg <- NA_character_
    if (platform == "GPL6947") pkg <- "illuminaHumanv3.db"
    if (platform %in% c("GPL10558", "GPL6884")) pkg <- "illuminaHumanv4.db"
    if (is.na(pkg)) return(data.frame())

    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
        msg("AnnotationDbi is not installed; skipping Bioconductor annotation package mapping.", level = "WARN")
        return(data.frame())
    }

    if (!requireNamespace(pkg, quietly = TRUE) && isTRUE(allow_install)) {
        stop("Runtime annotation-package installation is disabled. Restore the frozen environment before running Part2.", call. = FALSE)
    }

    if (!requireNamespace(pkg, quietly = TRUE)) {
        msg(pkg, " is not installed; skipping package mapping.", level = "WARN")
        return(data.frame())
    }

    db <- get(pkg, envir = asNamespace(pkg))
    keytypes <- AnnotationDbi::keytypes(db)
    kt <- intersect(c("PROBEID", "PROBE", "ID"), keytypes)[1]
    if (is.na(kt)) kt <- keytypes[1]

    cols <- intersect(c("SYMBOL", "GENENAME", "ENTREZID"), AnnotationDbi::columns(db))
    if (!"SYMBOL" %in% cols) {
        msg(pkg, " has no SYMBOL column; skipping package mapping.", level = "WARN")
        return(data.frame())
    }

    tab <- tryCatch(
        AnnotationDbi::select(
            db,
            keys = intersect(probe_ids, AnnotationDbi::keys(db, keytype = kt)),
            keytype = kt,
            columns = "SYMBOL"
        ),
        error = function(e) {
            msg("AnnotationDbi mapping failed for ", pkg, ": ", conditionMessage(e), level = "WARN")
            data.frame()
        }
    )

    if (!nrow(tab)) return(data.frame())
    names(tab) <- toupper(names(tab))
    if (!all(c(kt, "SYMBOL") %in% names(tab))) {
        ## After toupper, kt may be already uppercase.
        key_col <- names(tab)[1]
    } else {
        key_col <- kt
    }

    out <- data.frame(
        probe_id = as.character(tab[[key_col]]),
        gene_symbol = toupper(trimws(as.character(tab$SYMBOL))),
        source = paste0(pkg, ":", kt, ":SYMBOL"),
        stringsAsFactors = FALSE
    )
    out <- out[nzchar(out$gene_symbol) & !is.na(out$gene_symbol), , drop = FALSE]
    out <- out[!duplicated(paste(out$probe_id, out$gene_symbol)), , drop = FALSE]
    msg("Mapped ", length(unique(out$probe_id)), " probes using ", pkg)
    out
}

map_probe_table <- function(tab, probe_ids, source_prefix = "table") {
    tab <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
    id_col <- grep("^ID$|^ID_REF$|probe|Probe_Id|ProbeID", names(tab), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(id_col)) {
        pid <- rownames(tab)
    } else {
        pid <- as.character(tab[[id_col]])
    }

    sym_cols <- grep(
        "gene.?symbol|symbol|gene_assignment|gene.?name|GENE_SYMBOL|Gene Symbol|ILMN_Gene|ILMN_GeneSymbol",
        names(tab),
        ignore.case = TRUE,
        value = TRUE
    )

    rows <- list()
    for (cc in sym_cols) {
        sym <- extract_first_symbol(tab[[cc]])
        keep <- pid %in% probe_ids & nzchar(sym)
        if (any(keep)) {
            rows[[cc]] <- data.frame(
                probe_id = pid[keep],
                gene_symbol = sym[keep],
                source = paste0(source_prefix, ":", cc),
                stringsAsFactors = FALSE
            )
        }
    }
    bind_rows(rows)
}

extract_first_symbol <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    vapply(x, function(s) {
        z <- unlist(strsplit(s, "///|//|;|,|\\|", perl = TRUE))
        z <- trimws(z)
        z <- gsub("^.*gene symbol[:= ]+", "", z, ignore.case = TRUE)
        z <- gsub("\\s.*$", "", z)
        z <- gsub("[^A-Za-z0-9_.-]", "", z)
        z <- toupper(z[nzchar(z)])
        z <- z[!grepl("^LOC[0-9]+$", z)]
        if (length(z)) z[1] else ""
    }, character(1))
}

gse84437_extract_olfml2b_target <- function(expr_probe, mapping, platform, target_gene = "OLFML2B", target_entrez = "25903") {
    target_gene <- toupper(target_gene)
    target_entrez <- as.character(target_entrez)
    usable <- mapping$usable
    all_map <- mapping$all
    if (is.null(usable) || !is.data.frame(usable) || !nrow(usable)) {
        stop("GSE84437 mapping is empty; cannot audit OLFML2B probes.", call. = FALSE)
    }
    usable$probe_id <- as.character(usable$probe_id)
    usable$gene_symbol <- toupper(trimws(as.character(usable$gene_symbol)))
    target_rows <- usable[usable$gene_symbol == target_gene & usable$probe_id %in% rownames(expr_probe), , drop = FALSE]
    if (!nrow(target_rows)) {
        stop("GSE84437 has no unambiguous OLFML2B probe after platform mapping.", call. = FALSE)
    }
    probes <- unique(target_rows$probe_id)
    target_probe_expression <- expr_probe[probes, , drop = FALSE]
    target_mapping <- data.frame(
        probe_id = probes,
        symbol = target_gene,
        entrez_id = target_entrez,
        mapping_source = vapply(probes, function(p) paste(unique(target_rows$source[target_rows$probe_id == p]), collapse = ";"), character(1)),
        symbol_annotations = target_gene,
        entrez_annotations = target_entrez,
        target_match_reason = "symbol",
        ambiguous_symbol = FALSE,
        in_expression = probes %in% rownames(expr_probe),
        primary_eligible = TRUE,
        cohort = "GSE84437",
        platform = platform,
        stringsAsFactors = FALSE
    )
    target_probe_stats <- data.frame(
        probe_id = probes,
        symbol = target_gene,
        entrez_id = target_entrez,
        mean = rowMeans(target_probe_expression, na.rm = TRUE),
        sd = matrixStats::rowSds(target_probe_expression, na.rm = TRUE),
        iqr = matrixStats::rowIQRs(target_probe_expression, na.rm = TRUE),
        missing_fraction = rowMeans(!is.finite(target_probe_expression)),
        cohort = "GSE84437",
        platform = platform,
        stringsAsFactors = FALSE
    )
    evidence <- all_map
    if (!is.null(evidence) && is.data.frame(evidence) && nrow(evidence)) {
        evidence <- evidence[evidence$probe_id %in% probes, , drop = FALSE]
        evidence$cohort <- "GSE84437"
        evidence$platform <- platform
    } else {
        evidence <- target_mapping
    }
    list(
        target_mapping = target_mapping,
        target_mapping_evidence = evidence,
        target_probes = probes,
        target_probe_expression = target_probe_expression,
        target_probe_stats = target_probe_stats
    )
}

aggregate_probes_to_genes <- function(expr_probe, mapping) {
    if (is.null(mapping) || nrow(mapping) == 0L) {
        stop("No usable platform mapping. Cannot build gene-level GSE84437 matrix.", call. = FALSE)
    }
    mapping <- mapping[mapping$probe_id %in% rownames(expr_probe), , drop = FALSE]
    if (!nrow(mapping)) stop("Usable mapping has no overlap with probe expression matrix.", call. = FALSE)

    vars <- matrixStats::rowVars(expr_probe, na.rm = TRUE)
    names(vars) <- rownames(expr_probe)

    rows <- list()
    for (g in unique(mapping$gene_symbol)) {
        probes <- unique(mapping$probe_id[mapping$gene_symbol == g])
        probes <- probes[probes %in% rownames(expr_probe)]
        if (!length(probes)) next
        if (length(probes) == 1L) {
            rows[[g]] <- expr_probe[probes, ]
        } else {
            sel <- probes[which.max(vars[probes])]
            rows[[g]] <- expr_probe[sel, ]
        }
    }
    mat <- do.call(rbind, rows)
    rownames(mat) <- names(rows)
    mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
    mat
}

normalize_sample_order <- function(expr, clinical) {
    if (!"sample_id" %in% names(clinical)) clinical$sample_id <- rownames(clinical)
    common <- intersect(colnames(expr), clinical$sample_id)
    if (length(common) >= 20L) {
        expr <- expr[, common, drop = FALSE]
        clinical <- clinical[match(common, clinical$sample_id), , drop = FALSE]
    } else if (ncol(expr) == nrow(clinical)) {
        clinical$sample_id <- colnames(expr)
        rownames(clinical) <- clinical$sample_id
    } else {
        stop("Cannot align GSE84437 expression samples with clinical sample_id.", call. = FALSE)
    }
    list(expr = expr, clinical = clinical)
}

make_gse84437_audit <- function(accession, platform, expr_probe, expr_gene, clinical, mapping, transformed, sm_file) {
    os_keep <- is.finite(clinical$os_time_days) & clinical$os_time_days > 0 & clinical$os_event %in% c(0L, 1L)
    data.frame(
        cohort = accession,
        platform = platform,
        source_file = sm_file,
        n_probes = nrow(expr_probe),
        n_genes = nrow(expr_gene),
        n_samples_expression = ncol(expr_gene),
        n_samples_clinical = nrow(clinical),
        n_mapping_all = nrow(mapping$all),
        n_mapping_usable = nrow(mapping$usable),
        n_ambiguous_probe_rows = nrow(mapping$ambiguous_probes),
        transformed = transformed$transformed,
        transform_reason = transformed$reason,
        n_os_complete = sum(os_keep),
        os_events = sum(clinical$os_event[os_keep] == 1L, na.rm = TRUE),
        n_dfs_complete = 0L,
        dfs_events = 0L,
        n_recurrence_complete = 0L,
        recurrence_events = 0L,
        n_stage = if ("stage" %in% names(clinical)) sum(!is.na(clinical$stage)) else 0L,
        n_age = if ("age" %in% names(clinical)) sum(is.finite(clinical$age)) else 0L,
        recurrence_endpoint = NA_character_,
        publication_evaluable_os = sum(os_keep) >= 30L && sum(clinical$os_event[os_keep] == 1L, na.rm = TRUE) >= 20L,
        publication_evaluable_dfs = FALSE,
        stringsAsFactors = FALSE
    )
}

update_part2_geo_index <- function(index_path, accession, object_path, audit, target_mapping = data.frame()) {
    if (file.exists(index_path)) {
        index <- tryCatch(readRDS(index_path), error = function(e) NULL)
    } else {
        index <- NULL
    }
    if (is.null(index)) {
        index <- list(
            version = GSE84437_PART2_VERSION,
            generated_at = timestamp(),
            cohort_files = setNames(character(), character()),
            cohort_audit = data.frame(),
            endpoint_audit = data.frame(),
            target_mapping_audit = data.frame()
        )
    }

    cf <- index$cohort_files
    if (is.null(cf)) cf <- setNames(character(), character())
    cf <- cf[names(cf) != accession]
    cf <- c(cf, setNames(normalizePath(object_path, winslash = "/", mustWork = FALSE), accession))
    index$cohort_files <- cf

    old_audit <- index$cohort_audit
    if (is.null(old_audit)) old_audit <- data.frame()
    if (nrow(old_audit) && "cohort" %in% names(old_audit)) {
        old_audit <- old_audit[old_audit$cohort != accession, , drop = FALSE]
    }
    index$cohort_audit <- bind_rows(list(old_audit, audit))

    endpoint <- audit[, intersect(c(
        "cohort", "platform", "n_samples_expression", "recurrence_endpoint",
        "n_os_complete", "os_events", "n_recurrence_complete", "recurrence_events",
        "publication_evaluable_os", "publication_evaluable_dfs"
    ), names(audit)), drop = FALSE]
    old_ep <- index$endpoint_audit
    if (is.null(old_ep)) old_ep <- data.frame()
    if (nrow(old_ep) && "cohort" %in% names(old_ep)) {
        old_ep <- old_ep[old_ep$cohort != accession, , drop = FALSE]
    }
    index$endpoint_audit <- bind_rows(list(old_ep, endpoint))

    old_map <- index$target_mapping_audit
    if (is.null(old_map)) old_map <- data.frame()
    if (nrow(old_map) && "cohort" %in% names(old_map)) {
        old_map <- old_map[old_map$cohort != accession, , drop = FALSE]
    }
    if (!is.null(target_mapping) && is.data.frame(target_mapping) && nrow(target_mapping)) {
        index$target_mapping_audit <- bind_rows(list(old_map, target_mapping))
        if (exists("olfml2b_assert_target_mapping_clean", mode = "function")) {
            olfml2b_assert_target_mapping_clean(index$target_mapping_audit, target_gene = "OLFML2B")
        }
    } else {
        index$target_mapping_audit <- old_map
    }

    part2_table_dir <- dir_create(file.path(dirname(dirname(index_path)), "tables", "Part2"))
    write_csv(index$cohort_audit, file.path(part2_table_dir, "Part2_GEO_bulk_cohort_audit.csv"))
    write_csv(index$endpoint_audit, file.path(part2_table_dir, "Part2_GEO_endpoint_audit.csv"))
    write_csv(index$target_mapping_audit, file.path(part2_table_dir, "Part2_GEO_OLFML2B_mapping_audit.csv"))

    index$version <- paste(OLFML2B_PART2_PATCH_VERSION, GSE84437_PART2_VERSION, sep = "+")
    index$release_contract_version <- OLFML2B_PART2_PATCH_VERSION
    index$generated_at <- timestamp()

    saveRDS(index, index_path, compress = "xz")
    index
}

find_local_gpl_file <- function(platform, gpl_dir, cohort_dir) {
    if (is.na(platform) || !nzchar(platform)) return(NA_character_)
    files <- c(
        list.files(gpl_dir, pattern = paste0(platform, ".*"), recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
        list.files(cohort_dir, pattern = paste0(platform, ".*"), recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    )
    files <- unique(files[file.exists(files)])
    files <- files[grepl("\\.txt$|\\.txt\\.gz$|\\.soft$|\\.soft\\.gz$|\\.annot$|\\.annot\\.gz$", files, ignore.case = TRUE)]
    if (!length(files)) return(NA_character_)
    files[order(nchar(files), files)][1]
}

try_download_gpl <- function(platform, gpl_dir, timeout_sec = 600L) {
    dir_create(gpl_dir)
    dest <- file.path(gpl_dir, paste0(platform, "_family.soft.gz"))
    if (file.exists(dest) && file.info(dest)$size >= 10000L) {
        msg("Using existing GPL soft file; skipping download: ", dest)
        return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }
    tmp <- paste0(dest, ".partial")

    urls <- c(
        sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%snnn/%s/soft/%s_family.soft.gz",
                substr(platform, 1, nchar(platform) - 3), platform, platform),
        sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full", platform)
    )

    old_timeout <- getOption("timeout")
    options(timeout = max(as.integer(timeout_sec), old_timeout %||% 60L, 600L))
    on.exit(options(timeout = old_timeout), add = TRUE)

    for (u in urls) {
        if (file.exists(tmp)) unlink(tmp, force = TRUE)
        msg("Downloading GPL from: ", u)
        ok <- tryCatch({
            utils::download.file(
                url = u,
                destfile = tmp,
                mode = "wb",
                method = "libcurl",
                quiet = FALSE
            )
            TRUE
        }, error = function(e) {
            msg("GPL download attempt failed: ", conditionMessage(e), level = "WARN")
            FALSE
        })
        if (!ok || !file.exists(tmp) || file.info(tmp)$size < 10000L) next

        ## Validate that it contains a platform table and gives usable text.
        valid <- tryCatch({
            lines <- read_text_lines(tmp)
            any(grepl("!platform_table_begin", lines, ignore.case = TRUE)) &&
                any(grepl("!platform_table_end", lines, ignore.case = TRUE))
        }, error = function(e) FALSE)

        if (isTRUE(valid)) {
            if (file.exists(dest)) unlink(dest, force = TRUE)
            file.rename(tmp, dest)
            return(normalizePath(dest, winslash = "/", mustWork = TRUE))
        } else {
            msg("Downloaded GPL file did not contain a valid platform table; ignoring partial file.", level = "WARN")
            unlink(tmp, force = TRUE)
        }
    }

    NA_character_
}

read_gpl_table <- function(path) {
    lines <- read_text_lines(path)
    begin <- grep("!platform_table_begin", lines, ignore.case = TRUE)
    end <- grep("!platform_table_end", lines, ignore.case = TRUE)
    if (length(begin) && length(end) && end[1] > begin[1]) {
        lines <- lines[(begin[1] + 1L):(end[1] - 1L)]
    } else {
        lines <- lines[!grepl("^!", lines)]
    }
    if (length(lines) < 2L) return(NULL)
    tryCatch(
        read.delim(textConnection(lines), header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
    )
}

expand_characteristics <- function(pdata) {
    x <- as.data.frame(pdata, stringsAsFactors = FALSE, check.names = FALSE)
    names(x) <- make.unique(clean_names(names(x)))
    char_cols <- grep("characteristics", names(x), ignore.case = TRUE, value = TRUE)
    extracted <- list()

    if (length(char_cols)) {
        for (cc in char_cols) {
            vals <- as.character(x[[cc]])
            for (i in seq_along(vals)) {
                val <- vals[i]
                if (is.na(val) || !nzchar(val)) next
                pieces <- if (grepl(";", val, fixed = TRUE) && length(gregexpr(":", val, fixed = TRUE)[[1]]) > 1L) {
                    strsplit(val, ";", fixed = TRUE)[[1]]
                } else {
                    val
                }
                for (piece in pieces) {
                    pos <- regexpr("\\s*:\\s*", piece, perl = TRUE)
                    if (pos[1] < 1L) next
                    key <- clean_names(substr(piece, 1L, pos[1] - 1L))
                    value <- trimws(substr(piece, pos[1] + attr(pos, "match.length"), nchar(piece)))
                    if (!nzchar(key)) next
                    if (is.null(extracted[[key]])) extracted[[key]] <- rep(NA_character_, nrow(x))
                    if (is.na(extracted[[key]][i]) || !nzchar(extracted[[key]][i])) extracted[[key]][i] <- value
                }
            }
        }
    }

    if (length(extracted)) {
        extra <- as.data.frame(extracted, stringsAsFactors = FALSE, check.names = FALSE)
        dup <- intersect(names(extra), names(x))
        if (length(dup)) names(extra)[match(dup, names(extra))] <- paste0(dup, "_characteristic")
        x <- cbind(x, extra)
    }
    x
}

transform_expression_if_needed <- function(expr) {
    q <- suppressWarnings(quantile(expr, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
    do_log <- is.finite(q[6]) && (q[6] > 100 || q[5] > 50 || (q[6] - q[1]) > 100)
    if (do_log) {
        list(matrix = log2(expr + 1), transformed = TRUE, reason = "quantile_rule_log2_x_plus_1")
    } else {
        list(matrix = expr, transformed = FALSE, reason = "already_log_scale_or_no_transform_needed")
    }
}

stage_factor <- function(x) {
    y <- toupper(trimws(as.character(x)))
    y <- gsub("PATHOLOGIC|PATHOLOGICAL|AJCC|STAGE|[ :_-]", "", y)
    y <- gsub("A|B|C", "", y)
    out <- rep(NA_character_, length(y))
    out[grepl("^0", y)] <- "Stage 0"
    out[grepl("^I$|^1$", y)] <- "Stage I"
    out[grepl("^II$|^2$", y)] <- "Stage II"
    out[grepl("^III$|^3$", y)] <- "Stage III"
    out[grepl("^IV$|^4$", y)] <- "Stage IV"
    factor(out, levels = c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"))
}

stage_numeric <- function(x) {
    f <- stage_factor(x)
    as.numeric(f) - 1
}

sex_factor <- function(x) {
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_character_, length(y))
    out[y %in% c("female", "f", "woman")] <- "Female"
    out[y %in% c("male", "m", "man")] <- "Male"
    factor(out, levels = c("Female", "Male"))
}

clean_names <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    make.unique(x)
}

bind_rows <- function(xs) {
    if (is.null(xs)) return(data.frame())
    if (is.data.frame(xs)) return(xs)
    xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
    if (!length(xs)) return(data.frame())
    xs <- lapply(xs, function(x) {
        names(x) <- make.unique(names(x))
        x
    })
    cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
    aligned <- lapply(xs, function(x) {
        miss <- setdiff(cols, names(x))
        for (m in miss) x[[m]] <- NA
        x[, cols, drop = FALSE]
    })
    out <- do.call(rbind, aligned)
    rownames(out) <- NULL
    out
}

read_text_lines <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, warn = FALSE)
}

file_manifest <- function(path) {
    if (!dir.exists(path)) return(data.frame())
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    files <- files[!grepl("\\.partial$|\\.part$|\\.tmp$", files, ignore.case = TRUE)]
    info <- file.info(files)
    data.frame(
        path = normalizePath(files, winslash = "/", mustWork = FALSE),
        file = basename(files),
        size_bytes = info$size,
        modified = as.character(info$mtime),
        stringsAsFactors = FALSE
    )
}

dir_create <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
    dir_create(dirname(path))
    if (is.null(x)) x <- data.frame()
    write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    invisible(path)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required package is not installed: ", pkg, call. = FALSE)
    }
}

msg <- function(..., level = "INFO") {
    message("[", level, "] [OLFML2B-STAD Part2] ", paste0(..., collapse = ""))
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}


# ==============================================================================
# OLFML2B explicit Part2 entry: original Part2 + GSE84437 extension
# ==============================================================================

olfml2b_patch_gse84437_target <- function(ctx) {
  fp <- file.path(ctx$dirs$objects, "Part2_GSE84437.rds")
  if (!file.exists(fp)) return(invisible(FALSE))
  obj <- readRDS(fp)
  if (!is.null(obj$expression) && "OLFML2B" %in% rownames(obj$expression) && is.data.frame(obj$sample_metadata)) {
    meta <- obj$sample_metadata
    if (nrow(meta) != 433L || any(grepl("^GSM441", meta$sample_id)) ||
        !"analysis_role" %in% names(meta) || any(meta$analysis_role != "FORMAL_OS")) {
      stop("Refusing to patch GSE84437 target: object is not the provenance-partitioned 433-sample formal OS object.", call. = FALSE)
    }
    meta$olfml2b_expression <- as.numeric(obj$expression["OLFML2B", meta$sample_id])
    meta$olfml2b_z_within_subseries <- ave(
      meta$olfml2b_expression, meta$source_subseries,
      FUN = function(z) olfml2b_z(as.numeric(z))
    )
    meta$olfml2b_z <- olfml2b_z(meta$olfml2b_z_within_subseries)
    batch_means <- tapply(meta$olfml2b_z, meta$source_subseries, mean, na.rm = TRUE)
    if (abs(stats::sd(meta$olfml2b_z, na.rm = TRUE) - 1) > 1e-8 || any(abs(batch_means) > 1e-8)) {
      stop("GSE84437 target z-score is not standardized within the 433-sample formal cohort.", call. = FALSE)
    }
    obj$sample_metadata <- meta
    obj$target_gene <- "OLFML2B"
    obj$target_measured <- TRUE
    obj$target_rows <- "OLFML2B"
    if (is.null(obj$target_mapping) || !is.data.frame(obj$target_mapping) || !nrow(obj$target_mapping)) {
      target_info <- gse84437_extract_olfml2b_target(
        expr_probe = obj$probe_expression,
        mapping = list(usable = obj$platform_mapping, all = obj$platform_mapping_all),
        platform = obj$platform,
        target_gene = "OLFML2B",
        target_entrez = olfml2b_target_entrez_id()
      )
      obj$target_mapping <- target_info$target_mapping
      obj$target_mapping_evidence <- target_info$target_mapping_evidence
      obj$target_probes <- target_info$target_probes
      obj$target_probe_expression <- target_info$target_probe_expression
      obj$target_probe_stats <- target_info$target_probe_stats
    }
    obj$gene_annotation_factcheck <- olfml2b_gene_annotation()
    obj$release_contract_version <- OLFML2B_PART2_PATCH_VERSION
    olfml2b_atomic_save_rds(obj, fp)
  }
  invisible(TRUE)
}


olfml2b_part2_extract_clinical_frame <- function(obj) {
  if (is.data.frame(obj)) return(obj)
  if (!is.list(obj)) return(data.frame())
  for (nm in c("sample_metadata", "clinical", "metadata", "pheno", "pdata")) {
    if (!is.null(obj[[nm]]) && is.data.frame(obj[[nm]]) && nrow(obj[[nm]]) > 0L) {
      return(obj[[nm]])
    }
  }
  frames <- list()
  collect_frames <- function(x) {
    if (is.data.frame(x) && nrow(x) > 0L) {
      frames[[length(frames) + 1L]] <<- x
    } else if (is.list(x)) {
      for (i in seq_along(x)) collect_frames(x[[i]])
    }
  }
  collect_frames(obj)
  if (!length(frames)) return(data.frame())
  scored <- vapply(frames, function(d) {
    sum(c("os_time_days", "os_event", "sample_id") %in% names(d)) * 1000L + nrow(d)
  }, numeric(1))
  frames[[which.max(scored)]]
}

olfml2b_part2_endpoint_counts <- function(clinical) {
  if (!is.data.frame(clinical) || !nrow(clinical)) {
    return(list(n_complete = 0L, n_events = 0L))
  }
  complete <- if (all(c("os_time_days", "os_event") %in% names(clinical))) {
    is.finite(clinical$os_time_days) & clinical$os_time_days > 0 & clinical$os_event %in% c(0L, 1L)
  } else rep(FALSE, nrow(clinical))
  n_complete <- sum(complete)
  n_events <- if ("os_event" %in% names(clinical)) sum(clinical$os_event[complete] == 1L, na.rm = TRUE) else 0L
  list(n_complete = as.integer(n_complete), n_events = as.integer(n_events))
}

olfml2b_part2_stage_association_audit <- function(index, ctx) {
  rows <- lapply(names(index$cohort_files), function(cohort) {
    fp <- as.character(index$cohort_files[[cohort]])
    obj <- if (file.exists(fp)) readRDS(fp) else list()
    d <- olfml2b_part2_extract_clinical_frame(obj)
    if (!nrow(d) || !all(c("stage", "olfml2b_expression") %in% names(d))) return(NULL)
    keep <- !is.na(d$stage) & is.finite(d$olfml2b_expression)
    stage_num <- olfml2b_stage_numeric(d$stage[keep])
    keep2 <- is.finite(stage_num)
    y <- d$olfml2b_expression[keep][keep2]
    s <- stage_num[keep2]
    if (length(y) < 20L || length(unique(s)) < 2L) return(NULL)
    ct <- suppressWarnings(stats::cor.test(y, s, method = "spearman", exact = FALSE))
    sf <- factor(s, levels = sort(unique(s)))
    kw <- stats::kruskal.test(y ~ sf)
    k <- nlevels(sf); n <- length(y)
    eps2 <- if (n > k) max(0, (as.numeric(kw$statistic) - k + 1) / (n - k)) else NA_real_
    data.frame(
      cohort = cohort, n = n, n_stage_groups = k,
      spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
      kruskal_wallis_H = unname(kw$statistic), kruskal_wallis_p = kw$p.value,
      epsilon_squared = eps2,
      interpretation = "exploratory_association_not_a_quality_gate",
      stringsAsFactors = FALSE
    )
  })
  tab <- olfml2b_bind_rows(rows)
  if (nrow(tab)) {
    tab$spearman_fdr_bh <- stats::p.adjust(tab$spearman_p, method = "BH")
    tab$kruskal_wallis_fdr_bh <- stats::p.adjust(tab$kruskal_wallis_p, method = "BH")
    tab$replicated_monotonic_stage_association <- FALSE
  }
  pdir <- olfml2b_part_paths(ctx, "Part2")
  olfml2b_atomic_write_csv(tab, file.path(pdir$tables, "Part2_OLFML2B_stage_association_audit.csv"))
  tab
}

olfml2b_part2_release_contract_audit <- function(index) {
  expected <- data.frame(
    cohort = c("GSE62254", "GSE15459", "GSE26253", "GSE84437"),
    expected_n = c(300L, 192L, 432L, 433L), stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(expected)), function(i) {
    co <- expected$cohort[i]
    fp <- index$cohort_files[[co]] %||% NA_character_
    obj <- if (!is.na(fp) && file.exists(fp)) readRDS(fp) else list()
    d <- olfml2b_part2_extract_clinical_frame(obj)
    n <- nrow(d)
    unique_samples <- if ("sample_id" %in% names(d)) length(unique(d$sample_id)) else 0L
    unique_patients <- if ("patient_id" %in% names(d)) length(unique(d$patient_id)) else NA_integer_
    extra_ok <- TRUE; detail <- "PASS"
    if (co == "GSE62254") {
      cnt <- olfml2b_part2_endpoint_counts(d)
      dfs_complete <- if (all(c("dfs_time_days", "dfs_event") %in% names(d))) {
        is.finite(d$dfs_time_days) & d$dfs_time_days > 0 & d$dfs_event %in% c(0L, 1L)
      } else rep(FALSE, n)
      nonmissing_text <- function(field) {
        if (!field %in% names(d)) return(character())
        x <- trimws(as.character(d[[field]]))
        x[!is.na(x) & nzchar(x)]
      }
      distribution <- function(field) {
        x <- nonmissing_text(field)
        if (!length(x)) return("none")
        z <- table(x)
        paste(paste0(names(z), "=", as.integer(z)), collapse = "/")
      }
      overall <- nonmissing_text("stage_overall")
      p_t <- nonmissing_text("stage_pT")
      p_n <- nonmissing_text("stage_pN")
      p_m <- nonmissing_text("stage_pM")
      subtype <- nonmissing_text("molecular_subtype")
      tnm_ok <- length(overall) >= 298L &&
        all(c("Stage I", "Stage II", "Stage III", "Stage IV") %in% unique(overall)) &&
        length(p_t) >= 298L && all(c("T2", "T3", "T4") %in% unique(p_t)) &&
        length(p_n) == 300L && all(c("N0", "N1", "N2", "N3") %in% unique(p_n)) &&
        length(p_m) >= 297L && all(c("M0", "M1") %in% unique(p_m))
      subtype_ok <- length(subtype) == 300L &&
        all(c("MSS_TP53_inactive", "MSS_TP53_active", "MSI", "EMT") %in% unique(subtype))
      extra_ok <- n == 300L && cnt$n_complete == 300L && cnt$n_events == 152L &&
        sum(dfs_complete) == 282L && sum(d$dfs_event[dfs_complete] == 1L, na.rm = TRUE) == 125L &&
        tnm_ok && subtype_ok
      detail <- paste0(
        "OS_complete=", cnt$n_complete, ";OS_events=", cnt$n_events,
        ";DFS_complete=", sum(dfs_complete), ";DFS_events=", sum(d$dfs_event[dfs_complete] == 1L, na.rm = TRUE),
        ";stage=", distribution("stage_overall"),
        ";pT=", distribution("stage_pT"),
        ";pN=", distribution("stage_pN"),
        ";pM=", distribution("stage_pM"),
        ";subtype=", distribution("molecular_subtype"),
        ";row_matched_identity_gate=", tnm_ok && subtype_ok
      )
    } else if (co == "GSE15459") {
      cnt <- olfml2b_part2_endpoint_counts(d)
      extra_ok <- n == 192L && !any(grepl("EXCLUDED", d$sample_title, ignore.case = TRUE)) &&
        cnt$n_complete == 191L && cnt$n_events == 95L
      detail <- paste0("retained=", n, ";OS_complete=", cnt$n_complete, ";events=", cnt$n_events)
    } else if (co == "GSE26253") {
      extra_ok <- n == 432L && unique_patients == 432L
      detail <- paste0("unique_patient_id=", unique_patients)
    } else if (co == "GSE84437") {
      cnt <- olfml2b_part2_endpoint_counts(d)
      role_ok <- "analysis_role" %in% names(d) && all(d$analysis_role == "FORMAL_OS")
      zsd <- if ("olfml2b_z" %in% names(d)) stats::sd(d$olfml2b_z, na.rm = TRUE) else NA_real_
      batch_centered <- all(c("olfml2b_z", "source_subseries") %in% names(d)) &&
        all(abs(tapply(d$olfml2b_z, d$source_subseries, mean, na.rm = TRUE)) < 1e-8)
      sex_levels <- if ("sex" %in% names(d)) unique(as.character(stats::na.omit(d$sex))) else character()
      partition_ok <- is.data.frame(obj$superseries_partition_audit) && nrow(obj$superseries_partition_audit) == 483L &&
        sum(obj$superseries_partition_audit$analysis_role == "FORMAL_OS") == 433L &&
        sum(obj$superseries_partition_audit$analysis_role == "CONTEXT_ONLY") == 50L
      batch_ok <- is.data.frame(obj$batch_audit) && nrow(obj$batch_audit) == 3L
      context_path <- index$context_files[["GSE147163"]] %||% NA_character_
      context_ok <- !is.na(context_path) && file.exists(context_path)
      extra_ok <- n == 433L && !any(grepl("^GSM441", d$sample_id)) && role_ok &&
        cnt$n_complete == 431L && cnt$n_events == 207L && is.finite(zsd) && abs(zsd - 1) < 1e-8 &&
        all(c("Female", "Male") %in% sex_levels) && batch_centered && partition_ok && batch_ok && context_ok
      detail <- paste0("formal=", role_ok, ";OS_complete=", cnt$n_complete, ";events=", cnt$n_events,
                       ";z_sd=", signif(zsd, 6), ";sex=", paste(sort(sex_levels), collapse = "/"),
                       ";batch_centered=", batch_centered, ";partition=", partition_ok,
                       ";batch=", batch_ok, ";context_object=", context_ok)
    }
    version_ok <- identical(as.character(obj$release_contract_version %||% NA_character_),
                            OLFML2B_PART2_PATCH_VERSION)
    data.frame(cohort = co, expected_n = expected$expected_n[i], observed_n = n,
               unique_samples = unique_samples, unique_patients = unique_patients,
               version_ok = version_ok, scientific_contract_ok = extra_ok,
               status = ifelse(n == expected$expected_n[i] && unique_samples == n && version_ok && extra_ok, "PASS", "FAIL"),
               detail = detail, stringsAsFactors = FALSE)
  })
  tab <- olfml2b_bind_rows(rows)
  if (any(tab$status != "PASS")) {
    stop("Part2 release contract failed: ", paste(tab$cohort[tab$status != "PASS"], tab$detail[tab$status != "PASS"], sep = "=", collapse = "; "), call. = FALSE)
  }
  tab
}

# ==============================================================================
# Clinical-variable and TNM derivation contract (v2.5)
#
# The public cohorts do not expose the same clinical variables.  This layer
# records structural absence explicitly, retains pT/pN/pM as separate fields,
# and prevents a completely absent variable from being handed to an imputer.
# Overall stage is only promoted to the primary analysis field when it was
# reported by the source or can be derived from a complete, frozen AJCC-7
# gastric TNM tuple.  A T/N-only, M0-assumed stage is exported as sensitivity
# data and is never silently promoted to the primary stage variable.
# ==============================================================================

olfml2b_part2_nonmissing <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) return(is.finite(as.numeric(x)))
  !is.na(x) & nzchar(trimws(as.character(x)))
}

olfml2b_part2_pick_best_column <- function(d, candidates) {
  hit <- candidates[candidates %in% names(d)]
  if (!length(hit)) return(NA_character_)
  score <- vapply(hit, function(nm) sum(olfml2b_part2_nonmissing(d[[nm]])), integer(1))
  hit[which.max(score)]
}

olfml2b_part2_t_component <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\s+", "", y)
  out <- rep(NA_character_, length(y))
  out[grepl("TIS", y)] <- "TIS"
  ok <- grepl("T[0-4](?:[AB])?", y, perl = TRUE) & is.na(out)
  out[ok] <- sub(".*?(T[0-4](?:[AB])?).*", "\\1", y[ok], perl = TRUE)
  bare <- is.na(out) & grepl("^[0-4](?:[AB])?$", y)
  out[bare] <- paste0("T", y[bare])
  out
}

olfml2b_part2_n_component <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\s+", "", y)
  out <- rep(NA_character_, length(y))
  ok <- grepl("N[0-3](?:[AB])?", y, perl = TRUE)
  out[ok] <- sub(".*?(N[0-3](?:[AB])?).*", "\\1", y[ok], perl = TRUE)
  bare <- is.na(out) & grepl("^[0-3](?:[AB])?$", y)
  out[bare] <- paste0("N", y[bare])
  sub("^N3[AB]$", "N3", out)
}

olfml2b_part2_m_component <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\s+", "", y)
  out <- rep(NA_character_, length(y))
  ok <- grepl("M[01]", y, perl = TRUE)
  out[ok] <- sub(".*?(M[01]).*", "\\1", y[ok], perl = TRUE)
  bare <- is.na(out) & y %in% c("0", "1")
  out[bare] <- paste0("M", y[bare])
  out
}

olfml2b_part2_overall_stage <- function(x) {
  y <- gsub("[^A-Z0-9]", "", toupper(as.character(x)))
  y <- sub("^(AJCC|PATHOLOGIC|PATHOLOGICAL|CLINICAL)", "", y)
  y <- sub("^STAGE", "", y)
  out <- rep(NA_character_, length(y))
  out[grepl("^(IV|4)", y)] <- "Stage IV"
  out[is.na(out) & grepl("^(III|3)", y)] <- "Stage III"
  out[is.na(out) & grepl("^(II|2)", y)] <- "Stage II"
  out[is.na(out) & grepl("^(I|1)", y)] <- "Stage I"
  out[is.na(out) & grepl("^(0|IS)$", y)] <- "Stage 0"
  out
}

olfml2b_part2_derive_gastric_ajcc7 <- function(t_stage, n_stage, m_stage,
                                              assume_m0 = FALSE) {
  t <- olfml2b_part2_t_component(t_stage)
  n <- olfml2b_part2_n_component(n_stage)
  m <- olfml2b_part2_m_component(m_stage)
  stage <- rep(NA_character_, length(t))
  status <- rep("NOT_DERIVABLE_INCOMPLETE_TNM", length(t))
  assumption <- rep(NA_character_, length(t))

  for (i in seq_along(t)) {
    mi <- m[i]
    if (is.na(mi) && isTRUE(assume_m0) && !is.na(t[i]) && !is.na(n[i])) {
      mi <- "M0"
      assumption[i] <- "M0_ASSUMED_FOR_SENSITIVITY_ONLY"
    }
    if (is.na(t[i]) || is.na(n[i]) || is.na(mi)) next
    if (mi == "M1") {
      stage[i] <- "Stage IV"
    } else if (t[i] == "TIS" && n[i] == "N0") {
      stage[i] <- "Stage 0"
    } else {
      tc <- sub("[AB]$", "", t[i])
      nc <- sub("[AB]$", "", n[i])
      if (tc == "T1" && nc %in% c("N0", "N1")) stage[i] <- "Stage I"
      if (tc == "T1" && nc %in% c("N2", "N3")) stage[i] <- "Stage II"
      if (tc == "T2" && nc == "N0") stage[i] <- "Stage I"
      if (tc == "T2" && nc %in% c("N1", "N2")) stage[i] <- "Stage II"
      if (tc == "T2" && nc == "N3") stage[i] <- "Stage III"
      if (tc == "T3" && nc %in% c("N0", "N1")) stage[i] <- "Stage II"
      if (tc == "T3" && nc %in% c("N2", "N3")) stage[i] <- "Stage III"
      if (t[i] == "T4A" && nc == "N0") stage[i] <- "Stage II"
      if (t[i] == "T4A" && nc %in% c("N1", "N2", "N3")) stage[i] <- "Stage III"
      if (t[i] == "T4B") stage[i] <- "Stage III"
      # T4 without a/b and N0 is ambiguous (IIB vs IIIB in AJCC-7), so it is
      # intentionally not assigned.  For N1-N3 the coarse group is III for
      # both T4a and T4b and can be retained.
      if (t[i] == "T4" && nc %in% c("N1", "N2", "N3")) stage[i] <- "Stage III"
    }
    if (!is.na(stage[i])) {
      status[i] <- if (is.na(assumption[i])) {
        "DERIVED_AJCC7_GASTRIC_COMPLETE_TNM"
      } else {
        "DERIVED_AJCC7_GASTRIC_TN_ASSUMED_M0_SENSITIVITY"
      }
    } else if (t[i] == "T4" && n[i] == "N0" && mi == "M0") {
      status[i] <- "NOT_DERIVABLE_AMBIGUOUS_T4_SUBCATEGORY"
    } else {
      status[i] <- "NOT_DERIVABLE_UNMAPPED_TNM"
    }
  }
  data.frame(stage = stage, status = status, assumption = assumption,
             stringsAsFactors = FALSE)
}

olfml2b_part2_normalize_clinical_contract <- function(d, cohort) {
  d <- as.data.frame(d, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- make.unique(olfml2b_clean_names(names(d)), sep = "_")
  n <- nrow(d)
  get_best <- function(candidates, default = NA_character_) {
    col <- olfml2b_part2_pick_best_column(d, candidates)
    if (is.na(col)) rep(default, n) else d[[col]]
  }

  observed_overall <- olfml2b_part2_overall_stage(get_best(c(
    "stage_overall", "stage", "pathologic_stage", "pathological_stage", "ajcc_stage"
  )))
  p_t <- olfml2b_part2_t_component(get_best(c(
    "stage_pt", "ptstage", "ptstage_ch1", "p_t_stage", "t_stage", "tstage", "pathologic_t", "ajcc_pathologic_t"
  )))
  p_n <- olfml2b_part2_n_component(get_best(c(
    "stage_pn", "pnstage", "pnstage_ch1", "p_n_stage", "n_stage", "nstage", "pathologic_n", "ajcc_pathologic_n"
  )))
  p_m <- olfml2b_part2_m_component(get_best(c(
    "stage_pm", "pmstage", "pmstage_ch1", "p_m_stage", "m_stage", "mstage", "pathologic_m", "ajcc_pathologic_m"
  )))

  complete_tnm <- olfml2b_part2_derive_gastric_ajcc7(p_t, p_n, p_m, assume_m0 = FALSE)
  tn_m0_sensitivity <- olfml2b_part2_derive_gastric_ajcc7(p_t, p_n, p_m, assume_m0 = TRUE)
  primary_stage <- observed_overall
  use_complete_tnm <- is.na(primary_stage) & !is.na(complete_tnm$stage)
  primary_stage[use_complete_tnm] <- complete_tnm$stage[use_complete_tnm]
  sensitivity_stage <- primary_stage
  use_sensitivity <- is.na(sensitivity_stage) & !is.na(tn_m0_sensitivity$stage)
  sensitivity_stage[use_sensitivity] <- tn_m0_sensitivity$stage[use_sensitivity]

  origin <- rep("NOT_AVAILABLE_OR_NOT_DERIVABLE", n)
  origin[!is.na(observed_overall)] <- "SOURCE_REPORTED_OVERALL_STAGE"
  origin[use_complete_tnm] <- complete_tnm$status[use_complete_tnm]
  sensitivity_origin <- origin
  sensitivity_origin[use_sensitivity] <- tn_m0_sensitivity$status[use_sensitivity]

  # Remove fields created by an earlier contract pass before recreating them.
  # This makes finalisation idempotent when a user deliberately reruns Part2.
  generated_clean <- intersect(c(
    "stage_overall", "stage_pt", "stage_pn", "stage_pm",
    "stage_derived_ajcc7", "stage_derived_ajcc7_m0_sensitivity",
    "stage_analysis_primary", "stage_analysis_sensitivity",
    "stage_derivation_status", "stage_sensitivity_derivation_status",
    "stage_derivation_standard", "stage_derivation_assumption"
  ), names(d))
  if (length(generated_clean)) d[generated_clean] <- NULL

  d$stage_overall <- observed_overall
  d$stage_pT <- p_t
  d$stage_pN <- p_n
  d$stage_pM <- p_m
  d$stage_derived_ajcc7 <- complete_tnm$stage
  d$stage_derived_ajcc7_m0_sensitivity <- tn_m0_sensitivity$stage
  d$stage_analysis_primary <- primary_stage
  d$stage_analysis_sensitivity <- sensitivity_stage
  d$stage_derivation_status <- origin
  d$stage_sensitivity_derivation_status <- sensitivity_origin
  d$stage_derivation_standard <- "AJCC_7_GASTRIC_COARSE_GROUP"
  d$stage_derivation_assumption <- tn_m0_sensitivity$assumption

  # Backward-compatible aliases point only to the primary, non-assumption
  # analysis stage.  GSE84437 therefore remains NA in `stage`; its pT and pN
  # fields and explicitly labelled M0 sensitivity stage remain available.
  d$stage <- factor(primary_stage,
                    levels = c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"))
  d$stage_numeric <- match(primary_stage,
                           c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV")) - 1L
  d$t_stage <- p_t
  d$n_stage <- p_n
  d$m_stage <- p_m
  d
}

olfml2b_part2_field_source <- function(cohort, variable) {
  if (cohort == "GSE26253") {
    if (variable %in% c("age", "sex", "stage_pT", "stage_pN", "stage_pM"))
      return("GEO_SAMPLE_METADATA_NOT_DEPOSITED")
    if (variable == "stage_overall") return("GEO_PATHOLOGICAL_STAGE")
  }
  if (cohort == "GSE84437") {
    if (variable %in% c("age", "sex", "stage_pT", "stage_pN"))
      return("GEO_GSE84426_GSE84433_SAMPLE_CHARACTERISTICS")
    if (variable %in% c("stage_overall", "stage_pM"))
      return("GEO_SAMPLE_METADATA_NOT_DEPOSITED")
    if (variable == "stage_analysis_sensitivity")
      return("DERIVED_AJCC7_FROM_PT_PN_WITH_EXPLICIT_M0_SENSITIVITY_ASSUMPTION")
  }
  "HARMONIZED_PUBLIC_SOURCE"
}

olfml2b_part2_field_contract <- function(d, cohort) {
  variables <- c("age", "sex", "stage_overall", "stage_pT", "stage_pN", "stage_pM",
                 "stage_analysis_primary", "stage_analysis_sensitivity")
  rows <- lapply(variables, function(v) {
    x <- if (v %in% names(d)) d[[v]] else rep(NA, nrow(d))
    keep <- olfml2b_part2_nonmissing(x)
    n_obs <- sum(keep)
    n_total <- nrow(d)
    miss <- if (n_total) 1 - n_obs / n_total else NA_real_
    unique_n <- if (n_obs) length(unique(as.character(x[keep]))) else 0L
    known_structural <-
      (cohort == "GSE26253" && v %in% c("age", "sex", "stage_pT", "stage_pN", "stage_pM")) ||
      (cohort == "GSE84437" && v %in% c("stage_overall", "stage_pM", "stage_analysis_primary"))
    status <- if (n_obs == n_total && n_total > 0L) {
      "OBSERVED_COMPLETE"
    } else if (n_obs > 0L) {
      if (grepl("stage_analysis", v, fixed = TRUE)) "DERIVED_OR_OBSERVED_PARTIAL" else "PARTIALLY_MISSING"
    } else if (known_structural) {
      "STRUCTURALLY_ABSENT_OR_NOT_DERIVABLE"
    } else {
      "NOT_OBSERVED_REQUIRES_SOURCE_REVIEW"
    }
    imputation_eligible <- v %in% c("age", "sex", "stage_overall", "stage_pT", "stage_pN", "stage_pM") &&
      n_obs >= 30L && n_obs < n_total && is.finite(miss) && miss <= 0.20 &&
      unique_n >= 2L && !known_structural
    model_eligible <- n_obs >= 30L && unique_n >= 2L
    data.frame(
      clinical_contract_version = OLFML2B_PART2_CLINICAL_CONTRACT_VERSION,
      cohort = cohort, variable = v, n_total = n_total, n_observed = n_obs,
      missing_fraction = miss, n_unique_observed = unique_n,
      availability_status = status,
      source = olfml2b_part2_field_source(cohort, v),
      imputation_eligible = imputation_eligible,
      model_eligible = model_eligible,
      excessive_missingness_no_imputation = is.finite(miss) && miss > 0.20,
      stringsAsFactors = FALSE
    )
  })
  olfml2b_bind_rows(rows)
}

olfml2b_part2_finalize_clinical_contract <- function(index, ctx) {
  contracts <- list()
  derivation <- list()
  pdir <- olfml2b_part_paths(ctx, "Part2")
  for (cohort in names(index$cohort_files)) {
    fp <- as.character(index$cohort_files[[cohort]])
    if (!file.exists(fp)) next
    obj <- readRDS(fp)
    d <- olfml2b_part2_extract_clinical_frame(obj)
    if (!nrow(d)) next
    d <- olfml2b_part2_normalize_clinical_contract(d, cohort)
    contract <- olfml2b_part2_field_contract(d, cohort)
    obj$sample_metadata <- d
    obj$clinical <- d
    obj$clinical_contract_version <- OLFML2B_PART2_CLINICAL_CONTRACT_VERSION
    obj$clinical_field_contract <- contract
    obj$tnm_stage_derivation_standard <- "AJCC_7_GASTRIC_COARSE_GROUP"
    obj$tnm_stage_derivation_policy <- paste(
      "Primary stage requires source-reported overall stage or complete T/N/M.",
      "T/N-only M0 assumption is sensitivity-only and never promoted to primary stage."
    )
    olfml2b_atomic_save_rds(obj, fp)
    olfml2b_atomic_write_csv(d, file.path(ctx$dirs$derived_clinical,
                                        paste0(cohort, "_sample_metadata_harmonized.csv")))
    olfml2b_atomic_write_csv(d, file.path(pdir$tables,
                                        paste0("Part2_", cohort, "_sample_metadata.csv")))
    contracts[[cohort]] <- contract
    st <- as.data.frame(table(d$stage_sensitivity_derivation_status, useNA = "ifany"),
                        stringsAsFactors = FALSE)
    names(st) <- c("stage_derivation_status", "n")
    st$cohort <- cohort
    st$standard <- "AJCC_7_GASTRIC_COARSE_GROUP"
    derivation[[cohort]] <- st[, c("cohort", "standard", "stage_derivation_status", "n")]
  }
  contract_all <- olfml2b_bind_rows(contracts)
  derivation_all <- olfml2b_bind_rows(derivation)
  olfml2b_atomic_write_csv(contract_all,
                          file.path(pdir$tables, "Part2_clinical_field_contract.csv"))
  olfml2b_atomic_write_csv(derivation_all,
                          file.path(pdir$tables, "Part2_TNM_stage_derivation_audit.csv"))

  g262 <- contract_all[contract_all$cohort == "GSE26253" &
                         contract_all$variable %in% c("age", "sex"), , drop = FALSE]
  if (nrow(g262) != 2L || any(g262$n_observed != 0L) ||
      any(!grepl("STRUCTURALLY_ABSENT", g262$availability_status))) {
    stop("Part2 clinical contract failed: GSE26253 age/sex must be explicitly structural, not imputable.",
         call. = FALSE)
  }
  g844 <- contract_all[contract_all$cohort == "GSE84437", , drop = FALSE]
  get_n <- function(v) g844$n_observed[match(v, g844$variable)]
  if (!isTRUE(get_n("stage_overall") == 0L) || get_n("stage_pT") < 300L ||
      get_n("stage_pN") < 300L || get_n("stage_analysis_sensitivity") < 300L) {
    stop("Part2 clinical contract failed: GSE84437 must retain pT/pN and an explicitly labelled TN-assuming-M0 sensitivity stage while overall stage remains absent.",
         call. = FALSE)
  }
  list(contract = contract_all, derivation = derivation_all)
}

run_olfml2b_part2_geo <- function(ctx = NULL, include_gse84437 = TRUE) {
  if (is.null(ctx)) ctx <- olfml2b_load_context()
  olfml2b_log("INFO", "Part2 repair version: ", OLFML2B_PART2_PATCH_VERSION,
              " | forcing cohort-object regeneration before survival analysis",
              log_file = ctx$log_file %||% NULL)
  ctx$contract <- olfml2b_default_contract()
  idx <- run_part2(ctx)
  if (isTRUE(include_gse84437)) {
    olfml2b_build_gse84437_part2_object(root = ctx$dirs$root, overwrite = TRUE)
    olfml2b_patch_gse84437_target(ctx)
    idx <- readRDS(file.path(ctx$dirs$objects, "Part2_GEO_bulk_index.rds"))
    if (is.data.frame(idx$target_mapping_audit) && nrow(idx$target_mapping_audit)) {
      olfml2b_assert_target_mapping_clean(idx$target_mapping_audit, target_gene = "OLFML2B")
    }
  }
  idx$release_contract_version <- OLFML2B_PART2_PATCH_VERSION
  idx$version <- paste(OLFML2B_PART2_PATCH_VERSION, GSE84437_PART2_VERSION, sep = "+")
  idx$context_files <- c(GSE147163 = file.path(ctx$dirs$objects, "Part2_GSE147163_context_only.rds"))
  clinical_contract <- olfml2b_part2_finalize_clinical_contract(idx, ctx)
  idx$clinical_contract_version <- OLFML2B_PART2_CLINICAL_CONTRACT_VERSION
  idx$clinical_field_contract <- clinical_contract$contract
  idx$tnm_stage_derivation_audit <- clinical_contract$derivation
  release_audit <- olfml2b_part2_release_contract_audit(idx)
  stage_audit <- olfml2b_part2_stage_association_audit(idx, ctx)
  idx$release_contract_audit <- release_audit
  idx$stage_association_audit <- stage_audit
  olfml2b_atomic_save_rds(idx, file.path(ctx$dirs$objects, "Part2_GEO_bulk_index.rds"))
  pdir <- olfml2b_part_paths(ctx, "Part2")
  olfml2b_atomic_write_csv(release_audit, file.path(pdir$tables, "Part2_release_contract_audit.csv"))
  idx_tbl <- data.frame(cohort = names(idx$cohort_files), file = as.character(idx$cohort_files),
                        exists = file.exists(as.character(idx$cohort_files)), stringsAsFactors = FALSE)
  g15459_path <- idx$cohort_files[["GSE15459"]] %||% file.path(ctx$dirs$objects, "GSE15459_harmonized.rds")
  if (!file.exists(g15459_path)) stop("Part2 final gate failed: GSE15459 cohort object is missing: ", g15459_path, call. = FALSE)
  g15459 <- readRDS(g15459_path)
  g15459_clinical <- olfml2b_part2_extract_clinical_frame(g15459)
  g15459_counts <- olfml2b_part2_endpoint_counts(g15459_clinical)
  g15459_complete <- g15459_counts$n_complete
  g15459_events <- g15459_counts$n_events
  if (g15459_complete != 191L || g15459_events != 95L) {
    stop("Part2 final gate failed after object save: GSE15459 complete OS n=", g15459_complete,
         ", events=", g15459_events, ". Refusing to run Part3 on a stale or damaged clinical object. ",
         "The gate inspected sample_metadata/clinical-compatible metadata fields in the cohort object.", call. = FALSE)
  }
  olfml2b_log("INFO", "Part2 final GSE15459 endpoint gate PASS | complete OS n=", g15459_complete,
              " | events=", g15459_events, " | metadata_field=sample_metadata/clinical-compatible", log_file = ctx$log_file %||% NULL)
  olfml2b_atomic_write_csv(idx_tbl, file.path(pdir$tables, "Part2_OLFML2B_final_index_audit.csv"))
  invisible(idx)
}

# ==============================================================================
# Part2 v1.1.0: tumor-only GEO role and source-endpoint display contract
# ==============================================================================
OLFML2B_PART2_PATCH_VERSION <- "v1.1.0_20260722_GEO_TUMOR_ONLY_CONTEXT_AND_SOURCE_ENDPOINT_FIX"
OLFML2B_PART2_SAMPLE_CONTEXT_VERSION <- "v1.0.0_20260722_THREE_LEVEL_SAMPLE_CONTEXT"

olfml2b_part2_geo_role_contract <- function() {
  data.frame(
    cohort = c("GSE62254", "GSE15459", "GSE26253", "GSE84437", "GSE147163"),
    expected_context = rep("TUMOR_ONLY", 5L),
    analysis_role = c("FORMAL_PROGNOSTIC", "FORMAL_PROGNOSTIC", "FORMAL_PROGNOSTIC", "FORMAL_PROGNOSTIC", "CONTEXT_ONLY"),
    canonical_endpoint_family = c("RECURRENCE", "OS", "RECURRENCE", "OS", "NONE"),
    source_endpoint_label = c("DFS", "OS", "RFS", "OS", "NONE_CONTEXT_ONLY"),
    normal_comparison_allowed = FALSE,
    stringsAsFactors = FALSE
  )
}

olfml2b_part2_apply_sample_context <- function(d, cohort) {
  d <- as.data.frame(d, stringsAsFactors = FALSE, check.names = FALSE)
  role <- olfml2b_part2_geo_role_contract()
  rr <- role[role$cohort == cohort, , drop = FALSE]
  if (!nrow(rr)) stop("No frozen GEO sample-context role for ", cohort, call. = FALSE)
  n <- nrow(d)
  d$sample_context <- rep("Tumor", n)
  d$is_tumor <- rep(TRUE, n)
  d$is_normal <- rep(FALSE, n)
  d$expected_sample_context <- rr$expected_context[1L]
  d$normal_comparison_allowed <- FALSE
  d$canonical_endpoint_family <- rr$canonical_endpoint_family[1L]
  d$source_endpoint_label <- rr$source_endpoint_label[1L]
  if (cohort == "GSE62254") d$recurrence_endpoint <- "DFS"
  if (cohort == "GSE26253") d$recurrence_endpoint <- "RFS"
  d$sample_context_contract_version <- OLFML2B_PART2_SAMPLE_CONTEXT_VERSION
  d
}

olfml2b_part2_sample_context_audit <- function(index, ctx) {
  role <- olfml2b_part2_geo_role_contract()
  files <- as.list(index$cohort_files)
  context_fp <- file.path(ctx$dirs$objects, "Part2_GSE147163_context_only.rds")
  if (file.exists(context_fp)) files[["GSE147163"]] <- context_fp
  rows <- list()
  for (cohort in names(files)) {
    fp <- as.character(files[[cohort]])
    if (!file.exists(fp)) next
    obj <- readRDS(fp)
    d <- olfml2b_part2_extract_clinical_frame(obj)
    rows[[cohort]] <- data.frame(
      cohort = cohort, n_total = nrow(d),
      n_tumor = sum(d$sample_context == "Tumor", na.rm = TRUE),
      n_normal = sum(d$sample_context == "Normal", na.rm = TRUE),
      n_unknown = sum(is.na(d$sample_context) | d$sample_context == "Unknown", na.rm = TRUE),
      expected_context = unique(as.character(d$expected_sample_context))[1L],
      analysis_role = role$analysis_role[match(cohort, role$cohort)],
      canonical_endpoint_family = unique(as.character(d$canonical_endpoint_family))[1L],
      source_endpoint_label = unique(as.character(d$source_endpoint_label))[1L],
      normal_comparison_allowed = any(d$normal_comparison_allowed %in% TRUE, na.rm = TRUE),
      context_contract_status = ifelse(
        nrow(d) > 0L && all(d$sample_context == "Tumor") && !any(d$normal_comparison_allowed %in% TRUE, na.rm = TRUE),
        "PASS_TUMOR_ONLY", "FAIL_CONTEXT_CONTRACT"
      ),
      stringsAsFactors = FALSE
    )
  }
  out <- olfml2b_bind_rows(rows)
  pdir <- olfml2b_part_paths(ctx, "Part2")
  olfml2b_atomic_write_csv(out, file.path(pdir$tables, "Part2_GEO_sample_context_audit.csv"))
  required <- c("GSE62254", "GSE15459", "GSE26253", "GSE84437", "GSE147163")
  if (!all(required %in% out$cohort) || any(out$context_contract_status != "PASS_TUMOR_ONLY") || any(out$n_normal > 0L) || any(out$n_unknown > 0L)) {
    stop("Part2 GEO tumor-only sample-context contract failed.", call. = FALSE)
  }
  out
}

.olfml2b_part2_v110_core <- run_olfml2b_part2_geo
run_olfml2b_part2_geo <- function(ctx = NULL, include_gse84437 = TRUE) {
  if (is.null(ctx)) ctx <- olfml2b_load_context()
  idx <- .olfml2b_part2_v110_core(ctx = ctx, include_gse84437 = include_gse84437)
  for (cohort in names(idx$cohort_files)) {
    fp <- as.character(idx$cohort_files[[cohort]])
    if (!file.exists(fp)) next
    obj <- readRDS(fp)
    d <- olfml2b_part2_extract_clinical_frame(obj)
    d <- olfml2b_part2_apply_sample_context(d, cohort)
    obj$sample_metadata <- d
    obj$clinical <- d
    obj$sample_context_contract_version <- OLFML2B_PART2_SAMPLE_CONTEXT_VERSION
    obj$expected_sample_context <- "TUMOR_ONLY"
    obj$normal_comparison_allowed <- FALSE
    obj$source_endpoint_label <- unique(d$source_endpoint_label)[1L]
    obj$canonical_endpoint_family <- unique(d$canonical_endpoint_family)[1L]
    obj$release_contract_version <- OLFML2B_PART2_PATCH_VERSION
    olfml2b_atomic_save_rds(obj, fp)
    pdir <- olfml2b_part_paths(ctx, "Part2")
    olfml2b_atomic_write_csv(d, file.path(pdir$tables, paste0("Part2_", cohort, "_sample_metadata.csv")))
    olfml2b_atomic_write_csv(d, file.path(ctx$dirs$derived_clinical, paste0(cohort, "_sample_metadata_harmonized.csv")))
  }
  context_fp <- file.path(ctx$dirs$objects, "Part2_GSE147163_context_only.rds")
  if (file.exists(context_fp)) {
    obj <- readRDS(context_fp)
    d <- olfml2b_part2_extract_clinical_frame(obj)
    d <- olfml2b_part2_apply_sample_context(d, "GSE147163")
    obj$sample_metadata <- d; obj$clinical <- d
    obj$sample_context_contract_version <- OLFML2B_PART2_SAMPLE_CONTEXT_VERSION
    obj$normal_comparison_allowed <- FALSE
    olfml2b_atomic_save_rds(obj, context_fp)
    pdir <- olfml2b_part_paths(ctx, "Part2")
    olfml2b_atomic_write_csv(d, file.path(pdir$tables, "Part2_GSE147163_context_only_sample_metadata.csv"))
  }
  idx <- readRDS(file.path(ctx$dirs$objects, "Part2_GEO_bulk_index.rds"))
  idx$release_contract_version <- OLFML2B_PART2_PATCH_VERSION
  idx$sample_context_contract_version <- OLFML2B_PART2_SAMPLE_CONTEXT_VERSION
  idx$geo_role_contract <- olfml2b_part2_geo_role_contract()
  audit <- olfml2b_part2_sample_context_audit(idx, ctx)
  idx$sample_context_audit <- audit
  olfml2b_atomic_save_rds(idx, file.path(ctx$dirs$objects, "Part2_GEO_bulk_index.rds"))
  invisible(idx)
}
