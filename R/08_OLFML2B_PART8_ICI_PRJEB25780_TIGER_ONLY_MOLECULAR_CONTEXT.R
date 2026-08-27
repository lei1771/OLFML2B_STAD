############################################################
## OLFML2B-STAD | Part8 | PRJEB25780 anti-PD-1 molecular-context boundary analysis
## TIGER-only expression/response implementation; no ENA manifest dependency
## Version: v1.0.0_20260723_OLFML2B_DIRECT_FORK_FIRTH_CONTEXT_BOUNDARY
############################################################

PART8_ICI_TIGER_VERSION <- "v1.0.0_20260723_OLFML2B_DIRECT_FORK_FIRTH_CONTEXT_BOUNDARY"

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

.p8_timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

.p8_msg <- function(..., log_file = NULL) {
    txt <- paste0("[", .p8_timestamp(), "] [OLFML2B-P8] ", paste0(..., collapse = ""))
    message(txt)
    if (!is.null(log_file) && nzchar(log_file)) {
        cat(txt, "\n", file = log_file, append = TRUE)
    }
    invisible(txt)
}

.p8_stop <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

.p8_clean <- function(x) trimws(as.character(x))

.p8_write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, ".tmp")
    utils::write.csv(x, tmp, row.names = FALSE, na = "")
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
    invisible(path)
}

.p8_save_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, ".tmp")
    saveRDS(x, tmp, compress = "xz")
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
    invisible(path)
}

.p8_bind_rows <- function(x) {
    x <- x[!vapply(x, is.null, logical(1))]
    x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
    if (!length(x)) return(data.frame())
    all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    x <- lapply(x, function(d) {
        miss <- setdiff(all_names, names(d))
        for (m in miss) d[[m]] <- NA
        d[, all_names, drop = FALSE]
    })
    do.call(rbind, x)
}

.p8_safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

.p8_z <- function(x) {
    x <- as.numeric(x)
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
    (x - m) / s
}

.p8_detect_delim <- function(path) {
    if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) return(",")
    "\t"
}

.p8_read_table_any <- function(path, nrows = -1L) {
    sep <- .p8_detect_delim(path)
    utils::read.table(
        path,
        sep = sep,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        nrows = nrows
    )
}

.p8_read_any <- function(path) {
    if (grepl("\\.rds$|\\.Rds$|\\.RDS$", path)) {
        return(readRDS(path))
    }
    .p8_read_table_any(path)
}

.p8_is_expression_object <- function(obj) {
    is.data.frame(obj) &&
        any(toupper(names(obj)) %in% c("GENE_SYMBOL", "GENE", "GENE_SYMBOLS", "SYMBOL")) &&
        ncol(obj) >= 20L &&
        nrow(obj) >= 5000L
}

.p8_is_clinical_object <- function(obj) {
    is.data.frame(obj) &&
        all(c("sample_id", "patient_name") %in% names(obj)) &&
        any(grepl("response|RECIST|benefit", names(obj), ignore.case = TRUE)) &&
        ncol(obj) >= 5L
}

.p8_find_tiger_files <- function(tiger_dir) {
    if (!dir.exists(tiger_dir)) .p8_stop("TIGER directory not found: ", tiger_dir)

    preferred_expr <- file.path(tiger_dir, "STAD-PRJEB25780.Response.Rds")
    preferred_clin <- file.path(tiger_dir, "STAD-PRJEB25780.Response (1).Rds")

    expr_file <- if (file.exists(preferred_expr)) preferred_expr else NA_character_
    clin_file <- if (file.exists(preferred_clin)) preferred_clin else NA_character_

    candidates <- list.files(
        tiger_dir,
        pattern = "\\.rds$|\\.Rds$|\\.RDS$|\\.tsv$|\\.txt$|\\.csv$|\\.tsv\\.gz$|\\.txt\\.gz$|\\.csv\\.gz$",
        full.names = TRUE,
        recursive = TRUE
    )
    candidates <- candidates[!grepl("Part8|audit|summary|wilcoxon|logistic|fisher|input|figure|claim", basename(candidates), ignore.case = TRUE)]

    file_audit <- lapply(candidates, function(f) {
        obj <- tryCatch(.p8_read_any(f), error = function(e) e)
        if (inherits(obj, "error")) {
            return(data.frame(
                file = f, basename = basename(f), readable = FALSE,
                class = NA_character_, dim = NA_character_,
                is_expression = FALSE, is_clinical = FALSE,
                error = conditionMessage(obj), stringsAsFactors = FALSE
            ))
        }
        data.frame(
            file = f,
            basename = basename(f),
            readable = TRUE,
            class = paste(class(obj), collapse = ";"),
            dim = if (!is.null(dim(obj))) paste(dim(obj), collapse = " x ") else NA_character_,
            is_expression = .p8_is_expression_object(obj),
            is_clinical = .p8_is_clinical_object(obj),
            error = NA_character_,
            stringsAsFactors = FALSE
        )
    })
    file_audit <- .p8_bind_rows(file_audit)

    if (is.na(expr_file) && nrow(file_audit)) {
        hit <- file_audit$file[file_audit$is_expression %in% TRUE]
        if (length(hit)) expr_file <- hit[1]
    }
    if (is.na(clin_file) && nrow(file_audit)) {
        hit <- file_audit$file[file_audit$is_clinical %in% TRUE]
        if (length(hit)) clin_file <- hit[1]
    }

    if (is.na(expr_file) || !file.exists(expr_file)) .p8_stop("Could not identify TIGER expression file in: ", tiger_dir)
    if (is.na(clin_file) || !file.exists(clin_file)) .p8_stop("Could not identify TIGER clinical/response file in: ", tiger_dir)

    list(expression_file = expr_file, clinical_file = clin_file, file_audit = file_audit)
}

.p8_gene_column <- function(expr) {
    candidates <- c("GENE_SYMBOL", "gene", "Gene", "gene_symbol", "GeneSymbol", "symbol", "SYMBOL")
    hit <- candidates[candidates %in% names(expr)]
    if (!length(hit)) .p8_stop("No gene symbol column found. First columns: ", paste(head(names(expr), 20), collapse = ", "))
    hit[1]
}

.p8_prepare_expression <- function(expr_raw, transform = c("auto", "none", "log2"), log_file = NULL) {
    transform <- match.arg(transform)
    expr <- as.data.frame(expr_raw, stringsAsFactors = FALSE, check.names = FALSE)
    gene_col <- .p8_gene_column(expr)
    sample_cols <- setdiff(names(expr), gene_col)

    mat <- as.matrix(expr[, sample_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    genes <- toupper(.p8_clean(expr[[gene_col]]))

    ok_gene <- !is.na(genes) & nzchar(genes) & genes != "NA"
    mat <- mat[ok_gene, , drop = FALSE]
    genes <- genes[ok_gene]

    finite_vals <- as.numeric(mat[is.finite(mat)])
    transform_used <- "none"
    if (transform == "log2" || (transform == "auto" && length(finite_vals) && stats::quantile(finite_vals, 0.95, na.rm = TRUE) > 50)) {
        mat <- log2(pmax(mat, 0) + 1)
        transform_used <- "log2_x_plus_1"
    }

    # Collapse duplicate symbols by row mean. Exact OLFML2B row count is still audited separately.
    split_idx <- split(seq_along(genes), genes)
    collapsed <- do.call(rbind, lapply(split_idx, function(idx) {
        if (length(idx) == 1L) return(mat[idx, , drop = FALSE])
        colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    }))
    rownames(collapsed) <- names(split_idx)

    .p8_msg("Prepared expression matrix | genes=", nrow(collapsed), " | samples=", ncol(collapsed), " | transform=", transform_used, log_file = log_file)

    list(
        matrix = collapsed,
        gene_col = gene_col,
        sample_ids = colnames(collapsed),
        transform_used = transform_used,
        original_gene_symbols = genes
    )
}

.p8_clean_response <- function(clin) {
    clin$response_clean <- toupper(.p8_clean(clin$response %||% NA_character_))
    clin$response_NR_clean <- toupper(.p8_clean(clin$response_NR %||% NA_character_))

    clin$response_NR_clean <- ifelse(clin$response_clean %in% c("CR", "PR"), "R", clin$response_NR_clean)
    clin$response_NR_clean <- ifelse(clin$response_clean %in% c("SD", "PD"), "NR", clin$response_NR_clean)
    clin$response_NR_clean <- ifelse(clin$response_NR_clean %in% c("N", "NR", "NONRESPONDER", "NON-RESPONDER", "NON_RESPONDER"), "NR", clin$response_NR_clean)
    clin$response_NR_clean <- ifelse(clin$response_NR_clean %in% c("R", "RESPONDER"), "R", clin$response_NR_clean)

    clin$response_evaluable <- clin$response_NR_clean %in% c("R", "NR")
    clin$response_bin_R <- ifelse(clin$response_NR_clean == "R", 1L, ifelse(clin$response_NR_clean == "NR", 0L, NA_integer_))
    clin
}

.p8_classify_tiger_sample_type <- function(clin) {
    pn <- .p8_clean(clin$patient_name)
    clin$sample_type_tiger <- ifelse(grepl("-N$", pn, perl = TRUE), "normal_like", "tumor")
    clin$patient_id_tiger <- sub("-N$", "", pn, perl = TRUE)
    clin$patient_name_normal_like <- clin$sample_type_tiger == "normal_like"
    clin
}

.p8_clinical_completeness <- function(clin) {
    data.frame(
        column = names(clin),
        n_non_missing = vapply(clin, function(x) sum(!is.na(x) & nzchar(.p8_clean(x))), integer(1)),
        n_missing_or_blank = vapply(clin, function(x) sum(is.na(x) | !nzchar(.p8_clean(x))), integer(1)),
        stringsAsFactors = FALSE
    )
}

.p8_make_sample_audit <- function(expr_prep, clin, merged) {
    data.frame(
        audit_item = c(
            "part8_version",
            "data_source_policy",
            "expression_sample_count",
            "clinical_row_count",
            "expression_clinical_overlap",
            "tiger_tumor_rows",
            "tiger_normal_like_rows_excluded",
            "tumor_response_evaluable_rows",
            "tumor_responders_CR_PR",
            "tumor_nonresponders_SD_PD",
            "unique_tumor_type_values",
            "unique_dataset_group_values",
            "unique_seq_type_values",
            "expression_transform_used"
        ),
        value = c(
            PART8_ICI_TIGER_VERSION,
            "TIGER_ONLY; no ENA manifest is read or used",
            length(expr_prep$sample_ids),
            nrow(clin),
            sum(clin$sample_id %in% expr_prep$sample_ids),
            sum(merged$sample_type_tiger == "tumor", na.rm = TRUE),
            sum(merged$sample_type_tiger == "normal_like", na.rm = TRUE),
            sum(merged$sample_type_tiger == "tumor" & merged$response_evaluable, na.rm = TRUE),
            sum(merged$sample_type_tiger == "tumor" & merged$response_NR_clean == "R", na.rm = TRUE),
            sum(merged$sample_type_tiger == "tumor" & merged$response_NR_clean == "NR", na.rm = TRUE),
            paste(unique(.p8_clean(clin$tumor_type %||% NA_character_)), collapse = "; "),
            paste(unique(.p8_clean(clin$dataset_group %||% NA_character_)), collapse = "; "),
            paste(unique(.p8_clean(clin$seq_type %||% NA_character_)), collapse = "; "),
            expr_prep$transform_used
        ),
        stringsAsFactors = FALSE
    )
}

.p8_target_expression <- function(expr_prep, target_gene) {
    target <- toupper(target_gene)
    raw_count <- sum(expr_prep$original_gene_symbols == target, na.rm = TRUE)
    detected <- target %in% rownames(expr_prep$matrix)
    target_values <- if (detected) as.numeric(expr_prep$matrix[target, ]) else rep(NA_real_, ncol(expr_prep$matrix))
    out <- data.frame(
        sample_id = expr_prep$sample_ids,
        target_gene = target_gene,
        target_expr = target_values,
        stringsAsFactors = FALSE
    )
    mapping <- data.frame(
        target_gene = target_gene,
        detected_after_symbol_collapse = detected,
        exact_rows_before_collapse = raw_count,
        rows_after_collapse = ifelse(detected, 1L, 0L),
        interpretation = if (detected && raw_count == 1L) "EXACT_SINGLE_SYMBOL_ROW_NO_COLLAPSING_REQUIRED" else if (detected) "DETECTED_DUPLICATE_SYMBOLS_COLLAPSED_BY_MEAN" else "NOT_DETECTED",
        stringsAsFactors = FALSE
    )
    list(values = out, mapping = mapping)
}

.p8_module_gene_sets <- function() {
    list(
        CAF_Core = c("FAP", "PDGFRA", "PDGFRB", "ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM"),
        TGFb_Response = c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF", "INHBA", "THBS1", "TGFBI", "POSTN"),
        ECM_Remodeling = c("COL1A1", "COL1A2", "COL3A1", "FN1", "SPARC", "MMP2", "MMP9", "MMP11", "MMP14", "TIMP1", "LOX", "LOXL2", "ITGA5", "ITGB1"),
        EMT = c("VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "CDH2", "FN1", "ITGA5", "MMP2", "MMP9"),
        Angiogenesis_Endothelial = c("VEGFA", "KDR", "FLT1", "PECAM1", "VWF", "ANGPT2", "ENG", "CDH5", "ESAM", "CLDN5"),
        Myeloid_Macrophage = c("CD68", "CSF1R", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "FCGR3A", "MSR1", "TYROBP"),
        CD8_Cytotoxic = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "IFNG", "CXCL9", "CXCL10"),
        IFNG_Response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP5"),
        Antigen_Presentation = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "TAP1", "TAP2"),
        Immune_Checkpoint = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1"),
        Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM4", "MCM6", "CCNB1", "CDK1"),
        Inflamed_5gene_exploratory = c("CD274", "CXCL9", "CXCL10", "GZMB", "IFNG")
    )
}

.p8_score_modules <- function(expr_mat, gene_sets) {
    sample_ids <- colnames(expr_mat)
    scores <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
    coverage_rows <- list()

    for (nm in names(gene_sets)) {
        requested <- unique(toupper(gene_sets[[nm]]))
        present <- intersect(requested, rownames(expr_mat))
        missing <- setdiff(requested, present)
        coverage_rows[[nm]] <- data.frame(
            module = nm,
            n_requested_genes = length(requested),
            n_present_genes = length(present),
            coverage_fraction = ifelse(length(requested) > 0, length(present) / length(requested), NA_real_),
            present_genes = paste(present, collapse = ";"),
            missing_genes = paste(missing, collapse = ";"),
            status = ifelse(length(present) >= 3L, "OK", "LOW_COVERAGE"),
            stringsAsFactors = FALSE
        )
        if (length(present) >= 1L) {
            sub <- expr_mat[present, , drop = FALSE]
            z <- t(apply(sub, 1L, .p8_z))
            scores[[nm]] <- colMeans(z, na.rm = TRUE)
        } else {
            scores[[nm]] <- NA_real_
        }
    }

    # Pre-specified composite axes built only from base module scores.
    if (all(c("CAF_Core", "TGFb_Response") %in% names(scores))) {
        scores$CAF_TGFb_axis <- rowMeans(cbind(.p8_z(scores$CAF_Core), .p8_z(scores$TGFb_Response)), na.rm = TRUE)
    }
    if (all(c("CAF_Core", "TGFb_Response", "ECM_Remodeling") %in% names(scores))) {
        scores$CAF_TGFb_ECM_axis <- rowMeans(cbind(.p8_z(scores$CAF_Core), .p8_z(scores$TGFb_Response), .p8_z(scores$ECM_Remodeling)), na.rm = TRUE)
    }
    if (all(c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "EMT") %in% names(scores))) {
        scores$Stromal_Remodeling_axis <- rowMeans(cbind(.p8_z(scores$CAF_Core), .p8_z(scores$TGFb_Response), .p8_z(scores$ECM_Remodeling), .p8_z(scores$EMT)), na.rm = TRUE)
    }

    composite_cov <- data.frame(
        module = c("CAF_TGFb_axis", "CAF_TGFb_ECM_axis", "Stromal_Remodeling_axis"),
        n_requested_genes = NA_integer_, n_present_genes = NA_integer_, coverage_fraction = NA_real_,
        present_genes = "COMPOSITE_OF_PREDEFINED_MODULE_SCORES",
        missing_genes = NA_character_, status = "COMPOSITE",
        stringsAsFactors = FALSE
    )

    list(scores = scores, coverage = rbind(.p8_bind_rows(coverage_rows), composite_cov))
}

.p8_auc_rank <- function(feature, response_bin_R) {
    ok <- is.finite(feature) & response_bin_R %in% c(0L, 1L)
    feature <- feature[ok]
    y <- response_bin_R[ok]
    n1 <- sum(y == 1L)
    n0 <- sum(y == 0L)
    if (n1 < 1L || n0 < 1L) return(NA_real_)
    r <- rank(feature, ties.method = "average")
    auc <- (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    as.numeric(auc)
}

.p8_group_summary_one <- function(d, feature) {
    rows <- lapply(c("NR", "R"), function(g) {
        x <- d[[feature]][d$response_NR_clean == g]
        x <- x[is.finite(x)]
        data.frame(
            feature = feature,
            response_NR_clean = g,
            n = length(x),
            median = ifelse(length(x), median(x), NA_real_),
            mean = ifelse(length(x), mean(x), NA_real_),
            sd = ifelse(length(x) > 1L, stats::sd(x), NA_real_),
            q1 = ifelse(length(x), as.numeric(stats::quantile(x, 0.25)), NA_real_),
            q3 = ifelse(length(x), as.numeric(stats::quantile(x, 0.75)), NA_real_),
            stringsAsFactors = FALSE
        )
    })
    .p8_bind_rows(rows)
}

.p8_wilcox_one <- function(d, feature) {
    x <- d[[feature]]
    ok <- is.finite(x) & d$response_NR_clean %in% c("R", "NR")
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 10L || length(unique(dd$response_NR_clean)) < 2L) {
        return(data.frame(feature = feature, n = nrow(dd), statistic = NA_real_, p_value = NA_real_, median_R = NA_real_, median_NR = NA_real_, delta_median_R_minus_NR = NA_real_, status = "INSUFFICIENT_GROUPS", stringsAsFactors = FALSE))
    }
    wt <- tryCatch(suppressWarnings(stats::wilcox.test(dd[[feature]] ~ dd$response_NR_clean, exact = FALSE)), error = function(e) NULL)
    med_R <- median(dd[[feature]][dd$response_NR_clean == "R"], na.rm = TRUE)
    med_NR <- median(dd[[feature]][dd$response_NR_clean == "NR"], na.rm = TRUE)
    data.frame(
        feature = feature,
        n = nrow(dd),
        n_R = sum(dd$response_NR_clean == "R"),
        n_NR = sum(dd$response_NR_clean == "NR"),
        statistic = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
        p_value = if (!is.null(wt)) wt$p.value else NA_real_,
        median_R = med_R,
        median_NR = med_NR,
        delta_median_R_minus_NR = med_R - med_NR,
        direction = ifelse(is.finite(med_R - med_NR) && med_R > med_NR, "HIGHER_IN_R", ifelse(is.finite(med_R - med_NR) && med_R < med_NR, "HIGHER_IN_NR", "NO_MEDIAN_DIFFERENCE")),
        status = if (!is.null(wt)) "OK" else "TEST_FAILED",
        stringsAsFactors = FALSE
    )
}

.p8_logistic_one <- function(d, feature) {
    x <- d[[feature]]
    ok <- is.finite(x) & d$response_bin_R %in% c(0L, 1L)
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 10L || length(unique(dd$response_bin_R)) < 2L || stats::sd(dd[[feature]], na.rm = TRUE) <= 0) {
        return(data.frame(feature = feature, n = nrow(dd), estimate = NA_real_, se = NA_real_, z = NA_real_, p_value = NA_real_, OR_per_1SD = NA_real_, OR_lcl = NA_real_, OR_ucl = NA_real_, auc_for_R_high_value = NA_real_, status = "INSUFFICIENT_DATA", stringsAsFactors = FALSE))
    }
    dd$feature_z <- .p8_z(dd[[feature]])
    fit <- tryCatch(stats::glm(response_bin_R ~ feature_z, data = dd, family = stats::binomial()), error = function(e) NULL)
    if (is.null(fit)) {
        return(data.frame(feature = feature, n = nrow(dd), estimate = NA_real_, se = NA_real_, z = NA_real_, p_value = NA_real_, OR_per_1SD = NA_real_, OR_lcl = NA_real_, OR_ucl = NA_real_, auc_for_R_high_value = .p8_auc_rank(dd[[feature]], dd$response_bin_R), status = "MODEL_FAILED", stringsAsFactors = FALSE))
    }
    cf <- summary(fit)$coefficients
    if (!"feature_z" %in% rownames(cf)) {
        return(data.frame(feature = feature, n = nrow(dd), estimate = NA_real_, se = NA_real_, z = NA_real_, p_value = NA_real_, OR_per_1SD = NA_real_, OR_lcl = NA_real_, OR_ucl = NA_real_, auc_for_R_high_value = .p8_auc_rank(dd[[feature]], dd$response_bin_R), status = "COEFFICIENT_MISSING", stringsAsFactors = FALSE))
    }
    est <- cf["feature_z", "Estimate"]
    se <- cf["feature_z", "Std. Error"]
    data.frame(
        feature = feature,
        n = nrow(dd),
        n_R = sum(dd$response_bin_R == 1L),
        n_NR = sum(dd$response_bin_R == 0L),
        estimate = est,
        se = se,
        z = cf["feature_z", "z value"],
        p_value = cf["feature_z", "Pr(>|z|)"],
        OR_per_1SD = exp(est),
        OR_lcl = exp(est - 1.96 * se),
        OR_ucl = exp(est + 1.96 * se),
        auc_for_R_high_value = .p8_auc_rank(dd[[feature]], dd$response_bin_R),
        status = "OK",
        stringsAsFactors = FALSE
    )
}

.p8_fisher_one <- function(d, feature) {
    x <- d[[feature]]
    ok <- is.finite(x) & d$response_NR_clean %in% c("R", "NR")
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 10L || length(unique(dd$response_NR_clean)) < 2L || stats::sd(dd[[feature]], na.rm = TRUE) <= 0) {
        return(data.frame(feature = feature, n = nrow(dd), cutoff = NA_real_, low_NR = NA_integer_, low_R = NA_integer_, high_NR = NA_integer_, high_R = NA_integer_, odds_ratio = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, status = "INSUFFICIENT_DATA", stringsAsFactors = FALSE))
    }
    cutoff <- median(dd[[feature]], na.rm = TRUE)
    dd$feature_high <- ifelse(dd[[feature]] >= cutoff, "High", "Low")
    tab <- table(factor(dd$feature_high, levels = c("Low", "High")), factor(dd$response_NR_clean, levels = c("NR", "R")))
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    data.frame(
        feature = feature,
        n = nrow(dd),
        cutoff = cutoff,
        low_NR = unname(tab["Low", "NR"]),
        low_R = unname(tab["Low", "R"]),
        high_NR = unname(tab["High", "NR"]),
        high_R = unname(tab["High", "R"]),
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        p_value = if (!is.null(ft)) ft$p.value else NA_real_,
        conf_low = if (!is.null(ft)) ft$conf.int[1] else NA_real_,
        conf_high = if (!is.null(ft)) ft$conf.int[2] else NA_real_,
        status = if (!is.null(ft)) "OK" else "TEST_FAILED",
        stringsAsFactors = FALSE
    )
}

.p8_analyze_features <- function(d, features) {
    features <- features[features %in% names(d)]
    features <- features[vapply(features, function(f) any(is.finite(d[[f]])), logical(1))]
    summary <- .p8_bind_rows(lapply(features, function(f) .p8_group_summary_one(d, f)))
    wilcox <- .p8_bind_rows(lapply(features, function(f) .p8_wilcox_one(d, f)))
    logistic <- .p8_bind_rows(lapply(features, function(f) .p8_logistic_one(d, f)))
    fisher <- .p8_bind_rows(lapply(features, function(f) .p8_fisher_one(d, f)))
    if (nrow(wilcox)) wilcox$p_fdr <- stats::p.adjust(wilcox$p_value, method = "BH")
    if (nrow(logistic)) logistic$p_fdr <- stats::p.adjust(logistic$p_value, method = "BH")
    if (nrow(fisher)) fisher$p_fdr <- stats::p.adjust(fisher$p_value, method = "BH")
    list(summary = summary, wilcox = wilcox, logistic = logistic, fisher = fisher)
}

.p8_correlate_target_modules <- function(d, target_col, module_cols) {
    rows <- lapply(module_cols, function(m) {
        if (!m %in% names(d)) return(NULL)
        ok <- is.finite(d[[target_col]]) & is.finite(d[[m]])
        dd <- d[ok, , drop = FALSE]
        if (nrow(dd) < 10L || stats::sd(dd[[target_col]]) <= 0 || stats::sd(dd[[m]]) <= 0) {
            return(data.frame(module = m, n = nrow(dd), rho = NA_real_, p_value = NA_real_, status = "INSUFFICIENT_DATA", stringsAsFactors = FALSE))
        }
        ct <- tryCatch(suppressWarnings(stats::cor.test(dd[[target_col]], dd[[m]], method = "spearman", exact = FALSE)), error = function(e) NULL)
        data.frame(
            module = m,
            n = nrow(dd),
            rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
            p_value = if (!is.null(ct)) ct$p.value else NA_real_,
            status = if (!is.null(ct)) "OK" else "TEST_FAILED",
            stringsAsFactors = FALSE
        )
    })
    out <- .p8_bind_rows(rows)
    if (nrow(out)) out$p_fdr <- stats::p.adjust(out$p_value, method = "BH")
    out
}

.p8_combination_test <- function(d, target_col, axis_col = "CAF_TGFb_ECM_axis") {
    if (!all(c(target_col, axis_col) %in% names(d))) return(data.frame())
    ok <- is.finite(d[[target_col]]) & is.finite(d[[axis_col]]) & d$response_NR_clean %in% c("R", "NR")
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 10L || length(unique(dd$response_NR_clean)) < 2L) return(data.frame())
    dd$target_high <- dd[[target_col]] >= median(dd[[target_col]], na.rm = TRUE)
    dd$axis_high <- dd[[axis_col]] >= median(dd[[axis_col]], na.rm = TRUE)
    dd$combined_high <- ifelse(dd$target_high & dd$axis_high, "Both_high", "Other")
    tab <- table(factor(dd$combined_high, levels = c("Other", "Both_high")), factor(dd$response_NR_clean, levels = c("NR", "R")))
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    data.frame(
        target_col = target_col,
        axis_col = axis_col,
        n = nrow(dd),
        other_NR = unname(tab["Other", "NR"]),
        other_R = unname(tab["Other", "R"]),
        both_high_NR = unname(tab["Both_high", "NR"]),
        both_high_R = unname(tab["Both_high", "R"]),
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        p_value = if (!is.null(ft)) ft$p.value else NA_real_,
        conf_low = if (!is.null(ft)) ft$conf.int[1] else NA_real_,
        conf_high = if (!is.null(ft)) ft$conf.int[2] else NA_real_,
        status = if (!is.null(ft)) "OK" else "TEST_FAILED",
        stringsAsFactors = FALSE
    )
}

.p8_claim_ceiling <- function(sample_audit, target_mapping, wilcox, logistic, cor_tab) {
    target_w <- wilcox[wilcox$feature == "target_expr", , drop = FALSE]
    target_l <- logistic[logistic$feature == "target_expr", , drop = FALSE]
    target_positive <- nrow(target_w) && is.finite(target_w$p_value[1]) && target_w$p_value[1] < 0.05
    target_logit <- nrow(target_l) && is.finite(target_l$p_value[1]) && target_l$p_value[1] < 0.05

    stromal_features <- c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "EMT", "CAF_TGFb_axis", "CAF_TGFb_ECM_axis", "Stromal_Remodeling_axis")
    stromal_w <- wilcox[wilcox$feature %in% stromal_features & wilcox$status == "OK", , drop = FALSE]
    stromal_signal <- nrow(stromal_w) && any(is.finite(stromal_w$p_value) & stromal_w$p_value < 0.05, na.rm = TRUE)

    target_cor <- cor_tab[cor_tab$status == "OK" & cor_tab$module %in% stromal_features, , drop = FALSE]
    cor_signal <- nrow(target_cor) && any(is.finite(target_cor$p_value) & target_cor$p_value < 0.05, na.rm = TRUE)

    data.frame(
        claim = c(
            "TIGER-only data source",
            "TIGER internal tumor-only response analysis",
            "target gene mapping",
            "target single-gene anti-PD-1 response prediction",
            "stromal module anti-PD-1 response association",
            "target-stromal module association",
            "final Part8 claim ceiling"
        ),
        status = c(
            "PASS",
            ifelse(any(sample_audit$audit_item == "tumor_response_evaluable_rows" & as.integer(sample_audit$value) >= 30L), "PASS", "REVIEW"),
            ifelse(isTRUE(target_mapping$detected_after_symbol_collapse[1]), "PASS", "FAIL"),
            ifelse(target_positive || target_logit, "EXPLORATORY_SIGNAL", "NO_GO_AS_SINGLE_GENE_PREDICTOR"),
            ifelse(stromal_signal, "EXPLORATORY_SIGNAL", "NO_CONSISTENT_SIGNAL_YET"),
            ifelse(cor_signal, "SUPPORTS_TARGET_ASSOCIATED_STATE", "NO_STRONG_CORRELATION_SIGNAL"),
            ifelse(target_positive || target_logit, "EXPLORATORY_ONLY_NOT_CLINICAL_MODEL", "BOUNDARY_VALIDATION_NO_SINGLE_GENE_ICI_CLAIM")
        ),
        interpretation = c(
            "Part8 reads only TIGER expression and TIGER response annotation; no ENA manifest is used.",
            "Samples ending in -N in TIGER patient_name are excluded before response analysis.",
            target_mapping$interpretation[1],
            "OLFML2B is not promoted to an anti-PD-1 biomarker from this single cohort, irrespective of nominal significance.",
            "Pre-specified CAF/TGFb/ECM/immune modules are tested without post-hoc feature selection.",
            "Spearman correlations test whether OLFML2B tracks the stromal-remodeling state within the TIGER tumor subset.",
            "Use Part8 to define response-analysis boundaries and avoid overclaiming."
        ),
        stringsAsFactors = FALSE
    )
}

.p8_make_figures <- function(d, target_col, wilcox, cor_tab, fig_dir, target_gene) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    # Use the ORR endpoint labels introduced in v2.1 when available; keep the
    # legacy response_NR_clean fallback for downstream compatibility.
    response_x <- if ("response_ORR_group" %in% names(d)) "response_ORR_group" else "response_NR_clean"

    p_target <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[response_x]], y = .data[[target_col]])) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.75, size = 2) +
        ggplot2::labs(
            title = paste0("PRJEB25780 anti-PD-1 tumor-only subset | ", target_gene),
            subtitle = "TIGER-only analysis; samples ending with -N excluded; ORR = CR/PR vs SD/PD",
            x = "ORR response group",
            y = paste0(target_gene, " expression, TIGER processed")
        ) +
        ggplot2::theme_classic(base_size = 12)
    ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_OLFML2B_by_response_TIGER_only.png"), p_target, width = 5.5, height = 4.5, dpi = 300)
    ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_OLFML2B_by_response_TIGER_only.pdf"), p_target, width = 5.5, height = 4.5)

    if (is.data.frame(wilcox) && nrow(wilcox)) {
        w <- wilcox[wilcox$status == "OK" & wilcox$feature != target_col, , drop = FALSE]
        if ("endpoint" %in% names(w)) {
            w <- w[w$endpoint == "ORR", , drop = FALSE]
        }

        # v2.0 used delta_median_R_minus_NR; v2.1 endpoint-aware tables use
        # delta_median_positive_minus_negative. Standardize the plotting field
        # before entering ggplot so missing-column evaluation cannot abort.
        if (!"delta_median_R_minus_NR" %in% names(w) && "delta_median_positive_minus_negative" %in% names(w)) {
            w$delta_median_R_minus_NR <- w$delta_median_positive_minus_negative
        }

        if (nrow(w) && "delta_median_R_minus_NR" %in% names(w)) {
            w <- w[is.finite(w$delta_median_R_minus_NR), , drop = FALSE]
        }

        if (nrow(w)) {
            p_for_size <- if ("p_fdr" %in% names(w)) w$p_fdr else w$p_value
            w$neg_log10_p <- -log10(pmax(p_for_size, .Machine$double.xmin))
            w$neg_log10_p[!is.finite(w$neg_log10_p)] <- NA_real_

            p_mod <- ggplot2::ggplot(w, ggplot2::aes(x = stats::reorder(feature, delta_median_R_minus_NR), y = delta_median_R_minus_NR)) +
                ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
                ggplot2::geom_point(ggplot2::aes(size = neg_log10_p), alpha = 0.85) +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    title = "Pre-specified module differences by anti-PD-1 ORR",
                    subtitle = "Delta median = CR/PR responder minus SD/PD non-responder",
                    x = NULL,
                    y = "Delta median score (R - NR)",
                    size = "-log10(FDR/P)"
                ) +
                ggplot2::theme_classic(base_size = 11)
            ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_module_response_delta_TIGER_only.png"), p_mod, width = 7.2, height = 5.2, dpi = 300)
            ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_module_response_delta_TIGER_only.pdf"), p_mod, width = 7.2, height = 5.2)
        }
    }

    if (is.data.frame(cor_tab) && nrow(cor_tab)) {
        cc <- cor_tab[cor_tab$status == "OK", , drop = FALSE]
        if (nrow(cc)) {
            p_cor <- ggplot2::ggplot(cc, ggplot2::aes(x = stats::reorder(module, rho), y = rho)) +
                ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
                ggplot2::geom_point(ggplot2::aes(size = -log10(pmax(p_value, .Machine$double.xmin))), alpha = 0.85) +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    title = paste0(target_gene, " correlation with pre-specified modules"),
                    subtitle = "Spearman correlation within TIGER tumor-only subset",
                    x = NULL,
                    y = "Spearman rho",
                    size = "-log10(P)"
                ) +
                ggplot2::theme_classic(base_size = 11)
            ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_OLFML2B_module_correlations_TIGER_only.png"), p_cor, width = 7.2, height = 5.2, dpi = 300)
            ggplot2::ggsave(file.path(fig_dir, "Part8_Fig_OLFML2B_module_correlations_TIGER_only.pdf"), p_cor, width = 7.2, height = 5.2)
        }
    }
    invisible(TRUE)
}



############################################################
## v2.1 molecular-context helper functions
############################################################

.p8_standard_patient_id <- function(x) {
    x <- .p8_clean(x)
    x <- sub("_RNA$", "", x, perl = TRUE)
    x <- sub("-N$", "", x, perl = TRUE)
    x
}

.p8_pick_col <- function(df, patterns, required = FALSE) {
    cn <- names(df)
    hits <- character(0)
    for (p in patterns) {
        hits <- c(hits, grep(p, cn, ignore.case = TRUE, value = TRUE, perl = TRUE))
    }
    hits <- unique(hits)
    if (!length(hits)) {
        if (isTRUE(required)) .p8_stop("Required column not found. Patterns: ", paste(patterns, collapse = "; "), ". Columns: ", paste(cn, collapse = ", "))
        return(NA_character_)
    }
    hits[1]
}

.p8_parse_numeric <- function(x) {
    suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", as.character(x))))
}

.p8_find_molecular_file <- function(root) {
    cache_dir <- file.path(root, "data/cache/ICI_PRJEB25780")
    candidates <- c(
        file.path(cache_dir, "PRJEB25780_Kim2018_Bagaev2021_clinical_curated.csv"),
        file.path(root, "config/PRJEB25780_Kim2018_Bagaev2021_clinical_curated.csv"),
        file.path(cache_dir, "Bagaev2021_Kim_PRJEB25780_clinical_rows_raw.csv"),
        file.path(cache_dir, "Bagaev2021_mmc5_hits_Gastric.csv")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) hit[1] else NA_character_
}

.p8_add_response_endpoints <- function(d) {
    d$response_clean <- toupper(.p8_clean(d$response_clean %||% d$response %||% NA_character_))
    d$response_ORR_group <- ifelse(
        d$response_clean %in% c("CR", "PR"),
        "R",
        ifelse(d$response_clean %in% c("SD", "PD"), "NR", NA_character_)
    )
    d$response_ORR_bin_R <- ifelse(d$response_ORR_group == "R", 1L, ifelse(d$response_ORR_group == "NR", 0L, NA_integer_))
    d$response_DCR_group <- ifelse(
        d$response_clean %in% c("CR", "PR", "SD"),
        "DCR",
        ifelse(d$response_clean %in% c("PD"), "PD", NA_character_)
    )
    d$response_DCR_bin <- ifelse(d$response_DCR_group == "DCR", 1L, ifelse(d$response_DCR_group == "PD", 0L, NA_integer_))
    d
}

.p8_group_summary_endpoint_one <- function(d, feature, group_col, positive_label, negative_label, endpoint_name) {
    rows <- lapply(c(negative_label, positive_label), function(g) {
        x <- d[[feature]][d[[group_col]] == g]
        x <- x[is.finite(x)]
        data.frame(
            endpoint = endpoint_name,
            feature = feature,
            group = g,
            n = length(x),
            median = ifelse(length(x), median(x), NA_real_),
            mean = ifelse(length(x), mean(x), NA_real_),
            sd = ifelse(length(x) > 1L, stats::sd(x), NA_real_),
            q1 = ifelse(length(x), as.numeric(stats::quantile(x, 0.25)), NA_real_),
            q3 = ifelse(length(x), as.numeric(stats::quantile(x, 0.75)), NA_real_),
            stringsAsFactors = FALSE
        )
    })
    .p8_bind_rows(rows)
}

.p8_wilcox_endpoint_one <- function(d, feature, group_col, positive_label, negative_label, endpoint_name) {
    ok <- is.finite(d[[feature]]) & d[[group_col]] %in% c(positive_label, negative_label)
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 8L || length(unique(dd[[group_col]])) < 2L) {
        return(data.frame(endpoint = endpoint_name, feature = feature, n = nrow(dd), n_positive = sum(dd[[group_col]] == positive_label), n_negative = sum(dd[[group_col]] == negative_label), statistic = NA_real_, p_value = NA_real_, median_positive = NA_real_, median_negative = NA_real_, delta_median_positive_minus_negative = NA_real_, direction = NA_character_, status = "INSUFFICIENT_GROUPS", stringsAsFactors = FALSE))
    }
    dd[[group_col]] <- factor(dd[[group_col]], levels = c(negative_label, positive_label))
    wt <- tryCatch(suppressWarnings(stats::wilcox.test(dd[[feature]] ~ dd[[group_col]], exact = FALSE)), error = function(e) NULL)
    med_pos <- median(dd[[feature]][dd[[group_col]] == positive_label], na.rm = TRUE)
    med_neg <- median(dd[[feature]][dd[[group_col]] == negative_label], na.rm = TRUE)
    delta <- med_pos - med_neg
    data.frame(
        endpoint = endpoint_name,
        feature = feature,
        n = nrow(dd),
        n_positive = sum(dd[[group_col]] == positive_label),
        n_negative = sum(dd[[group_col]] == negative_label),
        statistic = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
        p_value = if (!is.null(wt)) wt$p.value else NA_real_,
        median_positive = med_pos,
        median_negative = med_neg,
        delta_median_positive_minus_negative = delta,
        direction = ifelse(is.finite(delta) && delta > 0, paste0("HIGHER_IN_", positive_label), ifelse(is.finite(delta) && delta < 0, paste0("HIGHER_IN_", negative_label), "NO_MEDIAN_DIFFERENCE")),
        status = if (!is.null(wt)) "OK" else "TEST_FAILED",
        stringsAsFactors = FALSE
    )
}

 .p8_logistic_endpoint_one <- function(d, feature, y_col, endpoint_name, positive_label = "positive") {
    ok <- is.finite(d[[feature]]) & d[[y_col]] %in% c(0L, 1L)
    dd <- d[ok, , drop = FALSE]
    base <- function(status, estimate = NA_real_, se = NA_real_, z = NA_real_, p_value = NA_real_, OR = NA_real_, lcl = NA_real_, ucl = NA_real_, note = NA_character_) {
        data.frame(
            endpoint = endpoint_name,
            feature = feature,
            n = nrow(dd),
            n_positive = sum(dd[[y_col]] == 1L, na.rm = TRUE),
            n_negative = sum(dd[[y_col]] == 0L, na.rm = TRUE),
            estimate = estimate,
            se = se,
            z = z,
            p_value = p_value,
            OR_per_1SD = OR,
            OR_lcl = lcl,
            OR_ucl = ucl,
            auc_for_positive_high_value = if (nrow(dd)) .p8_auc_rank(dd[[feature]], dd[[y_col]]) else NA_real_,
            positive_class = positive_label,
            status = status,
            inference_note = note,
            stringsAsFactors = FALSE
        )
    }

    if (nrow(dd) < 8L || length(unique(dd[[y_col]])) < 2L || stats::sd(dd[[feature]], na.rm = TRUE) <= 0) {
        return(base("INSUFFICIENT_DATA", note = "Not enough evaluable observations, outcome classes, or feature variation."))
    }

    dd$feature_z <- .p8_z(dd[[feature]])
    fit <- tryCatch(suppressWarnings(stats::glm(stats::as.formula(paste0(y_col, " ~ feature_z")), data = dd, family = stats::binomial())), error = function(e) NULL)
    if (is.null(fit)) {
        return(base("MODEL_FAILED", note = "Logistic model failed; use Wilcoxon/Fisher/AUC for interpretation."))
    }

    cf <- tryCatch(suppressWarnings(summary(fit)$coefficients), error = function(e) NULL)
    if (is.null(cf) || !"feature_z" %in% rownames(cf)) {
        return(base("COEFFICIENT_MISSING", note = "Logistic coefficient missing; use Wilcoxon/Fisher/AUC for interpretation."))
    }

    est <- cf["feature_z", "Estimate"]
    se <- cf["feature_z", "Std. Error"]
    zval <- cf["feature_z", "z value"]
    pval <- cf["feature_z", "Pr(>|z|)"]
    probs <- tryCatch(stats::fitted(fit), error = function(e) rep(NA_real_, nrow(dd)))
    unstable <- any(is.finite(probs) & (probs < 1e-6 | probs > 1 - 1e-6), na.rm = TRUE) ||
        !is.finite(est) || !is.finite(se) || se > 5 || abs(est) > 8

    if (isTRUE(unstable)) {
        return(base(
            "QUASI_COMPLETE_SEPARATION_OR_UNSTABLE",
            estimate = est,
            se = se,
            z = zval,
            p_value = NA_real_,
            OR = NA_real_,
            lcl = NA_real_,
            ucl = NA_real_,
            note = "Logistic OR/P suppressed because fitted probabilities or coefficients indicate quasi-complete separation/instability; interpret group counts, exact tests and AUC instead."
        ))
    }

    base(
        "OK",
        estimate = est,
        se = se,
        z = zval,
        p_value = pval,
        OR = exp(est),
        lcl = exp(est - 1.96 * se),
        ucl = exp(est + 1.96 * se),
        note = "Logistic OR per 1 SD is interpretable as exploratory effect size."
    )
}

.p8_fisher_endpoint_one <- function(d, feature, y_col, endpoint_name, positive_label = "positive") {
    ok <- is.finite(d[[feature]]) & d[[y_col]] %in% c(0L, 1L)
    dd <- d[ok, , drop = FALSE]
    if (nrow(dd) < 8L || length(unique(dd[[y_col]])) < 2L || stats::sd(dd[[feature]], na.rm = TRUE) <= 0) {
        return(data.frame(endpoint = endpoint_name, feature = feature, n = nrow(dd), cutoff = NA_real_, low_negative = NA_integer_, low_positive = NA_integer_, high_negative = NA_integer_, high_positive = NA_integer_, odds_ratio = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, positive_class = positive_label, status = "INSUFFICIENT_DATA", stringsAsFactors = FALSE))
    }
    cutoff <- median(dd[[feature]], na.rm = TRUE)
    dd$feature_high <- ifelse(dd[[feature]] >= cutoff, "High", "Low")
    tab <- table(factor(dd$feature_high, levels = c("Low", "High")), factor(dd[[y_col]], levels = c(0L, 1L)))
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    data.frame(
        endpoint = endpoint_name,
        feature = feature,
        n = nrow(dd),
        cutoff = cutoff,
        low_negative = unname(tab["Low", "0"]),
        low_positive = unname(tab["Low", "1"]),
        high_negative = unname(tab["High", "0"]),
        high_positive = unname(tab["High", "1"]),
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        p_value = if (!is.null(ft)) ft$p.value else NA_real_,
        conf_low = if (!is.null(ft)) ft$conf.int[1] else NA_real_,
        conf_high = if (!is.null(ft)) ft$conf.int[2] else NA_real_,
        positive_class = positive_label,
        status = if (!is.null(ft)) "OK" else "TEST_FAILED",
        stringsAsFactors = FALSE
    )
}

.p8_analyze_features_endpoint <- function(d, features, endpoint = c("ORR", "DCR")) {
    endpoint <- match.arg(endpoint)
    if (endpoint == "ORR") {
        group_col <- "response_ORR_group"; y_col <- "response_ORR_bin_R"; positive <- "R"; negative <- "NR"; positive_label <- "CR_PR_responder"
    } else {
        group_col <- "response_DCR_group"; y_col <- "response_DCR_bin"; positive <- "DCR"; negative <- "PD"; positive_label <- "CR_PR_SD_disease_control"
    }
    features <- features[features %in% names(d)]
    features <- features[vapply(features, function(f) any(is.finite(d[[f]])), logical(1))]
    summary <- .p8_bind_rows(lapply(features, function(f) .p8_group_summary_endpoint_one(d, f, group_col, positive, negative, endpoint)))
    wilcox <- .p8_bind_rows(lapply(features, function(f) .p8_wilcox_endpoint_one(d, f, group_col, positive, negative, endpoint)))
    logistic <- .p8_bind_rows(lapply(features, function(f) .p8_logistic_endpoint_one(d, f, y_col, endpoint, positive_label)))
    fisher <- .p8_bind_rows(lapply(features, function(f) .p8_fisher_endpoint_one(d, f, y_col, endpoint, positive_label)))
    if (nrow(wilcox)) wilcox$p_fdr <- stats::p.adjust(wilcox$p_value, method = "BH")
    if (nrow(logistic)) logistic$p_fdr <- stats::p.adjust(logistic$p_value, method = "BH")
    if (nrow(fisher)) fisher$p_fdr <- stats::p.adjust(fisher$p_value, method = "BH")
    list(summary = summary, wilcox = wilcox, logistic = logistic, fisher = fisher)
}

.p8_clean_molecular_context <- function(analysis_df, molecular_file = NA_character_) {
    if (is.na(molecular_file) || !file.exists(molecular_file)) {
        analysis_df$molecular_annotation_available <- FALSE
        return(list(
            data = analysis_df,
            merge_audit = data.frame(item = c("molecular_file", "status"), value = c(NA_character_, "NOT_AVAILABLE"), stringsAsFactors = FALSE),
            completeness = data.frame(),
            available = FALSE
        ))
    }
    mol <- utils::read.csv(molecular_file, stringsAsFactors = FALSE, check.names = FALSE)
    mol_pid_col <- .p8_pick_col(mol, c("^patient_id$", "patient", "patient_name", "sample", "case"), required = TRUE)
    mol$patient_id_part8 <- .p8_standard_patient_id(mol[[mol_pid_col]])
    mol <- mol[!duplicated(mol$patient_id_part8), , drop = FALSE]

    analysis_df$patient_id_part8 <- .p8_standard_patient_id(analysis_df$patient_id_tiger %||% analysis_df$patient_name)
    merged <- merge(analysis_df, mol, by = "patient_id_part8", all.x = TRUE, suffixes = c("", "_mol"))
    matched <- !is.na(merged[[mol_pid_col]])

    tcga_col <- .p8_pick_col(merged, c("^tcga_subtype$", "molecular.*subtype", "subtype"))
    msi_col <- .p8_pick_col(merged, c("^msi", "msi", "mmr"))
    ebv_col <- .p8_pick_col(merged, c("^ebv_status$", "ebv"))
    tmb_col <- .p8_pick_col(merged, c("^tmb_mut_mb$", "tmb", "mutation_load", "Total Mutation", "mutation"))
    pdl1_col <- .p8_pick_col(merged, c("^pdl1_cps$", "pd.?l1.*cps", "cps"))

    merged$MSI_H <- NA
    if (!is.na(msi_col)) {
        v <- toupper(.p8_clean(merged[[msi_col]]))
        merged$MSI_H <- ifelse(grepl("MSS|PROFICIENT|PMMR", v), FALSE, ifelse(grepl("MSI|DMMR|DEFIC|HIGH", v), TRUE, NA))
    }
    if (!is.na(tcga_col)) {
        v <- toupper(.p8_clean(merged[[tcga_col]]))
        tcga_msi <- ifelse(grepl("MSI", v), TRUE, ifelse(nzchar(v), FALSE, NA))
        merged$MSI_H <- ifelse(is.na(merged$MSI_H), tcga_msi, merged$MSI_H | (tcga_msi %in% TRUE))
    }

    merged$EBV_pos <- NA
    if (!is.na(ebv_col)) {
        v <- toupper(.p8_clean(merged[[ebv_col]]))
        merged$EBV_pos <- ifelse(grepl("NEG|NEGATIVE|0$", v), FALSE, ifelse(grepl("POS|POSITIVE|\\+|EBV_POS|^EBV$", v), TRUE, NA))
    }
    if (!is.na(tcga_col)) {
        v <- toupper(.p8_clean(merged[[tcga_col]]))
        tcga_ebv <- ifelse(grepl("EBV", v), TRUE, ifelse(nzchar(v), FALSE, NA))
        merged$EBV_pos <- ifelse(is.na(merged$EBV_pos), tcga_ebv, merged$EBV_pos | (tcga_ebv %in% TRUE))
    }

    merged$TMB_value <- NA_real_
    if (!is.na(tmb_col)) merged$TMB_value <- .p8_parse_numeric(merged[[tmb_col]])
    merged$PDL1_CPS <- NA_real_
    if (!is.na(pdl1_col)) merged$PDL1_CPS <- .p8_parse_numeric(merged[[pdl1_col]])

    merged$TMB_high_10 <- ifelse(is.finite(merged$TMB_value), merged$TMB_value >= 10, NA)
    merged$TMB_high_median <- ifelse(is.finite(merged$TMB_value), merged$TMB_value >= median(merged$TMB_value, na.rm = TRUE), NA)
    merged$PDL1_CPS_ge1 <- ifelse(is.finite(merged$PDL1_CPS), merged$PDL1_CPS >= 1, NA)
    merged$PDL1_CPS_ge5 <- ifelse(is.finite(merged$PDL1_CPS), merged$PDL1_CPS >= 5, NA)
    merged$PDL1_CPS_ge10 <- ifelse(is.finite(merged$PDL1_CPS), merged$PDL1_CPS >= 10, NA)

    merged$strong_ici_context_strict <- (
        (merged$MSI_H %in% TRUE) |
        (merged$EBV_pos %in% TRUE) |
        (merged$TMB_high_10 %in% TRUE)
    )
    merged$strong_ici_context_broad <- (
        (merged$MSI_H %in% TRUE) |
        (merged$EBV_pos %in% TRUE) |
        (merged$TMB_high_median %in% TRUE) |
        (merged$PDL1_CPS_ge5 %in% TRUE)
    )
    merged$MSS_EBVneg_TMBlow_like <- (
        (!is.na(merged$MSI_H) & !merged$MSI_H) &
        (!is.na(merged$EBV_pos) & !merged$EBV_pos) &
        (!is.na(merged$TMB_high_median) & !merged$TMB_high_median)
    )

    if ("Inflamed_5gene_exploratory" %in% names(merged)) {
        merged$Inflamed5_score_value <- suppressWarnings(as.numeric(merged$Inflamed_5gene_exploratory))
        merged$Inflamed5_high_median <- ifelse(is.finite(merged$Inflamed5_score_value), merged$Inflamed5_score_value >= median(merged$Inflamed5_score_value, na.rm = TRUE), NA)
        merged$Inflamed5_MSS_EBVneg_TMBlow_like <- (merged$MSS_EBVneg_TMBlow_like %in% TRUE) & (merged$Inflamed5_high_median %in% TRUE)
    } else {
        merged$Inflamed5_score_value <- NA_real_
        merged$Inflamed5_high_median <- NA
        merged$Inflamed5_MSS_EBVneg_TMBlow_like <- NA
    }

    marker_cols <- c("MSI_H", "EBV_pos", "TMB_value", "TMB_high_10", "TMB_high_median", "PDL1_CPS", "PDL1_CPS_ge1", "PDL1_CPS_ge5", "PDL1_CPS_ge10", "Inflamed5_score_value", "Inflamed5_high_median", "Inflamed5_MSS_EBVneg_TMBlow_like", "strong_ici_context_strict", "strong_ici_context_broad", "MSS_EBVneg_TMBlow_like")
    completeness <- data.frame(
        marker = marker_cols,
        n_non_missing = vapply(marker_cols, function(v) sum(!is.na(merged[[v]]) & (is.finite(merged[[v]]) | is.logical(merged[[v]]) | is.character(merged[[v]]))), integer(1)),
        stringsAsFactors = FALSE
    )
    merge_audit <- data.frame(
        item = c("molecular_file", "molecular_rows", "analysis_rows", "matched_by_patient_id", "unmatched_analysis_rows", "patient_id_column_in_molecular", "tcga_col", "msi_col", "ebv_col", "tmb_col", "pdl1_col", "inflamed5_signature_genes", "inflamed5_signature_status"),
        value = c(molecular_file, nrow(mol), nrow(analysis_df), sum(matched), sum(!matched), mol_pid_col, tcga_col, msi_col, ebv_col, tmb_col, pdl1_col, "CD274;CXCL9;CXCL10;GZMB;IFNG", ifelse("Inflamed_5gene_exploratory" %in% names(merged), "SCORED_EQUAL_WEIGHT_MEAN_Z_EXPLORATORY", "NOT_SCORED")),
        stringsAsFactors = FALSE
    )
    merged$molecular_annotation_available <- TRUE
    list(data = merged, merge_audit = merge_audit, completeness = completeness, available = TRUE)
}

.p8_fisher_binary_marker <- function(d, marker, y_col, endpoint) {
    dd <- d[!is.na(d[[marker]]) & d[[y_col]] %in% c(0L, 1L), , drop = FALSE]
    if (nrow(dd) < 8L || length(unique(dd[[marker]])) < 2L || length(unique(dd[[y_col]])) < 2L) {
        return(data.frame(
            endpoint = endpoint, marker = marker, n = nrow(dd),
            n_positive = sum(dd[[y_col]] == 1L, na.rm = TRUE),
            n_negative = sum(dd[[y_col]] == 0L, na.rm = TRUE),
            marker_FALSE_negative = NA_integer_, marker_FALSE_positive = NA_integer_,
            marker_TRUE_negative = NA_integer_, marker_TRUE_positive = NA_integer_,
            table_string = NA_character_, odds_ratio = NA_real_, p_value = NA_real_,
            zero_cell_flag = NA, status = "INSUFFICIENT_DATA",
            inference_note = "Not enough observations/classes for exact test.",
            stringsAsFactors = FALSE
        ))
    }
    mm <- as.logical(dd[[marker]])
    tab <- table(factor(mm, levels = c(FALSE, TRUE)), factor(dd[[y_col]], levels = c(0L, 1L)))
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    zero_cell <- any(tab == 0)
    data.frame(
        endpoint = endpoint,
        marker = marker,
        n = nrow(dd),
        n_positive = sum(dd[[y_col]] == 1L),
        n_negative = sum(dd[[y_col]] == 0L),
        marker_FALSE_negative = unname(tab["FALSE", "0"]),
        marker_FALSE_positive = unname(tab["FALSE", "1"]),
        marker_TRUE_negative = unname(tab["TRUE", "0"]),
        marker_TRUE_positive = unname(tab["TRUE", "1"]),
        table_string = paste(capture.output(print(tab)), collapse = " | "),
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        p_value = if (!is.null(ft)) ft$p.value else NA_real_,
        zero_cell_flag = zero_cell,
        status = if (!is.null(ft)) "OK" else "TEST_FAILED",
        inference_note = ifelse(zero_cell, "Zero cell present; report exact counts/P value and avoid overinterpreting OR magnitude.", "Exact Fisher test is interpretable as exploratory marker enrichment."),
        stringsAsFactors = FALSE
    )
}

.p8_numeric_marker_logistic <- function(d, marker, y_col, endpoint) {
    .p8_logistic_endpoint_one(d, marker, y_col, endpoint, positive_label = endpoint)
}

.p8_target_by_binary_context <- function(d, target_col, binary_markers) {
    .p8_bind_rows(lapply(binary_markers, function(m) {
        if (!m %in% names(d)) return(NULL)
        dd <- d[!is.na(d[[m]]) & is.finite(d[[target_col]]), , drop = FALSE]
        if (nrow(dd) < 8L || length(unique(dd[[m]])) < 2L) return(data.frame(context_marker = m, n = nrow(dd), status = "INSUFFICIENT_DATA", p_value = NA_real_, stringsAsFactors = FALSE))
        dd$grp <- factor(as.logical(dd[[m]]), levels = c(FALSE, TRUE))
        wt <- tryCatch(stats::wilcox.test(dd[[target_col]] ~ dd$grp, exact = FALSE), error = function(e) NULL)
        data.frame(context_marker = m, n = nrow(dd), n_FALSE = sum(dd$grp == FALSE), n_TRUE = sum(dd$grp == TRUE), median_FALSE = median(dd[[target_col]][dd$grp == FALSE], na.rm = TRUE), median_TRUE = median(dd[[target_col]][dd$grp == TRUE], na.rm = TRUE), delta_TRUE_minus_FALSE = median(dd[[target_col]][dd$grp == TRUE], na.rm = TRUE) - median(dd[[target_col]][dd$grp == FALSE], na.rm = TRUE), p_value = if (!is.null(wt)) wt$p.value else NA_real_, status = if (!is.null(wt)) "OK" else "TEST_FAILED", stringsAsFactors = FALSE)
    }))
}

.p8_target_numeric_correlations <- function(d, target_col, numeric_markers) {
    .p8_bind_rows(lapply(numeric_markers, function(m) {
        if (!m %in% names(d)) return(NULL)
        ok <- is.finite(d[[target_col]]) & is.finite(d[[m]])
        dd <- d[ok, , drop = FALSE]
        if (nrow(dd) < 8L || stats::sd(dd[[target_col]]) <= 0 || stats::sd(dd[[m]]) <= 0) return(data.frame(biomarker = m, n = nrow(dd), rho = NA_real_, p_value = NA_real_, status = "INSUFFICIENT_DATA", stringsAsFactors = FALSE))
        ct <- tryCatch(stats::cor.test(dd[[target_col]], dd[[m]], method = "spearman", exact = FALSE), error = function(e) NULL)
        data.frame(biomarker = m, n = nrow(dd), rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_, p_value = if (!is.null(ct)) ct$p.value else NA_real_, status = if (!is.null(ct)) "OK" else "TEST_FAILED", stringsAsFactors = FALSE)
    }))
}

.p8_subgroup_target_orr <- function(d, target_col) {
    if (!"strong_ici_context_strict" %in% names(d)) d$strong_ici_context_strict <- NA
    if (!"strong_ici_context_broad" %in% names(d)) d$strong_ici_context_broad <- NA
    if (!"MSS_EBVneg_TMBlow_like" %in% names(d)) d$MSS_EBVneg_TMBlow_like <- NA
    subsets <- list(
        all_tumor = rep(TRUE, nrow(d)),
        strong_ici_context_strict_TRUE = d$strong_ici_context_strict %in% TRUE,
        strong_ici_context_strict_FALSE = d$strong_ici_context_strict %in% FALSE,
        strong_ici_context_broad_TRUE = d$strong_ici_context_broad %in% TRUE,
        strong_ici_context_broad_FALSE = d$strong_ici_context_broad %in% FALSE,
        MSS_EBVneg_TMBlow_like_TRUE = d$MSS_EBVneg_TMBlow_like %in% TRUE
    )
    .p8_bind_rows(lapply(names(subsets), function(snm) {
        idx <- subsets[[snm]]
        dd <- d[idx & is.finite(d[[target_col]]) & d$response_ORR_bin_R %in% c(0L, 1L), , drop = FALSE]
        status <- if (nrow(dd) < 8L) "NOT_EVALUABLE_SMALL_N" else if (length(unique(dd$response_ORR_bin_R)) < 2L) "NOT_EVALUABLE_SINGLE_CLASS" else "EVALUABLE"
        if (status != "EVALUABLE") {
            return(data.frame(subset = snm, n = nrow(dd), n_R = sum(dd$response_ORR_bin_R == 1L, na.rm = TRUE), n_NR = sum(dd$response_ORR_bin_R == 0L, na.rm = TRUE), status = status, median_R = NA_real_, median_NR = NA_real_, wilcox_p = NA_real_, OR_per_1SD = NA_real_, logistic_p = NA_real_, AUC = NA_real_, stringsAsFactors = FALSE))
        }
        wt <- tryCatch(stats::wilcox.test(dd[[target_col]] ~ dd$response_ORR_bin_R, exact = FALSE), error = function(e) NULL)
        lr <- .p8_logistic_endpoint_one(dd, target_col, "response_ORR_bin_R", "ORR", "CR_PR_responder")
        data.frame(subset = snm, n = nrow(dd), n_R = sum(dd$response_ORR_bin_R == 1L), n_NR = sum(dd$response_ORR_bin_R == 0L), status = status, median_R = median(dd[[target_col]][dd$response_ORR_bin_R == 1L], na.rm = TRUE), median_NR = median(dd[[target_col]][dd$response_ORR_bin_R == 0L], na.rm = TRUE), wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_, OR_per_1SD = lr$OR_per_1SD[1], logistic_p = lr$p_value[1], AUC = lr$auc_for_positive_high_value[1], stringsAsFactors = FALSE)
    }))
}

.p8_one_context_marker_models <- function(d, target_col, binary_markers) {
    .p8_bind_rows(lapply(binary_markers, function(m) {
        if (!m %in% names(d)) return(NULL)
        dd <- d[is.finite(d[[target_col]]) & d$response_ORR_bin_R %in% c(0L, 1L) & !is.na(d[[m]]), , drop = FALSE]
        if (nrow(dd) < 15L || length(unique(dd$response_ORR_bin_R)) < 2L || length(unique(dd[[m]])) < 2L) return(NULL)
        dd$x_scaled <- .p8_z(dd[[target_col]])
        dd$z <- as.numeric(as.logical(dd[[m]]))
        fit <- tryCatch(suppressWarnings(stats::glm(response_ORR_bin_R ~ x_scaled + z, data = dd, family = stats::binomial())), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        cf <- tryCatch(suppressWarnings(summary(fit)$coefficients), error = function(e) NULL)
        if (is.null(cf) || !"x_scaled" %in% rownames(cf)) return(NULL)
        probs <- tryCatch(stats::fitted(fit), error = function(e) rep(NA_real_, nrow(dd)))
        unstable <- any(is.finite(probs) & (probs < 1e-6 | probs > 1 - 1e-6), na.rm = TRUE) ||
            any(!is.finite(cf[, "Estimate"]) | !is.finite(cf[, "Std. Error"]) | cf[, "Std. Error"] > 5 | abs(cf[, "Estimate"]) > 8, na.rm = TRUE)
        data.frame(
            model = paste0("ORR ~ ", target_col, "_per_1SD + ", m),
            n = nrow(dd),
            n_R = sum(dd$response_ORR_bin_R == 1L),
            n_NR = sum(dd$response_ORR_bin_R == 0L),
            context_marker = m,
            OLFML2B_OR_per_1SD = if (!unstable) exp(cf["x_scaled", "Estimate"]) else NA_real_,
            OLFML2B_p = if (!unstable) cf["x_scaled", "Pr(>|z|)"] else NA_real_,
            marker_OR = if (!unstable && "z" %in% rownames(cf)) exp(cf["z", "Estimate"]) else NA_real_,
            marker_p = if (!unstable && "z" %in% rownames(cf)) cf["z", "Pr(>|z|)"] else NA_real_,
            model_status = ifelse(unstable, "QUASI_COMPLETE_SEPARATION_OR_UNSTABLE", "OK"),
            warning = ifelse(
                unstable,
                "Exploratory model unstable due to separation/small n; OR and P values suppressed. Use Fisher counts/AUC instead.",
                "Exploratory one-context-marker model; not a definitive multivariable model due to small n."
            ),
            stringsAsFactors = FALSE
        )
    }))
}

.p8_final_molecular_claim <- function(orr_wilcox, orr_logistic, known_binary_orr, known_numeric_orr, subgroup_results, molecular_available) {
    target_w <- orr_wilcox[orr_wilcox$feature == "target_expr", , drop = FALSE]
    target_l <- orr_logistic[orr_logistic$feature == "target_expr", , drop = FALSE]
    target_supported <- (nrow(target_w) && is.finite(target_w$p_value[1]) && target_w$p_value[1] < 0.05) || (nrow(target_l) && is.finite(target_l$p_value[1]) && target_l$p_value[1] < 0.05)
    known_signal <- FALSE
    if (!is.null(known_binary_orr) && nrow(known_binary_orr)) known_signal <- known_signal || any(known_binary_orr$p_value < 0.05, na.rm = TRUE)
    if (!is.null(known_numeric_orr) && nrow(known_numeric_orr)) known_signal <- known_signal || any(known_numeric_orr$p_value < 0.05, na.rm = TRUE)
    target_population_evaluable <- FALSE
    if (!is.null(subgroup_results) && nrow(subgroup_results)) {
        ss <- subgroup_results[subgroup_results$subset == "MSS_EBVneg_TMBlow_like_TRUE", , drop = FALSE]
        target_population_evaluable <- nrow(ss) == 1L && ss$status == "EVALUABLE"
    }
    data.frame(
        claim_domain = c("molecular_annotation", "known_ici_biomarker_context", "OLFML2B_standalone_ORR_prediction", "intended_MSS_EBVneg_TMBlow_population", "recommended_claim_ceiling"),
        status = c(
            ifelse(molecular_available, "AVAILABLE_AND_MERGED", "NOT_AVAILABLE"),
            ifelse(known_signal, "KNOWN_MARKER_SIGNAL_PRESENT_OR_DIRECTIONAL", "KNOWN_MARKER_SIGNAL_WEAK_OR_NOT_DETECTED"),
            ifelse(target_supported, "NOMINAL_EXPLORATORY_ASSOCIATION_NOT_VALIDATED", "NO_NOMINAL_ASSOCIATION"),
            ifelse(target_population_evaluable, "EVALUABLE_EXPLORATORY", "UNDERPOWERED_OR_NOT_EVALUABLE"),
            "Do not claim OLFML2B as a standalone anti-PD-1 biomarker. Use PRJEB25780 only to test transportability and molecular-context boundaries around the established OLFML2B stromal/CAF-ECM association."
        ),
        stringsAsFactors = FALSE
    )
}


.p8_fmt_num <- function(x, digits = 4) {
    if (length(x) == 0L || is.na(x) || !is.finite(as.numeric(x))) return(NA_character_)
    formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

.p8_value_or_na <- function(df, row_filter, col) {
    if (!is.data.frame(df) || !nrow(df) || !col %in% names(df)) return(NA)
    hit <- df[row_filter, , drop = FALSE]
    if (!nrow(hit)) return(NA)
    hit[[col]][1]
}

.p8_publication_summary <- function(analysis_mol, orr_wilcox, orr_logistic, dcr_wilcox, dcr_logistic, known_binary_orr, known_numeric_orr, subgroup_results, claim, target_gene) {
    resp <- table(analysis_mol$response_clean, useNA = "ifany")
    orr <- table(analysis_mol$response_ORR_group, useNA = "ifany")
    dcr <- table(analysis_mol$response_DCR_group, useNA = "ifany")
    tw <- orr_wilcox[orr_wilcox$feature == "target_expr", , drop = FALSE]
    tl <- orr_logistic[orr_logistic$feature == "target_expr", , drop = FALSE]
    dw <- dcr_wilcox[dcr_wilcox$feature == "target_expr", , drop = FALSE]
    dl <- dcr_logistic[dcr_logistic$feature == "target_expr", , drop = FALSE]
    strict <- known_binary_orr[known_binary_orr$marker == "strong_ici_context_strict", , drop = FALSE]
    broad <- known_binary_orr[known_binary_orr$marker == "strong_ici_context_broad", , drop = FALSE]
    mss <- subgroup_results[subgroup_results$subset == "MSS_EBVneg_TMBlow_like_TRUE", , drop = FALSE]
    pdl1 <- known_numeric_orr[known_numeric_orr$feature == "PDL1_CPS", , drop = FALSE]
    inflamed5 <- known_numeric_orr[known_numeric_orr$feature %in% c("Inflamed5_score_value", "Inflamed_5gene_exploratory"), , drop = FALSE]
    tmb <- known_numeric_orr[known_numeric_orr$feature == "TMB_value", , drop = FALSE]
    data.frame(
        item = c(
            "part8_primary_dataset_policy",
            "tumor_only_evaluable_n",
            "response_counts_CR_PR_SD_PD",
            "ORR_counts_R_NR",
            "DCR_counts_DCR_PD",
            paste0(target_gene, "_ORR_wilcoxon_p"),
            paste0(target_gene, "_ORR_AUC"),
            paste0(target_gene, "_DCR_wilcoxon_p"),
            paste0(target_gene, "_DCR_AUC"),
            "strong_ici_context_strict_ORR_fisher_p",
            "strong_ici_context_strict_counts_TRUE_R_NR",
            "TMB_value_ORR_AUC",
            "PDL1_CPS_ORR_AUC",
            "Inflamed5_score_ORR_AUC",
            "MSS_EBVneg_TMBlow_like_subgroup",
            "final_claim_ceiling"
        ),
        value = c(
            "PRJEB25780/TIGER tumor-only cohort is the only Part8 primary ICI sequencing cohort; no additional external cohort is pursued in this part.",
            as.character(nrow(analysis_mol)),
            paste0("CR=", resp[["CR"]] %||% 0, "; PR=", resp[["PR"]] %||% 0, "; SD=", resp[["SD"]] %||% 0, "; PD=", resp[["PD"]] %||% 0),
            paste0("R=", orr[["R"]] %||% 0, "; NR=", orr[["NR"]] %||% 0),
            paste0("DCR=", dcr[["DCR"]] %||% 0, "; PD=", dcr[["PD"]] %||% 0),
            .p8_fmt_num(if (nrow(tw)) tw$p_value[1] else NA),
            .p8_fmt_num(if (nrow(tl)) tl$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(dw)) dw$p_value[1] else NA),
            .p8_fmt_num(if (nrow(dl)) dl$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(strict)) strict$p_value[1] else NA),
            if (nrow(strict)) paste0("TRUE_R=", strict$marker_TRUE_positive[1], "; TRUE_NR=", strict$marker_TRUE_negative[1]) else NA_character_,
            .p8_fmt_num(if (nrow(tmb)) tmb$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(pdl1)) pdl1$auc_for_positive_high_value[1] else NA),
            if (nrow(mss)) paste0("n=", mss$n[1], "; R=", mss$n_R[1], "; NR=", mss$n_NR[1], "; status=", mss$status[1]) else NA_character_,
            claim$status[claim$claim_domain == "recommended_claim_ceiling"][1]
        ),
        interpretation = c(
            "Part8 is now finalized as a single-primary-dataset molecular-context boundary analysis.",
            "Tumor-only samples are defined internally by TIGER patient_name after excluding -N samples.",
            "RECIST response composition in the analysis subset.",
            "Primary endpoint is ORR: CR/PR versus SD/PD.",
            "DCR is exploratory: CR/PR/SD versus PD.",
            paste0(target_gene, " single-gene ORR signal; should not be overclaimed if not significant."),
            "Rank-based AUC for high expression predicting ORR.",
            paste0(target_gene, " exploratory DCR signal."),
            "Rank-based AUC for high expression predicting DCR.",
            "Positive-control molecular-context signal; supports MSI/EBV/TMB/PD-L1-driven response biology.",
            "Counts in strong ICI context; exact counts are preferred over unstable logistic ORs.",
            "Continuous mutation-burden signal for ORR.",
            "Continuous PD-L1 CPS signal for ORR.",
            "Target MSS/EBV-negative/TMB-low-like population evaluability.",
            "Publication claim ceiling for Part8."
        ),
        stringsAsFactors = FALSE
    )
}

.p8_write_methods_results_text <- function(path, summary_tab, target_gene = "OLFML2B") {
    get <- function(item) summary_tab$value[summary_tab$item == item][1]
    txt <- c(
        "Part8 finalized design and interpretation",
        "",
        "Design: Part8 was finalized as a single-primary-dataset analysis based on the PRJEB25780/TIGER anti-PD-1 cohort. No additional external immunotherapy cohort is introduced into Part8. TIGER expression and response annotations are used as the expression-response backbone, and Kim/Bagaev molecular annotations are used only to interpret MSI/EBV/TMB/PD-L1 molecular context.",
        "",
        paste0("Analysis set: ", get("tumor_only_evaluable_n"), " tumor-only response-evaluable samples were retained after excluding TIGER samples with patient_name ending in -N. Response composition: ", get("response_counts_CR_PR_SD_PD"), ". Primary ORR endpoint: ", get("ORR_counts_R_NR"), ". Exploratory DCR endpoint: ", get("DCR_counts_DCR_PD"), "."),
        "",
        paste0("Target-gene boundary: ", target_gene, " did not show a robust standalone ORR/DCR predictive signal in this primary cohort. ORR Wilcoxon P=", get(paste0(target_gene, "_ORR_wilcoxon_p")), ", ORR AUC=", get(paste0(target_gene, "_ORR_AUC")), "; DCR Wilcoxon P=", get(paste0(target_gene, "_DCR_wilcoxon_p")), ", DCR AUC=", get(paste0(target_gene, "_DCR_AUC")), "."),
        "",
        paste0("Molecular context: established ICI-sensitive molecular contexts were strongly enriched for ORR. The strict strong-ICI context Fisher P=", get("strong_ici_context_strict_ORR_fisher_p"), " with ", get("strong_ici_context_strict_counts_TRUE_R_NR"), ". Mutation burden, PD-L1 CPS, and prespecified exploratory five-gene inflamed-expression score (an equal-weight mean-z score, not a validated clinical classifier) were summarized as context controls (TMB AUC=", get("TMB_value_ORR_AUC"), "; PD-L1 CPS AUC=", get("PDL1_CPS_ORR_AUC"), "; five-gene inflamed-score AUC=", get("Inflamed5_score_ORR_AUC"), ")."),
        "",
        paste0("Target population boundary: the MSS/EBV-negative/TMB-low-like subgroup was ", get("MSS_EBVneg_TMBlow_like_subgroup"), ", indicating that the intended non-strong-ICI-sensitive stromal-remodeling population is not evaluable in PRJEB25780."),
        "",
        paste0("Claim ceiling: ", get("final_claim_ceiling"))
    )
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(txt, con = path, useBytes = TRUE)
    invisible(path)
}

.p8_output_manifest <- function(out_tables, out_figures) {
    files <- c(list.files(out_tables, full.names = TRUE, recursive = FALSE), list.files(out_figures, full.names = TRUE, recursive = FALSE))
    data.frame(
        file = files,
        basename = basename(files),
        type = ifelse(grepl("\\.csv$", files, ignore.case = TRUE), "table_csv", ifelse(grepl("\\.png$|\\.pdf$", files, ignore.case = TRUE), "figure", "other")),
        size_bytes = file.info(files)$size,
        stringsAsFactors = FALSE
    )
}

run_part8_ici_prjeb25780_tiger_only <- function(
    root = "D:/OLFML2B_STAD",
    target_gene = "OLFML2B",
    tiger_dir = file.path(root, "data/cache/ICI_PRJEB25780/TIGER"),
    molecular_file = NA_character_,
    expression_transform = c("auto", "none", "log2"),
    make_figures = TRUE
) {
    expression_transform <- match.arg(expression_transform)

    out_tables <- file.path(root, "output/tables/Part8_Immunotherapy")
    out_figures <- file.path(root, "output/figures/Part8_Immunotherapy")
    out_objects <- file.path(root, "output/objects")
    out_reports <- file.path(root, "output/reports/Part8_Immunotherapy")
    derived_dir <- file.path(root, "data/derived/ICI_PRJEB25780")
    dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_objects, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_reports, recursive = TRUE, showWarnings = FALSE)
    dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(out_reports, "Part8_Immunotherapy_TIGER_only_molecular_context.log")
    if (file.exists(log_file)) unlink(log_file)
    .p8_msg("Starting Part8 | ", PART8_ICI_TIGER_VERSION, log_file = log_file)
    .p8_msg("Root: ", root, log_file = log_file)
    .p8_msg("TIGER directory: ", tiger_dir, log_file = log_file)
    .p8_msg("Policy: TIGER-only expression/response; no ENA manifest is read or used. Kim/Bagaev molecular annotation is optional context only.", log_file = log_file)

    files <- .p8_find_tiger_files(tiger_dir)
    .p8_write_csv(files$file_audit, file.path(out_tables, "00_TIGER_file_discovery_audit.csv"))
    .p8_msg("Expression file: ", files$expression_file, log_file = log_file)
    .p8_msg("Clinical file: ", files$clinical_file, log_file = log_file)

    expr_raw <- .p8_read_any(files$expression_file)
    clin_raw <- .p8_read_any(files$clinical_file)
    clin <- as.data.frame(clin_raw, stringsAsFactors = FALSE, check.names = FALSE)

    if (!all(c("sample_id", "patient_name", "response") %in% names(clin))) {
        .p8_stop("TIGER clinical table must contain sample_id, patient_name, and response columns. Found: ", paste(names(clin), collapse = ", "))
    }

    expr_prep <- .p8_prepare_expression(expr_raw, transform = expression_transform, log_file = log_file)
    target <- .p8_target_expression(expr_prep, target_gene)
    .p8_write_csv(target$mapping, file.path(out_tables, "01_target_gene_mapping_audit.csv"))
    if (!isTRUE(target$mapping$detected_after_symbol_collapse[1])) .p8_stop("Target gene not detected in TIGER expression matrix: ", target_gene)

    clin <- .p8_classify_tiger_sample_type(clin)
    clin <- .p8_clean_response(clin)
    clin <- .p8_add_response_endpoints(clin)
    clin$sample_id <- .p8_clean(clin$sample_id)

    target_expr <- target$values
    names(target_expr)[names(target_expr) == "target_expr"] <- paste0(target_gene, "_TIGER_expr")
    merged <- merge(clin, target_expr[, c("sample_id", paste0(target_gene, "_TIGER_expr")), drop = FALSE], by = "sample_id", all.x = TRUE, all.y = FALSE)
    names(merged)[names(merged) == paste0(target_gene, "_TIGER_expr")] <- "target_expr"

    gene_sets <- .p8_module_gene_sets()
    module <- .p8_score_modules(expr_prep$matrix, gene_sets)
    .p8_write_csv(module$coverage, file.path(out_tables, "02_module_gene_coverage.csv"))
    merged <- merge(merged, module$scores, by = "sample_id", all.x = TRUE, all.y = FALSE)

    completeness <- .p8_clinical_completeness(clin)
    .p8_write_csv(completeness, file.path(out_tables, "03_TIGER_clinical_column_completeness.csv"))

    sample_audit <- .p8_make_sample_audit(expr_prep, clin, merged)
    .p8_write_csv(sample_audit, file.path(out_tables, "04_TIGER_internal_sample_type_and_matching_audit.csv"))

    normal_like <- merged[merged$sample_type_tiger == "normal_like", , drop = FALSE]
    tumor_all <- merged[merged$sample_type_tiger == "tumor", , drop = FALSE]
    analysis_df <- tumor_all[
        is.finite(tumor_all$target_expr) & tumor_all$response_ORR_group %in% c("R", "NR"),
        , drop = FALSE
    ]
    analysis_df <- .p8_add_response_endpoints(analysis_df)

    .p8_write_csv(merged, file.path(derived_dir, "Part8_TIGER_all_78_with_sample_type_flags.csv"))
    .p8_write_csv(normal_like, file.path(out_tables, "05_TIGER_normal_like_samples_excluded_by_patient_name_suffix.csv"))
    .p8_write_csv(analysis_df, file.path(out_tables, "06_TIGER_tumor_only_response_analysis_input.csv"))

    .p8_msg("TIGER matched samples: ", nrow(merged), log_file = log_file)
    .p8_msg("TIGER normal-like excluded by -N suffix: ", nrow(normal_like), log_file = log_file)
    .p8_msg("TIGER tumor-only response-evaluable samples: ", nrow(analysis_df), log_file = log_file)
    .p8_msg("Response table: ", paste(capture.output(print(table(analysis_df$response_clean, useNA = "ifany"))), collapse = " | "), log_file = log_file)
    .p8_msg("ORR table: ", paste(capture.output(print(table(analysis_df$response_ORR_group, useNA = "ifany"))), collapse = " | "), log_file = log_file)
    .p8_msg("DCR table: ", paste(capture.output(print(table(analysis_df$response_DCR_group, useNA = "ifany"))), collapse = " | "), log_file = log_file)

    feature_cols <- c("target_expr", names(module$scores)[names(module$scores) != "sample_id"])

    orr_results <- .p8_analyze_features_endpoint(analysis_df, feature_cols, endpoint = "ORR")
    dcr_results <- .p8_analyze_features_endpoint(analysis_df, feature_cols, endpoint = "DCR")

    .p8_write_csv(orr_results$summary, file.path(out_tables, "07_ORR_response_group_summary_all_features.csv"))
    .p8_write_csv(orr_results$wilcox, file.path(out_tables, "08_ORR_response_wilcoxon_all_features.csv"))
    .p8_write_csv(orr_results$logistic, file.path(out_tables, "09_ORR_response_logistic_per_1SD_all_features.csv"))
    .p8_write_csv(orr_results$fisher, file.path(out_tables, "10_ORR_response_median_split_fisher_all_features.csv"))

    .p8_write_csv(dcr_results$summary, file.path(out_tables, "11_DCR_response_group_summary_all_features.csv"))
    .p8_write_csv(dcr_results$wilcox, file.path(out_tables, "12_DCR_response_wilcoxon_all_features.csv"))
    .p8_write_csv(dcr_results$logistic, file.path(out_tables, "13_DCR_response_logistic_per_1SD_all_features.csv"))
    .p8_write_csv(dcr_results$fisher, file.path(out_tables, "14_DCR_response_median_split_fisher_all_features.csv"))

    # Backward-compatible legacy ORR names for downstream scripts that still expect the old table names.
    .p8_write_csv(orr_results$summary, file.path(out_tables, "07_response_group_summary_all_features.csv"))
    .p8_write_csv(orr_results$wilcox, file.path(out_tables, "08_response_wilcoxon_all_features.csv"))
    .p8_write_csv(orr_results$logistic, file.path(out_tables, "09_response_logistic_per_1SD_all_features.csv"))
    .p8_write_csv(orr_results$fisher, file.path(out_tables, "10_response_median_split_fisher_all_features.csv"))

    if (is.na(molecular_file) || !nzchar(molecular_file)) molecular_file <- .p8_find_molecular_file(root)
    mol_ctx <- .p8_clean_molecular_context(analysis_df, molecular_file)
    analysis_mol <- mol_ctx$data
    .p8_write_csv(mol_ctx$merge_audit, file.path(out_tables, "15_molecular_annotation_merge_audit.csv"))
    .p8_write_csv(mol_ctx$completeness, file.path(out_tables, "16_known_ici_biomarker_completeness.csv"))
    .p8_write_csv(analysis_mol, file.path(out_tables, "14_TIGER_tumor_only_with_Bagaev_Kim_molecular_context.csv"))
    .p8_msg("Molecular annotation available: ", mol_ctx$available, log_file = log_file)
    if (isTRUE(mol_ctx$available)) {
        .p8_msg("Molecular matched samples: ", sum(.p8_standard_patient_id(analysis_df$patient_id_tiger %||% analysis_df$patient_name) %in% .p8_standard_patient_id(analysis_mol$patient_id)), " / ", nrow(analysis_df), log_file = log_file)
    }

    binary_markers <- c("MSI_H", "EBV_pos", "TMB_high_10", "TMB_high_median", "PDL1_CPS_ge1", "PDL1_CPS_ge5", "PDL1_CPS_ge10", "Inflamed5_high_median", "Inflamed5_MSS_EBVneg_TMBlow_like", "strong_ici_context_strict", "strong_ici_context_broad", "MSS_EBVneg_TMBlow_like")
    binary_markers <- binary_markers[binary_markers %in% names(analysis_mol)]
    numeric_markers <- c("TMB_value", "PDL1_CPS", "Inflamed5_score_value", "Inflamed_5gene_exploratory")
    numeric_markers <- numeric_markers[numeric_markers %in% names(analysis_mol)]

    known_binary_orr <- .p8_bind_rows(lapply(binary_markers, function(v) .p8_fisher_binary_marker(analysis_mol, v, "response_ORR_bin_R", "ORR")))
    known_binary_dcr <- .p8_bind_rows(lapply(binary_markers, function(v) .p8_fisher_binary_marker(analysis_mol, v, "response_DCR_bin", "DCR")))
    if (nrow(known_binary_orr)) known_binary_orr$p_fdr <- stats::p.adjust(known_binary_orr$p_value, method = "BH")
    if (nrow(known_binary_dcr)) known_binary_dcr$p_fdr <- stats::p.adjust(known_binary_dcr$p_value, method = "BH")
    .p8_write_csv(known_binary_orr, file.path(out_tables, "17_known_ici_binary_biomarkers_ORR_fisher.csv"))
    .p8_write_csv(known_binary_dcr, file.path(out_tables, "18_known_ici_binary_biomarkers_DCR_fisher.csv"))

    known_numeric_orr <- .p8_bind_rows(lapply(numeric_markers, function(v) .p8_numeric_marker_logistic(analysis_mol, v, "response_ORR_bin_R", "ORR")))
    known_numeric_dcr <- .p8_bind_rows(lapply(numeric_markers, function(v) .p8_numeric_marker_logistic(analysis_mol, v, "response_DCR_bin", "DCR")))
    if (nrow(known_numeric_orr)) known_numeric_orr$p_fdr <- stats::p.adjust(known_numeric_orr$p_value, method = "BH")
    if (nrow(known_numeric_dcr)) known_numeric_dcr$p_fdr <- stats::p.adjust(known_numeric_dcr$p_value, method = "BH")
    .p8_write_csv(known_numeric_orr, file.path(out_tables, "19_known_ici_numeric_biomarkers_ORR_logistic_auc.csv"))
    .p8_write_csv(known_numeric_dcr, file.path(out_tables, "20_known_ici_numeric_biomarkers_DCR_logistic_auc.csv"))

    target_by_context <- .p8_target_by_binary_context(analysis_mol, "target_expr", binary_markers)
    if (nrow(target_by_context)) target_by_context$p_fdr <- stats::p.adjust(target_by_context$p_value, method = "BH")
    .p8_write_csv(target_by_context, file.path(out_tables, "21_OLFML2B_by_known_ici_binary_context.csv"))

    target_cor_numeric <- .p8_target_numeric_correlations(analysis_mol, "target_expr", numeric_markers)
    if (nrow(target_cor_numeric)) target_cor_numeric$p_fdr <- stats::p.adjust(target_cor_numeric$p_value, method = "BH")
    .p8_write_csv(target_cor_numeric, file.path(out_tables, "22_OLFML2B_known_ici_numeric_correlations.csv"))

    module_cols <- setdiff(feature_cols, "target_expr")
    cor_tab <- .p8_correlate_target_modules(analysis_df, "target_expr", module_cols)
    .p8_write_csv(cor_tab, file.path(out_tables, "23_OLFML2B_module_spearman_correlations.csv"))

    combo <- .p8_combination_test(analysis_df, "target_expr", "CAF_TGFb_ECM_axis")
    .p8_write_csv(combo, file.path(out_tables, "24_OLFML2B_high_plus_CAF_TGFb_ECM_high_ORR_fisher.csv"))

    subgroup_results <- .p8_subgroup_target_orr(analysis_mol, "target_expr")
    .p8_write_csv(subgroup_results, file.path(out_tables, "25_OLFML2B_ORR_in_molecular_context_subsets.csv"))

    context_models <- .p8_one_context_marker_models(analysis_mol, "target_expr", binary_markers)
    .p8_write_csv(context_models, file.path(out_tables, "26_OLFML2B_ORR_one_biomarker_context_models.csv"))

    claim <- .p8_final_molecular_claim(orr_results$wilcox, orr_results$logistic, known_binary_orr, known_numeric_orr, subgroup_results, mol_ctx$available)
    .p8_write_csv(claim, file.path(out_tables, "27_final_molecular_context_claim_ceiling.csv"))
    # Also keep the old expected name as a pointer to the final claim ceiling.
    .p8_write_csv(claim, file.path(out_tables, "13_claim_ceiling_go_nogo_summary.csv"))

    publication_summary <- .p8_publication_summary(
        analysis_mol = analysis_mol,
        orr_wilcox = orr_results$wilcox,
        orr_logistic = orr_results$logistic,
        dcr_wilcox = dcr_results$wilcox,
        dcr_logistic = dcr_results$logistic,
        known_binary_orr = known_binary_orr,
        known_numeric_orr = known_numeric_orr,
        subgroup_results = subgroup_results,
        claim = claim,
        target_gene = target_gene
    )
    .p8_write_csv(publication_summary, file.path(out_tables, "28_part8_primary_dataset_publication_summary.csv"))
    .p8_write_methods_results_text(file.path(out_reports, "29_part8_primary_dataset_methods_results_text.txt"), publication_summary, target_gene = target_gene)

    if (isTRUE(make_figures)) {
        .p8_make_figures(analysis_df, "target_expr", orr_results$wilcox, cor_tab, out_figures, target_gene)
    }
    .p8_write_csv(.p8_output_manifest(out_tables, out_figures), file.path(out_tables, "30_part8_output_manifest.csv"))

    out <- list(
        version = PART8_ICI_TIGER_VERSION,
        generated_at = .p8_timestamp(),
        root = root,
        tiger_dir = tiger_dir,
        expression_file = files$expression_file,
        clinical_file = files$clinical_file,
        molecular_file = molecular_file,
        sample_audit = sample_audit,
        target_mapping = target$mapping,
        module_coverage = module$coverage,
        analysis_input = analysis_df,
        molecular_context_input = analysis_mol,
        ORR_summary = orr_results$summary,
        ORR_wilcox = orr_results$wilcox,
        ORR_logistic = orr_results$logistic,
        ORR_fisher = orr_results$fisher,
        DCR_summary = dcr_results$summary,
        DCR_wilcox = dcr_results$wilcox,
        DCR_logistic = dcr_results$logistic,
        DCR_fisher = dcr_results$fisher,
        known_binary_ORR = known_binary_orr,
        known_numeric_ORR = known_numeric_orr,
        target_context_subgroups = subgroup_results,
        target_module_correlations = cor_tab,
        combination_test = combo,
        claim_ceiling = claim,
        publication_summary = publication_summary,
        output_tables = out_tables,
        output_figures = out_figures
    )

    .p8_save_rds(out, file.path(out_objects, "Part8_Immunotherapy_TIGER_only_molecular_context.rds"))
    .p8_save_rds(out, file.path(out_objects, "Part8_Immunotherapy_TIGER_only.rds"))
    .p8_msg("Part8 complete. Tables: ", out_tables, log_file = log_file)
    .p8_msg("Part8 complete. Figures: ", out_figures, log_file = log_file)
    invisible(out)
}

# Backward-compatible short alias for the new Part8.
run_part8 <- run_part8_ici_prjeb25780_tiger_only


# ==============================================================================
# OLFML2B-specific final methodology overrides
# ==============================================================================

OLFML2B_PART8_VERSION <- PART8_ICI_TIGER_VERSION

.p8_candidate_data_roots <- function(root, shared_data_root = Sys.getenv("OLFML2B_STAD_SHARED_DATA_ROOT", unset = "")) {
    roots <- c(file.path(root, "data"))
    if (nzchar(shared_data_root)) roots <- c(roots, shared_data_root)
    # Public-release portability rule: no discovery from unrelated project trees.
    # Optional shared storage remains available only through OLFML2B_STAD_SHARED_DATA_ROOT.
    roots <- unique(normalizePath(roots, winslash = "/", mustWork = FALSE))
    roots
}

.p8_resolve_tiger_dir <- function(root, tiger_dir = NULL, shared_data_root = Sys.getenv("OLFML2B_STAD_SHARED_DATA_ROOT", unset = "")) {
    candidates <- character()
    if (!is.null(tiger_dir) && length(tiger_dir) && !is.na(tiger_dir[1]) && nzchar(tiger_dir[1])) {
        candidates <- c(candidates, tiger_dir[1])
    }
    roots <- .p8_candidate_data_roots(root, shared_data_root)
    candidates <- c(candidates, file.path(roots, "cache", "ICI_PRJEB25780", "TIGER"))
    candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
    hit <- candidates[dir.exists(candidates)]
    if (!length(hit)) {
        .p8_stop(
            "PRJEB25780/TIGER directory was not found. Checked: ",
            paste(candidates, collapse = "; "),
            ". Place the TIGER files in the OLFML2B project data tree or set OLFML2B_STAD_SHARED_DATA_ROOT."
        )
    }
    hit[1]
}

# Search project-local first, then the configured immutable/shared STAD data tree.
.p8_find_molecular_file <- function(root) {
    roots <- .p8_candidate_data_roots(root)
    candidates <- unique(c(
        unlist(lapply(roots, function(dr) file.path(dr, "cache", "ICI_PRJEB25780", c(
            "PRJEB25780_Kim2018_Bagaev2021_clinical_curated.csv",
            "Bagaev2021_Kim_PRJEB25780_clinical_rows_raw.csv",
            "Bagaev2021_mmc5_hits_Gastric.csv"
        )))),
        file.path(root, "config", "PRJEB25780_Kim2018_Bagaev2021_clinical_curated.csv")
    ))
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) hit[1] else NA_character_
}

# Freeze Part8 signatures to the same names and genes used by OLFML2B Part4.
# The five-gene inflamed score is an exploratory context control only.
.p8_module_gene_sets <- function() {
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
        Proliferation = c("MKI67","TOP2A","PCNA","MCM2","MCM5","MCM6","UBE2C","CCNB1","CCNB2","CDK1","AURKA","BUB1"),
        Inflamed_5gene_exploratory = c("CD274","CXCL9","CXCL10","GZMB","IFNG")
    )
}

.p8_firth_fit <- function(formula, data) {
    if (!requireNamespace("logistf", quietly = TRUE)) {
        return(structure(list(message = "Package 'logistf' is required for prespecified Firth logistic regression."), class = "p8_missing_logistf"))
    }
    tryCatch(
        suppressWarnings(logistf::logistf(formula, data = data, pl = TRUE, firth = TRUE)),
        error = function(e) e
    )
}

# Continuous one-feature model: Firth penalization is mandatory because the
# cohort contains only 45 tumors and 12 ORR responders.
.p8_logistic_endpoint_one <- function(d, feature, y_col, endpoint_name, positive_label = "positive") {
    ok <- is.finite(d[[feature]]) & d[[y_col]] %in% c(0L, 1L)
    dd <- d[ok, , drop = FALSE]
    auc <- if (nrow(dd)) .p8_auc_rank(dd[[feature]], dd[[y_col]]) else NA_real_
    base <- function(status, estimate = NA_real_, se = NA_real_, z = NA_real_, p_value = NA_real_,
                     OR = NA_real_, lcl = NA_real_, ucl = NA_real_, note = NA_character_) {
        data.frame(
            endpoint = endpoint_name, feature = feature, n = nrow(dd),
            n_positive = sum(dd[[y_col]] == 1L, na.rm = TRUE),
            n_negative = sum(dd[[y_col]] == 0L, na.rm = TRUE),
            estimate = estimate, se = se, z = z, p_value = p_value,
            OR_per_1SD = OR, OR_lcl = lcl, OR_ucl = ucl,
            auc_for_positive_high_value = auc, positive_class = positive_label,
            model_method = "Firth_penalized_logistic_profile_likelihood",
            status = status, inference_note = note, stringsAsFactors = FALSE
        )
    }
    class_counts <- table(factor(dd[[y_col]], levels = c(0L, 1L)))
    if (nrow(dd) < 8L || min(class_counts) < 3L || stats::sd(dd[[feature]], na.rm = TRUE) <= 0) {
        return(base("INSUFFICIENT_DATA", note = "Requires >=8 evaluable tumors and >=3 observations in each response class."))
    }
    dd$feature_z <- .p8_z(dd[[feature]])
    fit <- .p8_firth_fit(stats::as.formula(paste0(y_col, " ~ feature_z")), dd)
    if (inherits(fit, "p8_missing_logistf")) return(base("MISSING_LOGISTF_DEPENDENCY", note = fit$message))
    if (inherits(fit, "error")) return(base("FIRTH_MODEL_FAILED", note = conditionMessage(fit)))
    if (!"feature_z" %in% names(fit$coefficients)) return(base("COEFFICIENT_MISSING", note = "feature_z was not estimable."))
    est <- as.numeric(fit$coefficients["feature_z"])
    se <- tryCatch(sqrt(as.numeric(fit$var["feature_z", "feature_z"])), error = function(e) NA_real_)
    pval <- tryCatch(as.numeric(fit$prob["feature_z"]), error = function(e) NA_real_)
    low <- tryCatch(as.numeric(fit$ci.lower["feature_z"]), error = function(e) NA_real_)
    high <- tryCatch(as.numeric(fit$ci.upper["feature_z"]), error = function(e) NA_real_)
    stable <- is.finite(est) && is.finite(se) && se > 0 && abs(est) <= 8 && se <= 5
    if (!stable) return(base("FIRTH_UNSTABLE", estimate = est, se = se, note = "Penalized estimate remained unstable; OR, CI and P are suppressed."))
    base(
        "OK_FIRTH", estimate = est, se = se, z = est / se, p_value = pval,
        OR = exp(est), lcl = exp(low), ucl = exp(high),
        note = "Exploratory single-cohort association; not a treatment-selection model."
    )
}

# One-context-marker sensitivity models are also fitted with Firth penalization.
.p8_one_context_marker_models <- function(d, target_col, binary_markers) {
    .p8_bind_rows(lapply(binary_markers, function(m) {
        if (!m %in% names(d)) return(NULL)
        dd <- d[is.finite(d[[target_col]]) & d$response_ORR_bin_R %in% c(0L, 1L) & !is.na(d[[m]]), , drop = FALSE]
        if (nrow(dd) < 15L || min(table(factor(dd$response_ORR_bin_R, levels = c(0L, 1L)))) < 3L || length(unique(dd[[m]])) < 2L) {
            return(data.frame(model = paste0("ORR ~ ", target_col, "_per_1SD + ", m), n = nrow(dd), n_R = sum(dd$response_ORR_bin_R == 1L), n_NR = sum(dd$response_ORR_bin_R == 0L), context_marker = m, OLFML2B_OR_per_1SD = NA_real_, OLFML2B_p = NA_real_, marker_OR = NA_real_, marker_p = NA_real_, model_method = "Firth_penalized_logistic_profile_likelihood", model_status = "INSUFFICIENT_DATA", warning = "Exploratory sensitivity model not evaluable.", stringsAsFactors = FALSE))
        }
        dd$x_scaled <- .p8_z(dd[[target_col]])
        dd$z <- as.numeric(as.logical(dd[[m]]))
        fit <- .p8_firth_fit(response_ORR_bin_R ~ x_scaled + z, dd)
        if (inherits(fit, "p8_missing_logistf")) {
            return(data.frame(model = paste0("ORR ~ ", target_col, "_per_1SD + ", m), n = nrow(dd), n_R = sum(dd$response_ORR_bin_R == 1L), n_NR = sum(dd$response_ORR_bin_R == 0L), context_marker = m, OLFML2B_OR_per_1SD = NA_real_, OLFML2B_p = NA_real_, marker_OR = NA_real_, marker_p = NA_real_, model_method = "Firth_penalized_logistic_profile_likelihood", model_status = "MISSING_LOGISTF_DEPENDENCY", warning = fit$message, stringsAsFactors = FALSE))
        }
        if (inherits(fit, "error") || !all(c("x_scaled", "z") %in% names(fit$coefficients))) return(NULL)
        data.frame(
            model = paste0("ORR ~ ", target_col, "_per_1SD + ", m),
            n = nrow(dd), n_R = sum(dd$response_ORR_bin_R == 1L), n_NR = sum(dd$response_ORR_bin_R == 0L),
            context_marker = m,
            OLFML2B_OR_per_1SD = exp(as.numeric(fit$coefficients["x_scaled"])),
            OLFML2B_p = as.numeric(fit$prob["x_scaled"]),
            marker_OR = exp(as.numeric(fit$coefficients["z"])), marker_p = as.numeric(fit$prob["z"]),
            model_method = "Firth_penalized_logistic_profile_likelihood",
            model_status = "OK_FIRTH",
            warning = "Exploratory one-context-marker sensitivity model; not a validated multivariable prediction model.",
            stringsAsFactors = FALSE
        )
    }))
}

.p8_patient_unit_audit <- function(analysis_df) {
    pid <- .p8_standard_patient_id(analysis_df$patient_id_tiger %||% analysis_df$patient_name)
    sid <- .p8_clean(analysis_df$sample_id)
    data.frame(
        criterion = c("analysis_rows", "unique_sample_ids", "unique_patient_ids", "duplicated_sample_ids", "duplicated_patient_ids", "inferential_unit"),
        value = c(nrow(analysis_df), length(unique(sid)), length(unique(pid)), sum(duplicated(sid)), sum(duplicated(pid)), "one tumor sample per patient"),
        status = c("INFO", ifelse(anyDuplicated(sid) == 0L, "PASS", "FAIL"), ifelse(anyDuplicated(pid) == 0L, "PASS", "FAIL"), ifelse(anyDuplicated(sid) == 0L, "PASS", "FAIL"), ifelse(anyDuplicated(pid) == 0L, "PASS", "FAIL"), ifelse(anyDuplicated(sid) == 0L && anyDuplicated(pid) == 0L, "PASS", "FAIL")),
        stringsAsFactors = FALSE
    )
}

.p8_structural_qc <- function(root, out, target_gene = "OLFML2B") {
    tab_dir <- file.path(root, "output", "tables", "Part8_Immunotherapy")
    read_if <- function(nm) {
        fp <- file.path(tab_dir, nm)
        if (file.exists(fp)) utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
    }
    target <- read_if("01_target_gene_mapping_audit.csv")
    unit <- read_if("06a_patient_level_inferential_unit_audit.csv")
    orr_models <- read_if("09_ORR_response_logistic_per_1SD_all_features.csv")
    dcr_models <- read_if("13_DCR_response_logistic_per_1SD_all_features.csv")
    mol <- read_if("15_molecular_annotation_merge_audit.csv")
    pub <- read_if("28_part8_primary_dataset_publication_summary.csv")
    exact_target <- nrow(target) == 1L && isTRUE(target$detected_after_symbol_collapse[1]) && target$exact_rows_before_collapse[1] >= 1L
    unit_ok <- nrow(unit) && all(unit$status[unit$criterion %in% c("unique_sample_ids", "unique_patient_ids", "inferential_unit")] == "PASS")
    model_rows <- rbind(orr_models, dcr_models)
    firth_ok <- nrow(model_rows) && "model_method" %in% names(model_rows) && all(model_rows$model_method == "Firth_penalized_logistic_profile_likelihood") && !any(model_rows$status == "MISSING_LOGISTF_DEPENDENCY")
    fdr_ok <- nrow(orr_models) && nrow(dcr_models) && "p_fdr" %in% names(orr_models) && "p_fdr" %in% names(dcr_models)
    data.frame(
        criterion = c("TIGER_only_input_policy", "exact_OLFML2B_symbol_mapping", "tumor_only_patient_level_unit", "Firth_logistic_contract", "BH_within_family_adjustment", "molecular_context_optional_not_primary", "exploratory_inflamed5_not_validated_classifier", "result_independent_completion", "claim_ceiling", "publication_summary_available"),
        status = c("PASS", ifelse(exact_target, "PASS", "FAIL"), ifelse(unit_ok, "PASS", "FAIL"), ifelse(firth_ok, "PASS", "FAIL"), ifelse(fdr_ok, "PASS", "FAIL"), ifelse(nrow(mol), "PASS", "NOT_AVAILABLE_OPTIONAL"), "PASS", "PASS", "EXPLORATORY_CONTEXT_ONLY", ifelse(nrow(pub), "PASS", "FAIL")),
        interpretation = c(
            "Only TIGER expression/response data are used; ENA run counts are not analytic samples.",
            "Exact OLFML2B symbol mapping is required; no similarly named gene or pseudogene substitution is allowed.",
            "Normal-like -N samples are excluded and the inferential unit is one tumor sample per patient.",
            "Small-sample response models use Firth profile-likelihood logistic regression.",
            "Multiplicity is controlled separately within prespecified ORR and DCR feature families.",
            "Kim/Bagaev MSI, EBV, TMB and PD-L1 annotations are contextual covariates/positive controls, not a second validation cohort.",
            "The equal-weight five-gene inflamed score is exploratory and must not be named or interpreted as a validated classifier.",
            "Pipeline completion is determined by data and method contracts, never by significance or direction.",
            "No standalone anti-PD-1 biomarker, treatment-selection, clinical-utility, mediation or causal claim.",
            "A compact publication summary and methods/results boundary text are generated."
        ), stringsAsFactors = FALSE
    )
}


# Apply the same minimum-coverage rule used by current OLFML2B Part4.
.p8_score_modules <- function(expr_mat, gene_sets) {
    sample_ids <- colnames(expr_mat)
    scores <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
    coverage_rows <- list()
    for (nm in names(gene_sets)) {
        requested <- unique(toupper(gene_sets[[nm]]))
        present <- intersect(requested, rownames(expr_mat))
        missing <- setdiff(requested, present)
        frac <- length(present) / max(length(requested), 1L)
        eligible <- if (identical(nm, "Inflamed_5gene_exploratory")) {
            length(present) == length(requested)
        } else {
            length(present) >= 5L && frac >= 0.50
        }
        coverage_rows[[nm]] <- data.frame(
            module = nm, n_requested_genes = length(requested), n_present_genes = length(present),
            coverage_fraction = frac, present_genes = paste(present, collapse = ";"),
            missing_genes = paste(missing, collapse = ";"),
            status = ifelse(eligible, "PASS_SCORE_ELIGIBLE", "FAIL_NOT_SCORED_LOW_COVERAGE"),
            stringsAsFactors = FALSE
        )
        if (eligible) {
            sub <- expr_mat[present, , drop = FALSE]
            z <- t(apply(sub, 1L, .p8_z))
            sc <- colMeans(z, na.rm = TRUE)
            sc[!is.finite(sc)] <- NA_real_
            scores[[nm]] <- sc
        } else {
            scores[[nm]] <- NA_real_
        }
    }
    add_axis <- function(nms) {
        mat <- do.call(cbind, lapply(nms, function(nm) .p8_z(scores[[nm]])))
        out <- rowMeans(mat, na.rm = TRUE)
        out[rowSums(is.finite(mat)) < length(nms)] <- NA_real_
        out
    }
    scores$CAF_TGFb_axis <- add_axis(c("CAF_Core", "TGFb_Response"))
    scores$CAF_TGFb_ECM_axis <- add_axis(c("CAF_Core", "TGFb_Response", "ECM_Remodeling"))
    scores$Stromal_Remodeling_axis <- add_axis(c("CAF_Core", "TGFb_Response", "ECM_Remodeling", "EMT"))
    composite_cov <- data.frame(
        module = c("CAF_TGFb_axis", "CAF_TGFb_ECM_axis", "Stromal_Remodeling_axis"),
        n_requested_genes = NA_integer_, n_present_genes = NA_integer_, coverage_fraction = NA_real_,
        present_genes = "COMPOSITE_OF_ELIGIBLE_FROZEN_MODULE_SCORES", missing_genes = NA_character_,
        status = "COMPOSITE_REQUIRES_ALL_COMPONENTS", stringsAsFactors = FALSE
    )
    list(scores = scores, coverage = rbind(.p8_bind_rows(coverage_rows), composite_cov))
}

# Replace result-assuming publication prose with a neutral, result-driven summary.
.p8_publication_summary <- function(analysis_mol, orr_wilcox, orr_logistic, dcr_wilcox, dcr_logistic,
                                    known_binary_orr, known_numeric_orr, subgroup_results, claim, target_gene) {
    resp <- table(analysis_mol$response_clean, useNA = "ifany")
    orr <- table(analysis_mol$response_ORR_group, useNA = "ifany")
    dcr <- table(analysis_mol$response_DCR_group, useNA = "ifany")
    tw <- orr_wilcox[orr_wilcox$feature == "target_expr", , drop = FALSE]
    tl <- orr_logistic[orr_logistic$feature == "target_expr", , drop = FALSE]
    dw <- dcr_wilcox[dcr_wilcox$feature == "target_expr", , drop = FALSE]
    dl <- dcr_logistic[dcr_logistic$feature == "target_expr", , drop = FALSE]
    strict <- known_binary_orr[known_binary_orr$marker == "strong_ici_context_strict", , drop = FALSE]
    mss <- subgroup_results[subgroup_results$subset == "MSS_EBVneg_TMBlow_like_TRUE", , drop = FALSE]
    pdl1 <- known_numeric_orr[known_numeric_orr$feature == "PDL1_CPS", , drop = FALSE]
    inflamed5 <- known_numeric_orr[known_numeric_orr$feature %in% c("Inflamed5_score_value", "Inflamed_5gene_exploratory"), , drop = FALSE]
    tmb <- known_numeric_orr[known_numeric_orr$feature == "TMB_value", , drop = FALSE]
    data.frame(
        item = c(
            "part8_primary_dataset_policy", "tumor_only_evaluable_n", "response_counts_CR_PR_SD_PD",
            "ORR_counts_R_NR", "DCR_counts_DCR_PD", paste0(target_gene, "_ORR_wilcoxon_p"),
            paste0(target_gene, "_ORR_Firth_p"), paste0(target_gene, "_ORR_AUC"),
            paste0(target_gene, "_DCR_wilcoxon_p"), paste0(target_gene, "_DCR_Firth_p"),
            paste0(target_gene, "_DCR_AUC"), "strong_ici_context_strict_ORR_fisher_p",
            "strong_ici_context_strict_counts_TRUE_R_NR", "TMB_value_ORR_AUC", "PDL1_CPS_ORR_AUC",
            "Inflamed5_score_ORR_AUC", "MSS_EBVneg_TMBlow_like_subgroup", "final_claim_ceiling"
        ),
        value = c(
            "PRJEB25780/TIGER tumor-only is the only Part8 ICI cohort; Kim/Bagaev annotations are molecular context, not an external validation cohort.",
            as.character(nrow(analysis_mol)),
            paste0("CR=", resp[["CR"]] %||% 0, "; PR=", resp[["PR"]] %||% 0, "; SD=", resp[["SD"]] %||% 0, "; PD=", resp[["PD"]] %||% 0),
            paste0("R=", orr[["R"]] %||% 0, "; NR=", orr[["NR"]] %||% 0),
            paste0("DCR=", dcr[["DCR"]] %||% 0, "; PD=", dcr[["PD"]] %||% 0),
            .p8_fmt_num(if (nrow(tw)) tw$p_value[1] else NA),
            .p8_fmt_num(if (nrow(tl)) tl$p_value[1] else NA),
            .p8_fmt_num(if (nrow(tl)) tl$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(dw)) dw$p_value[1] else NA),
            .p8_fmt_num(if (nrow(dl)) dl$p_value[1] else NA),
            .p8_fmt_num(if (nrow(dl)) dl$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(strict)) strict$p_value[1] else NA),
            if (nrow(strict)) paste0("TRUE_R=", strict$marker_TRUE_positive[1], "; TRUE_NR=", strict$marker_TRUE_negative[1]) else NA_character_,
            .p8_fmt_num(if (nrow(tmb)) tmb$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(pdl1)) pdl1$auc_for_positive_high_value[1] else NA),
            .p8_fmt_num(if (nrow(inflamed5)) inflamed5$auc_for_positive_high_value[1] else NA),
            if (nrow(mss)) paste0("n=", mss$n[1], "; R=", mss$n_R[1], "; NR=", mss$n_NR[1], "; status=", mss$status[1]) else NA_character_,
            claim$status[claim$claim_domain == "recommended_claim_ceiling"][1]
        ),
        interpretation = c(
            "Single-cohort transportability/boundary analysis only.",
            "Tumor-only response-evaluable patient-level units.",
            "RECIST composition.", "Primary endpoint: ORR (CR/PR versus SD/PD).",
            "Exploratory endpoint: DCR (CR/PR/SD versus PD).",
            "Unadjusted rank-based group comparison for OLFML2B.",
            "Firth-penalized continuous OLFML2B association; exploratory.",
            "Rank-based discrimination; not internally validated predictive performance.",
            "Exploratory DCR rank comparison.", "Exploratory DCR Firth association.",
            "Exploratory DCR rank AUC.", "Known molecular-context positive-control comparison.",
            "Exact counts are prioritized over unstable effect estimates.",
            "Mutation burden context control.", "PD-L1 CPS context control.",
            "Equal-weight mean-z five-gene inflamed score; not a validated classifier.",
            "Evaluability of the non-strong-ICI molecular subgroup.", "Maximum permissible Part8 claim."
        ), stringsAsFactors = FALSE
    )
}

.p8_write_methods_results_text <- function(path, summary_tab, target_gene = "OLFML2B") {
    get <- function(item) {
        z <- summary_tab$value[summary_tab$item == item]
        if (length(z)) z[1] else NA_character_
    }
    txt <- c(
        "Part8 design and interpretation boundary",
        "",
        "Design: PRJEB25780/TIGER is analyzed as one anti-PD-1 gastric-cancer cohort. TIGER expression and response annotations define the analytic dataset. Kim/Bagaev MSI, EBV, mutation-burden and PD-L1 annotations are optional molecular-context variables and do not constitute an independent validation cohort.",
        "",
        paste0("Analysis set: ", get("tumor_only_evaluable_n"), " tumor-only patient-level units after excluding patient_name ending in -N. Response composition: ", get("response_counts_CR_PR_SD_PD"), ". ORR: ", get("ORR_counts_R_NR"), "; exploratory DCR: ", get("DCR_counts_DCR_PD"), "."),
        "",
        paste0("OLFML2B association: ORR Wilcoxon P=", get(paste0(target_gene, "_ORR_wilcoxon_p")), ", Firth P=", get(paste0(target_gene, "_ORR_Firth_p")), ", rank AUC=", get(paste0(target_gene, "_ORR_AUC")), "; DCR Wilcoxon P=", get(paste0(target_gene, "_DCR_wilcoxon_p")), ", Firth P=", get(paste0(target_gene, "_DCR_Firth_p")), ", rank AUC=", get(paste0(target_gene, "_DCR_AUC")), ". These are exploratory single-cohort associations, not validated predictive performance."),
        "",
        paste0("Molecular context: strict MSI-H/EBV-positive/TMB-high context Fisher P=", get("strong_ici_context_strict_ORR_fisher_p"), " with ", get("strong_ici_context_strict_counts_TRUE_R_NR"), ". TMB AUC=", get("TMB_value_ORR_AUC"), "; PD-L1 CPS AUC=", get("PDL1_CPS_ORR_AUC"), "; exploratory five-gene inflamed-score AUC=", get("Inflamed5_score_ORR_AUC"), ". The five-gene score is an equal-weight mean-z context control and is not a validated classifier."),
        "",
        paste0("MSS/EBV-negative/TMB-low-like subgroup: ", get("MSS_EBVneg_TMBlow_like_subgroup"), "."),
        "",
        paste0("Claim ceiling: ", get("final_claim_ceiling"))
    )
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(txt, con = path, useBytes = TRUE)
    invisible(path)
}

# Preserve the transformed analysis body as the core, then expose an
# OLFML2B-specific public interface with shared-data resolution and final QC.
.olfml2b_part8_core <- run_part8_ici_prjeb25780_tiger_only

run_olfml2b_part8_ici_prjeb25780 <- function(
    root = Sys.getenv("OLFML2B_STAD_ROOT", unset = "D:/OLFML2B_STAD"),
    target_gene = "OLFML2B",
    tiger_dir = NULL,
    shared_data_root = Sys.getenv("OLFML2B_STAD_SHARED_DATA_ROOT", unset = ""),
    molecular_file = NA_character_,
    expression_transform = c("auto", "none", "log2"),
    make_figures = TRUE,
    require_firth = TRUE
) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    if (!identical(toupper(target_gene), "OLFML2B")) .p8_stop("This Part8 release is gene-locked to OLFML2B.")
    if (isTRUE(require_firth) && !requireNamespace("logistf", quietly = TRUE)) {
        .p8_stop("Package 'logistf' is required. Install it once with install.packages('logistf'), then rerun Part8.")
    }
    if (nzchar(shared_data_root)) Sys.setenv(OLFML2B_STAD_SHARED_DATA_ROOT = shared_data_root)
    tiger_dir <- .p8_resolve_tiger_dir(root, tiger_dir, shared_data_root)
    out <- .olfml2b_part8_core(
        root = root, target_gene = "OLFML2B", tiger_dir = tiger_dir,
        molecular_file = molecular_file, expression_transform = expression_transform,
        make_figures = make_figures
    )
    unit_audit <- .p8_patient_unit_audit(out$analysis_input)
    .p8_write_csv(unit_audit, file.path(root, "output", "tables", "Part8_Immunotherapy", "06a_patient_level_inferential_unit_audit.csv"))
    if (any(unit_audit$status == "FAIL")) .p8_stop("Part8 patient-level inferential-unit audit failed.")
    qc <- .p8_structural_qc(root, out, target_gene = "OLFML2B")
    .p8_write_csv(qc, file.path(root, "output", "tables", "Part8_Immunotherapy", "31_part8_structural_qc.csv"))
    .p8_write_csv(.p8_output_manifest(file.path(root, "output", "tables", "Part8_Immunotherapy"), file.path(root, "output", "figures", "Part8_Immunotherapy")), file.path(root, "output", "tables", "Part8_Immunotherapy", "30_part8_output_manifest.csv"))
    if (any(qc$status == "FAIL")) .p8_stop("Part8 structural QC failed: ", paste(qc$criterion[qc$status == "FAIL"], collapse = "; "))
    out$version <- OLFML2B_PART8_VERSION
    out$shared_data_root <- shared_data_root
    out$patient_unit_audit <- unit_audit
    out$structural_qc <- qc
    .p8_save_rds(out, file.path(root, "output", "objects", "Part8_OLFML2B_Immunotherapy_TIGER_molecular_context.rds"))
    invisible(out)
}

run_part8_ici_prjeb25780_tiger_only <- run_olfml2b_part8_ici_prjeb25780
run_part8 <- run_olfml2b_part8_ici_prjeb25780
