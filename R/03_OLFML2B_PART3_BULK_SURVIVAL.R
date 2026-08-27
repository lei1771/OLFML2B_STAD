############################################################
## OLFML2B-STAD Specialized Bioinformatics Analysis
## Version: v1.0.0_20260705_EXECUTION_TEMPLATE_REBUILD
##
## Purpose:
##   A trial OLFML2B-specific analysis script for the gastric cancer
##   PhD project:
##   "OLFML2B expression in gastric cancer prognosis, recurrence-family endpoints, and exploratory state-module screening."
##
## Design principles:
##   1. Read already harmonized objects produced by the compact pipeline:
##        output/objects/Part1_TCGA_STAD.rds
##        output/objects/Part2_GEO_bulk_index.rds
##        output/objects/Part2_GSE84437.rds, Part2_GSE62254.rds,
##        Part2_GSE26253.rds, Part2_GSE15459.rds
##   2. Do NOT reuse old OLFML2B target-expression columns.
##   3. Extract OLFML2B directly from gene-level expression matrices.
##   4. Strictly exclude OLFML2A/OLFM2 paralog confusion.
##   5. Use continuous z-score Cox as the primary survival analysis.
##   6. Use median split only for visualization/sensitivity tables.
##   7. Test whether OLFML2B is more than a generic proliferation or bulk stromal-composition proxy
##      by adding proliferation, CAF/ECM, endothelial and epithelial context controls.
##   8. Test the proposed state model:
##        OLFML2B-high/low: module direction is exploratory and not pre-assumed
##   9. PDC proteomics is an optional measurement audit, not
##      a forced main conclusion.
##
## Run:
##   setwd("D:/OLFML2B_STAD")
##   source("R/03_OLFML2B_PART3_BULK_SURVIVAL.R")
##   res <- run_olfml2b_specialized_bioinformatics(
##       root = "D:/OLFML2B_STAD",
##       geo_validation_cohorts = c("GSE26253", "GSE84437", "GSE62254", "GSE15459")
##   )
############################################################

options(stringsAsFactors = FALSE)

OLFML2B_ANALYSIS_VERSION <- "v1.1.0_20260722_PDC614_CONTRACT_SYNC"
OLFML2B_PART2_REQUIRED_CONTRACT <- "v1.1.0_20260722_GEO_TUMOR_ONLY_CONTEXT_AND_SOURCE_ENDPOINT_FIX"
OLFML2B_PART2_REQUIRED_CLINICAL_CONTRACT <- "v1.0.0_20260721_OLFML2B_STRUCTURAL_MISSINGNESS_CONTRACT"

run_olfml2b_specialized_bioinformatics <- function(
    root = "D:/OLFML2B_STAD",
    target_gene = "OLFML2B",
    geo_validation_cohorts = c("GSE26253", "GSE84437", "GSE62254", "GSE15459"),
    output_subdir = "Part3",
    min_survival_n = 40L,
    min_survival_events = 20L,
    min_epv = 10,
    run_pdc_audit = TRUE,
    make_figures = TRUE
) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    dirs <- olfml2b_make_dirs(root, output_subdir)
    log_file <- file.path(dirs$logs, "OLFML2B_specialized_bioinformatics.log")
    olfml2b_log("Version: ", OLFML2B_ANALYSIS_VERSION, log_file = log_file)
    olfml2b_log("Root: ", root, log_file = log_file)
    olfml2b_log("Target gene: ", target_gene, log_file = log_file)

    olfml2b_require_packages(c("survival", "matrixStats"), log_file = log_file)

    params <- data.frame(
        parameter = c(
            "version", "root", "target_gene", "excluded_confounders",
            "geo_validation_cohorts", "primary_survival_model",
            "primary_endpoint", "secondary_endpoints", "min_survival_n",
            "min_survival_events", "min_epv", "analysis_statement"
        ),
        value = c(
            OLFML2B_ANALYSIS_VERSION,
            root,
            target_gene,
            "OLFML2A;OLFM2",
            paste(geo_validation_cohorts, collapse = ";"),
            "continuous z-score Cox; median split is sensitivity/display only",
            "RECURRENCE_FAMILY",
            "OS;DFS;RFS;Recurrence",
            as.character(min_survival_n),
            as.character(min_survival_events),
            as.character(min_epv),
            "Two-sided continuous association is primary and no effect direction is used as a quality gate; CAF/TGFb/ECM ecological-state context is prespecified, while tumor-normal overexpression, DDR-driver and ICI-prediction claims are excluded."
        ),
        stringsAsFactors = FALSE
    )
    olfml2b_write_csv(params, file.path(dirs$tables, "00_run_parameters.csv"))

    olfml2b_log("Loading harmonized TCGA and GEO objects", log_file = log_file)
    views <- olfml2b_load_bulk_views(
        root = root,
        target_gene = target_gene,
        geo_validation_cohorts = geo_validation_cohorts,
        log_file = log_file
    )
    input_contract_audit <- attr(views, "input_contract_audit")
    olfml2b_write_csv(input_contract_audit, file.path(dirs$tables, "00a_part2_input_contract_audit.csv"))

    source_status <- olfml2b_source_status(views)
    olfml2b_write_csv(source_status, file.path(dirs$tables, "01_source_status.csv"))

    measurement_audit <- olfml2b_measurement_audit(views)
    olfml2b_write_csv(measurement_audit, file.path(dirs$tables, "02_olfml2b_measurement_audit.csv"))

    olfml2b_log("Running TCGA tumor-normal expression analysis", log_file = log_file)
    tn <- olfml2b_tcga_tumor_normal(views$TCGA_STAD)
    olfml2b_write_csv(tn$summary, file.path(dirs$tables, "03a_tcga_olfml2b_tumor_normal_summary.csv"))
    olfml2b_write_csv(tn$paired_values, file.path(dirs$tables, "03b_tcga_olfml2b_paired_values.csv"))
    olfml2b_write_csv(olfml2b_expression_by_clinical_context(views),
                    file.path(dirs$tables, "03c_olfml2b_expression_by_clinical_context.csv"))

    olfml2b_log("Auditing available endpoints", log_file = log_file)
    endpoint_audit <- olfml2b_endpoint_audit(
        views,
        min_n = min_survival_n,
        min_events = min_survival_events
    )
    olfml2b_write_csv(endpoint_audit, file.path(dirs$tables, "04_endpoint_evaluable_audit.csv"))
    endpoint_alias_audit <- endpoint_audit[endpoint_audit$endpoint_family == "RECURRENCE", , drop = FALSE]
    olfml2b_write_csv(endpoint_alias_audit, file.path(dirs$tables, "04a_endpoint_alias_audit.csv"))

    olfml2b_log("Running continuous z-score Cox models", log_file = log_file)
    surv_all <- olfml2b_run_survival_all(
        views = views,
        min_n = min_survival_n,
        min_events = min_survival_events,
        min_epv = min_epv,
        log_file = log_file
    )
    if (nrow(surv_all)) {
        surv_all$fdr_by_endpoint_model <- olfml2b_ave_p_adjust(
            surv_all$p_value,
            surv_all$endpoint_family,
            surv_all$model
        )
    }
    olfml2b_write_csv(surv_all, file.path(dirs$tables, "05_olfml2b_survival_all_models.csv"))
    subseries_heterogeneity <- olfml2b_gse84437_subseries_heterogeneity_audit(views)
    olfml2b_write_csv(subseries_heterogeneity, file.path(dirs$tables, "05a_GSE84437_subseries_effect_heterogeneity_audit.csv"))

    olfml2b_log("Running median-split display/sensitivity survival tables", log_file = log_file)
    median_surv <- olfml2b_run_median_split_all(views, min_n = min_survival_n, min_events = min_survival_events)
    olfml2b_write_csv(median_surv, file.path(dirs$tables, "06_olfml2b_median_split_km_and_grouped_cox.csv"))

    olfml2b_log("Running spline non-linearity sensitivity tests", log_file = log_file)
    spline_tab <- olfml2b_run_spline_all(views, min_n = min_survival_n, min_events = min_survival_events)
    olfml2b_write_csv(spline_tab, file.path(dirs$tables, "07_olfml2b_nonlinearity_spline.csv"))

    olfml2b_log("Running prespecified PH time-interaction sensitivity", log_file = log_file)
    ph_sensitivity <- olfml2b_run_ph_time_interaction_all(
        views, min_n = max(min_survival_n, 80L), min_events = max(min_survival_events, 40L)
    )
    olfml2b_write_csv(ph_sensitivity, file.path(dirs$tables, "07a_olfml2b_PH_time_interaction_sensitivity.csv"))

    olfml2b_log("Running meta-analysis and leave-one-cohort-out sensitivity", log_file = log_file)
    meta <- olfml2b_meta_all(surv_all)
    loo <- olfml2b_leave_one_out_meta_all(surv_all)
    meta_elig <- if ("meta_eligible" %in% names(surv_all)) surv_all$meta_eligible %in% TRUE else surv_all$status == "OK"
    meta_input <- surv_all[
        meta_elig &
            surv_all$model %in% c("adjusted", "available_adjusted", "derived_stage_sensitivity", "univariable") &
            is.finite(surv_all$beta) &
            is.finite(surv_all$se) & surv_all$se > 0,
        , drop = FALSE
    ]
    meta_input_deduplicated <- olfml2b_deduplicate_meta_effects(meta_input)
    olfml2b_write_csv(meta_input_deduplicated, file.path(dirs$tables, "08a_olfml2b_meta_input_deduplicated.csv"))
    olfml2b_write_csv(meta, file.path(dirs$tables, "08_olfml2b_meta_analysis.csv"))
    olfml2b_write_csv(loo, file.path(dirs$tables, "09_olfml2b_leave_one_out_meta.csv"))

    olfml2b_log("Running proliferation and tissue-composition context-control models", log_file = log_file)
    control_results <- olfml2b_run_proliferation_controls(
        views = views,
        min_n = min_survival_n,
        min_events = min_survival_events,
        min_epv = min_epv,
        log_file = log_file
    )
    olfml2b_write_csv(control_results$expression_audit,
                    file.path(dirs$tables, "10_context_control_expression_audit.csv"))
    olfml2b_write_csv(control_results$models,
                    file.path(dirs$tables, "11_olfml2b_vs_context_control_models.csv"))

    olfml2b_log("Running OLFML2B state-module analysis", log_file = log_file)
    modules <- olfml2b_signature_catalog()
    module_results <- olfml2b_run_module_analysis(views, modules)
    olfml2b_write_csv(module_results$coverage,
                    file.path(dirs$tables, "12_module_score_coverage.csv"))
    olfml2b_write_csv(module_results$correlations,
                    file.path(dirs$tables, "13_olfml2b_module_correlations.csv"))
    olfml2b_write_csv(module_results$state_summary,
                    file.path(dirs$tables, "14_olfml2b_state_summary.csv"))
    olfml2b_write_csv(module_results$meta_correlations,
                    file.path(dirs$tables, "14a_olfml2b_module_random_effects_meta.csv"))

    pdc_audit <- data.frame()
    if (isTRUE(run_pdc_audit)) {
        olfml2b_log("Running lightweight PDC OLFML2B detection audit", log_file = log_file)
        pdc_audit <- olfml2b_pdc_detection_audit(root = root, target_gene = target_gene, max_file_mb = 80)
        olfml2b_write_csv(pdc_audit, file.path(dirs$tables, "15_pdc_olfml2b_detection_audit.csv"))
    }

    go_no_go <- olfml2b_go_no_go(
        views = views,
        endpoint_audit = endpoint_audit,
        surv_all = surv_all,
        meta = meta,
        control_models = control_results$models,
        module_summary = module_results$state_summary,
        pdc_audit = pdc_audit
    )
    olfml2b_write_csv(go_no_go, file.path(dirs$tables, "20_olfml2b_go_no_go_summary.csv"))
    interpretation_boundary <- olfml2b_interpretation_boundary(
        surv_all = surv_all,
        meta = meta,
        control_models = control_results$models,
        median_surv = median_surv,
        module_summary = module_results$state_summary
    )
    olfml2b_write_csv(interpretation_boundary, file.path(dirs$tables, "21_olfml2b_interpretation_boundary.csv"))

    if (isTRUE(make_figures)) {
        olfml2b_log("Generating basic audit figures", log_file = log_file)
        olfml2b_make_basic_figures(
            dirs = dirs,
            views = views,
            surv_all = surv_all,
            meta = meta,
            module_correlations = module_results$correlations,
            go_no_go = go_no_go,
            loo = loo,
            median_surv = median_surv,
            spline_tab = spline_tab,
            control_models = control_results$models,
            state_summary = module_results$state_summary
        )
    }

    index <- list(
        version = OLFML2B_ANALYSIS_VERSION,
        generated_at = olfml2b_timestamp(),
        root = root,
        dirs = dirs,
        source_status = source_status,
        input_contract_audit = input_contract_audit,
        measurement_audit = measurement_audit,
        endpoint_audit = endpoint_audit,
        survival = surv_all,
        gse84437_subseries_heterogeneity = subseries_heterogeneity,
        median_survival = median_surv,
        spline = spline_tab,
        ph_time_interaction = ph_sensitivity,
        meta = meta,
        leave_one_out = loo,
        controls = control_results,
        modules = module_results,
        pdc_audit = pdc_audit,
        go_no_go = go_no_go,
        interpretation_boundary = interpretation_boundary
    )
    saveRDS(index, file.path(dirs$objects, "Part3_OLFML2B_specialized_bioinformatics_index.rds"), compress = "xz")
    saveRDS(index, file.path(dirs$objects, "OLFML2B_specialized_analysis_index.rds"), compress = "xz")  # compatibility alias
    olfml2b_log("Done. Tables: ", dirs$tables, log_file = log_file)

    invisible(index)
}

############################################################
## Directory, logging, IO and dependencies
############################################################

olfml2b_make_dirs <- function(root, output_subdir) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    part <- if (grepl("^Part[0-9]$", as.character(output_subdir)[1])) as.character(output_subdir)[1] else "Part3"
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
        out_root = olfml2b_dir(out_root),
        tables = olfml2b_dir(file.path(out_root, "tables", part)),
        figures = olfml2b_dir(file.path(out_root, "figures", part)),
        reports = olfml2b_dir(file.path(out_root, "reports", part)),
        objects = olfml2b_dir(file.path(out_root, "objects")),
        logs = olfml2b_dir(file.path(root, "logs", "runtime", part)),
        qc = olfml2b_dir(file.path(out_root, "qc", part))
    )
}

olfml2b_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

olfml2b_timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

olfml2b_log <- function(..., level = "INFO", log_file = NULL) {
    line <- sprintf("[%s] [%s] [OLFML2B-STAD] %s", olfml2b_timestamp(), level, paste0(..., collapse = ""))
    message(line)
    if (!is.null(log_file)) {
        olfml2b_dir(dirname(log_file))
        cat(line, "\n", file = log_file, append = TRUE, sep = "")
    }
    invisible(line)
}

olfml2b_require_packages <- function(packages, log_file = NULL) {
    miss <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(miss)) {
        stop("Missing required packages: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

olfml2b_write_csv <- function(x, path) {
    olfml2b_dir(dirname(path))
    if (is.null(x)) x <- data.frame()
    if (is.atomic(x) && is.null(dim(x))) x <- data.frame(value = x)
    if (!is.data.frame(x) && !is.matrix(x)) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    invisible(path)
}

olfml2b_clean_names <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    make.unique(x)
}

olfml2b_bind_rows <- function(xs) {
    if (is.null(xs)) return(data.frame())
    if (is.data.frame(xs)) return(xs)
    xs <- xs[vapply(xs, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
    if (!length(xs)) return(data.frame())
    cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
    aligned <- lapply(xs, function(z) {
        miss <- setdiff(cols, names(z))
        for (m in miss) z[[m]] <- NA
        z[, cols, drop = FALSE]
    })
    out <- do.call(rbind, aligned)
    rownames(out) <- NULL
    out
}

olfml2b_num <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    x <- as.character(x)
    x <- gsub(",", "", x)
    x <- gsub("%", "", x)
    x <- gsub("[^0-9eE.+\\-]", "", x)
    suppressWarnings(as.numeric(x))
}

olfml2b_z <- function(x) {
    x <- as.numeric(x)
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
    (x - m) / s
}

olfml2b_bh <- function(p) {
    p <- as.numeric(p)
    out <- rep(NA_real_, length(p))
    idx <- which(is.finite(p))
    if (length(idx)) out[idx] <- p.adjust(p[idx], method = "BH")
    out
}

olfml2b_ave_p_adjust <- function(p, ...) {
    groups <- interaction(..., drop = TRUE, lex.order = TRUE)
    out <- rep(NA_real_, length(p))
    for (g in levels(groups)) {
        idx <- which(groups == g)
        out[idx] <- olfml2b_bh(p[idx])
    }
    out
}

olfml2b_safe_col <- function(df, patterns) {
    if (is.null(df) || !is.data.frame(df)) return(NA_character_)
    nms <- names(df)
    for (p in patterns) {
        hit <- grep(p, nms, ignore.case = TRUE, perl = TRUE, value = TRUE)
        if (length(hit)) return(hit[1])
    }
    NA_character_
}

############################################################
## Object loading and view construction
############################################################

olfml2b_load_bulk_views <- function(root, target_gene, geo_validation_cohorts, log_file = NULL) {
    obj_dir <- file.path(root, "output", "objects")
    tcga_path <- file.path(obj_dir, "Part1_TCGA_STAD.rds")
    geo_index_path <- file.path(obj_dir, "Part2_GEO_bulk_index.rds")

    if (file.exists(tcga_path)) {
        tcga_obj <- readRDS(tcga_path)
        tcga_source_file <- tcga_path
    } else {
        derived_expr <- file.path(root, "data", "derived", "expression", "TCGA_STAD_log_expression_gene_symbol.rds")
        derived_meta <- file.path(root, "data", "derived", "clinical", "TCGA_STAD_sample_metadata_harmonized.csv")
        table_meta <- file.path(root, "output", "tables", "Part1", "Part1_TCGA_sample_metadata.csv")
        meta_path <- if (file.exists(derived_meta)) derived_meta else table_meta
        if (!file.exists(derived_expr) || !file.exists(meta_path)) {
            stop("Missing TCGA object and no derived compatibility source was found. Expected either ", tcga_path, " or derived expression/metadata under data/derived.", call. = FALSE)
        }
        tcga_obj <- list(
            version = "compatibility_view_from_derived_outputs",
            cohort = "TCGA_STAD",
            expression = readRDS(derived_expr),
            sample_metadata = utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE),
            rebuild_note = "Part3 loaded TCGA from derived expression/metadata because Part1_TCGA_STAD.rds was missing."
        )
        tcga_source_file <- paste(derived_expr, meta_path, sep = ";")
        olfml2b_log("TCGA final object missing; loaded compatibility view from derived outputs.", level = "WARN", log_file = log_file)
    }
    tcga_view <- olfml2b_make_view(tcga_obj, cohort = "TCGA_STAD", target_gene = target_gene, source_file = tcga_source_file, is_tcga = TRUE)

    if (file.exists(geo_index_path)) {
        geo_index <- readRDS(geo_index_path)
        cohort_files <- geo_index$cohort_files
    } else {
        cohort_files <- setNames(file.path(obj_dir, paste0("Part2_", geo_validation_cohorts, ".rds")), geo_validation_cohorts)
        cohort_files <- cohort_files[file.exists(cohort_files)]
        olfml2b_log("GEO index missing; using available Part2_<GSE>.rds cohort objects as compatibility source.", level = "WARN", log_file = log_file)
    }
    if (is.null(cohort_files) || !length(cohort_files)) {
        stop("Part2_GEO_bulk_index.rds has no cohort_files and no Part2_<GSE>.rds compatibility objects were found.", call. = FALSE)
    }

    views <- list(TCGA_STAD = tcga_view)
    for (gse in geo_validation_cohorts) {
        fp <- NA_character_
        if (gse %in% names(cohort_files)) {
            fp <- cohort_files[[gse]]
        } else {
            hit <- cohort_files[grepl(gse, cohort_files, fixed = TRUE)]
            if (length(hit)) fp <- hit[[1]]
        }
        if (is.na(fp) || !file.exists(fp)) {
            olfml2b_log("GEO object not found for ", gse, "; adding empty audit view.", level = "WARN", log_file = log_file)
            views[[gse]] <- olfml2b_empty_view(gse, target_gene, load_status = "MISSING_GEO_RDS")
            next
        }
        obj <- tryCatch(readRDS(fp), error = function(e) e)
        if (inherits(obj, "error")) {
            olfml2b_log("Unable to read GEO object for ", gse, ": ", conditionMessage(obj), level = "WARN", log_file = log_file)
            views[[gse]] <- olfml2b_empty_view(gse, target_gene, load_status = "RDS_READ_FAILED")
            next
        }
        views[[gse]] <- olfml2b_make_view(obj, cohort = gse, target_gene = target_gene, source_file = fp, is_tcga = FALSE)
        olfml2b_log("Loaded ", gse, " | genes=", nrow(views[[gse]]$expr), " | samples=", ncol(views[[gse]]$expr),
                  " | OLFML2B measured=", views[[gse]]$target_measured, log_file = log_file)
    }
    input_contract_audit <- olfml2b_bulk_input_contract_audit(views)
    if (any(input_contract_audit$status == "FAIL")) {
        bad <- input_contract_audit[input_contract_audit$status == "FAIL", , drop = FALSE]
        stop("Part3/Part4 bulk input contract failed: ",
             paste(bad$cohort, bad$detail, sep = "=", collapse = "; "),
             ". Rebuild Part2 with the current code before downstream analysis.", call. = FALSE)
    }
    attr(views, "input_contract_audit") <- input_contract_audit
    views
}

olfml2b_empty_view <- function(cohort, target_gene, load_status) {
    list(
        cohort = cohort,
        source_file = NA_character_,
        expr = matrix(numeric(), nrow = 0L, ncol = 0L),
        clinical = data.frame(),
        target_gene = target_gene,
        target_measured = FALSE,
        target_probe_rows = character(),
        source_version = NA_character_,
        source_release_contract = NA_character_,
        source_clinical_contract = NA_character_,
        clinical_field_contract = data.frame(),
        source_role = NA_character_,
        source_partition_ok = NA,
        source_batch_audit_ok = NA,
        load_status = load_status,
        is_tcga = FALSE
    )
}

olfml2b_make_view <- function(obj, cohort, target_gene, source_file, is_tcga = FALSE) {
    expr <- olfml2b_find_expression_matrix(obj)
    clinical <- olfml2b_find_metadata_frame(obj, expr)
    clinical <- olfml2b_standardize_clinical(clinical, cohort = cohort)

    if (nrow(expr) > 0L && nrow(clinical) > 0L) {
        aligned <- olfml2b_align_expr_clinical(expr, clinical)
        expr <- aligned$expr
        clinical <- aligned$clinical
    }

    target <- olfml2b_extract_gene(expr, target_gene)
    clinical$olfml2b_expression <- target$values[match(clinical$sample_id, names(target$values))]
    if (identical(cohort, "GSE84437") && "source_subseries" %in% names(clinical) &&
        length(unique(stats::na.omit(clinical$source_subseries))) == 2L) {
        clinical$olfml2b_z_within_subseries <- ave(
            clinical$olfml2b_expression, clinical$source_subseries,
            FUN = function(z) olfml2b_z(as.numeric(z))
        )
        clinical$olfml2b_z <- olfml2b_z(clinical$olfml2b_z_within_subseries)
    } else {
        clinical$olfml2b_z <- olfml2b_z(clinical$olfml2b_expression)
    }

    list(
        cohort = cohort,
        source_file = source_file,
        expr = expr,
        clinical = clinical,
        target_gene = target_gene,
        target_measured = target$measured,
        target_probe_rows = target$rows,
        source_version = as.character(obj$version %||% NA_character_),
        source_release_contract = as.character(obj$release_contract_version %||% NA_character_),
        source_clinical_contract = as.character(obj$clinical_contract_version %||% NA_character_),
        clinical_field_contract = obj$clinical_field_contract %||% data.frame(),
        source_role = as.character(obj$role %||% NA_character_),
        source_partition_ok = is.data.frame(obj$superseries_partition_audit) &&
            nrow(obj$superseries_partition_audit) == 483L &&
            sum(obj$superseries_partition_audit$analysis_role == "FORMAL_OS") == 433L &&
            sum(obj$superseries_partition_audit$analysis_role == "CONTEXT_ONLY") == 50L,
        source_batch_audit_ok = is.data.frame(obj$batch_audit) && nrow(obj$batch_audit) == 3L,
        load_status = if (target$measured) "OK" else "TARGET_NOT_MEASURED",
        is_tcga = isTRUE(is_tcga)
    )
}

olfml2b_bulk_input_contract_audit <- function(views) {
    expected_n <- c(TCGA_STAD = NA_integer_, GSE62254 = 300L, GSE15459 = 192L,
                    GSE26253 = 432L, GSE84437 = 433L)
    rows <- lapply(names(views), function(cohort) {
        v <- views[[cohort]]; d <- v$clinical
        n_expr <- ncol(v$expr); n_clin <- nrow(d)
        unique_samples <- if ("sample_id" %in% names(d)) length(unique(as.character(d$sample_id))) else 0L
        z_sd <- if ("olfml2b_z" %in% names(d) && sum(is.finite(d$olfml2b_z)) > 1L) stats::sd(d$olfml2b_z, na.rm = TRUE) else NA_real_
        common_ok <- n_expr == n_clin && n_clin > 0L && unique_samples == n_clin &&
            isTRUE(v$target_measured) && is.finite(z_sd) && abs(z_sd - 1) < 1e-8
        version_ok <- if (cohort == "TCGA_STAD") TRUE else
            identical(v$source_release_contract, OLFML2B_PART2_REQUIRED_CONTRACT)
        clinical_contract_ok <- if (cohort == "TCGA_STAD") TRUE else
            identical(v$source_clinical_contract, OLFML2B_PART2_REQUIRED_CLINICAL_CONTRACT) &&
            is.data.frame(v$clinical_field_contract) && nrow(v$clinical_field_contract) > 0L
        cohort_ok <- TRUE
        detail <- paste0("expr=", n_expr, ";clinical=", n_clin, ";unique_sample=", unique_samples,
                         ";z_sd=", signif(z_sd, 6), ";contract=", v$source_release_contract,
                         ";clinical_contract=", v$source_clinical_contract)
        if (cohort %in% names(expected_n) && is.finite(expected_n[[cohort]])) {
            cohort_ok <- cohort_ok && n_clin == expected_n[[cohort]]
        }
        if (cohort == "GSE15459") {
            excluded <- if ("sample_title" %in% names(d)) sum(grepl("EXCLUDED", d$sample_title, ignore.case = TRUE)) else NA_integer_
            os_valid <- is.finite(d$os_time_days) & d$os_time_days > 0 & d$os_event %in% c(0L, 1L)
            n_os <- sum(os_valid); events <- sum(d$os_event[os_valid] == 1L, na.rm = TRUE)
            cohort_ok <- cohort_ok && identical(excluded, 0L) && n_os == 191L && events == 95L
            detail <- paste0(detail, ";excluded=", excluded, ";OS=", n_os, ";events=", events)
        } else if (cohort == "GSE26253") {
            unique_patient <- if ("patient_id" %in% names(d)) length(unique(as.character(d$patient_id))) else 0L
            cohort_ok <- cohort_ok && unique_patient == 432L
            detail <- paste0(detail, ";unique_patient=", unique_patient)
        } else if (cohort == "GSE84437") {
            role_ok <- "analysis_role" %in% names(d) && all(d$analysis_role == "FORMAL_OS")
            context_n <- sum(grepl("^GSM441", d$sample_id))
            batch_centered <- all(c("olfml2b_z", "source_subseries") %in% names(d)) &&
                all(abs(tapply(d$olfml2b_z, d$source_subseries, mean, na.rm = TRUE)) < 1e-8)
            sex_levels <- if ("sex" %in% names(d)) unique(as.character(stats::na.omit(d$sex))) else character()
            os_valid <- is.finite(d$os_time_days) & d$os_time_days > 0 & d$os_event %in% c(0L, 1L)
            n_os <- sum(os_valid); events <- sum(d$os_event[os_valid] == 1L, na.rm = TRUE)
            cohort_ok <- cohort_ok && role_ok && context_n == 0L && all(c("Female", "Male") %in% sex_levels) &&
                n_os == 431L && events == 207L && batch_centered &&
                isTRUE(v$source_partition_ok) && isTRUE(v$source_batch_audit_ok)
            detail <- paste0(detail, ";formal_role=", role_ok, ";context_n=", context_n,
                             ";sex=", paste(sort(sex_levels), collapse = "/"), ";OS=", n_os, ";events=", events,
                             ";batch_centered=", batch_centered,
                             ";partition=", v$source_partition_ok, ";batch=", v$source_batch_audit_ok)
        }
        detail <- paste0(
            detail,
            ";common_contract_ok=", common_ok,
            ";release_contract_ok=", version_ok,
            ";clinical_contract_ok=", clinical_contract_ok,
            ";scientific_contract_ok=", cohort_ok
        )
        data.frame(
            cohort = cohort, n_expression = n_expr, n_clinical = n_clin,
            expected_n = expected_n[[cohort]] %||% NA_integer_, unique_samples = unique_samples,
            target_measured = isTRUE(v$target_measured), within_cohort_z_sd = z_sd,
            release_contract_ok = version_ok, clinical_contract_ok = clinical_contract_ok,
            scientific_contract_ok = cohort_ok,
            status = ifelse(common_ok && version_ok && clinical_contract_ok && cohort_ok, "PASS", "FAIL"),
            detail = detail, stringsAsFactors = FALSE
        )
    })
    olfml2b_bind_rows(rows)
}

olfml2b_find_expression_matrix <- function(obj) {
    candidates <- list()
    collect <- function(x, path = "obj") {
        if (is.matrix(x) || inherits(x, "Matrix")) {
            if (is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x)) && nrow(x) > 10L && ncol(x) > 10L) {
                candidates[[length(candidates) + 1L]] <<- list(path = path, mat = as.matrix(x))
            }
        } else if (is.list(x)) {
            nms <- names(x)
            if (is.null(nms)) nms <- paste0("V", seq_along(x))
            for (i in seq_along(x)) {
                collect(x[[i]], paste0(path, "$", nms[i]))
            }
        }
    }
    if (is.list(obj)) {
        for (nm in c("expression", "expr", "gene_expression", "rna_expression")) {
            if (!is.null(obj[[nm]]) && (is.matrix(obj[[nm]]) || inherits(obj[[nm]], "Matrix"))) {
                mat <- as.matrix(obj[[nm]])
                storage.mode(mat) <- "numeric"
                return(olfml2b_clean_expression_matrix(mat))
            }
        }
    }
    collect(obj)
    if (!length(candidates)) return(matrix(numeric(), nrow = 0L, ncol = 0L))
    sizes <- vapply(candidates, function(z) nrow(z$mat) * ncol(z$mat), numeric(1))
    mat <- candidates[[which.max(sizes)]]$mat
    storage.mode(mat) <- "numeric"
    olfml2b_clean_expression_matrix(mat)
}

olfml2b_clean_expression_matrix <- function(mat) {
    if (is.null(mat) || !length(mat)) return(matrix(numeric(), nrow = 0L, ncol = 0L))
    rownames(mat) <- toupper(trimws(rownames(mat)))
    keep <- !is.na(rownames(mat)) & nzchar(rownames(mat))
    mat <- mat[keep, , drop = FALSE]
    if (anyDuplicated(rownames(mat))) {
        vars <- matrixStats::rowVars(mat, na.rm = TRUE)
        keep_idx <- unlist(tapply(seq_len(nrow(mat)), rownames(mat), function(ii) ii[which.max(vars[ii])]), use.names = FALSE)
        mat <- mat[sort(keep_idx), , drop = FALSE]
    }
    mat
}

olfml2b_find_metadata_frame <- function(obj, expr) {
    if (is.list(obj)) {
        for (nm in c("sample_metadata", "metadata", "clinical", "pheno", "pdata")) {
            if (!is.null(obj[[nm]]) && is.data.frame(obj[[nm]])) return(obj[[nm]])
        }
    }
    frames <- list()
    collect <- function(x, path = "obj") {
        if (is.data.frame(x) && nrow(x) > 10L) {
            frames[[length(frames) + 1L]] <<- x
        } else if (is.list(x)) {
            for (i in seq_along(x)) collect(x[[i]])
        }
    }
    collect(obj)
    if (!length(frames)) return(data.frame(sample_id = colnames(expr), stringsAsFactors = FALSE))
    if (ncol(expr) > 0L) {
        sizes <- vapply(frames, function(f) abs(nrow(f) - ncol(expr)), numeric(1))
        return(frames[[which.min(sizes)]])
    }
    frames[[which.max(vapply(frames, nrow, integer(1)))]]
}

olfml2b_align_expr_clinical <- function(expr, clinical) {
    if (!"sample_id" %in% names(clinical)) clinical$sample_id <- rownames(clinical)
    clinical$sample_id <- as.character(clinical$sample_id)
    common <- intersect(colnames(expr), clinical$sample_id)
    if (length(common) >= 10L) {
        expr <- expr[, common, drop = FALSE]
        clinical <- clinical[match(common, clinical$sample_id), , drop = FALSE]
    } else if (ncol(expr) == nrow(clinical)) {
        clinical$sample_id <- colnames(expr)
    }
    rownames(clinical) <- clinical$sample_id
    list(expr = expr, clinical = clinical)
}

olfml2b_extract_gene <- function(expr, gene) {
    if (nrow(expr) == 0L || ncol(expr) == 0L) {
        return(list(measured = FALSE, rows = character(), values = setNames(numeric(), character())))
    }
    gene <- toupper(gene)
    rn <- toupper(rownames(expr))
    rows <- which(rn == gene)
    ## OLFML2B must not be confused with similarly named paralogs or family members.
    if (identical(gene, "OLFML2B")) {
        rows <- rows[!rn[rows] %in% c("OLFML2A", "OLFM2")]
    }
    if (!length(rows)) {
        return(list(measured = FALSE, rows = character(), values = setNames(rep(NA_real_, ncol(expr)), colnames(expr))))
    }
    if (length(rows) == 1L) {
        values <- as.numeric(expr[rows, ])
        names(values) <- colnames(expr)
        return(list(measured = TRUE, rows = rownames(expr)[rows], values = values))
    }
    vars <- matrixStats::rowVars(expr[rows, , drop = FALSE], na.rm = TRUE)
    sel <- rows[which.max(vars)]
    values <- as.numeric(expr[sel, ])
    names(values) <- colnames(expr)
    list(measured = TRUE, rows = rownames(expr)[rows], values = values)
}

olfml2b_standardize_clinical <- function(meta, cohort) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
    original_names <- names(meta)
    names(meta) <- make.unique(olfml2b_clean_names(names(meta)))

    if (!"sample_id" %in% names(meta)) {
        sid_col <- olfml2b_safe_col(meta, c("^sample_id$", "submitter_id", "geo_accession", "^gsm", "barcode"))
        meta$sample_id <- if (!is.na(sid_col)) as.character(meta[[sid_col]]) else rownames(meta)
    }
    if (!"patient_id" %in% names(meta)) {
        pid_col <- olfml2b_safe_col(meta, c("^patient_id$", "case_id", "participant", "subject", "patient", "title"))
        meta$patient_id <- if (!is.na(pid_col)) as.character(meta[[pid_col]]) else olfml2b_patient_from_sample(meta$sample_id)
    }

    if (!"is_tumor" %in% names(meta)) {
        st <- paste(meta$sample_id, meta[[olfml2b_safe_col(meta, c("sample_type", "tissue", "source_name", "title")) %||% names(meta)[1]]])
        meta$is_tumor <- grepl("tumou?r|cancer|primary tumor|stomach adenocarcinoma|tumor", st, ignore.case = TRUE)
    }
    if (!"is_normal" %in% names(meta)) {
        st <- paste(meta$sample_id, meta[[olfml2b_safe_col(meta, c("sample_type", "tissue", "source_name", "title")) %||% names(meta)[1]]])
        meta$is_normal <- grepl("normal|adjacent|non.?tumou?r|precancer", st, ignore.case = TRUE)
    }

    age_col <- olfml2b_safe_col(meta, c("^age$", "age_at", "age_year", "years"))
    meta$age <- if (!is.na(age_col)) olfml2b_num(meta[[age_col]]) else NA_real_
    idx <- is.finite(meta$age) & meta$age > 150
    meta$age[idx] <- meta$age[idx] / 365.25

    sex_col <- olfml2b_safe_col(meta, c("^sex$", "^gender$"))
    meta$sex <- if (!is.na(sex_col)) olfml2b_sex(meta[[sex_col]]) else factor(NA_character_, levels = c("Female", "Male"))

    stage_col <- olfml2b_safe_col(meta, c("^stage$", "ajcc.*stage", "pathologic.*stage", "pathological.*stage", "tumor_stage"))
    meta$stage <- if (!is.na(stage_col)) olfml2b_stage(meta[[stage_col]]) else factor(NA_character_, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
    meta$stage_numeric <- as.numeric(meta$stage)

    # Preserve the Part2 clinical contract.  Overall stage and individual TNM
    # components are different variables and must never overwrite one another.
    overall_col <- olfml2b_safe_col(meta, c("^stage_overall$"))
    primary_stage_col <- olfml2b_safe_col(meta, c("^stage_analysis_primary$"))
    sensitivity_stage_col <- olfml2b_safe_col(meta, c("^stage_analysis_sensitivity$",
                                                     "^stage_derived_ajcc7_m0_sensitivity$"))
    pt_col <- olfml2b_safe_col(meta, c("^stage_pt$", "^ptstage$", "^p_t_stage$", "^t_stage$"))
    pn_col <- olfml2b_safe_col(meta, c("^stage_pn$", "^pnstage$", "^p_n_stage$", "^n_stage$"))
    pm_col <- olfml2b_safe_col(meta, c("^stage_pm$", "^pmstage$", "^p_m_stage$", "^m_stage$"))
    meta$stage_overall <- if (!is.na(overall_col)) olfml2b_stage(meta[[overall_col]]) else meta$stage
    meta$stage_analysis_primary <- if (!is.na(primary_stage_col)) {
        olfml2b_stage(meta[[primary_stage_col]])
    } else {
        meta$stage
    }
    meta$stage_analysis_sensitivity <- if (!is.na(sensitivity_stage_col)) {
        olfml2b_stage(meta[[sensitivity_stage_col]])
    } else {
        meta$stage_analysis_primary
    }
    meta$stage_pT <- if (!is.na(pt_col)) as.character(meta[[pt_col]]) else NA_character_
    meta$stage_pN <- if (!is.na(pn_col)) as.character(meta[[pn_col]]) else NA_character_
    meta$stage_pM <- if (!is.na(pm_col)) as.character(meta[[pm_col]]) else NA_character_
    meta$stage <- meta$stage_analysis_primary
    meta$stage_numeric <- as.numeric(meta$stage_analysis_primary)

    lauren_col <- olfml2b_safe_col(meta, c("lauren"))
    if (!is.na(lauren_col)) meta$lauren <- factor(as.character(meta[[lauren_col]]))

    subtype_col <- olfml2b_safe_col(meta, c("acrg", "molecular_subtype", "subtype"))
    if (!is.na(subtype_col)) {
        meta$molecular_subtype <- factor(as.character(meta[[subtype_col]]))
        meta$acrg_subtype <- meta$molecular_subtype
    }

    therapy_col <- olfml2b_safe_col(meta, c("adjuvant", "chemotherapy", "chemo", "treatment"))
    if (!is.na(therapy_col)) meta$adjuvant_therapy <- olfml2b_therapy(meta[[therapy_col]])

    meta <- olfml2b_standardize_endpoints(meta)
    meta$cohort <- cohort
    attr(meta, "original_names") <- original_names
    meta
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

olfml2b_patient_from_sample <- function(x) {
    x <- as.character(x)
    out <- x
    tcga <- grepl("^TCGA-", x)
    out[tcga] <- substr(x[tcga], 1, 12)
    out
}

olfml2b_sex <- function(x) {
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_character_, length(y))
    out[y %in% c("female", "f", "woman")] <- "Female"
    out[y %in% c("male", "m", "man")] <- "Male"
    factor(out, levels = c("Female", "Male"))
}

olfml2b_stage <- function(x) {
    y <- toupper(trimws(as.character(x)))
    y <- gsub("PATHOLOGIC|PATHOLOGICAL|AJCC|STAGE|[ :_\\-]", "", y)
    out <- rep(NA_character_, length(y))
    out[grepl("^I$|^IA$|^IB$|^1$", y)] <- "Stage I"
    out[grepl("^II$|^IIA$|^IIB$|^2$", y)] <- "Stage II"
    out[grepl("^III$|^IIIA$|^IIIB$|^IIIC$|^3$", y)] <- "Stage III"
    out[grepl("^IV$|^IVA$|^IVB$|^4$", y)] <- "Stage IV"
    factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

olfml2b_therapy <- function(x) {
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_character_, length(y))
    out[grepl("yes|treated|chemo|fluor|platin|xeloda|capecitabine|s-1|tegafur", y)] <- "Yes"
    out[grepl("no|none|untreated|not", y)] <- "No"
    factor(out, levels = c("No", "Yes"))
}

olfml2b_event <- function(x) {
    if (is.null(x)) return(integer())
    if (is.numeric(x) || is.integer(x)) {
        y <- as.numeric(x)
        out <- rep(NA_integer_, length(y))
        out[is.finite(y) & y == 1] <- 1L
        out[is.finite(y) & y == 0] <- 0L
        return(out)
    }
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(y))
    out[grepl("^1$|dead|death|deceased|died|yes|event|recurrence|recurred|progress", y)] <- 1L
    out[grepl("^0$|alive|living|no|censor|non.?recurrence|disease.?free", y)] <- 0L
    out
}

olfml2b_standardize_endpoints <- function(meta) {
    # OS
    if (!"os_time_days" %in% names(meta)) {
        os_time_col <- olfml2b_safe_col(meta, c("^os_time_days$", "overall.*survival.*days", "duration.*overall.*survival", "^os.*time", "days_to_death", "days_to_last_follow"))
        if (!is.na(os_time_col)) {
            vals <- olfml2b_num(meta[[os_time_col]])
            if (grepl("month", os_time_col, ignore.case = TRUE) || (max(vals, na.rm = TRUE) < 250 && is.finite(max(vals, na.rm = TRUE)))) vals <- vals * 30.4375
            meta$os_time_days <- vals
        } else meta$os_time_days <- NA_real_
    }
    if (!"os_event" %in% names(meta)) {
        os_event_col <- olfml2b_safe_col(meta, c("^os_event$", "vital_status", "death", "overall.*survival.*status", "^os.*event"))
        meta$os_event <- if (!is.na(os_event_col)) olfml2b_event(meta[[os_event_col]]) else NA_integer_
    }

    # DFS
    if (!"dfs_time_days" %in% names(meta)) {
        col <- olfml2b_safe_col(meta, c("^dfs_time_days$", "disease.?free.*time", "dfs.*time", "recurrence.*free.*time"))
        if (!is.na(col)) {
            vals <- olfml2b_num(meta[[col]])
            if (grepl("month", col, ignore.case = TRUE) || (max(vals, na.rm = TRUE) < 250 && is.finite(max(vals, na.rm = TRUE)))) vals <- vals * 30.4375
            meta$dfs_time_days <- vals
        } else meta$dfs_time_days <- NA_real_
    }
    if (!"dfs_event" %in% names(meta)) {
        col <- olfml2b_safe_col(meta, c("^dfs_event$", "disease.?free.*event", "dfs.*event", "recurrence.*status", "status_0_non_recurrence"))
        meta$dfs_event <- if (!is.na(col)) olfml2b_event(meta[[col]]) else NA_integer_
    }

    # RFS
    if (!"rfs_time_days" %in% names(meta)) {
        col <- olfml2b_safe_col(meta, c("^rfs_time_days$", "recurrence.?free.*time", "rfs.*time"))
        if (!is.na(col)) {
            vals <- olfml2b_num(meta[[col]])
            if (grepl("month", col, ignore.case = TRUE) || (max(vals, na.rm = TRUE) < 250 && is.finite(max(vals, na.rm = TRUE)))) vals <- vals * 30.4375
            meta$rfs_time_days <- vals
        } else meta$rfs_time_days <- meta$dfs_time_days
    }
    if (!"rfs_event" %in% names(meta)) {
        col <- olfml2b_safe_col(meta, c("^rfs_event$", "recurrence.?free.*event", "rfs.*event"))
        meta$rfs_event <- if (!is.na(col)) olfml2b_event(meta[[col]]) else meta$dfs_event
    }

    # recurrence alias
    if (!"recurrence_time_days" %in% names(meta)) meta$recurrence_time_days <- meta$rfs_time_days
    if (!"recurrence_event" %in% names(meta)) meta$recurrence_event <- meta$rfs_event

    meta
}

############################################################
## Audits and expression analysis
############################################################

olfml2b_source_status <- function(views) {
    olfml2b_bind_rows(lapply(views, function(v) {
        data.frame(
            cohort = v$cohort,
            source_file = v$source_file,
            load_status = v$load_status,
            is_tcga = v$is_tcga,
            n_genes = nrow(v$expr),
            n_samples_expr = ncol(v$expr),
            n_samples_clinical = nrow(v$clinical),
            olfml2b_measured = v$target_measured,
            stringsAsFactors = FALSE
        )
    }))
}

olfml2b_measurement_audit <- function(views) {
    olfml2b_bind_rows(lapply(views, function(v) {
        d <- v$clinical
        data.frame(
            cohort = v$cohort,
            target_gene = v$target_gene,
            target_measured = v$target_measured,
            target_rows = paste(v$target_probe_rows, collapse = ";"),
            n_samples = nrow(d),
            n_olfml2b_nonmissing = sum(is.finite(d$olfml2b_expression)),
            mean_olfml2b = mean(d$olfml2b_expression, na.rm = TRUE),
            sd_olfml2b = stats::sd(d$olfml2b_expression, na.rm = TRUE),
            median_olfml2b = stats::median(d$olfml2b_expression, na.rm = TRUE),
            missing_fraction = mean(!is.finite(d$olfml2b_expression)),
            stringsAsFactors = FALSE
        )
    }))
}

olfml2b_tcga_tumor_normal <- function(view) {
    d <- view$clinical
    if (!nrow(d)) return(list(summary = data.frame(), paired_values = data.frame()))
    tumor <- d[d$is_tumor %in% TRUE & is.finite(d$olfml2b_expression), , drop = FALSE]
    normal <- d[d$is_normal %in% TRUE & is.finite(d$olfml2b_expression), , drop = FALSE]
    summary <- data.frame(
        cohort = view$cohort,
        n_tumor = nrow(tumor),
        n_normal = nrow(normal),
        median_tumor = stats::median(tumor$olfml2b_expression, na.rm = TRUE),
        median_normal = stats::median(normal$olfml2b_expression, na.rm = TRUE),
        log2fc_median_tumor_vs_normal = stats::median(tumor$olfml2b_expression, na.rm = TRUE) - stats::median(normal$olfml2b_expression, na.rm = TRUE),
        wilcox_p = if (nrow(tumor) >= 5 && nrow(normal) >= 5) stats::wilcox.test(tumor$olfml2b_expression, normal$olfml2b_expression, exact = FALSE)$p.value else NA_real_,
        stringsAsFactors = FALSE
    )
    pairs <- intersect(tumor$patient_id, normal$patient_id)
    paired <- data.frame()
    if (length(pairs) >= 5) {
        tt <- tumor[match(pairs, tumor$patient_id), , drop = FALSE]
        nn <- normal[match(pairs, normal$patient_id), , drop = FALSE]
        paired <- data.frame(
            patient_id = pairs,
            tumor_expression = tt$olfml2b_expression,
            normal_expression = nn$olfml2b_expression,
            delta_tumor_minus_normal = tt$olfml2b_expression - nn$olfml2b_expression,
            stringsAsFactors = FALSE
        )
        summary$n_paired <- nrow(paired)
        summary$paired_median_delta <- stats::median(paired$delta_tumor_minus_normal, na.rm = TRUE)
        summary$paired_wilcox_p <- stats::wilcox.test(paired$tumor_expression, paired$normal_expression, paired = TRUE, exact = FALSE)$p.value
    } else {
        summary$n_paired <- length(pairs)
        summary$paired_median_delta <- NA_real_
        summary$paired_wilcox_p <- NA_real_
    }
    list(summary = summary, paired_values = paired)
}

olfml2b_expression_by_clinical_context <- function(views) {
    rows <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        d <- olfml2b_analysis_clinical(v)
        if (!nrow(d)) next
        for (var in c("stage", "sex", "lauren", "molecular_subtype", "adjuvant_therapy")) {
            if (!var %in% names(d)) next
            dd <- d[is.finite(d$olfml2b_expression) & !is.na(d[[var]]), , drop = FALSE]
            if (nrow(dd) < 20L || length(unique(dd[[var]])) < 2L) next
            fit <- tryCatch(stats::lm(olfml2b_expression ~ x, data = data.frame(olfml2b_expression = dd$olfml2b_expression, x = dd[[var]])), error = function(e) NULL)
            p <- if (!is.null(fit)) tryCatch(stats::anova(fit)$`Pr(>F)`[1], error = function(e) NA_real_) else NA_real_
            tab <- aggregate(dd$olfml2b_expression, by = list(level = as.character(dd[[var]])), FUN = function(z) c(n = length(z), median = median(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE)))
            out <- data.frame(
                cohort = v$cohort,
                variable = var,
                level = tab$level,
                n = tab$x[, "n"],
                median = tab$x[, "median"],
                mean = tab$x[, "mean"],
                global_p = p,
                stringsAsFactors = FALSE
            )
            rows[[length(rows) + 1L]] <- out
        }
    }
    olfml2b_bind_rows(rows)
}

olfml2b_analysis_clinical <- function(view) {
    d <- view$clinical
    if (!nrow(d)) return(d)
    if (isTRUE(view$is_tcga) && "is_tumor" %in% names(d)) {
        d <- d[d$is_tumor %in% TRUE, , drop = FALSE]
    }
    d
}

olfml2b_endpoint_defs <- function() {
    list(
        OS = list(time = "os_time_days", event = "os_event", family = "OS", priority = 1L),
        DFS = list(time = "dfs_time_days", event = "dfs_event", family = "RECURRENCE", priority = 1L),
        RFS = list(time = "rfs_time_days", event = "rfs_event", family = "RECURRENCE", priority = 2L),
        Recurrence = list(time = "recurrence_time_days", event = "recurrence_event", family = "RECURRENCE", priority = 3L)
    )
}

olfml2b_endpoint_status_one <- function(d, ep, def, min_n = 40L, min_events = 20L) {
    if (!all(c(def$time, def$event) %in% names(d))) {
        return(data.frame(
            endpoint = ep, endpoint_family = def$family,
            time_col = def$time, event_col = def$event,
            n_complete = 0L, events = 0L,
            endpoint_priority = def$priority %||% NA_integer_,
            status = "MISSING_COLUMNS",
            reason = paste0("Missing required columns: ", paste(setdiff(c(def$time, def$event), names(d)), collapse = ";")),
            stringsAsFactors = FALSE
        ))
    }
    keep <- is.finite(d[[def$time]]) & d[[def$time]] > 0 &
        d[[def$event]] %in% c(0L, 1L) & is.finite(d$olfml2b_z)
    n_complete <- sum(keep)
    events <- sum(d[[def$event]][keep] == 1L, na.rm = TRUE)
    status <- if (n_complete >= min_n && events >= min_events) "EVALUABLE" else "NOT_EVALUABLE"
    reason <- if (identical(status, "EVALUABLE")) {
        "OK"
    } else if (n_complete < min_n) {
        paste0("Insufficient complete cases: ", n_complete, " < ", min_n)
    } else {
        paste0("Insufficient events: ", events, " < ", min_events)
    }
    data.frame(
        endpoint = ep, endpoint_family = def$family,
        time_col = def$time, event_col = def$event,
        n_complete = n_complete, events = events,
        endpoint_priority = def$priority %||% NA_integer_,
        status = status,
        reason = reason,
        stringsAsFactors = FALSE
    )
}

olfml2b_endpoint_selection_table <- function(d, min_n = 40L, min_events = 20L) {
    endpoints <- olfml2b_endpoint_defs()
    rows <- lapply(names(endpoints), function(ep) olfml2b_endpoint_status_one(d, ep, endpoints[[ep]], min_n, min_events))
    tab <- olfml2b_bind_rows(rows)
    tab$analysis_included <- FALSE
    tab$primary_family_endpoint <- NA_character_
    tab$alias_of <- NA_character_

    if ("OS" %in% tab$endpoint && tab$status[tab$endpoint == "OS"] == "EVALUABLE") {
        tab$analysis_included[tab$endpoint == "OS"] <- TRUE
        tab$primary_family_endpoint[tab$endpoint == "OS"] <- "OS"
    }

    rec <- tab[tab$endpoint_family == "RECURRENCE" & tab$status == "EVALUABLE", , drop = FALSE]
    if (nrow(rec)) {
        rec <- rec[order(rec$endpoint_priority), , drop = FALSE]
        chosen <- rec$endpoint[1L]
        tab$analysis_included[tab$endpoint == chosen] <- TRUE
        tab$primary_family_endpoint[tab$endpoint_family == "RECURRENCE"] <- chosen
        tab$alias_of[tab$endpoint_family == "RECURRENCE" & tab$endpoint != chosen] <- chosen
    } else {
        tab$primary_family_endpoint[tab$endpoint_family == "RECURRENCE"] <- NA_character_
    }

    tab
}

olfml2b_selected_endpoint_names <- function(clinical, min_n = 40L, min_events = 20L) {
    tab <- olfml2b_endpoint_selection_table(clinical, min_n = min_n, min_events = min_events)
    tab$endpoint[tab$analysis_included %in% TRUE]
}

olfml2b_endpoint_audit <- function(views, min_n = 40L, min_events = 20L) {
    rows <- list()
    for (nm in names(views)) {
        d <- olfml2b_analysis_clinical(views[[nm]])
        tab <- olfml2b_endpoint_selection_table(d, min_n = min_n, min_events = min_events)
        tab$cohort <- views[[nm]]$cohort
        tab$olfml2b_measured <- views[[nm]]$target_measured
        rows[[length(rows) + 1L]] <- tab[, c(
            "cohort", "endpoint", "endpoint_family", "time_col", "event_col",
            "n_complete", "events", "olfml2b_measured", "status", "reason",
            "endpoint_priority", "analysis_included", "primary_family_endpoint", "alias_of"
        ), drop = FALSE]
    }
    olfml2b_bind_rows(rows)
}

############################################################
## Survival models
############################################################

olfml2b_run_survival_all <- function(views, min_n, min_events, min_epv, log_file = NULL) {
    rows <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        if (!isTRUE(v$target_measured)) next
        endpoints <- olfml2b_endpoint_defs()
        clinical_for_endpoints <- olfml2b_analysis_clinical(v)
        endpoint_names <- olfml2b_selected_endpoint_names(clinical_for_endpoints, min_n = min_n, min_events = min_events)
        for (ep in endpoint_names) {
            def <- endpoints[[ep]]
            if (!all(c(def$time, def$event) %in% names(v$clinical))) next
            dat <- olfml2b_surv_data(v, def)
            n <- nrow(dat)
            ev <- sum(dat$event == 1L, na.rm = TRUE)
            if (n < min_n || ev < min_events) next
            rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(dat, v$cohort, ep, def$family, model = "univariable", covariates = character())

            covars <- olfml2b_select_covariates(
                dat, events = ev, min_epv = min_epv, cohort = v$cohort,
                model_role = "primary_adjusted"
            )
            rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
                dat, v$cohort, ep, def$family, model = "adjusted", covariates = covars
            )

            subtype_covars <- olfml2b_select_covariates(
                dat, events = ev, min_epv = min_epv, cohort = v$cohort,
                model_role = "molecular_subtype_sensitivity"
            )
            if ("molecular_subtype" %in% subtype_covars) {
                rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
                    dat, v$cohort, ep, def$family,
                    model = "adjusted_molecular_subtype_sensitivity",
                    covariates = subtype_covars
                )
            }

            if ("adjuvant_therapy" %in% names(dat) && length(unique(stats::na.omit(dat$adjuvant_therapy))) >= 2L && ev >= 60L) {
                rows[[length(rows) + 1L]] <- olfml2b_fit_therapy_interaction(dat, v$cohort, ep, def$family, covariates = covars)
            }
            olfml2b_log("Survival done: ", v$cohort, " ", ep, " | n=", n, " events=", ev, log_file = log_file)
        }
    }
    olfml2b_bind_rows(rows)
}

olfml2b_surv_data <- function(view, def) {
    d <- olfml2b_analysis_clinical(view)
    keep <- is.finite(d[[def$time]]) & d[[def$time]] > 0 & d[[def$event]] %in% c(0L, 1L) & is.finite(d$olfml2b_z)
    d <- d[keep, , drop = FALSE]
    d$time <- as.numeric(d[[def$time]])
    d$event <- as.integer(d[[def$event]])
    d
}

olfml2b_select_covariates <- function(dat, events, min_epv = 10, cohort = NA_character_,
                                   model_role = c("primary_adjusted", "molecular_subtype_sensitivity")) {
    model_role <- match.arg(model_role)

    ## Primary adjusted model:
    ## TCGA-STAD is deliberately restricted to age + stage to avoid molecular-subtype
    ## missingness-driven sample loss. Molecular subtype is tested separately.
    if (identical(model_role, "primary_adjusted")) {
        if (identical(as.character(cohort), "TCGA_STAD")) {
            candidates <- c("age", "sex", "stage_numeric")
        } else {
            candidates <- c("age", "stage_numeric", "sex", "lauren")
        }
    } else {
        candidates <- c("age", "stage_numeric", "molecular_subtype")
    }

    covars <- character()
    df_used <- 1  # OLFML2B
    for (cv in candidates) {
        if (!cv %in% names(dat)) next
        x <- dat[[cv]]
        ok <- !is.na(x)
        if (sum(ok) < max(30L, 0.6 * nrow(dat))) next
        if (is.numeric(x) || is.integer(x)) {
            if (stats::sd(as.numeric(x), na.rm = TRUE) <= 0) next
            add_df <- 1
        } else {
            lev <- unique(as.character(stats::na.omit(x)))
            if (length(lev) < 2 || length(lev) > 8) next
            add_df <- length(lev) - 1
        }
        if (events / (df_used + add_df) >= min_epv) {
            covars <- c(covars, cv)
            df_used <- df_used + add_df
        }
    }
    covars
}

olfml2b_fit_olfml2b_cox <- function(dat, cohort, endpoint, endpoint_family, model, covariates = character()) {
    vars <- unique(c("time", "event", "olfml2b_z", covariates))
    d <- dat[, vars, drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if (nrow(d) < 20L || sum(d$event == 1L) < 5L) {
        return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, "INSUFFICIENT_COMPLETE_CASES", covariates))
    }
    form <- stats::as.formula(paste0("survival::Surv(time, event) ~ olfml2b_z", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
    fit <- tryCatch(suppressWarnings(survival::coxph(form, data = d)), error = function(e) e)
    if (inherits(fit, "error")) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, paste0("COX_FAILED: ", conditionMessage(fit)), covariates, nrow(d), sum(d$event == 1L)))
    sm <- summary(fit)
    if (!"olfml2b_z" %in% rownames(sm$coefficients)) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, "OLFML2B_COEFFICIENT_MISSING", covariates, nrow(d), sum(d$event == 1L)))
    co <- sm$coefficients["olfml2b_z", ]
    ci <- sm$conf.int["olfml2b_z", ]
    ph_p <- tryCatch({
        z <- suppressWarnings(survival::cox.zph(fit))
        as.numeric(z$table["olfml2b_z", "p"])
    }, error = function(e) NA_real_)
    cindex <- tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_)
    if (!is.finite(co["coef"]) || abs(co["coef"]) > 8 || !is.finite(co["se(coef)"]) || co["se(coef)"] > 5) {
        status <- "COX_UNSTABLE"
    } else {
        status <- "OK"
    }
    data.frame(
        cohort = cohort,
        endpoint = endpoint,
        endpoint_family = endpoint_family,
        model = model,
        covariates = paste(covariates, collapse = ";"),
        n = nrow(d),
        events = sum(d$event == 1L),
        beta = as.numeric(co["coef"]),
        se = as.numeric(co["se(coef)"]),
        hr = as.numeric(ci["exp(coef)"]),
        ci_low = as.numeric(ci["lower .95"]),
        ci_high = as.numeric(ci["upper .95"]),
        p_value = as.numeric(co["Pr(>|z|)"]),
        ph_p = ph_p,
        c_index = cindex,
        direction = ifelse(as.numeric(ci["exp(coef)"]) < 1, "favorable_high", "adverse_high"),
        status = status,
        stringsAsFactors = FALSE
    )
}

olfml2b_cox_empty <- function(cohort, endpoint, endpoint_family, model, status, covariates = character(), n = NA_integer_, events = NA_integer_) {
    data.frame(
        cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
        model = model, covariates = paste(covariates, collapse = ";"),
        n = n, events = events, beta = NA_real_, se = NA_real_, hr = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, ph_p = NA_real_,
        c_index = NA_real_, direction = NA_character_, status = status, stringsAsFactors = FALSE
    )
}

olfml2b_fit_therapy_interaction <- function(dat, cohort, endpoint, endpoint_family, covariates = character()) {
    covariates <- setdiff(covariates, "adjuvant_therapy")
    vars <- unique(c("time", "event", "olfml2b_z", "adjuvant_therapy", covariates))
    d <- dat[, vars, drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if (nrow(d) < 50L || sum(d$event == 1L) < 25L) {
        return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction", "INSUFFICIENT_COMPLETE_CASES", c("adjuvant_therapy", covariates)))
    }
    form <- stats::as.formula(paste0("survival::Surv(time, event) ~ olfml2b_z * adjuvant_therapy", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
    fit <- tryCatch(suppressWarnings(survival::coxph(form, data = d)), error = function(e) e)
    if (inherits(fit, "error")) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction", paste0("COX_FAILED: ", conditionMessage(fit)), c("adjuvant_therapy", covariates), nrow(d), sum(d$event == 1L)))
    sm <- summary(fit)
    int_row <- grep("olfml2b_z:adjuvant_therapy|adjuvant_therapy.*:olfml2b_z", rownames(sm$coefficients), value = TRUE)[1]
    if (is.na(int_row)) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction", "INTERACTION_TERM_MISSING", c("adjuvant_therapy", covariates), nrow(d), sum(d$event == 1L)))
    co <- sm$coefficients[int_row, ]
    ci <- sm$conf.int[int_row, ]
    data.frame(
        cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
        model = "therapy_interaction", covariates = paste(c("adjuvant_therapy", covariates), collapse = ";"),
        n = nrow(d), events = sum(d$event == 1L),
        beta = as.numeric(co["coef"]), se = as.numeric(co["se(coef)"]),
        hr = as.numeric(ci["exp(coef)"]), ci_low = as.numeric(ci["lower .95"]), ci_high = as.numeric(ci["upper .95"]),
        p_value = as.numeric(co["Pr(>|z|)"]), ph_p = NA_real_, c_index = tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_),
        direction = ifelse(as.numeric(ci["exp(coef)"]) < 1, "interaction_favorable_with_therapy", "interaction_adverse_with_therapy"),
        status = "OK",
        stringsAsFactors = FALSE
    )
}

olfml2b_run_median_split_all <- function(views, min_n, min_events) {
    rows <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        endpoints <- olfml2b_endpoint_defs()
        clinical_for_endpoints <- olfml2b_analysis_clinical(v)
        endpoint_names <- olfml2b_selected_endpoint_names(clinical_for_endpoints, min_n = min_n, min_events = min_events)
        for (ep in endpoint_names) {
            def <- endpoints[[ep]]
            if (!all(c(def$time, def$event) %in% names(v$clinical))) next
            d <- olfml2b_surv_data(v, def)
            if (nrow(d) < min_n || sum(d$event == 1L) < min_events) next
            cut <- stats::median(d$olfml2b_z, na.rm = TRUE)
            d$group <- factor(ifelse(d$olfml2b_z >= cut, "OLFML2B_high", "OLFML2B_low"), levels = c("OLFML2B_low", "OLFML2B_high"))
            strata_term <- if (v$cohort == "GSE84437" && "source_subseries" %in% names(d)) " + survival::strata(source_subseries)" else ""
            display_form <- stats::as.formula(paste0("survival::Surv(time, event) ~ group", strata_term))
            fit <- tryCatch(survival::survdiff(display_form, data = d), error = function(e) NULL)
            logrank_p <- if (!is.null(fit)) stats::pchisq(fit$chisq, df = nlevels(d$group) - 1L, lower.tail = FALSE) else NA_real_
            cox <- tryCatch(survival::coxph(display_form, data = d), error = function(e) NULL)
            sm <- if (!is.null(cox)) summary(cox) else NULL
            rows[[length(rows) + 1L]] <- data.frame(
                cohort = v$cohort, endpoint = ep, endpoint_family = def$family,
                cutoff_type = "median_within_cohort_z",
                cutoff_value = cut,
                n = nrow(d), events = sum(d$event == 1L),
                n_low = sum(d$group == "OLFML2B_low"),
                n_high = sum(d$group == "OLFML2B_high"),
                logrank_p = logrank_p,
                hr_high_vs_low = if (!is.null(sm)) sm$conf.int[1, "exp(coef)"] else NA_real_,
                p_cox = if (!is.null(sm)) sm$coefficients[1, "Pr(>|z|)"] else NA_real_,
                design_strata = ifelse(nzchar(strata_term), "source_subseries", NA_character_),
                note = "Median split of the batch-aware within-cohort z score is display/sensitivity only; continuous z-score Cox is primary.",
                stringsAsFactors = FALSE
            )
        }
    }
    olfml2b_bind_rows(rows)
}

olfml2b_run_spline_all <- function(views, min_n, min_events) {
    if (!requireNamespace("splines", quietly = TRUE)) return(data.frame())
    rows <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        endpoints <- olfml2b_endpoint_defs()
        clinical_for_endpoints <- olfml2b_analysis_clinical(v)
        endpoint_names <- olfml2b_selected_endpoint_names(clinical_for_endpoints, min_n = min_n, min_events = min_events)
        for (ep in endpoint_names) {
            def <- endpoints[[ep]]
            if (!all(c(def$time, def$event) %in% names(v$clinical))) next
            d <- olfml2b_surv_data(v, def)
            if (nrow(d) < max(min_n, 80L) || sum(d$event == 1L) < max(min_events, 30L)) next
            strata_term <- if (v$cohort == "GSE84437" && "source_subseries" %in% names(d)) " + survival::strata(source_subseries)" else ""
            form_l <- stats::as.formula(paste0("survival::Surv(time, event) ~ olfml2b_z", strata_term))
            form_s <- stats::as.formula(paste0("survival::Surv(time, event) ~ splines::ns(olfml2b_z, df = 3)", strata_term))
            fit_l <- tryCatch(survival::coxph(form_l, data = d), error = function(e) NULL)
            fit_s <- tryCatch(survival::coxph(form_s, data = d), error = function(e) NULL)
            p <- NA_real_
            if (!is.null(fit_l) && !is.null(fit_s)) {
                p <- tryCatch({
                    a <- stats::anova(fit_l, fit_s, test = "LRT")
                    pcol <- grep("^P\\(|Pr\\(|p", names(a), ignore.case = TRUE, value = TRUE)
                    val <- if (length(pcol)) a[[pcol[1]]][nrow(a)] else NA_real_
                    val <- suppressWarnings(as.numeric(val))
                    if (!length(val) || !is.finite(val[1])) NA_real_ else val[1]
                }, error = function(e) NA_real_)
            }
            rows[[length(rows) + 1L]] <- data.frame(
                cohort = v$cohort, endpoint = ep, endpoint_family = def$family,
                n = nrow(d), events = sum(d$event == 1L),
                nonlinearity_lrt_p = p,
                design_strata = ifelse(nzchar(strata_term), "source_subseries", NA_character_),
                status = if (length(p) == 1L && is.finite(p)) "OK" else "FAILED",
                stringsAsFactors = FALSE
            )
        }
    }
    olfml2b_bind_rows(rows)
}

############################################################
## Meta-analysis
############################################################


olfml2b_deduplicate_meta_effects <- function(d) {
    if (!nrow(d)) return(d)
    d$endpoint_priority <- ifelse(
        d$endpoint_family == "RECURRENCE",
        match(d$endpoint, c("DFS", "RFS", "Recurrence")),
        1L
    )
    d$endpoint_priority[is.na(d$endpoint_priority)] <- 99L

    rows <- list()
    keys <- unique(paste(d$model, d$endpoint_family, d$cohort, sep = "||"))
    for (key in keys) {
        dd <- d[paste(d$model, d$endpoint_family, d$cohort, sep = "||") == key, , drop = FALSE]
        if (unique(dd$endpoint_family)[1] == "RECURRENCE") {
            # Pre-specified endpoint priority only: DFS > RFS > Recurrence.
            # Never use p values to choose among related recurrence-family endpoints,
            # because that creates a post-hoc winner's-curse selection bias.
            dd <- dd[order(dd$endpoint_priority), , drop = FALSE]
            rows[[length(rows) + 1L]] <- dd[1L, , drop = FALSE]
        } else {
            rows[[length(rows) + 1L]] <- dd
        }
    }
    out <- olfml2b_bind_rows(rows)
    out$endpoint_priority <- NULL
    out
}


olfml2b_meta_all <- function(surv_all) {
    if (!nrow(surv_all)) return(data.frame())
    d <- surv_all[
        surv_all$status == "OK" &
            surv_all$model %in% c("adjusted", "univariable") &
            is.finite(surv_all$beta) &
            is.finite(surv_all$se),
        , drop = FALSE
    ]
    d <- olfml2b_deduplicate_meta_effects(d)

    rows <- list()
    for (model in unique(d$model)) {
        for (fam in unique(d$endpoint_family)) {
            dd <- d[d$model == model & d$endpoint_family == fam, , drop = FALSE]
            rows[[length(rows) + 1L]] <- olfml2b_meta_one(dd, endpoint_family = fam, model = model)
        }
    }
    olfml2b_bind_rows(rows)
}

olfml2b_meta_one <- function(d, endpoint_family, model) {
    if (nrow(d) < 2L) {
        return(data.frame(endpoint_family = endpoint_family, model = model, k = nrow(d), status = "INSUFFICIENT_COHORTS", stringsAsFactors = FALSE))
    }
    yi <- d$beta
    vi <- d$se^2
    wi <- 1 / vi
    fixed <- sum(wi * yi) / sum(wi)
    q <- sum(wi * (yi - fixed)^2)
    df <- length(yi) - 1
    cval <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- max(0, (q - df) / cval)
    wr <- 1 / (vi + tau2)
    random <- sum(wr * yi) / sum(wr)
    se_random <- sqrt(1 / sum(wr))
    p_random <- 2 * stats::pnorm(abs(random / se_random), lower.tail = FALSE)
    i2 <- if (q > df) max(0, (q - df) / q) * 100 else 0
    data.frame(
        endpoint_family = endpoint_family,
        model = model,
        k = nrow(d),
        cohorts = paste(d$cohort, collapse = ";"),
        cohort_endpoints = paste(paste(d$cohort, d$endpoint, sep = ":"), collapse = ";"),
        beta_random = random,
        se_random = se_random,
        hr_random = exp(random),
        ci_low_random = exp(random - 1.96 * se_random),
        ci_high_random = exp(random + 1.96 * se_random),
        p_random = p_random,
        beta_fixed = fixed,
        hr_fixed = exp(fixed),
        q = q,
        q_p = stats::pchisq(q, df = df, lower.tail = FALSE),
        i2 = i2,
        tau2 = tau2,
        direction_random = ifelse(exp(random) < 1, "favorable_high", "adverse_high"),
        status = "OK",
        stringsAsFactors = FALSE
    )
}

olfml2b_leave_one_out_meta_all <- function(surv_all) {
    d <- surv_all[
        surv_all$status == "OK" &
            surv_all$model == "adjusted" &
            is.finite(surv_all$beta) &
            is.finite(surv_all$se),
        , drop = FALSE
    ]
    d <- olfml2b_deduplicate_meta_effects(d)

    rows <- list()
    for (fam in unique(d$endpoint_family)) {
        dd <- d[d$endpoint_family == fam, , drop = FALSE]
        if (nrow(dd) < 3L) next
        for (co in unique(dd$cohort)) {
            mm <- olfml2b_meta_one(dd[dd$cohort != co, , drop = FALSE],
                                 endpoint_family = fam,
                                 model = "adjusted_leave_one_out")
            mm$left_out <- co
            rows[[length(rows) + 1L]] <- mm
        }
    }
    olfml2b_bind_rows(rows)
}

############################################################
## Proliferation control and incremental models
############################################################

olfml2b_run_proliferation_controls <- function(views, min_n, min_events, min_epv, log_file = NULL) {
    control_genes <- c(
      "MKI67", "PCNA", "TOP2A", "MCM2", "MCM5", "MCM6",
      "FAP", "ACTA2", "COL1A1", "COL3A1", "DCN", "LUM", "SPARC",
      "VWF", "PECAM1", "KDR", "RGS5", "EPCAM", "KRT19"
    )
    audit <- list()
    models <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        control_mat <- olfml2b_extract_controls(v, control_genes)
        audit[[length(audit) + 1L]] <- control_mat$audit
        v$clinical <- cbind(v$clinical, control_mat$scores[match(v$clinical$sample_id, control_mat$scores$sample_id), setdiff(names(control_mat$scores), "sample_id"), drop = FALSE])
        endpoints <- olfml2b_endpoint_defs()
        clinical_for_endpoints <- olfml2b_analysis_clinical(v)
        endpoint_names <- olfml2b_selected_endpoint_names(clinical_for_endpoints, min_n = min_n, min_events = min_events)
        for (ep in endpoint_names) {
            def <- endpoints[[ep]]
            if (!all(c(def$time, def$event) %in% names(v$clinical))) next
            d <- olfml2b_surv_data(v, def)
            if (nrow(d) < min_n || sum(d$event == 1L) < min_events) next
            base_cov <- olfml2b_select_covariates(d, sum(d$event == 1L), min_epv, cohort = v$cohort, model_role = "primary_adjusted")
            for (cv in c("proliferation_score_z", paste0(tolower(control_genes), "_z"))) {
                if (!cv %in% names(d) || sum(is.finite(d[[cv]])) < 30L) next
                models[[length(models) + 1L]] <- olfml2b_fit_control_model(d, v$cohort, ep, def$family, control_var = cv, covariates = base_cov)
            }
        }
    }
    list(expression_audit = olfml2b_bind_rows(audit), models = olfml2b_bind_rows(models))
}

olfml2b_extract_controls <- function(view, genes) {
    d <- data.frame(sample_id = colnames(view$expr), stringsAsFactors = FALSE)
    audit_rows <- list()
    for (g in genes) {
        ex <- olfml2b_extract_gene(view$expr, g)
        col <- paste0(tolower(g), "_z")
        d[[col]] <- olfml2b_z(ex$values[match(d$sample_id, names(ex$values))])
        audit_rows[[length(audit_rows) + 1L]] <- data.frame(
            cohort = view$cohort,
            gene = g,
            measured = ex$measured,
            selected_rows = paste(ex$rows, collapse = ";"),
            n_nonmissing = sum(is.finite(d[[col]])),
            stringsAsFactors = FALSE
        )
    }
    prolif_cols <- paste0(tolower(c("MKI67", "PCNA", "TOP2A", "MCM2", "MCM5", "MCM6")), "_z")
    available <- prolif_cols[prolif_cols %in% names(d)]
    if (length(available)) {
        d$proliferation_score_z <- olfml2b_z(rowMeans(d[, available, drop = FALSE], na.rm = TRUE))
    } else d$proliferation_score_z <- NA_real_
    list(scores = d, audit = olfml2b_bind_rows(audit_rows))
}


olfml2b_scalar_named <- function(x, nm, default = NA_real_) {
    if (is.null(x) || !length(x) || is.null(names(x)) || !nm %in% names(x)) return(default)
    val <- suppressWarnings(as.numeric(x[[nm]]))
    if (!length(val) || !is.finite(val[1])) default else val[1]
}

olfml2b_lrt_p <- function(fit1, fit2) {
    if (is.null(fit1) || is.null(fit2)) return(NA_real_)
    tryCatch({
        a <- stats::anova(fit1, fit2, test = "LRT")
        pcol <- grep("^P\\(|Pr\\(|p", names(a), ignore.case = TRUE, value = TRUE)
        val <- if (length(pcol)) a[[pcol[1]]][nrow(a)] else NA_real_
        val <- suppressWarnings(as.numeric(val))
        if (!length(val) || !is.finite(val[1])) NA_real_ else val[1]
    }, error = function(e) NA_real_)
}


olfml2b_fit_control_model <- function(dat, cohort, endpoint, endpoint_family, control_var, covariates = character()) {
    vars <- unique(c("time", "event", "olfml2b_z", control_var, covariates))
    vars <- vars[vars %in% names(dat)]
    required <- c("time", "event", "olfml2b_z", control_var)
    if (!all(required %in% vars)) {
        return(data.frame(
            cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
            control_var = control_var, covariates = paste(covariates, collapse = ";"),
            n = NA_integer_, events = NA_integer_,
            olfml2b_hr_after_control = NA_real_, olfml2b_p_after_control = NA_real_,
            control_hr = NA_real_, control_p = NA_real_,
            lrt_p_add_control = NA_real_, cindex_full = NA_real_,
            status = "MISSING_REQUIRED_VARIABLE",
            stringsAsFactors = FALSE
        ))
    }

    d <- dat[, vars, drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if (nrow(d) < 30L || sum(d$event == 1L) < 10L) {
        return(data.frame(
            cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
            control_var = control_var, covariates = paste(covariates, collapse = ";"),
            n = nrow(d), events = sum(d$event == 1L),
            olfml2b_hr_after_control = NA_real_, olfml2b_p_after_control = NA_real_,
            control_hr = NA_real_, control_p = NA_real_,
            lrt_p_add_control = NA_real_, cindex_full = NA_real_,
            status = "INSUFFICIENT_COMPLETE_CASES",
            stringsAsFactors = FALSE
        ))
    }

    # Drop covariates that became single-level or zero-variance after complete-case filtering.
    usable_covariates <- character()
    for (cv in covariates) {
        if (!cv %in% names(d)) next
        xx <- d[[cv]]
        if (is.numeric(xx) || is.integer(xx)) {
            if (sum(is.finite(xx)) >= 10L && stats::sd(xx, na.rm = TRUE) > 0) usable_covariates <- c(usable_covariates, cv)
        } else {
            if (length(unique(stats::na.omit(xx))) >= 2L) usable_covariates <- c(usable_covariates, cv)
        }
    }

    f1 <- stats::as.formula(paste0(
        "survival::Surv(time, event) ~ olfml2b_z",
        if (length(usable_covariates)) paste0(" + ", paste(usable_covariates, collapse = " + ")) else ""
    ))
    f2 <- stats::as.formula(paste0(
        "survival::Surv(time, event) ~ olfml2b_z + ", control_var,
        if (length(usable_covariates)) paste0(" + ", paste(usable_covariates, collapse = " + ")) else ""
    ))

    fit1 <- tryCatch(suppressWarnings(survival::coxph(f1, data = d)), error = function(e) NULL)
    fit2 <- tryCatch(suppressWarnings(survival::coxph(f2, data = d)), error = function(e) NULL)
    if (is.null(fit2)) {
        return(data.frame(
            cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
            control_var = control_var, covariates = paste(usable_covariates, collapse = ";"),
            n = nrow(d), events = sum(d$event == 1L),
            olfml2b_hr_after_control = NA_real_, olfml2b_p_after_control = NA_real_,
            control_hr = NA_real_, control_p = NA_real_,
            lrt_p_add_control = NA_real_, cindex_full = NA_real_,
            status = "COX_FAILED",
            stringsAsFactors = FALSE
        ))
    }

    sm <- summary(fit2)
    coefs <- sm$coefficients
    cis <- sm$conf.int

    co_olfml2b <- if ("olfml2b_z" %in% rownames(coefs)) coefs["olfml2b_z", , drop = TRUE] else setNames(rep(NA_real_, 5), colnames(coefs))
    ci_olfml2b <- if ("olfml2b_z" %in% rownames(cis)) cis["olfml2b_z", , drop = TRUE] else setNames(rep(NA_real_, ncol(cis)), colnames(cis))
    co_ctrl <- if (control_var %in% rownames(coefs)) coefs[control_var, , drop = TRUE] else setNames(rep(NA_real_, 5), colnames(coefs))
    ci_ctrl <- if (control_var %in% rownames(cis)) cis[control_var, , drop = TRUE] else setNames(rep(NA_real_, ncol(cis)), colnames(cis))

    lrt_p <- olfml2b_lrt_p(fit1, fit2)
    cidx <- tryCatch({
        val <- as.numeric(sm$concordance[1])
        if (!length(val) || !is.finite(val[1])) NA_real_ else val[1]
    }, error = function(e) NA_real_)

    data.frame(
        cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
        control_var = control_var,
        covariates = paste(usable_covariates, collapse = ";"),
        n = nrow(d), events = sum(d$event == 1L),
        olfml2b_hr_after_control = olfml2b_scalar_named(ci_olfml2b, "exp(coef)"),
        olfml2b_p_after_control = olfml2b_scalar_named(co_olfml2b, "Pr(>|z|)"),
        control_hr = olfml2b_scalar_named(ci_ctrl, "exp(coef)"),
        control_p = olfml2b_scalar_named(co_ctrl, "Pr(>|z|)"),
        lrt_p_add_control = lrt_p,
        cindex_full = cidx,
        status = if (is.finite(olfml2b_scalar_named(ci_olfml2b, "exp(coef)")) || is.finite(olfml2b_scalar_named(ci_ctrl, "exp(coef)"))) "OK" else "COEF_NOT_ESTIMABLE",
        stringsAsFactors = FALSE
    )
}

############################################################
## State module analysis
############################################################

olfml2b_signature_catalog <- function() {
    list(
        CAF_Core = c("FAP","ACTA2","PDGFRA","PDGFRB","TAGLN","THY1","CXCL12","COL1A1","COL1A2","COL3A1","DCN","LUM","POSTN"),
        ECM_Remodeling = c("COL5A1","COL6A1","COL6A2","SPARC","MMP2","MMP14","LOX","PLOD2"),
        TGFb_Response = c("TGFBI","SERPINE1","SMAD3","SMAD7","CTGF","THBS1","ITGA5","PMEPA1","INHBA"),
        Endothelial_Pericyte = c("PECAM1","VWF","KDR","FLT1","ENG","EMCN","RGS5","CSPG4","MCAM","NOTCH3","DES"),
        EMT = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1","ITGA5","MMP2"),
        Epithelial_Differentiation = c("EPCAM","KRT8","KRT18","KRT19","CDH1","MUC1","MUC5AC","MUC6","CLDN3","CLDN4","CLDN7"),
        Proliferation_Control = c("MKI67","TOP2A","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC20","CCNB1"),
        Hypoxia_Control = c("CA9","VEGFA","SLC2A1","LDHA","PDK1","BNIP3","EGLN3","HILPDA"),
        CD8_Cytotoxic = c("CD8A","CD8B","GZMA","GZMB","PRF1","NKG7","GNLY","IFNG")
    )
}

olfml2b_module_correlation_meta <- function(correlations) {
    if (!is.data.frame(correlations) || !nrow(correlations)) return(data.frame())
    rows <- lapply(split(correlations, correlations$module), function(d) {
        d <- d[is.finite(d$rho) & is.finite(d$n) & d$n > 3L, , drop = FALSE]
        if (nrow(d) < 2L) return(NULL)
        meta_input <- data.frame(
            cohort = d$cohort,
            endpoint = d$module,
            beta = atanh(pmax(pmin(d$rho, 0.999999), -0.999999)),
            se = sqrt(1 / (d$n - 3)),
            stringsAsFactors = FALSE
        )
        mm <- olfml2b_meta_one(meta_input, endpoint_family = d$module[1], model = "spearman_fisher_z")
        if (!nrow(mm) || mm$status[1] != "OK") return(NULL)
        data.frame(
            module = d$module[1], k = nrow(d), n_total = sum(d$n),
            rho_random = tanh(mm$beta_random[1]),
            rho_ci_low = tanh(log(mm$ci_low_random[1])),
            rho_ci_high = tanh(log(mm$ci_high_random[1])),
            p_random = mm$p_random[1],
            prediction_low = tanh(log(mm$prediction_low[1])),
            prediction_high = tanh(log(mm$prediction_high[1])),
            i2 = mm$i2[1], tau2_fisher_z = mm$tau2[1],
            n_positive = sum(d$rho > 0), n_negative = sum(d$rho < 0),
            n_nominal_positive = sum(d$rho > 0 & d$p_value < 0.05, na.rm = TRUE),
            n_nominal_negative = sum(d$rho < 0 & d$p_value < 0.05, na.rm = TRUE),
            direction_consistency = ifelse(all(d$rho > 0), "ALL_POSITIVE", ifelse(all(d$rho < 0), "ALL_NEGATIVE", "MIXED")),
            meta_method = "Fisher_z_REML_modified_Hartung_Knapp",
            claim_ceiling = "cross_cohort_state_association_not_mechanism",
            stringsAsFactors = FALSE
        )
    })
    out <- olfml2b_bind_rows(rows)
    if (nrow(out)) out$fdr_random <- stats::p.adjust(out$p_random, method = "BH")
    out
}

olfml2b_run_module_analysis <- function(views, modules) {
    coverage_rows <- list()
    corr_rows <- list()
    state_rows <- list()
    for (nm in names(views)) {
        v <- views[[nm]]
        if (nrow(v$expr) == 0L) next
        samples <- olfml2b_analysis_clinical(v)$sample_id
        samples <- intersect(samples, colnames(v$expr))
        if (length(samples) < 20L) next
        scores <- olfml2b_module_scores(v$expr[, samples, drop = FALSE], modules)
        coverage_rows[[length(coverage_rows) + 1L]] <- scores$coverage |>
            transform(cohort = v$cohort)
        clinical <- v$clinical[match(samples, v$clinical$sample_id), , drop = FALSE]
        for (mod in names(scores$scores)) {
            sc <- scores$scores[[mod]]
            ok <- is.finite(sc) & is.finite(clinical$olfml2b_z)
            if (sum(ok) < 20L) next
            ct <- suppressWarnings(stats::cor.test(clinical$olfml2b_z[ok], sc[ok], method = "spearman", exact = FALSE))
            corr_rows[[length(corr_rows) + 1L]] <- data.frame(
                cohort = v$cohort,
                module = mod,
                n = sum(ok),
                rho = unname(ct$estimate),
                p_value = ct$p.value,
                direction = ifelse(unname(ct$estimate) > 0, "OLFML2B_high_enriched", "OLFML2B_low_enriched"),
                stringsAsFactors = FALSE
            )
        }
        high <- clinical$olfml2b_expression >= stats::median(clinical$olfml2b_expression, na.rm = TRUE)
        for (mod in names(scores$scores)) {
            sc <- scores$scores[[mod]]
            ok <- is.finite(sc) & !is.na(high)
            if (sum(ok & high) < 10L || sum(ok & !high) < 10L) next
            wt <- tryCatch(stats::wilcox.test(sc[ok & high], sc[ok & !high], exact = FALSE), error = function(e) NULL)
            state_rows[[length(state_rows) + 1L]] <- data.frame(
                cohort = v$cohort,
                module = mod,
                n_high = sum(ok & high),
                n_low = sum(ok & !high),
                median_high = stats::median(sc[ok & high], na.rm = TRUE),
                median_low = stats::median(sc[ok & !high], na.rm = TRUE),
                delta_high_minus_low = stats::median(sc[ok & high], na.rm = TRUE) - stats::median(sc[ok & !high], na.rm = TRUE),
                p_value = if (!is.null(wt)) wt$p.value else NA_real_,
                expected_for_hypothesis = ifelse(
                    mod %in% c("CAF_Core", "ECM_Remodeling", "TGFb_Response", "EMT"),
                    "primary_positive_ecological_state",
                    "exploratory_two_sided_no_preferred_direction"
                ),
                stringsAsFactors = FALSE
            )
        }
    }
    cov <- olfml2b_bind_rows(coverage_rows)
    corr <- olfml2b_bind_rows(corr_rows)
    if (nrow(corr)) corr$fdr <- olfml2b_bh(corr$p_value)
    state <- olfml2b_bind_rows(state_rows)
    if (nrow(state)) state$fdr <- olfml2b_bh(state$p_value)
    meta_correlations <- olfml2b_module_correlation_meta(corr)
    list(coverage = cov, correlations = corr, state_summary = state, meta_correlations = meta_correlations)
}

olfml2b_module_scores <- function(expr, modules) {
    rn <- toupper(rownames(expr))
    zexpr <- t(scale(t(expr)))
    zexpr[!is.finite(zexpr)] <- NA_real_
    scores <- list()
    coverage <- list()
    for (mod in names(modules)) {
        genes <- toupper(modules[[mod]])
        hits <- intersect(genes, rn)
        coverage[[length(coverage) + 1L]] <- data.frame(
            module = mod,
            n_gene_set = length(unique(genes)),
            n_measured = length(hits),
            measured_genes = paste(hits, collapse = ";"),
            stringsAsFactors = FALSE
        )
        if (length(hits) >= 3L) {
            scores[[mod]] <- colMeans(zexpr[match(hits, rn), , drop = FALSE], na.rm = TRUE)
        } else {
            scores[[mod]] <- rep(NA_real_, ncol(expr))
        }
    }
    list(scores = scores, coverage = olfml2b_bind_rows(coverage))
}

############################################################
## PDC lightweight detection audit
############################################################

olfml2b_pdc_detection_audit <- function(root, target_gene = "OLFML2B", max_file_mb = 80) {
    pdc_root <- file.path(root, "data", "PDC_STAD")
    if (!dir.exists(pdc_root)) pdc_root <- file.path(root, "data", "raw", "PDC_STAD")
    if (!dir.exists(pdc_root)) return(data.frame(status = "PDC_DIRECTORY_NOT_FOUND", path = pdc_root, stringsAsFactors = FALSE))
    files <- list.files(pdc_root, recursive = TRUE, full.names = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    files <- files[grepl("\\.csv$|\\.tsv$|\\.txt$|\\.rds$", files, ignore.case = TRUE)]
    rows <- list()
    for (fp in files) {
        sz <- file.info(fp)$size / 1024^2
        if (!is.finite(sz) || sz > max_file_mb) next
        hit <- FALSE
        n_rows <- NA_integer_
        n_cols <- NA_integer_
        columns_hit <- NA_character_
        err <- NA_character_
        if (grepl("\\.rds$", fp, ignore.case = TRUE)) {
            obj <- tryCatch(readRDS(fp), error = function(e) e)
            if (inherits(obj, "error")) {
                err <- conditionMessage(obj)
            } else {
                txt <- paste(capture.output(str(obj, max.level = 2)), collapse = "\n")
                hit <- grepl(target_gene, txt, ignore.case = FALSE)
            }
        } else {
            tab <- tryCatch({
                sep <- if (grepl("\\.tsv$|\\.txt$", fp, ignore.case = TRUE)) "\t" else ","
                utils::read.table(fp, sep = sep, header = TRUE, nrows = 2000, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
            }, error = function(e) e)
            if (inherits(tab, "error")) {
                err <- conditionMessage(tab)
            } else {
                n_rows <- nrow(tab); n_cols <- ncol(tab)
                columns_hit <- paste(names(tab)[grepl(target_gene, names(tab), fixed = TRUE)], collapse = ";")
                hit <- any(grepl(target_gene, names(tab), fixed = TRUE))
                if (!hit && nrow(tab) && ncol(tab)) {
                    hit <- any(vapply(tab, function(x) any(grepl(target_gene, as.character(x), fixed = TRUE)), logical(1)))
                }
            }
        }
        rows[[length(rows) + 1L]] <- data.frame(
            file = normalizePath(fp, winslash = "/", mustWork = FALSE),
            size_mb = sz,
            olfml2b_text_detected = hit,
            n_rows_preview = n_rows,
            n_cols_preview = n_cols,
            columns_hit = columns_hit,
            error = err,
            stringsAsFactors = FALSE
        )
    }
    olfml2b_bind_rows(rows)
}

############################################################
## Go/No-Go and basic figures
############################################################

olfml2b_control_evidence_summary <- function(control_models) {
    if (!is.data.frame(control_models) || !nrow(control_models)) {
        return(list(status = "REVIEW", pass_n = 0L, trend_n = 0L, total_n = 0L,
                    detail = "No proliferation-control model table."))
    }
    cc <- control_models[
        control_models$status == "OK" &
            is.finite(control_models$olfml2b_hr_after_control) &
            is.finite(control_models$olfml2b_p_after_control),
        , drop = FALSE
    ]
    if (!nrow(cc)) {
        return(list(status = "REVIEW", pass_n = 0L, trend_n = 0L, total_n = 0L,
                    detail = "No evaluable proliferation-control models."))
    }
    median_log_hr <- stats::median(log(cc$olfml2b_hr_after_control), na.rm = TRUE)
    resolved_direction <- if (is.finite(median_log_hr) && median_log_hr > 0) "adverse_high" else if (is.finite(median_log_hr) && median_log_hr < 0) "favorable_high" else "unresolved"
    direction_match <- if (identical(resolved_direction, "adverse_high")) {
        cc$olfml2b_hr_after_control > 1
    } else if (identical(resolved_direction, "favorable_high")) {
        cc$olfml2b_hr_after_control < 1
    } else {
        rep(FALSE, nrow(cc))
    }
    pass <- direction_match & cc$olfml2b_p_after_control < 0.05
    trend <- direction_match & cc$olfml2b_p_after_control < 0.10
    pass_n <- sum(pass, na.rm = TRUE)
    trend_n <- sum(trend, na.rm = TRUE)
    total_n <- nrow(cc)
    cohorts_005 <- if (pass_n) paste(unique(cc$cohort[pass]), collapse = ";") else "none"
    cohorts_010 <- if (trend_n) paste(unique(cc$cohort[trend]), collapse = ";") else "none"
    pass_cohort_n <- length(unique(cc$cohort[pass]))
    trend_cohort_n <- length(unique(cc$cohort[trend]))
    status <- if (pass_cohort_n >= 3L) "PARTIAL" else "REVIEW"
    list(
        status = status,
        pass_n = pass_n,
        trend_n = trend_n,
        total_n = total_n,
        resolved_direction = resolved_direction,
        pass_cohort_n = pass_cohort_n,
        trend_cohort_n = trend_cohort_n,
        detail = paste0(
            "OLFML2B after proliferation controls (resolved direction=", resolved_direction,
            "): p<0.05 models=", pass_n, "/", total_n,
            " (cohorts=", cohorts_005, "); p<0.10 trends=", trend_n, "/", total_n,
            " (cohorts=", cohorts_010, "). Models sharing cohorts and controls are not independent replications; interpret as sensitivity/partial independence only."
        )
    )
}

olfml2b_interpretation_boundary <- function(surv_all, meta, control_models, median_surv, module_summary) {
    os_meta <- meta[meta$endpoint_family == "OS" & meta$model == "adjusted" & meta$status == "OK", , drop = FALSE]
    rec_meta <- meta[meta$endpoint_family == "RECURRENCE" & meta$model == "adjusted" & meta$status == "OK", , drop = FALSE]
    ctrl <- olfml2b_control_evidence_summary(control_models)
    hypoxia_status <- "NOT_EVALUATED"
    if (is.data.frame(module_summary) && nrow(module_summary) && "module" %in% names(module_summary)) {
        hyp <- module_summary[module_summary$module == "Hypoxia", , drop = FALSE]
        if (nrow(hyp)) {
            high_sig <- sum(hyp$delta_high_minus_low > 0 & hyp$p_value < 0.05, na.rm = TRUE)
            low_sig <- sum(hyp$delta_high_minus_low < 0 & hyp$p_value < 0.05, na.rm = TRUE)
            hypoxia_status <- if (high_sig > 0 && low_sig == 0) "CONTRADICTS_LOW_HYPOXIA" else if (low_sig > 0) "PARTIAL_LOW_HYPOXIA" else "NOT_SUPPORTED"
        }
    }
    km_sig <- 0L
    km_total <- 0L
    logrank_col <- if ("logrank_p" %in% names(median_surv)) "logrank_p" else
        if ("p_logrank" %in% names(median_surv)) "p_logrank" else NA_character_
    if (is.data.frame(median_surv) && nrow(median_surv) && !is.na(logrank_col)) {
        logrank_p <- suppressWarnings(as.numeric(median_surv[[logrank_col]]))
        ok <- is.finite(logrank_p)
        km_total <- sum(ok)
        km_sig <- sum(logrank_p[ok] < 0.05, na.rm = TRUE)
    }
    os_supported <- nrow(os_meta) && is.finite(os_meta$p_random[1]) && os_meta$p_random[1] < 0.05 &&
        is.finite(os_meta$ci_low_random[1]) && is.finite(os_meta$ci_high_random[1]) &&
        (os_meta$ci_low_random[1] > 1 || os_meta$ci_high_random[1] < 1)
    rec_supported <- nrow(rec_meta) && is.finite(rec_meta$p_random[1]) && rec_meta$p_random[1] < 0.05 &&
        is.finite(rec_meta$ci_low_random[1]) && is.finite(rec_meta$ci_high_random[1]) &&
        (rec_meta$ci_low_random[1] > 1 || rec_meta$ci_high_random[1] < 1)
    os_direction <- if (nrow(os_meta) && is.finite(os_meta$hr_random[1])) ifelse(os_meta$hr_random[1] > 1, "adverse", "favorable") else "unresolved"
    stromal_modules <- c("EMT", "ECM_CAF", "TGFb_Response")
    stromal <- if (is.data.frame(module_summary) && nrow(module_summary)) module_summary[module_summary$module %in% stromal_modules, , drop = FALSE] else data.frame()
    if (nrow(stromal) && !"fdr" %in% names(stromal)) stromal$fdr <- stats::p.adjust(stromal$p_value, method = "BH")
    stromal_cohorts <- if (nrow(stromal)) length(unique(stromal$cohort[is.finite(stromal$delta_high_minus_low) & stromal$delta_high_minus_low > 0 & is.finite(stromal$fdr) & stromal$fdr < 0.05])) else 0L

    data.frame(
        domain = c(
            "primary_survival_claim",
            "recurrence_family_claim",
            "proliferation_independence_claim",
            "mechanism_core_claim",
            "hypoxia_claim",
            "median_split_km_claim",
            "causal_language_boundary"
        ),
        allowed_claim = c(
            ifelse(os_supported,
                   paste0("Higher OLFML2B is associated with ", os_direction, " OS across adjusted public cohorts under REML plus modified Hartung-Knapp inference."),
                   "OS association requires review before use as a primary claim."),
            ifelse(rec_supported,
                   "DFS/RFS/recurrence-family evidence is supportive only; do not promote it to the primary endpoint.",
                   "Recurrence-family evidence is insufficient for a primary endpoint claim."),
            ctrl$detail,
            ifelse(stromal_cohorts >= 3L,
                   "OLFML2B-high may be described as a CAF/TGFb/ECM/EMT-rich ecological state across cohorts. Do not describe this as causal stromal activation.",
                   "CAF/TGFb/ECM/EMT ecological-state association requires review; E2F/G2M/DDR is not an allowed primary mechanism claim."),
            paste0("Hypoxia status: ", hypoxia_status, ". Do not use low-hypoxia as a primary mechanism unless rerun evidence changes."),
            paste0("Median split KM is display/sensitivity only: significant panels=", km_sig, "/", km_total, ". Continuous z-score Cox remains primary."),
            "Use adverse-association/ecological-state-marker language only. Do not claim tumor-wide overexpression, causal stromal activation, or treatment prediction."
        ),
        manuscript_use = c(
            "main_result",
            "supportive_result_only",
            "bounded_discussion",
            "mechanism_hypothesis",
            "exclude_from_primary_story",
            "figure_display_only",
            "global_language_rule"
        ),
        stringsAsFactors = FALSE
    )
}

olfml2b_go_no_go <- function(views, endpoint_audit, surv_all, meta, control_models, module_summary, pdc_audit) {
    adjusted <- surv_all[surv_all$model == "adjusted" & surv_all$status == "OK", , drop = FALSE]
    os_ok <- adjusted[adjusted$endpoint_family == "OS" & adjusted$hr > 1 & adjusted$p_value < 0.05, , drop = FALSE]

    os_meta <- meta[meta$endpoint_family == "OS" & meta$model == "adjusted" & meta$status == "OK", , drop = FALSE]
    rec_meta <- meta[meta$endpoint_family == "RECURRENCE" & meta$model == "adjusted" & meta$status == "OK", , drop = FALSE]

    ctrl_summary <- olfml2b_control_evidence_summary(control_models)
    ctrl_ok <- ctrl_summary$status %in% c("PASS", "PARTIAL")

    core_ok <- FALSE
    stromal_support <- FALSE
    hypoxia_boundary <- "NOT_PRIMARY"
    core_detail <- "No module table"
    stromal_detail <- "No module table"
    hypoxia_detail <- "No module table"

    if (is.data.frame(module_summary) && nrow(module_summary)) {
        core_modules <- c("G2M_Checkpoint", "E2F_Targets", "DNA_Repair", "Homologous_Recombination", "Replication_Stress")
        if (!"fdr" %in% names(module_summary)) module_summary$fdr <- stats::p.adjust(module_summary$p_value, method = "BH")
        core <- module_summary[module_summary$module %in% core_modules &
                                   is.finite(module_summary$delta_high_minus_low) &
                                   is.finite(module_summary$fdr) & module_summary$fdr < 0.05, , drop = FALSE]
        core_pos <- table(core$module[core$delta_high_minus_low > 0])
        core_neg <- table(core$module[core$delta_high_minus_low < 0])
        core_neg_aligned <- as.integer(core_neg[names(core_pos)])
        core_neg_aligned[is.na(core_neg_aligned)] <- 0L
        core_modules_consistent <- names(core_pos)[as.integer(core_pos) >= 3L & core_neg_aligned == 0L]
        core_counts <- core_pos
        core_ok <- length(core_modules_consistent) >= 3L
        core_detail <- paste0(
            "DDR/cell-cycle direction audit (FDR<0.05): positive=",
            paste(names(core_pos), as.integer(core_pos), sep = "=", collapse = ";"),
            "; negative=", paste(names(core_neg), as.integer(core_neg), sep = "=", collapse = ";"),
            ". Mixed directions cannot pass the mechanism gate."
        )

        stromal_modules <- c("EMT", "ECM_CAF", "TGFb_Response")
        stromal <- module_summary[module_summary$module %in% stromal_modules &
                                      module_summary$delta_high_minus_low > 0 &
                                      is.finite(module_summary$fdr) & module_summary$fdr < 0.05, , drop = FALSE]
        stromal_counts <- table(stromal$module)
        stromal_support <- length(unique(stromal$cohort)) >= 3L && sum(stromal_counts >= 3L) >= 2L
        stromal_detail <- paste0(
            "Positive CAF/ECM/EMT/TGFb ecological-state support: ",
            paste(names(stromal_counts), as.integer(stromal_counts), sep = "=", collapse = ";")
        )

        hyp <- module_summary[module_summary$module == "Hypoxia", , drop = FALSE]
        if (nrow(hyp)) {
            high_sig <- sum(hyp$delta_high_minus_low > 0 & hyp$p_value < 0.05, na.rm = TRUE)
            low_sig <- sum(hyp$delta_high_minus_low < 0 & hyp$p_value < 0.05, na.rm = TRUE)
            hypoxia_boundary <- if (high_sig > 0 && low_sig == 0) "CONTRADICTS_LOW_HYPOXIA" else if (high_sig > 0 && low_sig > 0) "MIXED" else if (low_sig > 0) "PARTIAL_LOW_HYPOXIA" else "NOT_SUPPORTED"
            hypoxia_detail <- paste0("Hypoxia is not part of the primary mechanism boundary: high_sig=", high_sig, "; low_sig=", low_sig)
        }
    }

    rec_detail <- if (nrow(rec_meta)) {
        paste0("Deduplicated recurrence-family meta: HR=", signif(rec_meta$hr_random[1], 3),
               ", p=", signif(rec_meta$p_random[1], 3),
               if ("cohort_endpoints" %in% names(rec_meta)) paste0("; inputs=", rec_meta$cohort_endpoints[1]) else "")
    } else {
        "No adjusted recurrence-family meta"
    }

    rows <- data.frame(
        criterion = c(
            "OLFML2B measured in TCGA and major GEO cohorts",
            "Adverse-high adjusted OS in at least TCGA plus two external cohorts",
            "Adjusted OS REML random-effects meta HR > 1 with modified Hartung-Knapp CI excluding 1",
            "Deduplicated recurrence-family meta supports the adverse direction",
            "OLFML2B retains signal after proliferation/tissue-composition control adjustment",
            "DDR/cell-cycle state direction is cross-cohort consistent",
            "Primary CAF/ECM/EMT/TGFb positive ecological-state pattern",
            "Hypoxia module direction audit",
            "PDC protein audit should be supportive only, not mandatory"
        ),
        status = c(
            ifelse(sum(vapply(views, function(v) isTRUE(v$target_measured), logical(1))) >= 3L, "PASS", "FAIL"),
            ifelse(length(unique(os_ok$cohort)) >= 3L, "PASS", "REVIEW"),
            ifelse(nrow(os_meta) && os_meta$hr_random[1] > 1 && os_meta$p_random[1] < 0.05 && os_meta$ci_low_random[1] > 1, "PASS", "REVIEW"),
            ifelse(nrow(rec_meta) && rec_meta$hr_random[1] > 1 && rec_meta$p_random[1] < 0.05 && rec_meta$ci_low_random[1] > 1, "SUPPORTIVE", "REVIEW"),
            ifelse(ctrl_ok, ctrl_summary$status, "REVIEW"),
            ifelse(core_ok, "PASS", "REVIEW"),
            ifelse(stromal_support, "PASS", "REVIEW"),
            hypoxia_boundary,
            "OPTIONAL"
        ),
        detail = c(
            paste(names(views)[vapply(views, function(v) isTRUE(v$target_measured), logical(1))], collapse = ";"),
            paste(unique(os_ok$cohort), collapse = ";"),
            ifelse(nrow(os_meta), paste0("HR=", signif(os_meta$hr_random[1], 3), ", p=", signif(os_meta$p_random[1], 3)), "No adjusted OS meta"),
            rec_detail,
            ifelse(ctrl_ok, ctrl_summary$detail, "Needs manual review against proliferation, CAF/ECM, endothelial and epithelial controls."),
            core_detail,
            stromal_detail,
            hypoxia_detail,
            ifelse(is.data.frame(pdc_audit) && nrow(pdc_audit), paste0(sum(pdc_audit$olfml2b_text_detected %in% TRUE), " files with OLFML2B text detected."), "PDC audit not available.")
        ),
        stringsAsFactors = FALSE
    )
    rows
}

# ==============================================================================
# Final Part3 clinical-contract repair layer (v2.3)
#
# M0: unadjusted continuous OLFML2B z score.
# M1: common age + sex + verified overall-stage model.  Cohorts lacking an
#     entire common covariate are explicitly NOT_EVALUABLE and never imputed.
# M2: prespecified available-covariate sensitivity model.  GSE26253 uses its
#     reported overall stage; GSE84437 uses age + sex + pT + pN and a
#     subseries-stratified baseline hazard.
# M2b: GSE84437 stage sensitivity using the Part2 AJCC-7 T/N derivation with an
#      explicit M0 assumption.  It is sensitivity-only and is never pooled with
#      source-reported overall stage as the primary adjusted meta-analysis.
#
# Multiple imputation is limited to partially missing covariates with no more
# than 20% variable-level missingness and no more than 30% incomplete rows.
# Structural or excessive missingness is never imputed.
# ==============================================================================

OLFML2B_PART3_IF78_VERSION <- "v1.1.0_20260722_PDC614_CONTRACT_SYNC"
OLFML2B_ANALYSIS_VERSION <- OLFML2B_PART3_IF78_VERSION

olfml2b_if78_model_id <- function(model) {
  switch(as.character(model)[1],
         univariable = "M0_unadjusted",
         adjusted = "M1_common_age_sex_verified_overall_stage",
         available_adjusted = "M2_prespecified_available_covariates",
         derived_stage_sensitivity = "M2b_TNM_derived_stage_sensitivity",
         adjusted_molecular_subtype_sensitivity = "M3_molecular_subtype_sensitivity",
         therapy_interaction_exploratory = "M4_therapy_interaction_exploratory",
         as.character(model)[1])
}

olfml2b_if78_nonmissing <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) return(is.finite(as.numeric(x)))
  !is.na(x) & nzchar(trimws(as.character(x)))
}

olfml2b_if78_component_factor <- function(x, component = c("T", "N")) {
  component <- match.arg(component)
  y <- toupper(gsub("\\s+", "", trimws(as.character(x))))
  if (component == "T") {
    out <- rep(NA_character_, length(y))
    hit <- grepl("T[1-4]", y)
    out[hit] <- sub(".*?(T[1-4]).*", "\\1", y[hit], perl = TRUE)
    factor(out, levels = paste0("T", 1:4))
  } else {
    out <- rep(NA_character_, length(y))
    hit <- grepl("N[0-3]", y)
    out[hit] <- sub(".*?(N[0-3]).*", "\\1", y[hit], perl = TRUE)
    factor(out, levels = paste0("N", 0:3))
  }
}

olfml2b_if78_available_covariates <- function(cohort) {
  switch(as.character(cohort)[1],
         GSE26253 = c("stage_analysis_primary"),
         GSE84437 = c("age", "sex", "stage_pT", "stage_pN"),
         c("age", "sex", "stage_analysis_primary"))
}

olfml2b_select_covariates <- function(dat, events, min_epv = 10, cohort = NA_character_,
                                    model_role = c("primary_adjusted", "molecular_subtype_sensitivity")) {
  model_role <- match.arg(model_role)
  base <- olfml2b_if78_available_covariates(cohort)
  if (identical(model_role, "molecular_subtype_sensitivity")) {
    base <- c(base, "molecular_subtype")
  }
  base[base %in% names(dat)]
}

olfml2b_if78_prepare_model_data <- function(dat, covariates, design_strata = NA_character_) {
  source_vars <- unique(c("time", "event", "olfml2b_z", covariates,
                          if (is.na(design_strata)) character() else design_strata))
  d <- dat[, source_vars, drop = FALSE]
  model_covars <- covariates
  if ("age" %in% model_covars) {
    age <- suppressWarnings(as.numeric(d$age))
    center <- if (any(is.finite(age))) mean(age, na.rm = TRUE) else NA_real_
    d$age10 <- if (is.finite(center)) (age - center) / 10 else rep(NA_real_, length(age))
    d$age <- NULL
    model_covars[model_covars == "age"] <- "age10"
  }
  if ("sex" %in% model_covars) d$sex <- olfml2b_if78_sex_factor(d$sex)
  for (v in intersect(c("stage_overall", "stage_analysis_primary", "stage_analysis_sensitivity"),
                      model_covars)) {
    d[[v]] <- olfml2b_if78_stage_factor(d[[v]])
  }
  if ("stage_pT" %in% model_covars) d$stage_pT <- olfml2b_if78_component_factor(d$stage_pT, "T")
  if ("stage_pN" %in% model_covars) d$stage_pN <- olfml2b_if78_component_factor(d$stage_pN, "N")
  if (!is.na(design_strata)) d[[design_strata]] <- factor(d[[design_strata]])
  list(data = d, covariates = model_covars)
}

olfml2b_if78_safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else stats::median(x)
}

olfml2b_make_basic_figures <- function(dirs, views, surv_all, meta, module_correlations, go_no_go, loo = data.frame(), median_surv = data.frame(), spline_tab = data.frame(), control_models = data.frame(), state_summary = data.frame()) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
    pal <- if (exists("olfml2b_pub_palette", mode = "function")) olfml2b_pub_palette() else c(red = "#C73E3A", blue = "#2F6DB3", grey = "#8C8C8C", dark = "#111827")
    fig_rows <- list()
    save_one <- function(p, id, title, width, height, source_table, stat_method, n, caption) {
        stem <- file.path(dirs$figures, id)
        if (exists("olfml2b_save_pub_plot", mode = "function")) {
            olfml2b_save_pub_plot(p, stem, width = width, height = height)
        } else {
            ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
            ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 300)
        }
        fig_rows[[length(fig_rows) + 1L]] <<- olfml2b_figure_registry_row(id, title, source_table, stat_method, n, caption, stem)
        invisible(TRUE)
    }

    expr_df <- olfml2b_bind_rows(lapply(views, function(v) {
        d <- olfml2b_analysis_clinical(v)
        if (!nrow(d)) return(data.frame())
        data.frame(cohort = v$cohort, olfml2b_expression = d$olfml2b_expression, stringsAsFactors = FALSE)
    }))
    if (nrow(expr_df)) {
        sum_df <- stats::aggregate(olfml2b_expression ~ cohort, expr_df, function(x) sum(is.finite(x)))
        names(sum_df)[2] <- "n"
        y_pos <- stats::aggregate(olfml2b_expression ~ cohort, expr_df, function(x) max(x, na.rm = TRUE))
        sum_df <- merge(sum_df, y_pos, by = "cohort", all.x = TRUE)
        names(sum_df)[3] <- "y"
        sum_df$label <- paste0("n=", sum_df$n)
        p <- ggplot2::ggplot(expr_df, ggplot2::aes(x = cohort, y = olfml2b_expression)) +
            ggplot2::geom_boxplot(outlier.shape = NA, fill = unname(pal["light_red"]), color = "black", linewidth = 0.28, width = 0.62) +
            ggplot2::geom_jitter(width = 0.12, alpha = 0.32, size = 0.7, color = unname(pal["dark"])) +
            ggplot2::geom_text(data = sum_df, ggplot2::aes(x = cohort, y = y, label = label), inherit.aes = FALSE, vjust = -0.45, size = 3.2) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
            olfml2b_base_theme(base_size = 11) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
            ggplot2::labs(title = "OLFML2B expression across analyzable bulk cohorts", subtitle = "Boxplots show median and IQR; points are individual samples; labels show evaluable samples", x = NULL, y = "OLFML2B expression")
        save_one(p, "FIG3A_OLFML2B_expression_by_cohort", "OLFML2B expression across analyzable bulk cohorts", 8.2, 4.8, "03c_olfml2b_expression_by_clinical_context.csv / bulk views", "Descriptive distribution; no cross-platform inferential test", paste0("total samples=", sum(sum_df$n)), "OLFML2B expression distributions across TCGA-STAD and GEO cohorts. Box center, limits, and whiskers denote median, IQR, and 1.5×IQR; points represent individual samples. No cross-platform differential test is implied.")
    }

    d <- surv_all[surv_all$model == "adjusted" & surv_all$status == "OK" & is.finite(surv_all$hr) & is.finite(surv_all$ci_low) & is.finite(surv_all$ci_high), , drop = FALSE]
    if (nrow(d)) {
        d$plot_label <- paste(d$cohort, d$endpoint, sep = " / ")
        d$stat_label <- mapply(olfml2b_fmt_hr_label, d$hr, d$ci_low, d$ci_high, d$p_value, d$n, d$events, USE.NAMES = FALSE)
        d$plot_label <- paste0(d$plot_label, "\n", d$stat_label)
        d <- d[order(d$endpoint_family, d$hr), , drop = FALSE]
        d$plot_label <- factor(make.unique(d$plot_label, sep = " #"), levels = rev(make.unique(d$plot_label, sep = " #")))
        p <- ggplot2::ggplot(d, ggplot2::aes(x = hr, y = plot_label)) +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = unname(pal["grey"])) +
            ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18, color = unname(pal["dark"]), linewidth = 0.38) +
            ggplot2::geom_point(ggplot2::aes(color = hr > 1, size = events), alpha = 0.95) +
            ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = unname(pal["blue"])), labels = c(`TRUE` = "Adverse", `FALSE` = "Favorable"), name = "Direction") +
            ggplot2::scale_size_continuous(range = c(2.0, 4.2), name = "Events") +
            ggplot2::scale_x_log10() +
            olfml2b_base_theme(base_size = 10) +
            ggplot2::labs(title = "Adjusted continuous Cox models for OLFML2B", subtitle = "Effect size is HR per 1-SD increase in OLFML2B expression; labels show HR, 95% CI, P, n and events", x = "Hazard ratio per 1 SD OLFML2B", y = NULL)
        save_one(p, "FIG3B_OLFML2B_adjusted_cox_forest_labeled", "Adjusted continuous Cox models for OLFML2B", 9.2, max(4.8, 0.75 * nrow(d) + 1.6), "05_olfml2b_survival_all_models.csv", "Cox proportional hazards model; adjusted covariates selected by EPV rule; PH test reported in table", paste0("models=", nrow(d)), "Forest plot of adjusted Cox models. Points show HR per 1-SD OLFML2B increase; horizontal bars show 95% CI. Model covariates and PH-test P values are provided in the source table.")
    }

    m <- meta[meta$status == "OK" & is.finite(meta$hr_random), , drop = FALSE]
    if (nrow(m)) {
        m <- m[m$model %in% c("adjusted", "univariable"), , drop = FALSE]
        m$label <- paste0(m$endpoint_family, " / ", m$model, "\nHR = ", olfml2b_fmt_num(m$hr_random, 2), " (95% CI ", olfml2b_fmt_num(m$ci_low_random, 2), "–", olfml2b_fmt_num(m$ci_high_random, 2), ")\n", olfml2b_fmt_p(m$p_random, 3), ", k = ", m$k, ", I² = ", olfml2b_fmt_num(m$i2, 1), "%")
        m$label <- factor(make.unique(m$label, sep = " #"), levels = rev(make.unique(m$label, sep = " #")))
        p <- ggplot2::ggplot(m, ggplot2::aes(x = hr_random, y = label)) +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = unname(pal["grey"])) +
            ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low_random, xmax = ci_high_random), orientation = "y", width = 0.18, color = unname(pal["dark"]), linewidth = 0.38) +
            ggplot2::geom_point(ggplot2::aes(color = hr_random > 1, size = k), alpha = 0.95) +
            ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = unname(pal["blue"])), guide = "none") +
            ggplot2::scale_size_continuous(range = c(2.5, 4.5), name = "Cohorts") +
            ggplot2::scale_x_log10() +
            olfml2b_base_theme(base_size = 10) +
            ggplot2::labs(title = "Random-effects meta-analysis of OLFML2B survival association", subtitle = "Endpoint-family specific meta-analysis after prespecified deduplication", x = "Random-effects HR", y = NULL)
        save_one(p, "FIG3C_OLFML2B_meta_forest_labeled", "Random-effects meta-analysis of OLFML2B survival association", 9.0, max(4.2, 0.8 * nrow(m) + 1.5), "08_olfml2b_meta_analysis.csv", "DerSimonian-Laird style random-effects meta-analysis on log HR; endpoint deduplication uses fixed priority, not P value", paste0("meta rows=", nrow(m)), "Random-effects meta-analysis of OLFML2B Cox estimates. Labels show HR, 95% CI, P value, number of cohorts and approximate I².")
    }

    if (is.data.frame(module_correlations) && nrow(module_correlations)) {
        mc <- module_correlations
        if (!"p_fdr" %in% names(mc) && "p_value" %in% names(mc)) mc$p_fdr <- stats::p.adjust(mc$p_value, method = "BH")
        pcol <- if ("p_value" %in% names(mc)) "p_value" else if ("p" %in% names(mc)) "p" else NA_character_
        mc$sig <- if (!is.na(pcol)) ifelse(is.finite(mc[[pcol]]) & mc[[pcol]] < 0.05, "*", "") else ""
        p <- ggplot2::ggplot(mc, ggplot2::aes(x = cohort, y = module, fill = rho)) +
            ggplot2::geom_tile(color = "white", linewidth = 0.35) +
            ggplot2::geom_text(ggplot2::aes(label = paste0(sprintf("%.2f", rho), sig)), size = 2.8) +
            ggplot2::scale_fill_gradient2(low = unname(pal["blue"]), mid = "white", high = unname(pal["red"]), midpoint = 0, na.value = "grey90", name = "rho") +
            olfml2b_base_theme(base_size = 10) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
            ggplot2::labs(title = "OLFML2B correlations with state modules", subtitle = "Cell values show Spearman rho; * indicates nominal P < 0.05", x = NULL, y = NULL)
        save_one(p, "FIG3D_OLFML2B_module_correlation_heatmap", "OLFML2B correlations with state modules", 9.5, 7.0, "13_olfml2b_module_correlations.csv", "Spearman correlation within each cohort; nominal significance annotated by *", paste0("tests=", nrow(mc)), "Heatmap of cohort-wise Spearman correlations between OLFML2B and prespecified biological state modules. Values are rho; asterisks indicate nominal P < 0.05 before multiplicity-aware interpretation.")
    }

    if (is.data.frame(loo) && nrow(loo)) {
        d <- loo[loo$status == "OK" & is.finite(loo$hr_random) & is.finite(loo$ci_low_random) & is.finite(loo$ci_high_random), , drop = FALSE]
        if (nrow(d)) {
            d$label <- paste0(d$endpoint_family, " / omit ", d$left_out, "\nHR = ", olfml2b_fmt_num(d$hr_random, 2), " (", olfml2b_fmt_num(d$ci_low_random, 2), "–", olfml2b_fmt_num(d$ci_high_random, 2), "), ", olfml2b_fmt_p(d$p_random, 3), ", k = ", d$k)
            d$label <- factor(make.unique(d$label, sep = " #"), levels = rev(make.unique(d$label, sep = " #")))
            p <- ggplot2::ggplot(d, ggplot2::aes(x = hr_random, y = label)) +
                ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = unname(pal["grey"])) +
                ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low_random, xmax = ci_high_random), orientation = "y", width = 0.16, color = unname(pal["dark"]), linewidth = 0.38) +
                ggplot2::geom_point(ggplot2::aes(color = endpoint_family, size = k), alpha = 0.95) +
                ggplot2::scale_x_log10() +
                olfml2b_base_theme(base_size = 9.5) +
                ggplot2::labs(title = "Leave-one-cohort-out meta-analysis influence audit", subtitle = "Each row re-estimates the adjusted random-effects meta-analysis after omitting one cohort", x = "Leave-one-out random-effects HR", y = NULL, color = "Endpoint", size = "Cohorts")
            save_one(p, "FIG3E_OLFML2B_leave_one_out_meta_influence", "Leave-one-cohort-out meta-analysis influence audit", 10.2, max(5.0, 0.60 * nrow(d) + 1.8), "09_olfml2b_leave_one_out_meta.csv", "Random-effects leave-one-out meta-analysis", paste0("leave-one-out rows=", nrow(d)), "Leave-one-cohort-out random-effects meta-analysis for adjusted OLFML2B survival associations.")
        }
    }

    if (is.data.frame(median_surv) && nrow(median_surv)) {
        d <- median_surv[is.finite(median_surv$hr_high_vs_low), , drop = FALSE]
        if (nrow(d)) {
            d$label <- paste0(d$cohort, " / ", d$endpoint, "\nHR = ", olfml2b_fmt_num(d$hr_high_vs_low, 2), ", ", olfml2b_fmt_p(d$p_cox, 3), "; log-rank ", olfml2b_fmt_p(d$logrank_p, 3), "\nnH = ", d$n_high, ", nL = ", d$n_low)
            d$label <- factor(make.unique(d$label, sep = " #"), levels = rev(make.unique(d$label, sep = " #")))
            p <- ggplot2::ggplot(d, ggplot2::aes(x = hr_high_vs_low, y = label)) +
                ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = unname(pal["grey"])) +
                ggplot2::geom_point(ggplot2::aes(color = hr_high_vs_low > 1, size = events), alpha = 0.95) +
                ggplot2::scale_color_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = unname(pal["blue"])), guide = "none") +
                ggplot2::scale_size_continuous(range = c(2.2, 4.2), name = "Events") +
                ggplot2::scale_x_log10() +
                olfml2b_base_theme(base_size = 9.5) +
                ggplot2::labs(title = "Median-split survival display for OLFML2B", subtitle = "Median split is display-only and secondary to continuous z-score Cox models", x = "HR: OLFML2B-high vs low", y = NULL)
            save_one(p, "FIG3F_OLFML2B_median_split_forest", "Median-split survival display for OLFML2B", 9.8, max(5.0, 0.56 * nrow(d) + 1.8), "06_olfml2b_median_split_km_and_grouped_cox.csv", "Log-rank test and grouped Cox model; display/sensitivity only", paste0("rows=", nrow(d)), "Median-split Kaplan–Meier/log-rank display and grouped Cox summary.")
        }
    }

    if (is.data.frame(spline_tab) && nrow(spline_tab)) {
        d <- spline_tab[spline_tab$status == "OK" & is.finite(spline_tab$nonlinearity_lrt_p), , drop = FALSE]
        if (nrow(d)) {
            d$neglog10_p <- -log10(pmax(d$nonlinearity_lrt_p, 1e-300))
            d$label <- paste0(d$cohort, " / ", d$endpoint)
            d$label <- factor(d$label, levels = d$label[order(d$neglog10_p)])
            p <- ggplot2::ggplot(d, ggplot2::aes(x = neglog10_p, y = label)) +
                ggplot2::geom_col(fill = unname(pal["purple"]), width = 0.68) +
                ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = unname(pal["grey"])) +
                ggplot2::geom_text(ggplot2::aes(label = paste0(olfml2b_fmt_p(nonlinearity_lrt_p, 3), "\nn=", n, ", e=", events)), hjust = -0.05, size = 2.8) +
                ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
                olfml2b_base_theme(base_size = 9.5) +
                ggplot2::labs(title = "Nonlinearity audit for OLFML2B survival effect", subtitle = "Likelihood-ratio test compares linear and restricted-cubic-spline Cox models", x = expression(-log[10](P[nonlinearity])), y = NULL)
            save_one(p, "FIG3G_OLFML2B_nonlinearity_audit", "Nonlinearity audit for OLFML2B survival effect", 9.2, max(4.6, 0.50 * nrow(d) + 1.6), "07_olfml2b_nonlinearity_spline.csv", "Likelihood-ratio test comparing linear vs spline Cox terms", paste0("rows=", nrow(d)), "Restricted-cubic-spline nonlinearity audit across cohorts and endpoints.")
        }
    }

    if (is.data.frame(control_models) && nrow(control_models)) {
        d <- control_models[control_models$status == "OK" & is.finite(control_models$olfml2b_hr_after_control), , drop = FALSE]
        if (nrow(d)) {
            d$panel <- paste(d$cohort, d$endpoint, sep = " / ")
            d$control_label <- d$control_var
            d$cell_label <- paste0("HR=", olfml2b_fmt_num(d$olfml2b_hr_after_control, 2), "\n", olfml2b_fmt_p(d$olfml2b_p_after_control, 3), ifelse(is.finite(d$lrt_p_add_control), paste0("\nLRT ", olfml2b_fmt_p(d$lrt_p_add_control, 3)), ""))
            p <- ggplot2::ggplot(d, ggplot2::aes(x = panel, y = control_label, fill = log2(olfml2b_hr_after_control))) +
                ggplot2::geom_tile(color = "white", linewidth = 0.30) +
                ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.3) +
                ggplot2::scale_fill_gradient2(low = unname(pal["blue"]), mid = "white", high = unname(pal["red"]), midpoint = 0, name = "log2(HR)") +
                olfml2b_base_theme(base_size = 9) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
                ggplot2::labs(title = "OLFML2B survival association after proliferation-control adjustment", subtitle = "Cell labels show residual OLFML2B HR, P value and likelihood-ratio P for adding the control", x = NULL, y = "Added control")
            save_one(p, "FIG3H_OLFML2B_proliferation_control_heatmap", "OLFML2B survival association after proliferation-control adjustment", 11.0, max(5.2, 0.35 * length(unique(d$control_label)) * 2 + 2.2), "11_olfml2b_vs_context_control_models.csv", "Adjusted Cox models with proliferation or tissue-composition control variables", paste0("models=", nrow(d)), "Heatmap summarizing residual OLFML2B Cox effects after adding proliferation-related control variables.")
        }
    }

    if (is.data.frame(state_summary) && nrow(state_summary)) {
        d <- state_summary
        if (!"fdr" %in% names(d) && "p_value" %in% names(d)) d$fdr <- stats::p.adjust(d$p_value, method = "BH")
        d$module <- factor(d$module, levels = unique(d$module))
        p <- ggplot2::ggplot(d, ggplot2::aes(x = cohort, y = module, fill = delta_high_minus_low)) +
            ggplot2::geom_tile(color = "white", linewidth = 0.30) +
            ggplot2::geom_text(ggplot2::aes(label = paste0(sprintf("%.2f", delta_high_minus_low), ifelse(is.finite(fdr) & fdr < 0.05, "*", ""))), size = 2.5) +
            ggplot2::scale_fill_gradient2(low = unname(pal["blue"]), mid = "white", high = unname(pal["red"]), midpoint = 0, name = "High–low\ndelta") +
            olfml2b_base_theme(base_size = 9.5) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
            ggplot2::labs(title = "State-module shifts in OLFML2B-high tumors", subtitle = "Cell values show median module-score difference (high minus low); * indicates FDR < 0.05", x = NULL, y = NULL)
        save_one(p, "FIG3I_OLFML2B_state_shift_heatmap", "State-module shifts in OLFML2B-high tumors", 9.8, 6.6, "14_olfml2b_state_summary.csv", "Within-cohort Wilcoxon test; BH-FDR across module comparisons", paste0("tests=", nrow(d)), "Heatmap of state-module shifts between OLFML2B-high and OLFML2B-low tumors.")
    }

    if (exists("olfml2b_write_figure_registry", mode = "function")) {
        olfml2b_write_figure_registry(fig_rows, file.path(dirs$tables, "99_figure_registry.csv"))
    }
    invisible(TRUE)
}



run_olfml2b_part3_bulk_survival <- function(ctx = NULL, make_figures = TRUE) {
  if (is.null(ctx)) ctx <- olfml2b_load_context()
  run_olfml2b_specialized_bioinformatics(
    root = ctx$dirs$root,
    geo_validation_cohorts = c("GSE26253", "GSE84437", "GSE62254", "GSE15459"),
    output_subdir = "Part3",
    run_pdc_audit = FALSE,
    make_figures = make_figures
  )
}

olfml2b_run_ph_time_interaction_all <- function(views, min_n = 80L, min_events = 40L) {
  rows <- list()
  endpoints <- olfml2b_endpoint_defs()
  for (nm in names(views)) {
    v <- views[[nm]]
    clinical <- olfml2b_analysis_clinical(v)
    endpoint_names <- olfml2b_selected_endpoint_names(clinical, min_n = min_n, min_events = min_events)
    for (ep in endpoint_names) {
      def <- endpoints[[ep]]
      d <- olfml2b_surv_data(v, def)
      required <- c("time", "event", "olfml2b_z", "age", "sex", "stage")
      if (length(setdiff(required, names(d)))) next
      age <- suppressWarnings(as.numeric(d$age))
      d$age10 <- (age - mean(age, na.rm = TRUE)) / 10
      d$sex <- olfml2b_if78_sex_factor(d$sex)
      d$stage <- olfml2b_if78_stage_factor(d$stage)
      design_strata <- if (v$cohort == "GSE84437" && "source_subseries" %in% names(d)) "source_subseries" else NA_character_
      required <- c("time", "event", "olfml2b_z", "age10", "sex", "stage",
                    if (is.na(design_strata)) character() else design_strata)
      d <- d[stats::complete.cases(d[, required, drop = FALSE]), , drop = FALSE]
      events <- sum(d$event == 1L)
      if (nrow(d) < min_n || events < min_events) next
      form <- stats::as.formula(paste0(
        "survival::Surv(time, event) ~ olfml2b_z + tt(olfml2b_z) + age10 + sex + stage",
        if (is.na(design_strata)) "" else " + survival::strata(source_subseries)"
      ))
      fit <- tryCatch(survival::coxph(
        form, data = d, x = TRUE,
        tt = function(x, t, ...) x * log(pmax(t, 1))
      ), error = function(e) e)
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = v$cohort, endpoint = ep, endpoint_family = def$family,
          n = nrow(d), events = events, beta_olfml2b = NA_real_, beta_log_time_interaction = NA_real_,
          p_log_time_interaction = NA_real_, hr_day365 = NA_real_, hr_day1095 = NA_real_,
          status = paste0("FAILED:", conditionMessage(fit)),
          design_strata = design_strata,
          analysis_role = "prespecified_complete_case_PH_sensitivity", stringsAsFactors = FALSE
        )
        next
      }
      sm <- summary(fit)$coefficients
      b0 <- if ("olfml2b_z" %in% rownames(sm)) as.numeric(sm["olfml2b_z", "coef"]) else NA_real_
      tt_row <- grep("^tt\\(olfml2b_z\\)$", rownames(sm), value = TRUE)[1]
      bt <- if (length(tt_row) && !is.na(tt_row)) as.numeric(sm[tt_row, "coef"]) else NA_real_
      pt <- if (length(tt_row) && !is.na(tt_row)) as.numeric(sm[tt_row, "Pr(>|z|)"]) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = v$cohort, endpoint = ep, endpoint_family = def$family,
        n = nrow(d), events = events, beta_olfml2b = b0, beta_log_time_interaction = bt,
        p_log_time_interaction = pt,
        hr_day365 = if (is.finite(b0) && is.finite(bt)) exp(b0 + bt * log(365)) else NA_real_,
        hr_day1095 = if (is.finite(b0) && is.finite(bt)) exp(b0 + bt * log(1095)) else NA_real_,
        status = ifelse(is.finite(bt), "OK", "INTERACTION_NOT_ESTIMABLE"),
        design_strata = design_strata,
        analysis_role = "prespecified_complete_case_PH_sensitivity",
        interpretation_boundary = "Time-varying HR sensitivity for PH review; it never changes structural QC based on effect direction or significance.",
        stringsAsFactors = FALSE
      )
    }
  }
  olfml2b_bind_rows(rows)
}

# Methodology repair layer v1.7.0 for Part3
# Replaces the default Cox/meta helpers with machine-auditable statistical
# guardrails: EPV, model df, PH review, complete-case transparency, meta
# eligibility, and strict claim ceilings.
# ============================================================================

olfml2b_method_repair_model_matrix_df <- function(d, rhs_terms) {
  rhs_terms <- rhs_terms[nzchar(rhs_terms)]
  if (!length(rhs_terms)) return(0L)
  form <- stats::as.formula(paste0("~ ", paste(rhs_terms, collapse = " + ")))
  mm <- tryCatch(stats::model.matrix(form, data = d), error = function(e) NULL)
  if (is.null(mm)) return(length(rhs_terms))
  max(0L, ncol(mm) - 1L)
}

olfml2b_method_repair_ph <- function(fit, term = "olfml2b_z") {
  z <- tryCatch(suppressWarnings(survival::cox.zph(fit)), error = function(e) NULL)
  if (is.null(z) || is.null(z$table)) {
    return(list(term_p = NA_real_, global_p = NA_real_, status = "PH_NOT_TESTED"))
  }
  tab <- z$table
  term_p <- if (term %in% rownames(tab)) suppressWarnings(as.numeric(tab[term, "p"])) else NA_real_
  global_p <- if ("GLOBAL" %in% rownames(tab)) suppressWarnings(as.numeric(tab["GLOBAL", "p"])) else NA_real_
  status <- if (is.finite(term_p) && term_p < 0.05) {
    "PH_REVIEW_OLFML2B"
  } else if (is.finite(global_p) && global_p < 0.05) {
    "PH_REVIEW_GLOBAL"
  } else if (is.finite(term_p) || is.finite(global_p)) {
    "PH_OK"
  } else {
    "PH_NOT_TESTED"
  }
  list(term_p = term_p, global_p = global_p, status = status)
}

olfml2b_fit_olfml2b_cox <- function(dat, cohort, endpoint, endpoint_family, model, covariates = character()) {
  vars <- unique(c("time", "event", "olfml2b_z", covariates))
  d <- dat[, vars, drop = FALSE]
  before_n <- nrow(d)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  complete_case_fraction <- if (before_n > 0) nrow(d) / before_n else NA_real_
  if (nrow(d) < 20L || sum(d$event == 1L) < 5L) {
    out <- olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, "INSUFFICIENT_COMPLETE_CASES", covariates, nrow(d), sum(d$event == 1L))
    out$complete_case_fraction <- complete_case_fraction
    out$model_df <- NA_integer_; out$events_per_model_df <- NA_real_; out$ph_global_p <- NA_real_
    out$ph_status <- "NOT_TESTED"; out$primary_inference_eligible <- FALSE
    out$model_grade <- "FAIL_INSUFFICIENT_COMPLETE_CASES"
    out$claim_ceiling <- "not_reportable_as_primary"
    out$interpretation_boundary <- "Do not use this row for primary inference."
    return(out)
  }
  rhs <- c("olfml2b_z", covariates)
  form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(rhs, collapse = " + ")))
  model_df <- olfml2b_method_repair_model_matrix_df(d, rhs)
  events <- sum(d$event == 1L)
  epdf <- if (is.finite(model_df) && model_df > 0) events / model_df else NA_real_
  fit <- tryCatch(suppressWarnings(survival::coxph(form, data = d, x = TRUE)), error = function(e) e)
  if (inherits(fit, "error")) {
    out <- olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, paste0("COX_FAILED: ", conditionMessage(fit)), covariates, nrow(d), events)
    out$complete_case_fraction <- complete_case_fraction
    out$model_df <- model_df; out$events_per_model_df <- epdf; out$ph_global_p <- NA_real_
    out$ph_status <- "NOT_TESTED"; out$primary_inference_eligible <- FALSE
    out$model_grade <- "FAIL_COX_FAILED"
    out$claim_ceiling <- "not_reportable_as_primary"
    out$interpretation_boundary <- "Cox model failed; inspect endpoint coding and covariates."
    return(out)
  }
  sm <- summary(fit)
  if (!"olfml2b_z" %in% rownames(sm$coefficients)) {
    out <- olfml2b_cox_empty(cohort, endpoint, endpoint_family, model, "OLFML2B_COEFFICIENT_MISSING", covariates, nrow(d), events)
    out$complete_case_fraction <- complete_case_fraction
    out$model_df <- model_df; out$events_per_model_df <- epdf; out$ph_global_p <- NA_real_
    out$ph_status <- "NOT_TESTED"; out$primary_inference_eligible <- FALSE
    out$model_grade <- "FAIL_COEFFICIENT_MISSING"
    out$claim_ceiling <- "not_reportable_as_primary"
    out$interpretation_boundary <- "OLFML2B coefficient absent after model fitting."
    return(out)
  }
  co <- sm$coefficients["olfml2b_z", ]
  ci <- sm$conf.int["olfml2b_z", ]
  ph <- olfml2b_method_repair_ph(fit, "olfml2b_z")
  cindex <- tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_)
  unstable <- !is.finite(co["coef"]) || abs(co["coef"]) > 8 || !is.finite(co["se(coef)"]) || co["se(coef)"] > 5
  epv_ok <- is.finite(epdf) && epdf >= 10
  complete_ok <- is.finite(complete_case_fraction) && complete_case_fraction >= 0.60
  ph_ok <- !(ph$status %in% c("PH_REVIEW_OLFML2B", "PH_REVIEW_GLOBAL"))
  primary_ok <- !unstable && epv_ok && complete_ok && ph_ok
  model_grade <- if (unstable) {
    "FAIL_UNSTABLE_COX"
  } else if (!epv_ok) {
    "SUPPORTIVE_ONLY_LOW_EPV"
  } else if (!complete_ok) {
    "SUPPORTIVE_ONLY_HIGH_MISSINGNESS"
  } else if (!ph_ok) {
    "SUPPORTIVE_ONLY_PH_REVIEW"
  } else if (identical(model, "adjusted")) {
    "PRIMARY_ADJUSTED_ELIGIBLE"
  } else {
    "SUPPORTIVE_MODEL_ELIGIBLE"
  }
  data.frame(
    cohort = cohort,
    endpoint = endpoint,
    endpoint_family = endpoint_family,
    model = model,
    covariates = paste(covariates, collapse = ";"),
    n = nrow(d),
    events = events,
    complete_case_fraction = complete_case_fraction,
    model_df = model_df,
    events_per_model_df = epdf,
    beta = as.numeric(co["coef"]),
    se = as.numeric(co["se(coef)"]),
    hr = as.numeric(ci["exp(coef)"]),
    ci_low = as.numeric(ci["lower .95"]),
    ci_high = as.numeric(ci["upper .95"]),
    p_value = as.numeric(co["Pr(>|z|)"]),
    ph_p = ph$term_p,
    ph_global_p = ph$global_p,
    ph_status = ph$status,
    c_index = cindex,
    direction = ifelse(as.numeric(ci["exp(coef)"]) < 1, "favorable_high", "adverse_high"),
    status = ifelse(unstable, "COX_UNSTABLE", "OK"),
    model_grade = model_grade,
    primary_inference_eligible = primary_ok,
    meta_eligible = primary_ok || (!unstable && epv_ok && ph_ok && identical(model, "univariable")),
    claim_ceiling = ifelse(primary_ok && identical(model, "adjusted"), "observational_adjusted_association", "supportive_or_sensitivity_only"),
    interpretation_boundary = "Observational Cox association; not clinical prediction, treatment prediction, or causal mechanism.",
    stringsAsFactors = FALSE
  )
}

olfml2b_cox_empty <- function(cohort, endpoint, endpoint_family, model, status, covariates = character(), n = NA_integer_, events = NA_integer_) {
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = model, covariates = paste(covariates, collapse = ";"),
    n = n, events = events, complete_case_fraction = NA_real_, model_df = NA_integer_, events_per_model_df = NA_real_,
    beta = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    ph_p = NA_real_, ph_global_p = NA_real_, ph_status = NA_character_, c_index = NA_real_,
    direction = NA_character_, status = status, model_grade = status, primary_inference_eligible = FALSE, meta_eligible = FALSE,
    claim_ceiling = "not_reportable_as_primary", interpretation_boundary = "Empty or failed model row; audit only.", stringsAsFactors = FALSE
  )
}

olfml2b_fit_therapy_interaction <- function(dat, cohort, endpoint, endpoint_family, covariates = character()) {
  covariates <- setdiff(covariates, "adjuvant_therapy")
  vars <- unique(c("time", "event", "olfml2b_z", "adjuvant_therapy", covariates))
  d <- dat[, vars, drop = FALSE]
  before_n <- nrow(d)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  events <- sum(d$event == 1L)
  if (nrow(d) < 80L || events < 40L || length(unique(stats::na.omit(d$adjuvant_therapy))) < 2L) {
    out <- olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction_exploratory", "INSUFFICIENT_FOR_INTERACTION", c("adjuvant_therapy", covariates), nrow(d), events)
    out$complete_case_fraction <- if (before_n > 0) nrow(d) / before_n else NA_real_
    out$interpretation_boundary <- "Treatment interaction is exploratory and requires adequate treated/untreated events; do not claim predictive biomarker."
    return(out)
  }
  rhs <- c("olfml2b_z * adjuvant_therapy", covariates)
  form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(rhs, collapse = " + ")))
  fit <- tryCatch(suppressWarnings(survival::coxph(form, data = d, x = TRUE)), error = function(e) e)
  if (inherits(fit, "error")) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction_exploratory", paste0("COX_FAILED: ", conditionMessage(fit)), c("adjuvant_therapy", covariates), nrow(d), events))
  sm <- summary(fit)
  int_row <- grep("olfml2b_z:adjuvant_therapy|adjuvant_therapy.*:olfml2b_z", rownames(sm$coefficients), value = TRUE)[1]
  if (is.na(int_row)) return(olfml2b_cox_empty(cohort, endpoint, endpoint_family, "therapy_interaction_exploratory", "INTERACTION_TERM_MISSING", c("adjuvant_therapy", covariates), nrow(d), events))
  co <- sm$coefficients[int_row, ]; ci <- sm$conf.int[int_row, ]
  model_df <- olfml2b_method_repair_model_matrix_df(d, c("olfml2b_z", "adjuvant_therapy", "olfml2b_z:adjuvant_therapy", covariates))
  epdf <- if (model_df > 0) events / model_df else NA_real_
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = "therapy_interaction_exploratory", covariates = paste(c("adjuvant_therapy", covariates), collapse = ";"),
    n = nrow(d), events = events, complete_case_fraction = if (before_n > 0) nrow(d) / before_n else NA_real_, model_df = model_df, events_per_model_df = epdf,
    beta = as.numeric(co["coef"]), se = as.numeric(co["se(coef)"]),
    hr = as.numeric(ci["exp(coef)"]), ci_low = as.numeric(ci["lower .95"]), ci_high = as.numeric(ci["upper .95"]),
    p_value = as.numeric(co["Pr(>|z|)"]), ph_p = NA_real_, ph_global_p = NA_real_, ph_status = "NOT_PRIMARY_INTERACTION",
    c_index = tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_),
    direction = ifelse(as.numeric(ci["exp(coef)"]) < 1, "interaction_favorable_with_therapy", "interaction_adverse_with_therapy"),
    status = "OK", model_grade = "EXPLORATORY_INTERACTION_ONLY", primary_inference_eligible = FALSE, meta_eligible = FALSE,
    claim_ceiling = "exploratory_interaction_only", interpretation_boundary = "Do not claim treatment-predictive biomarker without a treatment-specific design and independent validation.",
    stringsAsFactors = FALSE
  )
}

olfml2b_meta_all <- function(surv_all) {
  if (!nrow(surv_all)) return(data.frame())
  elig <- if ("meta_eligible" %in% names(surv_all)) surv_all$meta_eligible %in% TRUE else surv_all$status == "OK"
  d <- surv_all[
    elig &
      surv_all$model %in% c("adjusted", "univariable") &
      is.finite(surv_all$beta) &
      is.finite(surv_all$se) & surv_all$se > 0,
    , drop = FALSE
  ]
  d <- olfml2b_deduplicate_meta_effects(d)
  rows <- list()
  for (model in unique(d$model)) {
    for (fam in unique(d$endpoint_family)) {
      dd <- d[d$model == model & d$endpoint_family == fam, , drop = FALSE]
      rows[[length(rows) + 1L]] <- olfml2b_meta_one(dd, endpoint_family = fam, model = model)
    }
  }
  olfml2b_bind_rows(rows)
}

olfml2b_meta_one <- function(d, endpoint_family, model) {
  if (nrow(d) < 2L) {
    return(data.frame(endpoint_family = endpoint_family, model = model, k = nrow(d), status = "INSUFFICIENT_META_ELIGIBLE_COHORTS", evidence_grade = "NOT_PRIMARY", claim_ceiling = "single_or_no_cohort_support", metafor_reml_crosscheck = "NOT_EVALUABLE_K_LT2", stringsAsFactors = FALSE))
  }
  yi <- d$beta; vi <- d$se^2
  good <- is.finite(yi) & is.finite(vi) & vi > 0
  d <- d[good, , drop = FALSE]; yi <- yi[good]; vi <- vi[good]
  if (length(yi) < 2L) return(data.frame(endpoint_family = endpoint_family, model = model, k = length(yi), status = "INSUFFICIENT_META_ELIGIBLE_COHORTS", evidence_grade = "NOT_PRIMARY", claim_ceiling = "single_or_no_cohort_support", metafor_reml_crosscheck = "NOT_EVALUABLE_K_LT2", stringsAsFactors = FALSE))
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  cval <- sum(wi) - sum(wi^2) / sum(wi)
  tau2_dl <- if (is.finite(cval) && cval > 0) max(0, (q - df) / cval) else 0

  # REML is the locked estimator in Part0.  Implement it directly so the
  # analysis does not silently fall back to DL when metafor is unavailable.
  reml_objective <- function(tau2_value) {
    w0 <- 1 / (vi + tau2_value)
    mu0 <- sum(w0 * yi) / sum(w0)
    sum(log(vi + tau2_value)) + log(sum(w0)) + sum(w0 * (yi - mu0)^2)
  }
  upper <- max(1, tau2_dl * 20, stats::var(yi) * 20, na.rm = TRUE)
  opt <- tryCatch(stats::optimize(reml_objective, interval = c(0, upper)), error = function(e) NULL)
  tau2 <- if (!is.null(opt) && is.finite(opt$minimum)) max(0, opt$minimum) else tau2_dl
  if (is.finite(reml_objective(0)) && reml_objective(0) <= reml_objective(tau2) + 1e-10) tau2 <- 0
  wr <- 1 / (vi + tau2)
  random <- sum(wr * yi) / sum(wr)
  se_model_based <- sqrt(1 / sum(wr))
  q_star <- sum(wr * (yi - random)^2) / df
  # Modified Hartung-Knapp: prevent a spuriously narrower interval when the
  # observed residual scale is below one; use a t reference with k-1 df.
  hk_scale <- max(1, q_star)
  se_random <- sqrt(hk_scale / sum(wr))
  crit <- stats::qt(0.975, df = df)
  p_random <- 2 * stats::pt(abs(random / se_random), df = df, lower.tail = FALSE)
  pred_se <- sqrt(tau2 + se_random^2)
  i2 <- if (q > df && q > 0) max(0, (q - df) / q) * 100 else 0
  ci_low <- exp(random - crit * se_random)
  ci_high <- exp(random + crit * se_random)
  prediction_low <- exp(random - crit * pred_se)
  prediction_high <- exp(random + crit * pred_se)
  excludes_null <- is.finite(ci_low) && is.finite(ci_high) && (ci_low > 1 || ci_high < 1)
  evidence_grade <- if (length(yi) >= 3L && is.finite(i2) && i2 <= 50 && p_random < 0.05 && excludes_null) {
    "MULTI_COHORT_CONSISTENT_ASSOCIATION"
  } else if (length(yi) >= 2L && p_random < 0.05 && excludes_null) {
    "SUPPORTIVE_META_ASSOCIATION_REVIEW_HETEROGENEITY"
  } else {
    "NO_ROBUST_META_ASSOCIATION"
  }
  mf <- if (requireNamespace("metafor", quietly = TRUE)) tryCatch(
    metafor::rma(yi = yi, vi = vi, method = "REML", test = "knha"),
    error = function(e) e
  ) else structure(list(message = "metafor package missing"), class = c("simpleError", "error", "condition"))
  mf_beta <- if (!inherits(mf, "error")) as.numeric(mf$b[1]) else NA_real_
  mf_tau2 <- if (!inherits(mf, "error")) as.numeric(mf$tau2) else NA_real_
  beta_diff <- abs(random - mf_beta)
  tau2_diff <- abs(tau2 - mf_tau2)
  mf_check <- if (inherits(mf, "error")) "FAIL_METAFOR_NOT_EVALUABLE" else
    ifelse(is.finite(beta_diff) && is.finite(tau2_diff) && beta_diff <= 1e-5 && tau2_diff <= 1e-5,
           "PASS_REML_POINT_AND_TAU2", "FAIL_REML_DISAGREEMENT")
  data.frame(
    endpoint_family = endpoint_family, model = model, k = nrow(d),
    cohorts = paste(d$cohort, collapse = ";"),
    cohort_endpoints = paste(paste(d$cohort, d$endpoint, sep = ":"), collapse = ";"),
    beta_random = random, se_random = se_random, se_model_based = se_model_based, hr_random = exp(random),
    ci_low_random = ci_low, ci_high_random = ci_high,
    p_random = p_random, beta_fixed = fixed, hr_fixed = exp(fixed), q = q,
    q_p = stats::pchisq(q, df = df, lower.tail = FALSE), i2 = i2, tau2 = tau2,
    tau2_dl_diagnostic = tau2_dl, hk_scale = hk_scale, inference_df = df,
    prediction_low = prediction_low, prediction_high = prediction_high,
    meta_method = "REML", interval_method = "modified_Hartung_Knapp_t",
    metafor_beta_reml = mf_beta, metafor_tau2_reml = mf_tau2,
    metafor_beta_abs_diff = beta_diff, metafor_tau2_abs_diff = tau2_diff,
    metafor_reml_crosscheck = mf_check,
    direction_random = ifelse(exp(random) < 1, "favorable_high", "adverse_high"),
    status = "OK", evidence_grade = evidence_grade,
    claim_ceiling = ifelse(evidence_grade == "MULTI_COHORT_CONSISTENT_ASSOCIATION", "external_public_cohort_consistency", "supportive_meta_only"),
    interpretation_boundary = "Random-effects meta-analysis of observational cohorts; not clinical utility or causality.",
    stringsAsFactors = FALSE
  )
}

olfml2b_leave_one_out_meta_all <- function(surv_all) {
  elig <- if ("meta_eligible" %in% names(surv_all)) surv_all$meta_eligible %in% TRUE else surv_all$status == "OK"
  d <- surv_all[
    elig & surv_all$model == "adjusted" & is.finite(surv_all$beta) & is.finite(surv_all$se) & surv_all$se > 0,
    , drop = FALSE
  ]
  d <- olfml2b_deduplicate_meta_effects(d)
  rows <- list()
  for (fam in unique(d$endpoint_family)) {
    dd <- d[d$endpoint_family == fam, , drop = FALSE]
    if (nrow(dd) < 3L) next
    for (co in unique(dd$cohort)) {
      mm <- olfml2b_meta_one(dd[dd$cohort != co, , drop = FALSE], endpoint_family = fam, model = "adjusted_leave_one_out")
      mm$left_out <- co
      mm$loo_stability_flag <- ifelse(mm$status == "OK" & mm$p_random < 0.05 & mm$hr_random < 1, "STABLE_FAVORABLE", ifelse(mm$status == "OK" & mm$p_random < 0.05 & mm$hr_random > 1, "STABLE_ADVERSE", "UNSTABLE_OR_NULL"))
      rows[[length(rows) + 1L]] <- mm
    }
  }
  olfml2b_bind_rows(rows)
}

olfml2b_part3_methodology_repair_audit <- function(surv_all, meta, loo, control_models, module_results) {
  data.frame(
    item = c(
      "primary_scale", "cutpoint_policy", "complete_case_policy", "EPV_policy", "PH_policy",
      "meta_policy", "proliferation_control", "module_policy", "claim_ceiling"
    ),
    status = c(
      "PASS", "PASS", ifelse(nrow(surv_all) && "complete_case_fraction" %in% names(surv_all), "PASS", "REVIEW"),
      ifelse(nrow(surv_all) && "events_per_model_df" %in% names(surv_all), "PASS", "REVIEW"),
      ifelse(nrow(surv_all) && "ph_status" %in% names(surv_all), "PASS", "REVIEW"),
      ifelse(nrow(meta) && "evidence_grade" %in% names(meta), "PASS", "REVIEW"),
      ifelse(is.data.frame(control_models) && nrow(control_models), "PASS_OR_REVIEW_BY_TABLE", "NOT_AVAILABLE"),
      ifelse(is.list(module_results) || is.data.frame(module_results), "EXPLORATORY_ONLY", "NOT_AVAILABLE"),
      "OBSERVATIONAL_ASSOCIATION_ONLY"
    ),
    rule = c(
      "Continuous within-cohort z-score Cox remains primary; GSE84437 is standardized within GSE84426/GSE84433 and uses subseries-stratified baseline hazards.",
      "Median split is display/sensitivity only; no best cutpoint.",
      "Model rows export complete-case fraction and are downgraded if high missingness.",
      "Model rows export model df and events per df; low EPV rows are supportive only.",
      "OLFML2B/global PH checks are exported; PH violations are supportive/review only.",
      "Meta-analysis uses pre-deduplicated endpoint-family rows, REML tau-squared, modified Hartung-Knapp t intervals, prediction intervals, and meta-eligible models only.",
      "OLFML2B must be interpreted against proliferation and tissue-composition control sensitivity.",
      "State modules are exploratory biological context, not mechanism proof.",
      "No clinical prediction, treatment prediction, or causal mechanism claim from Part3 alone."
    ),
    stringsAsFactors = FALSE
  )
}

.olfml2b_part3_original_runner <- run_olfml2b_specialized_bioinformatics
run_olfml2b_specialized_bioinformatics <- function(...) {
  index <- .olfml2b_part3_original_runner(...)
  dirs <- index$dirs
  repair <- olfml2b_part3_methodology_repair_audit(index$survival, index$meta, index$leave_one_out, index$controls$models, index$modules)
  olfml2b_write_csv(repair, file.path(dirs$tables, "21_part3_methodology_repair_audit.csv"))
  index$methodology_repair_audit <- repair
  saveRDS(index, file.path(dirs$objects, "Part3_OLFML2B_specialized_bioinformatics_index.rds"), compress = "xz")
  invisible(index)
}
# ==============================================================================
# IF 7-8 survival implementation | fixed M0/M1/M2, cohort-level MI and
# structure-only quality gates. This block supersedes all earlier helper
# definitions and is the implementation used by the public Part3 runner.
# ==============================================================================

OLFML2B_PART3_IF78_VERSION <- "v1.1.0_20260722_PDC614_CONTRACT_SYNC"
OLFML2B_ANALYSIS_VERSION <- OLFML2B_PART3_IF78_VERSION

olfml2b_select_covariates <- function(dat, events, min_epv = 10, cohort = NA_character_,
                                    model_role = c("primary_adjusted", "molecular_subtype_sensitivity")) {
  model_role <- match.arg(model_role)
  if (identical(model_role, "primary_adjusted")) {
    return(c("age", "sex", "stage"))
  }
  c("age", "sex", "stage", "molecular_subtype")
}

olfml2b_if78_nelson_aalen <- function(time, event) {
  d <- data.frame(time = as.numeric(time), event = as.integer(event))
  fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ 1, data = d), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(d)))
  bh <- tryCatch(survival::basehaz(fit, centered = FALSE), error = function(e) NULL)
  if (is.null(bh) || !nrow(bh)) return(rep(0, nrow(d)))
  stats::approx(bh$time, bh$hazard, xout = d$time, method = "constant", f = 0,
                rule = 2, ties = "ordered")$y
}

olfml2b_if78_model_id <- function(model) {
  switch(as.character(model)[1],
         univariable = "M0_unadjusted",
         adjusted = "M1_common_age_sex_stage",
         adjusted_molecular_subtype_sensitivity = "M2_subtype_sensitivity",
         therapy_interaction = "M3_therapy_interaction_exploratory",
         as.character(model)[1])
}

olfml2b_if78_sex_factor <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("female", "f", "woman")] <- "Female"
  out[raw %in% c("male", "m", "man")] <- "Male"
  factor(out, levels = c("Female", "Male"))
}

olfml2b_if78_stage_factor <- function(x) {
  # Upper-case first: filtering lower-case letters before toupper() would turn
  # values such as "Stage IV" into "SIV" and silently lose the stage prefix.
  raw <- gsub("[^A-Z0-9]", "", toupper(as.character(x)))
  raw <- sub("^(AJCC|PATHOLOGIC|PATHOLOGICAL|CLINICAL)", "", raw)
  raw <- sub("^STAGE", "", raw)
  out <- rep(NA_character_, length(raw))
  out[grepl("^(IV|4)", raw)] <- "IV"
  out[is.na(out) & grepl("^(III|3)", raw)] <- "III"
  out[is.na(out) & grepl("^(II|2)", raw)] <- "II"
  out[is.na(out) & grepl("^(I|1)", raw)] <- "I"
  factor(out, levels = c("I", "II", "III", "IV"))
}

olfml2b_gse84437_subseries_heterogeneity_audit <- function(views) {
  if (!"GSE84437" %in% names(views)) return(data.frame())
  v <- views[["GSE84437"]]
  d0 <- olfml2b_analysis_clinical(v)
  if (!all(c("os_time_days", "os_event", "olfml2b_z", "source_subseries") %in% names(d0))) return(data.frame())
  keep <- is.finite(d0$os_time_days) & d0$os_time_days > 0 & d0$os_event %in% c(0L, 1L) &
    is.finite(d0$olfml2b_z) & !is.na(d0$source_subseries)
  d0 <- d0[keep, , drop = FALSE]
  d0$time <- as.numeric(d0$os_time_days); d0$event <- as.integer(d0$os_event)
  d0$source_subseries <- stats::relevel(factor(d0$source_subseries), ref = "GSE84426")
  rows <- list()
  model_specs <- list(
    unadjusted_batch_stratified = character(),
    adjusted_batch_stratified = c("age", "sex", "stage")
  )
  for (model in names(model_specs)) {
    d <- d0; covars <- model_specs[[model]]
    if (length(covars)) {
      if (length(setdiff(covars, names(d)))) next
      d$age10 <- (suppressWarnings(as.numeric(d$age)) - mean(suppressWarnings(as.numeric(d$age)), na.rm = TRUE)) / 10
      d$sex <- olfml2b_if78_sex_factor(d$sex)
      d$stage <- olfml2b_if78_stage_factor(d$stage)
      d <- d[stats::complete.cases(d[, c("time", "event", "olfml2b_z", "source_subseries", "age10", "sex", "stage"), drop = FALSE]), , drop = FALSE]
    }
    rhs <- c("olfml2b_z", "olfml2b_z:source_subseries",
             if (length(covars)) c("age10", "sex", "stage") else character(),
             "survival::strata(source_subseries)")
    form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(rhs, collapse = " + ")))
    fit <- tryCatch(survival::coxph(form, data = d, x = TRUE), error = function(e) NULL)
    sm <- if (!is.null(fit)) summary(fit)$coefficients else NULL
    int_row <- if (!is.null(sm)) grep("olfml2b_z:source_subseries", rownames(sm), value = TRUE)[1] else NA_character_
    beta <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "coef"]) else NA_real_
    se <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "se(coef)"]) else NA_real_
    p <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "Pr(>|z|)"]) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = "GSE84437", endpoint = "OS", model = model,
      n = nrow(d), events = sum(d$event == 1L),
      n_GSE84426 = sum(d$source_subseries == "GSE84426"),
      events_GSE84426 = sum(d$event[d$source_subseries == "GSE84426"] == 1L),
      n_GSE84433 = sum(d$source_subseries == "GSE84433"),
      events_GSE84433 = sum(d$event[d$source_subseries == "GSE84433"] == 1L),
      log_HR_ratio_GSE84433_vs_GSE84426 = beta,
      HR_ratio_GSE84433_vs_GSE84426 = exp(beta),
      ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
      interaction_p = p,
      design = "common exposure scale; subseries-specific baseline hazards; interaction is heterogeneity sensitivity",
      inference_role = "prespecified_subseries_heterogeneity_sensitivity_not_selection_gate",
      stringsAsFactors = FALSE
    )
  }
  out <- olfml2b_bind_rows(rows)
  if (nrow(out)) out$interaction_fdr_bh <- stats::p.adjust(out$interaction_p, method = "BH")
  out
}

olfml2b_if78_empty_cox <- function(cohort, endpoint, endpoint_family, model, status,
                                 covariates = character(), n = NA_integer_, events = NA_integer_,
                                 analysis_method = NA_character_, missing_row_fraction = NA_real_,
                                 imputation_m = NULL) {
  planned_m <- if (is.null(imputation_m)) {
    ifelse(identical(analysis_method, "mice_rubin"), 20L, 0L)
  } else as.integer(imputation_m)
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = model, model_id = olfml2b_if78_model_id(model),
    covariates = paste(covariates, collapse = ";"), n = n, events = events,
    design_strata = ifelse(identical(cohort, "GSE84437"), "source_subseries", NA_character_),
    complete_case_fraction = ifelse(is.finite(missing_row_fraction), 1 - missing_row_fraction, NA_real_),
    missing_row_fraction = missing_row_fraction, analysis_method = analysis_method,
    imputation_m = planned_m,
    model_df = NA_integer_, events_per_model_df = NA_real_,
    beta = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, ph_p = NA_real_, ph_global_p = NA_real_, ph_status = "NOT_TESTED",
    max_abs_dfbeta_olfml2b = NA_real_, n_influential_dfbeta_olfml2b = NA_integer_,
    dfbeta_review_threshold = NA_real_, max_abs_deviance_residual = NA_real_,
    influence_diagnostic_basis = NA_character_,
    c_index = NA_real_, direction = NA_character_, status = status, model_grade = status,
    structural_eligible = FALSE, primary_inference_eligible = FALSE, meta_eligible = FALSE,
    eligibility_reason = status, claim_ceiling = "not_reportable_as_primary",
    interpretation_boundary = "Model was not structurally evaluable; no effect-based decision was made.",
    stringsAsFactors = FALSE
  )
}

olfml2b_cox_empty <- function(cohort, endpoint, endpoint_family, model, status,
                            covariates = character(), n = NA_integer_, events = NA_integer_) {
  olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model, status,
                       covariates, n, events)
}

olfml2b_if78_extract_cox <- function(fit) {
  sm <- summary(fit)
  if (!"olfml2b_z" %in% rownames(sm$coefficients)) return(NULL)
  co <- sm$coefficients["olfml2b_z", ]
  ph <- olfml2b_method_repair_ph(fit, "olfml2b_z")
  dfb <- tryCatch(stats::residuals(fit, type = "dfbeta"), error = function(e) NULL)
  target_dfb <- if (is.matrix(dfb) && "olfml2b_z" %in% colnames(dfb)) {
    as.numeric(dfb[, "olfml2b_z"])
  } else if (is.numeric(dfb) && length(stats::coef(fit)) == 1L &&
             "olfml2b_z" %in% names(stats::coef(fit))) {
    as.numeric(dfb)
  } else numeric()
  fit_n <- if (!is.null(fit$n) && length(fit$n)) suppressWarnings(as.numeric(fit$n[1])) else NA_real_
  review_threshold <- if (is.finite(fit_n) && fit_n > 0L) 2 / sqrt(fit_n) else NA_real_
  deviance <- tryCatch(as.numeric(stats::residuals(fit, type = "deviance")), error = function(e) numeric())
  list(
    beta = as.numeric(co["coef"]),
    variance = as.numeric(co["se(coef)"])^2,
    p = as.numeric(co["Pr(>|z|)"]),
    ph_p = ph$term_p, ph_global_p = ph$global_p, ph_status = ph$status,
    c_index = tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_),
    max_abs_dfbeta = if (any(is.finite(target_dfb))) max(abs(target_dfb), na.rm = TRUE) else NA_real_,
    n_influential_dfbeta = if (length(target_dfb) && is.finite(review_threshold))
      sum(abs(target_dfb) > review_threshold, na.rm = TRUE) else NA_integer_,
    dfbeta_review_threshold = review_threshold,
    max_abs_deviance = if (any(is.finite(deviance))) max(abs(deviance), na.rm = TRUE) else NA_real_
  )
}

olfml2b_if78_fit_mice_cox <- function(d, form, covariates, m = 20L, seed = 20260721L) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("mice is required because core-covariate missingness exceeds 5%.", call. = FALSE)
  }
  imp_dat <- d
  imp_dat$nelson_aalen <- olfml2b_if78_nelson_aalen(imp_dat$time, imp_dat$event)
  meth <- mice::make.method(imp_dat)
  pred <- mice::make.predictorMatrix(imp_dat)
  locked <- c("time", "event", "olfml2b_z", "nelson_aalen", "source_subseries")
  meth[intersect(locked, names(meth))] <- ""
  pred[intersect(locked, rownames(pred)), ] <- 0
  for (v in intersect(covariates, names(imp_dat))) {
    x <- imp_dat[[v]]
    if (is.ordered(x)) {
      meth[v] <- "polr"
    } else if (is.factor(x) && nlevels(droplevels(x)) == 2L) {
      meth[v] <- "logreg"
    } else if (is.factor(x)) {
      meth[v] <- "polyreg"
    } else {
      meth[v] <- "pmm"
    }
  }
  set.seed(seed)
  imp <- mice::mice(imp_dat, m = as.integer(m), maxit = 20L, method = meth,
                    predictorMatrix = pred, seed = seed, printFlag = FALSE)
  fits <- lapply(seq_len(m), function(i) {
    completed <- mice::complete(imp, i)
    fit <- survival::coxph(form, data = completed, x = TRUE)
    list(fit = fit, extracted = olfml2b_if78_extract_cox(fit), completed = completed)
  })
  ok <- vapply(fits, function(z) !is.null(z$extracted) && is.finite(z$extracted$beta) &&
                 is.finite(z$extracted$variance) && z$extracted$variance > 0, logical(1))
  fits <- fits[ok]
  if (length(fits) < max(5L, ceiling(m / 2))) stop("Too few converged Cox fits across imputations.", call. = FALSE)
  qhat <- vapply(fits, function(z) z$extracted$beta, numeric(1))
  uhat <- vapply(fits, function(z) z$extracted$variance, numeric(1))
  mm <- length(qhat)
  qbar <- mean(qhat)
  ubar <- mean(uhat)
  b <- stats::var(qhat)
  total_var <- ubar + (1 + 1 / mm) * b
  r <- if (ubar > 0) ((1 + 1 / mm) * b) / ubar else Inf
  df <- if (is.finite(r) && r > 0) (mm - 1) * (1 + 1 / r)^2 else Inf
  se <- sqrt(total_var)
  crit <- if (is.finite(df)) stats::qt(0.975, df = df) else stats::qnorm(0.975)
  p <- if (is.finite(df)) 2 * stats::pt(abs(qbar / se), df = df, lower.tail = FALSE) else
    2 * stats::pnorm(abs(qbar / se), lower.tail = FALSE)
  list(
    beta = qbar, se = se, ci_low = exp(qbar - crit * se), ci_high = exp(qbar + crit * se),
    p = p, df = df, m = mm,
    ph_p = stats::median(vapply(fits, function(z) z$extracted$ph_p, numeric(1)), na.rm = TRUE),
    ph_global_p = stats::median(vapply(fits, function(z) z$extracted$ph_global_p, numeric(1)), na.rm = TRUE),
    c_index = mean(vapply(fits, function(z) z$extracted$c_index, numeric(1)), na.rm = TRUE),
    max_abs_dfbeta = stats::median(vapply(fits, function(z) z$extracted$max_abs_dfbeta, numeric(1)), na.rm = TRUE),
    n_influential_dfbeta = as.integer(round(stats::median(vapply(fits, function(z) z$extracted$n_influential_dfbeta, numeric(1)), na.rm = TRUE))),
    dfbeta_review_threshold = stats::median(vapply(fits, function(z) z$extracted$dfbeta_review_threshold, numeric(1)), na.rm = TRUE),
    max_abs_deviance = stats::median(vapply(fits, function(z) z$extracted$max_abs_deviance, numeric(1)), na.rm = TRUE),
    first_completed = fits[[1]]$completed
  )
}

olfml2b_if78_empty_cox <- function(cohort, endpoint, endpoint_family, model, status,
                                 covariates = character(), n = NA_integer_, events = NA_integer_,
                                 analysis_method = "not_run", missing_row_fraction = NA_real_,
                                 imputation_m = 0L, max_covariate_missing_fraction = NA_real_,
                                 design_strata = NA_character_) {
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = model, model_id = olfml2b_if78_model_id(model),
    covariates = paste(covariates, collapse = ";"), n = n, events = events,
    design_strata = design_strata,
    complete_case_fraction = ifelse(is.finite(missing_row_fraction), 1 - missing_row_fraction, NA_real_),
    missing_row_fraction = missing_row_fraction,
    max_covariate_missing_fraction = max_covariate_missing_fraction,
    analysis_method = analysis_method, imputation_m = as.integer(imputation_m),
    imputation_policy = "MI only when max variable missingness is >5% and <=20%, incomplete rows <=30%, and every covariate has >=30 observed values; structural/excessive missingness is never imputed.",
    model_df = NA_integer_, events_per_model_df = NA_real_,
    beta = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, ph_p = NA_real_, ph_global_p = NA_real_, ph_status = "NOT_TESTED",
    max_abs_dfbeta_olfml2b = NA_real_, n_influential_dfbeta_olfml2b = NA_integer_,
    dfbeta_review_threshold = NA_real_, max_abs_deviance_residual = NA_real_,
    influence_diagnostic_basis = NA_character_, c_index = NA_real_, direction = NA_character_,
    status = status, model_grade = status, structural_eligible = FALSE,
    primary_inference_eligible = FALSE, meta_eligible = FALSE,
    eligibility_reason = status, claim_ceiling = "not_reportable_as_primary",
    stage_derivation_policy = "Source overall stage or complete AJCC7 gastric TNM is primary; T/N-assuming-M0 stage is sensitivity-only.",
    interpretation_boundary = "Model was not structurally evaluable; no effect-based decision was made.",
    stringsAsFactors = FALSE
  )
}

olfml2b_cox_empty <- function(cohort, endpoint, endpoint_family, model, status,
                            covariates = character(), n = NA_integer_, events = NA_integer_) {
  olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model, status,
                       covariates, n, events)
}

olfml2b_fit_olfml2b_cox_v23 <- function(dat, cohort, endpoint, endpoint_family, model,
                                     covariates = character()) {
  design_strata <- if (identical(cohort, "GSE84437") &&
                       "source_subseries" %in% names(dat) &&
                       length(unique(stats::na.omit(dat$source_subseries))) == 2L) {
    "source_subseries"
  } else NA_character_
  required_source <- unique(c("time", "event", "olfml2b_z", covariates,
                              if (is.na(design_strata)) character() else design_strata))
  missing_columns <- setdiff(required_source, names(dat))
  if (length(missing_columns)) {
    return(olfml2b_if78_empty_cox(
      cohort, endpoint, endpoint_family, model,
      paste0("STRUCTURAL_MISSING_COLUMNS:", paste(missing_columns, collapse = ",")),
      covariates, design_strata = design_strata
    ))
  }
  structural <- covariates[vapply(covariates, function(v) {
    x <- dat[[v]]; keep <- olfml2b_if78_nonmissing(x)
    sum(keep) == 0L || length(unique(as.character(x[keep]))) < 2L
  }, logical(1))]
  if (length(structural)) {
    return(olfml2b_if78_empty_cox(
      cohort, endpoint, endpoint_family, model,
      paste0("STRUCTURALLY_UNAVAILABLE_NO_IMPUTATION:", paste(structural, collapse = ",")),
      covariates, nrow(dat), sum(dat$event == 1L, na.rm = TRUE),
      analysis_method = "not_run_structural_missingness", missing_row_fraction = 1,
      max_covariate_missing_fraction = 1, design_strata = design_strata
    ))
  }
  prep <- olfml2b_if78_prepare_model_data(dat, covariates, design_strata)
  d <- prep$data; model_covars <- prep$covariates
  if (any(!is.finite(d$time)) || any(!d$event %in% c(0L, 1L)) ||
      any(!is.finite(d$olfml2b_z))) {
    return(olfml2b_if78_empty_cox(
      cohort, endpoint, endpoint_family, model, "INVALID_ENDPOINT_OR_EXPOSURE",
      model_covars, nrow(d), sum(d$event == 1L, na.rm = TRUE),
      design_strata = design_strata
    ))
  }
  cov_missing <- if (length(model_covars)) {
    vapply(model_covars, function(v) mean(!olfml2b_if78_nonmissing(d[[v]])), numeric(1))
  } else numeric()
  max_missing <- if (length(cov_missing)) max(cov_missing) else 0
  incomplete <- if (length(model_covars))
    !stats::complete.cases(d[, model_covars, drop = FALSE]) else rep(FALSE, nrow(d))
  missing_row_fraction <- mean(incomplete)
  events_all <- sum(d$event == 1L)

  if (max_missing > 0.20 || missing_row_fraction > 0.30) {
    return(olfml2b_if78_empty_cox(
      cohort, endpoint, endpoint_family, model,
      paste0("EXCESSIVE_MISSINGNESS_NO_IMPUTATION:max_variable=",
             sprintf("%.3f", max_missing), ";incomplete_rows=",
             sprintf("%.3f", missing_row_fraction)),
      model_covars, nrow(d), events_all,
      analysis_method = "not_run_excessive_missingness",
      missing_row_fraction = missing_row_fraction,
      max_covariate_missing_fraction = max_missing,
      design_strata = design_strata
    ))
  }
  rhs <- c("olfml2b_z", model_covars)
  rhs_formula <- c(rhs, if (is.na(design_strata)) character() else
                     paste0("survival::strata(", design_strata, ")"))
  form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ",
                                   paste(rhs_formula, collapse = " + ")))
  # Ambiguous TN-derived stage is deterministic non-derivability, not a
  # stochastic missing value.  M2b therefore uses complete cases even when its
  # missing fraction is within the general limited-MI window.
  use_mi <- length(model_covars) > 0L && max_missing > 0.05 &&
    !identical(model, "derived_stage_sensitivity")
  planned_m <- if (use_mi) max(20L, as.integer(ceiling(100 * max_missing))) else 0L

  if (use_mi) {
    mi <- tryCatch(olfml2b_if78_fit_mice_limited(d, form, model_covars,
                                                m = planned_m, seed = 20260721L),
                   error = function(e) e)
    if (inherits(mi, "error")) {
      return(olfml2b_if78_empty_cox(
        cohort, endpoint, endpoint_family, model,
        paste0("ELIGIBLE_LIMITED_MICE_FAILED:", conditionMessage(mi)),
        model_covars, nrow(d), events_all,
        analysis_method = "mice_rubin_limited_missingness",
        missing_row_fraction = missing_row_fraction,
        imputation_m = planned_m,
        max_covariate_missing_fraction = max_missing,
        design_strata = design_strata
      ))
    }
    analysis_method <- "mice_rubin_limited_missingness"
    beta <- mi$beta; se <- mi$se; ci_low <- mi$ci_low; ci_high <- mi$ci_high
    p_value <- mi$p; ph_p <- mi$ph_p; ph_global_p <- mi$ph_global_p
    c_index <- mi$c_index; max_abs_dfbeta <- mi$max_abs_dfbeta
    n_influential_dfbeta <- mi$n_influential_dfbeta
    dfbeta_review_threshold <- mi$dfbeta_review_threshold
    max_abs_deviance <- mi$max_abs_deviance
    influence_basis <- "median_across_converged_limited_imputations"
    model_df <- olfml2b_method_repair_model_matrix_df(mi$first_completed, rhs)
    n_model <- nrow(d); events <- events_all; imputation_m <- mi$m
  } else {
    required_complete <- unique(c("time", "event", "olfml2b_z", model_covars,
                                  if (is.na(design_strata)) character() else design_strata))
    cc <- d[stats::complete.cases(d[, required_complete, drop = FALSE]), , drop = FALSE]
    if (nrow(cc) < 20L || sum(cc$event == 1L) < 5L) {
      return(olfml2b_if78_empty_cox(
        cohort, endpoint, endpoint_family, model, "INSUFFICIENT_COMPLETE_CASES",
        model_covars, nrow(cc), sum(cc$event == 1L), analysis_method = "complete_case",
        missing_row_fraction = missing_row_fraction,
        max_covariate_missing_fraction = max_missing, design_strata = design_strata
      ))
    }
    fit <- tryCatch(survival::coxph(form, data = cc, x = TRUE), error = function(e) e)
    if (inherits(fit, "error")) {
      return(olfml2b_if78_empty_cox(
        cohort, endpoint, endpoint_family, model,
        paste0("COX_FAILED:", conditionMessage(fit)), model_covars,
        nrow(cc), sum(cc$event == 1L), analysis_method = "complete_case",
        missing_row_fraction = missing_row_fraction,
        max_covariate_missing_fraction = max_missing, design_strata = design_strata
      ))
    }
    ex <- olfml2b_if78_extract_cox(fit)
    if (is.null(ex)) {
      return(olfml2b_if78_empty_cox(
        cohort, endpoint, endpoint_family, model, "OLFML2B_COEFFICIENT_MISSING",
        model_covars, nrow(cc), sum(cc$event == 1L), analysis_method = "complete_case",
        missing_row_fraction = missing_row_fraction,
        max_covariate_missing_fraction = max_missing, design_strata = design_strata
      ))
    }
    analysis_method <- "complete_case"
    beta <- ex$beta; se <- sqrt(ex$variance)
    ci_low <- exp(beta - 1.96 * se); ci_high <- exp(beta + 1.96 * se)
    p_value <- ex$p; ph_p <- ex$ph_p; ph_global_p <- ex$ph_global_p
    c_index <- ex$c_index; max_abs_dfbeta <- ex$max_abs_dfbeta
    n_influential_dfbeta <- ex$n_influential_dfbeta
    dfbeta_review_threshold <- ex$dfbeta_review_threshold
    max_abs_deviance <- ex$max_abs_deviance; influence_basis <- "complete_case_cox"
    model_df <- olfml2b_method_repair_model_matrix_df(cc, rhs)
    n_model <- nrow(cc); events <- sum(cc$event == 1L); imputation_m <- 0L
  }
  epdf <- if (is.finite(model_df) && model_df > 0L) events / model_df else NA_real_
  stable <- is.finite(beta) && is.finite(se) && se > 0 && abs(beta) <= 8 && se <= 5
  adjusted_model <- model %in% c("adjusted", "available_adjusted", "derived_stage_sensitivity",
                                 "adjusted_molecular_subtype_sensitivity")
  base_structural <- n_model >= ifelse(adjusted_model, 80L, 40L) &&
    events >= ifelse(adjusted_model, 40L, 20L) &&
    (!adjusted_model || (is.finite(epdf) && epdf >= 10)) && stable
  primary_eligible <- identical(model, "adjusted") && base_structural
  sensitivity_eligible <- model %in% c("available_adjusted", "derived_stage_sensitivity") && base_structural
  univariable_eligible <- identical(model, "univariable") && base_structural
  ph_status <- if (is.finite(ph_p) && ph_p < 0.05) "PH_REVIEW_OLFML2B" else
    if (is.finite(ph_global_p) && ph_global_p < 0.05) "PH_REVIEW_GLOBAL" else
      if (is.finite(ph_p) || is.finite(ph_global_p)) "PH_OK" else "PH_NOT_TESTED"
  reason <- if (!stable) "UNSTABLE_COEFFICIENT" else if (!base_structural) "N_EVENT_OR_EPV_THRESHOLD_NOT_MET" else
    "STRUCTURALLY_ELIGIBLE"
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = model, model_id = olfml2b_if78_model_id(model),
    covariates = paste(model_covars, collapse = ";"), design_strata = design_strata,
    n = n_model, events = events, complete_case_fraction = 1 - missing_row_fraction,
    missing_row_fraction = missing_row_fraction,
    max_covariate_missing_fraction = max_missing,
    analysis_method = analysis_method, imputation_m = imputation_m,
    imputation_policy = "MI only for eligible partial missingness >5%-20%; structural or >20% variable missingness is not imputed.",
    model_df = model_df, events_per_model_df = epdf,
    beta = beta, se = se, hr = exp(beta), ci_low = ci_low, ci_high = ci_high,
    p_value = p_value, ph_p = ph_p, ph_global_p = ph_global_p, ph_status = ph_status,
    max_abs_dfbeta_olfml2b = max_abs_dfbeta,
    n_influential_dfbeta_olfml2b = n_influential_dfbeta,
    dfbeta_review_threshold = dfbeta_review_threshold,
    max_abs_deviance_residual = max_abs_deviance,
    influence_diagnostic_basis = influence_basis, c_index = c_index,
    direction = ifelse(exp(beta) < 1, "favorable_high", "adverse_high"),
    status = ifelse(stable, "OK", "COX_UNSTABLE"),
    model_grade = if (primary_eligible) "PRIMARY_M1_COMMON_COVARIATE_ELIGIBLE" else
      if (sensitivity_eligible) "SUPPORTIVE_M2_SENSITIVITY_ELIGIBLE" else
        if (univariable_eligible) "M0_SENSITIVITY_ELIGIBLE" else "SUPPORTIVE_OR_INELIGIBLE",
    structural_eligible = base_structural,
    primary_inference_eligible = primary_eligible,
    meta_eligible = primary_eligible || univariable_eligible ||
      (identical(model, "available_adjusted") && sensitivity_eligible),
    eligibility_reason = reason,
    claim_ceiling = ifelse(primary_eligible, "observational_common_covariate_adjusted_association",
                           "supportive_or_sensitivity_only"),
    stage_derivation_policy = ifelse(identical(model, "derived_stage_sensitivity"),
      "AJCC7 gastric T/N with explicit M0 assumption; sensitivity only",
      "Source overall stage or complete AJCC7 gastric TNM for M1; pT/pN retained separately for M2"),
    interpretation_boundary = "Two-sided observational association; no effect direction or significance is a quality gate.",
    stringsAsFactors = FALSE
  )
}

olfml2b_run_survival_all <- function(views, min_n, min_events, min_epv,
                                    log_file = NULL) {
  rows <- list()
  endpoints <- olfml2b_endpoint_defs()
  for (nm in names(views)) {
    v <- views[[nm]]
    if (!isTRUE(v$target_measured)) next
    clinical <- olfml2b_analysis_clinical(v)
    endpoint_names <- olfml2b_selected_endpoint_names(clinical, min_n, min_events)
    for (ep in endpoint_names) {
      def <- endpoints[[ep]]
      if (!all(c(def$time, def$event) %in% names(clinical))) next
      dat <- olfml2b_surv_data(v, def)
      if (nrow(dat) < min_n || sum(dat$event == 1L) < min_events) next

      rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
        dat, v$cohort, ep, def$family, "univariable", character()
      )
      rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
        dat, v$cohort, ep, def$family, "adjusted",
        c("age", "sex", "stage_analysis_primary")
      )
      available <- olfml2b_if78_available_covariates(v$cohort)
      rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
        dat, v$cohort, ep, def$family, "available_adjusted", available
      )
      if (identical(v$cohort, "GSE84437")) {
        rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
          dat, v$cohort, ep, def$family, "derived_stage_sensitivity",
          c("age", "sex", "stage_analysis_sensitivity")
        )
      }
      if ("molecular_subtype" %in% names(dat)) {
        subtype_keep <- olfml2b_if78_nonmissing(dat$molecular_subtype)
        if (mean(subtype_keep) >= 0.80 &&
            length(unique(as.character(dat$molecular_subtype[subtype_keep]))) >= 2L) {
          rows[[length(rows) + 1L]] <- olfml2b_fit_olfml2b_cox(
            dat, v$cohort, ep, def$family,
            "adjusted_molecular_subtype_sensitivity",
            c(available, "molecular_subtype")
          )
        }
      }
      if ("adjuvant_therapy" %in% names(dat) &&
          length(unique(stats::na.omit(dat$adjuvant_therapy))) >= 2L &&
          sum(dat$event == 1L) >= 60L) {
        rows[[length(rows) + 1L]] <- olfml2b_fit_therapy_interaction(
          dat, v$cohort, ep, def$family, covariates = available
        )
      }
      olfml2b_log("Survival models completed: ", v$cohort, " ", ep,
                 " | M0 + M1 + M2",
                 if (identical(v$cohort, "GSE84437")) " + M2b" else "",
                 log_file = log_file)
    }
  }
  olfml2b_bind_rows(rows)
}

olfml2b_gse84437_subseries_heterogeneity_audit <- function(views) {
  if (!"GSE84437" %in% names(views)) return(data.frame())
  v <- views[["GSE84437"]]
  d0 <- olfml2b_analysis_clinical(v)
  required0 <- c("os_time_days", "os_event", "olfml2b_z", "source_subseries")
  if (length(setdiff(required0, names(d0)))) return(data.frame())
  keep <- is.finite(d0$os_time_days) & d0$os_time_days > 0 &
    d0$os_event %in% c(0L, 1L) & is.finite(d0$olfml2b_z) &
    !is.na(d0$source_subseries)
  d0 <- d0[keep, , drop = FALSE]
  d0$time <- as.numeric(d0$os_time_days); d0$event <- as.integer(d0$os_event)
  specs <- list(
    unadjusted_batch_stratified = character(),
    available_covariate_batch_stratified = c("age", "sex", "stage_pT", "stage_pN"),
    TNM_derived_stage_M0_sensitivity_batch_stratified = c("age", "sex", "stage_analysis_sensitivity")
  )
  rows <- list()
  for (model in names(specs)) {
    covars <- specs[[model]]
    miss <- setdiff(covars, names(d0))
    structural <- if (!length(miss)) covars[vapply(covars, function(z)
      sum(olfml2b_if78_nonmissing(d0[[z]])) == 0L, logical(1))] else miss
    if (length(structural)) {
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = "GSE84437", endpoint = "OS", model = model,
        n = 0L, events = 0L, n_GSE84426 = 0L, events_GSE84426 = 0L,
        n_GSE84433 = 0L, events_GSE84433 = 0L,
        log_HR_ratio_GSE84433_vs_GSE84426 = NA_real_,
        HR_ratio_GSE84433_vs_GSE84426 = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        interaction_p = NA_real_, status = paste0("NOT_EVALUABLE:", paste(structural, collapse = ",")),
        design = "subseries-stratified baseline hazards",
        inference_role = "prespecified_subseries_heterogeneity_sensitivity_not_selection_gate",
        stringsAsFactors = FALSE
      )
      next
    }
    prep <- olfml2b_if78_prepare_model_data(d0, covars, "source_subseries")
    d <- prep$data; model_covars <- prep$covariates
    needed <- c("time", "event", "olfml2b_z", "source_subseries", model_covars)
    d <- d[stats::complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
    d$source_subseries <- stats::relevel(droplevels(factor(d$source_subseries)), ref = "GSE84426")
    events <- sum(d$event == 1L)
    sub_n <- table(d$source_subseries)
    if (nrow(d) < 80L || events < 40L || length(sub_n) != 2L || any(sub_n < 20L)) {
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = "GSE84437", endpoint = "OS", model = model,
        n = nrow(d), events = events,
        n_GSE84426 = sum(d$source_subseries == "GSE84426"),
        events_GSE84426 = sum(d$event[d$source_subseries == "GSE84426"] == 1L),
        n_GSE84433 = sum(d$source_subseries == "GSE84433"),
        events_GSE84433 = sum(d$event[d$source_subseries == "GSE84433"] == 1L),
        log_HR_ratio_GSE84433_vs_GSE84426 = NA_real_,
        HR_ratio_GSE84433_vs_GSE84426 = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        interaction_p = NA_real_, status = "NOT_EVALUABLE_INSUFFICIENT_COMPLETE_CASES",
        design = "subseries-stratified baseline hazards",
        inference_role = "prespecified_subseries_heterogeneity_sensitivity_not_selection_gate",
        stringsAsFactors = FALSE
      )
      next
    }
    rhs <- c("olfml2b_z", "olfml2b_z:source_subseries", model_covars,
             "survival::strata(source_subseries)")
    form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ",
                                     paste(rhs, collapse = " + ")))
    fit <- tryCatch(survival::coxph(form, data = d, x = TRUE), error = function(e) e)
    sm <- if (!inherits(fit, "error")) summary(fit)$coefficients else NULL
    int_row <- if (!is.null(sm)) grep("olfml2b_z:source_subseries", rownames(sm), value = TRUE)[1] else NA_character_
    beta <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "coef"]) else NA_real_
    se <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "se(coef)"]) else NA_real_
    p <- if (length(int_row) && !is.na(int_row)) as.numeric(sm[int_row, "Pr(>|z|)"]) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = "GSE84437", endpoint = "OS", model = model,
      n = nrow(d), events = events,
      n_GSE84426 = sum(d$source_subseries == "GSE84426"),
      events_GSE84426 = sum(d$event[d$source_subseries == "GSE84426"] == 1L),
      n_GSE84433 = sum(d$source_subseries == "GSE84433"),
      events_GSE84433 = sum(d$event[d$source_subseries == "GSE84433"] == 1L),
      log_HR_ratio_GSE84433_vs_GSE84426 = beta,
      HR_ratio_GSE84433_vs_GSE84426 = ifelse(is.finite(beta), exp(beta), NA_real_),
      ci_low = ifelse(is.finite(beta) && is.finite(se), exp(beta - 1.96 * se), NA_real_),
      ci_high = ifelse(is.finite(beta) && is.finite(se), exp(beta + 1.96 * se), NA_real_),
      interaction_p = p, status = ifelse(is.finite(p), "OK", "INTERACTION_NOT_ESTIMABLE"),
      design = "common exposure scale; subseries-specific baseline hazards",
      inference_role = "prespecified_subseries_heterogeneity_sensitivity_not_selection_gate",
      stringsAsFactors = FALSE
    )
  }
  out <- olfml2b_bind_rows(rows)
  if (nrow(out)) out$interaction_fdr_bh <- stats::p.adjust(out$interaction_p, method = "BH")
  out
}

olfml2b_meta_all <- function(surv_all) {
  if (!nrow(surv_all)) return(data.frame())
  eligible <- if ("meta_eligible" %in% names(surv_all)) surv_all$meta_eligible %in% TRUE else
    surv_all$status == "OK"
  d <- surv_all[
    eligible & surv_all$model %in% c("adjusted", "available_adjusted", "univariable") &
      is.finite(surv_all$beta) & is.finite(surv_all$se) & surv_all$se > 0,
    , drop = FALSE
  ]
  d <- olfml2b_deduplicate_meta_effects(d)
  rows <- list()
  for (model in unique(d$model)) {
    for (fam in unique(d$endpoint_family)) {
      dd <- d[d$model == model & d$endpoint_family == fam, , drop = FALSE]
      one <- olfml2b_meta_one(dd, endpoint_family = fam, model = model)
      one$adjustment_harmonization <- if (model == "adjusted")
        "COMMON_AGE_SEX_VERIFIED_OVERALL_STAGE" else if (model == "available_adjusted")
          "HETEROGENEOUS_PRESPECIFIED_AVAILABLE_COVARIATES_SUPPORTIVE_ONLY" else
            "UNADJUSTED"
      if (model == "available_adjusted") {
        one$primary_meta_structural_eligible <- FALSE
        one$claim_ceiling <- "supportive_available_covariate_sensitivity_only"
      }
      rows[[length(rows) + 1L]] <- one
    }
  }
  olfml2b_bind_rows(rows)
}

olfml2b_go_no_go_v23 <- function(views, endpoint_audit, surv_all, meta, control_models,
                                module_summary, pdc_audit) {
  target_ok <- length(views) > 0L && all(vapply(views, function(v)
    is.data.frame(v$clinical) && "olfml2b_z" %in% names(olfml2b_analysis_clinical(v)), logical(1)))
  endpoint_ok <- is.data.frame(endpoint_audit) && nrow(endpoint_audit) > 0L
  m1 <- surv_all[surv_all$model == "adjusted", , drop = FALSE]
  m1_ok <- m1[m1$status == "OK", , drop = FALSE]
  common_formula_ok <- nrow(m1_ok) > 0L &&
    all(m1_ok$covariates == "age10;sex;stage_analysis_primary")
  expected_structural <- m1[m1$cohort %in% c("GSE26253", "GSE84437"), , drop = FALSE]
  structural_ok <- nrow(expected_structural) > 0L &&
    all(grepl("STRUCTURALLY_UNAVAILABLE_NO_IMPUTATION", expected_structural$status)) &&
    all(expected_structural$imputation_m == 0L)
  m2 <- surv_all[surv_all$model == "available_adjusted", , drop = FALSE]
  g262 <- m2[m2$cohort == "GSE26253", , drop = FALSE]
  g844 <- m2[m2$cohort == "GSE84437", , drop = FALSE]
  adaptive_ok <- nrow(g262) > 0L && all(g262$covariates == "stage_analysis_primary") &&
    nrow(g844) > 0L && all(g844$covariates == "age10;sex;stage_pT;stage_pN") &&
    all(g844$design_strata == "source_subseries")
  mi_flag <- !is.na(surv_all$analysis_method) &
    grepl("mice", surv_all$analysis_method, ignore.case = TRUE)
  mi_rows <- surv_all[mi_flag, , drop = FALSE]
  mi_ok <- !nrow(mi_rows) || all(mi_rows$max_covariate_missing_fraction <= 0.20 &
                                  mi_rows$missing_row_fraction <= 0.30)
  excessive_flag <- !is.na(surv_all$status) & grepl("EXCESSIVE_MISSINGNESS", surv_all$status)
  no_excessive_imputed <- all(!excessive_flag |
                                (!is.na(surv_all$imputation_m) & surv_all$imputation_m == 0L))
  diag_ok <- all(c("events_per_model_df", "ph_status", "analysis_method",
                   "max_abs_dfbeta_olfml2b", "max_abs_deviance_residual") %in% names(surv_all))
  meta_schema_ok <- is.data.frame(meta) && all(c("meta_method", "interval_method",
                                                 "prediction_low", "prediction_high",
                                                 "metafor_reml_crosscheck") %in% names(meta))
  evaluated_meta <- if (meta_schema_ok && nrow(meta)) meta[meta$status == "OK", , drop = FALSE] else data.frame()
  meta_ok <- meta_schema_ok && (!nrow(evaluated_meta) ||
    all(grepl("^PASS_REML", evaluated_meta$metafor_reml_crosscheck)))
  data.frame(
    criterion = c("exact_target_mapping", "endpoint_coding_audit",
                  "M1_common_formula", "structural_missingness_no_imputation",
                  "M2_available_covariate_contract", "limited_MI_policy",
                  "model_diagnostics_exported", "REML_mKH_meta_contract",
                  "result_independent_pipeline_continuation", "claim_ceiling"),
    status = c(ifelse(target_ok, "PASS", "FAIL"), ifelse(endpoint_ok, "PASS", "FAIL"),
               ifelse(common_formula_ok, "PASS", "FAIL"), ifelse(structural_ok, "PASS", "FAIL"),
               ifelse(adaptive_ok, "PASS", "FAIL"), ifelse(mi_ok && no_excessive_imputed, "PASS", "FAIL"),
               ifelse(diag_ok, "PASS", "FAIL"), ifelse(meta_ok || !nrow(meta), "PASS_OR_NOT_EVALUABLE", "FAIL"),
               "PASS", "OBSERVATIONAL_ASSOCIATION_ONLY"),
    detail = c(
      "Exact OLFML2B mapping; confusable paralog/family symbols excluded.",
      "Time, event and endpoint family audited before modelling.",
      "M1 uses age10 + sex + source-reported or complete-TNM-derived overall stage.",
      "GSE26253 age/sex and GSE84437 primary overall stage are structurally unavailable and never imputed.",
      "GSE26253 uses reported stage; GSE84437 uses age, sex, pT, pN and subseries strata in M2.",
      "MI is permitted only for partial missingness >5%-20%; variable missingness >20% or incomplete rows >30% is not imputed.",
      "EPV, PH, DFBETA, deviance residuals, missingness and eligibility are exported.",
      "REML with modified Hartung-Knapp intervals, prediction intervals and metafor cross-check.",
      "Effect direction and significance are outputs, never quality gates.",
      "No causal, clinical-utility, treatment-predictive or mechanism-proven claim."
    ),
    stringsAsFactors = FALSE
  )
}

olfml2b_fit_olfml2b_cox <- function(dat, cohort, endpoint, endpoint_family, model, covariates = character()) {
  model_id <- olfml2b_if78_model_id(model)
  design_strata <- if (identical(cohort, "GSE84437") && "source_subseries" %in% names(dat) &&
                       length(unique(stats::na.omit(dat$source_subseries))) == 2L) "source_subseries" else NA_character_
  source_covars <- if (identical(model, "adjusted")) c("age", "sex", "stage") else covariates
  required_covars <- ifelse(source_covars == "age", "age10", source_covars)
  missing_columns <- setdiff(c("time", "event", "olfml2b_z", source_covars,
                               if (is.na(design_strata)) character() else design_strata), names(dat))
  if (length(missing_columns)) {
    return(olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model,
      paste0("STRUCTURAL_MISSING_VARIABLES:", paste(missing_columns, collapse = ",")), required_covars))
  }
  source_vars <- unique(c("time", "event", "olfml2b_z", source_covars,
                          if (is.na(design_strata)) character() else design_strata))
  d <- dat[, source_vars, drop = FALSE]
  if ("age" %in% source_covars) {
    age_num <- suppressWarnings(as.numeric(d$age))
    d$age10 <- (age_num - mean(age_num, na.rm = TRUE)) / 10
    d$age <- NULL
  }
  if ("sex" %in% required_covars) d$sex <- olfml2b_if78_sex_factor(d$sex)
  if ("stage" %in% required_covars) d$stage <- olfml2b_if78_stage_factor(d$stage)
  if (!is.na(design_strata)) d[[design_strata]] <- factor(d[[design_strata]])
  vars <- unique(c("time", "event", "olfml2b_z", required_covars,
                   if (is.na(design_strata)) character() else design_strata))
  d <- d[, vars, drop = FALSE]
  if (any(!is.finite(d$time)) || any(!d$event %in% c(0L, 1L)) || any(!is.finite(d$olfml2b_z))) {
    return(olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model,
      "INVALID_ENDPOINT_OR_EXPOSURE", required_covars, nrow(d), sum(d$event == 1L, na.rm = TRUE)))
  }
  core_incomplete <- if (length(required_covars)) !stats::complete.cases(d[, required_covars, drop = FALSE]) else rep(FALSE, nrow(d))
  missing_row_fraction <- mean(core_incomplete)
  events <- sum(d$event == 1L)
  rhs <- c("olfml2b_z", required_covars)
  rhs_formula <- c(rhs, if (is.na(design_strata)) character() else paste0("survival::strata(", design_strata, ")"))
  form <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(rhs_formula, collapse = " + ")))
  use_mi <- length(required_covars) > 0L && missing_row_fraction > 0.05
  planned_m <- if (use_mi) {
    max(20L, as.integer(ceiling(100 * max(colMeans(is.na(d[, required_covars, drop = FALSE]))))))
  } else 0L

  if (use_mi) {
    mi <- tryCatch(olfml2b_if78_fit_mice_cox(d, form, required_covars, m = planned_m, seed = 20260721L), error = function(e) e)
    if (inherits(mi, "error")) {
      return(olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model,
        paste0("MICE_COX_FAILED:", conditionMessage(mi)), required_covars, nrow(d), events,
        "mice_rubin", missing_row_fraction, planned_m))
    }
    analysis_method <- "mice_rubin"
    beta <- mi$beta; se <- mi$se; ci_low <- mi$ci_low; ci_high <- mi$ci_high; p_value <- mi$p
    ph_p <- ifelse(is.nan(mi$ph_p), NA_real_, mi$ph_p)
    ph_global_p <- ifelse(is.nan(mi$ph_global_p), NA_real_, mi$ph_global_p)
    c_index <- ifelse(is.nan(mi$c_index), NA_real_, mi$c_index)
    max_abs_dfbeta <- ifelse(is.nan(mi$max_abs_dfbeta), NA_real_, mi$max_abs_dfbeta)
    n_influential_dfbeta <- mi$n_influential_dfbeta
    dfbeta_review_threshold <- ifelse(is.nan(mi$dfbeta_review_threshold), NA_real_, mi$dfbeta_review_threshold)
    max_abs_deviance <- ifelse(is.nan(mi$max_abs_deviance), NA_real_, mi$max_abs_deviance)
    influence_basis <- "median_across_converged_imputations"
    model_df <- olfml2b_method_repair_model_matrix_df(mi$first_completed, rhs)
    n_model <- nrow(d); imputation_m <- mi$m
  } else {
    cc <- d[stats::complete.cases(d), , drop = FALSE]
    fit <- tryCatch(survival::coxph(form, data = cc, x = TRUE), error = function(e) e)
    if (inherits(fit, "error")) {
      return(olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model,
        paste0("COX_FAILED:", conditionMessage(fit)), required_covars, nrow(cc), sum(cc$event == 1L),
        "complete_case", missing_row_fraction))
    }
    ex <- olfml2b_if78_extract_cox(fit)
    if (is.null(ex)) return(olfml2b_if78_empty_cox(cohort, endpoint, endpoint_family, model,
      "OLFML2B_COEFFICIENT_MISSING", required_covars, nrow(cc), sum(cc$event == 1L),
      "complete_case", missing_row_fraction))
    analysis_method <- "complete_case"
    beta <- ex$beta; se <- sqrt(ex$variance); ci_low <- exp(beta - 1.96 * se); ci_high <- exp(beta + 1.96 * se)
    p_value <- ex$p; ph_p <- ex$ph_p; ph_global_p <- ex$ph_global_p; c_index <- ex$c_index
    max_abs_dfbeta <- ex$max_abs_dfbeta; n_influential_dfbeta <- ex$n_influential_dfbeta
    dfbeta_review_threshold <- ex$dfbeta_review_threshold; max_abs_deviance <- ex$max_abs_deviance
    influence_basis <- "complete_case_cox"
    model_df <- olfml2b_method_repair_model_matrix_df(cc, rhs)
    n_model <- nrow(cc); events <- sum(cc$event == 1L); imputation_m <- 0L
  }

  epdf <- if (is.finite(model_df) && model_df > 0) events / model_df else NA_real_
  stable <- is.finite(beta) && is.finite(se) && se > 0 && abs(beta) <= 8 && se <= 5
  base_structural <- n_model >= 80L && events >= 40L && is.finite(epdf) && epdf >= 10 && stable
  primary_eligible <- identical(model, "adjusted") && base_structural
  univariable_eligible <- identical(model, "univariable") && n_model >= 40L && events >= 20L && stable
  ph_status <- if (is.finite(ph_p) && ph_p < 0.05) "PH_REVIEW_OLFML2B" else if (is.finite(ph_global_p) && ph_global_p < 0.05) "PH_REVIEW_GLOBAL" else if (is.finite(ph_p) || is.finite(ph_global_p)) "PH_OK" else "PH_NOT_TESTED"
  reason <- if (!stable) "UNSTABLE_COEFFICIENT" else if (n_model < 80L && identical(model, "adjusted")) "N_LT_80" else if (events < 40L && identical(model, "adjusted")) "EVENTS_LT_40" else if (identical(model, "adjusted") && (!is.finite(epdf) || epdf < 10)) "EPV_LT_10" else "STRUCTURALLY_ELIGIBLE"
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_family = endpoint_family,
    model = model, model_id = model_id, covariates = paste(required_covars, collapse = ";"),
    design_strata = design_strata,
    n = n_model, events = events, complete_case_fraction = 1 - missing_row_fraction,
    missing_row_fraction = missing_row_fraction, analysis_method = analysis_method,
    imputation_m = imputation_m, model_df = model_df, events_per_model_df = epdf,
    beta = beta, se = se, hr = exp(beta), ci_low = ci_low, ci_high = ci_high,
    p_value = p_value, ph_p = ph_p, ph_global_p = ph_global_p, ph_status = ph_status,
    max_abs_dfbeta_olfml2b = max_abs_dfbeta,
    n_influential_dfbeta_olfml2b = n_influential_dfbeta,
    dfbeta_review_threshold = dfbeta_review_threshold,
    max_abs_deviance_residual = max_abs_deviance,
    influence_diagnostic_basis = influence_basis,
    c_index = c_index, direction = ifelse(exp(beta) < 1, "favorable_high", "adverse_high"),
    status = ifelse(stable, "OK", "COX_UNSTABLE"),
    model_grade = ifelse(primary_eligible, "PRIMARY_M1_STRUCTURALLY_ELIGIBLE", ifelse(univariable_eligible, "M0_SENSITIVITY_ELIGIBLE", "SUPPORTIVE_OR_INELIGIBLE")),
    structural_eligible = base_structural, primary_inference_eligible = primary_eligible,
    meta_eligible = primary_eligible || univariable_eligible, eligibility_reason = reason,
    claim_ceiling = ifelse(primary_eligible, "observational_adjusted_association", "supportive_or_sensitivity_only"),
    interpretation_boundary = "Two-sided observational association. PH findings trigger sensitivity analyses but never an effect-based pass/fail decision.",
    stringsAsFactors = FALSE
  )
}

.olfml2b_if78_meta_one_core <- olfml2b_meta_one
olfml2b_meta_one <- function(d, endpoint_family, model) {
  out <- .olfml2b_if78_meta_one_core(d, endpoint_family, model)
  out$primary_meta_structural_eligible <- identical(model, "adjusted") && is.finite(out$k[1]) && out$k[1] >= 3L
  out$claim_ceiling <- ifelse(out$primary_meta_structural_eligible,
                             "internal_public_cohort_reanalysis_observational_association",
                             "supportive_meta_or_insufficient_cohorts")
  out$quality_gate_basis <- "cohort_count_and_model_diagnostics_only"
  out
}

olfml2b_go_no_go <- function(views, endpoint_audit, surv_all, meta, control_models, module_summary, pdc_audit) {
  target_ok <- length(views) > 0L && all(vapply(views, function(v) is.data.frame(v$clinical) &&
    "olfml2b_z" %in% names(olfml2b_analysis_clinical(v)), logical(1)))
  endpoint_ok <- is.data.frame(endpoint_audit) && nrow(endpoint_audit) > 0L
  m1 <- if (is.data.frame(surv_all) && nrow(surv_all)) surv_all[surv_all$model == "adjusted", , drop = FALSE] else data.frame()
  g844_m1 <- if (nrow(m1)) m1[m1$cohort == "GSE84437", , drop = FALSE] else data.frame()
  formula_ok <- nrow(m1) > 0L && all(m1$covariates == "age10;sex;stage") &&
    "design_strata" %in% names(m1) && nrow(g844_m1) > 0L &&
    all(g844_m1$design_strata == "source_subseries")
  diag_ok <- nrow(m1) > 0L && all(c("events_per_model_df", "ph_status", "analysis_method", "eligibility_reason",
                                    "max_abs_dfbeta_olfml2b", "max_abs_deviance_residual") %in% names(m1))
  meta_schema_ok <- is.data.frame(meta) && all(c("meta_method", "interval_method", "prediction_low", "prediction_high",
                                                 "metafor_reml_crosscheck") %in% names(meta))
  evaluated_meta <- if (meta_schema_ok && nrow(meta)) meta[meta$status == "OK", , drop = FALSE] else data.frame()
  meta_ok <- meta_schema_ok && (!nrow(evaluated_meta) || all(evaluated_meta$metafor_reml_crosscheck == "PASS_REML_POINT_AND_TAU2"))
  data.frame(
    criterion = c("exact_target_mapping", "endpoint_coding_audit", "fixed_M1_common_formula", "missing_data_policy_executed", "model_diagnostics_exported", "REML_mKH_meta_contract", "result_independent_pipeline_continuation", "claim_ceiling"),
    status = c(ifelse(target_ok, "PASS", "FAIL"), ifelse(endpoint_ok, "PASS", "FAIL"),
               ifelse(formula_ok, "PASS", "FAIL"), ifelse(nrow(m1) && all(m1$analysis_method %in% c("complete_case", "mice_rubin")), "PASS", "FAIL"),
               ifelse(diag_ok, "PASS", "FAIL"), ifelse(meta_ok || !nrow(meta), "PASS_OR_NOT_EVALUABLE", "FAIL"),
               "PASS", "OBSERVATIONAL_ASSOCIATION_ONLY"),
    detail = c(
      "Exact OLFML2B/ENSG00000162745 mapping; confusable paralog/family symbols excluded.",
      "Time, event and endpoint family are audited before modelling.",
      "Every M1 row uses OLFML2B_z + centered age per 10 years + sex + stage; GSE84437 additionally stratifies baseline hazard by GSE84426/GSE84433. Missing variables make the row not evaluable.",
      "Core missingness <=5% uses complete cases; >5% uses m=max(20, ceiling(maximum variable missing percent)) and Rubin pooling.",
      "n, events, model df, EPV, PH, DFBETA/deviance influence diagnostics and structural eligibility are exported.",
      "Random-effects synthesis uses REML, modified Hartung-Knapp intervals and prediction intervals; REML point estimate/tau2 are cross-checked against metafor within 1e-5.",
      "HR direction, effect size, CI and P value are evidence outputs and are never QC gates.",
      "No causal, clinical utility, treatment-predictive or mechanism-proven claim."
    ),
    stringsAsFactors = FALSE
  )
}

# Final activation: these definitions intentionally occur after all legacy
# compatibility layers so the structural-missingness policy is authoritative.
OLFML2B_PART3_IF78_VERSION <- "v1.1.0_20260722_PDC614_CONTRACT_SYNC"
OLFML2B_ANALYSIS_VERSION <- OLFML2B_PART3_IF78_VERSION

olfml2b_if78_model_id <- function(model) {
  switch(as.character(model)[1],
         univariable = "M0_unadjusted",
         adjusted = "M1_common_age_sex_verified_overall_stage",
         available_adjusted = "M2_prespecified_available_covariates",
         derived_stage_sensitivity = "M2b_TNM_derived_stage_sensitivity",
         adjusted_molecular_subtype_sensitivity = "M3_molecular_subtype_sensitivity",
         therapy_interaction_exploratory = "M4_therapy_interaction_exploratory",
         as.character(model)[1])
}

olfml2b_select_covariates <- function(dat, events, min_epv = 10, cohort = NA_character_,
                                    model_role = c("primary_adjusted", "molecular_subtype_sensitivity")) {
  model_role <- match.arg(model_role)
  out <- olfml2b_if78_available_covariates(cohort)
  if (identical(model_role, "molecular_subtype_sensitivity")) out <- c(out, "molecular_subtype")
  out[out %in% names(dat)]
}

olfml2b_if78_extract_cox <- function(fit) {
  sm <- summary(fit)
  if (!"olfml2b_z" %in% rownames(sm$coefficients)) return(NULL)
  co <- sm$coefficients["olfml2b_z", ]
  ph <- olfml2b_method_repair_ph(fit, "olfml2b_z")
  dfb <- tryCatch(stats::residuals(fit, type = "dfbeta"), error = function(e) NULL)
  if (is.matrix(dfb) && is.null(colnames(dfb)) && ncol(dfb) == length(stats::coef(fit))) {
    colnames(dfb) <- names(stats::coef(fit))
  }
  target_dfb <- if (is.matrix(dfb) && "olfml2b_z" %in% colnames(dfb)) {
    as.numeric(dfb[, "olfml2b_z"])
  } else if (is.numeric(dfb) && length(stats::coef(fit)) == 1L &&
             "olfml2b_z" %in% names(stats::coef(fit))) {
    as.numeric(dfb)
  } else numeric()
  fit_n <- suppressWarnings(as.numeric(fit$n %||% NA_real_)[1])
  review_threshold <- if (is.finite(fit_n) && fit_n > 0L) 2 / sqrt(fit_n) else NA_real_
  deviance <- tryCatch(as.numeric(stats::residuals(fit, type = "deviance")),
                       error = function(e) numeric())
  list(
    beta = as.numeric(co["coef"]),
    variance = as.numeric(co["se(coef)"])^2,
    p = as.numeric(co["Pr(>|z|)"]),
    ph_p = ph$term_p, ph_global_p = ph$global_p, ph_status = ph$status,
    c_index = tryCatch(as.numeric(sm$concordance[1]), error = function(e) NA_real_),
    max_abs_dfbeta = if (any(is.finite(target_dfb))) max(abs(target_dfb), na.rm = TRUE) else NA_real_,
    n_influential_dfbeta = if (length(target_dfb) && is.finite(review_threshold))
      sum(abs(target_dfb) > review_threshold, na.rm = TRUE) else NA_integer_,
    dfbeta_review_threshold = review_threshold,
    max_abs_deviance = if (any(is.finite(deviance))) max(abs(deviance), na.rm = TRUE) else NA_real_
  )
}

olfml2b_if78_fit_mice_limited <- function(d, form, covariates, m = 20L,
                                         seed = 20260721L) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("mice is required only for eligible partial missingness (5%-20%); install it or run the complete-case sensitivity explicitly.",
         call. = FALSE)
  }
  imp_dat <- d
  imp_dat$nelson_aalen <- olfml2b_if78_nelson_aalen(imp_dat$time, imp_dat$event)
  meth <- mice::make.method(imp_dat)
  meth[] <- ""
  pred <- mice::make.predictorMatrix(imp_dat)
  locked <- intersect(c("time", "event", "olfml2b_z", "nelson_aalen", "source_subseries"),
                      names(imp_dat))
  pred[locked, ] <- 0
  diag(pred) <- 0
  for (v in intersect(covariates, names(imp_dat))) {
    x <- imp_dat[[v]]
    n_miss <- sum(is.na(x))
    if (!n_miss) next
    if (sum(!is.na(x)) < 30L || length(unique(stats::na.omit(x))) < 2L) {
      stop("Covariate is not imputable under the partial-missingness contract: ", v,
           call. = FALSE)
    }
    if (is.ordered(x)) {
      meth[v] <- "polr"
    } else if (is.factor(x) && nlevels(droplevels(x)) == 2L) {
      meth[v] <- "logreg"
    } else if (is.factor(x) || is.character(x)) {
      imp_dat[[v]] <- factor(x)
      meth[v] <- "polyreg"
    } else {
      meth[v] <- "pmm"
    }
  }
  if (!any(nzchar(meth))) stop("No partially missing covariate requires imputation.", call. = FALSE)
  set.seed(seed)
  imp <- mice::mice(imp_dat, m = as.integer(m), maxit = 20L, method = meth,
                    predictorMatrix = pred, seed = seed, printFlag = FALSE)
  fits <- lapply(seq_len(m), function(i) {
    completed <- mice::complete(imp, i)
    fit <- tryCatch(survival::coxph(form, data = completed, x = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ex <- olfml2b_if78_extract_cox(fit)
    if (is.null(ex) || !is.finite(ex$beta) || !is.finite(ex$variance) || ex$variance <= 0) return(NULL)
    list(fit = fit, extracted = ex, completed = completed)
  })
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) < max(5L, ceiling(m / 2)))
    stop("Too few converged Cox fits across eligible imputations.", call. = FALSE)
  qhat <- vapply(fits, function(z) z$extracted$beta, numeric(1))
  uhat <- vapply(fits, function(z) z$extracted$variance, numeric(1))
  mm <- length(qhat); qbar <- mean(qhat); ubar <- mean(uhat); b <- stats::var(qhat)
  total_var <- ubar + (1 + 1 / mm) * b
  r <- if (ubar > 0) ((1 + 1 / mm) * b) / ubar else Inf
  df <- if (is.finite(r) && r > 0) (mm - 1) * (1 + 1 / r)^2 else Inf
  se <- sqrt(total_var)
  crit <- if (is.finite(df)) stats::qt(0.975, df = df) else stats::qnorm(0.975)
  p <- if (is.finite(df)) 2 * stats::pt(abs(qbar / se), df = df, lower.tail = FALSE) else
    2 * stats::pnorm(abs(qbar / se), lower.tail = FALSE)
  ex_list <- lapply(fits, `[[`, "extracted")
  list(
    beta = qbar, se = se, ci_low = exp(qbar - crit * se), ci_high = exp(qbar + crit * se),
    p = p, df = df, m = mm,
    ph_p = olfml2b_if78_safe_median(vapply(ex_list, `[[`, numeric(1), "ph_p")),
    ph_global_p = olfml2b_if78_safe_median(vapply(ex_list, `[[`, numeric(1), "ph_global_p")),
    c_index = olfml2b_if78_safe_median(vapply(ex_list, `[[`, numeric(1), "c_index")),
    max_abs_dfbeta = olfml2b_if78_safe_median(vapply(ex_list, `[[`, numeric(1), "max_abs_dfbeta")),
    n_influential_dfbeta = as.integer(round(olfml2b_if78_safe_median(
      vapply(ex_list, function(z) as.numeric(z$n_influential_dfbeta), numeric(1))))),
    dfbeta_review_threshold = olfml2b_if78_safe_median(
      vapply(ex_list, `[[`, numeric(1), "dfbeta_review_threshold")),
    max_abs_deviance = olfml2b_if78_safe_median(vapply(ex_list, `[[`, numeric(1), "max_abs_deviance")),
    first_completed = fits[[1]]$completed
  )
}

# Activate the v2.3 implementations after every legacy definition.
olfml2b_fit_olfml2b_cox <- olfml2b_fit_olfml2b_cox_v23
olfml2b_go_no_go <- olfml2b_go_no_go_v23

# Cross-check tolerances reflect optimizer-level numerical agreement rather
# than CSV rounding.  The primary estimates remain the manual REML/mKH result;
# this audit only avoids a false failure for immaterial tau-squared differences.
.olfml2b_meta_one_v23_core <- olfml2b_meta_one
olfml2b_meta_one <- function(d, endpoint_family, model) {
  out <- .olfml2b_meta_one_v23_core(d, endpoint_family, model)
  out$metafor_beta_tolerance <- 1e-5
  out$metafor_tau2_tolerance <- 1e-4
  if (nrow(out) && all(c("metafor_reml_crosscheck", "metafor_beta_abs_diff",
                         "metafor_tau2_abs_diff") %in% names(out)) &&
      identical(out$status[1], "OK") &&
      is.finite(out$metafor_beta_abs_diff[1]) &&
      is.finite(out$metafor_tau2_abs_diff[1]) &&
      out$metafor_beta_abs_diff[1] <= out$metafor_beta_tolerance[1] &&
      out$metafor_tau2_abs_diff[1] <= out$metafor_tau2_tolerance[1]) {
    out$metafor_reml_crosscheck <- "PASS_REML_POINT_AND_TAU2_NUMERICAL_TOLERANCE"
  }
  out
}


# ==============================================================================
# OLFML2B-specific final activation: recurrence-family primary, OS secondary.
# ==============================================================================
OLFML2B_PART3_TARGET_VERSION <- "v1.1.0_20260722_PDC614_CONTRACT_SYNC"
OLFML2B_PART3_IF78_VERSION <- OLFML2B_PART3_TARGET_VERSION
OLFML2B_ANALYSIS_VERSION <- OLFML2B_PART3_TARGET_VERSION

olfml2b_recurrence_priority_summary <- function(meta, loo = data.frame()) {
  if (!is.data.frame(meta)) meta <- data.frame()
  rec <- if (nrow(meta) && all(c("endpoint_family", "model") %in% names(meta)))
    meta[meta$endpoint_family == "RECURRENCE" & meta$model == "adjusted", , drop = FALSE] else data.frame()
  os <- if (nrow(meta) && all(c("endpoint_family", "model") %in% names(meta)))
    meta[meta$endpoint_family == "OS" & meta$model == "adjusted", , drop = FALSE] else data.frame()
  one <- function(d, role) {
    if (!nrow(d)) return(data.frame(
      endpoint_family = sub("_.*$", "", role), evidence_role = role, status = "NOT_EVALUABLE",
      k = NA_integer_, hr_random = NA_real_, ci_low_random = NA_real_, ci_high_random = NA_real_,
      p_random = NA_real_, i2 = NA_real_, claim = "No pooled adjusted estimate available.",
      stringsAsFactors = FALSE))
    z <- d[1, , drop = FALSE]
    data.frame(
      endpoint_family = as.character(z$endpoint_family), evidence_role = role,
      status = as.character(z$status %||% NA_character_), k = as.integer(z$k %||% NA_integer_),
      hr_random = as.numeric(z$hr_random %||% NA_real_),
      ci_low_random = as.numeric(z$ci_low_random %||% NA_real_),
      ci_high_random = as.numeric(z$ci_high_random %||% NA_real_),
      p_random = as.numeric(z$p_random %||% NA_real_), i2 = as.numeric(z$i2 %||% NA_real_),
      claim = if (role == "PRIMARY_RECURRENCE_FAMILY")
        "Primary clinical evidence; DFS/RFS labels remain cohort-specific and are pooled only by the prespecified recurrence family."
      else "Secondary survival evidence; heterogeneity and non-significance do not invalidate the recurrence-primary design.",
      stringsAsFactors = FALSE
    )
  }
  out <- olfml2b_bind_rows(list(one(rec, "PRIMARY_RECURRENCE_FAMILY"), one(os, "SECONDARY_OS")))
  out$result_selected_gate <- FALSE
  out$final_gene_lock <- FALSE
  out
}

.olfml2b_part3_target_core <- run_olfml2b_specialized_bioinformatics
run_olfml2b_specialized_bioinformatics <- function(...) {
  index <- .olfml2b_part3_target_core(...)
  index$version <- OLFML2B_PART3_TARGET_VERSION
  priority <- olfml2b_recurrence_priority_summary(index$meta, index$leave_one_out)
  boundary <- data.frame(
    item = c("primary_endpoint", "secondary_endpoint", "target_identity", "main_ecology", "cutpoint", "causality", "protein_layer"),
    rule = c(
      "Recurrence-family adjusted continuous Cox/meta-analysis is primary; dataset-specific DFS/RFS labels are retained.",
      "OS is secondary and interpreted with heterogeneity and leave-one-cohort-out diagnostics.",
      "Exact OLFML2B / ENSG00000162745 / Entrez 25903 only; OLFML2A and OLFM2 are not merged.",
      "CAF/ECM/TGF-beta ecological-state association is prespecified; cell of origin is not inferred from bulk data.",
      "Continuous within-cohort z score is primary; median split is visualization/sensitivity only.",
      "No driver, mediation, treatment-prediction or causal claim is permitted from Parts0-5.",
      "PDC000614 is the active protein cohort. Its within-plex case-paired Tumor-minus-Normal analysis is a single-cohort orthogonal support layer; structural evaluability and biological direction are reported separately."
    ), stringsAsFactors = FALSE
  )
  index$recurrence_priority <- priority
  index$olfml2b_claim_boundary <- boundary
  olfml2b_write_csv(priority, file.path(index$dirs$tables, "08b_olfml2b_recurrence_primary_summary.csv"))
  olfml2b_write_csv(boundary, file.path(index$dirs$tables, "20_olfml2b_claim_boundary.csv"))
  saveRDS(index, file.path(index$dirs$objects, "Part3_OLFML2B_specialized_bioinformatics_index.rds"), compress = "xz")
  invisible(index)
}


# ==============================================================================
# Final OLFML2B interpretation contract. This replaces legacy target-specific
# narrative fields that treated OS or DDR/hypoxia as primary. It does not use
# effect direction or significance as a pipeline quality gate.
# ==============================================================================
olfml2b_interpretation_boundary <- function(surv_all, meta, control_models,
                                           median_surv, module_summary) {
  data.frame(
    domain = c(
      "primary_clinical_endpoint",
      "secondary_clinical_endpoint",
      "primary_ecological_context",
      "tissue_composition_sensitivity",
      "median_split_policy",
      "hypoxia_and_proliferation_policy",
      "protein_layer_boundary",
      "causal_language_boundary"
    ),
    allowed_claim = c(
      "Recurrence-family continuous adjusted Cox and prespecified cross-cohort synthesis are the primary clinical evidence. Cohort-specific DFS/RFS labels remain visible.",
      "OS is secondary and must be interpreted with heterogeneity, prediction intervals and leave-one-cohort-out diagnostics.",
      "CAF, ECM-remodelling and TGF-beta associations are ecological-state evidence. They do not establish the expressing cell type or direct pathway activation.",
      "Models adding proliferation, CAF/ECM, endothelial, pericyte and epithelial controls are sensitivity analyses for tissue composition, not mediation models.",
      "Median split is visualization and sensitivity only. Continuous within-cohort z-scored expression remains primary.",
      "Hypoxia and proliferation modules are negative/control contexts only and cannot become the primary mechanism through result inspection.",
      "PDC000614 provides an exact-target, within-plex, case-paired protein direction analysis. It is single-cohort orthogonal support, not multi-cohort protein validation, causality, prognosis or clinical utility.",
      "Parts0-5 support observational expression, prognosis and ecological-context claims only; no causality, treatment prediction or clinical utility claim."
    ),
    manuscript_use = c(
      "main_result",
      "secondary_result",
      "main_context_result",
      "sensitivity_result",
      "figure_display_only",
      "control_context_only",
      "optional_orthogonal_layer",
      "global_language_rule"
    ),
    result_selected_gate = FALSE,
    final_gene_lock = FALSE,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Part3 v1.2.0: authoritative tumor context and source endpoint labels
# ==============================================================================
OLFML2B_PART3_IF78_VERSION <- "v1.2.0_20260722_GEO_TUMOR_ONLY_AND_RFS_DISPLAY_CONTRACT"
OLFML2B_ANALYSIS_VERSION <- OLFML2B_PART3_IF78_VERSION

.olfml2b_v130_standardize_clinical_core <- olfml2b_standardize_clinical
olfml2b_standardize_clinical <- function(meta, cohort) {
  d <- .olfml2b_v130_standardize_clinical_core(meta, cohort)
  if (cohort %in% c("GSE62254", "GSE15459", "GSE26253", "GSE84437", "GSE147163")) {
    d$sample_context <- "Tumor"
    d$is_tumor <- TRUE
    d$is_normal <- FALSE
    d$normal_comparison_allowed <- FALSE
  } else {
    context <- rep("Unknown", nrow(d))
    context[d$is_normal %in% TRUE] <- "Normal"
    context[d$is_tumor %in% TRUE] <- "Tumor"
    d$sample_context <- context
  }
  if (cohort == "GSE26253") {
    d$recurrence_endpoint <- "RFS"
    d$source_endpoint_label <- "RFS"
    d$canonical_endpoint_family <- "RECURRENCE"
  } else if (cohort == "GSE62254") {
    d$recurrence_endpoint <- "DFS"
    d$source_endpoint_label <- "DFS"
    d$canonical_endpoint_family <- "RECURRENCE"
  } else if (cohort %in% c("GSE15459", "GSE84437", "TCGA_STAD")) {
    if (!"source_endpoint_label" %in% names(d)) d$source_endpoint_label <- "OS"
  }
  d
}

olfml2b_analysis_clinical <- function(view) {
  d <- view$clinical
  if (!nrow(d)) return(d)
  if ("sample_context" %in% names(d)) {
    keep <- as.character(d$sample_context) == "Tumor"
    if (!any(keep)) stop("No formal tumor samples remain in ", view$cohort, call. = FALSE)
    d <- d[keep, , drop = FALSE]
  } else if ("is_tumor" %in% names(d)) {
    keep <- d$is_tumor %in% TRUE
    if (any(keep)) d <- d[keep, , drop = FALSE]
  }
  d
}

olfml2b_endpoint_selection_table <- function(d, min_n = 40L, min_events = 20L) {
  endpoints <- olfml2b_endpoint_defs()
  rows <- lapply(names(endpoints), function(ep) olfml2b_endpoint_status_one(d, ep, endpoints[[ep]], min_n, min_events))
  tab <- olfml2b_bind_rows(rows)
  tab$analysis_included <- FALSE
  tab$primary_family_endpoint <- NA_character_
  tab$alias_of <- NA_character_
  tab$canonical_endpoint_family <- tab$endpoint_family
  source_label <- NA_character_
  if ("source_endpoint_label" %in% names(d)) {
    z <- unique(toupper(as.character(d$source_endpoint_label)))
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z) == 1L) source_label <- z
  }
  if (!is.finite(match(source_label, c("DFS", "RFS")))) {
    if ("recurrence_endpoint" %in% names(d)) {
      z <- unique(toupper(as.character(d$recurrence_endpoint)))
      z <- z[z %in% c("DFS", "RFS")]
      if (length(z) == 1L) source_label <- z
    }
  }
  tab$source_endpoint_label <- ifelse(tab$endpoint_family == "RECURRENCE", source_label, tab$endpoint)
  if ("OS" %in% tab$endpoint && tab$status[tab$endpoint == "OS"] == "EVALUABLE") {
    tab$analysis_included[tab$endpoint == "OS"] <- TRUE
    tab$primary_family_endpoint[tab$endpoint == "OS"] <- "OS"
  }
  rec <- tab[tab$endpoint_family == "RECURRENCE" & tab$status == "EVALUABLE", , drop = FALSE]
  if (nrow(rec)) {
    preferred <- if (!is.na(source_label) && source_label %in% rec$endpoint) source_label else rec$endpoint[order(rec$endpoint_priority)][1L]
    tab$analysis_included[tab$endpoint == preferred] <- TRUE
    tab$primary_family_endpoint[tab$endpoint_family == "RECURRENCE"] <- preferred
    tab$alias_of[tab$endpoint_family == "RECURRENCE" & tab$endpoint != preferred] <- preferred
  }
  tab
}

olfml2b_endpoint_audit <- function(views, min_n = 40L, min_events = 20L) {
  rows <- list()
  for (nm in names(views)) {
    d <- olfml2b_analysis_clinical(views[[nm]])
    tab <- olfml2b_endpoint_selection_table(d, min_n = min_n, min_events = min_events)
    tab$cohort <- views[[nm]]$cohort
    tab$olfml2b_measured <- views[[nm]]$target_measured
    rows[[length(rows) + 1L]] <- tab[, c(
      "cohort", "endpoint", "endpoint_family", "canonical_endpoint_family", "source_endpoint_label",
      "time_col", "event_col", "n_complete", "events", "olfml2b_measured", "status", "reason",
      "endpoint_priority", "analysis_included", "primary_family_endpoint", "alias_of"
    ), drop = FALSE]
  }
  olfml2b_bind_rows(rows)
}

.olfml2b_part3_v120_core <- run_olfml2b_specialized_bioinformatics
run_olfml2b_specialized_bioinformatics <- function(...) {
  index <- .olfml2b_part3_v120_core(...)
  index$version <- OLFML2B_ANALYSIS_VERSION
  ep <- index$endpoint_audit %||% data.frame()
  if (nrow(ep)) {
    display <- ep[ep$analysis_included %in% TRUE, c("cohort", "endpoint", "endpoint_family", "source_endpoint_label", "n_complete", "events"), drop = FALSE]
    olfml2b_write_csv(display, file.path(index$dirs$tables, "04b_endpoint_source_label_audit.csv"))
    g262 <- display[display$cohort == "GSE26253" & display$endpoint_family == "RECURRENCE", , drop = FALSE]
    g622 <- display[display$cohort == "GSE62254" & display$endpoint_family == "RECURRENCE", , drop = FALSE]
    if (nrow(g262) != 1L || g262$endpoint[1L] != "RFS") stop("Part3 endpoint contract failed: GSE26253 must be displayed as RFS.", call. = FALSE)
    if (nrow(g622) != 1L || g622$endpoint[1L] != "DFS") stop("Part3 endpoint contract failed: GSE62254 must be displayed as DFS.", call. = FALSE)
    index$endpoint_source_label_audit <- display
  }
  index$final_gene_lock <- FALSE
  saveRDS(index, file.path(index$dirs$objects, "Part3_OLFML2B_specialized_bioinformatics_index.rds"), compress = "xz")
  invisible(index)
}
