# =============================================================================
# OLFML2B-STAD Part9 publication outputs
# Main figures + supplementary figures + exact figure-source tables
# + curated supplementary tables
#
# v8.10.0 Figure 2 ecology-row gutter repair over the v8.9 layout repairs:
#   The complete Figure 2 top row is released from the unrelated programme-label
#   gutter required by Figure 2c, allowing Figure 2a to use the left-side space
#   Figure 2 panel membership, 55/45 top-row split and three-facet axes remain unchanged
#   The complete Figure 4 top row is released from the unrelated y-axis gutter
#   required by the association rows, allowing the atlas to use the left space
#   Figure 4 panel membership, 75/25 top-row split and spatial coordinates remain unchanged
#   Figure 5b and Figure 5d are released from the unrelated long-label gutter
#   inherited from Figure 5a, so their panels use the available left-side space
#   Figure 5d restores the prior unboxed two-line annotation grammar while
#   retaining the inset safe-zone position that avoids a bottom white band
#   The ggplot2 version preflight uses only exported utils APIs
#   (packageVersion and compareVersion)
#   Multi-section spatial coordinates are centred and isotropically normalized
#   for display, preserving within-section shape without claiming absolute
#   between-section size comparability
#   Main panel layout = 3 / 4 / 5 / 6 / 4 (22 panels total)
#   Figure 1 contains only raw expression evidence: unpaired, paired and cross-cohort stage context
#   The complete single-cell LODO matrix and the protein orientation-sentinel
#   audit are promoted to the main figures; duplicated supplementary displays are removed
#   Protein plex direction and leave-one-plex-out sensitivity use two aligned
#   effect-size columns, preventing unlike estimands from overplotting
#   Bootstrap attenuation uncertainty is a readable status matrix, not an interval forest
#   Attenuation, survival LOCO and technical diagnostics are supplementary
#   Figure 4 uses an atlas-plus-section-summary top row, a compact association row,
#   and an enlarged joint-spatial/source-boundary row
#   Cross-row patchwork alignment is released where unrelated y-axis label widths
#   created false gutters in Fig2a, Fig4a, Fig5b and Fig5d
#   Spatial point sizes and draw order are tuned from the actual 600-dpi v8.5 render
#   rather than from nominal ggplot dimensions
#   Fig4a section labels are anchored to each observed tissue cloud rather than strips
#   Long supplementary subtitles and legends are wrapped or stacked inside their panels
#   Inline statistics and legends are kept inside their owning panels
#   Figure 3b uses compact cohort labels; complete test details remain in source tables/captions
#   Figure 5b uses a compact 31% share of the protein row while Fig5c receives 69%
#   Figure 5d legends occupy a preflighted empty upper-left region instead of a bottom band
#   All figures are exported at final two-column width rather than oversized canvases
#   Only two main panels may use line/interval/forest grammar: Fig2c and Fig5a
#   Radar, polar, slope, dumbbell, lollipop, significance-only and main-panel boxplots are banned
#   Standard visual forms may repeat when the scientific estimand repeats; novelty of geometry
#   is never an acceptance criterion
#   Main panels contain no plot title or subtitle; panel letters, axes and captions carry meaning
#   Patient/sample is the display and inference unit wherever upstream data permit
#   ZIP packaging is performed inside Part9 and verified against the archive contents
#   Optional manuscript-media identity auditing detects stale embedded figure files
#   Wide-table reads inspect up to 1,000,000 rows before type assignment
#   Any parse problem is written to Audit and stops the run
# =============================================================================
# v8.10.1 namespace-isolation repair:
#   Part2 defines a project-local single-argument row-binding helper in the
#   global environment. Because .GlobalEnv precedes attached packages, that
#   helper can mask dplyr::bind_rows during a full Part0->Part9 run.
#   All Part9 row-binding calls are therefore explicitly namespaced.
#   No statistical, data-selection, panel-contract, or figure-design rule changed.
# =============================================================================
# v8.10.2 publication-table stage-test reconciliation:
#   The frozen Part3 clinical-context table stores omnibus one-way ANOVA F-test
#   global_p values for stage, whereas Main Figure 1c and the manuscript use the
#   prespecified omnibus Kruskal-Wallis stage test followed by BH-FDR.
#   Patient/stage counts are identical across the two publication paths.
#   Therefore, only Table S3 stage rows are canonicalised at Part9 export time:
#     - stage global_p <- the exact Figure 1c Kruskal-Wallis P value;
#     - all non-stage global_p values remain byte-for-value unchanged;
#     - an explicit method/source note is added to Table S3.
#   Part1-Part8 frozen outputs are not modified. No figure statistic is changed.
# =============================================================================
# v8.10.3 GSE84437 publication-facing provenance repair:
#   The frozen Part2 GEO cohort audit correctly contains 433 formal-source
#   samples and 431 complete OS observations, but older publication-facing
#   wording could collapse these two distinct quantities into "formal OS, 431".
#   In addition, the frozen cohort-audit row may leave the target-probe fields
#   blank even though Part2 exported an independent exact OLFML2B mapping audit.
#   Therefore Table S2 is reconciled at Part9 export time only:
#     - 433 = formal 2016 source-subseries patients/samples;
#     - 431 = OS-evaluable patients; 207 OS events;
#     - 50 = context-only superseries samples excluded from formal inference;
#     - exact OLFML2B probe IDs/mapping sources are joined from the frozen
#       Part2_GEO_OLFML2B_mapping_audit.csv;
#     - no Part2 object or analysis result is modified.
# =============================================================================
# v8.10.4 supplementary-table presentation and S12 reader-view repair:
#   * Adds a frozen presentation contract for Tables S1-S12 containing:
#       title, reader-facing note, missing-value rule, and abbreviations.
#   * Writes one sidecar NOTE file per supplementary table so the final ESM
#       workbook can be assembled directly from Part9 without editorial guessing.
#   * Rebuilds Table S12 from heterogeneous internal claim-contract schemas into
#       a non-sparse reader-facing six-column boundary table while preserving
#       every source statement and source identifier.
#   * Does not change any statistical result, figure, Part1-Part8 output, or
#       existing numerical table value.
# =============================================================================
# v8.10.5 Table S12 publication-language normalization:
#   QA of the v8.10.4 reader view showed five residual machine-facing values
#   (e.g. FALSE, NOT_AVAILABLE, NO_NOMINAL_ASSOCIATION) in the publication
#   statement column. These values are valid source-contract codes, but are not
#   suitable as reader-facing prose.
#   v8.10.5 therefore:
#     - preserves the original code in boundary_class;
#     - translates only the publication_statement field into literal,
#       non-causal prose that does not exceed the frozen source meaning;
#     - retains source_input for every row;
#     - adds a fail-closed QA guard against raw Boolean/status-code statements.
#   No statistical value, source contract, figure, or Part1-Part8 output changes.
# =============================================================================

PART9_VERSION <- paste0(
  "v8.10.5_20260827_TABLE_S12_",
  "PUBLICATION_LANGUAGE_NORMALIZATION"
)
PART9_SEED <- 20260825L

options(stringsAsFactors = FALSE, scipen = 6)
set.seed(PART9_SEED)

auto_install <- identical(toupper(Sys.getenv("OLFML2B_AUTO_INSTALL", "TRUE")), "TRUE")
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "forcats", "stringr",
  "tibble", "purrr", "scales", "patchwork", "ggrepel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L && auto_install) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

if (!"free" %in% getNamespaceExports("patchwork")) {
  stop(
    "Part9 v8.10.0 requires a patchwork release exporting free(). ",
    "Update patchwork and rerun the standalone launcher."
  )
}
ggplot2_version_string <- as.character(utils::packageVersion("ggplot2"))
if (utils::compareVersion(ggplot2_version_string, "3.5.0") < 0L) {
  stop(
    "Part9 v8.10.0 requires ggplot2 >= 3.5.0 for an explicit inside-legend position. ",
    "Update ggplot2 and rerun the standalone launcher."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(stringr)
  library(tibble)
  library(purrr)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

ROOT <- normalizePath(
  Sys.getenv("OLFML2B_ROOT", unset = "D:/OLFML2B_STAD"),
  winslash = "/", mustWork = FALSE
)
TABLE_ROOT <- normalizePath(
  Sys.getenv("OLFML2B_TABLES", unset = file.path(ROOT, "output", "tables")),
  winslash = "/", mustWork = FALSE
)

OUT_ID <- "Part9_Publication_v8_10_0"
FIG_ROOT <- normalizePath(
  Sys.getenv(
    "OLFML2B_PART9_FIG_OUT",
    unset = file.path(ROOT, "output", "figures", OUT_ID)
  ),
  winslash = "/", mustWork = FALSE
)
TAB_ROOT <- normalizePath(
  Sys.getenv(
    "OLFML2B_PART9_TABLE_OUT",
    unset = file.path(ROOT, "output", "tables", OUT_ID)
  ),
  winslash = "/", mustWork = FALSE
)

MAIN_FIG_DIR <- file.path(FIG_ROOT, "Main_Figures")
SUPP_FIG_DIR <- file.path(FIG_ROOT, "Supplementary_Figures")
MAIN_SOURCE_DIR <- file.path(TAB_ROOT, "Main_Figure_Source_Tables")
SUPP_SOURCE_DIR <- file.path(TAB_ROOT, "Supplementary_Figure_Source_Tables")
SUPP_TABLE_DIR <- file.path(TAB_ROOT, "Supplementary_Tables")
AUDIT_DIR <- file.path(TAB_ROOT, "Audit")
PACKAGE_PARENT <- normalizePath(
  Sys.getenv(
    "OLFML2B_PART9_PACKAGE_OUT",
    unset = file.path(ROOT, "output", "publication")
  ),
  winslash = "/", mustWork = FALSE
)
PACKAGE_DIR <- file.path(PACKAGE_PARENT, OUT_ID)
PACKAGE_ZIP <- file.path(PACKAGE_PARENT, paste0(OUT_ID, ".zip"))

clean_version_outputs <- identical(
  toupper(Sys.getenv("OLFML2B_PART9_CLEAN_OUTPUT", unset = "TRUE")), "TRUE"
)
version_output_dirs <- c(
  MAIN_FIG_DIR, SUPP_FIG_DIR, MAIN_SOURCE_DIR, SUPP_SOURCE_DIR,
  SUPP_TABLE_DIR, AUDIT_DIR
)
if (clean_version_outputs) {
  allowed_leaf_names <- c(
    "Main_Figures", "Supplementary_Figures", "Main_Figure_Source_Tables",
    "Supplementary_Figure_Source_Tables", "Supplementary_Tables", "Audit"
  )
  for (path in version_output_dirs) {
    if (!basename(path) %in% allowed_leaf_names) {
      stop("Refusing to clean an unexpected Part9 output directory: ", path)
    }
    if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  }
}

for (path in c(
  MAIN_FIG_DIR, SUPP_FIG_DIR, MAIN_SOURCE_DIR, SUPP_SOURCE_DIR,
  SUPP_TABLE_DIR, AUDIT_DIR
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
dir.create(PACKAGE_PARENT, recursive = TRUE, showWarnings = FALSE)

COL <- list(
  red = "#C94B40", orange = "#E69F00", blue = "#3B6FB6",
  purple = "#7057A3", teal = "#178F8B", gold = "#D8A62A",
  grey = "#6B6B6B", light = "#E5E7EB", black = "#222222",
  green = "#26976F"
)

# Figures are authored at the intended two-column print width.  Font sizes are
# therefore the final sizes and are not silently reduced by manuscript scaling.
FINAL_MAIN_WIDTH_IN <- 7.05
FINAL_SUPP_WIDTH_IN <- 7.05
COHORT_COLORS <- c(
  TCGA_STAD = "#C94B40", GSE62254 = "#E69F00", GSE15459 = "#CC79A7",
  GSE26253 = "#3B6FB6", GSE84437 = "#26976F",
  GSE134520 = "#178F8B", GSE150290 = "#C94B40",
  GSE167297 = "#3B6FB6", GSE183904 = "#7057A3"
)

feature_label <- function(x) {
  dplyr::recode(
    x,
    CAF_Core = "CAF core", CAF_ECM = "CAF-ECM",
    TGFb_Response = "TGF-β response",
    ECM_Remodeling = "ECM remodelling",
    Checkpoint_Exhaustion = "Checkpoint exhaustion",
    IFNg_Response = "IFN-γ response",
    CD8_Cytotoxic = "CD8 cytotoxicity",
    Proliferation_Control = "Proliferation control",
    Proliferation = "Proliferation",
    CAF_TGFb_axis = "CAF/TGF-β axis",
    CAF_TGFb_ECM_axis = "CAF/TGF-β/ECM axis",
    Immune_Exclusion_Index = "Immune exclusion index",
    Suppressive_TME_Index = "Suppressive TME index",
    target_expr = "OLFML2B",
    Inflamed_5gene_exploratory = "Inflamed 5-gene",
    Stromal_Remodeling_axis = "Stromal remodelling axis",
    M2_Macrophage = "M2 macrophage",
    Myeloid_Macrophage = "Myeloid/macrophage",
    Myofibroblast = "Myofibroblast", Fibroblast = "Fibroblast",
    Smooth_Muscle = "Smooth muscle", Pericyte = "Pericyte",
    Epithelial = "Epithelial", Endothelial = "Endothelial",
    Epithelial_Differentiation = "Epithelial differentiation",
    Endothelial_Angiogenic = "Endothelial/angiogenic",
    Inflammatory_Fibroblast = "Inflammatory fibroblast",
    T_NK = "T/NK", B_Plasma = "B/plasma",
    .default = stringr::str_replace_all(x, "_", " ")
  )
}

cohort_label <- function(x) dplyr::recode(x, TCGA_STAD = "TCGA-STAD", .default = x)

fmt_p <- function(p) {
  superscript <- function(x) {
    chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(x))
  }
  vapply(p, function(value) {
    if (is.na(value)) return("NA")
    if (value <= 0) return("< 1 × 10⁻³⁰⁰")
    if (value < 1e-4) {
      exponent <- floor(log10(value))
      mantissa <- value / (10^exponent)
      return(paste0(formatC(mantissa, digits = 1, format = "f"),
                    " × 10", superscript(exponent)))
    }
    formatC(value, digits = 3, format = "f")
  }, character(1))
}

pairwise_auc <- function(positive_values, negative_values) {
  positive_values <- positive_values[is.finite(positive_values)]
  negative_values <- negative_values[is.finite(negative_values)]
  if (length(positive_values) < 1L || length(negative_values) < 1L) return(NA_real_)
  differences <- outer(positive_values, negative_values, FUN = "-")
  mean(differences > 0) + 0.5 * mean(differences == 0)
}

bootstrap_two_group_effect <- function(positive_values, negative_values,
                                       B = 2000L, seed = PART9_SEED) {
  positive_values <- positive_values[is.finite(positive_values)]
  negative_values <- negative_values[is.finite(negative_values)]
  observed_auc <- pairwise_auc(positive_values, negative_values)
  if (length(positive_values) < 2L || length(negative_values) < 2L || !is.finite(observed_auc)) {
    return(tibble(auc = observed_auc, auc_low = NA_real_, auc_high = NA_real_,
                  rank_biserial = 2 * observed_auc - 1,
                  rb_low = NA_real_, rb_high = NA_real_))
  }
  set.seed(seed)
  boot_auc <- replicate(B, {
    pairwise_auc(
      sample(positive_values, length(positive_values), replace = TRUE),
      sample(negative_values, length(negative_values), replace = TRUE)
    )
  })
  auc_ci <- unname(quantile(boot_auc, c(0.025, 0.975), na.rm = TRUE, type = 6))
  tibble(
    auc = observed_auc, auc_low = auc_ci[1], auc_high = auc_ci[2],
    rank_biserial = 2 * observed_auc - 1,
    rb_low = 2 * auc_ci[1] - 1, rb_high = 2 * auc_ci[2] - 1
  )
}

bootstrap_median_ci <- function(values, B = 5000L, seed = PART9_SEED) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(c(NA_real_, NA_real_))
  set.seed(seed)
  draws <- replicate(B, median(sample(values, length(values), replace = TRUE)))
  unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE, type = 6))
}

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.16), hjust = 0),
      plot.subtitle = element_text(colour = COL$grey, size = rel(0.86), hjust = 0),
      plot.tag = element_text(face = "bold", size = rel(1.10)),
      axis.title = element_text(colour = COL$black),
      axis.text = element_text(colour = COL$black),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = COL$black),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.82)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

# Applied to complete main-figure compositions.  Panel-specific statistical prose
# belongs in the external legend/caption, not inside the rendered panel.
main_panel_sanitizer <- theme(
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank(),
  plot.tag = element_text(face = "bold", size = 9.0, colour = COL$black)
)

READ_PARSE_AUDIT <- tibble()
PARSE_AUDIT_PATH <- file.path(AUDIT_DIR, "08_READR_PARSE_PROBLEMS_v8_10_0.csv")

read_contract <- function(
    part, file, required_columns,
    selected_columns = NULL, col_types = NULL,
    parse_policy = c("stop", "warn")
) {
  parse_policy <- match.arg(parse_policy)
  path <- file.path(TABLE_ROOT, part, file)
  if (!file.exists(path)) stop("Required frozen table is absent: ", path)
  dat <- suppressWarnings(suppressMessages({
    if (is.null(selected_columns)) {
      readr::read_csv(
        path, show_col_types = FALSE, progress = FALSE,
        guess_max = 1000000L, col_types = col_types
      )
    } else {
      readr::read_csv(
        path,
        col_select = tidyselect::all_of(selected_columns),
        show_col_types = FALSE, progress = FALSE,
        guess_max = 1000000L, col_types = col_types
      )
    }
  }))
  parse_issues <- readr::problems(dat)
  if (nrow(parse_issues) > 0L) {
    READ_PARSE_AUDIT <<- dplyr::bind_rows(
      READ_PARSE_AUDIT,
      parse_issues %>% mutate(source_file = path, .before = 1)
    )
    readr::write_csv(
      READ_PARSE_AUDIT %>% mutate(status = "REVIEW_SOURCE_CELL"),
      PARSE_AUDIT_PATH,
      na = ""
    )
    parse_message <- paste0(
      "Parsing audit: ", nrow(parse_issues), " issue(s) in ", path,
      ". See ", PARSE_AUDIT_PATH, "."
    )
    if (identical(parse_policy, "stop")) stop(parse_message)
    warning(parse_message)
  }
  absent <- setdiff(required_columns, names(dat))
  if (length(absent) > 0L) {
    stop("Schema mismatch in ", path, "; absent columns: ", paste(absent, collapse = ", "))
  }
  dat
}

assert_nonempty <- function(dat, label) {
  if (!is.data.frame(dat) || nrow(dat) < 1L) {
    stop("Panel preflight failed: ", label, " has zero eligible rows after filtering.")
  }
  invisible(dat)
}

source_registry <- tibble()

write_panel_source <- function(
    dat, branch = c("MAIN", "SUPPLEMENTARY"), file_name, panel_id,
    source_inputs, transformation, inference_unit, claim_ceiling
) {
  branch <- match.arg(branch)
  assert_nonempty(dat, panel_id)
  out_dir <- if (branch == "MAIN") MAIN_SOURCE_DIR else SUPP_SOURCE_DIR
  path <- file.path(out_dir, file_name)
  readr::write_csv(dat, path, na = "")
  source_registry <<- dplyr::bind_rows(
    source_registry,
    tibble(
      branch = branch, panel_id = panel_id, source_table_file = basename(path),
      source_inputs = source_inputs, transformation = transformation,
      inference_unit = inference_unit, claim_ceiling = claim_ceiling,
      n_rows = nrow(dat), n_columns = ncol(dat)
    )
  )
  invisible(dat)
}

save_figure_bundle <- function(plot_object, branch = c("MAIN", "SUPPLEMENTARY"),
                               stem, width, height) {
  branch <- match.arg(branch)
  out_dir <- if (branch == "MAIN") MAIN_FIG_DIR else SUPP_FIG_DIR
  target <- file.path(out_dir, stem)
  ggsave(
    paste0(target, ".png"), plot_object, width = width, height = height,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  ggsave(
    paste0(target, ".tiff"), plot_object, width = width, height = height,
    units = "in", dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(
    paste0(target, ".pdf"), plot_object, width = width, height = height,
    units = "in", device = pdf_device, bg = "white", limitsize = FALSE
  )
  invisible(target)
}

main_panel_ids <- c(
  paste0("Fig1", letters[1:3]), paste0("Fig2", letters[1:4]),
  paste0("Fig3", letters[1:5]), paste0("Fig4", letters[1:6]),
  paste0("Fig5", letters[1:4])
)
main_figure_ids <- rep(paste("Figure", 1:5), c(3, 4, 5, 6, 4))
figure_questions <- c(
  "Figure 1" = "How does OLFML2B differ between tumour and normal tissue, and how does it vary across pathological stage?",
  "Figure 2" = "Which bulk-tissue ecological state travels with OLFML2B across cohorts?",
  "Figure 3" = "Which cells carry OLFML2B and which programmes accompany that localization?",
  "Figure 4" = "Is the spatial pattern reproduced across sections and distinguishable from competing sources?",
  "Figure 5" = "How strong and transportable are the clinical and orthogonal protein associations?"
)

panel_contract_main <- tibble(
  branch = "MAIN",
  figure_id = main_figure_ids,
  panel_id = main_panel_ids,
  scientific_question = unname(figure_questions[main_figure_ids]),
  claim = c(
    "Tumour-normal expression difference", "Within-patient tumour-normal increase",
    "Cross-cohort pathological-stage context",
    "Continuous sample-level stromal association", "Cohort-by-program replication",
    "Pooled ecological effect with heterogeneity", "Composite stromal-state shift",
    "Dataset-specific cellular localization", "Paired patient pseudobulk increase",
    "Paired target-programme coupling", "OLFML2B-positive-cell programme context",
    "Cellular localization survives dataset omission", "Target detection across all spatial sections",
    "Section-level target-detection range", "Same-spot and neighbour concordance", "Target-context threshold overlap",
    "Representative target-high and CAF-ECM-high spatial overlap", "Competing-source disambiguation",
    "Adjusted prognosis and pooled uncertainty", "Case-paired protein direction",
    "Plex direction and omission stability", "Protein-assay orientation and biological plausibility"
  ),
  inference_unit = c(
    "Sample for unpaired comparison", "Patient", "Patient within cohort",
    "Sample within cohort", "Patient/sample within cohort", "Cohort", "Patient/sample within cohort",
    "Official sample/patient within dataset", "Patient", "Patient", "Official sample/patient",
    "Dataset omission", "Section; patient for inference elsewhere", "Section", "Patient after section collapse",
    "Patient after section collapse", "Section for display", "Patient after section collapse",
    "Patient within cohort; cohort for meta-analysis", "Case", "Analytical plex or case after plex omission", "Case-paired sentinel"
  ),
  contrast = c(
    "Tumour versus normal", "Tumour minus normal", "Stage I-IV within cohort",
    "OLFML2B decile versus programme score", "Across cohorts and programmes", "Positive versus null association",
    "OLFML2B-high minus low display split", "Across annotated compartments", "Tumour minus adjacent/normal",
    "Within-patient target change versus programme change", "OLFML2B-positive minus negative cells",
    "All datasets versus one omitted", "Across observed sections", "Across evaluable sections", "Same spot versus neighbour",
    "Target-high with context-high", "Upper-quartile overlap", "Standardized competing-source coefficients",
    "HR per 1-SD OLFML2B", "Tumour minus normal protein", "Across plexes and plex omissions", "Prespecified sentinel orientation"
  ),
  estimand = c(
    "Rank-biserial expression effect", "Median paired delta", "Within-cohort stage distribution and median",
    "Decile median", "Spearman rho", "Random-effects pooled rho", "Median programme-score difference",
    "Detection fraction", "Median paired fraction delta", "Spearman rho of paired changes",
    "Median programme-score difference", "Top-three support after omission", "Within-section percentile",
    "Positive spot fraction and analysed spots",
    "Patient-level median rho", "Patient-level median overlap fraction", "High-high overlap class",
    "Median standardized ridge coefficient", "Hazard ratio", "Paired protein delta", "Plex-specific and post-omission paired delta",
    "Sentinel median paired delta"
  ),
  uncertainty = c(
    "Raw samples and frozen Wilcoxon test", "Bootstrap CI in caption/source", "Raw patients and global Kruskal-Wallis FDR",
    "Formal continuous tests upstream", "Cohort replication", "95% CI and prediction interval",
    "Cross-cohort spread", "Dataset replication", "Raw patients and cohort median",
    "Two cohorts shown separately", "Across-dataset consistency", "Exhaustive eligible omissions",
    "All observed sections", "All evaluable sections; descriptive", "Bootstrap CI retained in source", "Bootstrap CI retained in source",
    "Display threshold only", "Bootstrap CI retained in source", "95% CI; PI only when identified",
    "Raw cases; HL and bootstrap CI in caption", "All plexes and leave-one-plex-out estimates", "All evaluable prespecified sentinels"
  ),
  panel_role = c(
    "Primary phenotype", "Patient-level confirmation", "Cross-cohort clinical context",
    "Raw-shape evidence", "Replication", "Synthesis", "State summary",
    "Localization", "Patient-level confirmation", "Ecological bridge", "Cell-state extension",
    "Robustness", "Spatial atlas", "Cross-section summary", "Spatial association", "Spatial overlap",
    "Spatial anchor", "Source boundary", "Clinical synthesis", "Orthogonal assay",
    "Batch robustness and influence", "Assay-orientation audit"
  ),
  analysis_branch = c(
    "Part1/Part3", "Part3", "Part1/Part2", "Part4", "Part4", "Part4", "Part4",
    "Part6", "Part6", "Part6", "Part6", "Part6", "Part7", "Part7", "Part7", "Part7",
    "Part7", "Part7", "Part3", "Part5", "Part5", "Part5"
  ),
  source_table = c(
    "MF1a_TCGA_tumour_normal_distribution.csv", "MF1b_TCGA_patient_paired_values.csv",
    "MF1c_stage_patient_values.csv",
    "MF2a_sample_level_decile_display.csv", "MF2b_cross_cohort_ecology.csv",
    "MF2c_pooled_ecology_nested_intervals.csv", "MF2d_composite_state_shift.csv",
    "MF3a_dataset_localization.csv", "MF3b_patient_paired_pseudobulk.csv",
    "MF3c_paired_change_program_coupling.csv", "MF3d_positive_cell_program_context.csv",
    "MF3e_localization_LODO_matrix.csv", "MF4a_all_section_spatial_atlas.csv",
    "MF4b_section_detection.csv", "MF4c_same_spot_neighbor_concordance.csv",
    "MF4d_threshold_overlap.csv", "MF4e_representative_joint_spatial.csv",
    "MF4f_competing_source_ridge.csv",
    "MF5a_adjusted_survival_and_meta_forest.csv", "MF5b_protein_paired_distribution.csv",
    "MF5c_two_column_plex_robustness.csv", "MF5d_orientation_sentinels.csv"
  ),
  visual_form = c(
    "sample violin and jitter", "signed histogram and patient rug", "faceted patient scatter and median",
    "faceted decile scatter", "bubble matrix", "nested interval blocks", "cohort swarm",
    "faceted detection bubbles", "signed patient swarm", "two-cohort dot rows", "faceted mirrored bars",
    "bubble omission matrix", "section-faceted spatial atlas", "ranked section detection points",
    "direct-labelled scatter", "context overlap points", "four-class bivariate spot map",
    "diverging coefficient bars",
    "grouped meta forest", "case violin and jitter", "faceted plex robustness scatter", "sentinel direction audit"
  ),
  claim_ceiling = c(
    "Unpaired expression phenotype", "Paired expression phenotype", "Cross-sectional stage context; not progression",
    "Descriptive shape; continuous tests upstream", "Association, not causality",
    "Transportability limited by heterogeneity", "Display split; not a threshold", "Dominant localization, not exclusivity",
    "Patient-level abundance, not cell-level inference", "Ecological coupling, not mechanism",
    "State association, not regulation", "Dataset robustness, not universal exclusivity",
    "Observed spatial distribution", "Cross-section detection; not patient-level independence",
    "Spatial association, not interaction", "Descriptive overlap, not cell identity",
    "Representative overlap only", "Exploratory source disambiguation", "Observational prognosis; no clinical utility",
    "Single-cohort orthogonal direction", "Plex-sensitive robustness boundary", "Orientation audit, not target validation"
  ),
  line_interval_forest = main_panel_ids %in% c("Fig2c", "Fig5a"),
  matrix_like = main_panel_ids %in% c("Fig2b", "Fig3e"),
  significance_only = FALSE,
  exotic_geometry = FALSE
)

supplement_contract <- tibble(
  branch = "SUPPLEMENTARY",
  figure_id = paste("Figure S", 1:6, sep = ""),
  panel_id = paste0("FigS", 1:6),
  scientific_question = c(
    "Clinical and survival diagnostics", "Bulk TME sensitivity and collinearity",
    "Single-cell orthogonal state", "Spatial sensitivity",
    "Protein pairing eligibility", "Exploratory ICI boundary"
  ),
  claim = scientific_question,
  inference_unit = "As specified per panel source table",
  contrast = "Diagnostic or exploratory",
  estimand = "Panel-specific",
  uncertainty = "Panel-specific",
  panel_role = "Qualifier",
  analysis_branch = paste0("Supplementary S", 1:6),
  source_table = NA_character_, visual_form = "multipanel",
  claim_ceiling = "Does not independently raise the main-text claim ceiling",
  line_interval_forest = FALSE, matrix_like = FALSE,
  significance_only = FALSE, exotic_geometry = FALSE
)
panel_contract <- dplyr::bind_rows(panel_contract_main, supplement_contract)
main_budget <- panel_contract %>% filter(branch == "MAIN")

expected_layout <- tibble(
  figure_id = paste("Figure", 1:5), expected_panels = c(3L, 4L, 5L, 6L, 4L)
)
observed_layout <- main_budget %>% count(figure_id, name = "observed_panels") %>%
  right_join(expected_layout, by = "figure_id")
if (nrow(main_budget) != 22L || any(observed_layout$observed_panels != observed_layout$expected_panels)) {
  stop("Main-figure panel contract must be exactly 3/4/5/6/4 (22 panels).")
}
if (!setequal(main_budget$panel_id[main_budget$line_interval_forest], c("Fig2c", "Fig5a"))) {
  stop("Only Fig2c and Fig5a may use line/interval/forest grammar.")
}
if (sum(main_budget$matrix_like) > 2L) stop("Main-figure matrix budget exceeded.")
if (any(main_budget$significance_only)) stop("Significance-only main panels are prohibited.")
if (any(main_budget$exotic_geometry)) stop("Exotic main-panel geometry is prohibited.")
readr::write_csv(panel_contract, file.path(AUDIT_DIR, "01_PANEL_CONTRACT_v8_10_0.csv"), na = "")

visual_budget_audit <- tibble(
  metric = c(
    "TOTAL_MAIN_PANELS", "LINE_INTERVAL_FOREST_TOTAL", "MATRIX_LIKE_MAIN_PANELS",
    "SIGNIFICANCE_ONLY_MAIN_PANELS", "EXOTIC_GEOMETRY_MAIN_PANELS"
  ),
  observed = c(
    nrow(main_budget), sum(main_budget$line_interval_forest), sum(main_budget$matrix_like),
    sum(main_budget$significance_only), sum(main_budget$exotic_geometry)
  ),
  required_or_maximum = c(22L, 2L, 2L, 0L, 0L),
  rule = c("EQUAL", "EQUAL", "MAXIMUM", "EQUAL", "EQUAL")
) %>% mutate(
  status = case_when(
    rule == "EQUAL" & observed == required_or_maximum ~ "PASS",
    rule == "MAXIMUM" & observed <= required_or_maximum ~ "PASS",
    TRUE ~ "FAIL"
  )
)
if (any(visual_budget_audit$status != "PASS")) stop("Main visual-grammar contract failed.")
readr::write_csv(
  visual_budget_audit,
  file.path(AUDIT_DIR, "01C_VISUAL_GRAMMAR_AUDIT_v8_10_0.csv")
)

geometry_signature_audit <- main_budget %>%
  transmute(
    panel_id, visual_form,
    perceptual_family = case_when(
      line_interval_forest ~ "LINE_INTERVAL_FOREST",
      matrix_like ~ "MATRIX",
      str_detect(visual_form, "spatial") ~ "SPATIAL",
      str_detect(visual_form, "scatter|swarm|raincloud|density|rug") ~ "POINT_DISTRIBUTION",
      str_detect(visual_form, "bar|mosaic|tile|card") ~ "BLOCK",
      TRUE ~ "OTHER_STANDARD"
    ),
    acceptance_basis = "Scientific role and estimand, not uniqueness of geometry",
    status = "PASS"
  )
readr::write_csv(
  geometry_signature_audit,
  file.path(AUDIT_DIR, "01D_GEOMETRY_SIGNATURE_AUDIT_v8_10_0.csv")
)

palette_contract <- tribble(
  ~semantic_role, ~hex, ~allowed_use,
  "Target / tumour / positive direction", COL$red,
  "OLFML2B, tumour increase, positive-direction effect or pooled summary",
  "Reference / inverse / contrast", COL$blue,
  "Normal tissue, negative change or inverse-direction effect; not model class",
  "Stromal localization", COL$purple,
  "Fibroblast/CAF source support and stromal emphasis",
  "Orthogonal context", COL$teal,
  "Technical or orthogonal contextual evidence",
  "Neutral evidence", COL$grey,
  "Cohort observations, uncertainty and non-highlighted elements",
  "Summary emphasis", COL$gold,
  "Median or prespecified summary marker only"
)
readr::write_csv(
  palette_contract,
  file.path(AUDIT_DIR, "00B_SEMANTIC_PALETTE_CONTRACT_v8_10_0.csv")
)

final_size_contract <- tibble(
  branch = c(rep("MAIN", 5), rep("SUPPLEMENTARY", 6)),
  figure = c(paste("Figure", 1:5), paste0("Figure S", 1:6)),
  width_in = c(rep(FINAL_MAIN_WIDTH_IN, 5), rep(FINAL_SUPP_WIDTH_IN, 6)),
  authored_at_final_width = TRUE,
  minimum_body_text_pt = 6.4,
  panel_tag_pt = 9.0,
  rule = "No manuscript-side downscaling is required for two-column placement"
)
readr::write_csv(
  final_size_contract,
  file.path(AUDIT_DIR, "00C_FINAL_SIZE_CONTRACT_v8_10_0.csv")
)

figure_narrative_contract <- tribble(
  ~figure, ~scientific_question, ~panel_sequence, ~claim_ceiling,
  "Figure 1", "How does OLFML2B differ between tumour and normal tissue, and how does it vary across pathological stage?",
  "Unpaired tumour-normal distribution -> paired patient delta -> cross-cohort stage context",
  "Expression phenotype and cross-sectional clinical context; stage is not interpreted as progression",
  "Figure 2", "Which bulk-tissue ecological state travels with OLFML2B across cohorts?",
  "Sample-level shape -> cohort map -> pooled uncertainty -> composite states",
  "Ecological association and shared information; not causality or mediation",
  "Figure 3", "Which cellular compartments carry OLFML2B and does detection increase within patients?",
  "Dataset localization -> paired pseudobulk -> programme coupling -> state context -> LODO",
  "Dominant localization and associated state; not cell exclusivity or regulation",
  "Figure 4", "Is the stromal localization locally coherent and distinguishable from competing sources?",
  "All-section atlas -> section detection -> local/neighbour concordance -> overlap -> representative joint map -> source disambiguation",
  "Spatial association; not cell-cell interaction or mechanism",
  "Figure 5", "How strong, transportable and orthogonally supported is the prognostic association?",
  "Adjusted cohort/meta forest -> case protein direction -> plex stability -> orientation sentinels",
  "Observational prognosis with small-k uncertainty; no clinical utility claim"
)
readr::write_csv(
  figure_narrative_contract,
  file.path(AUDIT_DIR, "00_FIGURE_NARRATIVE_CONTRACT_v8_10_0.csv")
)

panel_rebalance_audit <- tribble(
  ~action, ~panel_id, ~from_branch, ~to_branch, ~decision_basis,
  "DELETE", "Old Fig1a", "MAIN", "REMOVED", "The evidence-path schematic occupied a full row without adding measured data",
  "PROMOTE", "Old FigS1a", "SUPPLEMENTARY", "MAIN Fig1c", "Cross-cohort patient-level stage expression completes the Figure 1 clinical context",
  "MOVE", "Old Fig1d", "MAIN Figure 1", "MAIN Fig4e", "All spatial evidence is consolidated in Figure 4",
  "PROMOTE", "Old FigS2a", "SUPPLEMENTARY", "MAIN Fig2a", "Continuous sample-level shape precedes cohort/meta summaries",
  "DEMOTE", "Old Fig2d", "MAIN", "SUPPLEMENTARY FigS2", "Seventeen of eighteen attenuation intervals cross zero",
  "PROMOTE", "Old FigS3a", "SUPPLEMENTARY", "MAIN Fig3d", "Across-dataset positive-cell programme context extends localization",
  "PROMOTE", "Complete old FigS3a", "SUPPLEMENTARY", "MAIN Fig3e", "Full-compartment LODO matrix is more informative than a single-row solid status strip",
  "PROMOTE", "Old FigS4c", "SUPPLEMENTARY", "MAIN Fig4d", "Threshold overlap directly addresses target-context spatial placement",
  "PROMOTE", "Old FigS4d", "SUPPLEMENTARY", "MAIN Fig4b", "The ranked section-detection summary fills the atlas aspect mismatch and quantifies cross-section breadth",
  "PROMOTE", "Old FigS5b", "SUPPLEMENTARY", "MAIN Fig5d", "Orientation sentinels qualify the orthogonal protein assay",
  "MERGE", "Old Fig5c/Fig5d", "MAIN", "MAIN Fig5c", "Plex-specific and post-omission effects are separated into aligned columns rather than overlapping or using solid colour strips",
  "REBALANCE", "Fig5b/c/d", "MAIN Figure 5", "MAIN Figure 5", "The single-group distribution is narrowed, the plex panel is widened, and sentinel legends move into a validated empty region",
  "DEMOTE", "Old Fig5b", "MAIN", "SUPPLEMENTARY FigS1", "Small-k survival LOCO is a fragility diagnostic",
  "STANDARDIZE", "All main figures", "MAIN", "MAIN", "Plot titles/subtitles removed; only two interval/forest panels retained"
)
readr::write_csv(
  panel_rebalance_audit,
  file.path(AUDIT_DIR, "01B_PANEL_REBALANCE_DECISIONS_v8_10_0.csv")
)

render_qa_repair <- tribble(
  ~figure_panel, ~observed_problem, ~v8_7_repair, ~information_loss,
  "All main panels", "Poster-like titles and subtitles dominated the marks", "Titles/subtitles stripped at composition level", "Captions must be regenerated externally",
  "Figure 1", "The hollow-node evidence path remained a low-information full-width row", "Delete the schematic; pair unpaired and paired expression above a full-width cross-cohort stage panel", "Stage remains descriptive and is not interpreted as progression",
  "Figure 2", "Repeated interval grammar and weak attenuation panel", "Only pooled ecology retains interval grammar", "Attenuation remains in FigS2",
  "Figure 3a/e", "Fixed empty cell-type rows and nearly invisible non-evaluable marks created low-density whitespace", "Free-y localization facets plus explicit non-evaluable crosses", "Absence remains represented in the LODO panel",
  "Figure 4", "The full-width 5-by-2 atlas was aspect-mismatched while the representative joint map was undersized", "Atlas paired with ranked section detection; joint map enlarged in a taller final row", "All ten sections remain in the main figure",
  "Figure 4a", "The atlas top row inherited the wide y-axis gutter required by the association row below", "Release the complete 75/25 top-row patch on the left so the atlas uses the otherwise empty strip", "All ten sections, section labels and the detection summary remain unchanged",
  "Figure 4d/f", "Two adjacent solid bar summaries produced repetitive visual mass", "Overlap is point-only; only source disambiguation retains diverging bars", "Uncertainty values remain in source tables and captions",
  "Figure 5d", "Repeated per-gene sample-size labels and a tall stacked legend dominated the sentinel panel", "One panel-level n statement, compact two-line guides and no per-gene n labels", "No sentinel rows removed",
  "Figure 2a", "The complete 2a/2b top row inherited the long programme-label gutter required by Figure 2c despite freeing only the inner 2a plot", "Release the complete 55/45 top-row patch on the left at panel level", "No observations, cohorts, programme axes, facets or legends removed",
  "Figure 2b/c", "The Fig2b colourbar labels were compressed and Fig2c text overlapped the rightmost intervals", "Stack Fig2b guides and add a true blank annotation column to Fig2c", "Full CI, PI, I² and FDR values remain in panel sources and captions",
  "Figure 3a", "The five-break cell-count legend forced the Cells title beyond the left panel boundary", "Use three labelled cell-count anchors and vertically stack the two compact guides", "Bubble area continues to encode analysed cell count",
  "Figure 3b", "Comparison prose and FDR on the x axis collided with an excessively long y title", "Use cohort/n/q x labels and the compact y title Paired Δ positive fraction", "Contrast type and complete statistics remain in the exported source table/caption",
  "Figure 4a/e", "Spatial points were blotchy at 600 dpi and facet strips remained detached from tissue clouds", "Use 0.28-mm atlas points, data-anchored section labels and four ordered point-size layers in the joint map", "Every section, spot and joint class remains displayed",
  "Figure 4c", "The CAF-ECM label overlapped its point and neighbouring labels", "Use the frozen feature-specific label coordinates already present in the panel source", "No programmes or effect estimates removed",
  "Figure 5a", "Pooled P/I²/k strings extended beyond the right page edge", "Use compact pooled labels, log-scale right expansion and in-panel clipping", "All survival estimates and uncertainty intervals remain visible",
  "Figure 5b/d", "Both panels inherited the longer Figure 5a cohort-label gutter, leaving a large unused strip on the left", "Free the left panel boundary of Fig5b and Fig5d independently while retaining their own adjacent axis text", "Every paired case, gene label and sentinel remains represented",
  "Figure 5c", "Plex-specific medians and post-omission pooled medians overlapped in one 12-category x-axis display", "Separate the estimands into two aligned effect-size columns with plexes on the y axis and state the boundary symbol in the right facet strip", "All 12 plex-specific and 12 post-omission estimates remain displayed",
  "Figure 5b/c/d", "The v8.6 single-group Fig5b was too wide and the two-line Fig5d bottom legend created a large white band", "Use 31/69 middle-row widths and restore the prior unboxed two-line guide grammar inside a validated empty upper-left region", "No cases, plexes or sentinels removed; the run stops if the inset region contains data",
  "Supplementary text", "FigS1c labels covered estimates and long FigS2b/S4b/S6b text crossed panel boundaries", "Reserve an S1c detail column and wrap/stack long supplementary titles, subtitles and guides", "No diagnostic rows removed",
  "All main figures", "Oversized 10.6-11.2 inch canvases would force unreadable manuscript downscaling", "Author directly at 7.05-inch two-column width", "No statistical or graphical rows removed"
)
readr::write_csv(
  render_qa_repair,
  file.path(AUDIT_DIR, "01E_RENDER_QA_REPAIR_v8_10_0.csv")
)

method_semantic_audit <- tribble(
  ~audit_item, ~preflight_rule, ~status_if_run_completes,
  "SURVIVAL_PI_IDENTITY", "PI is plotted only when finite and distinct from the pooled CI", "PASS",
  "ECOLOGY_HETEROGENEITY", "I² and PI crossing are visible in Main Fig. 2c", "PASS",
  "PROTEIN_ESTIMAND_ALIGNMENT", "Wilcoxon is paired with Hodges-Lehmann pseudomedian/CI", "PASS",
  "PROTEIN_PLEX_SENSITIVITY", "Maximum leave-one-plex-out P is exported for the caption and all omissions are displayed", "PASS",
  "PLEX_ESTIMAND_SEPARATION", "Plex-specific and post-omission pooled medians occupy separate aligned facets and cannot overplot", "PASS",
  "ATTRITION_BRANCHING", "Tumour and normal counts cannot be interpreted as a single funnel", "PASS",
  "SINGLE_CELL_COLOUR_RANGE", "Symmetric colour limit covers every plotted effect", "PASS",
  "ICI_COMMON_EFFECT_SCALE", "Every response-screen x value is a rank-biserial effect", "PASS",
  "EXACT_ZERO_P_DISPLAY", "P = 0 is shown as a bound and not an infinite numeric estimate", "PASS",
  "MAIN_MATRIX_BUDGET", "Only Main Fig. 2b and Fig. 3e use matrix-like circle layouts", "PASS",
  "MAIN_INTERVAL_BUDGET", "Only Main Fig. 2c and Fig. 5a use line/interval/forest grammar", "PASS",
  "FINAL_SIZE_AUTHORING", "All main figures are authored at 7.05-inch two-column width", "PASS",
  "FIGURE4_DENSITY_BALANCE", "Figure 4 contains six panels in atlas-summary, association and spatial-source rows", "PASS",
  "INLINE_LABEL_CONTAINMENT", "Fig2c, Fig2d and Fig5a place inline statistical labels inside their own panel viewports", "PASS",
  "LEGEND_CONTAINMENT", "Fig3a uses three compact cell-count anchors and a two-row guide box", "PASS",
  "AXIS_LABEL_HIERARCHY", "Fig3b and Fig5b/c use short axis titles while detailed tests remain in captions/source tables", "PASS",
  "NESTED_PANEL_BALANCE", "Fig5b/Fig5c receive 31%/69% of the protein row and the row receives 0.80 relative height", "PASS",
  "SENTINEL_LEGEND_SAFE_ZONE", "Fig5d inset guides are allowed only when the upper-left safety region contains no sentinel estimate", "PASS",
  "EXPORTED_VERSION_CHECK_API", "ggplot2 compatibility preflight uses exported utils::packageVersion and utils::compareVersion only", "PASS",
  "CROSS_ROW_AXIS_RELEASE", "Fig2a, the Fig4 atlas row, Fig5b and Fig5d do not inherit unrelated long y-label gutters from other rows", "PASS",
  "FIGURE2_TOP_ROW_GUTTER_RELEASE", "The complete Fig2 55/45 ecology row is freed on the left while the three Fig2a facets remain internally aligned", "PASS",
  "FIGURE4_TOP_ROW_GUTTER_RELEASE", "The complete Fig4 75/25 atlas-detection row is freed on the left while all lower-row axes remain unchanged", "PASS",
  "FIGURE5_PRIOR_ANNOTATION_STYLE", "Fig5d uses the prior unboxed one-row colour plus one-row direction guide grammar without restoring a bottom legend band", "PASS",
  "SPATIAL_MARK_HIERARCHY", "Fig4a uses smaller points and Fig4e draws Neither, single-high and Joint-high classes in ordered layers", "PASS",
  "SUPPLEMENTARY_TEXT_CONTAINMENT", "FigS1c/S2b/S4b/S6b labels, titles, subtitles and guides remain inside their panels", "PASS",
  "PROMOTION_NON_DUPLICATION", "Promoted panels are removed or restricted to disjoint rows in supplements", "PASS",
  "STAGE_CLAIM_BOUNDARY", "Main Fig. 1c is labelled as cross-sectional context and never as progression", "PASS",
  "SENTINEL_POINT_ONLY", "Main Fig. 5d encodes direction with points and contains no area-proportional colour bars", "PASS",
  "TECHNICAL_RESIDUAL_LABEL", "Residual collapse is called attenuation and not robustness", "PASS"
)
readr::write_csv(
  method_semantic_audit,
  file.path(AUDIT_DIR, "01F_METHOD_SEMANTIC_AUDIT_v8_10_0.csv")
)

# =============================================================================
# MAIN FIGURE 1. Expression phenotype and cross-cohort stage context
# =============================================================================

tcga_meta <- read_contract(
  "Part1", "Part1_TCGA_sample_metadata.csv",
  c("patient_id", "sample_id", "sample_type", "olfml2b_expression",
    "olfml2b_z", "is_tumor", "is_normal", "stage")
)
tn_summary <- read_contract(
  "Part3", "03a_tcga_olfml2b_tumor_normal_summary.csv",
  c("n_tumor", "n_normal", "median_tumor", "median_normal", "wilcox_p",
    "n_paired", "paired_median_delta", "paired_wilcox_p")
)
paired_values <- read_contract(
  "Part3", "03b_tcga_olfml2b_paired_values.csv",
  c("patient_id", "tumor_expression", "normal_expression", "delta_tumor_minus_normal")
)

mf1a_expression <- tcga_meta %>%
  filter(is_tumor %in% TRUE | is_normal %in% TRUE, is.finite(olfml2b_expression)) %>%
  transmute(
    patient_id, sample_id,
    tissue = factor(if_else(is_tumor %in% TRUE, "Tumour", "Normal"),
                    levels = c("Normal", "Tumour")),
    olfml2b_expression
  )
write_panel_source(
  mf1a_expression, "MAIN", "MF1a_TCGA_tumour_normal_distribution.csv", "Fig1a",
  "Part1/Part1_TCGA_sample_metadata.csv; Part3/03a_tcga_olfml2b_tumor_normal_summary.csv",
  "Tumour and normal samples retained; no resampling or thresholding",
  "Sample for unpaired comparison", "Differential abundance only"
)
tn <- tn_summary[1, ]
unpaired_auc <- pairwise_auc(
  mf1a_expression$olfml2b_expression[mf1a_expression$tissue == "Tumour"],
  mf1a_expression$olfml2b_expression[mf1a_expression$tissue == "Normal"]
)
tn_display_stats <- tn_summary %>%
  mutate(
    probability_tumour_greater_than_normal = unpaired_auc,
    rank_biserial_effect = 2 * unpaired_auc - 1
  )
write_panel_source(
  tn_display_stats, "MAIN", "MF1a_TCGA_tumour_normal_statistics.csv", "Fig1a-statistics",
  "Part3/03a_tcga_olfml2b_tumor_normal_summary.csv; Part1/Part1_TCGA_sample_metadata.csv",
  "Frozen tests retained; probability-of-superiority and rank-biserial effect added from plotted values",
  "Sample/patient as specified upstream", "Expression difference only"
)
p1a <- ggplot(mf1a_expression, aes(x = tissue, y = olfml2b_expression, fill = tissue)) +
  geom_violin(trim = FALSE, width = 0.82, alpha = 0.58, colour = NA) +
  geom_point(
    position = position_jitter(width = 0.075, height = 0, seed = PART9_SEED),
    size = 0.65, alpha = 0.22, colour = COL$black
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23,
    size = 3.0, fill = "white", colour = COL$black, stroke = 0.55
  ) +
  scale_fill_manual(values = c(Normal = COL$blue, Tumour = COL$red), guide = "none") +
  labs(
    title = "Tumour-normal expression",
    subtitle = paste0(
      "TCGA-STAD; normal n=", tn$n_normal, ", tumour n=", tn$n_tumor,
      "; medians ", formatC(tn$median_normal, digits = 2, format = "f"),
      " vs ", formatC(tn$median_tumor, digits = 2, format = "f"),
      "\nrank-biserial=", formatC(2 * unpaired_auc - 1, digits = 2, format = "f"),
      "; Wilcoxon P = ", fmt_p(tn$wilcox_p)
    ),
    x = NULL, y = "OLFML2B expression"
  ) + theme_pub() +
  theme(plot.subtitle = element_text(colour = COL$grey, size = 6.4, lineheight = 1.08))

mf1b_paired <- paired_values %>%
  filter(is.finite(delta_tumor_minus_normal)) %>%
  arrange(delta_tumor_minus_normal) %>%
  mutate(
    case_order = row_number(),
    direction = factor(
      if_else(delta_tumor_minus_normal >= 0, "Increase", "Decrease"),
      levels = c("Decrease", "Increase")
    )
  )
write_panel_source(
  mf1b_paired, "MAIN", "MF1b_TCGA_patient_paired_values.csv", "Fig1b",
  "Part3/03b_tcga_olfml2b_paired_values.csv; Part3/03a_tcga_olfml2b_tumor_normal_summary.csv",
  "All eligible paired patients retained; signed histogram uses frozen paired deltas",
  "Patient", "Paired directional change"
)
paired_delta_range <- diff(range(mf1b_paired$delta_tumor_minus_normal, na.rm = TRUE))
paired_binwidth <- if (is.finite(paired_delta_range) && paired_delta_range > 0) {
  paired_delta_range / 12
} else {
  0.25
}
paired_n_increase <- sum(mf1b_paired$direction == "Increase")
paired_n_decrease <- sum(mf1b_paired$direction == "Decrease")
paired_median_ci <- bootstrap_median_ci(mf1b_paired$delta_tumor_minus_normal)
p1b <- ggplot(mf1b_paired, aes(x = delta_tumor_minus_normal, fill = direction)) +
  geom_histogram(
    binwidth = paired_binwidth, boundary = 0, closed = "left",
    colour = "white", linewidth = 0.35, alpha = 0.92
  ) +
  geom_rug(aes(colour = direction), sides = "b", alpha = 0.72, linewidth = 0.35) +
  annotate(
    "point", x = tn$paired_median_delta, y = 0,
    shape = 23, size = 3.0, fill = COL$gold, colour = "white", stroke = 0.55
  ) +
  scale_fill_manual(values = c(Decrease = COL$blue, Increase = COL$red), guide = "none") +
  scale_colour_manual(values = c(Decrease = COL$blue, Increase = COL$red), guide = "none") +
  labs(
    title = "Paired change distribution",
    subtitle = paste0(
      "n=", tn$n_paired, "; median Δ=",
      formatC(tn$paired_median_delta, digits = 2, format = "f"),
      " [bootstrap 95% CI ", formatC(paired_median_ci[1], digits = 2, format = "f"),
      ", ", formatC(paired_median_ci[2], digits = 2, format = "f"), "]\n",
      paired_n_increase, " increased / ", paired_n_decrease, " decreased; P = ",
      fmt_p(tn$paired_wilcox_p)
    ),
    x = "Tumour minus normal OLFML2B expression", y = "Paired patients", fill = NULL
  ) + theme_pub() +
  theme(plot.subtitle = element_text(colour = COL$grey, size = 6.4, lineheight = 1.08))

extract_stage_source <- function(dat, cohort) {
  if ("is_tumor" %in% names(dat)) {
    keep <- is.na(dat$is_tumor) | dat$is_tumor %in% TRUE
    dat <- dat[keep, , drop = FALSE]
  }
  stage_col <- intersect(c("stage_numeric", "stage", "stage_analysis_primary"), names(dat))[1]
  if (is.na(stage_col)) stop("No stage column for ", cohort)
  s <- toupper(trimws(as.character(dat[[stage_col]])))
  stage_numeric <- suppressWarnings(as.integer(s))
  stage_numeric[is.na(stage_numeric) & str_detect(s, "IV")] <- 4L
  stage_numeric[is.na(stage_numeric) & str_detect(s, "III")] <- 3L
  stage_numeric[is.na(stage_numeric) & str_detect(s, "II")] <- 2L
  stage_numeric[
    is.na(stage_numeric) & str_detect(s, "(^|[^A-Z])I([^A-Z]|$)")
  ] <- 1L
  tibble(
    cohort = cohort,
    patient_id = if ("patient_id" %in% names(dat)) dat$patient_id else dat$sample_id,
    stage = stage_numeric,
    olfml2b_z = dat$olfml2b_z
  ) %>% filter(stage %in% 1:4, is.finite(olfml2b_z))
}

gse62254_meta <- read_contract(
  "Part2", "Part2_GSE62254_sample_metadata.csv", c("patient_id", "stage_numeric", "olfml2b_z")
)
gse15459_meta <- read_contract(
  "Part2", "Part2_GSE15459_sample_metadata.csv", c("patient_id", "stage_numeric", "olfml2b_z")
)
gse26253_meta <- read_contract(
  "Part2", "Part2_GSE26253_sample_metadata.csv", c("patient_id", "stage_numeric", "olfml2b_z")
)
mf1c_stage <- dplyr::bind_rows(
  extract_stage_source(tcga_meta, "TCGA_STAD"),
  extract_stage_source(gse62254_meta, "GSE62254"),
  extract_stage_source(gse15459_meta, "GSE15459"),
  extract_stage_source(gse26253_meta, "GSE26253")
) %>% mutate(
  cohort = factor(cohort_label(cohort), levels = c("TCGA-STAD", "GSE62254", "GSE15459", "GSE26253")),
  stage_label = factor(stage, levels = 1:4, labels = c("I", "II", "III", "IV"))
)
stage_global <- mf1c_stage %>%
  group_by(cohort) %>%
  summarise(
    kw_p = tryCatch(
      stats::kruskal.test(olfml2b_z ~ stage_label)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(kw_fdr = p.adjust(kw_p, method = "BH"))
stage_counts <- mf1c_stage %>%
  count(cohort, stage_label, name = "stage_n")
mf1c_stage <- mf1c_stage %>%
  left_join(stage_global, by = "cohort") %>%
  left_join(stage_counts, by = c("cohort", "stage_label")) %>%
  mutate(
    cohort_panel = paste0(as.character(cohort), "\nKW P = ", fmt_p(kw_p),
                          "; FDR = ", fmt_p(kw_fdr)),
    stage_axis_label = paste0(as.character(stage_label), "\nn=", stage_n)
  )
stage_panel_levels <- mf1c_stage %>%
  distinct(cohort, cohort_panel) %>%
  arrange(cohort) %>%
  pull(cohort_panel)
stage_axis_levels <- mf1c_stage %>%
  distinct(cohort, stage, stage_axis_label) %>%
  arrange(cohort, stage) %>%
  pull(stage_axis_label)
mf1c_stage <- mf1c_stage %>% mutate(
  cohort_panel = factor(cohort_panel, levels = stage_panel_levels),
  stage_axis_label = factor(stage_axis_label, levels = unique(stage_axis_levels))
)
write_panel_source(
  mf1c_stage, "MAIN", "MF1c_stage_patient_values.csv", "Fig1c",
  paste(
    "Part1/Part1_TCGA_sample_metadata.csv;",
    "Part2/Part2_GSE62254_sample_metadata.csv;",
    "Part2/Part2_GSE15459_sample_metadata.csv;",
    "Part2/Part2_GSE26253_sample_metadata.csv"
  ),
  "Stage I-IV patients retained; within-cohort OLFML2B z scores; all plotted patients written; promoted from Supplementary Fig. S1a",
  "Patient", "Cross-sectional stage context; not evidence of progression"
)
p1c <- ggplot(mf1c_stage, aes(x = stage_axis_label, y = olfml2b_z)) +
  geom_hline(yintercept = 0, colour = COL$light, linewidth = 0.45) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = PART9_SEED),
    size = 0.62, alpha = 0.20, colour = COL$grey
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23, size = 2.75,
    fill = COL$gold, colour = "white", stroke = 0.55
  ) +
  facet_wrap(~cohort_panel, nrow = 1, scales = "free_x") +
  labs(
    title = "Stage context across cohorts",
    subtitle = "Patients and group medians; global Kruskal-Wallis tests are descriptive; no monotonic trend is assumed",
    x = "Pathological stage", y = "Within-cohort OLFML2B z score"
  ) + theme_pub(7.2) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.05)))

fig1_top <- (p1a | p1b) + plot_layout(widths = c(0.44, 0.56))
fig1 <- fig1_top / p1c +
  plot_layout(heights = c(0.88, 1.0)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) &
  main_panel_sanitizer
save_figure_bundle(
  fig1, "MAIN", "MAIN_FIGURE_01_CNS_v8_10_0",
  FINAL_MAIN_WIDTH_IN, 5.25
)

# =============================================================================
# MAIN FIGURE 2. Cross-cohort stromal ecology
# =============================================================================

tme_corr <- read_contract(
  "Part4", "09_olfml2b_tme_correlations_by_cohort.csv",
  c("cohort", "signature", "n", "rho", "p_value", "fdr_global", "analysis_population")
)
tme_meta <- read_contract(
  "Part4", "10_olfml2b_tme_meta_correlations.csv",
  c("signature", "k", "n_total", "rho_meta", "rho_lcl", "rho_ucl",
    "prediction_low", "prediction_high", "I2_approx", "fdr")
)
tme_shift <- read_contract(
  "Part4", "11_olfml2b_high_low_tme_shift.csv",
  c("cohort", "signature", "n", "high_minus_low", "fdr", "analysis_population")
)
tme_boot <- read_contract(
  "Part4", "12c_tme_bootstrap_stability_audit.csv",
  c("cohort", "endpoint", "attenuation_axis", "attenuation_bootstrap_lcl",
    "attenuation_bootstrap_ucl", "attenuation_bootstrap_median", "bootstrap_B")
)
tme_scores <- read_contract(
  "Part4", "08_sample_level_tme_scores_and_axes.csv",
  c("cohort", "OLFML2B_z", "CAF_Core", "ECM_Remodeling", "TGFb_Response", "analysis_population"),
  selected_columns = c(
    "cohort", "OLFML2B_z", "CAF_Core", "ECM_Remodeling",
    "TGFb_Response", "analysis_population"
  ),
  col_types = readr::cols(
    cohort = readr::col_character(), OLFML2B_z = readr::col_double(),
    CAF_Core = readr::col_double(), ECM_Remodeling = readr::col_double(),
    TGFb_Response = readr::col_double(), analysis_population = readr::col_character()
  )
)

mf2a <- tme_scores %>%
  filter(analysis_population == "TUMOR_ONLY", is.finite(OLFML2B_z)) %>%
  select(cohort, OLFML2B_z, CAF_Core, ECM_Remodeling, TGFb_Response) %>%
  pivot_longer(c(CAF_Core, ECM_Remodeling, TGFb_Response), names_to = "axis", values_to = "score") %>%
  filter(is.finite(score)) %>%
  group_by(cohort, axis) %>%
  mutate(decile = ntile(OLFML2B_z, 10L)) %>%
  group_by(cohort, axis, decile) %>%
  summarise(x = median(OLFML2B_z), y = median(score), n = n(), .groups = "drop") %>%
  mutate(axis_label = feature_label(axis))
write_panel_source(
  mf2a, "MAIN", "MF2a_sample_level_decile_display.csv", "Fig2a",
  "Part4/08_sample_level_tme_scores_and_axes.csv",
  "Within-cohort decile medians generated for display; formal tests remain continuous",
  "Sample within cohort", "Descriptive shape; not a binned inferential model"
)
p2a <- ggplot(mf2a, aes(x = x, y = y, colour = cohort)) +
  geom_point(size = 1.85, alpha = 0.78) +
  facet_wrap(~axis_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = COHORT_COLORS) +
  labs(x = "OLFML2B within-cohort z score", y = "Programme score", colour = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_pub(7.0) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.8),
    axis.title.y = element_text(margin = margin(r = 0)),
    axis.text.y = element_text(margin = margin(r = 1.0)),
    panel.spacing.x = grid::unit(2.2, "mm"),
    plot.margin = margin(4, 4, 5, 1)
  )

main_programs <- c(
  "CAF_Core", "ECM_Remodeling", "TGFb_Response",
  "Checkpoint_Exhaustion", "CD8_Cytotoxic", "Proliferation_Control"
)
main_cohorts <- c("TCGA_STAD", "GSE84437", "GSE62254", "GSE15459", "GSE26253")
mf2b <- tme_corr %>%
  filter(
    analysis_population == "TUMOR_ONLY", signature %in% main_programs,
    cohort %in% main_cohorts, is.finite(rho)
  ) %>%
  mutate(
    cohort = factor(cohort_label(cohort), levels = cohort_label(main_cohorts)),
    program = factor(feature_label(signature), levels = rev(feature_label(main_programs))),
    mark = if_else(!is.na(fdr_global) & fdr_global < 0.05, "*", "")
  )
write_panel_source(
  mf2b, "MAIN", "MF2b_cross_cohort_ecology.csv", "Fig2b",
  "Part4/09_olfml2b_tme_correlations_by_cohort.csv",
  "Prespecified six-program display; tumour samples only; no value re-estimation",
  "Patient/sample within cohort; cohort for replication", "Ecological association; not causality"
)
p2b <- ggplot(mf2b, aes(x = cohort, y = program)) +
  geom_point(aes(size = abs(rho), colour = rho), alpha = 0.92) +
  geom_text(data = mf2b %>% filter(abs(rho) >= 0.28), aes(label = mark),
            colour = "white", size = 2.2, fontface = "bold") +
  geom_text(data = mf2b %>% filter(abs(rho) < 0.28), aes(label = mark),
            colour = COL$black, size = 2.2, fontface = "bold") +
  scale_size_continuous(
    range = c(1.8, 6.2), limits = c(0, 1), breaks = c(0.25, 0.50, 0.75)
  ) +
  scale_colour_gradient2(low = COL$blue, mid = "white", high = COL$red, midpoint = 0, limits = c(-1, 1)) +
  labs(
    title = "Cross-cohort ecological map",
    subtitle = "Circle area is |Spearman ρ|; * global FDR<0.05",
    x = NULL, y = NULL, size = "|ρ|", colour = "ρ"
  ) +
  guides(
    size = guide_legend(direction = "horizontal", nrow = 1, title.position = "left"),
    colour = guide_colourbar(
      direction = "horizontal", title.position = "left",
      barwidth = grid::unit(18, "mm"), barheight = grid::unit(2.2, "mm")
    )
  ) +
  theme_pub(7.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_line(colour = "#F0F0F0"),
    legend.position = "bottom", legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = grid::unit(0.2, "mm"),
    legend.text = element_text(size = 5.4)
  )

mf2c_meta <- tme_meta %>%
  filter(signature %in% main_programs) %>%
  transmute(
    signature, program = feature_label(signature), k, n_total,
    pooled_rho = rho_meta, ci_low = rho_lcl, ci_high = rho_ucl,
    prediction_low, prediction_high, i2 = I2_approx, fdr
  )
mf2c_cohort <- tme_corr %>%
  filter(
    analysis_population == "TUMOR_ONLY", signature %in% main_programs,
    cohort %in% main_cohorts, is.finite(rho)
  ) %>%
  transmute(signature, program = feature_label(signature), cohort, rho)
mf2c_source <- dplyr::bind_rows(
  mf2c_cohort %>% mutate(row_type = "cohort"),
  mf2c_meta %>% mutate(row_type = "meta")
)
write_panel_source(
  mf2c_source, "MAIN", "MF2c_pooled_ecology_nested_intervals.csv", "Fig2c",
  "Part4/09_olfml2b_tme_correlations_by_cohort.csv; Part4/10_olfml2b_tme_meta_correlations.csv",
  "Cohort estimates displayed with frozen random-effects summary, CI and prediction interval",
  "Cohort", "Cross-cohort ecological association; heterogeneity retained"
)
mf2c_meta <- mf2c_meta %>% mutate(
  program = factor(program, levels = rev(feature_label(main_programs))),
  pi_crosses_zero = prediction_low <= 0 & prediction_high >= 0,
  summary_label = paste0(
    "ρ ", formatC(pooled_rho, digits = 2, format = "f"),
    " [", formatC(ci_low, digits = 2, format = "f"), ", ",
    formatC(ci_high, digits = 2, format = "f"), "] · I² ", round(i2), "%"
  )
)
mf2c_cohort <- mf2c_cohort %>% mutate(
  program = factor(program, levels = rev(feature_label(main_programs)))
)
p2c <- ggplot() +
  geom_tile(
    data = mf2c_meta,
    aes(
      x = (prediction_low + prediction_high) / 2,
      y = program,
      width = pmax(prediction_high - prediction_low, .Machine$double.eps)
    ),
    height = 0.58, fill = "#E3E5E8", colour = NA
  ) +
  geom_tile(
    data = mf2c_meta,
    aes(
      x = (ci_low + ci_high) / 2,
      y = program,
      width = pmax(ci_high - ci_low, .Machine$double.eps)
    ),
    height = 0.22, fill = COL$grey, colour = NA
  ) +
  geom_point(
    data = mf2c_cohort, aes(x = rho, y = program),
    position = position_jitter(height = 0.10, width = 0, seed = PART9_SEED),
    size = 1.45, colour = COL$grey, alpha = 0.65
  ) +
  geom_point(
    data = mf2c_meta, aes(x = pooled_rho, y = program),
    shape = 23, size = 3.3, fill = COL$red, colour = "white", stroke = 0.45
  ) +
  geom_text(
    data = mf2c_meta,
    aes(x = 2.02, y = program, label = summary_label),
    hjust = 1, size = 1.48, colour = COL$black
  ) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  labs(
    title = "Heterogeneous pooled ecological effects",
    subtitle = paste0(
      "Dark blocks: 95% CI; pale blocks: prediction interval; ",
      sum(mf2c_meta$pi_crosses_zero), "/", nrow(mf2c_meta),
      " PIs cross 0; I² range ", round(min(mf2c_meta$i2)), "-",
      round(max(mf2c_meta$i2)), "%"
    ),
    x = "Spearman ρ", y = NULL
  ) +
  scale_x_continuous(
    limits = c(-0.75, 2.08),
    breaks = c(-0.5, 0, 0.5, 1.0),
    expand = expansion(mult = 0.01)
  ) +
  theme_pub(7.2)

composite_programs <- c("CAF_TGFb_axis", "Immune_Exclusion_Index", "Suppressive_TME_Index")
mf2d <- tme_shift %>%
  filter(
    analysis_population == "TUMOR_ONLY", signature %in% composite_programs,
    cohort %in% main_cohorts, is.finite(high_minus_low)
  ) %>%
  mutate(
    cohort_label = cohort_label(cohort),
    program = factor(feature_label(signature), levels = rev(feature_label(composite_programs)))
  ) %>%
  group_by(program) %>%
  mutate(
    cross_cohort_median = median(high_minus_low, na.rm = TRUE),
    support_count = sum(is.finite(fdr) & fdr < 0.05 & high_minus_low > 0),
    eligible_count = sum(is.finite(high_minus_low))
  ) %>%
  ungroup()
write_panel_source(
  mf2d, "MAIN", "MF2d_composite_state_shift.csv", "Fig2d",
  "Part4/11_olfml2b_high_low_tme_shift.csv",
  "Only prespecified composite axes not displayed in Fig2a; cohort estimates plus cross-cohort median",
  "Patient/sample within cohort", "Display split only; not a prognostic threshold"
)
p2d <- ggplot(mf2d, aes(x = high_minus_low, y = program)) +
  annotate("rect", xmin = -0.018, xmax = 0.018, ymin = -Inf, ymax = Inf,
           fill = "#EFEFEF", colour = NA) +
  geom_point(aes(colour = cohort), size = 1.9, alpha = 0.80,
             position = position_jitter(height = 0.08, width = 0, seed = PART9_SEED)) +
  geom_point(
    data = mf2d %>% distinct(program, cross_cohort_median),
    aes(x = cross_cohort_median, y = program), inherit.aes = FALSE,
    shape = 23, size = 3.4, fill = COL$black, colour = "white"
  ) +
  geom_text(
    data = mf2d %>% distinct(program, support_count, eligible_count) %>%
      mutate(label_x = max(mf2d$high_minus_low, na.rm = TRUE) + 0.012),
    aes(x = label_x, y = program, label = paste0(support_count, "/", eligible_count, " q<0.05")),
    inherit.aes = FALSE, hjust = 0, size = 1.65, colour = COL$grey
  ) +
  scale_colour_manual(values = COHORT_COLORS) +
  labs(
    title = "Composite state shift",
    subtitle = "OLFML2B-high minus low median; diamonds are cross-cohort medians",
    x = "Programme-score difference", y = NULL, colour = NULL
  ) +
  guides(colour = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  theme_pub(7.2) +
  theme(plot.margin = margin(5, 4, 5, 5)) +
  coord_cartesian(clip = "on")

primary_axes <- c("CAF_Core", "ECM_Remodeling", "TGFb_Response")
sf2c_attenuation <- tme_boot %>%
  filter(
    attenuation_axis %in% primary_axes,
    is.finite(attenuation_bootstrap_median),
    is.finite(attenuation_bootstrap_lcl),
    is.finite(attenuation_bootstrap_ucl)
  ) %>%
  mutate(
    cohort_endpoint = paste(cohort_label(cohort), endpoint),
    axis = factor(feature_label(attenuation_axis), levels = feature_label(primary_axes)),
    interval_crosses_zero = attenuation_bootstrap_lcl <= 0 & attenuation_bootstrap_ucl >= 0,
    interval_status = factor(
      if_else(interval_crosses_zero, "95% interval crosses 0", "95% interval excludes 0"),
      levels = c("95% interval crosses 0", "95% interval excludes 0")
    ),
    attenuation_label = paste0(
      if_else(
        abs(attenuation_bootstrap_median) >= 1000,
        paste0(if_else(attenuation_bootstrap_median < 0, "−", ""), "≥1000"),
        formatC(attenuation_bootstrap_median, digits = 0, format = "f")
      ),
      "%"
    )
  )
write_panel_source(
  sf2c_attenuation, "SUPPLEMENTARY", "SF2c_bootstrap_attenuation.csv", "FigS2c",
  "Part4/12c_tme_bootstrap_stability_audit.csv",
  "All eligible patient-resampling intervals for the three primary axes",
  "Patient resampling within cohort/endpoint",
  "Uncertainty of shared prognostic information; not causal mediation"
)
p_s2c_attenuation <- ggplot(
  sf2c_attenuation,
  aes(x = cohort_endpoint, y = axis, fill = interval_status)
) +
  geom_tile(colour = "white", linewidth = 0.58) +
  geom_text(aes(label = attenuation_label), size = 2.25, colour = COL$black) +
  scale_fill_manual(
    values = c(
      "95% interval crosses 0" = "#E1E3E6",
      "95% interval excludes 0" = "#D95A4E"
    )
  ) +
  labs(
    title = "Uncertain prognostic attenuation",
    subtitle = paste0(
      sum(sf2c_attenuation$interval_crosses_zero), "/", nrow(sf2c_attenuation),
      " patient-bootstrap intervals cross 0; cells show median attenuation; not mediation"
    ),
    x = "Cohort and endpoint",
    y = NULL, fill = NULL
  ) +
  theme_pub(6.7) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom", legend.text = element_text(size = 5.8)
  )

top_row_2_core <- (p2a | p2b) +
  plot_layout(widths = c(0.55, 0.45))
top_row_2 <- patchwork::free(
  top_row_2_core, side = "l", type = "panel"
)
bottom_row_2 <- (p2c | p2d) + plot_layout(widths = c(0.56, 0.44))
fig2 <- top_row_2 / bottom_row_2 +
  plot_layout(heights = c(1.05, 0.90)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) &
  main_panel_sanitizer
save_figure_bundle(
  fig2, "MAIN", "MAIN_FIGURE_02_CNS_v8_10_0",
  FINAL_MAIN_WIDTH_IN, 5.35
)

# =============================================================================
# MAIN FIGURE 3. Single-cell localization and paired patient evidence
# =============================================================================

sc_localization <- read_contract(
  "Part6", "54_strict_localization_by_dataset.csv",
  c("cohort", "marker_celltype", "n_units", "total_cells", "total_positive",
    "pooled_positive_fraction", "eligible", "threshold_top3_consistent",
    "high_confidence_top3_consistent")
)
sc_paired <- read_contract(
  "Part6", "41_primary_paired_patient_tumor_normal_delta.csv",
  c("cohort", "patient", "tumor_fraction", "normal_fraction", "delta_tumor_minus_normal")
)
sc_paired_test <- read_contract(
  "Part6", "42_primary_paired_patient_tumor_normal_test.csv",
  c("cohort", "n_paired_patients", "median_delta_tumor_minus_normal", "wilcox_p",
    "inference_unit", "wilcox_p_fdr")
)
sc_bridge <- read_contract(
  "Part6", "65_fulltx_paired_change_program_bridge_exact.csv",
  c("cohort", "contrast", "program", "n_paired_patients", "rho_delta", "p_value",
    "formal_significance_eligible", "fdr_within_contrast")
)
sc_evidence <- read_contract(
  "Part6", "70_final_cell_source_evidence_matrix.csv",
  c("marker_celltype", "n_datasets_eligible", "n_datasets_top3",
    "leave_one_dataset_out_robust", "evidence_grade", "interpretation")
)
sc_program <- read_contract(
  "Part6", "32_OLFML2B_positive_marker_program_tests.csv",
  c("cohort", "marker_program", "n_samples", "median_delta_positive_minus_negative", "wilcox_p_fdr")
)
sc_loo <- read_contract(
  "Part6", "66_strict_localization_leave_one_dataset_out.csv",
  c("marker_celltype", "omitted_dataset", "n_datasets_remaining", "n_remaining_top3",
    "replicated_after_omission")
)

cell_order <- c(
  "Fibroblast", "Myofibroblast", "Pericyte", "Myeloid",
  "Epithelial", "Endothelial", "T_NK", "B_Plasma"
)
mf3a <- sc_localization %>%
  filter(marker_celltype %in% cell_order, eligible %in% TRUE, is.finite(pooled_positive_fraction)) %>%
  left_join(
    sc_evidence %>% select(marker_celltype, evidence_grade),
    by = "marker_celltype"
  ) %>%
  mutate(
    celltype = factor(feature_label(marker_celltype), levels = rev(feature_label(cell_order))),
    cohort = factor(cohort, levels = c("GSE134520", "GSE150290", "GSE167297", "GSE183904")),
    top3 = threshold_top3_consistent %in% TRUE | high_confidence_top3_consistent %in% TRUE
  )
write_panel_source(
  mf3a, "MAIN", "MF3a_dataset_localization.csv", "Fig3a",
  "Part6/54_strict_localization_by_dataset.csv; Part6/70_final_cell_source_evidence_matrix.csv",
  "Eligible dataset-celltype rows only; no pooling across datasets in the display",
  "Official sample/patient within dataset", "Dominant localization, not cell exclusivity"
)
p3a <- ggplot(mf3a, aes(x = pooled_positive_fraction, y = celltype)) +
  geom_point(aes(size = total_cells, colour = top3), alpha = 0.88) +
  facet_wrap(~cohort, ncol = 2, scales = "free_y") +
  scale_x_log10(labels = label_percent(accuracy = 0.01)) +
  scale_colour_manual(
    values = c(`FALSE` = "#A8A8A8", `TRUE` = COL$purple),
    labels = c(`FALSE` = "Other", `TRUE` = "Top-three")
  ) +
  scale_size_continuous(
    range = c(1.4, 3.8),
    breaks = c(10000, 30000, 50000),
    labels = c("10k", "30k", "50k")
  ) +
  labs(
    title = "Dataset-level OLFML2B detection",
    subtitle = "Position is detection fraction; area is analysed cells, not expression; purple denotes top-three support",
    x = "OLFML2B-positive fraction (log scale)", y = NULL,
    colour = "Top-three", size = "Cells"
  ) +
  guides(
    size = guide_legend(order = 1, nrow = 1, title.position = "left"),
    colour = guide_legend(
      order = 2, nrow = 1, title.position = "left",
      override.aes = list(size = 2.2)
    )
  ) +
  theme_pub(7.0) +
  scale_y_discrete(drop = TRUE) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.spacing = grid::unit(0.4, "mm"),
    legend.margin = margin(t = 1, r = 0, b = 0, l = 0),
    plot.margin = margin(4, 4, 5, 5)
  )

mf3b <- sc_paired %>%
  left_join(sc_paired_test, by = "cohort") %>%
  mutate(
    comparison = recode(cohort, GSE150290 = "tumour-adjacent", GSE183904 = "tumour-normal",
                        .default = "paired contrast"),
    cohort_label = paste0(
      cohort, "\nn=", n_paired_patients,
      "; q=", vapply(wilcox_p_fdr, fmt_p, character(1))
    )
  )
write_panel_source(
  mf3b, "MAIN", "MF3b_patient_paired_pseudobulk.csv", "Fig3b",
  "Part6/41_primary_paired_patient_tumor_normal_delta.csv; Part6/42_primary_paired_patient_tumor_normal_test.csv",
  "Every paired patient retained; frozen cohort-level test appended",
  "Patient", "Paired pseudobulk abundance; no cell-level inference"
)
p3b <- ggplot(mf3b, aes(x = cohort_label, y = delta_tumor_minus_normal)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.002, ymax = 0.002,
           fill = "#EFEFEF", colour = NA) +
  geom_point(
    aes(colour = delta_tumor_minus_normal >= 0),
    position = position_jitter(width = 0.08, height = 0, seed = PART9_SEED),
    size = 1.7, alpha = 0.82
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23,
    fill = "white", colour = COL$black, size = 3.2, stroke = 0.65
  ) +
  scale_colour_manual(values = c(`FALSE` = COL$blue, `TRUE` = COL$red), guide = "none") +
  labs(
    title = "Paired patient pseudobulk",
    subtitle = "Points are patients; open diamonds denote cohort medians",
    x = NULL, y = "Paired Δ positive fraction"
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.005, base = 10),
    breaks = c(-0.10, -0.03, -0.01, 0, 0.01, 0.03, 0.10, 0.30)
  ) +
  theme_pub(7.2) +
  theme(
    axis.title.y = element_text(margin = margin(r = 0.5)),
    axis.text.x = element_text(size = 6.1, lineheight = 0.92),
    plot.margin = margin(4, 4, 5, 2)
  )

bridge_programs <- c(
  "CAF_ECM", "ECM_Remodeling", "TGFb_Response", "Myofibroblast",
  "CD8_Cytotoxic", "Epithelial_Differentiation"
)
mf3c_long <- sc_bridge %>%
  filter(
    cohort %in% c("GSE150290", "GSE183904"),
    program %in% bridge_programs, is.finite(rho_delta)
  ) %>%
  select(cohort, program, n_paired_patients, rho_delta, p_value, fdr_within_contrast)
mf3c <- mf3c_long %>%
  select(cohort, program, rho_delta) %>%
  pivot_wider(names_from = cohort, values_from = rho_delta) %>%
  filter(is.finite(GSE150290), is.finite(GSE183904)) %>%
  mutate(
    program_label = feature_label(program),
    direction = case_when(
      GSE150290 > 0 & GSE183904 > 0 ~ "Positive in both",
      GSE150290 < 0 & GSE183904 < 0 ~ "Negative in both",
      TRUE ~ "Direction discordant"
    )
  ) %>%
  left_join(
    mf3c_long %>%
      group_by(program) %>%
      summarise(
        n_cohorts_fdr_lt_005 = sum(is.finite(fdr_within_contrast) & fdr_within_contrast < 0.05),
        .groups = "drop"
      ),
    by = "program"
  )
coupling_n_150290 <- mf3c_long %>%
  filter(cohort == "GSE150290") %>% distinct(n_paired_patients) %>% pull(n_paired_patients)
coupling_n_183904 <- mf3c_long %>%
  filter(cohort == "GSE183904") %>% distinct(n_paired_patients) %>% pull(n_paired_patients)
if (length(coupling_n_150290) != 1L || length(coupling_n_183904) != 1L) {
  stop("Fig3c preflight failed: paired-patient n must be unique within each cohort.")
}
write_panel_source(
  mf3c, "MAIN", "MF3c_paired_change_program_coupling.csv", "Fig3c",
  "Part6/65_fulltx_paired_change_program_bridge_exact.csv",
  "Prespecified programmes with evaluable paired-change correlations in both paired cohorts",
  "Patient", "Within-patient ecological coupling; not mechanism"
)
mf3c_plot <- mf3c %>%
  mutate(program_order = fct_reorder(program_label, (GSE150290 + GSE183904) / 2)) %>%
  select(program, program_label, program_order, n_cohorts_fdr_lt_005, GSE150290, GSE183904) %>%
  pivot_longer(c(GSE150290, GSE183904), names_to = "cohort", values_to = "rho") %>%
  mutate(cohort = factor(cohort, levels = c("GSE150290", "GSE183904")))
p3c <- ggplot(mf3c_plot, aes(x = rho, y = program_order, colour = cohort, shape = cohort)) +
  annotate("rect", xmin = -0.018, xmax = 0.018, ymin = -Inf, ymax = Inf,
           fill = "#EFEFEF", colour = NA) +
  geom_point(size = 2.45, alpha = 0.90, position = position_dodge(width = 0.34)) +
  geom_text(
    data = mf3c %>% mutate(program_order = fct_reorder(program_label, (GSE150290 + GSE183904) / 2)),
    aes(x = pmax(GSE150290, GSE183904) + 0.035, y = program_order,
        label = paste0(n_cohorts_fdr_lt_005, "/2")),
    inherit.aes = FALSE, hjust = 0, size = 1.60, colour = COL$grey
  ) +
  scale_colour_manual(values = c(GSE150290 = COL$red, GSE183904 = COL$purple)) +
  scale_shape_manual(values = c(GSE150290 = 16, GSE183904 = 17)) +
  labs(
    title = "Paired target-programme coupling",
    subtitle = paste0(
      "Circles: GSE150290 n=", coupling_n_150290,
      "; triangles: GSE183904 n=", coupling_n_183904,
      "; cohort estimates are not pooled"
    ),
    x = "Within-patient change correlation (Spearman ρ)", y = NULL,
    colour = NULL, shape = NULL
  ) + theme_pub(7.2) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  theme(plot.margin = margin(6, 14, 6, 8), legend.position = "bottom") +
  coord_cartesian(clip = "off")

program_context <- c(
  "CAF_ECM", "ECM_Remodeling", "TGFb_Response", "Myofibroblast",
  "Myeloid_Macrophage", "CD8_Cytotoxic", "Epithelial_Differentiation"
)
mf3d <- sc_program %>%
  filter(cohort != "POOLED_PRIMARY_GC", marker_program %in% program_context) %>%
  mutate(
    cohort = factor(cohort, levels = c("GSE134520", "GSE150290", "GSE167297", "GSE183904")),
    program = factor(feature_label(marker_program), levels = rev(feature_label(program_context))),
    direction = median_delta_positive_minus_negative >= 0,
    mark = if_else(!is.na(wilcox_p_fdr) & wilcox_p_fdr < 0.05, "*", ""),
    label_text = mark,
    label_hjust = if_else(median_delta_positive_minus_negative >= 0, -0.08, 1.08)
  ) %>%
  group_by(cohort) %>%
  mutate(
    cohort_panel = if (n_distinct(n_samples) == 1L) {
      paste0(as.character(first(cohort)), "\nn=", first(n_samples))
    } else {
      as.character(first(cohort))
    }
  ) %>%
  ungroup() %>%
  mutate(
    cohort_panel = factor(
      cohort_panel,
      levels = unique(cohort_panel[order(as.integer(cohort))])
    )
  )
write_panel_source(
  mf3d, "MAIN", "MF3d_positive_cell_program_context.csv", "Fig3d",
  "Part6/32_OLFML2B_positive_marker_program_tests.csv",
  "Sample-level positive-minus-negative programme contrasts; pooled row excluded",
  "Official sample/patient", "State association; not cell-level DE or regulation"
)
p3d_limit <- ceiling(max(abs(mf3d$median_delta_positive_minus_negative), na.rm = TRUE) * 5) / 5
p3d <- ggplot(mf3d, aes(x = median_delta_positive_minus_negative, y = program, fill = direction)) +
  annotate("rect", xmin = -0.018, xmax = 0.018, ymin = -Inf, ymax = Inf,
           fill = "#EFEFEF", colour = NA) +
  geom_col(width = 0.64, alpha = 0.86) +
  geom_text(aes(label = label_text, hjust = label_hjust), size = 1.6, colour = COL$black) +
  facet_wrap(~cohort_panel, ncol = 2) +
  scale_fill_manual(values = c(`FALSE` = COL$blue, `TRUE` = COL$purple), guide = "none") +
  scale_x_continuous(limits = c(-p3d_limit * 1.18, p3d_limit * 1.18)) +
  labs(
    title = "OLFML2B-positive-cell programme context",
    subtitle = "Sample-level positive-minus-negative differences; * within-family FDR<0.05",
    x = "Median programme-score difference", y = NULL
  ) + theme_pub(6.4) + coord_cartesian(clip = "off")

loo_display_levels <- c(
  "All datasets", "GSE134520", "GSE150290", "GSE167297", "GSE183904"
)
mf3e_observed <- sc_loo %>%
  filter(n_datasets_remaining > 0) %>%
  mutate(
    omitted_dataset = if_else(omitted_dataset == "NONE", "All datasets", omitted_dataset),
    omitted_dataset = as.character(omitted_dataset)
  )
mf3e <- tidyr::expand_grid(
  marker_celltype = cell_order,
  omitted_dataset = loo_display_levels
) %>%
  left_join(mf3e_observed, by = c("marker_celltype", "omitted_dataset")) %>%
  mutate(
    evaluable = is.finite(n_datasets_remaining),
    display_size = if_else(evaluable, as.numeric(n_datasets_remaining), NA_real_),
    replication_status = factor(
      case_when(
        !evaluable ~ "Not evaluable",
        replicated_after_omission %in% TRUE ~ "Replicated",
        TRUE ~ "Not replicated"
      ),
      levels = c("Not evaluable", "Not replicated", "Replicated")
    ),
    ratio_label = if_else(
      evaluable,
      paste0(n_remaining_top3, "/", n_datasets_remaining),
      ""
    ),
    omitted_dataset = factor(omitted_dataset, levels = loo_display_levels),
    celltype = factor(
      feature_label(marker_celltype),
      levels = rev(feature_label(cell_order))
    )
  )
write_panel_source(
  mf3e, "MAIN", "MF3e_localization_LODO_matrix.csv", "Fig3e",
  "Part6/66_strict_localization_leave_one_dataset_out.csv",
  "Every compartment-by-omission combination retained; non-evaluable cells are shown explicitly as neutral crosses",
  "Dataset omission", "Localization robustness; not universal cell exclusivity"
)
p3e <- ggplot(mf3e, aes(x = omitted_dataset, y = celltype)) +
  geom_point(
    data = mf3e %>% filter(evaluable),
    aes(size = display_size, fill = replication_status),
    shape = 21, colour = "white", stroke = 0.42, alpha = 0.94
  ) +
  geom_text(
    data = mf3e %>% filter(!evaluable),
    aes(label = "×"), size = 2.35, colour = "#B6BAC1", fontface = "bold"
  ) +
  geom_text(
    data = mf3e %>% filter(evaluable),
    aes(label = ratio_label), size = 1.65, colour = COL$black
  ) +
  scale_fill_manual(
    values = c(
      "Not evaluable" = "#F3F4F6",
      "Not replicated" = "#D7D9DD",
      "Replicated" = "#A78BD4"
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(1.8, 6.2), limits = c(1, 4), breaks = 1:4, guide = "none"
  ) +
  labs(x = "Dataset omitted", y = NULL) +
  theme_pub(6.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

top_row_3 <- (p3a | p3b) + plot_layout(widths = c(0.52, 0.48))
bottom_row_3 <- (p3c | p3d) + plot_layout(widths = c(0.46, 0.54))
fig3 <- top_row_3 / bottom_row_3 / p3e +
  plot_layout(heights = c(1.0, 0.98, 0.64)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) &
  main_panel_sanitizer
save_figure_bundle(
  fig3, "MAIN", "MAIN_FIGURE_03_CNS_v8_10_0",
  FINAL_MAIN_WIDTH_IN, 6.85
)

# =============================================================================
# MAIN FIGURE 4. Spatial localization and source disambiguation
# =============================================================================

sp_detection <- read_contract(
  "Part7", "07_spatial_target_detection_by_section.csv",
  c("sample_id", "patient_id", "n_spots", "n_positive", "positive_fraction")
)
sp_spots <- read_contract(
  "Part7", "20_spot_level_continuous_scores.csv",
  c("spot_id", "sample_id", "patient_id", "pxl_row", "pxl_col", "in_tissue",
    "OLFML2B", "Fibroblast", "CAF_ECM")
)
sp_assoc <- read_contract(
  "Part7", "22_patient_level_continuous_associations.csv",
  c("effect_source", "feature", "n_sections", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "wilcoxon_fdr")
)
sp_ridge <- read_contract(
  "Part7", "41_multivariable_source_disambiguation_patient.csv",
  c("effect_source", "feature", "n_sections", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "exact_signflip_fdr", "wilcoxon_fdr")
)
sp_residual <- read_contract(
  "Part7", "27_technical_burden_residual_sensitivity_by_patient.csv",
  c("effect_source", "feature", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "wilcoxon_fdr")
)
sp_overlap <- read_contract(
  "Part7", "73_spatial_target_context_overlap_patient_summary.csv",
  c("effect_source", "feature", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "exact_signflip_fdr", "wilcoxon_fdr")
)

representative_section <- sp_detection %>%
  filter(is.finite(positive_fraction), n_positive > 0) %>%
  mutate(distance_to_median = abs(positive_fraction - median(positive_fraction))) %>%
  arrange(distance_to_median, sample_id) %>%
  slice(1) %>%
  pull(sample_id)
if (length(representative_section) != 1L) stop("Representative spatial section selection failed.")

mf4e_joint <- sp_spots %>%
  filter(
    sample_id == representative_section, in_tissue %in% TRUE,
    is.finite(OLFML2B), is.finite(CAF_ECM), is.finite(pxl_row), is.finite(pxl_col)
  ) %>%
  mutate(
    target_percentile = percent_rank(OLFML2B),
    context_percentile = percent_rank(CAF_ECM),
    target_high = target_percentile >= 0.75,
    context_high = context_percentile >= 0.75,
    joint_class = factor(
      case_when(
        target_high & context_high ~ "Joint high",
        target_high ~ "OLFML2B high",
        context_high ~ "CAF-ECM high",
        TRUE ~ "Neither high"
      ),
      levels = c("Neither high", "OLFML2B high", "CAF-ECM high", "Joint high")
    )
  ) %>%
  select(spot_id, sample_id, patient_id, pxl_row, pxl_col, target_percentile,
         context_percentile, target_high, context_high, joint_class)
write_panel_source(
  mf4e_joint, "MAIN", "MF4e_representative_joint_spatial.csv", "Fig4e",
  "Part7/20_spot_level_continuous_scores.csv",
  "Within-section upper quartiles define target-high, CAF-ECM-high and joint-high display classes",
  "Section for display", "Representative overlap; threshold is not a cell-identity rule"
)
p4e_joint <- ggplot() +
  geom_point(
    data = mf4e_joint %>% filter(joint_class == "Neither high"),
    aes(x = pxl_col, y = -pxl_row, colour = joint_class),
    size = 0.23, alpha = 0.55
  ) +
  geom_point(
    data = mf4e_joint %>% filter(joint_class == "CAF-ECM high"),
    aes(x = pxl_col, y = -pxl_row, colour = joint_class),
    size = 0.31, alpha = 0.74
  ) +
  geom_point(
    data = mf4e_joint %>% filter(joint_class == "OLFML2B high"),
    aes(x = pxl_col, y = -pxl_row, colour = joint_class),
    size = 0.35, alpha = 0.86
  ) +
  geom_point(
    data = mf4e_joint %>% filter(joint_class == "Joint high"),
    aes(x = pxl_col, y = -pxl_row, colour = joint_class),
    size = 0.42, alpha = 0.96
  ) +
  annotate(
    "label",
    x = min(mf4e_joint$pxl_col, na.rm = TRUE),
    y = max(-mf4e_joint$pxl_row, na.rm = TRUE),
    label = str_extract(representative_section, "^[0-9]{2}_[0-9]{5}"),
    hjust = 0, vjust = 1, size = 2.05, linewidth = 0,
    fill = scales::alpha("white", 0.78), colour = COL$black
  ) +
  scale_colour_manual(values = c(
    "Neither high" = "#E5E7EB", "OLFML2B high" = COL$red,
    "CAF-ECM high" = COL$purple, "Joint high" = "#32194F"
  ), labels = c(
    "Neither high" = "Neither", "OLFML2B high" = "OLFML2B",
    "CAF-ECM high" = "CAF-ECM", "Joint high" = "Joint"
  ), name = NULL, drop = FALSE) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  coord_equal() +
  theme_void(base_family = "sans") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.6),
        plot.margin = margin(6, 8, 6, 8))

mf4a <- sp_spots %>%
  filter(in_tissue %in% TRUE, is.finite(OLFML2B), is.finite(pxl_row), is.finite(pxl_col)) %>%
  group_by(sample_id) %>%
  mutate(
    within_section_percentile = percent_rank(OLFML2B),
    section_x_mid = mean(range(pxl_col, na.rm = TRUE)),
    section_y_mid = mean(range(pxl_row, na.rm = TRUE)),
    section_span = pmax(
      diff(range(pxl_col, na.rm = TRUE)),
      diff(range(pxl_row, na.rm = TRUE))
    ),
    display_x = (pxl_col - section_x_mid) / section_span,
    display_y = -(pxl_row - section_y_mid) / section_span,
    section_label = str_extract(sample_id, "^[0-9]{2}_[0-9]{5}")
  ) %>%
  ungroup() %>%
  select(
    spot_id, sample_id, section_label, patient_id,
    pxl_row, pxl_col, display_x, display_y, section_span,
    OLFML2B, within_section_percentile
  )
if (
  any(!is.finite(mf4a$section_span)) || any(mf4a$section_span <= 0) ||
  any(is.na(mf4a$section_label))
) {
  stop("Fig4a preflight failed: invalid section span or display label.")
}
write_panel_source(
  mf4a, "MAIN", "MF4a_all_section_spatial_atlas.csv", "Fig4a",
  "Part7/20_spot_level_continuous_scores.csv",
  "Every in-tissue spot retained; coordinates are centred and divided by one isotropic within-section span for display only; target is ranked within section",
  "Section for display; patient for inference elsewhere",
  "Observed cross-section distribution; normalized coordinates preserve shape but not absolute between-section size"
)
mf4a_labels <- mf4a %>%
  group_by(section_label) %>%
  summarise(
    label_x = mean(range(display_x, na.rm = TRUE)),
    label_y = max(display_y, na.rm = TRUE) + 0.018,
    .groups = "drop"
  )
p4a <- ggplot(mf4a, aes(x = display_x, y = display_y, colour = within_section_percentile)) +
  geom_point(size = 0.28, alpha = 0.86) +
  geom_text(
    data = mf4a_labels,
    aes(x = label_x, y = label_y, label = section_label),
    inherit.aes = FALSE, hjust = 0.5, vjust = 0,
    size = 1.95, fontface = "bold", colour = COL$black
  ) +
  facet_wrap(~section_label, ncol = 5, scales = "fixed") +
  scale_colour_gradient(low = "#ECEFF1", high = COL$red, limits = c(0, 1),
                        labels = label_percent(), name = "OLFML2B percentile") +
  guides(
    colour = guide_colourbar(
      direction = "horizontal", title.position = "left",
      barwidth = grid::unit(24, "mm"), barheight = grid::unit(2.2, "mm")
    )
  ) +
  coord_equal(xlim = c(-0.54, 0.54), ylim = c(-0.56, 0.60), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(
    strip.text = element_blank(),
    panel.spacing.x = grid::unit(0.8, "mm"),
    panel.spacing.y = grid::unit(1.0, "mm"),
    legend.position = "bottom",
    legend.title = element_text(size = 5.7),
    legend.text = element_text(size = 5.3),
    plot.margin = margin(1, 4, 2, 1)
  )

spatial_features <- c(
  "CAF_ECM", "ECM_Remodeling", "Fibroblast", "Myofibroblast",
  "TGFb_Response", "CD8_Cytotoxic", "Epithelial"
)
mf4b <- sp_assoc %>%
  filter(
    effect_source %in% c("same_spot_spearman", "neighbor_mean_spearman"),
    feature %in% spatial_features
  ) %>%
  select(effect_source, feature, n_patients, median_effect, bootstrap_ci_low,
         bootstrap_ci_high, wilcoxon_fdr) %>%
  pivot_wider(
    names_from = effect_source,
    values_from = c(n_patients, median_effect, bootstrap_ci_low, bootstrap_ci_high, wilcoxon_fdr)
  ) %>%
  mutate(
    feature_label = feature_label(feature),
    label_x = case_when(
      feature == "ECM_Remodeling" ~ 0.245,
      feature == "Fibroblast" ~ 0.215,
      feature == "CAF_ECM" ~ 0.325,
      feature == "Myofibroblast" ~ 0.205,
      feature == "TGFb_Response" ~ 0.055,
      feature == "CD8_Cytotoxic" ~ 0.030,
      TRUE ~ -0.055
    ),
    label_y = case_when(
      feature == "ECM_Remodeling" ~ 0.305,
      feature == "Fibroblast" ~ 0.275,
      feature == "CAF_ECM" ~ 0.238,
      feature == "Myofibroblast" ~ 0.205,
      feature == "TGFb_Response" ~ 0.175,
      feature == "CD8_Cytotoxic" ~ 0.120,
      TRUE ~ -0.073
    ),
    label_hjust = if_else(feature == "CAF_ECM", 0, 0.5)
  )
write_panel_source(
  mf4b, "MAIN", "MF4c_same_spot_neighbor_concordance.csv", "Fig4c",
  "Part7/22_patient_level_continuous_associations.csv",
  "Patient-level same-spot and neighbour summaries joined by prespecified feature",
  "Patient after within-patient section collapse", "Spatial association, not cell-cell interaction"
)
p4b <- ggplot(
  mf4b,
  aes(x = median_effect_same_spot_spearman, y = median_effect_neighbor_mean_spearman)
) +
  geom_point(aes(colour = feature %in% c("CD8_Cytotoxic", "Epithelial")), size = 2.7) +
  geom_text(
    data = mf4b,
    aes(x = label_x, y = label_y, label = feature_label, hjust = label_hjust),
    inherit.aes = FALSE, size = 1.92, colour = COL$black
  ) +
  scale_colour_manual(values = c(`FALSE` = COL$red, `TRUE` = COL$blue), guide = "none") +
  scale_x_continuous(limits = c(-0.09, 0.36), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-0.09, 0.325), expand = expansion(mult = 0.02)) +
  coord_equal(clip = "off") +
  labs(
    title = "Local and neighbouring context",
    subtitle = "Each point is one programme; bootstrap CIs remain in the source table",
    x = "Same-spot Spearman ρ", y = "Neighbour Spearman ρ"
  ) + theme_pub(7.3)

mf4c <- sp_overlap %>%
  filter(effect_source == "same_spot_target_context_overlap", is.finite(median_effect)) %>%
  mutate(
    context = fct_reorder(feature_label(feature), median_effect),
    detail_label = paste0("n=", n_patients, "; exact q=", fmt_p(exact_signflip_fdr))
  )
write_panel_source(
  mf4c, "MAIN", "MF4d_threshold_overlap.csv", "Fig4d",
  "Part7/73_spatial_target_context_overlap_patient_summary.csv",
  "All frozen target-context overlap summaries retained; uncertainty remains in the source table",
  "Patient after within-patient section collapse", "Descriptive overlap; not cell identity or interaction"
)
p4c <- ggplot(mf4c, aes(x = median_effect, y = context)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.65,
             fill = COL$teal, colour = "white", alpha = 0.92) +
  geom_text(
    aes(label = label_percent(accuracy = 1)(median_effect)),
    hjust = -0.35, size = 1.85, colour = COL$grey
  ) +
  scale_x_continuous(
    limits = c(0, max(mf4c$median_effect, na.rm = TRUE) * 1.20),
    labels = label_percent(accuracy = 1)
  ) +
  labs(x = "Median target-high/context-high overlap", y = NULL) +
  theme_pub(7.0) + coord_cartesian(clip = "off")

mf4d <- sp_ridge %>%
  filter(effect_source == "ridge_competing_source", is.finite(median_effect)) %>%
  mutate(
    feature_label = fct_reorder(feature_label(feature), median_effect),
    interval_excludes_zero = bootstrap_ci_low > 0 | bootstrap_ci_high < 0,
    exact_supported = !is.na(exact_signflip_fdr) & exact_signflip_fdr < 0.05
  )
write_panel_source(
  mf4d, "MAIN", "MF4f_competing_source_ridge.csv", "Fig4f",
  "Part7/41_multivariable_source_disambiguation_patient.csv",
  "All frozen competing-source coefficients retained",
  "Patient after within-patient section collapse", "Source disambiguation is associative"
)
p4d <- ggplot(mf4d, aes(x = median_effect, y = feature_label, fill = median_effect > 0)) +
  annotate("rect", xmin = -0.012, xmax = 0.012, ymin = -Inf, ymax = Inf,
           fill = "#EFEFEF", colour = NA) +
  geom_col(width = 0.62, alpha = 0.84) +
  scale_fill_manual(values = c(`FALSE` = COL$blue, `TRUE` = COL$red), guide = "none") +
  labs(
    title = "Competing-source disambiguation",
    subtitle = paste0(
      "n=", max(mf4d$n_patients, na.rm = TRUE),
      " patients; bootstrap intervals remain in the source table; minimum exact sign-flip FDR = ",
      fmt_p(min(mf4d$exact_signflip_fdr, na.rm = TRUE)), "; exploratory"
    ),
    x = "Standardized ridge coefficient", y = NULL
  ) + theme_pub(7.1) + coord_cartesian(clip = "off")

sf4d_technical <- mf4b %>%
  select(feature, raw_same_spot = median_effect_same_spot_spearman) %>%
  inner_join(
    sp_residual %>%
      filter(
        effect_source == "technical_residual_spearman",
        feature %in% spatial_features,
        is.finite(median_effect)
      ) %>%
      transmute(
        feature,
        residual_effect = median_effect,
        residual_ci_low = bootstrap_ci_low,
        residual_ci_high = bootstrap_ci_high,
        residual_fdr = wilcoxon_fdr
      ),
    by = "feature"
  ) %>%
  filter(is.finite(raw_same_spot), is.finite(residual_effect)) %>%
  mutate(
    feature_label = fct_reorder(feature_label(feature), residual_effect),
    retained_direction = sign(raw_same_spot) == sign(residual_effect)
  )
write_panel_source(
  sf4d_technical, "SUPPLEMENTARY", "SF4c_raw_vs_technical_residual.csv", "FigS4c",
  "Part7/22_patient_level_continuous_associations.csv; Part7/27_technical_burden_residual_sensitivity_by_patient.csv",
  "Raw same-spot patient-level medians joined to the prespecified technical-residual sensitivity estimate",
  "Patient after within-patient section collapse",
  "Measured technical-burden sensitivity; substantial attenuation limits robustness claims"
)
p_s4d_technical <- ggplot(sf4d_technical, aes(y = feature_label)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_segment(
    aes(x = raw_same_spot, xend = residual_effect, yend = feature_label),
    colour = "#A9A9A9", linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
  ) +
  geom_point(aes(x = raw_same_spot, shape = "Raw same-spot"), size = 2.25,
             fill = "white", colour = COL$grey, stroke = 0.55) +
  geom_point(aes(x = residual_effect, shape = "Technical residual",
                 colour = retained_direction), size = 2.35) +
  scale_colour_manual(
    values = c(`FALSE` = COL$blue, `TRUE` = COL$red),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Raw same-spot" = 21, "Technical residual" = 16)
  ) +
  labs(
    title = "Attenuation after technical residualization",
    subtitle = "Open circles: raw same-spot estimate; filled arrowheads: technical-residual estimate",
    x = "Patient-level median Spearman ρ", y = NULL, shape = NULL
  ) +
  theme_pub(7.4) +
  theme(legend.position = "bottom")

mf4b_detection <- sp_detection %>%
  filter(is.finite(positive_fraction), n_spots > 0) %>%
  mutate(
    section_label = str_extract(sample_id, "^[0-9]{2}_[0-9]{5}"),
    section = fct_reorder(section_label, positive_fraction)
  )
write_panel_source(
  mf4b_detection, "MAIN", "MF4b_section_detection.csv", "Fig4b",
  "Part7/07_spatial_target_detection_by_section.csv",
  "Every evaluable spatial section retained and ordered by target-positive fraction",
  "Section", "Cross-section target detection; not patient-level independence"
)
p4b_detection <- ggplot(
  mf4b_detection,
  aes(x = positive_fraction, y = section, size = n_spots, colour = positive_fraction)
) +
  geom_point(alpha = 0.92) +
  scale_x_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, max(mf4b_detection$positive_fraction) * 1.12)
  ) +
  scale_colour_gradient(low = "#E7A098", high = COL$red, guide = "none") +
  scale_size_continuous(range = c(1.8, 4.2), guide = "none") +
  labs(
    title = "Target detection across sections",
    subtitle = "Point position is positive fraction; area is analysed spots",
    x = "OLFML2B-positive fraction", y = NULL
  ) +
  theme_pub(6.4) +
  theme(
    axis.text.y = element_text(size = 5.4),
    plot.margin = margin(4, 4, 4, 2)
  ) +
  coord_cartesian(clip = "off")

top_row_4_core <- (p4a | p4b_detection) +
  plot_layout(widths = c(0.75, 0.25))
top_row_4 <- patchwork::free(
  top_row_4_core, side = "l", type = "panel"
)
middle_row_4 <- (p4b | p4c) + plot_layout(widths = c(0.40, 0.60))
bottom_row_4 <- (p4e_joint | p4d) + plot_layout(widths = c(0.44, 0.56))
fig4 <- top_row_4 / middle_row_4 / bottom_row_4 +
  plot_layout(heights = c(0.93, 0.84, 1.03)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) &
  main_panel_sanitizer
save_figure_bundle(
  fig4, "MAIN", "MAIN_FIGURE_04_CNS_v8_10_0",
  FINAL_MAIN_WIDTH_IN, 7.05
)

# =============================================================================
# MAIN FIGURE 5. Prognostic synthesis and orthogonal protein support
# =============================================================================

survival_models <- read_contract(
  "Part3", "05_olfml2b_survival_all_models.csv",
  c("cohort", "endpoint", "endpoint_family", "model", "model_id", "n", "events",
    "hr", "ci_low", "ci_high", "p_value", "ph_status", "status",
    "primary_inference_eligible", "meta_eligible")
)
survival_meta <- read_contract(
  "Part3", "08_olfml2b_meta_analysis.csv",
  c("endpoint_family", "model", "k", "cohorts", "hr_random", "ci_low_random",
    "ci_high_random", "p_random", "prediction_low", "prediction_high", "i2", "status")
)
survival_loo <- read_contract(
  "Part3", "09_olfml2b_leave_one_out_meta.csv",
  c("endpoint_family", "model", "hr_random", "ci_low_random", "ci_high_random",
    "p_random", "i2", "left_out", "status")
)
protein_pairs <- read_contract(
  "Part5", "09_OLFML2B_case_paired_deltas.csv",
  c("value_type", "case_id", "n_same_plex_pairs", "paired_delta_tumor_minus_normal", "pairing_scope")
)
protein_test <- read_contract(
  "Part5", "10_OLFML2B_case_paired_direction_tests.csv",
  c("value_type", "n_paired_cases", "median_delta", "bootstrap_median_ci_low",
    "bootstrap_median_ci_high", "positive_fraction", "wilcox_signed_rank_p", "exact_sign_p")
)
protein_loo <- read_contract(
  "Part5", "10a_OLFML2B_leave_one_plex_out_sensitivity.csv",
  c("omitted_plex", "n_remaining_cases", "median_delta", "positive_fraction",
    "wilcox_signed_rank_p", "exact_sign_p", "direction")
)
protein_plex <- read_contract(
  "Part5", "10b_OLFML2B_plex_level_direction_summary.csv",
  c("analytical_sample", "n_cases", "median_delta", "positive_fraction", "plex_direction")
)
protein_sentinel <- read_contract(
  "Part5", "12_sentinel_direction_summary.csv",
  c("gene", "group", "expected_direction", "n_pairs", "median_delta",
    "positive_fraction", "direction_concordant", "wilcox_fdr")
)

mf5a_individual <- survival_models %>%
  filter(
    status == "OK",
    (primary_inference_eligible %in% TRUE & model == "adjusted") |
      (cohort %in% c("GSE84437", "GSE26253") & model == "available_adjusted")
  ) %>%
  distinct(cohort, endpoint, .keep_all = TRUE) %>%
  mutate(
    endpoint_group = if_else(endpoint_family == "OS", "Overall survival", "DFS/RFS"),
    row_label = paste0(cohort_label(cohort), " ", endpoint, " (", events, "/", n, ")"),
    model_class = if_else(model == "adjusted", "Common adjusted", "Available-covariate"),
    row_type = "Cohort",
    prediction_low = NA_real_, prediction_high = NA_real_,
    meta_p = NA_real_, meta_i2 = NA_real_, meta_k = NA_integer_
  )
mf5a_meta <- survival_meta %>%
  filter(
    status == "OK",
    (endpoint_family == "OS" & model == "adjusted") |
      (endpoint_family == "RECURRENCE" & model == "available_adjusted")
  ) %>%
  transmute(
    endpoint_family, model, cohorts,
    endpoint_group = if_else(endpoint_family == "OS", "Overall survival", "DFS/RFS"),
    row_label = paste0("Pooled random effects (k=", k, ")"),
    row_type = "Pooled",
    model_class = "Pooled random effects",
    hr = hr_random, ci_low = ci_low_random, ci_high = ci_high_random,
    prediction_low, prediction_high,
    meta_p = p_random, meta_i2 = i2, meta_k = k,
    n = NA_integer_, events = NA_integer_, p_value = p_random,
    meta_eligible = TRUE
  ) %>%
  mutate(
    pi_separately_identified =
      is.finite(prediction_low) & is.finite(prediction_high) &
      !(dplyr::near(prediction_low, ci_low) & dplyr::near(prediction_high, ci_high)),
    pi_plot_low = if_else(pi_separately_identified, prediction_low, NA_real_),
    pi_plot_high = if_else(pi_separately_identified, prediction_high, NA_real_)
  )
mf5a_pool_members <- mf5a_meta %>%
  select(endpoint_group, cohorts) %>%
  tidyr::separate_rows(cohorts, sep = ";") %>%
  transmute(endpoint_group, cohort = cohorts, included_in_endpoint_pool = TRUE)
mf5a_individual <- mf5a_individual %>%
  left_join(mf5a_pool_members, by = c("endpoint_group", "cohort")) %>%
  mutate(
    included_in_endpoint_pool = included_in_endpoint_pool %in% TRUE,
    pool_display = factor(
      if_else(included_in_endpoint_pool, "Included in endpoint pool", "Shown, not pooled"),
      levels = c("Included in endpoint pool", "Shown, not pooled")
    )
  )
mf5a_order <- dplyr::bind_rows(
  mf5a_individual %>%
    select(endpoint_group, row_label) %>%
    mutate(group_order = if_else(endpoint_group == "Overall survival", 1L, 2L),
           pooled_order = 0L),
  mf5a_meta %>%
    select(endpoint_group, row_label) %>%
    mutate(group_order = if_else(endpoint_group == "Overall survival", 1L, 2L),
           pooled_order = 1L)
) %>%
  arrange(group_order, pooled_order, row_label)
mf5a_forest <- dplyr::bind_rows(mf5a_individual, mf5a_meta) %>%
  mutate(
    endpoint_group = factor(endpoint_group, levels = c("Overall survival", "DFS/RFS")),
    row_label = factor(row_label, levels = rev(unique(mf5a_order$row_label)))
  )
write_panel_source(
  mf5a_forest, "MAIN", "MF5a_adjusted_survival_and_meta_forest.csv", "Fig5a",
  "Part3/05_olfml2b_survival_all_models.csv; Part3/08_olfml2b_meta_analysis.csv",
  "Selected adjusted cohort estimates integrated with frozen pooled estimates; PI displayed only when separately identified from CI",
  "Patient within cohort; cohort for meta-analysis",
  "Observational prognosis with explicit small-k uncertainty; no clinical utility claim"
)
p5a <- ggplot(mf5a_forest, aes(x = hr, y = row_label)) +
  geom_vline(xintercept = 1, colour = "#B8B8B8", linewidth = 0.45) +
  geom_tile(
    data = mf5a_forest %>% filter(row_type == "Pooled", pi_separately_identified),
    aes(
      x = (pi_plot_low + pi_plot_high) / 2,
      width = pmax(pi_plot_high - pi_plot_low, .Machine$double.eps)
    ),
    height = 0.54, fill = "#E4E6E9", colour = NA
  ) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = row_label),
    colour = COL$grey, linewidth = 0.62
  ) +
  geom_point(
    data = mf5a_forest %>% filter(row_type == "Cohort"),
    aes(fill = pool_display, shape = pool_display), size = 2.45,
    colour = COL$grey, stroke = 0.58
  ) +
  geom_point(
    data = mf5a_forest %>% filter(row_type == "Pooled"),
    shape = 23, size = 3.5, fill = COL$red, colour = "white", stroke = 0.55
  ) +
  geom_text(
    data = mf5a_forest %>% filter(row_type == "Pooled"),
    aes(
      x = ci_high * 1.02,
      label = paste0("P=", vapply(meta_p, fmt_p, character(1)),
                     " · I² ", round(meta_i2), "% · k=", meta_k)
    ),
    hjust = 0, size = 1.82
  ) +
  facet_wrap(vars(endpoint_group), ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Included in endpoint pool" = COL$grey,
    "Shown, not pooled" = "white"
  )) +
  scale_shape_manual(values = c(
    "Included in endpoint pool" = 21,
    "Shown, not pooled" = 24
  )) +
  scale_x_log10(
    breaks = c(0.5, 1, 1.5, 2, 3),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.24))
  ) +
  labs(
    title = "Adjusted prognosis and pooled uncertainty",
    subtitle = if (any(mf5a_meta$pi_separately_identified)) {
      "Cohort and pooled 95% CIs; pale blocks are separately identified prediction intervals"
    } else {
      "Cohort and pooled 95% CIs; no separately identified prediction interval was exported for k≤3"
    },
    x = "Hazard ratio per 1-SD OLFML2B (log scale)", y = NULL, fill = NULL, shape = NULL
  ) +
  theme_pub(6.9) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#F3F4F6", colour = NA),
    strip.text = element_text(face = "bold", hjust = 0),
    panel.spacing.y = grid::unit(0.55, "lines"),
    plot.margin = margin(4, 8, 2, 8)
  ) +
  coord_cartesian(clip = "on")

mf5b <- protein_pairs %>%
  filter(value_type == "Log_Ratio", is.finite(paired_delta_tumor_minus_normal)) %>%
  arrange(paired_delta_tumor_minus_normal) %>%
  mutate(case_order = row_number(), direction = paired_delta_tumor_minus_normal >= 0)
protein_stat <- protein_test %>% filter(value_type == "Log_Ratio")
write_panel_source(
  mf5b, "MAIN", "MF5b_protein_paired_distribution.csv", "Fig5b",
  "Part5/09_OLFML2B_case_paired_deltas.csv; Part5/10_OLFML2B_case_paired_direction_tests.csv",
  "Same-plex Log_Ratio pairs only; every eligible case shown in the raw distribution",
  "Case", "Single-cohort orthogonal protein direction"
)
if (nrow(protein_stat) != 1L) {
  stop("Fig5b preflight failed: exactly one Log_Ratio protein-test row is required.")
}
protein_hl <- tryCatch(
  suppressWarnings(stats::wilcox.test(
    mf5b$paired_delta_tumor_minus_normal, mu = 0,
    conf.int = TRUE, conf.level = 0.95, exact = FALSE
  )),
  error = function(e) NULL
)
protein_hl_estimate <- if (!is.null(protein_hl) && length(protein_hl$estimate) > 0L) {
  unname(protein_hl$estimate[[1]])
} else NA_real_
protein_hl_ci <- if (!is.null(protein_hl) && length(protein_hl$conf.int) >= 2L) {
  unname(protein_hl$conf.int[1:2])
} else c(NA_real_, NA_real_)
protein_loo_max_p <- max(protein_loo$wilcox_signed_rank_p, na.rm = TRUE)
if (!is.finite(protein_loo_max_p)) {
  stop("Fig5b preflight failed: no finite leave-one-plex-out Wilcoxon P value.")
}
protein_plex_status <- if (protein_loo_max_p >= 0.05) "plex-sensitive" else "plex-robust"
protein_stat <- protein_stat %>%
  mutate(
    hodges_lehmann_pseudomedian = protein_hl_estimate,
    hodges_lehmann_ci_low = protein_hl_ci[1],
    hodges_lehmann_ci_high = protein_hl_ci[2],
    leave_one_plex_out_max_wilcox_p = protein_loo_max_p,
    plex_sensitive = protein_loo_max_p >= 0.05
  )
write_panel_source(
  protein_stat, "MAIN", "MF5b_protein_statistics.csv", "Fig5b-statistics",
  "Part5/10_OLFML2B_case_paired_direction_tests.csv; Part5/10a_OLFML2B_leave_one_plex_out_sensitivity.csv",
  "Frozen direction test retained; Hodges-Lehmann pseudomedian/CI and maximum leave-one-plex-out P added",
  "Case", "Single-cohort protein direction with plex-sensitivity qualification"
)
ps <- protein_stat[1, ]
protein_n_positive <- sum(mf5b$paired_delta_tumor_minus_normal > 0)
protein_n_negative <- sum(mf5b$paired_delta_tumor_minus_normal < 0)
p5b <- ggplot(mf5b, aes(x = factor("Paired cases"), y = paired_delta_tumor_minus_normal)) +
  geom_violin(fill = COL$teal, colour = NA, alpha = 0.28, width = 0.78, trim = FALSE) +
  geom_point(
    aes(colour = direction),
    position = position_jitter(width = 0.09, height = 0, seed = PART9_SEED),
    size = 1.35, alpha = 0.82
  ) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.2,
               fill = COL$gold, colour = "white", stroke = 0.55) +
  scale_colour_manual(values = c(`FALSE` = COL$blue, `TRUE` = COL$red), guide = "none") +
  labs(
    title = "Same-plex protein direction",
    subtitle = paste0(
      "n=", ps$n_paired_cases, "; ↑", protein_n_positive, " / ↓", protein_n_negative,
      "; Hodges-Lehmann=", formatC(ps$hodges_lehmann_pseudomedian, digits = 2, format = "f"),
      " [95% CI ", formatC(ps$hodges_lehmann_ci_low, digits = 2, format = "f"),
      ", ", formatC(ps$hodges_lehmann_ci_high, digits = 2, format = "f"),
      "]\nWilcoxon P = ", fmt_p(ps$wilcox_signed_rank_p),
      "; leave-one-plex-out max P = ", fmt_p(ps$leave_one_plex_out_max_wilcox_p),
      " (", protein_plex_status, ")"
    ),
    x = NULL, y = "Protein Δ (tumour–normal)"
  ) +
  scale_x_discrete(expand = expansion(add = 0.28)) +
  theme_pub(7.0) +
  theme(
    axis.title.y = element_text(margin = margin(r = 0)),
    axis.text.y = element_text(margin = margin(r = 1.0)),
    plot.margin = margin(1, 4, 2, 1)
  )

mf5d <- protein_sentinel %>%
  filter(is.finite(median_delta)) %>%
  mutate(
    gene = fct_reorder(gene, median_delta),
    expected = factor(
      expected_direction,
      levels = c("NEGATIVE", "NO_FIXED_EXPECTATION", "POSITIVE")
    )
  )
write_panel_source(
  mf5d, "MAIN", "MF5d_orientation_sentinels.csv", "Fig5d",
  "Part5/12_sentinel_direction_summary.csv",
  "All evaluable prespecified assay-orientation sentinels retained; promoted from Supplementary Fig. S5b",
  "Case-paired sentinel", "Direction and biological-plausibility audit; not target-gene validation"
)
sentinel_n_text <- if (n_distinct(mf5d$n_pairs) == 1L) {
  paste0("n=", first(mf5d$n_pairs), " paired cases")
} else {
  paste0("n=", min(mf5d$n_pairs), "–", max(mf5d$n_pairs), " paired cases")
}
sentinel_top_gene <- tail(levels(mf5d$gene), 1)
sentinel_x_range <- range(mf5d$median_delta, na.rm = TRUE)
sentinel_x_span <- diff(sentinel_x_range)
if (!all(is.finite(sentinel_x_range)) || !is.finite(sentinel_x_span) || sentinel_x_span <= 0) {
  stop("Fig5d preflight failed: sentinel effect range is invalid.")
}
sentinel_legend_top_rows <- tail(
  levels(mf5d$gene),
  max(6L, ceiling(nlevels(mf5d$gene) * 0.38))
)
sentinel_legend_x_cut <- sentinel_x_range[1] + 0.62 * sentinel_x_span
sentinel_legend_collision <- mf5d %>%
  filter(
    as.character(gene) %in% sentinel_legend_top_rows,
    median_delta <= sentinel_legend_x_cut
  )
if (nrow(sentinel_legend_collision) > 0L) {
  stop(
    "Fig5d inset-legend safety preflight failed: the upper-left legend region contains ",
    nrow(sentinel_legend_collision), " sentinel estimate(s)."
  )
}
p5d <- ggplot(mf5d, aes(x = median_delta, y = gene)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.42) +
  geom_point(
    aes(colour = group, shape = expected, size = n_pairs),
    alpha = 0.92, stroke = 0.72
  ) +
  annotate(
    "text", x = Inf, y = sentinel_top_gene, label = sentinel_n_text,
    hjust = 1.04, vjust = -0.75, size = 1.75, colour = COL$grey
  ) +
  scale_colour_manual(values = c(
    CAF_ECM = COL$red,
    Gastric_Normal = COL$gold,
    Immune_Control = COL$green,
    TGFb = COL$blue,
    Tumor_Epithelial = COL$purple
  ), labels = c(
    CAF_ECM = "CAF/ECM", Gastric_Normal = "Gastric normal",
    Immune_Control = "Immune control", TGFb = "TGF-β",
    Tumor_Epithelial = "Tumour epithelial"
  )) +
  scale_shape_manual(values = c(
    NEGATIVE = 25,
    NO_FIXED_EXPECTATION = 21,
    POSITIVE = 24
  ), labels = c(
    NEGATIVE = "Negative", NO_FIXED_EXPECTATION = "No fixed", POSITIVE = "Positive"
  )) +
  scale_size_continuous(range = c(1.9, 3.5), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  labs(
    x = "Median tumour-normal protein difference", y = NULL,
    colour = "Sentinel group", shape = "Expected"
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    shape = guide_legend(order = 2, nrow = 1)
  ) +
  coord_cartesian(clip = "off") +
  theme_pub(6.4) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.015, 0.985),
    legend.justification.inside = c(0, 1),
    legend.box = "vertical", legend.box.just = "left",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0, "mm"),
    legend.spacing.y = grid::unit(0.25, "mm"),
    legend.spacing.x = grid::unit(0.45, "mm"),
    legend.title = element_text(size = 5.1),
    legend.text = element_text(size = 5.0),
    legend.key.width = grid::unit(0.70, "lines"),
    legend.key.height = grid::unit(0.58, "lines"),
    plot.margin = margin(2, 6, 2, 1)
  )

mf5c_plex <- protein_plex %>%
  filter(is.finite(median_delta), n_cases > 0) %>%
  transmute(
    plex_code = str_extract(analytical_sample, "^[0-9]{2}CPTAC"),
    display_layer = "Plex-specific median",
    median_delta,
    nominal_boundary = FALSE,
    supporting_n = n_cases
  )
mf5c_loo <- protein_loo %>%
  filter(is.finite(median_delta), is.finite(wilcox_signed_rank_p)) %>%
  transmute(
    plex_code = str_extract(omitted_plex, "^[0-9]{2}CPTAC"),
    display_layer = "Pooled median after omission",
    median_delta,
    nominal_boundary = wilcox_signed_rank_p >= 0.05,
    supporting_n = n_remaining_cases
  )
if (
  any(is.na(mf5c_plex$plex_code)) || anyDuplicated(mf5c_plex$plex_code) ||
  any(is.na(mf5c_loo$plex_code)) || anyDuplicated(mf5c_loo$plex_code)
) {
  stop("Fig5c preflight failed: plex codes must be complete and unique within each layer.")
}
if (!setequal(mf5c_plex$plex_code, mf5c_loo$plex_code)) {
  stop("Fig5c preflight failed: plex-specific and post-omission plex sets do not match.")
}
plex_levels <- sort(unique(mf5c_plex$plex_code))
mf5c <- dplyr::bind_rows(mf5c_plex, mf5c_loo) %>%
  mutate(
    plex = factor(plex_code, levels = rev(plex_levels)),
    display_layer = factor(
      display_layer,
      levels = c("Plex-specific median", "Pooled median after omission")
    ),
    boundary_status = factor(
      if_else(
        display_layer == "Pooled median after omission" & nominal_boundary,
        "P≥0.05 after omission",
        "Estimate"
      ),
      levels = c("Estimate", "P≥0.05 after omission")
    )
  )
write_panel_source(
  mf5c, "MAIN", "MF5c_two_column_plex_robustness.csv", "Fig5c",
  "Part5/10b_OLFML2B_plex_level_direction_summary.csv; Part5/10a_OLFML2B_leave_one_plex_out_sensitivity.csv",
  "Plex-specific and leave-one-plex-out pooled medians retained in separate aligned effect-size columns; post-omission nominal boundary is encoded by point shape",
  "Analytical plex or case after plex omission", "Batch-direction stability with explicit plex sensitivity"
)
protein_stability_limit <- max(abs(mf5c$median_delta), na.rm = TRUE)
p5c <- ggplot(
  mf5c,
  aes(x = median_delta, y = plex, shape = boundary_status, fill = median_delta >= 0)
) +
  annotate(
    "rect",
    xmin = -protein_stability_limit * 0.045,
    xmax = protein_stability_limit * 0.045,
    ymin = -Inf, ymax = Inf,
    fill = "#F0F0F0", colour = NA
  ) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.42) +
  geom_point(
    size = 2.25, stroke = 0.68, colour = COL$black
  ) +
  scale_shape_manual(values = c(
    "Estimate" = 21,
    "P≥0.05 after omission" = 25
  ), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = COL$blue, `TRUE` = COL$red), guide = "none") +
  scale_x_continuous(
    limits = c(-protein_stability_limit * 1.12, protein_stability_limit * 1.12),
    breaks = scales::pretty_breaks(n = 3),
    expand = expansion(mult = 0.01)
  ) +
  facet_grid(
    cols = vars(display_layer),
    labeller = as_labeller(c(
      "Plex-specific median" = "Plex-specific median",
      "Pooled median after omission" = "Pooled after omission\n▼ Wilcoxon P≥0.05"
    ))
  ) +
  labs(x = "Protein Δ (tumour–normal)", y = NULL, shape = NULL) +
  theme_pub(6.1) +
  theme(
    axis.text.y = element_text(size = 5.35),
    axis.text.x = element_text(size = 5.35),
    strip.text = element_text(size = 5.9, face = "bold", lineheight = 0.94),
    panel.spacing.x = grid::unit(3.0, "mm"),
    legend.position = "none",
    plot.margin = margin(1, 4, 2, 2)
  )

sf1f <- survival_loo %>%
  filter(status == "OK", endpoint_family == "OS") %>%
  mutate(
    label = paste0("Omit ", cohort_label(left_out)),
    ci_crosses_one = ci_low_random <= 1 & ci_high_random >= 1
  )
full_os_meta <- mf5a_meta %>% filter(endpoint_group == "Overall survival")
if (nrow(full_os_meta) != 1L) {
  stop("FigS1f preflight failed: exactly one full OS meta row is required.")
}
sf1f <- sf1f %>%
  mutate(
    pooled_hr_change_percent = 100 * (hr_random / full_os_meta$hr - 1),
    label = fct_reorder(label, pooled_hr_change_percent),
    inference_status = factor(
      if_else(ci_crosses_one, "95% CI crosses 1", "95% CI excludes 1"),
      levels = c("95% CI crosses 1", "95% CI excludes 1")
    ),
    detail_label = paste0(
      "HR ", formatC(hr_random, digits = 2, format = "f"),
      " [", formatC(ci_low_random, digits = 2, format = "f"),
      ", ", formatC(ci_high_random, digits = 2, format = "f"), "]"
    ),
    detail_x = if_else(pooled_hr_change_percent >= 0,
                       pooled_hr_change_percent - 0.10, 0.12),
    detail_hjust = if_else(pooled_hr_change_percent >= 0, 1, 0),
    detail_colour = if_else(pooled_hr_change_percent >= 0, "white", COL$black)
  )
loo_status_subtitle <- if (all(sf1f$ci_crosses_one)) {
  paste0("All ", nrow(sf1f), " leave-one-out 95% CIs cross 1\nFull pooled HR=",
         formatC(full_os_meta$hr, digits = 2, format = "f"),
         "; bars show relative HR change")
} else {
  paste0(sum(sf1f$ci_crosses_one), "/", nrow(sf1f),
         " leave-one-out 95% CIs cross 1\nFull pooled HR=",
         formatC(full_os_meta$hr, digits = 2, format = "f"),
         "; bars show relative HR change")
}
write_panel_source(
  sf1f, "SUPPLEMENTARY", "SF1d_leave_one_cohort_out_OS.csv", "FigS1d",
  "Part3/09_olfml2b_leave_one_out_meta.csv; Part3/08_olfml2b_meta_analysis.csv",
  "All frozen leave-one-cohort-out OS estimates expressed as percentage change from the full pooled HR; CI-crossing status retained",
  "Cohort omission", "Small-k influence analysis; not independent validation"
)
p_s1f <- ggplot(sf1f, aes(x = pooled_hr_change_percent, y = label)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_col(aes(fill = inference_status), width = 0.62, alpha = 0.92) +
  geom_text(
    aes(
      x = detail_x, label = detail_label,
      hjust = detail_hjust, colour = detail_colour
    ),
    size = 2.05
  ) +
  scale_fill_manual(
    values = c("95% CI crosses 1" = COL$blue, "95% CI excludes 1" = COL$red),
    guide = "none"
  ) +
  scale_colour_identity() +
  labs(
    title = "OS leave-one-cohort-out influence",
    subtitle = loo_status_subtitle,
    x = "Change in pooled HR after cohort omission (%)", y = NULL, fill = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  theme_pub(7.0) +
  coord_cartesian(clip = "off")

middle_row_5_core <- (p5b | p5c) +
  plot_layout(widths = c(0.31, 0.69))
middle_row_5 <- patchwork::free(
  middle_row_5_core, side = "l", type = "panel"
)
bottom_row_5 <- patchwork::free(p5d, side = "l", type = "panel")
fig5 <- p5a / middle_row_5 / bottom_row_5 +
  plot_layout(heights = c(0.88, 0.80, 0.86)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) &
  main_panel_sanitizer
save_figure_bundle(
  fig5, "MAIN", "MAIN_FIGURE_05_CNS_v8_10_0",
  FINAL_MAIN_WIDTH_IN, 6.80
)

# =============================================================================
# SUPPLEMENTARY FIGURE S1. Clinical and survival diagnostics
# Non-stage clinical screens and survival-model diagnostics remain supplementary qualifiers.
# =============================================================================

clinical_context <- read_contract(
  "Part3", "03c_olfml2b_expression_by_clinical_context.csv",
  c("cohort", "variable", "level", "n", "median", "mean", "global_p")
)
spline_audit <- read_contract(
  "Part3", "07_olfml2b_nonlinearity_spline.csv",
  c("cohort", "endpoint", "n", "events", "nonlinearity_lrt_p", "status")
)
ph_time <- read_contract(
  "Part3", "07a_olfml2b_PH_time_interaction_sensitivity.csv",
  c("cohort", "endpoint", "n", "events", "p_log_time_interaction",
    "hr_day365", "hr_day1095", "status")
)
sf1a <- clinical_context %>%
  mutate(variable_lower = str_to_lower(variable)) %>%
  filter(variable_lower %in% c("sex", "lauren", "lauren type", "molecular_subtype", "molecular subtype")) %>%
  group_by(cohort, variable_lower) %>%
  summarise(
    global_p = {
      observed_p <- global_p[!is.na(global_p)]
      if (length(observed_p) > 0L) observed_p[[1]] else NA_real_
    },
    n = sum(n, na.rm = TRUE), .groups = "drop"
  ) %>%
  mutate(
    variable_label = recode(
      variable_lower, sex = "Sex", lauren = "Lauren type", `lauren type` = "Lauren type",
      molecular_subtype = "Molecular subtype", `molecular subtype` = "Molecular subtype"
    ),
    cohort = cohort_label(cohort)
  ) %>%
  select(cohort, variable_label, global_p, n) %>%
  tidyr::complete(
    cohort = c("GSE15459", "GSE62254", "GSE84437", "TCGA-STAD"),
    variable_label = c("Sex", "Lauren type", "Molecular subtype")
  ) %>%
  mutate(
    strength_raw = case_when(
      is.na(global_p) ~ NA_real_,
      global_p <= 0 ~ 300,
      TRUE ~ -log10(global_p)
    ),
    strength_display = pmin(strength_raw, 20),
    cell_label = case_when(
      is.na(global_p) ~ "×",
      global_p <= 0 ~ "P<1×10⁻³⁰⁰",
      global_p < 0.001 ~ "***",
      global_p < 0.01 ~ "**",
      global_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
write_panel_source(
  sf1a, "SUPPLEMENTARY", "SF1a_nonstage_clinical_screen.csv", "FigS1a",
  "Part3/03c_olfml2b_expression_by_clinical_context.csv",
  "Stage excluded; full cohort-variable grid completed; exact zero P is displayed as a bound and colour is capped",
  "Patient/sample", "Exploratory global-test screen; no direction or effect size is implied"
)
p_s1a <- ggplot(sf1a, aes(x = cohort, y = variable_label)) +
  geom_tile(aes(fill = strength_display), colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 2.0, fontface = "bold") +
  scale_fill_gradient(low = "#EFEFF1", high = COL$purple, na.value = "#F7F7F7") +
  labs(
    title = "Non-stage global-test screen",
    subtitle = "Stage is excluded; × = unavailable; colour caps -log10(P) at 20",
    x = NULL, y = NULL, fill = "Capped -log10(P)"
  ) + theme_pub(7.2) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

sf1b <- spline_audit %>%
  filter(status == "OK", is.finite(nonlinearity_lrt_p)) %>%
  mutate(
    label = paste0(cohort_label(cohort), " ", endpoint),
    neglogp = -log10(pmax(nonlinearity_lrt_p, .Machine$double.xmin)),
    label = fct_reorder(label, neglogp)
  )
write_panel_source(
  sf1b, "SUPPLEMENTARY", "SF1b_nonlinearity_audit.csv", "FigS1b",
  "Part3/07_olfml2b_nonlinearity_spline.csv",
  "All successful prespecified spline tests retained",
  "Patient within cohort/endpoint", "Shape diagnostic; not a threshold search"
)
p_s1b <- ggplot(sf1b, aes(x = neglogp, y = label)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = COL$grey, linewidth = 0.45) +
  geom_point(colour = COL$blue, size = 2.2) +
  labs(
    title = "Nonlinearity-test summary",
    subtitle = "Dashed line: nominal P = 0.05",
    x = "-log10(P for nonlinearity)", y = NULL
  ) + theme_pub(7.3)

sf1c <- ph_time %>%
  filter(status == "OK", is.finite(hr_day365), is.finite(hr_day1095)) %>%
  mutate(
    label = paste0(cohort_label(cohort), " ", endpoint),
    label = fct_reorder(label, hr_day365),
    label_x = pmax(hr_day365, hr_day1095) + 0.025,
    detail_label = paste0(
      "n/events=", n, "/", events,
      "; P-int=", vapply(p_log_time_interaction, fmt_p, character(1))
    )
  )
write_panel_source(
  sf1c, "SUPPLEMENTARY", "SF1c_time_varying_HR.csv", "FigS1c",
  "Part3/07a_olfml2b_PH_time_interaction_sensitivity.csv",
  "Frozen day-365 and day-1095 estimates for successful models",
  "Patient within cohort/endpoint", "PH sensitivity; does not replace primary Cox estimate"
)
p_s1c <- ggplot(sf1c, aes(y = label)) +
  geom_vline(xintercept = 1, colour = "#BDBDBD", linewidth = 0.45) +
  geom_segment(aes(x = hr_day365, xend = hr_day1095, yend = label),
               colour = "#C6C8CC", linewidth = 0.55) +
  geom_point(aes(x = hr_day365), colour = COL$blue, shape = 21,
             fill = "white", size = 2.4, position = position_nudge(y = -0.10)) +
  geom_point(aes(x = hr_day1095), colour = COL$red, shape = 17,
             size = 2.5, position = position_nudge(y = 0.10)) +
  geom_text(
    aes(x = label_x, label = detail_label),
    hjust = 0, nudge_y = -0.13, size = 1.85
  ) +
  labs(
    title = "Time-varying HR sensitivity",
    subtitle = "Open circles: day 365; triangles: day 1095; time-specific CIs are unavailable upstream",
    x = "Estimated time-specific HR", y = NULL
  ) +
  scale_x_continuous(
    limits = c(0.97, max(sf1c$label_x, na.rm = TRUE) + 0.40),
    breaks = c(1.0, 1.2, 1.4, 1.6)
  ) +
  theme_pub(7.2) +
  coord_cartesian(clip = "on")

fig_s1 <- (p_s1a | p_s1b) / (p_s1c | p_s1f) +
  plot_layout(heights = c(1.0, 1.0)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s1, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S1_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 5.70
)

# =============================================================================
# SUPPLEMENTARY FIGURE S2. Bulk TME sensitivity and collinearity
# =============================================================================

tme_collin <- read_contract(
  "Part4", "12d_OLFML2B_TME_collinearity_audit.csv",
  c("cohort", "endpoint", "attenuation_axis", "olfml2b_tme_spearman_rho",
    "vif_OLFML2B", "condition_index", "coefficient_sign_reversal")
)

secondary_shift_programs <- c("Treg", "M2_Macrophage", "IFNg_Response")
sf2a_secondary <- tme_shift %>%
  filter(
    analysis_population == "TUMOR_ONLY", signature %in% secondary_shift_programs,
    cohort %in% main_cohorts
  ) %>%
  mutate(
    cohort = factor(cohort_label(cohort), levels = cohort_label(main_cohorts)),
    program = factor(feature_label(signature), levels = rev(feature_label(secondary_shift_programs))),
    mark = if_else(!is.na(fdr) & fdr < 0.05, "*", "")
  )
write_panel_source(
  sf2a_secondary, "SUPPLEMENTARY", "SF2a_secondary_state_shift.csv", "FigS2a",
  "Part4/11_olfml2b_high_low_tme_shift.csv",
  "Only secondary programmes not used in Main Fig. 2a or 2c",
  "Patient/sample within cohort", "Display split only"
)
sf2b_plot <- sf2a_secondary %>%
  tidyr::complete(cohort, program) %>%
  mutate(
    display_label = if_else(
      is.finite(high_minus_low),
      paste0(formatC(high_minus_low, digits = 2, format = "f"), coalesce(mark, "")),
      "×"
    )
  )
p_s2a <- ggplot(sf2b_plot, aes(x = cohort, y = program, fill = high_minus_low)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = display_label), size = 2.2) +
  scale_fill_gradient2(
    low = COL$blue, mid = "white", high = COL$red, midpoint = 0,
    na.value = "#F2F2F2"
  ) +
  labs(
    title = "Secondary programme state shift",
    subtitle = "OLFML2B-high minus low median; × = unavailable; rows are excluded from main figures",
    x = NULL, y = NULL, fill = "Difference"
  ) + theme_pub(7.2) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

sf2b_collin <- tme_collin %>%
  filter(is.finite(olfml2b_tme_spearman_rho), is.finite(vif_OLFML2B), is.finite(condition_index)) %>%
  mutate(
    axis_label = feature_label(attenuation_axis),
    axis_class = if_else(attenuation_axis %in% primary_axes, "Primary axis", "Composite/exploratory"),
    reversal = coefficient_sign_reversal %in% TRUE,
    label = if_else(reversal, paste0(cohort_label(cohort), " ", endpoint), "")
  )
sf2d_labels <- sf2b_collin %>%
  filter(reversal) %>%
  group_by(cohort, endpoint) %>%
  slice_max(order_by = vif_OLFML2B, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste0(cohort_label(cohort), " ", endpoint))
write_panel_source(
  sf2b_collin, "SUPPLEMENTARY", "SF2b_collinearity_audit.csv", "FigS2b",
  "Part4/12d_OLFML2B_TME_collinearity_audit.csv",
  "All frozen collinearity rows retained",
  "Patient within cohort/endpoint", "Model diagnostic"
)
sf2d_long <- sf2b_collin %>%
  select(cohort, endpoint, attenuation_axis, axis_label, axis_class, reversal,
         olfml2b_tme_spearman_rho, vif_OLFML2B, condition_index) %>%
  pivot_longer(
    c(vif_OLFML2B, condition_index),
    names_to = "diagnostic", values_to = "diagnostic_value"
  ) %>%
  mutate(
    diagnostic = recode(
      diagnostic, vif_OLFML2B = "OLFML2B VIF", condition_index = "Condition index"
    ),
    diagnostic = factor(diagnostic, levels = c("OLFML2B VIF", "Condition index"))
  )
sf2d_thresholds <- tibble(
  diagnostic = factor(c("OLFML2B VIF", "Condition index"),
                      levels = c("OLFML2B VIF", "Condition index")),
  threshold = c(5, 10)
)
sf2d_labels_long <- sf2d_labels %>%
  transmute(
    cohort, endpoint, attenuation_axis, axis_label, axis_class, reversal,
    olfml2b_tme_spearman_rho, diagnostic_value = vif_OLFML2B,
    diagnostic = factor("OLFML2B VIF", levels = c("OLFML2B VIF", "Condition index")),
    label
  )
p_s2b <- ggplot(
  sf2d_long,
  aes(x = olfml2b_tme_spearman_rho, y = diagnostic_value)
) +
  geom_hline(
    data = sf2d_thresholds, aes(yintercept = threshold),
    linetype = 2, colour = "#D0D0D0", linewidth = 0.42
  ) +
  geom_point(aes(colour = axis_class, shape = reversal), size = 1.95, alpha = 0.80) +
  ggrepel::geom_text_repel(
    data = sf2d_labels_long,
    aes(label = label), size = 2.1, seed = PART9_SEED,
    box.padding = 0.55, point.padding = 0.35,
    force = 2.5, max.overlaps = Inf, show.legend = FALSE
  ) +
  facet_wrap(~diagnostic, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("Primary axis" = COL$purple,
                                 "Composite/exploratory" = COL$teal)) +
  labs(
    title = "Collinearity and sign reversal",
    subtitle = stringr::str_wrap(
      "Dashed thresholds: VIF=5; condition index=10; labels mark reversing endpoints",
      width = 48
    ),
    x = "Spearman ρ: OLFML2B versus TME axis", y = NULL,
    colour = NULL, shape = "Crossed HR=1"
  ) +
  guides(
    colour = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  ) +
  theme_pub(6.5) +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    legend.box.just = "left", legend.spacing.y = grid::unit(0.2, "mm"),
    legend.text = element_text(size = 5.5)
  )

fig_s2 <- (p_s2a | p_s2b) / p_s2c_attenuation +
  plot_layout(heights = c(1, 0.92)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s2, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S2_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 5.75
)

# =============================================================================
# SUPPLEMENTARY FIGURE S3. Single-cell extensions and sensitivity
# =============================================================================

sc_state <- read_contract(
  "Part6", "69_within_compartment_program_tests.csv",
  c("cohort", "tissue_or_condition", "marker_celltype", "annotation_confidence_scope",
    "marker_program", "orthogonality_status", "n_units", "median_delta",
    "formal_evaluable", "wilcox_fdr_within_family")
)

state_programs <- c(
  "Endothelial_Angiogenic", "Epithelial_Differentiation", "Inflammatory_Fibroblast",
  "Myeloid_Macrophage", "Pericyte", "Proliferation", "TGFb_Response"
)
sf3b_state <- sc_state %>%
  filter(
    annotation_confidence_scope == "ALL_ASSIGNED",
    orthogonality_status == "ORTHOGONAL_TO_COMPARTMENT_CLASSIFIER",
    formal_evaluable %in% TRUE, marker_program %in% state_programs
  ) %>%
  mutate(
    column = paste0(cohort, "\n", feature_label(marker_celltype)),
    program = factor(feature_label(marker_program), levels = rev(feature_label(state_programs))),
    mark = if_else(!is.na(wilcox_fdr_within_family) & wilcox_fdr_within_family < 0.05, "*", "")
  )
write_panel_source(
  sf3b_state, "SUPPLEMENTARY", "SF3b_within_compartment_orthogonal_state.csv", "FigS3b",
  "Part6/69_within_compartment_program_tests.csv",
  "ALL_ASSIGNED, orthogonal, formally evaluable programmes only",
  "Patient/official sample", "State association within compartment; not mechanism"
)
p_s3b <- ggplot(sf3b_state, aes(x = median_delta, y = program, fill = median_delta >= 0)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.42) +
  geom_col(width = 0.66, alpha = 0.88) +
  geom_point(
    data = sf3b_state %>% filter(mark == "*"),
    aes(x = median_delta, y = program), inherit.aes = FALSE,
    shape = 8, size = 1.55, colour = COL$black
  ) +
  facet_wrap(~column, ncol = 3) +
  scale_fill_manual(
    values = c(`FALSE` = COL$blue, `TRUE` = COL$red),
    labels = c(`FALSE` = "Negative difference", `TRUE` = "Positive difference")
  ) +
  labs(
    title = "Within-compartment orthogonal state",
    subtitle = "Identity-coupled programmes excluded\nBlack star denotes within-family FDR<0.05",
    x = "Median positive-minus-negative programme difference", y = NULL, fill = NULL
  ) + theme_pub(6.4) + theme(legend.position = "bottom")

fig_s3 <- p_s3b +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s3, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S3_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 4.10
)

# =============================================================================
# SUPPLEMENTARY FIGURE S4. Spatial sensitivity
# Global autocorrelation, distance and technical residualization are retained as
# diagnostic qualifiers; threshold overlap is promoted to Main Fig. 4c.
# =============================================================================

sp_moran <- read_contract(
  "Part7", "51_global_moransI_patient_summary.csv",
  c("effect_source", "feature", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "wilcoxon_fdr")
)
sp_distance <- read_contract(
  "Part7", "71_spatial_distance_patient_summary.csv",
  c("effect_source", "feature", "n_patients", "median_effect",
    "bootstrap_ci_low", "bootstrap_ci_high", "exact_signflip_fdr", "wilcoxon_fdr")
)
sf4a <- sp_moran %>%
  filter(is.finite(median_effect)) %>%
  mutate(
    label_text = paste0(
      feature_label(feature), " | q=", vapply(wilcoxon_fdr, fmt_p, character(1))
    ),
    label = fct_reorder(label_text, median_effect)
  )
write_panel_source(
  sf4a, "SUPPLEMENTARY", "SF4a_global_MoransI.csv", "FigS4a",
  "Part7/51_global_moransI_patient_summary.csv",
  "All frozen patient-level Moran summaries retained",
  "Patient after within-patient section collapse", "Global spatial autocorrelation"
)
p_s4a <- ggplot(sf4a, aes(x = median_effect, y = label)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_tile(
    aes(
      x = (bootstrap_ci_low + bootstrap_ci_high) / 2,
      width = pmax(bootstrap_ci_high - bootstrap_ci_low, .Machine$double.eps)
    ),
    height = 0.52, fill = "#E0D9EB", colour = NA
  ) +
  geom_point(colour = COL$purple, size = 2.1) +
  labs(
    title = "Global spatial autocorrelation",
    subtitle = "Patient-level Moran's I summaries; retained as a spatial diagnostic",
    x = "Median Moran's I (bootstrap 95% CI)", y = NULL
  ) + theme_pub(7.0)

sf4b <- sp_distance %>%
  filter(is.finite(median_effect)) %>%
  mutate(
    proximity = -1000 * median_effect,
    proximity_low = -1000 * bootstrap_ci_high,
    proximity_high = -1000 * bootstrap_ci_low,
    significant = !is.na(exact_signflip_fdr) & exact_signflip_fdr < 0.05,
    label = paste0(feature_label(feature), " | exact q=", vapply(exact_signflip_fdr, fmt_p, character(1)))
  )
write_panel_source(
  sf4b, "SUPPLEMENTARY", "SF4b_distance_permutation.csv", "FigS4b",
  "Part7/71_spatial_distance_patient_summary.csv",
  "Sign reversed and multiplied by 1000 so positive means closer than matched null",
  "Patient after within-patient section collapse", "Distance sensitivity; not physical cell-cell contact"
)
p_s4b <- ggplot(sf4b, aes(x = proximity, y = fct_reorder(label, proximity))) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_tile(
    aes(
      x = (proximity_low + proximity_high) / 2,
      width = pmax(proximity_high - proximity_low, .Machine$double.eps)
    ),
    height = 0.52, fill = "#DCE7F4", colour = NA
  ) +
  geom_point(aes(colour = significant), size = 2.1) +
  scale_colour_manual(values = c(`FALSE` = COL$grey, `TRUE` = COL$blue), guide = "none") +
  labs(
    title = "Exploratory distance-\npermutation estimates",
    subtitle = stringr::str_wrap(
      paste0(
        "Positive means closer than matched null; ", sum(sf4b$significant), "/",
        nrow(sf4b), " exact sign-flip FDR<0.05"
      ),
      width = 48
    ),
    x = "Proximity index (10^-3 normalized spot units)", y = NULL
  ) +
  theme_pub(7.1) +
  theme(plot.title = element_text(lineheight = 0.92))

fig_s4 <- (p_s4a | p_s4b) / p_s4d_technical +
  plot_layout(heights = c(1, 0.82)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s4, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S4_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 5.95
)

# =============================================================================
# SUPPLEMENTARY FIGURE S5. Protein pairing audit
# Orientation sentinels and the combined plex-stability point display are promoted to Main Fig. 5.
# =============================================================================

protein_attrition <- read_contract(
  "Part5", "09a2_OLFML2B_pair_attrition_summary.csv",
  c("step", "n_cases", "interpretation")
)
sf5a <- protein_attrition %>%
  mutate(
    step_label = recode(
      step,
      cases_with_any_mapped_primary_sample = "Any mapped primary sample",
      cases_with_mapped_tumor = "Mapped tumour",
      cases_with_mapped_normal = "Mapped normal",
      mapped_tumor_normal_cases = "Mapped tumour + normal",
      cases_with_finite_target_tumor = "Finite OLFML2B tumour",
      cases_with_finite_target_normal = "Finite OLFML2B normal",
      cases_with_finite_target_both_tissues = "Finite OLFML2B in both tissues",
      included_same_plex_target_pairs = "Included same-plex pairs",
      .default = str_replace_all(step, "_", " ")
    ),
    step_label = factor(step_label, levels = rev(unique(step_label)))
  )
write_panel_source(
  sf5a, "SUPPLEMENTARY", "SF5a_pair_attrition.csv", "FigS5a",
  "Part5/09a2_OLFML2B_pair_attrition_summary.csv",
  "Frozen attrition counts retained; plotting coordinates encode tumour and normal branches before intersection",
  "Case", "Pairing audit"
)
get_attrition_n <- function(step_id) {
  value <- sf5a %>% filter(step == step_id) %>% pull(n_cases)
  if (length(value) != 1L || !is.finite(value)) {
    stop("FigS5a preflight failed for attrition step: ", step_id)
  }
  value
}
sf5a_nodes <- tribble(
  ~node, ~label, ~x, ~y,
  "root", "Any mapped\nprimary", 0, 0,
  "mapped_t", "Mapped tumour", 1.4, 0.85,
  "finite_t", "Finite tumour\nOLFML2B", 3.0, 0.85,
  "mapped_n", "Mapped normal", 1.4, -0.85,
  "finite_n", "Finite normal\nOLFML2B", 3.0, -0.85,
  "mapped_both", "Mapped tumour\n+ normal", 4.6, 0.58,
  "finite_both", "Finite OLFML2B\nin both tissues", 4.6, -0.58,
  "same_plex", "Included same-plex\npairs", 6.2, -0.58
) %>%
  mutate(n_cases = c(
    get_attrition_n("cases_with_any_mapped_primary_sample"),
    get_attrition_n("cases_with_mapped_tumor"),
    get_attrition_n("cases_with_finite_target_tumor"),
    get_attrition_n("cases_with_mapped_normal"),
    get_attrition_n("cases_with_finite_target_normal"),
    get_attrition_n("mapped_tumor_normal_cases"),
    get_attrition_n("cases_with_finite_target_both_tissues"),
    get_attrition_n("included_same_plex_target_pairs")
  ), display = paste0(label, "\nn=", n_cases))
sf5a_edges <- tribble(
  ~from, ~to,
  "root", "mapped_t", "mapped_t", "finite_t",
  "root", "mapped_n", "mapped_n", "finite_n",
  "mapped_t", "mapped_both", "mapped_n", "mapped_both",
  "finite_t", "finite_both", "finite_n", "finite_both",
  "mapped_both", "finite_both", "finite_both", "same_plex"
) %>%
  left_join(sf5a_nodes %>% select(from = node, x, y), by = "from") %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sf5a_nodes %>% select(to = node, x, y), by = "to") %>%
  rename(x_to = x, y_to = y)
p_s5a <- ggplot() +
  geom_segment(
    data = sf5a_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    colour = "#B9BBC0", linewidth = 0.65
  ) +
  geom_label(
    data = sf5a_nodes, aes(x = x, y = y, label = display),
    fill = "white", colour = COL$black, linewidth = 0.25,
    label.padding = grid::unit(0.12, "lines"), size = 1.85
  ) +
  annotate("text", x = 2.2, y = 1.42, label = "Tumour branch",
           colour = COL$purple, fontface = "bold", size = 2.05) +
  annotate("text", x = 2.2, y = -1.42, label = "Normal branch",
           colour = COL$purple, fontface = "bold", size = 2.05) +
  labs(
    title = "Branched case-pair eligibility",
    subtitle = "Tumour and normal availability are parallel branches, not a single decreasing funnel"
  ) +
  coord_cartesian(xlim = c(-0.45, 6.65), ylim = c(-1.70, 1.70), clip = "off") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 8.3, hjust = 0),
    plot.subtitle = element_text(colour = COL$grey, size = 6.2, hjust = 0),
    plot.margin = margin(8, 10, 8, 10)
  )

fig_s5 <- p_s5a +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s5, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S5_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 3.85
)

# =============================================================================
# SUPPLEMENTARY FIGURE S6. Exploratory immunotherapy boundary
# Entire ICI module is supplementary because it is single-cohort and exploratory.
# =============================================================================

ici_input <- read_contract(
  "Part8_Immunotherapy", "06_TIGER_tumor_only_response_analysis_input.csv",
  c("sample_id", "response_ORR_group", "target_expr")
)
ici_wilcox <- read_contract(
  "Part8_Immunotherapy", "08_ORR_response_wilcoxon_all_features.csv",
  c("feature", "n", "n_positive", "n_negative", "p_value", "p_fdr",
    "delta_median_positive_minus_negative", "status")
)
ici_auc <- read_contract(
  "Part8_Immunotherapy", "09_ORR_response_logistic_per_1SD_all_features.csv",
  c("feature", "n", "auc_for_positive_high_value", "status")
)
ici_corr <- read_contract(
  "Part8_Immunotherapy", "23_OLFML2B_module_spearman_correlations.csv",
  c("module", "n", "rho", "p_value", "p_fdr", "status")
)

ici_bootstrap_B <- suppressWarnings(as.integer(
  Sys.getenv("OLFML2B_PART9_ICI_BOOTSTRAP_B", unset = "2000")
))
if (!is.finite(ici_bootstrap_B) || ici_bootstrap_B < 500L) ici_bootstrap_B <- 2000L
ici_effects <- purrr::map_dfr(seq_along(unique(ici_wilcox$feature)), function(i) {
  current_feature <- unique(ici_wilcox$feature)[i]
  if (!current_feature %in% names(ici_input)) {
    return(tibble(
      feature = current_feature, auc = NA_real_, auc_low = NA_real_, auc_high = NA_real_,
      rank_biserial = NA_real_, rb_low = NA_real_, rb_high = NA_real_
    ))
  }
  positive_values <- ici_input[[current_feature]][
    ici_input$response_ORR_group == "R" & is.finite(ici_input[[current_feature]])
  ]
  negative_values <- ici_input[[current_feature]][
    ici_input$response_ORR_group == "NR" & is.finite(ici_input[[current_feature]])
  ]
  bootstrap_two_group_effect(
    positive_values, negative_values,
    B = ici_bootstrap_B, seed = PART9_SEED + i
  ) %>% mutate(feature = current_feature, .before = 1)
})

sf6a <- ici_input %>%
  filter(response_ORR_group %in% c("R", "NR"), is.finite(target_expr)) %>%
  mutate(response = factor(response_ORR_group, levels = c("NR", "R")))
sf6a_counts <- sf6a %>% count(response, name = "n")
sf6a <- sf6a %>%
  left_join(sf6a_counts, by = "response") %>%
  mutate(
    response_axis = factor(
      paste0(as.character(response), "\nn=", n),
      levels = paste0(c("NR", "R"), "\nn=", sf6a_counts$n[match(c("NR", "R"), sf6a_counts$response)])
    )
  )
sf6a_target_effect <- ici_effects %>% filter(feature == "target_expr")
sf6a_target_test <- ici_wilcox %>% filter(feature == "target_expr", status == "OK")
if (nrow(sf6a_target_effect) != 1L || nrow(sf6a_target_test) != 1L) {
  stop("FigS6a preflight failed: one target effect and one target test are required.")
}
write_panel_source(
  sf6a, "SUPPLEMENTARY", "SF6a_OLFML2B_by_ORR.csv", "FigS6a",
  "Part8_Immunotherapy/06_TIGER_tumor_only_response_analysis_input.csv",
  "Tumour-only ORR-evaluable patients; no model refitting",
  "Patient", "Single-cohort exploratory transportability"
)
p_s6a <- ggplot(sf6a, aes(x = response_axis, y = target_expr, fill = response)) +
  geom_boxplot(width = 0.34, outlier.shape = NA, alpha = 0.22, linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = PART9_SEED),
             size = 1.25, alpha = 0.75) +
  scale_fill_manual(values = c(NR = "#BDBDBD", R = COL$red), guide = "none") +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(
    title = "OLFML2B by objective response",
    subtitle = paste0(
      "PRJEB25780/TIGER; tumour-only; rank-biserial = ",
      formatC(sf6a_target_effect$rank_biserial, digits = 2, format = "f"),
      " [", formatC(sf6a_target_effect$rb_low, digits = 2, format = "f"), ", ",
      formatC(sf6a_target_effect$rb_high, digits = 2, format = "f"),
      "]; FDR = ", fmt_p(sf6a_target_test$p_fdr)
    ),
    x = NULL, y = "OLFML2B expression"
  ) + theme_pub(7.4)

sf6b <- ici_auc %>%
  left_join(ici_effects, by = "feature") %>%
  filter(is.finite(auc), is.finite(auc_low), is.finite(auc_high)) %>%
  mutate(
    label = fct_reorder(feature_label(feature), auc),
    target = feature == "target_expr"
  )
write_panel_source(
  sf6b, "SUPPLEMENTARY", "SF6b_directional_AUC.csv", "FigS6b",
  "Part8_Immunotherapy/09_ORR_response_logistic_per_1SD_all_features.csv",
  "Rank AUC recalculated from patient-level values with stratified bootstrap CI; unstable Firth coefficients remain suppressed",
  "Patient", "Descriptive discrimination; not validated prediction"
)
p_s6b <- ggplot(sf6b, aes(x = auc, y = label, colour = target)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#BDBDBD", linewidth = 0.45) +
  geom_segment(aes(x = auc_low, xend = auc_high, yend = label), linewidth = 0.55) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(`FALSE` = COL$grey, `TRUE` = COL$blue), guide = "none") +
  labs(
    title = "Directional response discrimination",
    subtitle = stringr::str_wrap(
      paste0("All Firth fits were unstable; stratified bootstrap B=", ici_bootstrap_B,
             "; AUC remains descriptive"),
      width = 48
    ),
    x = "Directional AUC", y = NULL
  ) + theme_pub(6.8)

sf6c <- ici_corr %>%
  filter(status == "OK", is.finite(rho)) %>%
  mutate(
    label = fct_reorder(feature_label(module), rho),
    mark = if_else(!is.na(p_fdr) & p_fdr < 0.05, "*", "")
  )
write_panel_source(
  sf6c, "SUPPLEMENTARY", "SF6c_OLFML2B_module_correlations.csv", "FigS6c",
  "Part8_Immunotherapy/23_OLFML2B_module_spearman_correlations.csv",
  "All successful frozen module correlations retained",
  "Patient", "Ecological context in one ICI cohort"
)
p_s6c <- ggplot(sf6c, aes(x = rho, y = label)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_point(aes(colour = rho), size = 2.15, alpha = 0.88) +
  geom_text(aes(label = mark), nudge_x = 0.035, size = 2.0, fontface = "bold") +
  scale_colour_gradient2(low = COL$blue, mid = "white", high = COL$red, midpoint = 0) +
  labs(
    title = "OLFML2B-module ecological context",
    subtitle = "Spearman correlations in the 45-patient anti-PD-1 cohort; * FDR<0.05",
    x = "Spearman ρ", y = NULL, colour = "ρ"
  ) + theme_pub(6.8)

sf6d <- ici_wilcox %>%
  left_join(ici_effects, by = "feature") %>%
  filter(status == "OK", is.finite(rank_biserial)) %>%
  mutate(
    neglogfdr = -log10(pmax(p_fdr, .Machine$double.xmin)),
    label = if_else(
      feature %in% c("target_expr", "Inflamed_5gene_exploratory", "IFNg_Response",
                     "CAF_Core", "ECM_Remodeling", "TGFb_Response"),
      feature_label(feature), ""
    ),
    direction = if_else(rank_biserial >= 0,
                        "Higher in responders", "Higher in nonresponders")
  )
write_panel_source(
  sf6d, "SUPPLEMENTARY", "SF6d_response_program_screen.csv", "FigS6d",
  "Part8_Immunotherapy/08_ORR_response_wilcoxon_all_features.csv",
  "All successful response-group comparisons converted to scale-free rank-biserial effects with bootstrap intervals",
  "Patient", "Exploratory multiplicity-controlled screen"
)
p_s6d <- ggplot(sf6d, aes(x = rank_biserial, y = neglogfdr)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.45) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = COL$grey, linewidth = 0.45) +
  geom_point(aes(colour = direction), size = 2.0, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.05, seed = PART9_SEED,
                           max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c("Higher in responders" = COL$red,
                                 "Higher in nonresponders" = COL$blue)) +
  labs(
    title = "Programme response-group screen",
    subtitle = "Scale-free rank-biserial effects; dashed line: FDR=0.05",
    x = "Rank-biserial effect (positive = higher in responders)",
    y = "-log10(FDR)", colour = NULL
  ) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  theme_pub(7.0) + theme(legend.position = "bottom")

fig_s6 <- (p_s6a | p_s6b) / (p_s6c | p_s6d) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_figure_bundle(
  fig_s6, "SUPPLEMENTARY", "SUPPLEMENTARY_FIGURE_S6_CNS_v8_10_0",
  FINAL_SUPP_WIDTH_IN, 6.20
)

# =============================================================================
# CURATED SUPPLEMENTARY TABLES
# These are publication-facing numerical tables, distinct from exact panel-source tables.
# =============================================================================

supp_table_contract <- tribble(
  ~table_id, ~title, ~part, ~file, ~scope,
  "Table S1", "TCGA-STAD cohort audit", "Part1", "Part1_TCGA_cohort_audit.csv",
  "Discovery-cohort composition and expression audit",
  "Table S2", "GEO bulk cohort audit", "Part2", "Part2_GEO_bulk_cohort_audit.csv",
  "External-cohort platforms, endpoints, completeness, exact OLFML2B mapping provenance, and explicit GSE84437 formal-source-versus-OS-evaluable counts",
  "Table S3", "Clinical-context expression summary", "Part3", "03c_olfml2b_expression_by_clinical_context.csv",
  "Clinical-context levels; stage rows use the canonical Figure 1c omnibus Kruskal-Wallis P value, while non-stage rows retain frozen Part3 global tests",
  "Table S4", "All survival models", "Part3", "05_olfml2b_survival_all_models.csv",
  "Univariable, adjusted, available-covariate and sensitivity Cox models",
  "Table S5", "Survival meta-analysis", "Part3", "08_olfml2b_meta_analysis.csv",
  "Random-effects, prediction intervals, heterogeneity and claim ceilings",
  "Table S6", "TME meta-correlations", "Part4", "10_olfml2b_tme_meta_correlations.csv",
  "Complete cross-cohort ecological programme synthesis",
  "Table S7", "TME-adjusted survival models", "Part4", "12_tme_attenuation_survival_models.csv",
  "Same-patient nested models, attenuation and diagnostics",
  "Table S8", "Protein direction tests", "Part5", "10_OLFML2B_case_paired_direction_tests.csv",
  "Same-plex case-paired protein inference",
  "Table S9", "Single-cell source evidence matrix", "Part6", "70_final_cell_source_evidence_matrix.csv",
  "Cross-dataset cell-source evidence and interpretation",
  "Table S10", "Single-cell/spatial concordance", "Part7", "80_Part6_Part7_concordance_matrix.csv",
  "Cross-modality source concordance without causal interpretation",
  "Table S11", "Immunotherapy publication summary", "Part8_Immunotherapy", "28_part8_primary_dataset_publication_summary.csv",
  "Single-cohort exploratory ICI results and claim ceiling"
)

# Publication-facing Table S3 stage-test reconciliation.
# IMPORTANT: this is deliberately confined to Part9 export.  The frozen Part3
# source table is not edited.  Figure 1c already computes the canonical stage
# omnibus test from the same patient-level Stage I-IV source rows.
canonical_stage_p <- stage_global %>%
  transmute(
    cohort_key = str_replace_all(as.character(cohort), "_", "-"),
    canonical_stage_global_p = as.numeric(kw_p),
    canonical_stage_fdr = as.numeric(kw_fdr)
  )

# Publication-facing GEO target-mapping provenance.
# Part2 writes this independently from the cohort-audit table; joining it here
# makes Table S2 self-contained without mutating any frozen Part2 output.
geo_target_mapping_audit <- read_contract(
  "Part2", "Part2_GEO_OLFML2B_mapping_audit.csv",
  c("probe_id", "symbol", "primary_eligible", "cohort")
)

geo_target_mapping_summary <- geo_target_mapping_audit %>%
  mutate(
    cohort = as.character(cohort),
    probe_id = as.character(probe_id),
    symbol = toupper(str_trim(as.character(symbol))),
    primary_eligible = as.logical(primary_eligible)
  ) %>%
  filter(
    symbol == "OLFML2B",
    primary_eligible %in% TRUE,
    !is.na(probe_id),
    nzchar(str_trim(probe_id))
  ) %>%
  group_by(cohort) %>%
  summarise(
    publication_olfml2b_probe_n = n_distinct(probe_id),
    publication_olfml2b_probe_ids = paste(sort(unique(probe_id)), collapse = ";"),
    publication_olfml2b_mapping_sources = if (
      "mapping_source" %in% names(cur_data_all())
    ) {
      paste(
        sort(unique(
          as.character(cur_data_all()$mapping_source)[
            !is.na(cur_data_all()$mapping_source) &
              nzchar(str_trim(as.character(cur_data_all()$mapping_source)))
          ]
        )),
        collapse = ";"
      )
    } else {
      NA_character_
    },
    publication_olfml2b_mapping_status = "PASS_EXACT_PRIMARY_ELIGIBLE_OLFML2B",
    .groups = "drop"
  )

supp_table_manifest <- pmap_dfr(
  supp_table_contract,
  function(table_id, title, part, file, scope) {
    dat <- read_contract(part, file, character())

    if (identical(table_id, "Table S2")) {
      required_s2 <- c(
        "cohort", "n_samples_expression", "n_samples_clinical",
        "n_os_complete", "os_events"
      )
      absent_s2 <- setdiff(required_s2, names(dat))
      if (length(absent_s2) > 0L) {
        stop(
          "Table S2 provenance reconciliation failed; absent columns: ",
          paste(absent_s2, collapse = ", ")
        )
      }

      dat_before_s2 <- dat

      dat <- dat %>%
        mutate(cohort = as.character(cohort)) %>%
        left_join(geo_target_mapping_summary, by = "cohort")

      # Where the legacy cohort-audit fields are blank, recover the exact target
      # probe count/IDs from the frozen Part2 mapping audit. Existing nonmissing
      # values are preserved.
      if ("n_olfml2b_probes" %in% names(dat)) {
        old_n <- suppressWarnings(as.numeric(dat$n_olfml2b_probes))
        mapped_n <- suppressWarnings(as.numeric(dat$publication_olfml2b_probe_n))
        dat$n_olfml2b_probes <- ifelse(
          is.finite(old_n), old_n, mapped_n
        )
      }
      if ("olfml2b_probes" %in% names(dat)) {
        old_probe <- as.character(dat$olfml2b_probes)
        mapped_probe <- as.character(dat$publication_olfml2b_probe_ids)
        use_old <- !is.na(old_probe) & nzchar(str_trim(old_probe))
        dat$olfml2b_probes <- ifelse(use_old, old_probe, mapped_probe)
      }

      # Reader-facing endpoint/source definitions. These columns clarify that
      # source-subseries membership and endpoint evaluability are distinct.
      dat$formal_source_subseries_n <- NA_integer_
      dat$os_evaluable_n <- suppressWarnings(as.integer(dat$n_os_complete))
      dat$os_evaluable_events <- suppressWarnings(as.integer(dat$os_events))
      dat$context_only_excluded_n <- NA_integer_
      dat$formal_source_endpoint_note <- NA_character_

      g844 <- dat$cohort == "GSE84437"
      if (sum(g844) != 1L) {
        stop(
          "Table S2 provenance reconciliation failed: expected exactly one ",
          "GSE84437 row; observed ", sum(g844), "."
        )
      }

      # These fields were explicitly frozen by Part2's formal GSE84437 contract.
      if (!all(c("n_formal_os", "n_context_only_excluded") %in% names(dat))) {
        stop(
          "Table S2 provenance reconciliation failed: GSE84437 formal-source ",
          "audit fields n_formal_os/n_context_only_excluded are absent."
        )
      }

      dat$formal_source_subseries_n[g844] <- suppressWarnings(
        as.integer(dat$n_formal_os[g844])
      )
      dat$context_only_excluded_n[g844] <- suppressWarnings(
        as.integer(dat$n_context_only_excluded[g844])
      )
      dat$formal_source_endpoint_note[g844] <- paste0(
        "Formal 2016 source subseries: ",
        dat$formal_source_subseries_n[g844],
        " patients/samples; complete OS: ",
        dat$os_evaluable_n[g844],
        " patients (", dat$os_evaluable_events[g844],
        " events); context-only superseries samples excluded: ",
        dat$context_only_excluded_n[g844], "."
      )

      # Hard publication guard for the frozen GSE84437 contract.
      g844_guard <- tibble(
        criterion = c(
          "formal_source_subseries_n",
          "os_evaluable_n",
          "os_events",
          "context_only_excluded_n",
          "expression_rows_equal_formal_source",
          "clinical_rows_equal_formal_source",
          "exact_olfml2b_mapping_available",
          "legacy_probe_fields_publication_filled"
        ),
        observed = c(
          as.character(dat$formal_source_subseries_n[g844]),
          as.character(dat$os_evaluable_n[g844]),
          as.character(dat$os_evaluable_events[g844]),
          as.character(dat$context_only_excluded_n[g844]),
          as.character(dat$n_samples_expression[g844]),
          as.character(dat$n_samples_clinical[g844]),
          as.character(dat$publication_olfml2b_probe_ids[g844]),
          if (
            "olfml2b_probes" %in% names(dat)
          ) as.character(dat$olfml2b_probes[g844]) else NA_character_
        ),
        expected = c(
          "433", "431", "207", "50", "433", "433",
          "non-empty exact primary-eligible OLFML2B probe ID(s)",
          "non-empty"
        ),
        pass = c(
          dat$formal_source_subseries_n[g844] == 433L,
          dat$os_evaluable_n[g844] == 431L,
          dat$os_evaluable_events[g844] == 207L,
          dat$context_only_excluded_n[g844] == 50L,
          suppressWarnings(as.integer(dat$n_samples_expression[g844])) == 433L,
          suppressWarnings(as.integer(dat$n_samples_clinical[g844])) == 433L,
          !is.na(dat$publication_olfml2b_probe_ids[g844]) &
            nzchar(str_trim(dat$publication_olfml2b_probe_ids[g844])) &
            dat$publication_olfml2b_mapping_status[g844] ==
              "PASS_EXACT_PRIMARY_ELIGIBLE_OLFML2B",
          if ("olfml2b_probes" %in% names(dat)) {
            !is.na(dat$olfml2b_probes[g844]) &
              nzchar(str_trim(as.character(dat$olfml2b_probes[g844])))
          } else {
            FALSE
          }
        )
      )

      if (any(!g844_guard$pass)) {
        stop(
          "Table S2 GSE84437 publication contract failed: ",
          paste(g844_guard$criterion[!g844_guard$pass], collapse = ", ")
        )
      }

      # Guard: no pre-existing non-GSE84437 endpoint/sample-count field may be
      # changed by the publication-facing reconciliation.
      invariant_cols <- intersect(
        c(
          "n_samples_expression", "n_samples_clinical", "n_os_complete",
          "os_events", "n_formal_os", "n_context_only_excluded"
        ),
        names(dat_before_s2)
      )
      for (cc in invariant_cols) {
        before <- dat_before_s2[[cc]]
        after <- dat[[cc]]
        if (!identical(before, after)) {
          stop(
            "Table S2 reconciliation contract failed: frozen column '", cc,
            "' changed."
          )
        }
      }

      readr::write_csv(
        g844_guard,
        file.path(
          AUDIT_DIR,
          "02B_TABLE_S2_GSE84437_PROVENANCE_RECONCILIATION_v8_10_3.csv"
        ),
        na = ""
      )
    }

    if (identical(table_id, "Table S3")) {
      required_s3 <- c("cohort", "variable", "level", "n", "median", "mean", "global_p")
      absent_s3 <- setdiff(required_s3, names(dat))
      if (length(absent_s3) > 0L) {
        stop(
          "Table S3 stage canonicalisation failed; absent columns: ",
          paste(absent_s3, collapse = ", ")
        )
      }

      dat_before <- dat
      dat <- dat %>%
        mutate(
          .row_id_s3 = row_number(),
          .cohort_key_s3 = str_replace_all(as.character(cohort), "_", "-"),
          .variable_lower_s3 = str_to_lower(str_trim(as.character(variable)))
        ) %>%
        left_join(canonical_stage_p, by = c(".cohort_key_s3" = "cohort_key"))

      stage_idx <- dat$.variable_lower_s3 == "stage"

      if (!any(stage_idx)) {
        stop("Table S3 stage canonicalisation failed: no stage rows were found.")
      }
      if (any(!is.finite(dat$canonical_stage_global_p[stage_idx]))) {
        bad <- unique(dat$.cohort_key_s3[
          stage_idx & !is.finite(dat$canonical_stage_global_p)
        ])
        stop(
          "Table S3 stage canonicalisation failed: Figure 1c Kruskal-Wallis P ",
          "value unavailable for cohort(s): ", paste(bad, collapse = ", ")
        )
      }

      # ONLY stage global_p is replaced.  Non-stage values are untouched.
      dat$global_p[stage_idx] <- dat$canonical_stage_global_p[stage_idx]

      # Reader-facing provenance note.  This does not alter any non-stage statistic.
      dat$global_test_method <- ifelse(
        stage_idx,
        "Kruskal-Wallis omnibus test across Stage I-IV",
        "Frozen Part3 clinical-context global test"
      )
      dat$global_p_source <- ifelse(
        stage_idx,
        "Part9 Figure 1c canonical stage test; same Stage I-IV patient set",
        "Part3/03c_olfml2b_expression_by_clinical_context.csv; unchanged"
      )

      # Guard 1: all non-stage global_p values must remain exactly unchanged.
      nonstage_idx_before <- str_to_lower(
        str_trim(as.character(dat_before$variable))
      ) != "stage"
      if (!identical(
        as.numeric(dat$global_p[!stage_idx]),
        as.numeric(dat_before$global_p[nonstage_idx_before])
      )) {
        stop(
          "Table S3 canonicalisation contract failed: a non-stage global_p ",
          "value changed."
        )
      }

      # Guard 2: every stage row within a cohort must equal Figure 1c kw_p.
      stage_check <- dat %>%
        filter(stage_idx) %>%
        transmute(
          cohort,
          level,
          n,
          exported_global_p = as.numeric(global_p),
          figure1c_kw_p = as.numeric(canonical_stage_global_p),
          figure1c_kw_fdr = as.numeric(canonical_stage_fdr),
          exact_match = is.finite(exported_global_p) &
            is.finite(figure1c_kw_p) &
            abs(exported_global_p - figure1c_kw_p) <=
              1e-12 * pmax(1, abs(figure1c_kw_p))
        )

      if (any(!stage_check$exact_match)) {
        stop(
          "Table S3 canonicalisation contract failed: exported stage global_p ",
          "does not equal Figure 1c Kruskal-Wallis P."
        )
      }

      readr::write_csv(
        stage_check,
        file.path(
          AUDIT_DIR,
          "02A_TABLE_S3_STAGE_GLOBAL_P_CANONICALISATION_v8_10_2.csv"
        ),
        na = ""
      )

      dat <- dat %>%
        select(
          -.row_id_s3,
          -.cohort_key_s3,
          -.variable_lower_s3,
          -canonical_stage_global_p,
          -canonical_stage_fdr
        )
    }

    suffix <- str_replace_all(table_id, " ", "_")
    out_name <- paste0(
      suffix, "_",
      str_replace_all(title, "[^A-Za-z0-9]+", "_"),
      ".csv"
    )
    out_path <- file.path(SUPP_TABLE_DIR, out_name)
    readr::write_csv(dat, out_path, na = "")

    tibble(
      table_id, title, output_file = out_name,
      source_input = if (
        identical(table_id, "Table S2")
      ) {
        paste0(
          part, "/", file,
          " + Part2/Part2_GEO_OLFML2B_mapping_audit.csv ",
          "(publication-facing GSE84437 count/mapping reconciliation only)"
        )
      } else if (
        identical(table_id, "Table S3")
      ) {
        paste0(
          part, "/", file,
          " + Part9 Figure 1c stage_global (Kruskal-Wallis canonicalisation only)"
        )
      } else {
        paste(part, file, sep = "/")
      },
      n_rows = nrow(dat), n_columns = ncol(dat), scope
    )
  }
)

claim_tables <- list(
  "Part3/21_olfml2b_interpretation_boundary.csv" = read_contract(
    "Part3", "21_olfml2b_interpretation_boundary.csv", character()
  ),
  "Part4/21_part4_methodology_claim_limits.csv" = read_contract(
    "Part4", "21_part4_methodology_claim_limits.csv", character()
  ),
  "Part6/103_OLFML2B_Part6_unified_claim_boundary.csv" = read_contract(
    "Part6", "103_OLFML2B_Part6_unified_claim_boundary.csv", character()
  ),
  "Part7/100_Part7_spatial_claim_boundary.csv" = read_contract(
    "Part7", "100_Part7_spatial_claim_boundary.csv", character()
  ),
  "Part8_Immunotherapy/27_final_molecular_context_claim_ceiling.csv" = read_contract(
    "Part8_Immunotherapy", "27_final_molecular_context_claim_ceiling.csv", character()
  )
)
# Reader-facing normalization of heterogeneous internal claim contracts.
# The source files have different schemas by design; the old raw row-bind created
# a sparse 13-column table.  This normalization does not infer new claims.
# It maps source-native fields into a common publication-facing vocabulary.
claim_boundary_reader <- dplyr::bind_rows(
  claim_tables[["Part3/21_olfml2b_interpretation_boundary.csv"]] %>%
    transmute(
      evidence_layer = "Bulk expression / survival / ecological context",
      domain = as.character(domain),
      boundary_class = as.character(manuscript_use),
      publication_statement = as.character(allowed_claim),
      prohibited_or_unsupported_claim = NA_character_,
      source_input = "Part3/21_olfml2b_interpretation_boundary.csv"
    ),
  claim_tables[["Part4/21_part4_methodology_claim_limits.csv"]] %>%
    transmute(
      evidence_layer = "Bulk tumour-microenvironment programmes",
      domain = as.character(analysis_block),
      boundary_class = "ALLOWED_WITH_EXPLICIT_PROHIBITION",
      publication_statement = as.character(allowed_claim),
      prohibited_or_unsupported_claim = as.character(prohibited_claim),
      source_input = "Part4/21_part4_methodology_claim_limits.csv"
    ),
  claim_tables[["Part6/103_OLFML2B_Part6_unified_claim_boundary.csv"]] %>%
    transmute(
      evidence_layer = "Single-cell transcriptomics",
      domain = "single_cell_claim_boundary",
      boundary_class = as.character(level),
      publication_statement = if_else(
        toupper(as.character(level)) == "FORBIDDEN",
        NA_character_,
        as.character(statement)
      ),
      prohibited_or_unsupported_claim = if_else(
        toupper(as.character(level)) == "FORBIDDEN",
        as.character(statement),
        NA_character_
      ),
      source_input = "Part6/103_OLFML2B_Part6_unified_claim_boundary.csv"
    ),
  claim_tables[["Part7/100_Part7_spatial_claim_boundary.csv"]] %>%
    transmute(
      evidence_layer = "Spatial transcriptomics",
      domain = as.character(domain),
      boundary_class = case_when(
        as.character(domain) == "final_gene_lock" ~ "INTERNAL_CONTRACT_FLAG",
        str_detect(str_to_lower(as.character(decision)), "not supported") ~
          "NOT_SUPPORTED",
        TRUE ~ "BOUNDARY_OR_ANALYSIS_RULE"
      ),
      publication_statement = case_when(
        as.character(domain) == "final_gene_lock" ~
          paste0(
            "Internal source-contract flag retained for provenance: ",
            "final_gene_lock = ", as.character(decision), "."
          ),
        str_detect(str_to_lower(as.character(decision)), "not supported") ~
          NA_character_,
        TRUE ~ as.character(decision)
      ),
      prohibited_or_unsupported_claim = if_else(
        str_detect(str_to_lower(as.character(decision)), "not supported"),
        as.character(decision),
        NA_character_
      ),
      source_input = "Part7/100_Part7_spatial_claim_boundary.csv"
    ),
  claim_tables[["Part8_Immunotherapy/27_final_molecular_context_claim_ceiling.csv"]] %>%
    transmute(
      evidence_layer = "Anti-PD-1 response boundary",
      domain = as.character(claim_domain),
      boundary_class = as.character(status),
      publication_statement = case_when(
        as.character(claim_domain) == "molecular_annotation" &
          as.character(status) == "NOT_AVAILABLE" ~
            paste(
              "Molecular annotations required for the intended contextual",
              "stratification were not available."
            ),
        as.character(claim_domain) == "known_ici_biomarker_context" &
          as.character(status) == "KNOWN_MARKER_SIGNAL_WEAK_OR_NOT_DETECTED" ~
            paste(
              "Established ICI-marker signals were weak or not detected in",
              "this cohort; independence from standard biomarkers was not established."
            ),
        as.character(claim_domain) == "OLFML2B_standalone_ORR_prediction" &
          as.character(status) == "NO_NOMINAL_ASSOCIATION" ~
            paste(
              "No nominal standalone association between OLFML2B and objective",
              "response was observed in this cohort."
            ),
        as.character(claim_domain) == "intended_MSS_EBVneg_TMBlow_population" &
          as.character(status) == "UNDERPOWERED_OR_NOT_EVALUABLE" ~
            paste(
              "The intended MSS/EBV-negative/TMB-low subgroup analysis was",
              "underpowered or not evaluable."
            ),
        as.character(claim_domain) == "recommended_claim_ceiling" ~
          as.character(status),
        TRUE ~ paste0("Source-contract status: ", as.character(status), ".")
      ),
      prohibited_or_unsupported_claim = case_when(
        as.character(claim_domain) == "OLFML2B_standalone_ORR_prediction" ~
          "Do not claim OLFML2B as a validated standalone anti-PD-1 response biomarker.",
        as.character(claim_domain) == "known_ici_biomarker_context" ~
          "Do not claim independence from established ICI biomarkers.",
        as.character(claim_domain) == "intended_MSS_EBVneg_TMBlow_population" ~
          "Do not make a subgroup-specific predictive claim from the intended MSS/EBV-negative/TMB-low analysis.",
        TRUE ~ NA_character_
      ),
      source_input =
        "Part8_Immunotherapy/27_final_molecular_context_claim_ceiling.csv"
    )
) %>%
  mutate(
    across(
      c(
        evidence_layer, domain, boundary_class, publication_statement,
        prohibited_or_unsupported_claim, source_input
      ),
      ~na_if(str_trim(as.character(.x)), "")
    )
  )

if (nrow(claim_boundary_reader) != sum(vapply(claim_tables, nrow, integer(1)))) {
  stop(
    "Table S12 reader-view contract failed: normalized row count differs from ",
    "the sum of frozen source claim-contract rows."
  )
}
if (any(is.na(claim_boundary_reader$source_input))) {
  stop("Table S12 reader-view contract failed: source_input is missing.")
}

# Fail closed if a machine-facing Boolean/status code remains as the complete
# reader-facing publication statement.  boundary_class may retain source codes.
machine_only_publication_values <- c(
  "TRUE", "FALSE", "NOT_AVAILABLE",
  "KNOWN_MARKER_SIGNAL_WEAK_OR_NOT_DETECTED",
  "NO_NOMINAL_ASSOCIATION", "UNDERPOWERED_OR_NOT_EVALUABLE"
)
bad_machine_statement <- !is.na(claim_boundary_reader$publication_statement) &
  str_trim(claim_boundary_reader$publication_statement) %in%
    machine_only_publication_values

if (any(bad_machine_statement)) {
  stop(
    "Table S12 publication-language QA failed; machine-facing statement(s) ",
    "remain in row(s): ",
    paste(which(bad_machine_statement), collapse = ", ")
  )
}

s12_language_audit <- claim_boundary_reader %>%
  transmute(
    row_id = row_number(),
    evidence_layer,
    domain,
    boundary_class,
    publication_statement,
    prohibited_or_unsupported_claim,
    source_input,
    publication_statement_is_machine_only =
      !is.na(publication_statement) &
      str_trim(publication_statement) %in% machine_only_publication_values,
    source_provenance_present =
      !is.na(source_input) & nzchar(str_trim(source_input))
  )

readr::write_csv(
  s12_language_audit,
  file.path(
    AUDIT_DIR,
    "02D_TABLE_S12_PUBLICATION_LANGUAGE_QA_v8_10_5.csv"
  ),
  na = ""
)

claim_out <- file.path(SUPP_TABLE_DIR, "Table_S12_Integrated_claim_boundaries.csv")
readr::write_csv(claim_boundary_reader, claim_out, na = "")

supp_table_manifest <- dplyr::bind_rows(
  supp_table_manifest,
  tibble(
    table_id = "Table S12", title = "Integrated claim boundaries",
    output_file = basename(claim_out),
    source_input = paste(names(claim_tables), collapse = "; "),
    n_rows = nrow(claim_boundary_reader), n_columns = ncol(claim_boundary_reader),
    scope = paste(
      "Reader-facing allowed claims, explicit limitations and interpretation",
      "ceilings across Parts 3-8; all source statements retain source provenance"
    )
  )
)

# ---------------------------------------------------------------------------
# Supplementary-table presentation contract: titles, notes, missing-value rules,
# and abbreviations.  These sidecars are publication metadata and do not alter
# the machine-readable numerical tables.
# ---------------------------------------------------------------------------
supp_table_presentation <- tribble(
  ~table_id, ~title, ~note, ~missing_value_rule, ~abbreviations,

  "Table S1", "TCGA-STAD cohort audit",
  paste(
    "Cohort composition and target-expression audit for TCGA-STAD.",
    "Counts refer to the frozen evaluable data used by the analysis pipeline;",
    "target identifiers and locus information are included for reproducibility."
  ),
  "Blank fields indicate information not applicable to this one-row TCGA audit.",
  "OS, overall survival; TCGA-STAD, The Cancer Genome Atlas stomach adenocarcinoma; TPM, transcripts per million; HGNC, HUGO Gene Nomenclature Committee.",

  "Table S2", "GEO bulk cohort audit",
  paste(
    "Platform, sample, endpoint, target-probe and provenance audit for the four",
    "external GEO bulk cohorts. For GSE84437, 433 denotes the formal source",
    "subseries, 431 denotes patients with complete OS data, 207 denotes OS events,",
    "and 50 context-only superseries samples were excluded from formal inference."
  ),
  paste(
    "Blank cells mean that the corresponding field was structurally unavailable,",
    "not applicable to that cohort, or not used for formal inference; blanks must",
    "not be interpreted as zero."
  ),
  "DFS, disease-free survival; GEO, Gene Expression Omnibus; OS, overall survival; RFS, recurrence-free survival.",

  "Table S3", "Clinical-context expression summary",
  paste(
    "OLFML2B expression summaries across available clinical strata.",
    "Stage global P values are omnibus Kruskal-Wallis tests across Stage I-IV and",
    "are repeated across stage rows only to retain a rectangular table.",
    "Non-stage global tests retain the frozen Part3 clinical-context tests."
  ),
  "Blank fields indicate unavailable clinical strata or non-applicable metadata.",
  "FDR, false discovery rate; NA, not available.",

  "Table S4", "All survival models",
  paste(
    "Complete univariable, common-adjustment, available-covariate and sensitivity",
    "Cox models for OLFML2B. Continuous effects are reported per one within-cohort",
    "standard-deviation increase in OLFML2B expression."
  ),
  paste(
    "Blank diagnostic or imputation fields indicate that the item was not",
    "applicable, not estimable, or intentionally not used for that model."
  ),
  "CI, confidence interval; DFS, disease-free survival; HR, hazard ratio; OS, overall survival; PH, proportional hazards; RFS, recurrence-free survival.",

  "Table S5", "Survival meta-analysis",
  paste(
    "Random-effects meta-analyses of structurally eligible cohort-specific Cox",
    "models, including heterogeneity, prediction intervals, and REML cross-checks."
  ),
  "Blank values denote quantities that were not estimable or not applicable for the specified small-k synthesis.",
  "CI, confidence interval; HR, hazard ratio; I2, I-squared heterogeneity statistic; REML, restricted maximum likelihood.",

  "Table S6", "TME meta-correlations",
  paste(
    "Cross-cohort random-effects synthesis of Spearman correlations between",
    "OLFML2B and prespecified tumour-microenvironment programme scores."
  ),
  "Blank values indicate non-estimable quantities for the relevant programme synthesis.",
  "CAF, cancer-associated fibroblast; ECM, extracellular matrix; FDR, false discovery rate; TME, tumour microenvironment.",

  "Table S7", "TME-adjusted survival models",
  paste(
    "Same-patient nested Cox models assessing statistical attenuation of the",
    "OLFML2B coefficient after adding prespecified ecological axes.",
    "Attenuation is a sensitivity analysis and is not a causal mediated proportion."
  ),
  paste(
    "Blank cells indicate unavailable, non-estimable or structurally non-applicable",
    "diagnostics; they are not zeros."
  ),
  "AIC, Akaike information criterion; CAF, cancer-associated fibroblast; CI, confidence interval; ECM, extracellular matrix; FDR, false discovery rate; HR, hazard ratio; PH, proportional hazards; TME, tumour microenvironment; VIF, variance inflation factor.",

  "Table S8", "Protein direction tests",
  paste(
    "Case-paired OLFML2B protein direction analyses in PDC000614, including",
    "bootstrap uncertainty, exact sign tests, analytical-plex direction and",
    "leave-one-plex-out sensitivity."
  ),
  "Blank fields indicate a statistic not estimable or not applicable to the stated value type.",
  "CI, confidence interval; PDC, Proteomic Data Commons.",

  "Table S9", "Single-cell source evidence matrix",
  paste(
    "Cross-dataset evidence matrix used to rank preferential cellular localisation",
    "of OLFML2B. Evidence grades integrate dataset eligibility, top-rank recurrence,",
    "leave-one-dataset-out stability, threshold consistency and paired contrasts."
  ),
  "Zero indicates no qualifying evidence under the frozen rule; blank would indicate not evaluated.",
  "CAF, cancer-associated fibroblast; ECM, extracellular matrix; LODO, leave-one-dataset-out.",

  "Table S10", "Single-cell/spatial concordance",
  paste(
    "Cross-modality concordance between single-cell source expectations and",
    "patient-level spatial same-spot, neighbourhood and competing-source analyses.",
    "The table supports spatial association only, not cell identity or causality."
  ),
  "Blank values indicate that the corresponding spatial summary was unavailable or not estimable.",
  "CAF, cancer-associated fibroblast; ECM, extracellular matrix.",

  "Table S11", "Immunotherapy publication summary",
  paste(
    "Patient-level exploratory anti-PD-1 results in PRJEB25780/TIGER.",
    "This single cohort defines the boundary of response-prediction claims and is",
    "not an independent validation of a predictive biomarker."
  ),
  paste(
    "Blank values denote molecular annotations or model quantities that were",
    "unavailable or not estimable; blanks are not negative results."
  ),
  "AUC, area under the receiver-operating-characteristic curve; CR, complete response; DCR, disease control rate; ICI, immune-checkpoint inhibitor; NR, non-responder; ORR, objective response rate; PD, progressive disease; PR, partial response; SD, stable disease.",

  "Table S12", "Integrated claim boundaries",
  paste(
    "Reader-facing synthesis of the explicit claim ceilings frozen across Parts 3-8.",
    "Rows preserve source provenance and distinguish supported publication language",
    "from prohibited, unsupported or not-evaluable claims."
  ),
  "Blank supported/prohibited cells mean that the source contract did not specify that side of the boundary.",
  "CAF, cancer-associated fibroblast; ECM, extracellular matrix; ICI, immune-checkpoint inhibitor; ORR, objective response rate; TME, tumour microenvironment."
)

if (nrow(supp_table_presentation) != 12L ||
    anyDuplicated(supp_table_presentation$table_id)) {
  stop("Supplementary-table presentation contract must contain exactly one row for Tables S1-S12.")
}

# Attach publication metadata to the manifest.
supp_table_manifest <- supp_table_manifest %>%
  left_join(supp_table_presentation, by = c("table_id", "title"))

if (any(is.na(supp_table_manifest$note)) ||
    any(is.na(supp_table_manifest$missing_value_rule)) ||
    any(is.na(supp_table_manifest$abbreviations))) {
  stop("Supplementary-table presentation metadata is incomplete after manifest join.")
}

presentation_out <- file.path(
  SUPP_TABLE_DIR,
  "00_Supplementary_Table_Titles_Notes_Abbreviations.csv"
)
readr::write_csv(supp_table_presentation, presentation_out, na = "")

readr::write_csv(
  supp_table_presentation,
  file.path(
    AUDIT_DIR,
    "02C_SUPPLEMENTARY_TABLE_PRESENTATION_CONTRACT_v8_10_4.csv"
  ),
  na = ""
)

# One human-readable sidecar note per table for deterministic ESM workbook assembly.
for (i in seq_len(nrow(supp_table_presentation))) {
  rr <- supp_table_presentation[i, ]
  note_lines <- c(
    paste0(rr$table_id, ". ", rr$title),
    "",
    paste0("Note: ", rr$note),
    paste0("Missing-value rule: ", rr$missing_value_rule),
    paste0("Abbreviations: ", rr$abbreviations)
  )
  note_file <- file.path(
    SUPP_TABLE_DIR,
    paste0(str_replace_all(rr$table_id, " ", "_"), "_NOTE.txt")
  )
  writeLines(note_lines, note_file, useBytes = TRUE)
}

readr::write_csv(
  supp_table_manifest,
  file.path(AUDIT_DIR, "02_SUPPLEMENTARY_TABLE_MANIFEST_v8_10_0.csv"),
  na = ""
)

# =============================================================================
# AUDITS, PROVENANCE, AND COMPLETION SENTINEL
# =============================================================================

source_registry <- source_registry %>%
  arrange(branch, panel_id, source_table_file)
readr::write_csv(
  source_registry,
  file.path(AUDIT_DIR, "03_FIGURE_SOURCE_TABLE_REGISTRY_v8_10_0.csv"),
  na = ""
)

main_panel_coverage <- tibble(panel_id = main_budget$panel_id) %>%
  left_join(
    source_registry %>%
      filter(branch == "MAIN", !str_detect(panel_id, "statistics")) %>%
      distinct(panel_id) %>% mutate(source_registered = TRUE),
    by = "panel_id"
  ) %>%
  mutate(
    source_registered = source_registered %in% TRUE,
    status = if_else(source_registered, "PASS", "FAIL")
  )
unexpected_main_sources <- setdiff(
  source_registry$panel_id[
    source_registry$branch == "MAIN" & !str_detect(source_registry$panel_id, "statistics")
  ],
  main_budget$panel_id
)
if (any(main_panel_coverage$status != "PASS") || length(unexpected_main_sources) > 0L) {
  stop(
    "Main panel/source-table coverage failed. Missing or unexpected IDs: ",
    paste(c(main_panel_coverage$panel_id[main_panel_coverage$status != "PASS"],
            unexpected_main_sources), collapse = ", ")
  )
}
readr::write_csv(
  main_panel_coverage,
  file.path(AUDIT_DIR, "03B_MAIN_PANEL_SOURCE_COVERAGE_v8_10_0.csv"),
  na = ""
)

nonduplication_audit <- tribble(
  ~domain, ~main_content, ~supplementary_content, ~status, ~reason,
  "Expression/clinical",
  "Unpaired expression, paired patient change and cross-cohort stage context",
  "Non-stage clinical screens and model diagnostics",
  "PASS", "The low-information evidence path is deleted and the promoted stage panel is removed from Fig. S1",
  "Survival",
  "Selected adjusted cohort estimates and full random-effects synthesis",
  "LOCO influence, spline and time-varying HR diagnostics",
  "PASS", "LOCO is supplementary and the main forest is not repeated",
  "Bulk TME",
  "Decile shape, six-program ecological map, pooled synthesis and composite axes",
  "Secondary programmes, collinearity and attenuation diagnostics",
  "PASS", "Promoted decile display is removed from Fig. S2; attenuation is supplementary",
  "Single cell",
  "Localization, patient pairing, coupling, positive-cell context and complete LODO matrix",
  "Within-compartment orthogonal states",
  "PASS", "The complete LODO matrix is shown once in Main Fig. 3e",
  "Spatial",
  "All-section atlas, section detection, concordance, threshold overlap, representative joint map and ridge source boundary",
  "Global Moran's I, distance and technical-residual attenuation",
  "PASS", "Section detection is promoted once to Main Fig. 4b and removed from Fig. S4",
  "Protein",
  "Case-level paired distribution, plex stability and orientation sentinels",
  "Pair eligibility",
  "PASS", "Protein orientation and plex displays occur once in Main Fig. 5",
  "Immunotherapy",
  "Not included in main figures",
  "Entire exploratory single-cohort module",
  "PASS", "Evidence strength is below the main-figure threshold"
)
readr::write_csv(
  nonduplication_audit,
  file.path(AUDIT_DIR, "04_MAIN_SUPPLEMENT_NONDUPLICATION_AUDIT_v8_10_0.csv")
)

input_rel <- unique(trimws(unlist(strsplit(source_registry$source_inputs, ";", fixed = TRUE))))
input_rel <- input_rel[nzchar(input_rel)]
input_abs <- file.path(TABLE_ROOT, input_rel)
input_manifest <- tibble(
  relative_path = input_rel,
  exists = file.exists(input_abs),
  size_bytes = ifelse(file.exists(input_abs), file.info(input_abs)$size, NA_real_),
  md5 = ifelse(file.exists(input_abs), unname(tools::md5sum(input_abs)), NA_character_)
)
if (any(!input_manifest$exists)) {
  stop("Figure-source provenance contains missing input files: ",
       paste(input_manifest$relative_path[!input_manifest$exists], collapse = ", "))
}
readr::write_csv(
  input_manifest,
  file.path(AUDIT_DIR, "05_FIGURE_INPUT_MD5_MANIFEST_v8_10_0.csv"),
  na = ""
)

expected_main <- unlist(lapply(1:5, function(i) {
  stem <- file.path(MAIN_FIG_DIR, sprintf("MAIN_FIGURE_%02d_CNS_v8_10_0", i))
  paste0(stem, c(".pdf", ".png", ".tiff"))
}))
expected_supp <- unlist(lapply(1:6, function(i) {
  stem <- file.path(SUPP_FIG_DIR, sprintf("SUPPLEMENTARY_FIGURE_S%d_CNS_v8_10_0", i))
  paste0(stem, c(".pdf", ".png", ".tiff"))
}))
expected_all <- c(expected_main, expected_supp)
export_audit <- tibble(
  file = basename(expected_all),
  branch = c(rep("MAIN", length(expected_main)), rep("SUPPLEMENTARY", length(expected_supp))),
  exists = file.exists(expected_all),
  size_bytes = ifelse(file.exists(expected_all), file.info(expected_all)$size, NA_real_),
  status = if_else(exists & size_bytes > 0, "PASS", "FAIL")
)
readr::write_csv(
  export_audit,
  file.path(AUDIT_DIR, "06_FIGURE_EXPORT_AUDIT_v8_10_0.csv"),
  na = ""
)
if (any(export_audit$status != "PASS")) {
  stop("One or more expected figure files are absent or empty. See export audit.")
}

figure_source_paths <- ifelse(
  source_registry$branch == "MAIN",
  file.path(MAIN_SOURCE_DIR, source_registry$source_table_file),
  file.path(SUPP_SOURCE_DIR, source_registry$source_table_file)
)
supplementary_table_paths <- file.path(SUPP_TABLE_DIR, supp_table_manifest$output_file)
table_export_audit <- dplyr::bind_rows(
  tibble(
    output_class = "FIGURE_SOURCE_TABLE",
    file = basename(figure_source_paths),
    exists = file.exists(figure_source_paths),
    size_bytes = ifelse(file.exists(figure_source_paths), file.info(figure_source_paths)$size, NA_real_)
  ),
  tibble(
    output_class = "SUPPLEMENTARY_TABLE",
    file = basename(supplementary_table_paths),
    exists = file.exists(supplementary_table_paths),
    size_bytes = ifelse(
      file.exists(supplementary_table_paths),
      file.info(supplementary_table_paths)$size,
      NA_real_
    )
  )
) %>% mutate(status = if_else(exists & size_bytes > 0, "PASS", "FAIL"))
readr::write_csv(
  table_export_audit,
  file.path(AUDIT_DIR, "07_TABLE_EXPORT_AUDIT_v8_10_0.csv"),
  na = ""
)
if (any(table_export_audit$status != "PASS")) {
  stop("One or more figure-source or supplementary tables are absent or empty. See table export audit.")
}

parse_audit_out <- if (nrow(READ_PARSE_AUDIT) == 0L) {
  tibble(
    source_file = NA_character_, row = NA_integer_, col = NA_integer_,
    expected = NA_character_, actual = NA_character_, file = NA_character_,
    status = "PASS_NO_PARSE_PROBLEMS"
  )
} else {
  READ_PARSE_AUDIT %>% mutate(status = "REVIEW_SOURCE_CELL")
}
readr::write_csv(
  parse_audit_out,
  PARSE_AUDIT_PATH,
  na = ""
)

manuscript_docx <- Sys.getenv("OLFML2B_MANUSCRIPT_DOCX", unset = "")
enforce_manuscript_audit <- identical(
  toupper(Sys.getenv("OLFML2B_ENFORCE_MANUSCRIPT_MEDIA_AUDIT", unset = "FALSE")),
  "TRUE"
)
if (nzchar(manuscript_docx) && file.exists(manuscript_docx)) {
  manuscript_docx <- normalizePath(manuscript_docx, winslash = "/", mustWork = TRUE)
  docx_listing <- utils::unzip(manuscript_docx, list = TRUE)
  media_names <- docx_listing$Name[str_detect(docx_listing$Name, "^word/media/[^/]+$")]
  if (length(media_names) < 1L) stop("Manuscript media audit: no embedded word/media files found.")
  media_dir <- tempfile("part9_docx_media_")
  dir.create(media_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(manuscript_docx, files = media_names, exdir = media_dir)
  media_paths <- file.path(media_dir, media_names)
  main_pngs <- expected_main[str_ends(expected_main, ".png")]
  main_md5 <- unname(tools::md5sum(main_pngs))
  media_md5 <- unname(tools::md5sum(media_paths))
  manuscript_media_audit <- tibble(
    embedded_media = basename(media_paths),
    embedded_md5 = media_md5,
    exact_main_match = media_md5 %in% main_md5,
    matched_main_figure = vapply(media_md5, function(value) {
      hit <- basename(main_pngs[main_md5 == value])
      if (length(hit) == 1L) hit else ""
    }, character(1)),
    status = if_else(exact_main_match, "PASS_EXACT_MATCH", "NOT_A_CURRENT_MAIN_PNG")
  )
  matched_main_n <- n_distinct(manuscript_media_audit$matched_main_figure[
    manuscript_media_audit$exact_main_match
  ])
  if (enforce_manuscript_audit && matched_main_n < 5L) {
    stop(
      "Manuscript media audit failed: only ", matched_main_n,
      "/5 current main PNGs are embedded exactly. Update the DOCX or disable enforcement for figure-only runs."
    )
  }
} else {
  manuscript_media_audit <- tibble(
    embedded_media = NA_character_, embedded_md5 = NA_character_,
    exact_main_match = NA, matched_main_figure = NA_character_,
    status = if_else(nzchar(manuscript_docx), "MANUSCRIPT_PATH_ABSENT", "NOT_REQUESTED")
  )
}
manuscript_audit_status <- if (all(is.na(manuscript_media_audit$exact_main_match))) {
  manuscript_media_audit$status[[1]]
} else {
  matched_main_n <- n_distinct(manuscript_media_audit$matched_main_figure[
    manuscript_media_audit$exact_main_match %in% TRUE
  ])
  if (matched_main_n >= 5L) "PASS_ALL_5_MAIN_PNGS" else paste0("REVIEW_", matched_main_n, "_OF_5_EXACT")
}
readr::write_csv(
  manuscript_media_audit,
  file.path(AUDIT_DIR, "08B_MANUSCRIPT_MEDIA_IDENTITY_AUDIT_v8_10_0.csv"),
  na = ""
)

run_summary <- c(
  paste0("PART9_VERSION=", PART9_VERSION),
  paste0("COMPLETED_AT=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("TABLE_ROOT=", TABLE_ROOT),
  paste0("FIG_ROOT=", FIG_ROOT),
  paste0("TAB_ROOT=", TAB_ROOT),
  paste0("REPRESENTATIVE_SPATIAL_SECTION=", representative_section),
  paste0("MAIN_FIGURES=5"),
  paste0("MAIN_PANELS=", nrow(main_budget)),
  paste0("MAIN_PANEL_LAYOUT=3/4/5/6/4"),
  paste0("FINAL_MAIN_WIDTH_IN=", FINAL_MAIN_WIDTH_IN),
  paste0("SUPPLEMENTARY_FIGURES=6"),
  paste0("MAIN_LINE_INTERVAL_FOREST_PANELS=", sum(main_budget$line_interval_forest)),
  paste0("MAIN_LINE_INTERVAL_FOREST_IDS=", paste(main_budget$panel_id[main_budget$line_interval_forest], collapse = ";")),
  paste0("MAIN_MATRIX_LIKE_PANELS=", sum(main_budget$matrix_like)),
  paste0("MAIN_SIGNIFICANCE_ONLY_PANELS=", sum(main_budget$significance_only)),
  paste0("MAIN_EXOTIC_GEOMETRY_PANELS=", sum(main_budget$exotic_geometry)),
  paste0("SEMANTIC_PALETTE_ROLES=", nrow(palette_contract)),
  paste0("CLAIM_DRIVEN_RESTRUCTURE=paired_phenotype;spatial_anchor;cell_state;cross_section_replication;protein_plex_boundary"),
  paste0("PROMOTED_WITHOUT_DUPLICATION=S1a_stage_to_Fig1c;S3a_to_Fig3e;S4c_to_Fig4d;S4d_detection_to_Fig4b;S5b_to_Fig5d"),
  paste0("DEMOTED_DIAGNOSTIC=attenuation;survival_LOCO;technical_residual"),
  paste0("FIGURE4_LAYOUT=ATLAS_PLUS_SECTION_SUMMARY_OVER_ASSOCIATION_PLUS_OVERLAP_OVER_JOINT_MAP_PLUS_RIDGE"),
  paste0("FIGURE5_LAYOUT=FOREST_OVER_PROTEIN_DISTRIBUTION_PLEX_POINTS_AND_FULL_WIDTH_SENTINELS"),
  paste0("RENDER_QA_REPAIRS=", nrow(render_qa_repair)),
  paste0("METHOD_SEMANTIC_AUDIT_ITEMS=", nrow(method_semantic_audit)),
  paste0("SURVIVAL_PI_FAIL_CLOSED=TRUE"),
  paste0("ICI_EFFECT_SCALE=RANK_BISERIAL"),
  paste0("ICI_BOOTSTRAP_B=", ici_bootstrap_B),
  paste0("MAIN_FIGURE_SOURCE_TABLES=", sum(source_registry$branch == "MAIN")),
  paste0("SUPPLEMENTARY_FIGURE_SOURCE_TABLES=", sum(source_registry$branch == "SUPPLEMENTARY")),
  paste0("SUPPLEMENTARY_TABLES=", nrow(supp_table_manifest)),
  paste0("READ_GUESS_MAX=1000000"),
  paste0("READ_PARSE_POLICY=FAIL_CLOSED"),
  paste0("FIG2A_EXPLICIT_INPUT_COLUMNS=cohort;OLFML2B_z;CAF_Core;ECM_Remodeling;TGFb_Response;analysis_population"),
  paste0("MANUSCRIPT_MEDIA_AUDIT=", manuscript_audit_status),
  paste0("PACKAGE_ZIP=", PACKAGE_ZIP),
  "STATUS=COMPLETE"
)
writeLines(run_summary, file.path(AUDIT_DIR, "99_PART9_COMPLETE.txt"), useBytes = TRUE)

session_text <- capture.output(sessionInfo())
writeLines(session_text, file.path(AUDIT_DIR, "Part9_sessionInfo_v8_10_0.txt"), useBytes = TRUE)

copy_tree <- function(source_dir, destination_dir) {
  source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE,
                             include.dirs = FALSE, all.files = FALSE)
  if (length(source_files) < 1L) stop("Package source directory is empty: ", source_dir)
  relative_files <- substring(source_files, nchar(source_dir) + 2L)
  destination_files <- file.path(destination_dir, relative_files)
  invisible(lapply(unique(dirname(destination_files)), dir.create,
                   recursive = TRUE, showWarnings = FALSE))
  copied <- file.copy(source_files, destination_files, overwrite = TRUE, copy.mode = TRUE)
  if (length(copied) != length(source_files) || any(!copied)) {
    stop("Package copy failed for: ", source_dir)
  }
  invisible(destination_files)
}

if (dir.exists(PACKAGE_DIR)) unlink(PACKAGE_DIR, recursive = TRUE, force = TRUE)
dir.create(PACKAGE_DIR, recursive = TRUE, showWarnings = FALSE)
copy_tree(MAIN_FIG_DIR, file.path(PACKAGE_DIR, "Main_Figures"))
copy_tree(SUPP_FIG_DIR, file.path(PACKAGE_DIR, "Supplementary_Figures"))
copy_tree(MAIN_SOURCE_DIR, file.path(PACKAGE_DIR, "Main_Figure_Source_Tables"))
copy_tree(SUPP_SOURCE_DIR, file.path(PACKAGE_DIR, "Supplementary_Figure_Source_Tables"))
copy_tree(SUPP_TABLE_DIR, file.path(PACKAGE_DIR, "Supplementary_Tables"))
copy_tree(AUDIT_DIR, file.path(PACKAGE_DIR, "Audit"))

canonical_r_files <- file.path(
  ROOT, "R",
  c(
    "09_OLFML2B_PART9_PUBLICATION_OUTPUTS_v8_10_0.R",
    "RUN_PART9_PUBLICATION_v8_10_0.R"
  )
)
canonical_r_files <- canonical_r_files[file.exists(canonical_r_files)]
if (length(canonical_r_files) > 0L) {
  dir.create(file.path(PACKAGE_DIR, "R"), recursive = TRUE, showWarnings = FALSE)
  copied_r_files <- file.copy(
    canonical_r_files,
    file.path(PACKAGE_DIR, "R", basename(canonical_r_files)),
    overwrite = TRUE, copy.mode = TRUE
  )
  if (length(copied_r_files) != length(canonical_r_files) || any(!copied_r_files)) {
    stop("Failed to include the canonical Part9 scripts in the publication package.")
  }
}

package_files_before_manifest <- list.files(
  PACKAGE_DIR, recursive = TRUE, full.names = TRUE, include.dirs = FALSE
)
package_manifest <- tibble(
  relative_path = substring(package_files_before_manifest, nchar(PACKAGE_DIR) + 2L),
  size_bytes = file.info(package_files_before_manifest)$size,
  md5 = unname(tools::md5sum(package_files_before_manifest)),
  status = if_else(size_bytes > 0, "PASS", "FAIL")
)
if (any(package_manifest$status != "PASS")) stop("Package contains an empty file before ZIP creation.")
readr::write_csv(
  package_manifest,
  file.path(AUDIT_DIR, "09_PACKAGE_FILE_MANIFEST_v8_10_0.csv"),
  na = ""
)
readr::write_csv(
  package_manifest,
  file.path(PACKAGE_DIR, "Audit", "09_PACKAGE_FILE_MANIFEST_v8_10_0.csv"),
  na = ""
)

build_package_zip <- function() {
  if (file.exists(PACKAGE_ZIP)) unlink(PACKAGE_ZIP, force = TRUE)
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  setwd(PACKAGE_PARENT)
  utils::zip(
    zipfile = basename(PACKAGE_ZIP),
    files = basename(PACKAGE_DIR),
    flags = "-r9"
  )
  if (!file.exists(PACKAGE_ZIP) || file.info(PACKAGE_ZIP)$size <= 0) {
    stop("Publication ZIP creation failed: ", PACKAGE_ZIP)
  }
  invisible(PACKAGE_ZIP)
}

build_package_zip()
final_package_files <- list.files(
  PACKAGE_DIR, recursive = TRUE, full.names = TRUE, include.dirs = FALSE
)
expected_archive_paths <- file.path(
  OUT_ID,
  substring(final_package_files, nchar(PACKAGE_DIR) + 2L)
) %>% str_replace_all("\\\\", "/")
archive_listing <- utils::unzip(PACKAGE_ZIP, list = TRUE)$Name %>%
  str_replace_all("\\\\", "/")
archive_content_audit <- tibble(
  expected_path = expected_archive_paths,
  present_in_zip = expected_archive_paths %in% archive_listing,
  status = if_else(present_in_zip, "PASS", "FAIL")
)
if (any(archive_content_audit$status != "PASS")) {
  stop(
    "ZIP content audit failed before finalization: ",
    paste(archive_content_audit$expected_path[archive_content_audit$status != "PASS"], collapse = ", ")
  )
}
readr::write_csv(
  archive_content_audit,
  file.path(AUDIT_DIR, "10_ARCHIVE_CONTENT_AUDIT_v8_10_0.csv"),
  na = ""
)
readr::write_csv(
  archive_content_audit,
  file.path(PACKAGE_DIR, "Audit", "10_ARCHIVE_CONTENT_AUDIT_v8_10_0.csv"),
  na = ""
)

# Rebuild once so the archive contains its own content audit, then verify again.
build_package_zip()
final_package_files <- list.files(
  PACKAGE_DIR, recursive = TRUE, full.names = TRUE, include.dirs = FALSE
)
expected_archive_paths <- file.path(
  OUT_ID,
  substring(final_package_files, nchar(PACKAGE_DIR) + 2L)
) %>% str_replace_all("\\\\", "/")
archive_listing <- utils::unzip(PACKAGE_ZIP, list = TRUE)$Name %>%
  str_replace_all("\\\\", "/")
missing_final_archive_paths <- setdiff(expected_archive_paths, archive_listing)
if (length(missing_final_archive_paths) > 0L) {
  stop("Final ZIP is incomplete: ", paste(missing_final_archive_paths, collapse = ", "))
}

message("Part9 v8.10.1 complete: ", FIG_ROOT)
message("Verified publication package: ", PACKAGE_ZIP)
