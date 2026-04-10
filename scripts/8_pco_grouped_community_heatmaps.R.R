# 08_pco_grouped_community_heatmaps.R
# Author: Rezvan Hatami
# Date: 11-03-2026
#
# Purpose:
# This script reproduces the community heatmap workflow separately from
# Script 07 and focuses only on grouped community heatmaps along PCO
# gradients.
# It rebuilds the species clustering tree, groups samples by PCO deciles,
# aggregates community composition within intervals, and writes the
# heatmaps and supporting outputs needed for interpretation.

rm(list = ls())

# ---- Setup: package management ------------------------------------------------
required_pkgs <- c(
  "dplyr", "readr", "tibble", "fs", "here", "vegan", "BiodiversityR"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "These packages must be installed before running the script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
set.seed(123)

# ---- Setup: path configuration ------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "08")

fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- Setup: read required inputs ----------------------------------------------
wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
pwangbug_file <- fs::path(output_dir, "02_pwangbug_scores.rds")

if (!fs::file_exists(wangbug_file)) {
  stop("Missing input file: ", wangbug_file)
}
if (!fs::file_exists(pwangbug_file)) {
  stop("Missing input file: ", pwangbug_file)
}

wangbug <- readRDS(wangbug_file)
pwangbug <- readRDS(pwangbug_file)

# ----------------------------------------------------------------------
# This chunk converts compatible columns to numeric so the community and
# PCO objects are ready for aggregation and heatmap plotting.
# ----------------------------------------------------------------------
prep_num <- function(df) {
  for (nm in names(df)) {
    suppressWarnings({
      converted <- as.numeric(df[[nm]])
    })
    if (!all(is.na(converted))) {
      df[[nm]] <- converted
    }
  }
  df
}

# ----------------------------------------------------------------------
# This chunk retrieves the tabasco function from BiodiversityR using the
# installed namespace available on the current system.
# ----------------------------------------------------------------------
get_tab <- function() {
  if ("tabasco" %in% getNamespaceExports("BiodiversityR")) {
    return(getExportedValue("BiodiversityR", "tabasco"))
  }
  get("tabasco", envir = asNamespace("BiodiversityR"))
}

# ----------------------------------------------------------------------
# This chunk creates the original decile break sequence used to group
# samples along each PCO axis.
# ----------------------------------------------------------------------
mk_breaks <- function(x) {
  as.numeric(stats::quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE, type = 7))
}

# ----------------------------------------------------------------------
# This chunk assigns each sample to its original PCO interval.
# ----------------------------------------------------------------------
mk_cut <- function(x) {
  breaks <- unique(mk_breaks(x))
  
  if (length(breaks) < 2) {
    return(factor(rep("all", length(x))))
  }
  
  cut(x, breaks = breaks)
}

# ----------------------------------------------------------------------
# This chunk aggregates community composition within grouped PCO
# intervals using the original averaging logic.
# ----------------------------------------------------------------------
agg_comm <- function(comm_df, breaks_factor) {
  keep <- which(!is.na(breaks_factor))
  
  if (length(keep) == 0) {
    return(NULL)
  }
  
  stats::aggregate(
    as.data.frame(comm_df[keep, , drop = FALSE]),
    by = list(breaks_factor[keep]),
    FUN = mean
  )
}

# ----------------------------------------------------------------------
# This chunk writes one tabasco heatmap using the original plotting
# logic and catches plotting errors cleanly.
# ----------------------------------------------------------------------
plot_tab <- function(comm_df, use_vals, sp_tree, lab_col, main_title, file_path) {
  tab_fun <- get_tab()
  
  tiff(
    file = file_path, width = 5, height = 5, units = "in", pointsize = 12,
    bg = "transparent", res = 1200, compression = "lzw"
  )
  
  tryCatch(
    {
      tab_fun(
        comm_df,
        use = use_vals,
        sp.ind = sp_tree,
        labCol = lab_col,
        main = main_title,
        cex.main = 0.2,
        cexRow = 0.6,
        cexCol = 0.5
      )
    },
    error = function(e) {
      plot.new()
      title(main_title)
      text(0.5, 0.5, labels = conditionMessage(e))
    }
  )
  
  dev.off()
}

# ----------------------------------------------------------------------
# 1. Prepare analysis objects and align row counts.
# ----------------------------------------------------------------------
wangbug <- prep_num(as.data.frame(wangbug))
pwangbug <- prep_num(as.data.frame(pwangbug))

shared_n <- min(nrow(wangbug), nrow(pwangbug))

if (shared_n == 0) {
  stop("At least one required input object is empty.")
}

if (nrow(wangbug) != nrow(pwangbug)) {
  warning(
    "wangbug and pwangbug have different row counts. The script will ",
    "align them by row order using the smallest common size."
  )
}

wangbug <- wangbug[seq_len(shared_n), , drop = FALSE]
pwangbug <- pwangbug[seq_len(shared_n), , drop = FALSE]

pco_axes <- paste0("pco", 1:6)
pco_axes <- pco_axes[pco_axes %in% names(pwangbug)]

if (length(pco_axes) == 0) {
  stop("No PCO axis columns were found in pwangbug.")
}

# ----------------------------------------------------------------------
# 2. Build the original species clustering tree after removing taxa with
#    maximum abundance below five.
# ----------------------------------------------------------------------
max_bug <- apply(wangbug, 2, max, na.rm = TRUE)
rare_taxa <- names(which(max_bug < 5))

if (length(rare_taxa) < ncol(wangbug)) {
  wangbug1 <- wangbug[, !names(wangbug) %in% rare_taxa, drop = FALSE]
} else {
  wangbug1 <- wangbug
}

sptree <- stats::hclust(
  vegan::vegdist(t((as.matrix(wangbug1))^0.25), method = "raup"),
  method = "average"
)

readr::write_csv(
  tibble::tibble(taxon = names(wangbug1)[sptree$order]),
  fs::path(output_dir, "08_taxon_order_from_sptree.csv")
)

# ----------------------------------------------------------------------
# 3. Reproduce the grouped community aggregation for PCO1 to PCO6 and
#    save the grouped community tables.
# ----------------------------------------------------------------------
interval_tbl <- tibble::tibble(sample_id = seq_len(nrow(pwangbug)))
agg_list <- list()

for (axis_name in pco_axes) {
  breaks_factor <- mk_cut(pwangbug[[axis_name]])
  interval_tbl[[paste0(axis_name, "_interval")]] <- as.character(breaks_factor)
  
  agg_df <- agg_comm(wangbug, breaks_factor)
  agg_list[[axis_name]] <- agg_df
  
  if (!is.null(agg_df)) {
    readr::write_csv(
      tibble::tibble(interval = agg_df$Group.1),
      fs::path(output_dir, paste0("08_", axis_name, "_interval_labels.csv"))
    )
    
    readr::write_csv(
      tibble::as_tibble(agg_df),
      fs::path(output_dir, paste0("08_", axis_name, "_aggregated_community.csv"))
    )
  }
}

readr::write_csv(interval_tbl, fs::path(output_dir, "08_pco_interval_assignments.csv"))
saveRDS(agg_list, fs::path(output_dir, "08_aggregated_community_tables.rds"))

# ----------------------------------------------------------------------
# 4. Reproduce the original ungrouped PCO1 heatmap.
# ----------------------------------------------------------------------
fig_manifest <- tibble::tibble(
  figure_file = character(),
  axis = character(),
  description = character()
)

if ("pco1" %in% names(pwangbug)) {
  fig_file <- fs::path(figure_dir, "08_pco1_ungrouped_heatmap.tif")
  
  plot_tab(
    comm_df = (as.matrix(wangbug))^0.25,
    use_vals = pwangbug[, "pco1"],
    sp_tree = sptree,
    lab_col = NULL,
    main_title = "PCO1",
    file_path = fig_file
  )
  
  fig_manifest <- dplyr::bind_rows(
    fig_manifest,
    tibble::tibble(
      figure_file = fs::path_file(fig_file),
      axis = "pco1",
      description = "Ungrouped heatmap ordered by PCO1 scores"
    )
  )
}

# ----------------------------------------------------------------------
# 5. Reproduce grouped heatmaps for PCO1 to PCO6 using the original
#    grouped tabasco workflow.
# ----------------------------------------------------------------------
for (axis_name in names(agg_list)) {
  agg_df <- agg_list[[axis_name]]
  
  if (is.null(agg_df) || nrow(agg_df) == 0 || ncol(agg_df) <= 1) {
    next
  }
  
  comm_mat <- as.data.frame(agg_df[, -1, drop = FALSE])
  int_labels <- agg_df$Group.1
  use_vals <- suppressWarnings(as.numeric(rownames(agg_df)))
  
  if (all(is.na(use_vals))) {
    use_vals <- seq_len(nrow(agg_df))
  }
  
  fig_file <- fs::path(figure_dir, paste0("08_", axis_name, "_grouped_heatmap.tif"))
  
  plot_tab(
    comm_df = (as.matrix(comm_mat))^0.25,
    use_vals = use_vals,
    sp_tree = sptree,
    lab_col = int_labels,
    main_title = toupper(axis_name),
    file_path = fig_file
  )
  
  fig_manifest <- dplyr::bind_rows(
    fig_manifest,
    tibble::tibble(
      figure_file = fs::path_file(fig_file),
      axis = axis_name,
      description = "Grouped heatmap across PCO score intervals"
    )
  )
}

readr::write_csv(fig_manifest, fs::path(output_dir, "08_figure_manifest.csv"))

cat("Script 08 completed successfully.\n")