# 07_biotic_indices.R
# Author: Rezvan Hatami
# Date: 10 March 2026
#
# Purpose:
# This script calculates biotic indices, relates diversity and selected
# taxa to PCO axes, and produces spatial plots and heat maps required
# for downstream interpretation of biological patterns.

rm(list = ls())

# ---- Setup: package management ------------------------------------------------
required_pkgs <- c(
  "dplyr", "readr", "tibble", "purrr", "tidyr", "broom",
  "vegan", "fs", "here"
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
figure_dir <- fs::path(project_dir, "figures", "07")

fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- Setup: read required inputs ----------------------------------------------
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
pwangbug_file <- fs::path(output_dir, "02_pwangbug_scores.rds")
diversity_file <- fs::path(output_dir, "01_diversity_raw.rds")

if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}
if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}
if (!fs::file_exists(wangbug_file)) {
  stop("Missing input file: ", wangbug_file)
}
if (!fs::file_exists(pwangbug_file)) {
  stop("Missing input file: ", pwangbug_file)
}
if (!fs::file_exists(diversity_file)) {
  stop("Missing input file: ", diversity_file)
}

script3_inputs <- readRDS(script3_inputs_file)
list2env(script3_inputs, envir = .GlobalEnv)

wangenv <- readRDS(wangenv_file)
wangbug <- readRDS(wangbug_file)
pwangbug <- readRDS(pwangbug_file)
diversity_tbl <- readRDS(diversity_file)

# ----------------------------------------------------------------------
# This chunk defines the common plotting key used across spatial panels.
# ----------------------------------------------------------------------
day_key <- tibble::tribble(
  ~day, ~col, ~pch, ~lty, ~label,
  1, "black", 19, 1, "Dec 2013",
  126, "gold2", 15, 2, "April 2014",
  260, "blue", 8, 3, "Aug 2014",
  336, "green3", 17, 4, "Nov 2014",
  518, "red", 18, 5, "May 2015"
)

# ----------------------------------------------------------------------
# This chunk returns the sampling days available in a plotting dataset.
# ----------------------------------------------------------------------
day_id <- function(data) {
  day_key$day[day_key$day %in% unique(data$day)]
}

# ----------------------------------------------------------------------
# This chunk rebuilds a plotting frame for a fitted model.
# ----------------------------------------------------------------------
mod_fr <- function(model, data) {
  rows <- rownames(stats::model.frame(model))
  out <- data[rows, , drop = FALSE]
  out$fitted <- as.numeric(stats::fitted(model))
  out$residual <- as.numeric(stats::residuals(model))
  out
}

# ----------------------------------------------------------------------
# This chunk draws a distance-by-day plot for one response variable.
# ----------------------------------------------------------------------
plot_day <- function(
    data, y_var, main, ylab, ylim = NULL, xlab = "Distance(km)",
    cex_points = 1.5, legend_cex = 0.8
) {
  days_present <- day_id(data)
  first_day <- days_present[1]
  first_row <- day_key[day_key$day == first_day, , drop = FALSE]
  first_data <- data[data$day == first_day, , drop = FALSE] |>
    dplyr::arrange(dist)
  
  y_vals <- suppressWarnings(as.numeric(data[[y_var]]))
  if (is.null(ylim)) {
    ylim <- range(y_vals, na.rm = TRUE)
  }
  
  plot(
    first_data$dist, first_data[[y_var]], type = "b", ylim = ylim,
    pch = first_row$pch, col = first_row$col, cex = cex_points,
    xlab = xlab, ylab = ylab, main = main, lwd = 3, lty = first_row$lty
  )
  
  for (current_day in days_present[-1]) {
    day_row <- day_key[day_key$day == current_day, , drop = FALSE]
    day_data <- data[data$day == current_day, , drop = FALSE] |>
      dplyr::arrange(dist)
    
    points(
      day_data$dist, day_data[[y_var]], type = "b", col = day_row$col,
      pch = day_row$pch, lwd = 3,
      cex = if (current_day == 518) cex_points + 0.5 else cex_points,
      lty = day_row$lty
    )
  }
  
  legend(
    "topright", inset = c(0, 0),
    legend = day_key$label[match(days_present, day_key$day)],
    lty = day_key$lty[match(days_present, day_key$day)],
    pch = day_key$pch[match(days_present, day_key$day)], lwd = 2,
    col = day_key$col[match(days_present, day_key$day)], ncol = 2,
    horiz = FALSE, cex = legend_cex, title = "months"
  )
  
  abline(v = 4, lty = 2)
}

# ----------------------------------------------------------------------
# This chunk fits a univariate linear model and extracts a compact
# summary table.
# ----------------------------------------------------------------------
fit_lm <- function(data, response, predictor) {
  fit <- stats::lm(stats::as.formula(paste(response, "~", predictor)), data = data)
  sm <- summary(fit)
  
  tibble::tibble(
    response = response,
    predictor = predictor,
    n = stats::nobs(fit),
    intercept = unname(stats::coef(fit)[1]),
    slope = unname(stats::coef(fit)[2]),
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    p_value = sm$coefficients[2, 4]
  )
}

# ----------------------------------------------------------------------
# This chunk builds the transformed community matrix used in heat maps.
# ----------------------------------------------------------------------
agg_tab <- function(comm, scores) {
  brks <- cut(scores, breaks = stats::quantile(scores, seq(0, 1, 0.1)))
  stats::aggregate(comm, by = list(brks), FUN = mean)
}

# ----------------------------------------------------------------------
# 1. Calculate Shannon diversity and save the combined table.
# ----------------------------------------------------------------------
shannon <- vegan::diversity(wangbug, index = "shannon", MARGIN = 1, base = exp(1))
shan_tbl <- cbind(wangenv, shannon = shannon)

readr::write_csv(as.data.frame(shan_tbl), fs::path(output_dir, "07_shannon.csv"))
saveRDS(as.data.frame(shan_tbl), fs::path(output_dir, "07_shannon.rds"))

tiff(
  file = fs::path(figure_dir, "07_shannon_spatial_pattern.tif"),
  width = 7, height = 5, units = "in", pointsize = 12,
  bg = "transparent", res = 800, compression = "lzw"
)

par(mar = c(5, 5, 4, 2), cex = 0.9)
plot_day(
  data = as.data.frame(shan_tbl), y_var = "shannon",
  main = "Shannon index plotted against spatial position",
  ylab = "Shannon index", ylim = c(0, 3.5), legend_cex = 0.6
)
dev.off()

# ----------------------------------------------------------------------
# 2. Relate diversity indices to the first six PCO axes.
# ----------------------------------------------------------------------
div_tbl <- cbind(as.data.frame(diversity_tbl), as.data.frame(pwangbug))

div_vars <- intersect(
  c("shannon", "signal", "taxarich", "abundance", "eptnum", "eptrich"),
  names(div_tbl)
)

pco_vars <- intersect(paste0("pco", 1:6), names(div_tbl))

div_sum <- purrr::map_dfr(
  pco_vars,
  function(pco_var) {
    purrr::map_dfr(div_vars, function(div_var) fit_lm(div_tbl, pco_var, div_var))
  }
)

readr::write_csv(div_sum, fs::path(output_dir, "07_diversity_pco_summary.csv"))

# ----------------------------------------------------------------------
# 3. Fit selected taxon relationships against the first four PCO axes.
# ----------------------------------------------------------------------
taxa_vars <- intersect(
  c("gripopterygidae", "oligochaeta", "simulidae", "chironominae"),
  colnames(wangbug)
)

taxa_tbl <- cbind(as.data.frame(wangenv), as.data.frame(wangbug), as.data.frame(pwangbug))

taxa_sum <- purrr::map_dfr(
  paste0("pco", 1:4),
  function(pco_var) {
    purrr::map_dfr(taxa_vars, function(taxon) fit_lm(taxa_tbl, pco_var, taxon))
  }
)

readr::write_csv(taxa_sum, fs::path(output_dir, "07_selected_taxa_pco_summary.csv"))

# ----------------------------------------------------------------------
# 4. Produce diversity-vs-PCO scatterplot figures for the first six axes.
# ----------------------------------------------------------------------
for (pco_var in pco_vars) {
  fig_file <- fs::path(figure_dir, paste0("07_", pco_var, "_diversity_relationships.tif"))
  
  tiff(
    file = fig_file, width = 12, height = 8, units = "in", pointsize = 12,
    bg = "transparent", res = 800, compression = "lzw"
  )
  
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2.5, 2), cex = 0.9, las = 1)
  
  for (div_var in div_vars) {
    fit <- stats::lm(stats::as.formula(paste(pco_var, "~", div_var)), data = div_tbl)
    
    plot(
      div_tbl[[div_var]], div_tbl[[pco_var]], pch = 19,
      xlab = div_var, ylab = pco_var,
      main = paste(pco_var, "against", div_var)
    )
    abline(fit, lwd = 2)
  }
  
  dev.off()
}

# ----------------------------------------------------------------------
# 5. Produce selected taxon-vs-PCO scatterplot figures.
# ----------------------------------------------------------------------
for (pco_var in paste0("pco", 1:4)) {
  fig_file <- fs::path(figure_dir, paste0("07_", pco_var, "_selected_taxa.tif"))
  
  tiff(
    file = fig_file, width = 10, height = 8, units = "in", pointsize = 12,
    bg = "transparent", res = 800, compression = "lzw"
  )
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 2), cex = 0.9, las = 1)
  
  for (taxon in taxa_vars) {
    fit <- stats::lm(stats::as.formula(paste(pco_var, "~", taxon)), data = taxa_tbl)
    
    plot(
      taxa_tbl[[taxon]], taxa_tbl[[pco_var]], pch = 19,
      xlab = taxon, ylab = pco_var,
      main = paste(pco_var, "against", taxon)
    )
    abline(fit, lwd = 2)
  }
  
  dev.off()
}

# ----------------------------------------------------------------------
# 6. Plot selected macroinvertebrate taxa against spatial position.
# ----------------------------------------------------------------------
spatial_taxa <- intersect(
  c(
    "baetidae", "caenidae", "coloburiscidae", "leptophlebiidae",
    "gripopterygidae", "hydroptilidae", "leptoceridae", "calocidae",
    "hydrobiosidae", "hydropsychidae", "ecnomidae", "ceratopogonidae",
    "chrysomelidae", "diamesinae", "empididae", "micronectidae",
    "orthocladiinae", "simulidae", "tanypodinae", "tipulidae",
    "oligochaeta", "chironominae", "psychodidae", "ancylidae",
    "anostraca", "hydrochidae", "hydrophilidae", "muscidae",
    "physidae", "podonominae", "scitidae", "staphylinidae",
    "culicidae", "curculionidae"
  ),
  colnames(wangbug)
)

bugenv_tbl <- cbind(as.data.frame(wangenv), as.data.frame(wangbug))

saveRDS(spatial_taxa, fs::path(output_dir, "07_spatial_taxa_used.rds"))

for (taxon in spatial_taxa) {
  ymax <- max(bugenv_tbl[[taxon]], na.rm = TRUE)
  fig_file <- fs::path(figure_dir, paste0("07_", taxon, "_spatial_pattern.tif"))
  
  tiff(
    file = fig_file, width = 7, height = 5, units = "in", pointsize = 12,
    bg = "transparent", res = 800, compression = "lzw"
  )
  
  par(mar = c(5, 5, 4, 2), cex = 0.9)
  plot_day(
    data = bugenv_tbl, y_var = taxon,
    main = paste(taxon, "plotted against spatial position"),
    ylab = taxon, ylim = c(0, ymax)
  )
  dev.off()
}

# ----------------------------------------------------------------------
# 7. Relate selected taxa to chromium as an exploratory check.
# ----------------------------------------------------------------------
if (all(c("oligochaeta", "cr") %in% names(bugenv_tbl))) {
  cor_cr <- stats::cor.test(bugenv_tbl$oligochaeta, bugenv_tbl$cr)
  
  cor_tbl <- tibble::tibble(
    response = "oligochaeta",
    predictor = "cr",
    estimate = unname(cor_cr$estimate),
    statistic = unname(cor_cr$statistic),
    p_value = cor_cr$p.value,
    method = cor_cr$method
  )
  
  readr::write_csv(cor_tbl, fs::path(output_dir, "07_oligochaeta_cr_correlation.csv"))
}




# ----------------------------------------------------------------------
# 8. Produce PCO heat maps using taxon aggregates across PCO classes.
# ----------------------------------------------------------------------
if (requireNamespace("labdsv", quietly = TRUE)) {
  max_bug <- apply(wangbug, 2, max, na.rm = TRUE)
  keep_taxa <- names(max_bug[max_bug >= 5])
  wangbug_heat <- wangbug[, keep_taxa, drop = FALSE]
  
  sptree <- stats::hclust(
    vegan::vegdist(t((wangbug_heat)^0.25), "raup"),
    method = "average"
  )
  
  pco_axes <- intersect(paste0("pco", 1:6), colnames(pwangbug))
  
  for (pco_var in pco_axes) {
    agg <- agg_tab(wangbug_heat, pwangbug[, pco_var])
    fig_file <- fs::path(figure_dir, paste0("07_heat_", pco_var, ".tif"))
    
    tiff(
      file = fig_file, width = 5, height = 5, units = "in", pointsize = 12,
      bg = "transparent", res = 1200, compression = "lzw"
    )
    
    labdsv:::tabasco(
      (as.data.frame(agg[, -1]))^0.25,
      use = as.numeric(rownames(agg)),
      sp.ind = sptree,
      labCol = agg$Group.1,
      main = toupper(pco_var),
      cex.main = 0.2,
      cexRow = 0.6,
      cexCol = 0.5
    )
    
    dev.off()
  }
}

# ----------------------------------------------------------------------
# 9. Save a compact manifest of key outputs from this script.
# ----------------------------------------------------------------------
manifest_tbl <- tibble::tribble(
  ~output_file, ~description,
  "07_shannon.csv", "Shannon index combined with environmental data",
  "07_shannon.rds", "Saved Shannon table for downstream use",
  "07_diversity_pco_summary.csv", "Linear model summaries for diversity indices against PCO1 to PCO6",
  "07_selected_taxa_pco_summary.csv", "Linear model summaries for selected taxa against PCO1 to PCO4",
  "07_spatial_taxa_used.rds", "Vector of taxa plotted against spatial position"
)

readr::write_csv(manifest_tbl, fs::path(output_dir, "07_output_manifest.csv"))

cat("Script 07 completed successfully.\n")