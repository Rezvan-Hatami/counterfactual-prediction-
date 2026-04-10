# 05_environmental_dbrda_models.R
# Author: Rezvan Hatami
# Date: 09 March 2026
#
# ===================================
# Objectives:
# Project root: Causality_Wang
# This script continues the pipeline after Script 04 and performs the
# environmental modelling stage of the analysis. It fits candidate
# environmental models for the main ordination axes, evaluates the
# environmental explanation of the community pattern, and tests whether
# spatiotemporal structure remains after conditioning on selected
# environmental predictors.
# ===================================
#
# ======================================================================
# Workflow in this script:
# 11. Fit candidate environmental models for individual PCO axes.
# 12. Keep the PCO1 environmental modelling block together with Figure 7,
#     so the main comparison between observed PCO1, the spatiotemporal
#     model prediction, and the environmental prediction remains in one place.
# 13. Fit dbRDA models using environmental predictors to test which
#     measured vaiables explain the community pattern.
# 14. Test whether any spatiotemporal structure remains after conditioning
#     on the selected environmental predictors.
# ==
# ======================================================================

rm(list = ls())

# ---- Setup: package management -------------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, vegan, permute, broom, fs, here, cli, MASS, lattice, car)

# ---- Setup: path configuration -------------------------------------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "05")
fs::dir_create(c(output_dir, figure_dir))

# ---- Setup: read saved analysis state from Script 03 ---------------------------------------------------------
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")
if (!fs::file_exists(script3_inputs_file)) cli::cli_abort("Missing input file: {.path {script3_inputs_file}}")

script3_inputs <- readr::read_rds(script3_inputs_file)
list2env(script3_inputs, envir = .GlobalEnv)

# ---- Setup: read auxiliary inputs used later in this script --------------------------------------------------
outliers_file <- fs::path(project_dir, "data", "wangenvOutliers.csv")
if (!fs::file_exists(outliers_file)) cli::cli_abort("Missing input file: {.path {outliers_file}}")

wangenvOutliers <- readr::read_csv(outliers_file, show_col_types = FALSE)
bugenv1 <- cbind(wangenvOutliers, pwangbug)

# ---- Setup: confirm that the Script 05 modelling environment is ready ----------------------------------------
cli::cli_alert_success("Script 05 modelling environment is ready.")

# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# 11. Rebuild the main space-time dbRDA model and prepare the complete-case
#     environmental analysis subset.
# ----------------------------------------------------------------------

# Refit the main space-time dbRDA model as the reference model for subsequent environmental analyses
wang.cap1 <- vegan::capscale(
  formula = wangbug.BC ~ dist + time + dist:time + eff + eff:time + eff:dist:time,
  data = wangenv,
  comm = wangbug,
  add = TRUE,
  na.action = na.omit
)

summary(wang.cap1)
wang.anova1 <- anova(wang.cap1, by = "term")
print(wang.anova1)

# Retain only complete cases across the environmental variables required for downstream analyses
env_vars_stage11 <- c("vel", "alk", "sed", "ph2", "turb", "toc", "cfpom", "temp", "no2", "nh3", "chla", "cond", "on")

wangenv1 <- wangenv |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(env_vars_stage11), ~ !is.na(.)))

# Match community and environmental datasets after environmental filtering
shared_sites <- intersect(rownames(wangbug), rownames(wangenv1))

wangenv1 <- wangenv1[shared_sites, , drop = FALSE]
wangbug1 <- wangbug[shared_sites, , drop = FALSE]

# Confirm that the filtered datasets are aligned correctly for analysis
stopifnot(nrow(wangbug1) == nrow(wangenv1))
stopifnot(identical(rownames(wangbug1), rownames(wangenv1)))


# ----------------------------------------------------------------------
# 12. Fit candidate environmental models for individual PCO axes.
#     Begin with PCO1.
# ----------------------------------------------------------------------

# Read the outlier-adjusted environmental table used for the COD check in the original workflow
outliers_file <- fs::path(project_dir, "data", "wangenvOutliers.csv")
if (!fs::file_exists(outliers_file)) cli::cli_abort("Missing input file: {.path {outliers_file}}")

wangenvOutliers <- readr::read_csv(outliers_file, show_col_types = FALSE)
bugenv1 <- cbind(wangenvOutliers, pwangbug)

# Define the candidate PCO1 models exactly as specified in the original script
pco1_model_specs <- tibble::tribble(
  ~model_id, ~data_name, ~formula,
  "pco1_temp_1",  "alldata", "pco1 ~ temp",
  "pco1_cond_1",  "alldata", "pco1 ~ cond",
  "pco1_temp_2",  "alldata", "pco1 ~ temp + cond + temp:cond",
  "pco1_cr_1",    "alldata", "pco1 ~ log(cr)",
  "pco1_cr_2",    "alldata", "pco1 ~ temp + cr + temp:cr",
  "pco1_cr_3",    "alldata", "pco1 ~ temp + cond + temp:cond + cr + temp:cr",
  "pco1_zn_1",    "alldata", "pco1 ~ log(zn)",
  "pco1_zn_2",    "alldata", "pco1 ~ temp + log(zn) + log(zn):temp",
  "pco1_zn_3",    "alldata", "pco1 ~ temp + cond + temp:cond + log(zn) + log(zn):temp",
  "pco1_cd_1",    "alldata", "pco1 ~ log(cd)",
  "pco1_cd_2",    "alldata", "pco1 ~ temp + log(cd) + log(cd):temp",
  "pco1_cd_3",    "alldata", "pco1 ~ temp + cond + temp:cond + cd + cd:temp",
  "pco1_cu_1",    "alldata", "pco1 ~ cu",
  "pco1_cu_2",    "alldata", "pco1 ~ temp + cu + cu:temp",
  "pco1_cu_3",    "alldata", "pco1 ~ temp + cond + temp:cond + cu + cu:temp",
  "pco1_chla_1",  "alldata", "pco1 ~ log(chla)",
  "pco1_chla_2",  "alldata", "pco1 ~ temp + log(chla) + log(chla):temp",
  "pco1_chla_3",  "alldata", "pco1 ~ temp + cond + temp:cond + log(chla) + temp:log(chla)",
  "pco1_do_1",    "alldata", "pco1 ~ do",
  "pco1_do_2",    "alldata", "pco1 ~ temp + do + do:temp",
  "pco1_do_3",    "alldata", "pco1 ~ temp + cond + temp:cond + do + do:temp",
  "pco1_cfpom_1", "alldata", "pco1 ~ cfpom",
  "pco1_cfpom_2", "alldata", "pco1 ~ cfpom + temp + temp:cfpom",
  "pco1_cfpom_3", "alldata", "pco1 ~ temp + cond + temp:cond + cfpom + cfpom:temp",
  "pco1_cod_1",   "bugenv1", "pco1 ~ cod",
  "pco1_cod_2",   "alldata", "pco1 ~ cod + temp + temp:cod",
  "pco1_cod_3",   "alldata", "pco1 ~ temp + cond + temp:cond + cod + cod:temp",
  "pco1_turb_1",  "alldata", "pco1 ~ turb",
  "pco1_turb_2",  "alldata", "pco1 ~ temp + turb + turb:temp",
  "pco1_turb_3",  "alldata", "pco1 ~ temp + cond + temp:cond + turb + turb:temp",
  "pco1_vel_1",   "alldata", "pco1 ~ vel",
  "pco1_vel_2",   "alldata", "pco1 ~ temp + vel + vel:temp",
  "pco1_vel_3",   "alldata", "pco1 ~ temp + cond + temp:cond + vel + vel:temp",
  "pco1_nh3_1",   "alldata", "pco1 ~ nh3",
  "pco1_nh3_2",   "alldata", "pco1 ~ temp + nh3 + nh3:temp",
  "pco1_nh3_3",   "alldata", "pco1 ~ temp + cond + temp:cond + nh3 + nh3:temp",
  "pco1_toc_1",   "alldata", "pco1 ~ toc",
  "pco1_toc_2",   "alldata", "pco1 ~ temp + toc + toc:temp",
  "pco1_toc_3",   "alldata", "pco1 ~ temp + cond + temp:cond + toc + toc:temp",
  "pco1_ph2_1",   "alldata", "pco1 ~ ph2",
  "pco1_ph2_2",   "alldata", "pco1 ~ temp + ph2 + ph2:temp",
  "pco1_ph2_3",   "alldata", "pco1 ~ temp + cond + temp:cond + ph2 + ph2:temp",
  "pco1_sed_1",   "alldata", "pco1 ~ sed",
  "pco1_sed_2",   "alldata", "pco1 ~ temp + sed + sed:temp",
  "pco1_on_1",    "alldata", "pco1 ~ on",
  "pco1_on_2",    "alldata", "pco1 ~ temp + on + on:temp"
)

# Fit all candidate models in a reproducible table-driven workflow
data_lookup <- list(alldata = alldata, bugenv1 = bugenv1)

pco1_candidate_models <- pco1_model_specs |>
  dplyr::mutate(
    model = purrr::map2(
      formula,
      data_name,
      ~ stats::lm(stats::as.formula(.x), data = data_lookup[[.y]])
    ),
    glance = purrr::map(model, broom::glance),
    tidy = purrr::map(model, broom::tidy)
  )

# Create a compact model comparison table for rapid screening
pco1_model_summary <- pco1_candidate_models |>
  dplyr::transmute(
    model_id = model_id,
    data_name = data_name,
    formula = formula,
    glance = glance
  ) |>
  tidyr::unnest(glance) |>
  dplyr::select(model_id, data_name, formula, r.squared, adj.r.squared, sigma, AIC, BIC, p.value) |>
  dplyr::arrange(dplyr::desc(adj.r.squared))

print(pco1_model_summary)

# Note: This model-fitting function defines the standard template used to evaluate each pre-specified candidate model in this step.
# Fit all candidate models in a reproducible table-driven workflow
data_lookup <- list(alldata = alldata, bugenv1 = bugenv1)

fit_candidate_model <- function(formula_text, data_name, data_lookup) {
  stats::lm(stats::as.formula(formula_text), data = data_lookup[[data_name]])
}

pco1_candidate_models <- pco1_model_specs |>
  dplyr::mutate(
    model = purrr::map2(formula, data_name, ~ fit_candidate_model(.x, .y, data_lookup)),
    glance = purrr::map(model, broom::glance)
  )

# Summarise all candidate models in one comparison table
pco1_model_summary <- pco1_candidate_models |>
  tidyr::unnest(glance) |>
  dplyr::select(model_id, data_name, formula, r.squared, adj.r.squared, sigma, AIC, BIC, p.value) |>
  dplyr::arrange(dplyr::desc(adj.r.squared))

print(pco1_model_summary)
# ----------------------------------------------------------------------
# 13. 13. The PCO1 environmental modelling block is presented together with Figure 7 
# because the figure requires the fitted environmental models and their predictions to be constructed.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# This chunk combines the environmental data and ordination scores so
# the observed PCO1 pattern can be plotted directly against distance and
# sampling occasion.
# ----------------------------------------------------------------------
bugenv <- cbind(wangenv, pwangbug)

# ----------------------------------------------------------------------
# This chunk defines a single plotting key so all figures use the same
# colours, symbols, line types, and labels for the five sampling times.
# ----------------------------------------------------------------------
day_key <- tibble::tribble(
  ~day, ~time, ~col, ~pch, ~lty, ~label,
  1, 1, "black", 19, 1, "Dec 13",
  126, 2, "gold2", 15, 2, "April 14",
  260, 3, "blue", 8, 3, "Aug 2014",
  336, 4, "green3", 17, 4, "Nov 2014",
  518, 5, "red", 18, 5, "May 2015"
)

# ----------------------------------------------------------------------
# This chunk fits a table of linear models from a specification table.
# ----------------------------------------------------------------------
fit_tbl <- function(specs, lookup) {
  specs |>
    dplyr::mutate(
      model = purrr::map2(
        formula, data_name,
        ~ stats::lm(stats::as.formula(.x), data = lookup[[.y]])
      ),
      glance = purrr::map(model, broom::glance)
    )
}

# ----------------------------------------------------------------------
# This chunk rebuilds the model frame and appends fitted values and
# residuals so plotting data stay aligned with the fitted model.
# ----------------------------------------------------------------------
mod_fr <- function(model, data) {
  rows <- rownames(stats::model.frame(model))
  out <- data[rows, , drop = FALSE]
  out$fitted <- stats::fitted(model)
  out$residual <- stats::residuals(model)
  out
}

# ----------------------------------------------------------------------
# This chunk returns the sampling days available in a plotting dataset.
# ----------------------------------------------------------------------
day_id <- function(data) {
  day_key$day[day_key$day %in% unique(data$day)]
}

# ----------------------------------------------------------------------
# This chunk draws the observed PCO1 pattern for each sampling day.
# ----------------------------------------------------------------------
plot_obs <- function(
    data, ylim, main, xlab = "", ylab = "PCO1", cex_points = 1.2
) {
  days_present <- day_id(data)
  first_day <- days_present[1]
  first_row <- day_key[day_key$day == first_day, , drop = FALSE]
  first_data <- data[data$day == first_day, , drop = FALSE] |>
    dplyr::arrange(dist)
  
  plot(
    first_data$dist, first_data$pco1, type = "b", ylim = ylim,
    pch = first_row$pch, col = first_row$col, cex = cex_points,
    xlab = xlab, ylab = ylab, main = main, lwd = 3, lty = first_row$lty
  )
  
  for (current_day in days_present[-1]) {
    day_row <- day_key[day_key$day == current_day, , drop = FALSE]
    day_data <- data[data$day == current_day, , drop = FALSE] |>
      dplyr::arrange(dist)
    
    points(
      day_data$dist, day_data$pco1, type = "b", col = day_row$col,
      pch = day_row$pch, cex = if (current_day == 518) 1.4 else cex_points,
      lwd = 3, lty = day_row$lty
    )
  }
  
  legend(
    "topright", inset = c(0, 0),
    legend = day_key$label[match(days_present, day_key$day)],
    lty = day_key$lty[match(days_present, day_key$day)],
    pch = day_key$pch[match(days_present, day_key$day)], lwd = 2,
    col = day_key$col[match(days_present, day_key$day)], ncol = 2,
    horiz = FALSE, cex = 0.8, title = "months"
  )
  
  abline(v = 4, lty = 2)
}

# ----------------------------------------------------------------------
# This chunk overlays observed values and modelled lines for each day.
# ----------------------------------------------------------------------
plot_fit <- function(
    point_data, point_var, line_data, line_var, ylim, main, xlab = "",
    ylab = "PCO1", point_cex = 1.2, legend_cex = 0.6
) {
  days_present <- intersect(day_id(point_data), day_id(line_data))
  first_day <- days_present[1]
  first_row <- day_key[day_key$day == first_day, , drop = FALSE]
  first_points <- point_data[point_data$day == first_day, , drop = FALSE] |>
    dplyr::arrange(dist)
  first_lines <- line_data[line_data$day == first_day, , drop = FALSE] |>
    dplyr::arrange(dist)
  
  plot(
    first_points$dist, first_points[[point_var]], type = "p",
    pch = first_row$pch, col = first_row$col, cex = point_cex,
    xlab = xlab, ylab = ylab, ylim = ylim, main = main
  )
  
  lines(
    first_lines$dist, first_lines[[line_var]], col = first_row$col,
    lty = first_row$lty, lwd = 3
  )
  
  for (current_day in days_present[-1]) {
    day_row <- day_key[day_key$day == current_day, , drop = FALSE]
    current_points <- point_data[
      point_data$day == current_day, , drop = FALSE
    ] |>
      dplyr::arrange(dist)
    current_lines <- line_data[
      line_data$day == current_day, , drop = FALSE
    ] |>
      dplyr::arrange(dist)
    
    points(
      current_points$dist, current_points[[point_var]], col = day_row$col,
      pch = day_row$pch, cex = if (current_day == 518) 1.4 else point_cex,
      lwd = 3
    )
    
    lines(
      current_lines$dist, current_lines[[line_var]], col = day_row$col,
      lty = day_row$lty, lwd = 3
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
# This chunk draws a day-series panel from one response variable.
# ----------------------------------------------------------------------
plot_day <- function(
    data, y_var, main, ylim, xlab = "Distance(km)", ylab = "PCO1",
    cex_points = 1.2
) {
  days_present <- day_id(data)
  first_day <- days_present[1]
  first_row <- day_key[day_key$day == first_day, , drop = FALSE]
  first_data <- data[data$day == first_day, , drop = FALSE] |>
    dplyr::arrange(dist)
  
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
      pch = day_row$pch, cex = if (current_day == 518) 1.4 else cex_points,
      lwd = 3, lty = day_row$lty
    )
  }
  
  legend(
    "topright", inset = c(0, 0),
    legend = day_key$label[match(days_present, day_key$day)],
    lty = day_key$lty[match(days_present, day_key$day)],
    pch = day_key$pch[match(days_present, day_key$day)], lwd = 2,
    col = day_key$col[match(days_present, day_key$day)], ncol = 2,
    horiz = FALSE, cex = 0.8, title = "months"
  )
  
  abline(v = 4, lty = 2)
}

# ----------------------------------------------------------------------
# This chunk fits the candidate PCO1 environmental models in a single
# table so the model set can be compared transparently.
# ----------------------------------------------------------------------
pco1_final_model_specs <- tibble::tribble(
  ~model_id, ~data_name, ~formula,
  "tempo1", "alldata", "pco1 ~ log(cd)",
  "tempo2", "alldata", "pco1 ~ log(chla)",
  "tempo3", "alldata", "pco1 ~ log(zn)",
  "tempo4", "alldata", "pco1 ~ temp",
  "tempo5", "alldata", "pco1 ~ temp + cod + temp:cod",
  "tempo6", "alldata", "pco1 ~ temp + cond + temp:cond",
  "tempo7", "alldata", "pco1 ~ temp + toc + temp:toc",
  "tempo8", "alldata", "pco1 ~ temp + log(zn) + temp:log(zn)",
  "tempo9", "alldata",
  "pco1 ~ temp + cond + temp:cond + cod + temp:cod",
  "tempo10", "alldata",
  "pco1 ~ temp + cod + temp:cod + log(zn) + temp:log(zn)",
  "tempo11", "alldata",
  "pco1 ~ temp + cod + temp:cod + log(zn) + temp:log(zn) + log(chla)",
  "tempo12", "alldata",
  "pco1 ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + temp:log(zn)",
  "tempo13", "alldata", "pco1 ~ temp + cod + temp:cod + log(zn)",
  "tempo14", "alldata",
  "pco1 ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + temp:log(zn) + log(chla) + log(cd)",
  "tempo15", "alldata",
  "pco1 ~ temp + cond + toc + temp:cond + log(zn) + temp:toc + log(chla) + log(cr)",
  "tempo16", "alldata", "pco1 ~ temp + cond + log(zn)"
)

data_lookup <- list(alldata = alldata)

pco1_final_models <- fit_tbl(pco1_final_model_specs, data_lookup)

pco1_final_model_summary <- pco1_final_models |>
  tidyr::unnest(glance) |>
  dplyr::select(
    model_id, data_name, formula, r.squared, adj.r.squared,
    sigma, AIC, BIC, p.value
  ) |>
  dplyr::arrange(dplyr::desc(adj.r.squared))

print(pco1_final_model_summary)

# ----------------------------------------------------------------------
# This chunk extracts the final preferred model and prepares fitted
# values that are later shown in Figure 7.
# ----------------------------------------------------------------------
tempo15 <- pco1_final_models$model[[which(
  pco1_final_models$model_id == "tempo15"
)]]
tempo16 <- pco1_final_models$model[[which(
  pco1_final_models$model_id == "tempo16"
)]]

anova(tempo15, tempo16)
summary(tempo15)
anova(tempo15)
car::avPlots(tempo15)

pco1_final_frame <- mod_fr(tempo15, alldata)
pco1_residual_against_pred1 <- stats::lm(
  residual ~ pred1, data = pco1_final_frame
)
summary(pco1_residual_against_pred1)

lattice::xyplot(
  fitted ~ dist, groups = day, data = pco1_final_frame,
  type = "l", auto.key = TRUE
)
lattice::xyplot(
  residual ~ dist, groups = day, data = pco1_final_frame,
  type = "l", auto.key = TRUE
)

# ----------------------------------------------------------------------
# This chunk keeps Figure 7 in a dedicated section because it is a
# standalone figure comparing the observed PCO1 pattern, the
# spatiotemporal predictions, and the final environmental model.
# ----------------------------------------------------------------------
tempo17 <- stats::lm(
  pco1 ~ temp + cond + toc + temp:cond + log(zn) + temp:toc +
    log(chla) + log(cr),
  data = bugenv
)

pco1_final_bugenv_frame <- mod_fr(tempo17, bugenv)

tiff(
  file = fs::path(figure_dir, "05_Figure_7.tif"),
  width = 10, height = 14, units = "in", pointsize = 12,
  bg = "transparent", res = 800, compression = "lzw"
)

par(
  mfrow = c(3, 1), mar = c(4.5, 4.5, 2.5, 3), cex = 1.5,
  cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8
)

# ----------------------------------------------------------------------
# This chunk produces panel A, which shows the observed PCO1 values
# against downstream distance before any model-based predictions.
# ----------------------------------------------------------------------
plot_obs(
  data = bugenv, ylim = c(-0.6, 1),
  main = "a) PCO1 plotted against spatial position",
  xlab = "", ylab = "PCO1", cex_points = 1.2
)

# ----------------------------------------------------------------------
# This chunk produces panel B, which overlays the observed PCO1 values
# with the spatiotemporal model predictions.
# ----------------------------------------------------------------------
plot_fit(
  point_data = bugenv, point_var = "pco1",
  line_data = alldata, line_var = "pred1", ylim = c(-0.6, 1),
  main = "b) PCO1 with spatiotemporal model predictions",
  xlab = "", ylab = "PCO1", point_cex = 1.2, legend_cex = 0.6
)

# ----------------------------------------------------------------------
# This chunk produces panel C, which overlays the observed PCO1 values
# with fitted values from the final environmental model.
# ----------------------------------------------------------------------
plot_fit(
  point_data = bugenv, point_var = "pco1",
  line_data = pco1_final_bugenv_frame, line_var = "fitted",
  ylim = c(-0.6, 1),
  main = paste0(
    "c) Prediction with zinc, chlorophyll A, temperature, TOC, conductivity,",
    "\ntemperature*TOC, and temperature*conductivity interaction"
  ),
  xlab = "Distance(km)", ylab = "PCO1", point_cex = 1.2,
  legend_cex = 0.6
)

dev.off()

# ----------------------------------------------------------------------
# This chunk defines the reduced display models used to show how the
# PCO1 prediction changes as temperature and conductivity are added.
# ----------------------------------------------------------------------
pco1_display_model_specs <- tibble::tribble(
  ~model_id, ~formula, ~title,
  "tempo18", "pco1 ~ temp", "a) Prediction with temperature",
  "tempo19", "pco1 ~ cond", "b) Prediction with conductivity",
  "tempo20", "pco1 ~ temp + cond",
  "c) Prediction with temperature and conductivity",
  "tempo21", "pco1 ~ temp + cond + temp:cond",
  "d) Prediction with temperature*conductivity interaction",
  "tempo22", "pco1 ~ tn + op + zn + temp + cond + temp:cond",
  paste0(
    "e) Prediction with temperature*conductivity interaction",
    "\n and other water quality parameters"
  )
)

pco1_display_models <- pco1_display_model_specs |>
  dplyr::mutate(
    model = purrr::map(
      formula, ~ stats::lm(stats::as.formula(.x), data = bugenv)
    ),
    frame = purrr::map(model, mod_fr, data = bugenv)
  )

tempo18 <- pco1_display_models$model[[which(
  pco1_display_models$model_id == "tempo18"
)]]
tempo19 <- pco1_display_models$model[[which(
  pco1_display_models$model_id == "tempo19"
)]]
tempo20 <- pco1_display_models$model[[which(
  pco1_display_models$model_id == "tempo20"
)]]
tempo21 <- pco1_display_models$model[[which(
  pco1_display_models$model_id == "tempo21"
)]]
tempo22 <- pco1_display_models$model[[which(
  pco1_display_models$model_id == "tempo22"
)]]

tiff(
  file = fs::path(
    figure_dir, "05_prediction_with_temperature_and_conductivity.tif"
  ),
  width = 14, height = 14, units = "in", pointsize = 12,
  bg = "transparent", res = 800, compression = "lzw"
)

par(
  mfrow = c(3, 2), mar = c(5, 5, 4, 2), cex = 0.9,
  cex.axis = 0.9, las = 1, cex.main = 1.2
)

for (i in seq_len(nrow(pco1_display_models))) {
  plot_day(
    data = pco1_display_models$frame[[i]], y_var = "fitted",
    main = pco1_display_models$title[[i]], ylim = c(-0.6, 0.6),
    xlab = "Distance(km)", ylab = "PCO1"
  )
}

plot_day(
  data = bugenv, y_var = "pco1",
  main = "f) PCO1 plotted against spatial position", ylim = c(-0.6, 0.68),
  xlab = "Distance(km)", ylab = "PCO1"
)

dev.off()

# ----------------------------------------------------------------------
# This chunk plots the residual pattern from the final environmental
# model across downstream distance and sampling day.
# ----------------------------------------------------------------------
tiff(
  file = fs::path(figure_dir, "05_pco1_environmental_model_residuals.tif"),
  width = 7, height = 5, units = "in", pointsize = 12,
  bg = "transparent", res = 800, compression = "lzw"
)

par(mar = c(5, 5, 4, 2), cex = 0.9, cex.axis = 0.9, las = 1)

plot_day(
  data = pco1_final_bugenv_frame, y_var = "residual", main = "Residuals",
  ylim = c(-0.6, 0.6), xlab = "Distance(km)", ylab = "PCO1"
)

dev.off()
# ----------------------------------------------------------------------
# 14. Fit dbRDA models using the selected environmental predictors.
# ----------------------------------------------------------------------

env_vars_stage14 <- c("cod", "temp", "zn", "chla", "cond")

wangenv1 <- wangenv |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(env_vars_stage14), ~ !is.na(.)))

shared_sites <- intersect(rownames(wangbug), rownames(wangenv1))

wangenv1 <- wangenv1[shared_sites, , drop = FALSE]
wangbug1 <- wangbug[shared_sites, , drop = FALSE]

stopifnot(nrow(wangbug1) == nrow(wangenv1))
stopifnot(identical(rownames(wangbug1), rownames(wangenv1)))

wangbug1.BC <- vegan::vegdist(sqrt(wangbug1))

dbrda_candidate_specs <- tibble::tribble(
  ~model_id, ~formula,
  "wang.cap2_1",  "wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + temp:log(zn) + toc + temp:toc + log(chla)",
  "wang.cap2_2",  "wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + toc + log(chla)",
  "wang.cap2_3",  "wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + log(chla)",
  "wang.cap2_4",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + temp:log(zn) + toc + temp:toc + log(chla)",
  "wang.cap2_5",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla)",
  "wang.cap2_6",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + nh3",
  "wang.cap2_7",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + no3",
  "wang.cap2_8",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + tp",
  "wang.cap2_9",  "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + tn",
  "wang.cap2_10", "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + sed",
  "wang.cap2_11", "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + ph2",
  "wang.cap2_12", "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + cfpom",
  "wang.cap2_13", "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + turb",
  "wang.cap2_14", "wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + log(cr)"
)

fit_dbrda_model <- function(formula_text, data, comm) {
  vegan::capscale(
    formula = stats::as.formula(formula_text),
    data = data,
    comm = comm,
    add = TRUE,
    na.action = na.omit
  )
}

wang.cap2_candidates <- dbrda_candidate_specs |>
  dplyr::mutate(
    model = purrr::map(formula, ~ fit_dbrda_model(.x, data = wangenv1, comm = wangbug1)),
    anova_overall = purrr::map(model, vegan::anova.cca),
    inertia = purrr::map_dbl(model, ~ .x$tot.chi),
    rank = purrr::map_int(model, ~ .x$CCA$qrank)
  )

wang.cap2_candidate_summary <- wang.cap2_candidates |>
  dplyr::transmute(
    model_id = model_id,
    formula = formula,
    inertia = inertia,
    rank = rank,
    p_value = purrr::map_dbl(anova_overall, ~ .x$`Pr(>F)`[1])
  ) |>
  dplyr::arrange(p_value)

print(wang.cap2_candidate_summary)

wang.cap2 <- wang.cap2_candidates$model[[which(wang.cap2_candidates$model_id == "wang.cap2_14")]]

summary(wang.cap2)
anova(wang.cap2)

wang.anova2 <- anova(wang.cap2, by = "term", permutations = permute::how(nperm = 9999))
print(wang.anova2)

# ----------------------------------------------------------------------
# 14. Fit dbRDA models using the selected environmental predictors.
# ----------------------------------------------------------------------

env_vars_stage14 <- c("cod", "temp", "zn", "chla", "cond", "toc", "cr")

wangenv1 <- wangenv |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(env_vars_stage14), ~ !is.na(.)))

shared_sites <- intersect(rownames(wangbug), rownames(wangenv1))

wangenv1 <- wangenv1[shared_sites, , drop = FALSE]
wangbug1 <- wangbug[shared_sites, , drop = FALSE]

stopifnot(nrow(wangbug1) == nrow(wangenv1))
stopifnot(identical(rownames(wangbug1), rownames(wangenv1)))

wangbug1.BC <- vegan::vegdist(sqrt(wangbug1))

wang.cap2.1 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + temp:log(zn) + toc + temp:toc + log(chla),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.1)
anova(wang.cap2.1)

wang.cap2.2 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + toc + log(chla),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.2)
anova(wang.cap2.2)

wang.cap2.3 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + cod + temp:cod + log(zn) + log(chla),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.3)
anova(wang.cap2.3)

wang.cap2.4 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + temp:log(zn) + toc + temp:toc + log(chla),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.4)
anova(wang.cap2.4)

wang.cap2.5 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.5)
anova(wang.cap2.5)

wang.cap2.6 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + nh3,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.6)
anova(wang.cap2.6)

wang.cap2.7 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + no3,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.7)
anova(wang.cap2.7)

wang.cap2.8 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + tp,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.8)
anova(wang.cap2.8)

wang.cap2.9 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + tn,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.9)
anova(wang.cap2.9)

wang.cap2.10 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + sed,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.10)
anova(wang.cap2.10)

wang.cap2.11 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + ph2,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.11)
anova(wang.cap2.11)

wang.cap2.12 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + cfpom,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.12)
anova(wang.cap2.12)

wang.cap2.13 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + turb,
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.13)
anova(wang.cap2.13)

wang.cap2.14 <- vegan::capscale(
  formula = wangbug1.BC ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla) + log(cr),
  data = wangenv1,
  comm = wangbug1,
  add = TRUE,
  na.action = na.omit
)
summary(wang.cap2.14)
anova(wang.cap2.14)

wang.cap2 <- wang.cap2.14

wang.anova2 <- anova(wang.cap2, by = "term", permutations = permute::how(nperm = 9999))
print(wang.anova2)

