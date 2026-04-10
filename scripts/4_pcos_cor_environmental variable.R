# 04_environmental_screening.R
# Author: Rezvan Hatami
# Date: 05 March 2026
#
# ===================================
# Objectives:
# Project root: Causality_Wang
# This script continues the pipeline after Script 03 and screens
# environmental variables before environmental dbRDA modelling.
# It uses exploratory visualisation and structured correlation summaries
# to identify candidate predictors for Script 05.
# ===================================
#
# ======================================================================
# Workflow in this script:
# 8. Screen environmental variables using grouped scatterplot matrices
#    and exploratory visual summaries.
# 9. Summarise variable-wise correlations with the main PCO axes.
# 10. Save the screening outputs required for Script 05.
# ======================================================================

rm(list = ls())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, GGally, broom, fs, here, cli)

project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "04")
fs::dir_create(c(output_dir, figure_dir))

input_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")
if (!fs::file_exists(input_file)) cli::cli_abort("Missing input file: {.path {input_file}}")

inputs <- readr::read_rds(input_file)
list2env(inputs, envir = .GlobalEnv)

wangalldata <- alldata

cli::cli_alert_success("Script 04 screening environment is ready.")

# ---- 8. Screen environmental variables against the main PCO axes using grouped scatterplot matrices ----------

panel_cor <- function(x, y, digits = 2, method = "spearman", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  keep <- complete.cases(x, y)
  
  if (sum(keep) < 3) {
    text(0.5, 0.5, "n < 3")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x[keep], y[keep], method = method))
  p <- suppressWarnings(cor.test(x[keep], y[keep], method = method)$p.value)
  
  r_txt <- paste0("r= ", formatC(r, digits = digits, format = "f"))
  p_txt <- if (is.na(p)) "p= NA" else if (p < 0.01) "p= <0.01" else paste0("p= ", formatC(p, digits = digits, format = "f"))
  
  text(0.5, 0.8, r_txt)
  text(0.5, 0.3, p_txt)
}

safe_log <- function(x) ifelse(is.na(x) | x <= 0, NA_real_, log(x))
safe_log1p <- function(x) ifelse(is.na(x) | x < 0, NA_real_, log1p(x))
safe_sqrt <- function(x) ifelse(is.na(x) | x < 0, NA_real_, sqrt(x))
safe_asin <- function(x) ifelse(is.na(x) | x < 0 | x > 1, NA_real_, asin(x))

alldata_pairs_log1p <- alldata |>
  dplyr::mutate(
    cd = safe_log1p(cd),
    cu = safe_log1p(cu),
    cr = safe_log1p(cr),
    zn = safe_log1p(zn),
    alk = safe_log1p(alk),
    temp = safe_log1p(temp),
    cond = safe_log1p(cond),
    do = safe_log1p(do),
    toc = safe_log1p(toc),
    cfpom = safe_log1p(cfpom),
    chla = safe_log1p(chla),
    turb = safe_log1p(turb),
    sed = safe_log1p(sed),
    vel = safe_log1p(vel),
    tp = safe_log1p(tp),
    nh3 = safe_log1p(nh3),
    no3 = safe_log1p(no3),
    no2 = safe_log1p(no2),
    on = safe_log1p(on),
    op = safe_log1p(op)
  )

alldata_pairs_pco_transformed <- alldata |>
  dplyr::mutate(
    cd = safe_log(cd),
    cu = safe_log(cu),
    cr = safe_log(cr),
    zn = safe_log(zn),
    cod = safe_sqrt(cod),
    cfpom = safe_log(cfpom),
    chla = safe_log(chla),
    ph2 = safe_log(ph2),
    turb = safe_log(turb),
    vel = safe_asin(vel)
  )

alldata_pairs_pred_transformed <- alldata |>
  dplyr::mutate(
    cd = safe_log(cd),
    cu = safe_log(cu),
    cr = safe_log(cr),
    zn = safe_log(zn),
    cfpom = safe_log(cfpom),
    chla = safe_log(chla),
    ph2 = safe_log(ph2),
    sed = safe_log(sed)
  )

save_pairs_plot <- function(data, formula_obj, filename, width = 14, height = 14, res = 800) {
  grDevices::tiff(
    file = fs::path(figure_dir, filename),
    width = width,
    height = height,
    units = "in",
    pointsize = 12,
    bg = "transparent",
    res = res,
    compression = "lzw"
  )
  par(mar = c(5, 5, 4, 2), cex = 0.9)
  pairs(formula_obj, data = data, upper.panel = panel_cor, pch = 20, na.action = na.omit)
  dev.off()
}

save_pairs_plot(
  data = alldata,
  formula_obj = ~ cd + cu + cr + zn + alk + temp + cond + do + toc + cfpom + chla + ph2 + turb + sed + vel + tp + nh3 + no3 + no2 + on + op + pco1 + pco2,
  filename = "04_pairs_pco1_pco2_raw_all.tif",
  width = 20,
  height = 14,
  res = 1200
)

save_pairs_plot(
  data = alldata_pairs_log1p,
  formula_obj = ~ cd + cu + cr + zn + alk + temp + cond + do + toc + cfpom + chla + ph2 + turb + sed + vel + tp + nh3 + no3 + no2 + on + op + pco1 + pco2,
  filename = "04_pairs_pco1_pco2_log1p_all.tif",
  width = 20,
  height = 14,
  res = 1200
)

save_pairs_plot(
  data = alldata,
  formula_obj = ~ ph1 + sali + do + pco1 + pco2,
  filename = "04_pairs_pco1_pco2_ph_sali_do_raw.tif"
)

save_pairs_plot(
  data = alldata_pairs_pco_transformed,
  formula_obj = ~ cd + cu + cr + zn + temp + cod + cfpom + chla + ph2 + turb + vel + pco1 + pco2,
  filename = "04_pairs_pco1_pco2_transformed.tif"
)

save_pairs_plot(
  data = alldata |>
    dplyr::mutate(
      ph1 = safe_log(ph1),
      sali = safe_log(sali),
      do = safe_sqrt(do)
    ),
  formula_obj = ~ ph1 + sali + do + pco1 + pco2,
  filename = "04_pairs_pco1_pco2_ph_sali_do_transformed.tif"
)

save_pairs_plot(
  data = alldata,
  formula_obj = ~ cd + cu + cr + zn + cod + temp + cfpom + chla + vel + ph2 + turb + sed + pred1 + pred2,
  filename = "04_pairs_pred1_pred2_raw.tif"
)

save_pairs_plot(
  data = alldata_pairs_pred_transformed,
  formula_obj = ~ cd + cu + cr + zn + cod + temp + cfpom + chla + vel + ph2 + turb + sed + pred1 + pred2,
  filename = "04_pairs_pred1_pred2_transformed.tif"
)

save_pairs_plot(
  data = alldata,
  formula_obj = ~ sali + do + ph1 + ph2 + pred1 + pred2,
  filename = "04_pairs_pred1_pred2_ph_sali_do.tif"
)



#####################################
# ---- 9. Generate scatterplot matrices for the first two PCO axes and selected environmental groups -----------

panel_cor <- function(x, y, digits = 2, method = "spearman", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  keep <- complete.cases(x, y)
  
  if (sum(keep) < 3) {
    text(0.5, 0.5, "n < 3")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x[keep], y[keep], method = method))
  p <- suppressWarnings(cor.test(x[keep], y[keep], method = method)$p.value)
  
  r_txt <- paste0("r= ", formatC(r, digits = digits, format = "f"))
  p_txt <- if (is.na(p)) "p= NA" else if (p < 0.01) "p= <0.01" else paste0("p= ", formatC(p, digits = digits, format = "f"))
  
  text(0.5, 0.8, r_txt)
  text(0.5, 0.3, p_txt)
}

alldata_screening <- alldata |>
  dplyr::mutate(
    log_vel = ifelse(is.na(vel) | vel <= 0, NA_real_, log(vel)),
    log_sali = ifelse(is.na(sali) | sali <= 0, NA_real_, log(sali)),
    log_chla = ifelse(is.na(chla) | chla <= 0, NA_real_, log(chla)),
    log_do = ifelse(is.na(do) | do <= 0, NA_real_, log(do))
  )

grDevices::tiff(
  file = fs::path(figure_dir, "04_first_two_pcos_and_environmental_variables.tif"),
  width = 14,
  height = 14,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800,
  compression = "lzw"
)
par(mar = c(5, 5, 4, 2), cex = 0.9)
pairs(
  ~ temp + cod + cfpom + chla + ph2 + turb + sed + vel + alk + cond + no3 + op + no2 + nh3 + tp + tn + on + pco1 + pco2,
  data = alldata,
  upper.panel = panel_cor,
  pch = 20,
  na.action = na.omit
)
dev.off()

grDevices::tiff(
  file = fs::path(figure_dir, "04_first_two_pcos_important_variables.tif"),
  width = 14,
  height = 14,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800,
  compression = "lzw"
)
par(mar = c(5, 5, 4, 2), cex = 0.9)
pairs(
  ~ log_vel + log_sali + temp + log_chla + log_do + pco1 + pco2,
  data = alldata_screening,
  upper.panel = panel_cor,
  pch = 20,
  na.action = na.omit
)
dev.off()

grDevices::tiff(
  file = fs::path(figure_dir, "04_first_two_pcos_nutrients.tif"),
  width = 14,
  height = 14,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800,
  compression = "lzw"
)
par(mar = c(5, 5, 4, 2), cex = 0.9)
pairs(
  ~ alk + vel + sed + ph2 + ph1 + turb + sali + do + cod + toc + cfpom + temp + no3 + no2 + nh3 + tp + tn + op + chla + cond + pco1 + pco2,
  data = alldata,
  upper.panel = panel_cor,
  pch = 20,
  na.action = na.omit
)
dev.off()


##########################################################


# ---- 10. Generate grouped environmental scatterplot matrices -----------------
wangalldata <- alldata
panel_cor <- function(x, y, digits = 2, method = "spearman", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  keep <- complete.cases(x, y)
  
  if (sum(keep) < 3) {
    text(0.5, 0.5, "n < 3")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x[keep], y[keep], method = method))
  p <- suppressWarnings(cor.test(x[keep], y[keep], method = method)$p.value)
  
  r_txt <- paste0("r= ", formatC(r, digits = digits, format = "f"))
  p_txt <- if (is.na(p)) "p= NA" else if (p < 0.01) "p= <0.01" else paste0("p= ", formatC(p, digits = digits, format = "f"))
  
  text(0.5, 0.8, r_txt)
  text(0.5, 0.3, p_txt)
}

safe_log <- function(x) ifelse(is.na(x) | x <= 0, NA_real_, log(x))

save_pairs_plot <- function(data, formula_obj, filename, width = 14, height = 14, res = 800) {
  grDevices::tiff(
    file = fs::path(figure_dir, filename),
    width = width,
    height = height,
    units = "in",
    pointsize = 12,
    bg = "transparent",
    res = res,
    compression = "lzw"
  )
  par(mar = c(5, 5, 4, 2), cex = 0.9)
  pairs(formula_obj, data = data, upper.panel = panel_cor, pch = 20, na.action = na.omit)
  dev.off()
}

alldata_groups <- alldata |>
  dplyr::mutate(
    log_vel = safe_log(vel),
    log_dflow = safe_log(dflow),
    log_depth = safe_log(depth),
    log_width = safe_log(width),
    log_chla = safe_log(chla),
    log_sed = safe_log(sed),
    log_ph2 = safe_log(ph2)
  )

grouped_pairs_specs <- list(
  list(
    data = alldata_groups,
    formula = ~ log_vel + log_dflow + log_depth + log_width + turb + log_chla + log_sed + log_ph2,
    filename = "04_velocity_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ temp + airtemp + eff + canop + veg30m + veg60m + veg90m,
    filename = "04_water_temperature_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ ph1 + ph2 + chla + alk + eff + dflow,
    filename = "04_ph_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ cod + dcod + eff + dflow + canop + veg30m + veg60m + veg90m + graz1000 + ind100 + dtoc + toc + rain1,
    filename = "04_cod_and_other_factors.tif"
  ),
  list(
    data = alldata,
    formula = ~ chla + eff + dchla + dflow + canop + no3 + tp + op + dtp + dop + dno3 + no2 + dno2 + turb + veg30m + veg60m + veg90m + vel + temp,
    filename = "04_chla_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ cfpom + dayflow + weekflow + monflow + eff + dtoc + dflow + canop + veg30m + veg60m + veg90m + graz1000 + ind100 + toc + rain1,
    filename = "04_cfpom_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ alk + dflow + dalk + dcond + cond + dtemp + temp,
    filename = "04_alkalinity_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ no3 + tn + nh3 + no2 + dtn + dnh3 + dno2 + dno3 + dflow + graz1000 + ind100 + eff,
    filename = "04_nitrate_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ op + tp + graz1000 + ind100 + dflow + dop + dtp + eff,
    filename = "04_phosphate_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ sali + graz1000 + ind100 + dflow + dsali + rain1 + eff,
    filename = "04_salinity_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ cr + cu + zn + cd + dzn + eff + dflow,
    filename = "04_heavy_metals_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ sb + dsb + eff + dflow,
    filename = "04_antimony_and_other_factors.tif"
  ),
  list(
    data = alldata,
    formula = ~ turb + vel + eff + dflow + rain1 + graz1000 + ind100,
    filename = "04_turbidity_and_other_factors.tif"
  ),
  list(
    data = wangalldata,
    formula = ~ sed + vel + eff + dflow,
    filename = "04_sediment_and_other_factors.tif"
  )
)

pairs_manifest_groups <- purrr::map_dfr(
  grouped_pairs_specs,
  function(spec) {
    save_pairs_plot(
      data = spec$data,
      formula_obj = spec$formula,
      filename = spec$filename
    )
    tibble::tibble(filename = spec$filename, status = "saved")
  }
)

readr::write_csv(pairs_manifest_groups, fs::path(output_dir, "04_grouped_pairs_manifest.csv"))




# ---- 11. Run targeted univariate screening checks for candidate environmental variables ------------------------
library(GGally)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)

vars_to_check <- c("zn", "cu", "sed", "alk", "cod", "toc", "temp", "chla", "tp")

plot_data <- alldata %>%
  select(pco1, all_of(vars_to_check)) %>%
  drop_na()

GGally::ggpairs(
  plot_data,
  lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size = 0.2)),
  upper = list(continuous = GGally::wrap("cor", size = 4)),
  diag = list(continuous = GGally::wrap("densityDiag"))
)

cor_results <- alldata %>%
  select(pco1, all_of(vars_to_check)) %>%
  pivot_longer(cols = -pco1, names_to = "variable", values_to = "value") %>%
  filter(is.finite(value), is.finite(pco1)) %>%
  group_by(variable) %>%
  group_modify(~ broom::tidy(cor.test(.x$value, .x$pco1))) %>%
  select(variable, estimate, statistic, p.value) %>%
  arrange(p.value)

print(cor_results)

facet_data <- alldata %>%
  select(pco1, all_of(vars_to_check)) %>%
  pivot_longer(cols = -pco1, names_to = "variable", values_to = "value") %>%
  filter(is.finite(value), is.finite(pco1))

ggplot(facet_data, aes(x = value, y = pco1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Environmental variables screened against PCO1",
    x = "Environmental value",
    y = "PCO1 score"
  )
