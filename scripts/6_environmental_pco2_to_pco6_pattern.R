# 06_environmental_pco2_to_pco6_pattern_modelling.R
# Author: Rezvan Hatami
# Date: 10 March 2026
#
# Purpose:
# This script models the environmental patterns of PCO2 to PCO6 after
# the main spatiotemporal analysis has been completed.
# It fits axis-specific environmental models, compares observed,
# spatiotemporal, and environmental patterns, and saves summaries,
# model objects, and publication-ready figures.
rm(list = ls())

# ---- Setup: package management ------------------------------------------------
required_pkgs <- c("dplyr", "readr", "tibble", "purrr", "fs", "here")

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
figure_dir <- fs::path(project_dir, "figures", "06")

fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- Setup: read required inputs ----------------------------------------------
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")

if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}

script3_inputs <- readRDS(script3_inputs_file)
list2env(script3_inputs, envir = .GlobalEnv)

bugenv <- cbind(wangenv, pwangbug)
envpcopred <- alldata
# ----------------------------------------------------------------------
# This chunk parses date columns when present so joins and plotting
# behave consistently across saved analysis objects.
# ----------------------------------------------------------------------
prep_date <- function(df) {
  if ("date" %in% names(df)) {
    d1 <- as.Date(as.character(df$date), format = "%Y-%m-%d")
    d2 <- as.Date(as.character(df$date), format = "%d-%m-%Y")
    df$date <- ifelse(!is.na(d1), d1, d2)
    df$date <- as.Date(df$date, origin = "1970-01-01")
  }
  df
}

# ----------------------------------------------------------------------
# This chunk safely log-transforms numeric vectors and returns NA for
# non-positive or non-numeric values.
# ----------------------------------------------------------------------
safe_log <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_num) | x_num <= 0, NA_real_, log(x_num))
}

# ----------------------------------------------------------------------
# This chunk converts compatible columns to numeric while leaving Date
# columns unchanged.
# ----------------------------------------------------------------------
prep_num <- function(df) {
  for (nm in names(df)) {
    if (inherits(df[[nm]], "Date")) next
    
    if (
      is.factor(df[[nm]]) || is.character(df[[nm]]) ||
      is.integer(df[[nm]]) || is.numeric(df[[nm]])
    ) {
      suppressWarnings({
        converted <- as.numeric(df[[nm]])
      })
      
      if (!all(is.na(converted))) {
        df[[nm]] <- converted
      }
    }
  }
  
  df
}

# ----------------------------------------------------------------------
# This chunk rebuilds the model frame and appends fitted values and
# residuals so plotting data stay aligned with each fitted model.
# ----------------------------------------------------------------------
mod_fr <- function(model, data) {
  rows <- rownames(stats::model.frame(model))
  out <- data[rows, , drop = FALSE]
  out$fitted <- as.numeric(stats::fitted(model))
  out$residual <- as.numeric(stats::residuals(model))
  out
}

# ----------------------------------------------------------------------
# This chunk fits one axis-specific linear model and returns both the
# fitted model and the aligned plotting frame.
# ----------------------------------------------------------------------
fit_ax <- function(data, response, formula_text, req_vars) {
  vars_needed <- unique(c(response, "dist", "day", req_vars))
  vars_missing <- vars_needed[!vars_needed %in% names(data)]
  
  if (length(vars_missing) > 0) {
    return(NULL)
  }
  
  model_df <- data[, vars_needed, drop = FALSE]
  
  for (nm in names(model_df)) {
    model_df[[nm]] <- suppressWarnings(as.numeric(model_df[[nm]]))
  }
  
  model_df <- stats::na.omit(model_df)
  
  if (nrow(model_df) == 0) {
    return(NULL)
  }
  
  fit <- stats::lm(stats::as.formula(formula_text), data = model_df)
  frame <- mod_fr(fit, model_df)
  
  list(model = fit, data = frame)
}

# ----------------------------------------------------------------------
# This chunk tests whether residuals from an environmental model remain
# associated with the corresponding spatiotemporal prediction.
# ----------------------------------------------------------------------
fit_res <- function(model_df, pred_df, axis_no) {
  pred_name <- paste0("pred", axis_no)
  pred_alt <- paste0("pco", axis_no, "pred")
  
  pred_var <- if (pred_name %in% names(pred_df)) {
    pred_name
  } else if (pred_alt %in% names(pred_df)) {
    pred_alt
  } else {
    NULL
  }
  
  if (is.null(pred_var)) {
    return(NULL)
  }
  
  join_vars <- intersect(c("day", "dist", "date"), names(model_df))
  join_vars <- intersect(join_vars, names(pred_df))
  
  if (length(join_vars) < 2) {
    return(NULL)
  }
  
  pred_tmp <- pred_df[, c(join_vars, pred_var), drop = FALSE]
  joined <- dplyr::left_join(model_df, pred_tmp, by = join_vars)
  joined[[pred_var]] <- suppressWarnings(as.numeric(joined[[pred_var]]))
  joined$residual <- suppressWarnings(as.numeric(joined$residual))
  joined <- joined[
    is.finite(joined$residual) & is.finite(joined[[pred_var]]), ,
    drop = FALSE
  ]
  
  if (nrow(joined) == 0) {
    return(NULL)
  }
  
  fit <- stats::lm(
    stats::as.formula(paste0("residual ~ ", pred_var)),
    data = joined
  )
  
  list(model = fit, data = joined, pred_var = pred_var)
}

# ----------------------------------------------------------------------
# This chunk plots one response variable along distance for each
# sampling day using a consistent visual key.
# ----------------------------------------------------------------------
plot_day <- function(
    df, y_var, ylab, main, ylim = NULL, days, labels, cols, pch, lty,
    leg_title = "Months", vline = 4
) {
  if (!all(c("dist", "day", y_var) %in% names(df))) {
    plot.new()
    title(main)
    text(0.5, 0.5, labels = paste("Missing required columns for", y_var))
    return(invisible(NULL))
  }
  
  plot_list <- vector("list", length(days))
  
  for (i in seq_along(days)) {
    tmp <- df[df$day == days[i], c("dist", y_var), drop = FALSE]
    tmp$dist <- suppressWarnings(as.numeric(tmp$dist))
    tmp[[y_var]] <- suppressWarnings(as.numeric(tmp[[y_var]]))
    tmp <- tmp[is.finite(tmp$dist) & is.finite(tmp[[y_var]]), , drop = FALSE]
    
    if (nrow(tmp) > 0) {
      tmp <- tmp[order(tmp$dist), , drop = FALSE]
      plot_list[[i]] <- tmp
    }
  }
  
  valid_idx <- which(!vapply(plot_list, is.null, logical(1)))
  
  if (length(valid_idx) == 0) {
    plot.new()
    title(main)
    text(0.5, 0.5, labels = "No data available for selected sampling days")
    return(invisible(NULL))
  }
  
  all_x <- unlist(lapply(plot_list[valid_idx], function(x) x$dist))
  all_y <- unlist(lapply(plot_list[valid_idx], function(x) x[[y_var]]))
  
  if (!any(is.finite(all_x)) || !any(is.finite(all_y))) {
    plot.new()
    title(main)
    text(0.5, 0.5, labels = "No finite values available for plotting")
    return(invisible(NULL))
  }
  
  if (is.null(ylim)) {
    ylim <- range(all_y, na.rm = TRUE)
  }
  
  xlim <- range(all_x, na.rm = TRUE)
  first_i <- valid_idx[1]
  first_df <- plot_list[[first_i]]
  
  plot(
    x = first_df$dist, y = first_df[[y_var]], type = "b", xlim = xlim,
    ylim = ylim, pch = pch[first_i], col = cols[first_i], cex = 1.2,
    xlab = "Distance (km)", ylab = ylab, main = main, lwd = 2,
    lty = lty[first_i]
  )
  
  if (length(valid_idx) > 1) {
    for (i in valid_idx[-1]) {
      points(
        x = plot_list[[i]]$dist, y = plot_list[[i]][[y_var]], type = "b",
        col = cols[i], pch = pch[i], cex = 1.2, lwd = 2, lty = lty[i]
      )
    }
  }
  
  legend(
    "topright", inset = c(0, 0), legend = labels[valid_idx],
    lty = lty[valid_idx], pch = pch[valid_idx], lwd = 2,
    col = cols[valid_idx], ncol = 2, horiz = FALSE, cex = 0.8,
    title = leg_title
  )
  
  abline(v = vline, lty = 2)
  
  invisible(valid_idx)
}

# ----------------------------------------------------------------------
# 1. Prepare analysis objects and derived transformed variables.
# ----------------------------------------------------------------------
bugenv <- prep_date(as.data.frame(bugenv))
envpcopred <- prep_date(as.data.frame(envpcopred))
alldata <- prep_date(as.data.frame(alldata))

bugenv <- prep_num(bugenv)
envpcopred <- prep_num(envpcopred)
alldata <- prep_num(alldata)

if ("zn" %in% names(alldata)) alldata$log_zn <- safe_log(alldata$zn)
if ("chla" %in% names(alldata)) alldata$log_chla <- safe_log(alldata$chla)
if ("cr" %in% names(alldata)) alldata$log_cr <- safe_log(alldata$cr)
if ("temp" %in% names(alldata)) alldata$log_temp <- safe_log(alldata$temp)

pref_days <- c(1, 126, 260, 336, 518)
pref_labels <- c("Dec 13", "Apr 14", "Aug 14", "Nov 14", "May 15")
pref_cols <- c("black", "gold2", "blue", "green3", "red")
pref_pch <- c(19, 15, 8, 17, 18)
pref_lty <- c(1, 2, 3, 4, 5)

avail_days <- sort(unique(stats::na.omit(as.numeric(bugenv$day))))
valid_idx <- which(pref_days %in% avail_days)

samp_days <- pref_days[valid_idx]
samp_labels <- pref_labels[valid_idx]
samp_cols <- pref_cols[valid_idx]
samp_pch <- pref_pch[valid_idx]
samp_lty <- pref_lty[valid_idx]

if (length(samp_days) == 0) {
  stop("None of the preferred sampling days were found in bugenv.")
}

# ----------------------------------------------------------------------
# 2. Define final environmental models for PCO2 to PCO6.
# ----------------------------------------------------------------------
axis_specs <- list(
  pco2 = list(
    axis_no = 2,
    formula_text = "pco2 ~ turb + temp + nh3 + temp:nh3 + ph2:temp + ph2",
    req_vars = c("temp", "nh3", "turb", "ph2"),
    fig_stub = "06_pco2"
  ),
  pco3 = list(
    axis_no = 3,
    formula_text = "pco3 ~ temp + toc + temp:toc + turb + temp:turb + log_chla + vel + log_temp",
    req_vars = c("temp", "toc", "turb", "log_chla", "vel", "log_temp"),
    fig_stub = "06_pco3"
  ),
  pco4 = list(
    axis_no = 4,
    formula_text = "pco4 ~ temp + toc + temp:toc + cod + temp:cod + cfpom + temp:cfpom + turb",
    req_vars = c("temp", "toc", "cod", "cfpom", "turb"),
    fig_stub = "06_pco4"
  ),
  pco5 = list(
    axis_no = 5,
    formula_text = "pco5 ~ ph2 + cod",
    req_vars = c("ph2", "cod"),
    fig_stub = "06_pco5"
  ),
  pco6 = list(
    axis_no = 6,
    formula_text = "pco6 ~ do + do:temp + temp",
    req_vars = c("do", "temp"),
    fig_stub = "06_pco6"
  )
)

# ----------------------------------------------------------------------
# 3. Fit environmental models and residual prediction tests for each
#    PCO axis from PCO2 to PCO6.
# ----------------------------------------------------------------------
axis_mods <- list()
axis_res <- list()
model_sum <- tibble::tibble()
resid_sum <- tibble::tibble()

for (axis_name in names(axis_specs)) {
  spec <- axis_specs[[axis_name]]
  
  axis_fit <- fit_ax(
    data = alldata, response = axis_name, formula_text = spec$formula_text,
    req_vars = spec$req_vars
  )
  
  if (is.null(axis_fit)) {
    message(
      "Skipping ", axis_name,
      " because required variables are missing or no complete cases exist."
    )
    axis_mods[[axis_name]] <- NULL
    axis_res[[axis_name]] <- NULL
    next
  }
  
  axis_mods[[axis_name]] <- axis_fit
  
  fit_obj <- axis_fit$model
  fit_df <- axis_fit$data
  
  model_sum <- dplyr::bind_rows(
    model_sum,
    tibble::tibble(
      axis = axis_name,
      axis_no = spec$axis_no,
      formula = paste(deparse(stats::formula(fit_obj)), collapse = " "),
      n = stats::nobs(fit_obj),
      r_squared = summary(fit_obj)$r.squared,
      adj_r_squared = summary(fit_obj)$adj.r.squared,
      sigma = summary(fit_obj)$sigma,
      aic = stats::AIC(fit_obj)
    )
  )
  
  axis_res[[axis_name]] <- fit_res(
    model_df = fit_df, pred_df = envpcopred, axis_no = spec$axis_no
  )
  
  if (!is.null(axis_res[[axis_name]])) {
    rp_fit <- axis_res[[axis_name]]$model
    rp_var <- axis_res[[axis_name]]$pred_var
    
    resid_sum <- dplyr::bind_rows(
      resid_sum,
      tibble::tibble(
        axis = axis_name,
        axis_no = spec$axis_no,
        predictor = rp_var,
        n = stats::nobs(rp_fit),
        slope = unname(stats::coef(rp_fit)[2]),
        intercept = unname(stats::coef(rp_fit)[1]),
        r_squared = summary(rp_fit)$r.squared,
        adj_r_squared = summary(rp_fit)$adj.r.squared,
        p_value = summary(rp_fit)$coefficients[2, 4]
      )
    )
  }
}

readr::write_csv(model_sum, fs::path(output_dir, "06_axis_model_summary.csv"))
readr::write_csv(resid_sum, fs::path(output_dir, "06_residual_prediction_summary.csv"))

saveRDS(
  lapply(axis_mods, function(x) if (is.null(x)) NULL else x$model),
  fs::path(output_dir, "06_axis_models.rds")
)

saveRDS(
  lapply(axis_mods, function(x) if (is.null(x)) NULL else x$data),
  fs::path(output_dir, "06_axis_model_data.rds")
)

saveRDS(
  lapply(axis_res, function(x) if (is.null(x)) NULL else x$model),
  fs::path(output_dir, "06_residual_prediction_models.rds")
)

# ----------------------------------------------------------------------
# 4. Produce one figure set per axis, including fitted and residual
#    plots, comparison plots, and residual prediction plots.
# ----------------------------------------------------------------------
fig_manifest <- tibble::tibble(
  figure_file = character(),
  axis = character(),
  description = character()
)

for (axis_name in names(axis_specs)) {
  spec <- axis_specs[[axis_name]]
  
  if (is.null(axis_mods[[axis_name]])) {
    next
  }
  
  axis_no <- spec$axis_no
  fit_df <- axis_mods[[axis_name]]$data
  pred_var <- paste0("pred", axis_no)
  
  if (!(pred_var %in% names(envpcopred))) {
    pred_var <- paste0("pco", axis_no, "pred")
  }
  
  obs_var <- axis_name
  
  fig_a <- fs::path(figure_dir, paste0(spec$fig_stub, "_fitted_and_residual.jpg"))
  
  jpeg(
    filename = fig_a, width = 10, height = 10, units = "in",
    pointsize = 12, bg = "transparent", res = 800
  )
  
  par(mfrow = c(2, 1), mar = c(5, 5, 4, 2), cex = 0.9,
      cex.axis = 0.9, las = 1)
  
  plot_day(
    df = fit_df, y_var = "fitted", ylab = toupper(axis_name),
    main = paste0("a) Environmental prediction for ", toupper(axis_name)),
    ylim = NULL, days = samp_days, labels = samp_labels, cols = samp_cols,
    pch = samp_pch, lty = samp_lty
  )
  
  plot_day(
    df = fit_df, y_var = "residual",
    ylab = paste("Residual", toupper(axis_name)),
    main = paste0("b) Residual pattern for ", toupper(axis_name)),
    ylim = NULL, days = samp_days, labels = samp_labels, cols = samp_cols,
    pch = samp_pch, lty = samp_lty
  )
  
  dev.off()
  
  fig_manifest <- dplyr::bind_rows(
    fig_manifest,
    tibble::tibble(
      figure_file = fs::path_file(fig_a),
      axis = axis_name,
      description = "Environmental fitted values and residual structure"
    )
  )
  
  fig_b <- fs::path(figure_dir, paste0(spec$fig_stub, "_comparison.jpg"))
  
  jpeg(
    filename = fig_b, width = 10, height = 14, units = "in",
    pointsize = 12, bg = "transparent", res = 800
  )
  
  par(mfrow = c(3, 1), mar = c(4.5, 4.5, 2.5, 3), cex = 1.0,
      cex.axis = 0.9, las = 1)
  
  plot_day(
    df = bugenv, y_var = obs_var, ylab = toupper(axis_name),
    main = paste0("a) Observed ", toupper(axis_name),
                  " plotted against spatial position"),
    ylim = NULL, days = samp_days, labels = samp_labels, cols = samp_cols,
    pch = samp_pch, lty = samp_lty
  )
  
  if (pred_var %in% names(envpcopred)) {
    plot_day(
      df = envpcopred, y_var = pred_var,
      ylab = paste("Predicted", toupper(axis_name)),
      main = paste0("b) Spatiotemporal prediction for ", toupper(axis_name)),
      ylim = NULL, days = samp_days, labels = samp_labels, cols = samp_cols,
      pch = samp_pch, lty = samp_lty
    )
  } else {
    plot.new()
    title(paste0("b) Spatiotemporal prediction unavailable for ",
                 toupper(axis_name)))
  }
  
  plot_day(
    df = fit_df, y_var = "fitted", ylab = paste("Fitted", toupper(axis_name)),
    main = paste0("c) Environmental fitted values for ", toupper(axis_name)),
    ylim = NULL, days = samp_days, labels = samp_labels, cols = samp_cols,
    pch = samp_pch, lty = samp_lty
  )
  
  dev.off()
  
  fig_manifest <- dplyr::bind_rows(
    fig_manifest,
    tibble::tibble(
      figure_file = fs::path_file(fig_b),
      axis = axis_name,
      description = "Observed, spatiotemporal, and environmental comparison"
    )
  )
  
  fig_c <- fs::path(figure_dir, paste0(spec$fig_stub, "_residual_vs_prediction.jpg"))
  
  jpeg(
    filename = fig_c, width = 8, height = 6, units = "in",
    pointsize = 12, bg = "transparent", res = 800
  )
  
  rp_res <- axis_res[[axis_name]]
  
  if (!is.null(rp_res)) {
    rp_df <- rp_res$data
    rp_fit <- rp_res$model
    rp_var <- rp_res$pred_var
    
    plot(
      rp_df[[rp_var]], rp_df$residual, pch = 19, cex = 1,
      xlab = paste("Spatiotemporal prediction (", rp_var, ")", sep = ""),
      ylab = paste("Residual", toupper(axis_name)),
      main = paste0("Residual association for ", toupper(axis_name))
    )
    abline(rp_fit, lwd = 2, lty = 2)
  } else {
    plot.new()
    title(paste0("Residual prediction test unavailable for ",
                 toupper(axis_name)))
  }
  
  dev.off()
  
  fig_manifest <- dplyr::bind_rows(
    fig_manifest,
    tibble::tibble(
      figure_file = fs::path_file(fig_c),
      axis = axis_name,
      description = "Residuals plotted against spatiotemporal predictions"
    )
  )
}

readr::write_csv(fig_manifest, fs::path(output_dir, "06_figure_manifest.csv"))

# ----------------------------------------------------------------------
# 5. Save the compact axis formula table used in this script.
# ----------------------------------------------------------------------
axis_formula_table <- tibble::tibble(
  axis = names(axis_specs),
  axis_no = vapply(axis_specs, function(x) x$axis_no, numeric(1)),
  formula = vapply(axis_specs, function(x) x$formula_text, character(1))
)

readr::write_csv(axis_formula_table, fs::path(output_dir, "06_axis_formulas.csv"))

cat("Script 06 completed successfully.\n")