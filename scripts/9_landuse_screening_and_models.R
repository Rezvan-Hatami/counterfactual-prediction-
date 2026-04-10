# 09_landuse_screening_and_models.R
# Author: Rezvan Hatami
# Date: 13-03-26
#
# Purpose:
# This script reproduces the original land-use screening workflow and
# scatterplot-matrix checks in a clean, self-sufficient, and reproducible form.
# It prints analytical results, and saves figure outputs.

rm(list = ls())

# ---- Setup: package management ------------------------------------------------
required_pkgs <- c("dplyr", "purrr", "fs", "here")

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
figure_dir <- fs::path(project_dir, "figures", "09")

fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- Setup: helper functions --------------------------------------------------
# This chunk converts compatible columns to numeric while leaving Date columns
# unchanged so modelling and plotting use consistent column types.
prep_num <- function(df) {
  for (nm in names(df)) {
    if (inherits(df[[nm]], "Date")) next
    suppressWarnings({
      converted <- as.numeric(df[[nm]])
    })
    if (!all(is.na(converted))) {
      df[[nm]] <- converted
    }
  }
  df
}

# This chunk applies log safely to positive values.
safe_log <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_num) | x_num <= 0, NA_real_, log(x_num))
}

# This chunk returns only variables that exist in the supplied dataset.
pick_vars <- function(data, vars) {
  vars[vars %in% names(data)]
}

# This chunk checks that required objects exist after the shared state file is
# loaded and stops with a precise message if anything is missing.
check_obj <- function(nms, env = parent.frame()) {
  missing <- nms[!vapply(nms, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(missing) > 0) {
    stop("Missing required objects in script state: ", paste(missing, collapse = ", "))
  }
}

# This chunk builds the upper-panel correlation display used in pairs plots and
# suppresses zero-variance warnings.
panel_cor <- function(x, y, digits = 2, method = "spearman", ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  keep <- complete.cases(x, y)
  if (sum(keep) < 3) {
    text(0.5, 0.5, "n < 3")
    return(invisible(NULL))
  }
  
  x2 <- x[keep]
  y2 <- y[keep]
  
  if (sd(x2) == 0 || sd(y2) == 0) {
    text(0.5, 0.5, "sd = 0")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x2, y2, method = method))
  p <- suppressWarnings(cor.test(x2, y2, method = method)$p.value)
  
  r_txt <- paste0("r = ", formatC(r, digits = digits, format = "f"))
  p_txt <- if (is.na(p)) {
    "p = NA"
  } else if (p < 0.01) {
    "p < 0.01"
  } else {
    paste0("p = ", formatC(p, digits = digits, format = "f"))
  }
  
  text(0.5, 0.8, r_txt)
  text(0.5, 0.4, p_txt)
}

# This chunk saves one scatterplot matrix from a variable specification.
save_pairs <- function(data, vars, file_name, width = 14, height = 14, res = 800) {
  vars <- pick_vars(data, vars)
  
  if (length(vars) < 2) {
    return(invisible(NULL))
  }
  
  tiff(
    file = fs::path(figure_dir, file_name),
    width = width,
    height = height,
    units = "in",
    pointsize = 12,
    bg = "transparent",
    res = res,
    compression = "lzw"
  )
  
  par(mar = c(0.2, 0.2, 0.2, 0.2), cex = 0.6, cex.axis = 0.8, las = 1)
  
  pairs(
    stats::as.formula(paste("~", paste(vars, collapse = "+"))),
    data = data,
    upper.panel = panel_cor,
    pch = 20,
    na.action = na.omit
  )
  
  dev.off()
  invisible(file_name)
}

# This chunk fits one land-use screening model and prints the selected model,
# ANOVA table, and summary
run_land_model <- function(label, formula_text, data) {
  vars_needed <- all.vars(stats::as.formula(formula_text))
  vars_needed <- pick_vars(data, vars_needed)
  work_df <- prep_num(data[, vars_needed, drop = FALSE])
  work_df <- stats::na.omit(work_df)
  
  if (length(vars_needed) == 0) {
    cat("\n============================================================\n")
    cat(label, "\n")
    cat("============================================================\n")
    cat("No matching variables were found for this model.\n")
    return(invisible(NULL))
  }
  
  fit <- stats::lm(stats::as.formula(formula_text), data = work_df)
  
  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  cat("\nModel formula:\n")
  print(stats::formula(fit))
  cat("\nANOVA:\n")
  print(anova(fit))
  cat("\nSUMMARY:\n")
  print(summary(fit, corr = TRUE))
  
  invisible(fit)
}

# This chunk creates the log-transformed environmental variables used in the revised
# causal-diagram scatterplot matrices.
mk_log_df <- function(df) {
  dplyr::mutate(
    df,
    cond = safe_log(cond),
    toc = safe_log(toc),
    chla = safe_log(chla),
    zn = safe_log(zn),
    temp = safe_log(temp),
    tp = safe_log(tp),
    no3 = safe_log(no3),
    alk = safe_log(alk),
    turb = safe_log(turb)
  )
}

# This chunk checks all supplied environmental and land-use combinations and
# prints the Spearman screening results directly to the console.
run_cor_grid <- function(data, env_vars, landuse_vars, label) {
  env_vars <- pick_vars(data, env_vars)
  landuse_vars <- pick_vars(data, landuse_vars)
  
  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  
  if (length(env_vars) == 0 || length(landuse_vars) == 0) {
    cat("No matching variables found.\n")
    return(invisible(NULL))
  }
  
  for (env_var in env_vars) {
    for (land_var in landuse_vars) {
      x <- suppressWarnings(as.numeric(data[[env_var]]))
      y <- suppressWarnings(as.numeric(data[[land_var]]))
      keep <- is.finite(x) & is.finite(y)
      
      if (sum(keep) < 3) next
      if (sd(x[keep]) == 0 || sd(y[keep]) == 0) next
      
      test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
      
      cat(
        env_var, " vs ", land_var,
        " | rho = ", formatC(unname(test$estimate), digits = 3, format = "f"),
        " | p = ", formatC(test$p.value, digits = 3, format = "f"),
        " | n = ", sum(keep), "\n",
        sep = ""
      )
    }
  }
  
  invisible(NULL)
}

# ---- Setup: read required inputs ----------------------------------------------
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")

if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}
if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}

script3_inputs <- readRDS(script3_inputs_file)
if (!is.list(script3_inputs)) {
  stop("The script state file must be a named list: ", script3_inputs_file)
}

local_env <- list2env(script3_inputs, parent = emptyenv())
check_obj(c("pwangbug"), env = local_env)

pwangbug <- get("pwangbug", envir = local_env)
wangenv <- readRDS(wangenv_file)

if (exists("envpcopred", envir = local_env, inherits = FALSE)) {
  envpcopred <- get("envpcopred", envir = local_env)
} else {
  envpcopred <- NULL
}

if (exists("alldata", envir = local_env, inherits = FALSE)) {
  alldata_state <- get("alldata", envir = local_env)
} else {
  alldata_state <- NULL
}

# ----------------------------------------------------------------------
# 1. Prepare the analysis tables and print the core structure.
# ----------------------------------------------------------------------
wangenv <- prep_num(as.data.frame(wangenv))
pwangbug <- prep_num(as.data.frame(pwangbug))

if (!is.null(envpcopred)) {
  envpcopred <- prep_num(as.data.frame(envpcopred))
  alldata <- cbind(envpcopred, pwangbug)
} else if (!is.null(alldata_state)) {
  alldata <- prep_num(as.data.frame(alldata_state))
} else {
  stop(
    "Neither envpcopred nor alldata was found in 03_environmental_script_inputs.rds."
  )
}

cat("\nDIMENSION OF alldata:\n")
print(dim(alldata))

cat("\nNAMES OF alldata:\n")
print(names(alldata))

cat("\nNAMES OF wangenv:\n")
print(names(wangenv))

landuse1 <- wangenv[c(
  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
  28, 35, 37, 38, 43, 45, 47, 48, 53, 54, 56, 57, 62, 64, 66, 67, 72,
  73, 75, 76, 81, 82, 84, 90, 91, 93, 94, 99, 100, 102, 103, 108, 109,
  111, 112, 117, 118, 120, 121, 126, 127, 129, 130, 135, 136, 138, 139,
  144, 177
)]

landuse2 <- log1p(landuse1)

cat("\nDIMENSION OF landuse1:\n")
print(dim(landuse1))

cat("\nNAMES OF landuse1:\n")
print(names(landuse1))

cat("\nDIMENSION OF landuse2:\n")
print(dim(landuse2))

cat("\nNAMES OF landuse2:\n")
print(names(landuse2))

# ----------------------------------------------------------------------
# 2. Fit the original land-use models in the original order and print
#    ANOVA and summary outputs to the console.
# ----------------------------------------------------------------------
land_models <- list()

land_models[["no2"]] <- run_land_model(
  "Model for no2",
  "log1p(no2) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["no3"]] <- run_land_model(
  "Model for no3",
  "log1p(no3) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["nh3_1"]] <- run_land_model(
  "Model for nh3",
  "log1p(nh3) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["nh3_2"]] <- run_land_model(
  "Repeated model for nh3",
  "log1p(nh3) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["op"]] <- run_land_model(
  "Model for op",
  "log1p(op) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["tp"]] <- run_land_model(
  "Model for tp",
  "log1p(tp) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["tn"]] <- run_land_model(
  "Model for tn",
  "log1p(tn) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["sali"]] <- run_land_model(
  "Model for sali",
  "log1p(sali) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["cond"]] <- run_land_model(
  "Model for cond",
  "log1p(cond) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["turb"]] <- run_land_model(
  "Model for turb",
  "log1p(turb) ~ eff + log1p(rec300) + log1p(treat100) + log1p(graz100) + log1p(treat60w)",
  alldata
)

land_models[["cod"]] <- run_land_model(
  "Model for cod",
  "sqrt(cod) ~ eff + log1p(rec300) + log1p(treat100) + log1p(treat100w) + log1p(graz100) + log1p(treat60w) + log1p(treat60)",
  alldata
)

land_models[["toc"]] <- run_land_model(
  "Model for toc",
  "log1p(toc) ~ eff + log1p(rec300) + log1p(treat100) + log1p(treat100w) + log1p(graz100) + log1p(treat60w) + log1p(treat60)",
  alldata
)

land_models[["zn"]] <- run_land_model(
  "Model for zn",
  "log1p(zn) ~ eff + log1p(res1000) + log1p(treat1000) + log1p(treat500) + log1p(road1000) + log1p(ind)",
  alldata
)

# ----------------------------------------------------------------------
# 3. Produce the original land-use scatterplot matrices in the original
#    order for metals and land-use variables.
# ----------------------------------------------------------------------
metal_pairs_specs <- list(
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz", "ind", "res", "for.", "rec", "road", "treat"), file = "09_pairs_metals_landuse_base.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz1000", "ind1000", "res1000", "for1000", "rec1000", "road1000", "treat1000"), file = "09_pairs_metals_landuse_1000.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz500", "ind500", "res500", "for500", "rec500", "road500", "treat500"), file = "09_pairs_metals_landuse_500.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz300", "ind300", "res300", "for300", "rec300", "road300", "treat300"), file = "09_pairs_metals_landuse_300.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz100", "ind100", "res100", "for100", "rec100", "road100", "treat100"), file = "09_pairs_metals_landuse_100.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz60", "ind60", "for60", "road60", "treat60"), file = "09_pairs_metals_landuse_60.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz500w", "ind500w", "res500w", "for500w", "rec500w", "road500w", "treat500w"), file = "09_pairs_metals_landuse_500w.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz400w", "ind400w", "res400w", "for400w", "rec400w", "road400w", "treat400w", "tra400w"), file = "09_pairs_metals_landuse_400w.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz300w", "ind300w", "res300w", "for300w", "rec300w", "road300w", "treat300w", "tra300w"), file = "09_pairs_metals_landuse_300w.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz200w", "ind200w", "res200w", "for200w", "rec200w", "road200w", "treat200w", "tra200w"), file = "09_pairs_metals_landuse_200w.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz100w", "ind100w", "res100w", "for100w", "rec100w", "road100w", "treat100w"), file = "09_pairs_metals_landuse_100w.tif"),
  list(data = wangenv, vars = c("cd", "cr", "cu", "zn", "graz60w", "ind60w", "res60w", "for60w", "rec60w", "road60w", "treat60w"), file = "09_pairs_metals_landuse_60w.tif")
)

purrr::walk(metal_pairs_specs, function(spec) {
  save_pairs(spec$data, spec$vars, spec$file)
})

# ----------------------------------------------------------------------
# 4. Produce the original land-use scatterplot matrices in the original
#    order for nutrients, chemistry, and physical variables.
# ----------------------------------------------------------------------
core_vars <- list(
  nutrients = c("no3", "no2", "nh3", "tn", "tp", "op", "ph2"),
  chemistry = c("alk", "toc", "cod", "sb", "cond", "temp", "sali", "do", "ph1", "ph2"),
  physical = c("turb", "vel", "chla", "cfpom", "sed")
)

landuse_specs <- list(
  base = c("graz", "res", "for.", "rec", "road", "treat"),
  `1000` = c("graz1000", "ind1000", "res1000", "for1000", "rec1000", "road1000", "treat1000"),
  `500` = c("graz500", "ind500", "res500", "for500", "rec500", "road500", "treat500"),
  `300` = c("graz300", "ind300", "res300", "for300", "rec300", "road300", "treat300"),
  `100` = c("graz100", "ind100", "res100", "for100", "rec100", "road100", "treat100"),
  `60` = c("graz60", "ind60", "res60", "for60", "road60", "treat60"),
  `500w` = c("graz500w", "ind500w", "res500w", "for500w", "rec500w", "road500w", "treat500w"),
  `400w` = c("graz400w", "ind400w", "res400w", "for400w", "rec400w", "road400w", "treat400w", "tra400w"),
  `300w` = c("graz300w", "ind300w", "res300w", "for300w", "rec300w", "road300w", "treat300w", "tra300w"),
  `200w` = c("graz200w", "ind200w", "res200w", "for200w", "rec200w", "road200w", "treat200w", "tra200w"),
  `100w` = c("graz100w", "ind100w", "res100w", "for100w", "rec100w", "road100w", "treat100w"),
  `60w` = c("graz60w", "ind60w", "res60w", "for60w", "rec60w", "road60w", "treat60w")
)

env_group_order <- c("nutrients", "chemistry", "physical")
scale_order <- c("base", "1000", "500", "300", "100", "60", "500w", "400w", "300w", "200w", "100w", "60w")

env_pairs_specs <- list()

for (group_name in env_group_order) {
  for (scale_name in scale_order) {
    env_data <- if (group_name == "nutrients" && scale_name == "base") alldata else wangenv
    env_vars <- core_vars[[group_name]]
    land_vars <- landuse_specs[[scale_name]]
    vars <- c(env_vars, land_vars)
    vars <- vars[vars %in% names(env_data)]
    
    file_name <- paste0("09_pairs_", group_name, "_landuse_", scale_name, ".tif")
    
    env_pairs_specs[[length(env_pairs_specs) + 1]] <- list(
      data = env_data,
      vars = vars,
      file = file_name
    )
  }
}

purrr::walk(env_pairs_specs, function(spec) {
  save_pairs(spec$data, spec$vars, spec$file)
})

# ----------------------------------------------------------------------
# 5. Produce the revised causal-diagram scatterplot matrices in the
#    original order for raw and log-transformed scales.
# ----------------------------------------------------------------------
wangenv_log <- mk_log_df(wangenv)

revised_core_vars <- c(
  "cond", "toc", "chla", "zn", "temp", "tp", "no3", "alk", "turb"
)

revised_landuse_specs <- list(
  base = c("graz", "res", "for.", "rec", "road", "treat"),
  `1000` = c("graz1000", "ind1000", "res1000", "for1000", "rec1000", "road1000", "treat1000"),
  `500` = c("graz500", "ind500", "res500", "for500", "rec500", "road500", "treat500"),
  `300` = c("graz300", "ind300", "res300", "for300", "rec300", "road300", "treat300"),
  `100` = c("graz100", "ind100", "res100", "for100", "rec100", "road100", "treat100"),
  `60` = c("graz60", "ind60", "res60", "for60", "road60", "treat60"),
  `500w` = c("graz500w", "ind500w", "res500w", "for500w", "rec500w", "road500w", "treat500w"),
  `400w` = c("graz400w", "ind400w", "res400w", "for400w", "rec400w", "road400w", "treat400w"),
  `300w` = c("graz300w", "ind300w", "res300w", "for300w", "rec300w", "road300w", "treat300w"),
  `200w` = c("graz200w", "ind200w", "res200w", "for200w", "rec200w", "road200w", "treat200w"),
  `100w` = c("graz100w", "ind100w", "res100w", "for100w", "rec100w", "road100w", "treat100w"),
  `60w` = c("graz60w", "ind60w", "res60w", "for60w", "rec60w", "road60w", "treat60w")
)

revised_scale_order <- c("base", "1000", "500", "300", "100", "60", "500w", "400w", "300w", "200w", "100w", "60w")
revised_pairs_specs <- list()

for (scale_name in revised_scale_order) {
  vars <- c(revised_core_vars, revised_landuse_specs[[scale_name]])
  
  revised_pairs_specs[[length(revised_pairs_specs) + 1]] <- list(
    data = wangenv,
    vars = vars[vars %in% names(wangenv)],
    file = paste0("09_pairs_revised_causal_raw_", scale_name, ".tif")
  )
  
  revised_pairs_specs[[length(revised_pairs_specs) + 1]] <- list(
    data = wangenv_log,
    vars = vars[vars %in% names(wangenv_log)],
    file = paste0("09_pairs_revised_causal_log_", scale_name, ".tif")
  )
}

purrr::walk(revised_pairs_specs, function(spec) {
  save_pairs(spec$data, spec$vars, spec$file)
})

# ----------------------------------------------------------------------
# 6. Check all land-use and environmental variable combinations that are
#    represented in the original script and print the results to console.
# ----------------------------------------------------------------------
screen_groups <- list(
  metal = c("cd", "cr", "cu", "zn"),
  nutrient = c("no3", "no2", "nh3", "tn", "tp", "op", "ph2"),
  chemistry = c("alk", "toc", "cod", "sb", "cond", "temp", "sali", "do", "ph1", "ph2"),
  physical = c("turb", "vel", "chla", "cfpom", "sed")
)

screen_landuse_specs <- list(
  base = c("graz", "ind", "res", "for.", "rec", "road", "treat"),
  `1000` = c("graz1000", "ind1000", "res1000", "for1000", "rec1000", "road1000", "treat1000", "tra1000", "un1000"),
  `500` = c("graz500", "ind500", "res500", "for500", "rec500", "road500", "treat500", "tra500", "un500"),
  `300` = c("graz300", "ind300", "res300", "for300", "rec300", "road300", "treat300", "tra300", "un300"),
  `100` = c("graz100", "ind100", "res100", "for100", "rec100", "road100", "treat100", "un100"),
  `60` = c("graz60", "ind60", "res60", "for60", "road60", "treat60", "un60"),
  `500w` = c("graz500w", "ind500w", "res500w", "for500w", "rec500w", "road500w", "treat500w", "un500w"),
  `400w` = c("graz400w", "ind400w", "res400w", "for400w", "rec400w", "road400w", "treat400w", "tra400w", "un400w"),
  `300w` = c("graz300w", "ind300w", "res300w", "for300w", "rec300w", "road300w", "treat300w", "tra300w", "un300w"),
  `200w` = c("graz200w", "ind200w", "res200w", "for200w", "rec200w", "road200w", "treat200w", "tra200w", "un200w"),
  `100w` = c("graz100w", "ind100w", "res100w", "for100w", "rec100w", "road100w", "treat100w", "un100w"),
  `60w` = c("graz60w", "ind60w", "res60w", "for60w", "rec60w", "road60w", "treat60w", "un60w")
)

screen_group_order <- c("metal", "nutrient", "chemistry", "physical")
screen_scale_order <- c("base", "1000", "500", "300", "100", "60", "500w", "400w", "300w", "200w", "100w", "60w")

for (scale_name in screen_scale_order) {
  for (group_name in screen_group_order) {
    run_cor_grid(
      wangenv,
      screen_groups[[group_name]],
      screen_landuse_specs[[scale_name]],
      paste(tools::toTitleCase(group_name), "screening:", scale_name)
    )
  }
}

cat("\nScript 09 completed successfully.\n")
