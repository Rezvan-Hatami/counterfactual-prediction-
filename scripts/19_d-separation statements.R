# 19_d_separation_statements.R
# Author: Rezvan Hatami
# Date: 2026-03-28

rm(list = ls())

# ----------------------------------------------------------------------
# 1. Package setup
# ----------------------------------------------------------------------
required_pkgs <- c(
  "fs",
  "here",
  "dagR",
  "gRbase",
  "gRain",
  "igraph",
  "vegan",
  "lattice"
)

optional_pkgs <- c(
  "ggm",
  "pcalg",
  "piecewiseSEM"
)

missing_required_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required_pkgs) > 0) {
  stop(
    "These packages must be installed before running the script: ",
    paste(missing_required_pkgs, collapse = ", ")
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

missing_optional_pkgs <- optional_pkgs[
  !vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_optional_pkgs) > 0) {
  cat(
    "\nOptional packages not installed: ",
    paste(missing_optional_pkgs, collapse = ", "),
    "\nThe script will continue because they are not required for the current workflow.\n",
    sep = ""
  )
}

options(stringsAsFactors = FALSE)
options(scipen = 1)
set.seed(123)

# ----------------------------------------------------------------------
# 2. Paths
# ----------------------------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
# script_output_dir <- fs::path(output_dir, "19")

fs::dir_create(output_dir)
# fs::dir_create(script_output_dir)

# ----------------------------------------------------------------------
# 3. Read inputs
# ----------------------------------------------------------------------
wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")

script3_objs <- readRDS(file.path(output_dir, "03_environmental_script_inputs.rds"))

wangbug.BC <- script3_objs$wangbug.BC
wangbug    <- script3_objs$wangbug
wangenv    <- script3_objs$wangenv
alldata    <- script3_objs$alldata

if (!fs::file_exists(wangbug_file)) {
  stop("Missing input file: ", wangbug_file)
}
if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}
if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}

wangbug <- readRDS(wangbug_file)
wangenv <- readRDS(wangenv_file)
script3_inputs <- readRDS(script3_inputs_file)
alldata <- script3_inputs$alldata

# ----------------------------------------------------------------------
# 4. Helper functions
# ----------------------------------------------------------------------
prepare_numeric_columns <- function(df) {
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

print_model_block <- function(model_name, model_obj) {
  cat("\n")
  cat("============================================================\n")
  cat(model_name, "\n")
  cat("============================================================\n")
  print(summary(model_obj))
  print(anova(model_obj))
}

run_dsep_test <- function(base_formula, test_formula, data, label, subset_expr = NULL) {
  if (!is.null(subset_expr)) {
    data <- subset_expr(data)
  }
  
  tempo <- lm(base_formula, data = data)
  independence.test1 <- lm(test_formula, data = data)
  
  cat("\n")
  cat("============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  print(summary(independence.test1))
  print(anova(independence.test1))
  print(anova(tempo, independence.test1))
  
  invisible(
    list(
      base_model = tempo,
      test_model = independence.test1
    )
  )
}

run_capscale_residual_check <- function(
    cap_formula,
    design_formula,
    x,
    xlab,
    ylab = "pco1 residuals",
    subset_cols = c("temp", "cond", "toc", "chla")
) {
  wang.cap3 <- vegan::capscale(
    formula = cap_formula,
    data = alldata,
    comm = wangbug,
    add = TRUE,
    na.action = na.omit
  )
  
  print(anova(wang.cap3))
  print(summary(wang.cap3))
  
  logchla <- log(alldata$chla)
  datalogchla <- cbind(logchla, alldata)
  datalogchla1 <- datalogchla[complete.cases(datalogchla[, subset_cols]), ]
  
  design <- scale(
    model.matrix(design_formula, data = alldata),
    center = FALSE,
    scale = FALSE
  )
  
  wangbug1 <- wangbug[9:40, , drop = FALSE]
  wangbug1.BC <- vegan::vegdist(sqrt(wangbug1))
  n <- dim(wangbug1)[1]
  p <- n - 1
  wangbug1.mds <- cmdscale(
    wangbug1.BC,
    k = p,
    eig = TRUE,
    add = TRUE,
    x.ret = FALSE
  )
  
  pco.predict <- qr.fitted(qr(design), wangbug1.mds$points)
  pco.resid <- wangbug1.mds$points - pco.predict
  
  plot(x, pco.resid[, 1], type = "p", xlab = xlab, ylab = ylab)
  print(lattice::xyplot(pco.resid[, 1] ~ x, xlab = xlab, ylab = ylab))
  
  independence.check <- lm(pco.resid[, 1] ~ x)
  print(summary(independence.check))
  print(anova(independence.check))
  
  invisible(
    list(
      capscale_model = wang.cap3,
      independence_check = independence.check,
      pco_resid = pco.resid,
      datalogchla1 = datalogchla1
    )
  )
}

# ----------------------------------------------------------------------
# 5. Prepare analysis objects
# ----------------------------------------------------------------------
wangbug <- as.data.frame(wangbug)
wangenv <- prepare_numeric_columns(as.data.frame(wangenv))
alldata <- prepare_numeric_columns(as.data.frame(alldata))

if (!"time" %in% names(alldata) && "day" %in% names(alldata)) {
  unique_days <- sort(unique(stats::na.omit(alldata$day)))
  alldata$time <- match(alldata$day, unique_days)
}

# ----------------------------------------------------------------------
# 6. DAG specification notes
# ----------------------------------------------------------------------
if ("dag.init" %in% getNamespaceExports("dagR")) {
  cat("\ndagR is available, but manual DAG reconstruction is skipped because the old DAG() constructor used in the original script is not available in the current dagR version.\n")
  cat("D-separation regression tests will continue.\n")
  
  # sink(fs::path(script_output_dir, "19_dag_notes.txt"), append = TRUE)
  # cat("The original script used DAG(...), but the installed dagR version does not export DAG().\n")
  # cat("Current dagR exports dag.init(), dag.draw(), and related functions instead.\n")
  # cat("DAG drawing and basis-set generation were skipped in this run.\n")
  # sink(NULL)
} else {
  cat("\ndagR helper functions are not available for DAG construction in this session.\n")
  cat("DAG drawing and basis-set generation were skipped.\n")
}

# ----------------------------------------------------------------------
# 7. d-separation checks
# ----------------------------------------------------------------------

# 2
run_dsep_test(
  ph2 ~ rain2,
  ph2 ~ rain2 + time,
  alldata,
  "#2 d-separation (ph2 des time | rain2)"
)

# 8
run_dsep_test(
  dayflow ~ dflow + rain3,
  dayflow ~ dflow + rain3 + time,
  alldata,
  "#8 d-separation (dayflow des time | rain3, dflow)"
)

# 14
run_dsep_test(
  vel ~ dist,
  vel ~ dist + time,
  alldata,
  "#14 d-separation (vel des time | dist)"
)

# 15
tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel), data = alldata)
print_model_block("#15 base turb model", tempo)
run_dsep_test(
  log(turb) ~ log(rain3) + log(dayflow) + log(vel),
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + time,
  alldata,
  "#15 d-separation (turb des time | rain3, dayflow, vel)"
)

# 16
run_dsep_test(
  res300 ~ dist,
  res300 ~ dist + time,
  alldata,
  "#16 d-separation (res300 des time | dist)"
)

# 17
tempo <- lm(temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3, data = alldata)
print_model_block("#17 base temp model", tempo)
run_dsep_test(
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3,
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3 + time,
  alldata,
  "#17 d-separation (temp des time | airtemp, rain3, dayflow, res300)"
)

# 19
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#19 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time,
  alldata,
  "#19 d-separation (chla des time | solar, dayflow, dchla, turb, temp, eff, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 20
run_dsep_test(
  log(zn) ~ eff,
  log(zn) ~ eff + time,
  alldata,
  "#20 d-separation (zn des time | eff)"
)

# 21
run_capscale_residual_check(
  wangbug.BC ~ time + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn)),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn),
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "time"],
  xlab = "time"
)

# 76
run_dsep_test(
  ph2 ~ rain2 + rain3 + dflow,
  ph2 ~ rain2 + rain3 + dflow + dayflow,
  alldata,
  "#76 d-separation (ph2 des dayflow | rain2, rain3, dflow)"
)

# 80
run_dsep_test(
  ph2 ~ rain2,
  ph2 ~ rain2 + dist,
  alldata,
  "#80 d-separation (ph2 des dist | rain2)"
)

# 81
run_dsep_test(
  ph2 ~ rain2 + dist,
  ph2 ~ rain2 + dist + vel,
  alldata,
  "#81 d-separation (ph2 des vel | rain2, dist)"
)

# 82
run_dsep_test(
  ph2 ~ rain2 + rain3 + dayflow + vel,
  ph2 ~ rain2 + rain3 + dayflow + vel + turb,
  alldata,
  "#82 d-separation (ph2 des turb | rain2, rain3, dayflow, vel)"
)

# 83
run_dsep_test(
  ph2 ~ rain2 + dist,
  ph2 ~ rain2 + dist + res300,
  alldata,
  "#83 d-separation (ph2 des res300 | rain2, dist)"
)

# 84
run_dsep_test(
  ph2 ~ rain2 + rain3 + airtemp + dayflow + res300,
  ph2 ~ rain2 + rain3 + airtemp + dayflow + res300 + temp,
  alldata,
  "#84 d-separation (ph2 des temp | rain2, rain3, airtemp, dayflow, res300)"
)

# 86
tempo <- lm(log(alk) ~ rain2 + rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time, data = alldata)
print_model_block("#86 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain2 + rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  log(alk) ~ rain2 + rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + ph2,
  alldata,
  "#86 d-separation (alk des ph2 | rain2, rain3, time, dalk, dayflow, dist, eff)"
)

# 87
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + rain2,
  log(cond) ~ tp + alk + time + dist:time:eff + rain2 + ph2,
  alldata,
  "#87 d-separation (cond des ph2 | rain2, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 88
run_dsep_test(
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + ph2,
  alldata,
  "#88 d-separation (no3 des ph2 | rain2, rain3, time, dno3, dayflow, dist, eff)"
)

# 89
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + rain2,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#89 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + rain2,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + rain2 + ph2,
  alldata,
  "#89 d-separation (chla des ph2 | rain2, solar, dayflow, dchla, turb, temp, eff, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 90
run_dsep_test(
  log(zn1) ~ rain2 + eff,
  log(zn1) ~ rain2 + eff + ph2,
  alldata,
  "#90 d-separation (zn des ph2 | rain2, eff)"
)

# 91
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + ph2,
  alldata,
  "#91 d-separation (toc des ph2 | rain2, time, rain3, dayflow, dtoc, dist, eff)"
)

# 92
run_capscale_residual_check(
  wangbug.BC ~ ph2 + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla)) + log(zn) + rain2,
  ~ cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + rain2,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "ph2"],
  xlab = "ph2"
)

# 208
run_dsep_test(
  dayflow ~ dflow + rain3,
  dayflow ~ dflow + rain3 + dist,
  alldata,
  "#208 d-separation (dayflow des dist | rain3, dflow)"
)

# 209
tempo <- lm(vel ~ rain3 + dflow + dist, data = alldata)
print_model_block("#209 base vel model", tempo)
run_dsep_test(
  vel ~ rain3 + dflow + dist,
  vel ~ rain3 + dflow + dist + dayflow,
  alldata,
  "#209 d-separation (vel des dayflow | rain3, dflow, dist)"
)

# 210
tempo <- lm(dayflow ~ rain3 + dflow + dist, data = alldata)
print_model_block("#210 base dayflow model", tempo)
run_dsep_test(
  dayflow ~ rain3 + dflow + dist,
  dayflow ~ rain3 + dflow + dist + res300,
  alldata,
  "#210 d-separation (dayflow des res300 | rain3, dflow, dist)"
)

# 212
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff + rain3 + dflow, data = alldata[complete.cases(alldata[, c("cond", "tp")]), ])
print_model_block("#212 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + rain3 + dflow,
  log(cond) ~ tp + alk + time + dist:time:eff + rain3 + dflow + dayflow,
  alldata,
  "#212 d-separation (cond des dayflow | rain3, dflow, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond", "tp")]), ]
)

# 213
tempo <- lm(log(zn) ~ rain3 + dflow + eff, data = alldata)
print_model_block("#213 base zn model", tempo)
run_dsep_test(
  log(zn) ~ rain3 + dflow + eff,
  log(zn) ~ rain3 + dflow + eff + dayflow,
  alldata,
  "#213 d-separation (zn des dayflow | rain3, dflow, eff)"
)

# 214
run_capscale_residual_check(
  wangbug.BC ~ dayflow + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + rain3 + dflow),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + rain3 + log(zn) + rain3 + dflow,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "dayflow"],
  xlab = "dayflow"
)

# 273
tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel), data = alldata)
print_model_block("#273 base turb model", tempo)
run_dsep_test(
  log(turb) ~ log(rain3) + log(dayflow) + log(vel),
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + dist,
  alldata,
  "#273 d-separation (turb des dist | rain3, dayflow, vel)"
)

# 274
tempo <- lm(temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3, data = alldata)
print_model_block("#274 base temp model", tempo)
run_dsep_test(
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3,
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3 + dist,
  alldata,
  "#274 d-separation (temp des dist | airtemp, rain3, dayflow, res300)"
)

# 275
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#275 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist,
  alldata,
  "#275 d-separation (tp des dist | alk, ph2, rain2, rain3, dtp, dayflow, eff)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 276
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#276 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist,
  alldata,
  "#276 d-separation (chla des dist | solar, dayflow, dchla, turb, temp, eff, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 277
tempo <- lm(log(zn) ~ eff, data = alldata)
print_model_block("#277 base zn model", tempo)
run_dsep_test(
  log(zn) ~ eff,
  log(zn) ~ dist + eff + dayflow,
  alldata,
  "#277 d-separation (zn des dist | eff)"
)

# 278
run_capscale_residual_check(
  wangbug.BC ~ dist + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla)) + log(zn),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn),
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "dist"],
  xlab = "dist"
)

# 279
tempo <- lm(vel ~ dist, data = alldata)
print_model_block("#279 base vel model", tempo)
run_dsep_test(
  vel ~ dist,
  vel ~ res300 + dist + dayflow,
  alldata,
  "#279 d-separation (vel des res300 | dist)"
)

# 280
tempo <- lm(temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3 + dist, data = alldata)
print_model_block("#280 base temp model", tempo)
run_dsep_test(
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3 + dist,
  temp ~ rain3 + airtemp + dayflow + res300 + dayflow:res300 + rain3:res300 + dayflow:res300:rain3 + dist + vel,
  alldata,
  "#280 d-separation (vel des temp | airtemp, rain3, dayflow, res300, dist)"
)

# 282
tempo <- lm(log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#282 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist,
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + vel,
  alldata,
  "#282 d-separation (alk des vel | rain3, time, dalk, dayflow, dist, eff)"
)

# 283
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#283 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist + vel,
  alldata,
  "#283 d-separation (tp des vel | dist, time, rain3, rain2, ph2, dayflow, dtp, eff, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 284
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff + dist, data = alldata[complete.cases(alldata[, c("cond")]), ])
print_model_block("#284 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + dist,
  log(cond) ~ tp + alk + time + dist:time:eff + dist + vel,
  alldata,
  "#284 d-separation (cond des vel | alk, tp, dist)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 285
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#285 base no3 model", tempo)
run_dsep_test(
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + vel,
  alldata,
  "#285 d-separation (no3 des vel | rain2, rain3, time, dno3, dayflow, dist, temp, eff)"
)

# 286
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#286 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist + vel,
  alldata,
  "#286 d-separation (chla des vel | solar, dayflow, dchla, turb, temp, eff, tp, no3, dist)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 287
run_dsep_test(
  log(zn1) ~ dist + eff,
  log(zn1) ~ dist + eff + vel,
  alldata,
  "#287 d-separation (zn des vel | dist, eff)"
)

# 288
run_dsep_test(
  log(toc) ~ dayflow + rain2 + rain3 + eff:I(log(dtoc / dayflow)):time + eff:I(log(dtoc / dayflow)):I(dist == 4.08):time + dist,
  log(toc) ~ dayflow + rain2 + rain3 + eff:I(log(dtoc / dayflow)):time + eff:I(log(dtoc / dayflow)):I(dist == 4.08):time + dist + vel,
  alldata,
  "#288 d-separation (toc des vel | rain2, rain3, time, dist, dayflow, dtoc, eff)"
)

# 289
run_capscale_residual_check(
  wangbug.BC ~ vel + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla)) + log(zn) + dist,
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + dist,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "vel"],
  xlab = "vel"
)

# 290
tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel) + dist, data = alldata)
print_model_block("#290 base turb model", tempo)
run_dsep_test(
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + dist,
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + dist + res300,
  alldata,
  "#290 d-separation (turb des res300 | rain3, dayflow, vel, dist)"
)

# 291
tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel) + airtemp + res300, data = alldata)
print_model_block("#291 base turb model", tempo)
run_dsep_test(
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + airtemp + res300,
  log(turb) ~ log(rain3) + log(dayflow) + log(vel) + airtemp + res300 + temp,
  alldata,
  "#291 d-separation (turb des temp | rain3, dayflow, vel, airtemp, res300)"
)

# 293
tempo <- lm(log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + vel, data = alldata)
print_model_block("#293 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + vel,
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + vel + turb,
  alldata,
  "#293 d-separation (alk des turb | rain3, dayflow, vel, time, dalk, dist, eff)"
)

# 294
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(log(dtp / dayflow)):time + eff:I(log(dtp / dayflow)):I(dist == 4.08):time + dist + vel, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#294 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(log(dtp / dayflow)):time + eff:I(log(dtp / dayflow)):I(dist == 4.08):time + dist + vel,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(log(dtp / dayflow)):time + eff:I(log(dtp / dayflow)):I(dist == 4.08):time + dist + vel + turb,
  alldata,
  "#294 d-separation (tp des turb | rain3, dayflow, vel, time, rain2, ph2, dist, dtp, eff, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 295
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff + vel + rain3 + dayflow, data = alldata[complete.cases(alldata[, c("cond")]), ])
print_model_block("#295 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + vel + rain3 + dayflow,
  log(cond) ~ tp + alk + time + dist:time:eff + vel + rain3 + dayflow + turb,
  alldata,
  "#295 d-separation (cond des turb | rain3, dayflow, vel, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 296
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#296 base no3 model", tempo)
run_dsep_test(
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + turb,
  alldata,
  "#296 d-separation (no3 des turb | rain3, dayflow, vel, time, rain2, dno3, dist, temp, eff)"
)

# 297
tempo <- lm(log(zn) ~ rain3 + dayflow + vel + eff, data = alldata)
print_model_block("#297 base zn model", tempo)
run_dsep_test(
  log(zn) ~ rain3 + dayflow + vel + eff,
  log(zn) ~ rain3 + dayflow + vel + eff + turb,
  alldata,
  "#297 d-separation (zn des turb | rain3, dayflow, vel, eff)"
)

# 298
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + vel, data = alldata)
print_model_block("#298 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + vel,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + vel + turb,
  alldata,
  "#298 d-separation (toc des turb | dayflow, rain3, vel, time, rain2, dtoc, dist, eff)"
)

# 299
run_capscale_residual_check(
  wangbug.BC ~ turb + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + rain3 + dayflow + vel),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + rain3 + log(zn) + rain3 + dflow + vel,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "turb"],
  xlab = "turb"
)

# 301
tempo <- lm(log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#301 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist,
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + res300,
  alldata,
  "#301 d-separation (alk des res300 | dist, time, rain, dalk, dayflow, eff)"
)

# 302
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#302 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist + res300,
  alldata,
  "#302 d-separation (tp des res300 | dist, time, dtp, dayflow, rain2, rain3, ph2, eff, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 303
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff + dist, data = alldata[complete.cases(alldata[, c("cond")]), ])
print_model_block("#303 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + dist,
  log(cond) ~ tp + alk + time + dist:time:eff + dist + res300,
  alldata,
  "#303 d-separation (cond des res300 | dist, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 304
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#304 base no3 model", tempo)
run_dsep_test(
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + res300,
  alldata,
  "#304 d-separation (no3 des res300 | dist, time, dno3, rain2, rain3, dayflow, temp, eff)"
)

# 305
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#305 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist + res300,
  alldata,
  "#305 d-separation (chla des res300 | dist, solar, dayflow, dchla, turb, temp, eff, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 306
run_dsep_test(
  log(zn1) ~ dist + eff,
  log(zn1) ~ dist + eff + res300,
  alldata,
  "#306 d-separation (zn des res300 | dist, eff)"
)

# 307
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + res300,
  alldata,
  "#307 d-separation (toc des res300 | dayflow, rain3, time, rain2, dtoc, dist, eff)"
)

# 308
run_capscale_residual_check(
  wangbug.BC ~ res300 + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + dist),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + dist,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "res300"],
  xlab = "res300"
)

# 310
tempo <- lm(log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + res300, data = alldata)
print_model_block("#310 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + res300,
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + res300 + temp,
  alldata,
  "#310 d-separation (alk des temp | airtemp, rain3, dayflow, res300, time, dalk, dist, eff)"
)

# 311
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist + airtemp, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#311 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist + airtemp,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist + airtemp + temp,
  alldata,
  "#311 d-separation (tp des temp | airtemp, rain3, dayflow, res300, time, rain2, ph2, dtp, eff, alk, dist)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 312
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff + airtemp + rain3 + dayflow + res300, data = alldata[complete.cases(alldata[, c("cond")]), ])
print_model_block("#312 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff + airtemp + rain3 + dayflow + res300,
  log(cond) ~ tp + alk + time + dist:time:eff + airtemp + rain3 + dayflow + res300 + temp,
  alldata,
  "#312 d-separation (cond des temp | airtemp, rain3, dayflow, res300, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 313
tempo <- lm(log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff, data = alldata)
print_model_block("#313 base zn model", tempo)
run_dsep_test(
  log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff,
  log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff + temp,
  alldata,
  "#313 d-separation (zn des temp | airtemp, rain3, dayflow, res300, dist, eff)"
)

# 314
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + airtemp + res300, data = alldata)
print_model_block("#314 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + airtemp + res300,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + airtemp + res300 + temp,
  alldata,
  "#314 d-separation (toc des temp | airtemp, dayflow, rain3, res300, time, rain2, dtoc, dist, eff)"
)

# 316
tempo <- lm(log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + dno3 + temp + rain2, data = alldata)
print_model_block("#316 base alk model", tempo)
run_dsep_test(
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + dno3 + temp + rain2,
  log(alk) ~ rain3 + dayflow + eff:I(log(dalk / dayflow)):time + eff:I(log(dalk / dayflow)):I(dist == 4.08):time + dist + dno3 + temp + rain2 + no3,
  alldata,
  "#316 d-separation (alk des no3 | time, dno3, dayflow, dist, temp, eff, rain3, rain2, dalk)"
)

# 317
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist + time + rain3 + dalk,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#317 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist + time + rain3 + dalk,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + dist + time + rain3 + dalk + alk,
  alldata,
  "#317 d-separation (chla des alk | time, rain3, dalk, dayflow, dist, eff, solar, dchla, turb, temp, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 318
tempo <- lm(log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff, data = alldata)
print_model_block("#318 base zn model", tempo)
run_dsep_test(
  log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff,
  log(zn1) ~ airtemp + rain3 + dayflow + res300 + dist + eff + alk,
  alldata,
  "#318 d-separation (zn des alk | time, rain3, dalk, dayflow, dist, eff)"
)

# 319
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dalk, data = alldata)
print_model_block("#319 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dalk,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dalk + alk,
  alldata,
  "#319 d-separation (toc des alk | rain3, rain2, time, dalk, dayflow, dtoc, dist, eff)"
)

# 320
run_capscale_residual_check(
  wangbug.BC ~ alk + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + eff + dist + dayflow + dalk + rain3 + log(zn)),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + eff + dist + dayflow + dalk + rain3 + log(zn),
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "alk"],
  xlab = "alk"
)

# 321
tempo <- lm(log(no3) ~ temp + ph2 + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + dtp + alk, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#321 base no3 model", tempo)
run_dsep_test(
  log(no3) ~ temp + ph2 + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + dtp + alk,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + dist + dtp + alk + tp,
  alldata,
  "#321 d-separation (no3 des tp | time, rain3, rain2, ph2, dno3, dayflow, dist, temp, eff, dtp, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 322
tempo <- lm(tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#322 base tp model", tempo)
run_dsep_test(
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(dtp / dayflow):time + eff:I(dtp / dayflow):I(dist == 4.08):time + dist,
  tp ~ alk + ph2 + rain2 + rain3 + eff:I(log(dtp / dayflow)):time + eff:I(log(dtp / dayflow)):I(dist == 4.08):time + dist + zn,
  alldata,
  "#322 d-separation (tp des zn | time, rain3, rain2, ph2, dayflow, dtp, dist, eff, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 323
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + ph2 + alk, data = alldata[complete.cases(alldata[, c("tp")]), ])
print_model_block("#323 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + ph2 + alk,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + ph2 + alk + tp,
  alldata,
  "#323 d-separation (toc des tp | time, rain3, rain2, ph2, dayflow, dtoc, dist, eff, dtp, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp")]), ]
)

# 324
run_capscale_residual_check(
  wangbug.BC ~ tp + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn1) + time + rain3 + rain2 + ph2 + dtp + dayflow + eff + alk + dist),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn1) + time + rain3 + rain2 + ph2 + dtp + dayflow + eff + alk + dist,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "tp"],
  xlab = "tp"
)

# 325
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + alk + tp, data = alldata[complete.cases(alldata[, c("tp", "cond")]), ])
print_model_block("#325 base no3 model", tempo)
run_dsep_test(
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + alk + tp,
  log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3) / dayflow):time + eff:I(log(dno3) / dayflow):I(dist == 4.08):time + alk + tp + cond,
  alldata,
  "#325 d-separation (no3 des cond | time, rain3, rain2, dno3, dayflow, dist, temp, eff, tp, alk)",
  subset_expr = function(d) d[complete.cases(d[, c("tp", "cond")]), ]
)

# 326
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + rain3 + alk,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp", "cond")]), ]
)
print_model_block("#326 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + rain3 + alk,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + rain3 + alk + cond,
  alldata,
  "#326 d-separation (chla des cond | time, dayflow, eff, alk, tp, solar, dchla, turb, temp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp", "cond")]), ]
)

# 327
tempo <- lm(log(cond) ~ tp + alk + time + dist:time:eff, data = alldata[complete.cases(alldata[, c("cond")]), ])
print_model_block("#327 base cond model", tempo)
run_dsep_test(
  log(cond) ~ tp + alk + time + dist:time:eff,
  log(cond) ~ tp + alk + time + dist:time:eff + zn,
  alldata,
  "#327 d-separation (cond des zn | eff, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("cond")]), ]
)

# 328
tempo <- lm(log(toc) ~ dayflow + rain2 + rain3 + eff:I(log(dtoc / dayflow)):time + eff:I(log(dtoc / dayflow)):I(dist == 4.08):time + tp + alk, data = alldata[complete.cases(alldata[, c("tp", "cond")]), ])
print_model_block("#328 base toc model", tempo)
run_dsep_test(
  log(toc) ~ dayflow + rain2 + rain3 + eff:I(log(dtoc / dayflow)):time + eff:I(log(dtoc / dayflow)):I(dist == 4.08):time + tp + alk,
  log(toc) ~ dayflow + rain2 + rain3 + eff:I(log(dtoc / dayflow)):time + eff:I(log(dtoc / dayflow)):I(dist == 4.08):time + tp + alk + cond,
  alldata,
  "#328 d-separation (toc des cond | time, rain3, rain2, dayflow, dtoc, dist, eff, alk, tp)",
  subset_expr = function(d) d[complete.cases(d[, c("tp", "cond")]), ]
)

# 329
tempo <- lm(log(zn1) ~ dist + eff + temp + time + dno3 + rain3 + rain2 + dayflow, data = alldata)
print_model_block("#329 base zn model", tempo)
run_dsep_test(
  log(zn1) ~ dist + eff + temp + time + dno3 + rain3 + rain2 + dayflow,
  log(zn1) ~ dist + eff + temp + time + dno3 + rain3 + rain2 + dayflow + no3,
  alldata,
  "#329 d-separation (zn des no3 | time, dno3, rain3, rain2, dayflow, dist, temp, eff)"
)

# 330
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dno3 + temp + no3, data = alldata)
print_model_block("#330 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dno3 + temp + no3,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + dno3 + temp,
  alldata,
  "#330 d-separation (toc des no3 | time, dno3, rain3, rain2, dayflow, dist, temp, eff, dtoc)"
)

# 331
run_capscale_residual_check(
  wangbug.BC ~ no3 + Condition(temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + dist + time + dno3 + dayflow + eff + rain3 + rain2),
  ~ temp + cond + temp:cond + toc + temp:toc + log(chla) + log(zn) + dist + time + dno3 + dayflow + eff + rain2 + rain3,
  x = cbind(log(alldata$chla), alldata)[complete.cases(cbind(log(alldata$chla), alldata)[, c("temp", "cond", "toc", "chla")]), "no3"],
  xlab = "no3"
)

# 332
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#332 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + zn,
  alldata,
  "#332 d-separation (chla des zn | eff, solar, dayflow, dchla, turb, temp, tp, no3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 333
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + dtoc + dist + rain2 + rain3,
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)
print_model_block("#333 base chla model", tempo)
run_dsep_test(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + dtoc + dist + rain2 + rain3,
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)) + time + dtoc + dist + rain2 + rain3 + toc,
  alldata,
  "#333 d-separation (chla des toc | time, dayflow, dtoc, dist, eff, solar, dchla, turb, temp, no3, tp, rain2, rain3)",
  subset_expr = function(d) d[complete.cases(d[, c("temp", "no3", "tp")]), ]
)

# 334
tempo <- lm(toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist, data = alldata)
print_model_block("#334 base toc model", tempo)
run_dsep_test(
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist,
  toc ~ dayflow + rain2 + rain3 + eff:I((dtoc / dayflow)):time + eff:I((dtoc / dayflow)):I(dist == 4.08):time + dist + zn,
  alldata,
  "#334 d-separation (toc des zn | dist, eff, time, dayflow, dtoc, rain2, rain3)"
)

# ----------------------------------------------------------------------
# 8. Fisher's C test
# ----------------------------------------------------------------------
C <- -2 * (
  log(0.2862) + log(0.2862) + log(0.281) + log(0.4666) + log(0.2584) +
    log(0.729) + log(0.9811) + log(0.9801) + log(0.8075) + log(0.08321) +
    log(0.768) + log(0.7981) + log(0.3135) + log(0.4091) + log(0.3185) +
    log(0.4001) + log(0.8362) + log(0.9707) + log(0.8159) + log(0.3239) +
    log(0.9299) + log(0.3545) + log(0.868) + log(0.1587) + log(0.1848) +
    log(0.859) + log(0.3215) + log(0.9016) + log(0.2587) + log(0.6589) +
    log(0.9717) + log(0.2118) + log(0.2745) + log(0.2654) + log(0.5033) +
    log(0.5644) + log(0.3253) + log(0.4007) + log(0.3457) + log(0.2598) +
    log(0.7497) + log(0.07954) + log(0.5078) + log(0.1927) + log(0.01745) +
    log(0.9612) + log(0.001129) + log(0.768) + log(0.197) + log(0.7472) +
    log(0.372) + log(0.4781) + log(0.04157) + log(0.4621) + log(0.6556) +
    log(0.9284) + log(0.3132) + log(0.2834) + log(0.3957) + log(0.7677) +
    log(0.4636) + log(0.1166) + log(0.159) + log(0.4545) + log(0.3448) +
    log(0.8076) + log(0.2935) + log(0.665) + log(0.2502) + log(0.9427) +
    log(0.3342) + log(0.365) + log(0.1427) + log(0.06621) + log(0.5555) +
    log(0.03393) + log(0.97) + log(0.1328) + log(0.5163) + log(0.9455)
)

cat("\n============================================================\n")
cat("Fisher's C statistic\n")
cat("============================================================\n")
print(C)
print(pchisq(C, df = 2 * 80, lower.tail = FALSE))

cat("\nScript 19 completed successfully.\n")
