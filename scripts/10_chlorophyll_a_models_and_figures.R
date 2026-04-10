# 10_chlorophyll_a_models_and_figures.R
# Author: Rezvan Hatami
# Date: 2026-03-15

rm(list = ls())

# ----------------------------------------------------------------------
# 1. Package setup
# ----------------------------------------------------------------------
required_pkgs <- c(
  "fs",
  "here",
  "lattice",
  "latticeExtra",
  "car"
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
options(scipen = 1)
set.seed(123)

# ----------------------------------------------------------------------
# 2. Paths
# ----------------------------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figures_dir <- fs::path(project_dir, "figures")
script_figure_dir <- fs::path(figures_dir, "10")

fs::dir_create(figures_dir)
fs::dir_create(script_figure_dir)

# ----------------------------------------------------------------------
# 3. Read inputs
# ----------------------------------------------------------------------
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
bugenv_file <- fs::path(output_dir, "02_bugenv_with_pco.rds")
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")

if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}
if (!fs::file_exists(bugenv_file)) {
  stop("Missing input file: ", bugenv_file)
}
if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}

wangenv <- readRDS(wangenv_file)
bugenv <- readRDS(bugenv_file)

script3_inputs <- readRDS(script3_inputs_file)
if (!is.list(script3_inputs) || !"alldata" %in% names(script3_inputs)) {
  stop("03_environmental_script_inputs.rds does not contain alldata.")
}
alldata <- as.data.frame(script3_inputs$alldata)

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

panel.cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(invisible(NULL))
  }
  
  if (sd(x) == 0 || sd(y) == 0) {
    text(0.5, 0.5, "sd = 0")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x, y, method = "spearman"))
  txt <- paste0("r= ", formatC(r, digits = digits, format = "f"))
  text(0.5, 0.8, txt)
  
  p <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  txt2 <- if (is.na(p)) {
    "p= NA"
  } else if (p < 0.01) {
    "p= <0.01"
  } else {
    paste0("p= ", formatC(p, digits = digits, format = "f"))
  }
  text(0.5, 0.3, txt2)
  
  invisible(NULL)
}

make_pairs_plot <- function(formula, data, upper.panel = panel.cor, pch = 20, ...) {
  pairs(
    formula,
    data = data,
    upper.panel = upper.panel,
    pch = pch,
    na.action = na.omit,
    ...
  )
}

print_model_block <- function(model_name, model_obj) {
  cat("\n")
  cat("============================================================\n")
  cat(model_name, "\n")
  cat("============================================================\n")
  print(summary(model_obj))
  print(anova(model_obj))
}

# ----------------------------------------------------------------------
# 5. Prepare analysis objects
# ----------------------------------------------------------------------
wangenv <- prepare_numeric_columns(as.data.frame(wangenv))
bugenv <- prepare_numeric_columns(as.data.frame(bugenv))
alldata <- prepare_numeric_columns(as.data.frame(alldata))

if (!"time" %in% names(alldata) && "day" %in% names(alldata)) {
  unique_days <- sort(unique(stats::na.omit(alldata$day)))
  alldata$time <- match(alldata$day, unique_days)
}

if (!"time" %in% names(bugenv) && "day" %in% names(bugenv)) {
  unique_days <- sort(unique(stats::na.omit(bugenv$day)))
  bugenv$time <- match(bugenv$day, unique_days)
}

# ----------------------------------------------------------------------
# 6. Testing other variables for causal diagram
# ----------------------------------------------------------------------

# Chlorophyll A exploratory pair plots
jpeg(
  filename = fs::path(script_figure_dir, "10_chla_scatterplot_checks.jpg"),
  width = 14,
  height = 14,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800
)

par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  dchla ~ eff + solar + temp + no3 + tp + turb + don + dtoc + canop +
    toc + log(chla) + chla,
  data = alldata
)

make_pairs_plot(
  log(I(dchla/dayflow)) ~ eff + solar + temp + no3 + tp + turb + don +
    dtoc + canop + toc + log(chla) + chla,
  data = alldata[alldata$eff == 1, ],
  upper.panel = NULL
)

pairs_data3 <- alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
if (exists("tempo") && length(residuals(tempo)) == nrow(pairs_data3)) {
  pairs_data3$tempo_residuals <- residuals(tempo)
  make_pairs_plot(
    ~ eff + solar + temp + no3 + tn + tp + log(turb) + don + dtoc +
      canop + toc + vel + I(log(dchla/dayflow)) + dayflow + log(chla) +
      tempo_residuals,
    data = pairs_data3
  )
}

dev.off()

# Tested chlorophyll A models
chla_formula_list_1_35 <- list(
  tempo.model1  = log(chla) ~ temp,
  tempo.model2  = log(chla) ~ on,
  tempo.model3  = log(chla) ~ eff,
  tempo.model4  = log(chla) ~ on + on:temp,
  tempo.model5  = log(chla) ~ no3,
  tempo.model6  = log(chla) ~ no3 + temp:no3,
  tempo.model7  = log(chla) ~ vel,
  tempo.model8  = log(chla) ~ vel + vel:temp,
  tempo.model9  = log(chla) ~ nh3,
  tempo.model10 = log(chla) ~ nh3 + temp:nh3,
  tempo.model11 = log(chla) ~ no2,
  tempo.model12 = log(chla) ~ no2 + no2:temp,
  tempo.model13 = log(chla) ~ tp,
  tempo.model14 = log(chla) ~ tp + tp:temp,
  tempo.model15 = log(chla) ~ op,
  tempo.model16 = log(chla) ~ op + op:temp,
  tempo.model17 = log(chla) ~ solar,
  tempo.model18 = log(chla) ~ solar + solar:temp,
  tempo.model19 = log(chla) ~ turb,
  tempo.model20 = log(chla) ~ turb + turb:temp,
  tempo.model21 = log(chla) ~ dchla,
  tempo.model22 = log(chla) ~ canop,
  tempo.model23 = log(chla) ~ canop + canop:temp,
  tempo.model24 = log(chla) ~ dayflow,
  tempo.model25 = log(chla) ~ don,
  tempo.model26 = log(chla) ~ dtoc,
  tempo.model27 = log(chla) ~ cfpom,
  tempo.model28 = log(chla) ~ toc,
  tempo.model29 = log(chla) ~ dtoc + dtoc:temp,
  tempo.model30 = log(chla) ~ don,
  tempo.model31 = log(chla) ~ don + don:temp,
  tempo.model32 = log(chla) ~ tn,
  tempo.model33 = log(chla) ~ tn + dtn:temp,
  tempo.model34 = log(chla) ~ cond,
  tempo.model35 = log(chla) ~ cond + cond:temp
)

tempo_candidate_models_1_35 <- lapply(chla_formula_list_1_35, function(fm) {
  lm(formula = fm, data = alldata)
})

for (model_name in names(tempo_candidate_models_1_35)) {
  print_model_block(model_name, tempo_candidate_models_1_35[[model_name]])
}

tempo <- tempo_candidate_models_1_35[["tempo.model35"]]

png(
  filename = fs::path(script_figure_dir, "10_chla_temp_plot.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(alldata$chla ~ alldata$temp)
dev.off()

cat("\n============================================================\n")
cat("cor.test(alldata$chla, alldata$temp)\n")
cat("============================================================\n")
print(cor.test(alldata$chla, alldata$temp))

if ("cod" %in% names(alldata)) {
  png(
    filename = fs::path(script_figure_dir, "10_logchla_cod_plot.png"),
    width = 2400,
    height = 1800,
    res = 300
  )
  plot(log(alldata$chla), alldata$cod)
  dev.off()
  
  cat("\n============================================================\n")
  cat("cor.test(log(alldata$chla), alldata$cod)\n")
  cat("============================================================\n")
  print(cor.test(log(alldata$chla), alldata$cod))
}

# Chlorophyll A plotted against spatial position
png(
  filename = fs::path(script_figure_dir, "10_chla_spatial_position.png"),
  width = 2400,
  height = 1800,
  res = 300
)

plot(
  bugenv[bugenv$day == 1, "dist"],
  bugenv[bugenv$day == 1, "chla"],
  type = "b",
  ylim = c(0, 40),
  pch = 19,
  xlab = "Distance(km)",
  ylab = "Chlorophyll A (mg/L)",
  main = "Chlorophyll A plotted against spatial position",
  cex = 1,
  lwd = 3,
  lty = 1
)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "chla"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "chla"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "chla"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "chla"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"), lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2, col = c("black", "gold2", "blue", "green3", "red"), ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()

# Multiple regression analysis for chlorophyll A as a function of time, distance and eff
tempo <- lm(
  log(chla) ~ dist + time + eff + eff:time + eff:dist + eff:time:dist +
    eff:time:I(dist^2),
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)

cat("\n============================================================\n")
cat("Spatiotemporal chlorophyll A model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(
  cbind(exp(fitted(tempo)), residuals(tempo)),
  alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)



print(xyplot(log(chla) ~ dist, groups = time, data = alldata, type = "b", lwd = 2, col = c("black", "gold2", "blue", "green3", "red")))
print(xyplot(dayflow ~ dist, groups = time, data = alldata, type = "b", lwd = 2, col = c("black", "gold2", "blue", "green3", "red")))
print(xyplot(turb ~ dist, groups = time, data = alldata, type = "b", lwd = 2, col = c("black", "gold2", "blue", "green3", "red")))

a <- xyplot(chla ~ dist, group = time, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], type = "b", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
b <- xyplot(exp(fitted(tempo)) ~ dist, group = time, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], type = "l", lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
c <- xyplot(toc ~ dist, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], panel = function(x, y) { panel.abline(v = 4, lty = 2) })

png(
  filename = fs::path(script_figure_dir, "10_chla_xyplot_overlay.png"),
  width = 2400,
  height = 1800,
  res = 300
)
print(a + as.layer(b) + as.layer(c))
dev.off()

# Chlorophyll A linear model with other environmental variables
tempo.model36 <- lm(log(chla) ~ solar + temp + tp + temp:tp + dchla + no3 + dtoc + turb, data = wangenv)
tempo.model37 <- lm(log(chla) ~ solar + temp2 + tp + no3 + dayflow + eff + eff:dist + eff:I(dist^2) + eff:dchla + eff:dchla:dist + eff:dchla:I(dist^2) + eff:dayflow + eff:dayflow:dist + eff:dayflow:I(dist^2) + eff:temp2 + eff:temp2:dist + eff:temp2:I(dist^2), data = alldata)
tempo.model38 <- lm(log(chla) ~ solar + temp2 + tp + no3 + dayflow + eff + eff:dchla + eff:temp2 + eff:temp2:dist, data = alldata)
tempo.model39 <- lm(log(chla) ~ solar + temp2 + tp + no3 + dayflow + eff:dchla, data = alldata)
tempo.model40 <- lm(log(chla) ~ solar + temp2 + tp + no3 + turb + dayflow, data = alldata)
tempo.model41 <- lm(log(chla) ~ solar + temp + tp + dchla + no3 + dtoc + turb + dist + temp:tp + dchla:eff, data = alldata)
tempo.model42 <- lm(log(chla) ~ solar + temp2 + tp + no3 + dtoc + turb + dayflow + eff:dist + eff:I(dchla/dayflow) + eff:I(dchla/dayflow):dist + eff:dist:temp2, data = alldata)
tempo.model43 <- lm(log(chla) ~ solar + temp + tp + no3 + turb + dayflow + solar:temp + eff:log(I(dchla/dayflow)), data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ])
tempo.model44 <- lm(log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla/dayflow)), data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ])

tempo_candidate_models_36_44 <- list(
  tempo.model36 = tempo.model36,
  tempo.model37 = tempo.model37,
  tempo.model38 = tempo.model38,
  tempo.model39 = tempo.model39,
  tempo.model40 = tempo.model40,
  tempo.model41 = tempo.model41,
  tempo.model42 = tempo.model42,
  tempo.model43 = tempo.model43,
  tempo.model44 = tempo.model44
)

for (model_name in names(tempo_candidate_models_36_44)) {
  print_model_block(model_name, tempo_candidate_models_36_44[[model_name]])
}

tempo <- tempo.model44
tempo1 <- cbind(
  cbind(exp(fitted(tempo)), exp(residuals(tempo))),
  alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)

attributes(tempo1)
names(tempo1[, 1:2]) <- c("fitted", "residual")
tempo1[, "2"]

cat("\n============================================================\n")
cat("Final chlorophyll A environmental model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

png(
  filename = fs::path(script_figure_dir, "10_chla_environmental_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)

plot(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "l", ylim = c(0, 40), pch = 15, col = "gold2", xlab = "Distance(km)", ylab = "Chlorophyll A (mg/L)", main = "Prediction with environmental variables", lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "l", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "l", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "l", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"), lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2, col = c("black", "gold2", "blue", "green3", "red"), ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()

tempo3 <- lm(tempo1[, 2] ~ tempo1$pred1)
cat("\n============================================================\n")
cat("tempo3\n")
cat("============================================================\n")
print(summary(tempo3))

cat("\n============================================================\n")
cat("summary(tempo)\n")
cat("============================================================\n")
print(summary(tempo))

cat("\n============================================================\n")
cat("anova(tempo)\n")
cat("============================================================\n")
print(anova(tempo))

cat("\n============================================================\n")
cat("sqrt(vif(tempo))\n")
cat("============================================================\n")
print(sqrt(vif(tempo)))

png(
  filename = fs::path(script_figure_dir, "10_chla_environmental_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "10_chla_environmental_avplots.png"),
  width = 2400,
  height = 1800,
  res = 300
)
avPlots(tempo)
dev.off()

tempo.check <- lm(
  residuals(tempo) ~ dist + time + eff + eff:time + eff:dist + eff:time:dist + eff:time:I(dist^2),
  data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
)

cat("\n============================================================\n")
cat("tempo.check\n")
cat("============================================================\n")
print(summary(tempo.check))
print(anova(tempo.check))

print(xyplot(residuals(tempo) ~ dist, groups = time, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], type = "b", pch = 19, lwd = 2, col = c("black", "gold2", "blue", "green3", "red")))

a <- xyplot(chla ~ dist, group = time, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], type = "b", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
b <- xyplot(exp(fitted(tempo)) ~ dist, group = time, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], type = "l", lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
c <- xyplot(toc ~ dist, data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ], panel = function(x, y) { panel.abline(v = 4, lty = 2) })

png(
  filename = fs::path(script_figure_dir, "10_chla_observed_fitted_toc_overlay.png"),
  width = 2400,
  height = 1800,
  res = 300
)
print(a + as.layer(b) + as.layer(c))
dev.off()


##  Calling data in to make sure the graph is using the right data
script3_inputs <- readRDS(fs::path(output_dir, "03_environmental_script_inputs.rds"))
alldata <- script3_inputs$alldata

alldata_chla_cc <- alldata[
  complete.cases(alldata[, c("temp", "no3", "tp")]),
]

tiff(
  file = fs::path(script_figure_dir, "10_chlorophyll_a_figure.tif"),
  width = 10, height = 14, units = "in",
  pointsize = 12, bg = "transparent",
  res = 800, compression = "lzw"
)

par(mfrow = c(3, 1), mar = c(4.5, 4.5, 2.5, 3), cex = 1.5, cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8)

# a) Chlorophyll A plotted against spatial position
plot(
  bugenv[bugenv$day == 1, "dist"],
  bugenv[bugenv$day == 1, "chla"],
  type = "b",
  ylim = c(0, 40),
  pch = 19,
  xlab = "",
  ylab = "Chlorophyll A (mg/L)",
  main = "a) Chlorophyll A plotted against spatial position",
  lwd = 3,
  lty = 1
)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "chla"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "chla"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "chla"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "chla"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0),
       legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

# b) Chlorophyll A with spatiotemporal prediction
tempo <- lm(
  log(chla) ~ dist + time + eff + eff:time + eff:dist + eff:time:dist + eff:time:I(dist^2),
  data = alldata_chla_cc
)
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(
  cbind(exp(fitted(tempo)), residuals(tempo)),
  alldata_chla_cc
)

plot(
  tempo1[tempo1$time == 1, "dist"],
  tempo1[tempo1$time == 1, "chla"],
  type = "p",
  pch = 19,
  col = "black",
  xlab = "",
  ylab = "Chlorophyll A (mg/L)",
  main = "b) Chlorophyll A with spatiotemporal prediction",
  ylim = c(0, 40),
  cex = 1,
  lwd = 3,
  lty = 1
)
points(tempo1[tempo1$time == 2, "dist"], tempo1[tempo1$time == 2, "chla"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1[tempo1$time == 3, "dist"], tempo1[tempo1$time == 3, "chla"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$time == 4, "dist"], tempo1[tempo1$time == 4, "chla"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$time == 5, "dist"], tempo1[tempo1$time == 5, "chla"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)

points(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", col = "black", pch = 19, cex = 1.5, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)

legend("topleft", inset = c(0, 0),
       legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

# c) Prediction with environmental variables
tempo <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla/dayflow)),
  data = alldata_chla_cc
)
print(summary(tempo))
print(anova(tempo))

tempo1 <- cbind(
  cbind(exp(fitted(tempo)), exp(residuals(tempo))),
  alldata_chla_cc
)
names(tempo1[, 1:2]) <- c("fitted", "residual")
tempo1[, "2"]

plot(
  tempo1[tempo1$day == 126, "dist"],
  tempo1[tempo1$day == 126, "1"],
  type = "b",
  ylim = c(0, 40),
  pch = 15,
  col = "gold2",
  xlab = "Distance(km)",
  ylab = "Chlorophyll A (mg/L)",
  main = "c) Prediction with environmental variables",
  lwd = 3,
  lty = 2
)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)

legend("topleft", inset = c(0, 0),
       legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()


# Chlorophyll A non-linear model
temp3 <- alldata$temp2 - 20
chlanlm <- 0.8 + 1.5 * (alldata$time == 4) + 2.8 * alldata$eff * (alldata$dchla / alldata$dayflow) * exp(-0.03 * 0.9^temp3 * alldata$dist / alldata$vel2)
chlanlm <- 0.8 + 1.5 * (alldata$time == 2) + 0.5 * (alldata$time == 3) + 1.5 * (alldata$time == 4) + 2.8 * alldata$eff * (alldata$dchla / alldata$dayflow) * exp(-0.03 * 1.1^temp3 * alldata$dist / alldata$vel2)
alldata2 <- cbind(alldata, chlanlm, temp3)

png(
  filename = fs::path(script_figure_dir, "10_chla_nonlinear_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)

par(mar = c(5, 5, 4, 2), cex = 0.9)
plot(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 1, "chla"], type = "p", pch = 19, col = "black", xlab = "Distance (km)", ylab = "Chla (mg/L)", main = "Chla with spatiotemporal prediction", ylim = c(0, 40), cex = 1, lwd = 3, lty = 1)
points(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 2, "chla"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 3, "chla"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 4, "chla"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 5, "chla"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(alldata2[alldata2$day == 1, "dist"], alldata2[alldata2$day == 1, "chlanlm"], type = "l", col = "black", pch = 15, lwd = 3, cex = 1.5, lty = 1)
points(alldata2[alldata2$day == 126, "dist"], alldata2[alldata2$day == 126, "chlanlm"], type = "l", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(alldata2[alldata2$day == 260, "dist"], alldata2[alldata2$day == 260, "chlanlm"], type = "l", col = "blue", pch = 8, cex = 1.5, lty = 3)
points(alldata2[alldata2$day == 336, "dist"], alldata2[alldata2$day == 336, "chlanlm"], type = "l", col = "green3", pch = 17, cex = 1.5, lty = 4)
points(alldata2[alldata2$day == 518, "dist"], alldata2[alldata2$day == 518, "chlanlm"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"), lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2, col = c("black", "gold2", "blue", "green3", "red"), ncol = 2, horiz = FALSE, cex = 0.9, title = "months")
abline(v = 4, lty = 2)

dev.off()

if (requireNamespace("caret", quietly = TRUE)) {
  lmfit <- caret::train(
    log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla/dayflow)),
    data = alldata[complete.cases(alldata[, c("temp", "no3", "tp")]), ]
  )
  
  cat("\n============================================================\n")
  cat("caret cross-validation model\n")
  cat("============================================================\n")
  print(lmfit)
}

cat("\nScript 10 completed successfully.\n")

