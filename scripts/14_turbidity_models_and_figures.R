# 14_turbidity_models_and_figures.R
# Author: Rezvan Hatami
# Date: 19 March 2026

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
script_figure_dir <- fs::path(figures_dir, "14")

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
  
  r <- cor(x, y, method = "spearman")
  txt <- paste("r= ", format(c(r, 0.123456789), digits = digits)[1], sep = "")
  text(0.5, 0.8, txt)
  
  p <- cor.test(x, y, method = "spearman")$p.value
  txt2 <- paste("p= ", format(c(p, 0.123456789), digits = digits)[1], sep = "")
  if (p < 0.01) {
    txt2 <- "p= <0.01"
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

if (!is.factor(alldata$time)) {
  alldata$time <- as.factor(alldata$time)
}

# ----------------------------------------------------------------------
# 6. Turbidity exploratory checks
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ time + dayflow + weekflow + monflow + eff + dturb + vel + rain2 + solar +
    log(rain2) + veg30m + res300 + depth + log(cond) + turb + log(turb),
  data = alldata
)

tempo <- lm(turb ~ rain3, data = alldata)
print_model_block("turb_rain3", tempo)

tempo <- lm(turb ~ vel, data = alldata)
print_model_block("turb_vel", tempo)

tempo <- lm(turb ~ dayflow, data = alldata)
print_model_block("turb_dayflow", tempo)

tempo <- lm(turb ~ eff, data = alldata)
print_model_block("turb_eff", tempo)

tempo <- lm(turb ~ dturb, data = alldata)
print_model_block("turb_dturb", tempo)

if ("veg30m" %in% names(alldata)) {
  tempo <- lm(turb ~ veg30m, data = alldata)
  print_model_block("turb_veg30m", tempo)
}

if ("depth" %in% names(alldata)) {
  tempo <- lm(turb ~ depth, data = alldata)
  print_model_block("turb_depth", tempo)
}

png(
  filename = fs::path(script_figure_dir, "14_turb_spatial_position.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(bugenv[bugenv$day == 1, "dist"], bugenv[bugenv$day == 1, "turb"], type = "b", ylim = c(0, 100), pch = 19, xlab = "Distance(km)", ylab = "Turbidity (NTU)", main = "Turbidity plotted against spatial position", lwd = 3, lty = 1)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "turb"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "turb"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "turb"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "turb"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)
dev.off()

# ----------------------------------------------------------------------
# 7. Turbidity spatiotemporal model
# ----------------------------------------------------------------------
tempo <- lm(log(turb) ~ dist + time + dist:time, data = alldata)

cat("\n============================================================\n")
cat("Turbidity spatiotemporal model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "14_turb_spatiotemporal_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$time == 1, "dist"], tempo1[tempo1$time == 1, "turb"], type = "p", pch = 19, col = "black", xlab = "Distance (km)", ylab = "Turbidity (NTU)", main = "Turbidity with spatiotemporal prediction", ylim = c(0, 100), cex = 1, lwd = 3, lty = 1)
points(tempo1[tempo1$time == 2, "dist"], tempo1[tempo1$time == 2, "turb"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1[tempo1$time == 3, "dist"], tempo1[tempo1$time == 3, "turb"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$time == 4, "dist"], tempo1[tempo1$time == 4, "turb"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$time == 5, "dist"], tempo1[tempo1$time == 5, "turb"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", col = "black", pch = 19, cex = 1.2, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)
dev.off()

cat("\n============================================================\n")
cat("Repeated turbidity spatiotemporal model output\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

png(
  filename = fs::path(script_figure_dir, "14_turb_spatiotemporal_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "14_turb_spatiotemporal_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "14_turb_spatiotemporal_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

par(mar = c(5, 5, 4, 2), cex = 0.9)

png(
  filename = fs::path(script_figure_dir, "14_turb_spatiotemporal_residuals_distance.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "2"], type = "b", ylim = c(-0.3, 0.4), pch = 19, cex = 1.5, xlab = "Distance(km)", ylab = "residuals", main = "prediction with space and time", lwd = 3, lty = 4)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "2"], type = "b", col = "gold2", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "2"], type = "b", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "2"], type = "b", col = "green3", pch = 17, cex = 1.5, lwd = 3)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "2"], type = "b", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
       lty = 1, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 5, horiz = FALSE, cex = 0.7, title = "months")
abline(v = 4, lty = 2)
dev.off()

# ----------------------------------------------------------------------
# 8. Turbidity environmental model
# ----------------------------------------------------------------------
tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel), data = alldata)

cat("\n============================================================\n")
cat("Turbidity environmental model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))
print(sqrt(vif(tempo)))

tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "14_turb_environmental_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", ylim = c(0, 100), pch = 19, col = "black", xlab = "Distance(km)", ylab = "Turbidity (NTU)", main = "Prediction with environmental parameters", lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "l", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "l", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "l", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "l", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)
dev.off()

png(
  filename = fs::path(script_figure_dir, "14_turb_environmental_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

tempo.check <- lm(residuals(tempo) ~ dist + time + dist:time, data = alldata)

cat("\n============================================================\n")
cat("Turbidity residual spatiotemporal check\n")
cat("============================================================\n")
print(summary(tempo.check))
print(anova(tempo.check))

par(mar = c(5, 5, 4, 2), cex = 0.9)

png(
  filename = fs::path(script_figure_dir, "14_turb_environmental_residuals_distance.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "2"], type = "b", ylim = c(-0.3, 0.4), pch = 19, cex = 1.5, xlab = "Distance(km)", ylab = "residuals", main = "prediction with space and time", lwd = 3, lty = 4)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "2"], type = "b", col = "gold2", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "2"], type = "b", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "2"], type = "b", col = "green3", pch = 17, cex = 1.5, lwd = 3)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "2"], type = "b", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
       lty = 1, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 5, horiz = FALSE, cex = 0.7, title = "months")
abline(v = 4, lty = 2)
dev.off()

tempo3 <- lm(tempo1[, 2] ~ tempo1$pred1)
cat("\n============================================================\n")
cat("tempo3\n")
cat("============================================================\n")
print(summary(tempo3))

tempo <- lm(tempo1[, 2] ~ dist + time + dist:time, data = alldata)

cat("\n============================================================\n")
cat("Residuals against time and space\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "14_turb_residual_time_space_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

# ----------------------------------------------------------------------
# 9. Turbidity figure
# ----------------------------------------------------------------------
tiff(
  file = fs::path(script_figure_dir, "14_turbidity_figure.tif"),
  width = 10,
  height = 14,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800,
  compression = "lzw"
)

par(mfrow = c(3, 1), mar = c(3, 4.5, 2.5, 3.5), cex = 1.5, cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8)

plot(bugenv[bugenv$day == 1, "dist"], bugenv[bugenv$day == 1, "turb"], type = "b", ylim = c(0, 100), pch = 19, xlab = "Distance(km)", ylab = "Turbidity (NTU)", main = "a) Turbidity plotted against spatial position", lwd = 3, lty = 1)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "turb"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "turb"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "turb"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "turb"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

tempo <- lm(log(turb) ~ dist + time + dist:time, data = alldata)
print(anova(tempo))
print(summary(tempo))
tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

plot(tempo1[tempo1$time == 1, "dist"], tempo1[tempo1$time == 1, "turb"], type = "p", pch = 19, col = "black", xlab = "Distance (km)", ylab = "Turbidity (NTU)", main = "b) Turbidity with spatiotemporal prediction", ylim = c(0, 100), cex = 1, lwd = 3, lty = 1)
points(tempo1[tempo1$time == 2, "dist"], tempo1[tempo1$time == 2, "turb"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1[tempo1$time == 3, "dist"], tempo1[tempo1$time == 3, "turb"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$time == 4, "dist"], tempo1[tempo1$time == 4, "turb"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$time == 5, "dist"], tempo1[tempo1$time == 5, "turb"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", col = "black", pch = 19, cex = 1.2, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

tempo <- lm(log(turb) ~ log(rain3) + log(dayflow) + log(vel), data = alldata)
print(anova(tempo))
print(summary(tempo))
tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "b", ylim = c(0, 100), pch = 19, col = "black", xlab = "Distance(km)", ylab = "Turbidity (NTU)", main = "c) Prediction with environmental parameters", lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "1"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "1"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "1"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "1"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()

cat("\nScript 14 completed successfully.\n")