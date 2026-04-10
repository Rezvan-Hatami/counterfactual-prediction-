# 18_other_structural_equations.R
# Author: Rezvan Hatami
# Date: 2026-03-27

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
script_figure_dir <- fs::path(figures_dir, "18")

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

# ----------------------------------------------------------------------
# 6. pH
# ----------------------------------------------------------------------
tempo <- lm(ph2 ~ time, data = alldata)

cat("\n============================================================\n")
cat("pH spatiotemporal model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(cbind(fitted(tempo), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "18_ph_spatiotemporal_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$time == 1, "dist"], tempo1[tempo1$time == 1, "ph2"], type = "p", pch = 19, col = "black", xlab = "Distance (km)", ylab = "pH", main = "pH with spatiotemporal prediction", ylim = c(6.5, 9), cex = 1, lwd = 3, lty = 1)
points(tempo1[tempo1$time == 2, "dist"], tempo1[tempo1$time == 2, "ph2"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1[tempo1$time == 3, "dist"], tempo1[tempo1$time == 3, "ph2"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$time == 4, "dist"], tempo1[tempo1$time == 4, "ph2"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$time == 5, "dist"], tempo1[tempo1$time == 5, "ph2"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
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

tempo <- lm(ph2 ~ rain2, data = alldata)

cat("\n============================================================\n")
cat("pH environmental model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(cbind(fitted(tempo), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "18_ph_environmental_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", ylim = c(6.5, 9), pch = 19, col = "black", xlab = "Distance(km)", ylab = "pH", main = "Prediction with environmental parameters", lwd = 3, lty = 1)
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

# ----------------------------------------------------------------------
# 7. Zinc
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ dzn + eff + dist + log(zn) + ph1,
  data = alldata
)

tempo <- lm(log(zn) ~ eff, data = alldata)
print_model_block("zn_eff", tempo)

tempo <- lm(log(zn) ~ dist, data = alldata)
print_model_block("zn_dist", tempo)

tempo <- lm(log(zn) ~ dzn, data = alldata)
print_model_block("zn_dzn", tempo)

tiff(
  file = fs::path(script_figure_dir, "18_figure_13_zinc_spatial.tif"),
  width = 10,
  height = 10,
  units = "in",
  pointsize = 12,
  bg = "transparent",
  res = 800,
  compression = "lzw"
)

par(mfrow = c(2, 1), mar = c(3, 4.5, 2.5, 3.5), cex = 1.5, cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8)

plot(alldata[alldata$day == 1, "dist"], alldata[alldata$day == 1, "zn"], type = "b", ylim = c(0, 120), pch = 19, xlab = "Distance(km)", ylab = "Zinc (mg/L)", main = "a) Zn plotted against spatial position", lwd = 3, lty = 1)
points(alldata[alldata$day == 126, "dist"], alldata[alldata$day == 126, "zn"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(alldata[alldata$day == 260, "dist"], alldata[alldata$day == 260, "zn"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(alldata[alldata$day == 336, "dist"], alldata[alldata$day == 336, "zn"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(alldata[alldata$day == 518, "dist"], alldata[alldata$day == 518, "zn"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

plot(alldata[alldata$day == 1, "dist"], alldata[alldata$day == 1, "zn1"], type = "b", ylim = c(0, 30), pch = 19, xlab = "Distance(km)", ylab = "Zinc (mg/L)", main = "b) Zn plotted against spatial position/n (without outliers)", lwd = 3, lty = 1)
points(alldata[alldata$day == 126, "dist"], alldata[alldata$day == 126, "zn1"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(alldata[alldata$day == 260, "dist"], alldata[alldata$day == 260, "zn1"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(alldata[alldata$day == 336, "dist"], alldata[alldata$day == 336, "zn1"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(alldata[alldata$day == 518, "dist"], alldata[alldata$day == 518, "zn1"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()

tempo <- lm(log(zn) ~ eff, data = alldata)

cat("\n============================================================\n")
cat("Zinc spatiotemporal model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_zinc_spatiotemporal_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_zinc_spatiotemporal_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_zinc_spatiotemporal_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "18_zinc_spatiotemporal_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$time == 1, "dist"], tempo1[tempo1$time == 1, "zn"], type = "p", pch = 19, col = "black", xlab = "Distance (km)", ylab = "Zn", main = "Zinc with spatiotemporal prediction", ylim = c(0, 20), cex = 1, lwd = 3, lty = 1)
points(tempo1[tempo1$time == 2, "dist"], tempo1[tempo1$time == 2, "zn"], col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$time == 3, "dist"], tempo1[tempo1$time == 3, "zn"], col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$time == 4, "dist"], tempo1[tempo1$time == 4, "zn"], col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$time == 5, "dist"], tempo1[tempo1$time == 5, "zn"], col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", col = "black", pch = 19, cex = 1.5, lwd = 3, lty = 1)
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

tempo <- lm(log(zn1) ~ eff, data = alldata)

cat("\n============================================================\n")
cat("Zinc environmental model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_zinc_environmental_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

tempo1 <- cbind(cbind(exp(fitted(tempo)), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

png(
  filename = fs::path(script_figure_dir, "18_zinc_environmental_prediction.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "1"], type = "l", ylim = c(0, 40), pch = 19, col = "black", xlab = "Distance(km)", ylab = "Zinc", main = "Prediction with environmental parameters", lwd = 3, lty = 1)
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

par(mar = c(5, 5, 4, 2), cex = 0.9)

png(
  filename = fs::path(script_figure_dir, "18_zinc_environmental_residuals_distance.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "2"], type = "b", ylim = c(-1, 20), pch = 19, cex = 1.5, xlab = "Distance(km)", ylab = "residuals", main = "prediction with space and time", lwd = 3, lty = 4)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "2"], type = "b", col = "gold2", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "2"], type = "b", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "2"], type = "b", col = "green3", pch = 17, cex = 1.5, lwd = 3)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "2"], type = "b", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
       lty = 1, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 5, horiz = FALSE, cex = 0.7, title = "months")
abline(v = 4, lty = 2)
dev.off()

print(xyplot(fitted(tempo) ~ dist, groups = day, data = tempo1, type = "l", auto.key = TRUE))
print(xyplot(residuals(tempo) ~ dist, groups = day, data = tempo1, type = "l", auto.key = TRUE))

tempo3 <- lm(tempo1[, 2] ~ tempo1$pred1)

cat("\n============================================================\n")
cat("tempo3\n")
cat("============================================================\n")
print(summary(tempo3))

tempo <- lm(tempo1[, 2] ~ as.factor(eff) + dist + time + as.factor(eff):time + as.factor(eff):dist, data = alldata)

cat("\n============================================================\n")
cat("Zinc residuals against time and space\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_zinc_residual_time_space_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

# ----------------------------------------------------------------------
# 8. Chlorophyll A in discharge
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ solar + airtemp + time + dchla,
  data = alldata
)

time1 <- 2 * pi * alldata$day / 365

tempo <- lm(dchla ~ solar, data = alldata)

cat("\n============================================================\n")
cat("Discharge chlorophyll A model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_dchla_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_dchla_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_dchla_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

# ----------------------------------------------------------------------
# 9. Air temperature
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ airtemp + solar + time,
  data = alldata
)

time1 <- 2 * pi * alldata$day / 365

cat("\n============================================================\n")
cat("time1\n")
cat("============================================================\n")
print(time1)

tempo <- lm(airtemp ~ time1, data = alldata)

cat("\n============================================================\n")
cat("Air temperature model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_airtemp_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_airtemp_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_airtemp_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

# ----------------------------------------------------------------------
# 10. Solar
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ solar + time,
  data = alldata
)

# ----------------------------------------------------------------------
# 11. Canopy cover
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ dist + canop + veg30m + log(veg30m),
  data = alldata
)

tempo <- lm(log(canop) ~ time + dist + time:dist, data = alldata)

cat("\n============================================================\n")
cat("Canopy spatiotemporal model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_canopy_spatiotemporal_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_canopy_spatiotemporal_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_canopy_spatiotemporal_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_canopy_spatial_position.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(bugenv[bugenv$day == 1, "dist"], bugenv[bugenv$day == 1, "canop"], type = "b", ylim = c(0, 100), pch = 19, xlab = "Distance(km)", ylab = "Alkalinity (mg/L)", main = "Canopy cover plotted against spatial position", lwd = 3, lty = 1)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "canop"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "canop"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "canop"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "canop"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)
dev.off()

tempo <- lm(canop ~ veg30m, data = alldata)

cat("\n============================================================\n")
cat("Canopy environmental model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_canopy_environmental_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_canopy_environmental_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_canopy_environmental_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

# ----------------------------------------------------------------------
# 12. Riparian vegetation
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ dist + veg30m,
  data = alldata
)

tempo <- lm(veg30m ~ dist, data = alldata)

cat("\n============================================================\n")
cat("Riparian vegetation model\n")
cat("============================================================\n")
print(summary(tempo))
print(anova(tempo))

png(
  filename = fs::path(script_figure_dir, "18_veg30m_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_veg30m_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_veg30m_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

# ----------------------------------------------------------------------
# 13. Creek flow rate
# ----------------------------------------------------------------------
par(mar = c(5, 5, 4, 2), cex = 0.9)

make_pairs_plot(
  ~ rain3 + dflow + dayflow,
  data = alldata
)

make_pairs_plot(
  ~ log(rain3) + log(rain2) + log(dflow) + log(dayflow),
  data = alldata
)

tempo <- lm(dayflow ~ rain3, data = alldata)

cat("\n============================================================\n")
cat("Simple creek flow regression\n")
cat("============================================================\n")
print(summary(tempo))

tempo <- lm(dayflow ~ dflow + rain3, data = alldata)

cat("\n============================================================\n")
cat("Creek flow rate model\n")
cat("============================================================\n")
print(anova(tempo))
print(summary(tempo))

png(
  filename = fs::path(script_figure_dir, "18_dayflow_diagnostics.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(tempo)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_dayflow_qqnorm.png"),
  width = 2400,
  height = 1800,
  res = 300
)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
dev.off()

png(
  filename = fs::path(script_figure_dir, "18_dayflow_residual_hist.png"),
  width = 2400,
  height = 1800,
  res = 300
)
hist(tempo$residuals)
dev.off()

tempo1 <- cbind(cbind(fitted(tempo), residuals(tempo)), alldata)
names(tempo1[, 1:2]) <- c("fitted", "residual")

if (exists("mydates") && exists("mydaterange")) {
  tiff(
    file = fs::path(script_figure_dir, "18_flow_rate_figure.tif"),
    width = 10,
    height = 4.6,
    units = "in",
    pointsize = 12,
    bg = "transparent",
    res = 800,
    compression = "lzw"
  )
  
  par(mfrow = c(1, 1), mar = c(3, 4.5, 2.5, 3.5), cex = 1.5, cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8)
  
  plot(mydates, bugenv[bugenv$dist == 0, "dayflow"], type = "b", ylim = c(10, 120), pch = 16, cex = 1.5, xlab = "Day", ylab = "Creek flow rate (ML/day)", main = "Observed and predicted values of \n creek flow rate against time ", lwd = 3, lty = 2, col = "red", xaxt = "n")
  axis.Date(side = 1, mydates, at = seq(mydaterange[1], mydaterange[2], by = "month"), format = "%b-%y")
  points(mydates, tempo1[tempo1$dist == 0, "1"], type = "l", ylim = c(0, 120), pch = 19, cex = 1.5, xlab = "day", ylab = "Creek flow rate (ML/day)", main = "Average daily Creek flow rate plotted against time", lwd = 3, lty = 1, xaxt = "n")
  legend("topright", inset = c(0, 0), legend = c("Observed values", "Predicted values"), lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2, col = c("red", "black"), ncol = 1, horiz = FALSE, cex = 0.8, title = "")
  dev.off()
}

# ----------------------------------------------------------------------
# 14. Additional checks
# ----------------------------------------------------------------------
tempo <- lm(res300 ~ dist, data = alldata)
print_model_block("res300_dist", tempo)

tempo <- lm(log(vel) ~ dist, data = alldata)
print_model_block("logvel_dist", tempo)

tempo <- lm(ph2 ~ rain2, data = alldata)
print_model_block("ph2_rain2", tempo)

par(mar = c(5, 5, 4, 2), cex = 4)

make_pairs_plot(
  ~ log(dalk) + log(dcond) + log(dtoc) + log(dno3) + log(dtp) +
    log(dchla) + log(dflow) + log(rain3) + log(rain2),
  data = alldata
)

cat("\nScript 18 completed successfully.\n")